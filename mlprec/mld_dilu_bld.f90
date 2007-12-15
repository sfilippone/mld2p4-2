!!$ 
!!$ 
!!$                                MLD2P4
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS v.2.0)
!!$  
!!$  (C) Copyright 2007  Alfredo Buttari      University of Rome Tor Vergata
!!$                      Pasqua D'Ambra       ICAR-CNR, Naples
!!$                      Daniela di Serafino  Second University of Naples
!!$                      Salvatore Filippone  University of Rome Tor Vergata       
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the MLD2P4 group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$
! File: mld_dilu_bld.f90.
!
! Subroutine: mld_dilu_bld.
! Version:    real.
!
!  This routine computes an incomplete LU (ILU) factorization of the local part
!  of the matrix stored into a. This factorization is used to build the 'base
!  preconditioner' (block-Jacobi preconditioner/solver, Additive Schwarz
!  preconditioner) corresponding to a certain level of a multilevel preconditioner.
!  The following factorizations are available:
!  - ILU(k), i.e. ILU factorization with fill-in level k,
!  - MILU(k), i.e. modified ILU factorization with fill-in level k,
!  - ILU(k,t), i.e. ILU with threshold (i.e. drop tolerance) t and k additional
!    entries in each row of the L and U factors with respect to the initial
!    sparsity pattern.
!  Note that the meaning of k in ILU(k,t) is different from that in ILU(k) and
!  MILU(k).
!
!  For details on the above factorizations see
!    Y. Saad, Iterative Methods for Sparse Linear Systems, Second Edition,
!    SIAM, 2003, Chapter 10.
!
!  Note that that this routine handles the ILU(0) factorization separately,
!  through mld_ilu_fct, for performance reasons.
!
!
! Arguments:
!    a       -  type(psb_dspmat_type), input.
!               The sparse matrix structure containing the local matrix to be
!               factorized. Note that if p%iprcparm(mld_n_ovr_) > 0, i.e. the
!               'base' Additive Schwarz preconditioner has overlap greater than
!               0, and p%iprcparm(mld_sub_ren_) = 0, i.e. a      reordering of the
!               matrix has not been performed (see mld_bjac_bld), then a contains
!               only the 'original' local part of the matrix to      be factorized,
!               i.e. the rows of the matrix held by the calling process according
!               to the initial data distribution.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor associated to a.
!    p       -  type(mld_dbaseprc_type), input/output.
!               The 'base preconditioner' data structure. In input, p%iprcparm
!               contains information on the type of factorization to be computed.
!               In output, p%av(mld_l_pr_) and p%av(mld_u_pr_) contain the
!               incomplete L and U factors (without their diagonals), and p%d
!               contains the diagonal of the incomplete U factor. For more
!               details on p see its description in mld_prec_type.f90.
!    info    -  integer, output.                                                             
!               Error code.
!    blck    -  type(psb_dspmat_type), input, optional.
!               The sparse matrix structure containing the remote rows of the
!               matrix to be factorized, that have been retrieved by mld_asmat_bld
!               to build an Additive Schwarz base preconditioner with overlap
!               greater than 0.      If the overlap is 0 or the matrix has been reordered
!               (see mld_bjac_bld), then blck does not contain any row.
!
subroutine mld_dilu_bld(a,desc_a,p,upd,info,blck)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_dilu_bld

  implicit none
                                                      
! Arguments                                                     
  type(psb_dspmat_type), intent(in), target   :: a
  type(mld_dbaseprc_type), intent(inout)      :: p
  type(psb_desc_type), intent(in)             :: desc_a
  character, intent(in)                       :: upd
  integer, intent(out)                        :: info
  type(psb_dspmat_type), intent(in), optional :: blck

  !     Local Variables                       
  integer   :: i, nztota, err_act, n_row, nrow_a
  character :: trans, unitd
  logical, parameter :: debugprt=.false., debug=.false., aggr_dump=.false.
  integer :: ictxt,np,me
  character(len=20)  :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  name='mld_dilu_bld'
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)

  trans = 'N'
  unitd = 'U'

  !
  ! Check the memory available to hold the incomplete L and U factors
  ! and allocate it if needed
  !

  if (allocated(p%av)) then 
    if (size(p%av) < mld_bp_ilu_avsz_) then 
      do i=1,size(p%av) 
        call psb_sp_free(p%av(i),info)
        if (info /= 0) then 
          ! Actually, we don't care here about this. Just let it go.
          ! return
        end if
      enddo
      deallocate(p%av,stat=info)
    endif
  end if
  if (.not.allocated(p%av)) then 
    allocate(p%av(mld_max_avsz_),stat=info)
    if (info /= 0) then
      call psb_errpush(4000,name)
      goto 9999
    end if
  endif
!!$  call psb_csprt(50+me,a,head='% (A)')    

  nrow_a = psb_cd_get_local_rows(desc_a)
  nztota = psb_sp_get_nnzeros(a)
  if (present(blck)) then 
    nztota = nztota + psb_sp_get_nnzeros(blck)
  end if
  if (debug) write(0,*)me,': out get_nnzeros',nztota,a%m,a%k
  if (debug) call psb_barrier(ictxt)

  n_row  = p%desc_data%matrix_data(psb_n_row_)
  p%av(mld_l_pr_)%m  = n_row
  p%av(mld_l_pr_)%k  = n_row
  p%av(mld_u_pr_)%m  = n_row
  p%av(mld_u_pr_)%k  = n_row
  call psb_sp_all(n_row,n_row,p%av(mld_l_pr_),nztota,info)
  if (info == 0) call psb_sp_all(n_row,n_row,p%av(mld_u_pr_),nztota,info)
  if(info/=0) then
    info=4010
    ch_err='psb_sp_all'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (allocated(p%d)) then 
    if (size(p%d) < n_row) then 
      deallocate(p%d)
    endif
  endif
  if (.not.allocated(p%d)) then 
    allocate(p%d(n_row),stat=info)
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

  endif

  select case(p%iprcparm(mld_sub_solve_))

  case (mld_ilu_t_)
  !
  ! ILU(k,t)
  !

    select case(p%iprcparm(mld_sub_fill_in_))

    case(:-1) 
      ! Error: fill-in <= -1
      call psb_errpush(30,name,i_err=(/3,p%iprcparm(mld_sub_fill_in_),0,0,0/))
      goto 9999

    case(0:)
      ! Fill-in >= 0
      call mld_ilut_fct(p%iprcparm(mld_sub_fill_in_),p%dprcparm(mld_fact_thrs_),&
           & p%iprcparm(mld_sub_solve_), a,&
           & p%av(mld_l_pr_),p%av(mld_u_pr_),p%d,info,blck=blck)
    end select
    if(info/=0) then
      info=4010
      ch_err='mld_ilut_fct'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  case(mld_ilu_n_,mld_milu_n_) 
    !
    ! ILU(k) and MILU(k)
    !
    select case(p%iprcparm(mld_sub_fill_in_))
    case(:-1) 
      ! Error: fill-in <= -1
      call psb_errpush(30,name,i_err=(/3,p%iprcparm(mld_sub_fill_in_),0,0,0/))
      goto 9999
    case(0)
      ! Fill-in 0
      
      ! Separate implementation of ILU(0) for better performance
      ! There seems to be a problem with the separate implementation of MILU(0),
      ! contained into mld_ilu_fct. This must be investigated. For the time being,
      ! resort to the implementation of MILU(k)      with k=0.
      if (p%iprcparm(mld_sub_solve_) == mld_ilu_n_) then 
        call mld_ilu_fct(p%iprcparm(mld_sub_solve_),a,p%av(mld_l_pr_),p%av(mld_u_pr_),&
             & p%d,info,blck=blck)
      else
        call mld_iluk_fct(p%iprcparm(mld_sub_fill_in_),p%iprcparm(mld_sub_solve_),&
             & a,p%av(mld_l_pr_),p%av(mld_u_pr_),p%d,info,blck=blck)
      endif
    case(1:)
      ! Fill-in >= 1
      ! The same routine implements both ILU(k) and MILU(k)
      call mld_iluk_fct(p%iprcparm(mld_sub_fill_in_),p%iprcparm(mld_sub_solve_),&
           & a,p%av(mld_l_pr_),p%av(mld_u_pr_),p%d,info,blck=blck)
    end select
    if (info/=0) then
      info=4010
      ch_err='mld_iluk_fct'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  case default
    ! If we end up here, something was wrong up in the call chain. 
    call psb_errpush(4000,name)
    goto 9999

  end select

  if (debugprt) then 
    !
    ! Print out the factors on file.
    !
    open(80+me)

    call psb_csprt(80+me,p%av(mld_l_pr_),head='% Local L factor')
    write(80+me,*) '% Diagonal: ',p%av(mld_l_pr_)%m
    do i=1,p%av(mld_l_pr_)%m
      write(80+me,*) i,i,p%d(i)
    enddo
    call psb_csprt(80+me,p%av(mld_u_pr_),head='% Local U factor')

    close(80+me)
  endif

  if (psb_sp_getifld(psb_upd_,p%av(mld_u_pr_),info) /= psb_upd_perm_) then
    call psb_sp_trim(p%av(mld_u_pr_),info)
  endif

  if (psb_sp_getifld(psb_upd_,p%av(mld_l_pr_),info) /= psb_upd_perm_) then
    call psb_sp_trim(p%av(mld_l_pr_),info)
  endif

  if (debug) write(0,*) me,'End of ilu_bld'

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_dilu_bld


