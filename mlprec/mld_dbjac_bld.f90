!!$
!!$ 
!!$                                MLD2P4
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS v.2.0)
!!$  
!!$  (C) Copyright 2007  Alfredo Buttari      University of Rome Tor Vergata
!!$                           Pasqua D'Ambra       ICAR-CNR, Naples
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
! File: mld_dbjac_bld.f90.
!
! Subroutine: mld_dbjac_bld.
! Version:    real.
!
!  This routine builds the local extended matrix associated to an Additive
!  Schwarz preconditioner and computes an LU or incomplete LU factorization
!  of this matrix, according to the value of p%iprcparm(iprcparm(sub_solve_),
!  set by the user through mld_dprecinit or mld_dprecset.
!  Alternatively, it splits the local matrix into its block-diagonal and
!  off block-diagonal parts, for the future application of multiple
!  block-Jacobi sweeps.
!
!  This routine is used by mld_dbaseprec_bld, to build a 'base' block-Jacobi or
!  Additive Schwarz (AS) preconditioner at any level of a multilevel preconditioner,
!  or a block-Jacobi or LU or ILU solver at the coarsest level of a multilevel
!  preconditioner. More precisely, the routine is used to perform one of the
!  following tasks: 
!
!  1. construction of a block-Jacobi or Additive Schwarz preconditioner associated
!     to a matrix A distributed among the processes (allowed at any level);
!
!  2. setup of block-Jacobi sweeps to compute an approximate solution of a
!     linear system
!                                    A*Y = X,
!
!     distributed among the processes (allowed only at the coarsest level);
!
!  3. LU factorization of a linear system
!
!                                    A*Y = X,
!
!     distributed among the processes (allowed only at the coarsest level);
!
!  4. LU or incomplete LU factorization of a linear system
!
!                                    A*Y = X,
!
!        replicated on the processes (allowed only at the coarsest level).
!
!  The following factorizations are available:
!  - ILU(k), i.e. ILU factorization with fill-in level k;
!  - MILU(k), i.e. modified ILU factorization with fill-in level k;
!  - ILU(k,t), i.e. ILU with threshold (i.e. drop tolerance) t and k additional
!    entries in each row of the L and U factors with respect to the initial
!    sparsity pattern;
!  - LU implemented in SuperLU version 3.0;
!  - LU implemented in UMFPACK version 4.4;
!  - distributed LU implemented in SuperLU_DIST version 2.0.
!
!  
! Arguments:
!    a       -  type(psb_dspmat_type), input.
!               The sparse matrix structure containing the local part of the
!               matrix to be preconditioned or factorized.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor associated to a.
!    p       -  type(mld_dbaseprec_type), input/output.
!               The 'base preconditioner' data structure containing the local 
!               part of the preconditioner or solver at the current level.
!    info    -  integer, output.
!               Error code.              
!  
subroutine mld_dbjac_bld(a,desc_a,p,upd,info)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_dbjac_bld

  implicit none
                                                                               
! Arguments
  type(psb_dspmat_type), intent(in), target :: a
  type(psb_desc_type), intent(in)           :: desc_a
  type(mld_dbaseprc_type), intent(inout)    :: p
  integer, intent(out)                      :: info
  character, intent(in)                     :: upd

  !      Local Variables                         
  integer  ::    i, k, m
  integer  ::    int_err(5)
  character ::        trans, unitd
  type(psb_dspmat_type) :: blck, atmp
  logical, parameter :: debugprt=.false., debug=.false., aggr_dump=.false.
  integer :: err_act, n_row, nrow_a,n_col
  integer :: ictxt,np,me
  character(len=20)      :: name
  character(len=5), parameter :: coofmt='COO', csrfmt='CSR'

  if(psb_get_errstatus().ne.0) return 
  info=0
  name='mld_dbjac_bld'
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)

  m = a%m
  if (m < 0) then
    info = 10
    int_err(1) = 1
    int_err(2) = m
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif
  trans = 'N'
  unitd = 'U'
  if (p%iprcparm(mld_n_ovr_) < 0) then
    info = 11
    int_err(1) = 1
    int_err(2) = p%iprcparm(mld_n_ovr_)
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  call psb_nullify_sp(blck)
  call psb_nullify_sp(atmp)


  if(debug) write(0,*)me,': calling mld_asmat_bld',&
       & p%iprcparm(mld_prec_type_),p%iprcparm(mld_n_ovr_)
  if (debug) call psb_barrier(ictxt)

  !
  ! Build the communication descriptor for the Additive Schwarz
  ! preconditioner and retrieves the remote pieces of the local
  ! extended matrix needed by that preconditioner. If the
  ! preconditioner is the block-Jacobi one, only a copy of the
  ! descriptor of the original matrix is made.
  !
  call mld_asmat_bld(p%iprcparm(mld_prec_type_),p%iprcparm(mld_n_ovr_),a,&
       & blck,desc_a,upd,p%desc_data,info,outfmt=csrfmt)

  if (debugprt) then 
    open(60+me)
    call psb_csprt(60+me,a,head='% A')
    close(60+me)
    open(70+me)
    call psb_csprt(70+me,blck,head='% BLCK')
    close(70+me)
  endif
  
  if(info/=0) then
    call psb_errpush(4010,name,a_err='mld_asmat_bld')
    goto 9999
  end if

  if (debug) write(0,*)me,': out of mld_asmat_bld'
  if (debug) call psb_barrier(ictxt)

  !
  ! Treat separately the case the local matrix has to be reordered
  ! and the case this is not required.
  !
  select case(p%iprcparm(mld_sub_ren_)) 

  !
  ! A reordering of the local matrix is required.
  !
  case (1:)

    !
    ! Reorder the rows and the columns of the local extended matrix,
    ! according to the value of p%iprcparm(sub_ren_). The reordered 
    ! matrix is stored into atmp, using the COO format.
    !
    call  mld_sp_renum(a,desc_a,blck,p,atmp,info)

    if (info/=0) then
      call psb_errpush(4010,name,a_err='mld_sp_renum')
      goto 9999
    end if

    !
    ! Clip into p%av(ap_nd_) the off block-diagonal part of the local
    ! matrix. The clipped matrix is then stored in CSR format.
    !
    call psb_sp_clip(atmp,p%av(mld_ap_nd_),info,&
         & jmin=atmp%m+1,rscale=.false.,cscale=.false.)

    call psb_spcnv(p%av(mld_ap_nd_),info,afmt='csr',dupl=psb_dupl_add_)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psb_spcnv csr 1')
      goto 9999
    end if

    k = psb_sp_get_nnzeros(p%av(mld_ap_nd_))
    call psb_sum(ictxt,k)

    if (k == 0) then
      ! 
      ! If the off diagonal part is emtpy, there is no point in doing
      ! multiple Jacobi sweeps. This is certain to happen when running
      ! on a single processor.
      !
      p%iprcparm(mld_smooth_sweeps_) = 1
    end if


    if (debugprt) then 
      call psb_barrier(ictxt)
      open(40+me) 
      call psb_csprt(40+me,atmp,head='% Local matrix')
      close(40+me)
    endif
    if (debug) write(0,*) me,' Factoring rows ',&
         &atmp%m,a%m,blck%m,atmp%ia2(atmp%m+1)-1

    ! 
    ! Compute a factorization of the diagonal block of the local matrix,
    ! according to the choice made by the user by setting p%iprcparm(sub_solve_)
    !
    select case(p%iprcparm(mld_sub_solve_))

    case(mld_ilu_n_,mld_milu_n_,mld_ilu_t_) 
    !
    ! ILU(k)/MILU(k)/ILU(k,t) factorization.
    !

      call psb_spcnv(atmp,info,afmt='csr',dupl=psb_dupl_add_)
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='psb_spcnv csr 2')
        goto 9999
      end if

      call mld_ilu_bld(atmp,p%desc_data,p,upd,info)

      if (info/=0) then
        call psb_errpush(4010,name,a_err='mld_ilu_bld')
        goto 9999
      end if


      if (debugprt) then 

        open(80+me)

        call psb_csprt(80+me,p%av(mld_l_pr_),head='% Local L factor')
        write(80+me,*) '% Diagonal: ',p%av(mld_l_pr_)%m
        do i=1,p%av(mld_l_pr_)%m
          write(80+me,*) i,i,p%d(i)
        enddo
        call psb_csprt(80+me,p%av(mld_u_pr_),head='% Local U factor')

        close(80+me)
      endif

    case(mld_slu_)
    !
    ! LU factorization through the SuperLU package.
    !

      call psb_spcnv(atmp,info,afmt='csr',dupl=psb_dupl_add_)
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='psb_spcnv csr 3')
        goto 9999
      end if

      call mld_slu_bld(atmp,p%desc_data,p,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='mld_slu_bld')
        goto 9999
      end if

    case(mld_umf_)
    !
    ! LU factorization through the UMFPACK package.
    !

      call psb_spcnv(atmp,info,afmt='csc',dupl=psb_dupl_add_)
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='psb_spcnv csc')
        goto 9999
      end if

      call mld_umf_bld(atmp,p%desc_data,p,info)
      if(debug) write(0,*)me,': Done mld_umf_bld ',info
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='mld_umf_bld')
        goto 9999
      end if

    case(mld_f_none_) 
    !
    ! Error: no factorization required.
    !
      info=4010
      call psb_errpush(info,name,a_err='Inconsistent prec  mld_f_none_')
      goto 9999

    case default
      info=4010
      call psb_errpush(info,name,a_err='Unknown mld_sub_solve_')
      goto 9999
    end select

    call psb_sp_free(atmp,info) 

    if(info/=0) then
      call psb_errpush(4010,name,a_err='psb_sp_free')
      goto 9999
    end if

  !
  ! No reordering of the local matrix is required
  !
  case(0)

    ! 
    ! Compute a factorization of the diagonal block of the local matrix,
    ! according to the choice made by the user by setting p%iprcparm(sub_solve_)
    !
    select case(p%iprcparm(mld_sub_solve_))


    case(mld_ilu_n_,mld_milu_n_,mld_ilu_t_) 
    !
    ! ILU(k)/MILU(k)/ILU(k,t) factorization.
    !

      !
      ! In case of multiple block-Jacobi sweeps, clip into p%av(ap_nd_)
      ! the off block-diagonal part of the local extended matrix. The
      ! clipped matrix is then stored in CSR format.
      !

      if (p%iprcparm(mld_smooth_sweeps_) > 1) then 
        n_row = psb_cd_get_local_rows(p%desc_data)
        n_col = psb_cd_get_local_cols(p%desc_data)
        nrow_a = a%m 
        ! The following is known to work 
        ! given that the output from CLIP is in COO. 
        call psb_sp_clip(a,p%av(mld_ap_nd_),info,&
             & jmin=nrow_a+1,rscale=.false.,cscale=.false.)
        call psb_sp_clip(blck,atmp,info,&
             & jmin=nrow_a+1,rscale=.false.,cscale=.false.)
        call psb_rwextd(n_row,p%av(mld_ap_nd_),info,b=atmp) 

        call psb_spcnv(p%av(mld_ap_nd_),info,afmt='csr',dupl=psb_dupl_add_)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_spcnv csr 4')
          goto 9999
        end if
        
        k = psb_sp_get_nnzeros(p%av(mld_ap_nd_))
        call psb_sum(ictxt,k)

        if (k == 0) then 
          !
          ! If the off block-diagonal part is emtpy, there is no point in doing
          ! multiple Jacobi sweeps. This is certain to happen when running
          ! on a single processor.
          !
          p%iprcparm(mld_smooth_sweeps_) = 1
        end if
        call psb_sp_free(atmp,info) 
      end if

      !
      ! Compute the incomplete LU factorization.
      !
      call mld_ilu_bld(a,desc_a,p,upd,info,blck=blck)

      if(info/=0) then
        call psb_errpush(4010,name,a_err='mld_ilu_bld')
        goto 9999
      end if

      if (debugprt) then 

        open(80+me)

        call psb_csprt(80+me,p%av(mld_l_pr_),head='% Local L factor')
        write(80+me,*) '% Diagonal: ',p%av(mld_l_pr_)%m
        do i=1,p%av(mld_l_pr_)%m
          write(80+me,*) i,i,p%d(i)
        enddo
        call psb_csprt(80+me,p%av(mld_u_pr_),head='% Local U factor')

        close(80+me)
      endif

    case(mld_slu_)
    !
    ! LU factorization through the SuperLU package.
    !

      call psb_spcnv(a,atmp,info,afmt='coo')
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='psb_spcnv')
        goto 9999
      end if

      n_row = psb_cd_get_local_rows(p%desc_data)
      n_col = psb_cd_get_local_cols(p%desc_data)
      call psb_rwextd(n_row,atmp,info,b=blck) 

      !
      ! In case of multiple block-Jacobi sweeps, clip into p%av(ap_nd_)
      ! the off block-diagonal part of the local extended matrix. The
      ! clipped matrix is then stored in CSR format.
      !
      if (p%iprcparm(mld_smooth_sweeps_) > 1) then 

        call psb_sp_clip(atmp,p%av(mld_ap_nd_),info,&
             & jmin=atmp%m+1,rscale=.false.,cscale=.false.)

        call psb_spcnv(p%av(mld_ap_nd_),info,afmt='csr',dupl=psb_dupl_add_)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_spcnv csr 6')
          goto 9999
        end if

        k = psb_sp_get_nnzeros(p%av(mld_ap_nd_))
        call psb_sum(ictxt,k)

        if (k == 0) then 
          !
          ! If the off block-diagonal part is emtpy, there is no point in doing
          ! multiple Jacobi sweeps. This is certain to happen when running
          ! on a single processor.
          !
          p%iprcparm(mld_smooth_sweeps_) = 1
        end if
      endif

      !
      ! Compute the LU factorization.
      !
      if (info == 0) call psb_spcnv(atmp,info,afmt='csr',dupl=psb_dupl_add_)
      if (info == 0) call mld_slu_bld(atmp,p%desc_data,p,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='mld_slu_bld')
        goto 9999
      end if

      call psb_sp_free(atmp,info) 
      if(info/=0) then
        call psb_errpush(4010,name,a_err='psb_sp_free')
        goto 9999
      end if

    case(mld_sludist_)
      !
      ! LU factorization through the SuperLU_DIST package. This works only
      ! when the matrix is distributed among the processes.
      ! NOTE: Should have NO overlap here!!!!   
      !

      call psb_spcnv(a,atmp,info,afmt='csr')
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='psb_spcnv')
        goto 9999
      end if
      
      n_row = psb_cd_get_local_rows(p%desc_data)
      n_col = psb_cd_get_local_cols(p%desc_data)
      if (info == 0) call mld_sludist_bld(atmp,p%desc_data,p,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='mld_slu_bld')
        goto 9999
      end if

      call psb_sp_free(atmp,info) 
      if(info/=0) then
        call psb_errpush(4010,name,a_err='psb_sp_free')
        goto 9999
      end if

    case(mld_umf_)
      !
      ! LU factorization through the UMFPACK package.
      !


      call psb_spcnv(a,atmp,info,afmt='coo')
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='psb_spcnv')
        goto 9999
      end if

      n_row = psb_cd_get_local_rows(p%desc_data)
      n_col = psb_cd_get_local_cols(p%desc_data)
      call psb_rwextd(n_row,atmp,info,b=blck) 

      !
      ! In case of multiple block-Jacobi sweeps, clip into p%av(ap_nd_)
      ! the off block-diagonal part of the local extended matrix. The
      ! clipped matrix is then stored in CSR format.
      !
      if (p%iprcparm(mld_smooth_sweeps_) > 1) then 

        call psb_sp_clip(atmp,p%av(mld_ap_nd_),info,&
             & jmin=atmp%m+1,rscale=.false.,cscale=.false.)

        call psb_spcnv(p%av(mld_ap_nd_),info,afmt='csr',dupl=psb_dupl_add_)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_spcnv csr 8')
          goto 9999
        end if

        k = psb_sp_get_nnzeros(p%av(mld_ap_nd_))
        call psb_sum(ictxt,k)

        if (k == 0) then 
          !
          ! If the off block-diagonal part is emtpy, there is no point in doing
          ! multiple Jacobi sweeps. This is certain to happen when running
          ! on a single processor.
          !
          p%iprcparm(mld_smooth_sweeps_) = 1
        end if
      endif

      !
      ! Compute the LU factorization.
      !
      if (info == 0) call psb_ipcoo2csc(atmp,info,clshr=.true.)
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='psb_ipcoo2csc')
        goto 9999
      end if
      call mld_umf_bld(atmp,p%desc_data,p,info)
      if(debug) write(0,*)me,': Done mld_umf_bld ',info
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='mld_umf_bld')
        goto 9999
      end if

      call psb_sp_free(atmp,info) 
      if(info/=0) then
        call psb_errpush(4010,name,a_err='psb_sp_free')
        goto 9999
      end if


    case(mld_f_none_) 
      !
      ! Error: no factorization required.
      !
      info=4010
      call psb_errpush(info,name,a_err='Inconsistent prec  mld_f_none_')
      goto 9999

    case default
      info=4010
      call psb_errpush(info,name,a_err='Unknown mld_sub_solve_')
      goto 9999
    end select

  case default
    info=4010
    call psb_errpush(info,name,a_err='Invalid renum_')
    goto 9999

  end select


  call psb_sp_free(blck,info)
  if(info/=0) then
    call psb_errpush(4010,name,a_err='psb_sp_free')
    goto 9999
  end if

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


end subroutine mld_dbjac_bld


