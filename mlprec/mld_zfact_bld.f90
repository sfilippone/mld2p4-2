!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010
!!$
!!$                      Salvatore Filippone  University of Rome Tor Vergata
!!$                      Alfredo Buttari      CNRS-IRIT, Toulouse
!!$                      Pasqua D'Ambra       ICAR-CNR, Naples
!!$                      Daniela di Serafino  Second University of Naples
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
! File: mld_zfact_bld.f90
!
! Subroutine: mld_zfact_bld
! Version:    complex
!
!  This routine computes an LU or incomplete LU (ILU) factorization of the
!  diagonal blocks of a distributed matrix, according to the value of
!  p%iprcparm(iprcparm(sub_solve_), set by the user through	mld_zprecinit
!  or mld_zprecset.
!  It may also compute an LU factorization of a distributed matrix, or split
!  a distributed matrix into its block-diagonal and off block-diagonal parts, 
!  for the future application of multiple block-Jacobi sweeps.
!
!  This routine is used by mld_as_bld, to build a 'base' block-Jacobi or
!  Additive Schwarz (AS) preconditioner at any level of a multilevel preconditioner,
!  or a block-Jacobi or LU or ILU solver at the coarsest level of a multilevel
!  preconditioner. For the AS preconditioners, the diagonal blocks to be factorized
!  are stored into the sparse matrix data structures a and blck, and blck contains
!  the remote rows needed to build the extended local matrix as required by the
!  AS preconditioner.
!
!  More precisely, the routine performs one of the following tasks: 
!
!  1. LU or ILU factorization of the diagonal blocks of the distributed matrix
!     for the construction of block-Jacobi or AS preconditioners (allowed at
!     any level of a multilevel preconditioner);
!
!  2. setup of block-Jacobi sweeps to compute an approximate solution of a
!     linear system
!                                    A*Y = X,
!     distributed among the processes (allowed only at the coarsest level);
!
!  3. LU factorization of the matrix of a linear system
!                                    A*Y = X,
!     distributed among the processes (allowed only at the coarsest level);
!
!  4. LU or incomplete LU factorization of the matrix of a linear system
!                                    A*Y = X,
!     replicated on the processes (allowed only at the coarsest level).
!
!  The following factorizations are available:
!  - ILU(k), i.e. ILU factorization with fill-in level k;
!  - MILU(k), i.e. modified ILU factorization with fill-in level k;
!  - ILU(k,t), i.e. ILU with threshold (i.e. drop tolerance) t and k additional
!    entries in each row of the L and U factors with respect to the initial
!    sparsity pattern;
!  - serial LU implemented in SuperLU version 3.0;
!  - serial LU implemented in UMFPACK version 4.4;
!  - distributed LU implemented in SuperLU_DIST version 2.0.
!
!  
! Arguments:
!    a       -  type(psb_zspmat_type), input.
!               The sparse matrix structure containing the local part of the
!               distributed matrix.
!    p       -  type(mld_zbaseprec_type), input/output.
!               The 'base preconditioner' data structure containing the local 
!               part of the preconditioner or solver at the current level.
!
!    info    -  integer, output.
!               Error code.              
!    upd     -  character, input.
!               If upd='F' then the preconditioner is built from scratch;
!               if upd=T' then the matrix to be preconditioned has the same
!               sparsity pattern of a matrix that has been previously
!               preconditioned, hence some information is reused in building
!               the new preconditioner.
!    blck    -  type(psb_zspmat_type), input, optional.
!               The sparse matrix structure containing the remote rows of the
!               distributed matrix, that have been retrieved by mld_as_bld
!               to build an Additive Schwarz base preconditioner with overlap
!               greater than 0. If the overlap is 0 blck is empty.
!  
subroutine mld_zfact_bld(a,p,upd,info,blck)

  use psb_sparse_mod
  use mld_z_inner_mod, mld_protect_name => mld_zfact_bld

  implicit none
                                                                               
! Arguments
  type(psb_zspmat_type), intent(in), target :: a
  type(mld_zbaseprec_type), intent(inout)    :: p
  integer, intent(out)                      :: info
  character, intent(in)                     :: upd
  type(psb_zspmat_type), intent(in), target, optional  :: blck

  !      Local Variables                         
  type(psb_zspmat_type), pointer :: blck_
  type(psb_zspmat_type)          :: atmp
  integer                        :: ictxt,np,me,err_act
  integer                        :: debug_level, debug_unit
  integer                        :: k, m, int_err(5), n_row, nrow_a, n_col
  character                      :: trans, unitd
  character(len=20)              :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=psb_success_
  name='mld_zfact_bld'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ictxt       = psb_cd_get_context(p%desc_data)
  call psb_info(ictxt, me, np)

  m = a%m
  if (m < 0) then
    info = psb_err_iarg_neg_
    int_err(1) = 1
    int_err(2) = m
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif
  trans = 'N'
  unitd = 'U'

  if (present(blck)) then 
    blck_ => blck
  else
    allocate(blck_,stat=info)
    if (info == psb_success_) call psb_sp_all(0,0,blck_,1,info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    blck_%fida            = 'COO'
    blck_%infoa(psb_nnz_) = 0
  end if
  call psb_nullify_sp(atmp)

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
    call  mld_sp_renum(a,blck_,p,atmp,info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_sp_renum')
      goto 9999
    end if

    !
    ! Clip into p%av(ap_nd_) the off block-diagonal part of the local
    ! matrix. The clipped matrix is then stored in CSR format.
    !
    if (p%iprcparm(mld_smoother_sweeps_) > 1) then 
      call psb_sp_clip(atmp,p%av(mld_ap_nd_),info,&
           & jmin=atmp%m+1,rscale=.false.,cscale=.false.)
      if (info == psb_success_) call psb_spcnv(p%av(mld_ap_nd_),info,&
           & afmt='csr',dupl=psb_dupl_add_)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_spcnv')
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
        p%iprcparm(mld_smoother_sweeps_) = 1
      end if
    end if
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' Factoring rows ',&
         & atmp%m,a%m,blck_%m,atmp%ia2(atmp%m+1)-1

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
      if (info == psb_success_) call mld_ilu_bld(atmp,p,upd,info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_ilu_bld')
        goto 9999
      end if

    case(mld_slu_)
      !
      ! LU factorization through the SuperLU package.
      !
      call psb_spcnv(atmp,info,afmt='csr',dupl=psb_dupl_add_)
      if (info == psb_success_) call mld_slu_bld(atmp,p%desc_data,p,info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_slu_bld')
        goto 9999
      end if

    case(mld_sludist_)
      !
      ! LU factorization through the SuperLU_DIST package. This works only
      ! when the matrix is distributed among the processes.
      ! NOTE: Should have NO overlap here!!!!   
      !
      call psb_spcnv(a,atmp,info,afmt='csr')
      if (info == psb_success_) call mld_sludist_bld(atmp,p%desc_data,p,info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_sludist_bld')
        goto 9999
      end if

    case(mld_umf_)
      !
      ! LU factorization through the UMFPACK package.
      !
      call psb_spcnv(atmp,info,afmt='csc',dupl=psb_dupl_add_)
      if (info == psb_success_) call mld_umf_bld(atmp,p%desc_data,p,info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_umf_bld')
        goto 9999
      end if

    case(mld_f_none_) 
      !
      ! Error: no factorization required.
      !
      info=psb_err_internal_error_
      call psb_errpush(info,name,a_err='Inconsistent prec  mld_f_none_')
      goto 9999

    case default
      info=psb_err_internal_error_
      call psb_errpush(info,name,a_err='Unknown mld_sub_solve_')
      goto 9999
    end select

    call psb_sp_free(atmp,info) 

    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_sp_free')
      goto 9999
    end if

    !
    ! No reordering of the local matrix is required
    !
  case(0)
    !
    ! In case of multiple block-Jacobi sweeps, clip into p%av(ap_nd_)
    ! the off block-diagonal part of the local extended matrix. The
    ! clipped matrix is then stored in CSR format.
    !
    
    if (p%iprcparm(mld_smoother_sweeps_) > 1) then 
      n_row = psb_cd_get_local_rows(p%desc_data)
      n_col = psb_cd_get_local_cols(p%desc_data)
      nrow_a = a%m 
      ! The following is known to work 
      ! given that the output from CLIP is in COO. 
      call psb_sp_clip(a,p%av(mld_ap_nd_),info,&
           & jmin=nrow_a+1,rscale=.false.,cscale=.false.)
      if (info == psb_success_) call psb_sp_clip(blck_,atmp,info,&
           & jmin=nrow_a+1,rscale=.false.,cscale=.false.)
      if (info == psb_success_) call psb_rwextd(n_row,p%av(mld_ap_nd_),info,b=atmp) 
      if (info == psb_success_) call psb_spcnv(p%av(mld_ap_nd_),info,&
           & afmt='csr',dupl=psb_dupl_add_)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='clip & psb_spcnv csr 4')
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
        p%iprcparm(mld_smoother_sweeps_) = 1
      end if
      call psb_sp_free(atmp,info) 
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_sp_free')
        goto 9999
      end if
    end if
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
      ! Compute the incomplete LU factorization.
      !
      call mld_ilu_bld(a,p,upd,info,blck=blck_)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_ilu_bld')
        goto 9999
      end if

    case(mld_slu_)
      !
      ! LU factorization through the SuperLU package.
      ! 
      n_row = psb_cd_get_local_rows(p%desc_data)
      n_col = psb_cd_get_local_cols(p%desc_data)
      call psb_spcnv(a,atmp,info,afmt='coo')
      if (info == psb_success_) call psb_rwextd(n_row,atmp,info,b=blck_) 

      !
      ! Compute the LU factorization.
      !
      if (info == psb_success_) call psb_spcnv(atmp,info,afmt='csr',dupl=psb_dupl_add_)
      if (info == psb_success_) call mld_slu_bld(atmp,p%desc_data,p,info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_slu_bld')
        goto 9999
      end if

      call psb_sp_free(atmp,info) 
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_sp_free')
        goto 9999
      end if

    case(mld_sludist_)
      !
      ! LU factorization through the SuperLU_DIST package. This works only
      ! when the matrix is distributed among the processes.
      ! NOTE: Should have NO overlap here!!!!   
      !
      call psb_spcnv(a,atmp,info,afmt='csr')
      if (info == psb_success_) call mld_sludist_bld(atmp,p%desc_data,p,info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_sludist_bld')
        goto 9999
      end if

      call psb_sp_free(atmp,info) 
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_sp_free')
        goto 9999
      end if

    case(mld_umf_)
      !
      ! LU factorization through the UMFPACK package.
      !

      call psb_spcnv(a,atmp,info,afmt='coo')
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_spcnv')
        goto 9999
      end if

      n_row = psb_cd_get_local_rows(p%desc_data)
      n_col = psb_cd_get_local_cols(p%desc_data)
      call psb_rwextd(n_row,atmp,info,b=blck_) 

      !
      ! Compute the LU factorization.
      !
      if (info == psb_success_) call psb_spcnv(atmp,info,afmt='csc',dupl=psb_dupl_add_)
      if (info == psb_success_) call mld_umf_bld(atmp,p%desc_data,p,info)
      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & ': Done mld_umf_bld ',info
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_umf_bld')
        goto 9999
      end if

      call psb_sp_free(atmp,info) 
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_sp_free')
        goto 9999
      end if

    case(mld_f_none_) 
      !
      ! Error: no factorization required.
      !
      info=psb_err_internal_error_
      call psb_errpush(info,name,a_err='Inconsistent prec  mld_f_none_')
      goto 9999

    case default
      info=psb_err_internal_error_
      call psb_errpush(info,name,a_err='Unknown mld_sub_solve_')
      goto 9999
    end select

  case default
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Invalid renum_')
    goto 9999
  end select
  
  if (.not.present(blck)) then 
    call psb_sp_free(blck_,info)
    if (info == psb_success_) deallocate(blck_)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_sp_free')
      goto 9999
    end if
  end if

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),'End '

  call psb_erractionrestore(err_act)

  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return


end subroutine mld_zfact_bld


