!!$
!!$ 
!!$                           MLD2P4  version 1.1
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 2.3.1)
!!$  
!!$  (C) Copyright 2008,2009
!!$
!!$                      Salvatore Filippone  University of Rome Tor Vergata
!!$                      Alfredo Buttari      University of Rome Tor Vergata
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
! File mld_zsub_solve.f90
!
! Subroutine: mld_zsub_solve
! Version: complex
!
!  This routine computes
!
!                       Y = beta*Y + alpha*op(K^(-1))*X,
!
!  where
!  - K is a factored matrix, as specified below,
!  - op(K^(-1)) is K^(-1) or its transpose, according to the value of the
!    argument trans,
!  - X and Y are vectors,
!  - alpha and beta are scalars.
!
!  Depending on K, alpha and beta (and on the communication descriptor desc_data
!  - see the arguments below), the above computation may correspond to one of
!  the following tasks:
!
!  1. approximate solution of a linear system
!
!                                    A*Y = X,
!                              
!     by using the L and U factors computed with an ILU factorization of A.
!     In this case K = L*U ~ A, alpha = 1 and beta = 0. The factors L and U
!     (and the matrix A) are either distributed and block-diagonal or replicated. 
!
!  2. Solution of a linear system
!
!                                    A*Y = X,
!
!     by using the L and U factors computed with a LU factorization of A. In this
!     case K = L*U = A, alpha = 1 and beta = 0. The LU factorization is performed
!     by one of the following auxiliary pakages:
!     a. UMFPACK,
!     b. SuperLU,
!     c. SuperLU_Dist.
!     In the cases a. and b., the factors L and U (and the matrix A) are either
!     distributed and block diagonal) or replicated; in the case c., L, U (and A)
!     are distributed.
!
!  This routine is used by mld_dsub_aply, to apply a 'base' block-Jacobi or
!  Additive Schwarz (AS) preconditioner at any level of a multilevel preconditioner,
!  or a block-Jacobi or LU or ILU solver at the coarsest level of a multilevel
!  preconditioner. 
!
!
! Arguments:
!
!   alpha      -  complex(psb_dpk_), input.
!                 The scalar alpha.
!   prec       -  type(mld_zbaseprec_type), input.
!                 The 'base preconditioner' data structure containing the local 
!                 part of the L and U factors of the matrix A.
!   x          -  complex(psb_dpk_), dimension(:), input.
!                 The local part of the vector X.
!   beta       -  complex(psb_dpk_), input.
!                 The scalar beta.
!   y          -  complex(psb_dpk_), dimension(:), input/output.
!                 The local part of the vector Y.
!   desc_data  -  type(psb_desc_type), input.
!                 The communication descriptor associated to the matrix to be
!                 preconditioned or 'inverted'.
!   trans      -  character(len=1), input.
!                 If trans='N','n' then op(K^(-1)) = K^(-1);
!                 if trans='T','t' then op(K^(-1)) = K^(-T) (transpose of K^(-1)).
!                 if trans='C','c' then op(K^(-1)) = K^(-C) (transpose conjugate of K^(-1)).
!                 If prec%iprcparm(mld_smoother_sweeps_) > 1, the value of trans provided
!                 in input is ignored.
!   work       -  complex(psb_dpk_), dimension (:), target.
!                 Workspace. Its size must be at least 4*psb_cd_get_local_cols(desc_data).
!   info       -  integer, output.
!                 Error code.
!  
subroutine mld_zsub_solve(alpha,prec,x,beta,y,desc_data,trans,work,info)

  use psb_sparse_mod
  use mld_inner_mod, mld_protect_name => mld_zsub_solve

  implicit none 

  ! Arguments
  type(psb_desc_type), intent(in)           :: desc_data
  type(mld_zbaseprec_type), intent(in)       :: prec
  complex(psb_dpk_),intent(in)            :: x(:)
  complex(psb_dpk_),intent(inout)         :: y(:)
  complex(psb_dpk_),intent(in)            :: alpha,beta
  character(len=1), intent(in)              :: trans
  complex(psb_dpk_),target, intent(inout) :: work(:)
  integer, intent(out)                      :: info

  ! Local variables
  integer :: n_row,n_col
  complex(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer :: ictxt,np,me,i, err_act
  character(len=20)   :: name
  character           :: trans_

  interface 
    subroutine mld_zumf_solve(flag,m,x,b,n,ptr,info)
      use psb_sparse_mod
      integer, intent(in)  :: flag,m,n,ptr
      integer, intent(out) :: info
      complex(psb_dpk_), intent(in)    :: b(*)
      complex(psb_dpk_), intent(inout) :: x(*)
    end subroutine mld_zumf_solve
  end interface

  name='mld_zsub_solve'
  info = psb_success_
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_data)
  call psb_info(ictxt, me, np)

  trans_ = psb_toupper(trans)
  select case(trans_)
  case('N')
  case('T','C')
  case default
    call psb_errpush(psb_err_iarg_invalid_i_,name)
    goto 9999
  end select


  n_row = psb_cd_get_local_rows(desc_data)
  n_col = psb_cd_get_local_cols(desc_data)

  if (n_col <= size(work)) then 
    ww => work(1:n_col)
    if ((4*n_col+n_col) <= size(work)) then 
      aux => work(n_col+1:)
    else
      allocate(aux(4*n_col),stat=info)
      if (info /= psb_success_) then 
        info=psb_err_alloc_request_
        call psb_errpush(info,name,i_err=(/4*n_col,0,0,0,0/),&
             & a_err='complex(psb_dpk_)')
        goto 9999      
      end if
    endif
  else
    allocate(ww(n_col),aux(4*n_col),stat=info)
    if (info /= psb_success_) then 
      info=psb_err_alloc_request_
      call psb_errpush(info,name,i_err=(/5*n_col,0,0,0,0/),&
           & a_err='complex(psb_dpk_)')
      goto 9999      
    end if
  endif


  select case(prec%iprcparm(mld_sub_solve_))
  case(mld_ilu_n_,mld_milu_n_,mld_ilu_t_) 
    !
    ! Apply a block-Jacobi preconditioner with ILU(k)/MILU(k)/ILU(k,t)
    ! factorization of the blocks (distributed matrix) or approximately
    ! solve a system through ILU(k)/MILU(k)/ILU(k,t) (replicated matrix).
    ! 

    select case(trans_)
    case('N')

      call psb_spsm(zone,prec%av(mld_l_pr_),x,zzero,ww,desc_data,info,&
           & trans=trans_,unit='L',diag=prec%d,choice=psb_none_,work=aux)
      if (info == psb_success_) call psb_spsm(alpha,prec%av(mld_u_pr_),ww,beta,y,desc_data,info,&
           & trans=trans_,unit='U',choice=psb_none_, work=aux)

    case('T')
      call psb_spsm(zone,prec%av(mld_u_pr_),x,zzero,ww,desc_data,info,&
           & trans=trans_,unit='L',diag=prec%d,choice=psb_none_, work=aux)
      if(info == psb_success_) call psb_spsm(alpha,prec%av(mld_l_pr_),ww,beta,y,desc_data,info,&
           & trans=trans_,unit='U',choice=psb_none_,work=aux)

    case('C')
      call psb_spsm(zone,prec%av(mld_u_pr_),x,zzero,ww,desc_data,info,&
           & trans=trans_,unit='L',diag=conjg(prec%d),choice=psb_none_, work=aux)
      if(info == psb_success_) call psb_spsm(alpha,prec%av(mld_l_pr_),ww,beta,y,desc_data,info,&
           & trans=trans_,unit='U',choice=psb_none_,work=aux)
    case default
      call psb_errpush(psb_err_internal_error_,name,a_err='Invalid TRANS in ILU subsolve')
      goto 9999
    end select

  case(mld_slu_)
    !
    ! Apply a block-Jacobi preconditioner with LU factorization of the
    ! blocks (distributed matrix) or approximately solve a local linear
    ! system through LU (replicated matrix). The SuperLU package is used 
    ! to apply the LU factorization in both cases.
    !

    ww(1:n_row) = x(1:n_row)

    select case(trans_)
    case('N')
      call mld_zslu_solve(0,n_row,1,ww,n_row,prec%iprcparm(mld_slu_ptr_),info)
    case('T')
      call mld_zslu_solve(1,n_row,1,ww,n_row,prec%iprcparm(mld_slu_ptr_),info)
    case('C')
      call mld_zslu_solve(2,n_row,1,ww,n_row,prec%iprcparm(mld_slu_ptr_),info)
    case default 
      call psb_errpush(psb_err_internal_error_,name,a_err='Invalid TRANS in SLU subsolve')
      goto 9999
    end select

    if (info == psb_success_) call psb_geaxpby(alpha,ww,beta,y,desc_data,info)

  case(mld_sludist_)
    !
    ! Solve a distributed linear system with the LU factorization.
    ! The SuperLU_DIST package is used.
    !

    ww(1:n_row) = x(1:n_row)

    select case(trans_)
    case('N')
      call mld_zsludist_solve(0,n_row,1,ww,n_row,prec%iprcparm(mld_slud_ptr_),info)
    case('T')
      call mld_zsludist_solve(1,n_row,1,ww,n_row,prec%iprcparm(mld_slud_ptr_),info)
    case('C')
      call mld_zsludist_solve(2,n_row,1,ww,n_row,prec%iprcparm(mld_slud_ptr_),info)
    case default 
      call psb_errpush(psb_err_internal_error_,name,a_err='Invalid TRANS in SLUDist subsolve')
      goto 9999
    end select

    if (info == psb_success_) call psb_geaxpby(alpha,ww,beta,y,desc_data,info)

  case (mld_umf_) 
    !
    ! Apply a block-Jacobi preconditioner with LU factorization of the
    ! blocks (distributed matrix) or approximately solve a local linear
    ! system through LU (replicated matrix). The UMFPACK package is used 
    ! to apply the LU factorization in both cases.
    !

    select case(trans_)
    case('N')
      call mld_zumf_solve(0,n_row,ww,x,n_row,prec%iprcparm(mld_umf_numptr_),info)
    case('T')
      call mld_zumf_solve(1,n_row,ww,x,n_row,prec%iprcparm(mld_umf_numptr_),info)
    case('C')
      call mld_zumf_solve(2,n_row,ww,x,n_row,prec%iprcparm(mld_umf_numptr_),info)
    case default 
      call psb_errpush(psb_err_internal_error_,name,a_err='Invalid TRANS in UMF subsolve')
      goto 9999
    end select

    if (info == psb_success_) call psb_geaxpby(alpha,ww,beta,y,desc_data,info)

  case default
    call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_solve_')
    goto 9999

  end select

  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='Error in subsolve ')
    goto 9999
  endif

  if (n_col <= size(work)) then 
    if ((4*n_col+n_col) <= size(work)) then 
    else
      deallocate(aux)
    endif
  else
    deallocate(ww,aux)
  endif

  call psb_erractionrestore(err_act)
  return

9999 continue

  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_zsub_solve

