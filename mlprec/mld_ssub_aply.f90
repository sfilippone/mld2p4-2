!!$
!!$ 
!!$                           MLD2P4  version 1.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 2.2)
!!$  
!!$  (C) Copyright 2008
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
! File mld_ssub_aply.f90
!
! Subroutine: mld_ssub_aply
! Version: real 
!
!  This routine computes
!
!                       Y = beta*Y + alpha*op(K^(-1))*X,
!
!  where
!  - K is a suitable matrix, as specified below,
!  - op(K^(-1)) is K^(-1) or its transpose, according to the value of the
!    argument trans,
!  - X and Y are vectors,
!  - alpha and beta are scalars.
!
!  Depending on K, alpha and beta (and on the communication descriptor desc_data
!  - see the arguments below), the above computation may correspond to one of
!  the following tasks:
!
!  1. Application of a block-Jacobi preconditioner associated to a matrix A
!     distributed among the processes. Here K is the preconditioner, op(K^(-1))
!     = K^(-1), alpha = 1 and beta = 0.
!
!  2. Application of block-Jacobi sweeps to compute an approximate solution of
!     a linear system
!                                    A*Y = X,
!
!     distributed among the processes (note that a single block-Jacobi sweep,
!     with null starting guess, corresponds to the application of a block-Jacobi
!     preconditioner). Here K^(-1) denotes the iteration matrix of the
!     block-Jacobi solver, op(K^(-1)) = K^(-1), alpha = 1 and beta = 0.
!
!  3. Solution, through the LU factorization, of a linear system
!
!                                    A*Y = X,
!
!     distributed among the processes. Here K = L*U = A, op(K^(-1)) = K^(-1),
!     alpha = 1 and beta = 0.
!
!  4. (Approximate) solution, through the LU or incomplete LU factorization, of
!     a linear system
!                                    A*Y = X,
!
!     replicated on the processes. Here K = L*U = A or K = L*U ~ A, op(K^(-1)) =
!     K^(-1), alpha = 1 and beta = 0.
!
!  The block-Jacobi preconditioner or solver and the L and U factors of the LU
!  or ILU factorizations have been built by the routine mld_fact_bld and stored
!  into the 'base preconditioner' data structure prec. See mld_fact_bld for more
!  details.
!
!  This routine is used by mld_as_aply, to apply a 'base' block-Jacobi or
!  Additive Schwarz (AS) preconditioner at any level of a multilevel preconditioner,
!  or a block-Jacobi or LU or ILU solver at the coarsest level of a multilevel
!  preconditioner. 
!
!  Tasks 1, 3 and 4 may be selected when prec%iprcparm(mld_smoother_sweeps_) = 1, 
!  while task 2 is selected when prec%iprcparm(mld_smoother_sweeps_) > 1. Furthermore
!  Tasks 1, 2 and 3 may be performed when the matrix A is
!  distributed among the processes (prec%iprcparm(mld_coarse_mat_) = mld_distr_mat_),
!  while task 4 may be performed when A is replicated on the processes
!  (prec%iprcparm(mld_coarse_mat_) = mld_repl_mat_). Note that the matrix A is
!  distributed among the processes at each level of the multilevel preconditioner,
!  except the coarsest one, where it may be either distributed or replicated on
!  the processes.  Tasks 2, 3 and 4 are performed only at the coarsest level.
!  Note also that this routine manages implicitly the fact that
!  the matrix is distributed or replicated, i.e. it does not make any explicit
!  reference to the value of prec%iprcparm(mld_coarse_mat_).
!
! Arguments:
!
!   alpha      -  real(psb_spk_), input.
!                 The scalar alpha.
!   prec       -  type(mld_sbaseprec_type), input.
!                 The 'base preconditioner' data structure containing the local 
!                 part of the preconditioner or solver.
!   x          -  real(psb_spk_), dimension(:), input.
!                 The local part of the vector X.
!   beta       -  real(psb_spk_), input.
!                 The scalar beta.
!   y          -  real(psb_spk_), dimension(:), input/output.
!                 The local part of the vector Y.
!   desc_data  -  type(psb_desc_type), input.
!                 The communication descriptor associated to the matrix to be
!                 preconditioned or 'inverted'.
!   trans      -  character(len=1), input.
!                 If trans='N','n' then op(K^(-1)) = K^(-1);
!                 if trans='T','t' then op(K^(-1)) = K^(-T) (transpose of K^(-1)).
!                 If prec%iprcparm(mld_smoother_sweeps_) > 1, the value of trans provided
!                 in input is ignored.
!   work       -  real(psb_spk_), dimension (:), target.
!                 Workspace. Its size must be at least 4*psb_cd_get_local_cols(desc_data).
!   info       -  integer, output.
!                 Error code.
!  
subroutine mld_ssub_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)

  use psb_base_mod
  use mld_inner_mod, mld_protect_name => mld_ssub_aply

  implicit none 

  ! Arguments
  type(psb_desc_type), intent(in)        :: desc_data
  type(mld_sbaseprec_type), intent(in)    :: prec
  real(psb_spk_),intent(in)            :: x(:)
  real(psb_spk_),intent(inout)         :: y(:)
  real(psb_spk_),intent(in)            :: alpha,beta
  character(len=1),intent(in)            :: trans
  real(psb_spk_),target, intent(inout) :: work(:)
  integer, intent(out)                   :: info

  ! Local variables
  integer :: n_row,n_col
  real(psb_spk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer :: ictxt,np,me,i, err_act
  character(len=20)   :: name
  character           :: trans_

  name='mld_ssub_aply'
  info = 0
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_data)
  call psb_info(ictxt, me, np)

  trans_ = psb_toupper(trans)
  select case(trans_)
  case('N')
  case('T','C')
  case default
    call psb_errpush(40,name)
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
      if (info /= 0) then 
        info=4025
        call psb_errpush(info,name,i_err=(/4*n_col,0,0,0,0/),&
             & a_err='real(psb_spk_)')
        goto 9999      
      end if
    endif
  else
    allocate(ww(n_col),aux(4*n_col),stat=info)
    if (info /= 0) then 
      info=4025
      call psb_errpush(info,name,i_err=(/5*n_col,0,0,0,0/),&
           & a_err='real(psb_spk_)')
      goto 9999      
    end if
  endif

  if (prec%iprcparm(mld_smoother_sweeps_) == 1) then 
    
    call mld_sub_solve(alpha,prec,x,beta,y,desc_data,trans_,aux,info) 
    
    if (info /= 0) then
      call psb_errpush(4001,name,a_err='Error in sub_aply Jacobi Sweeps = 1')
      goto 9999
    endif
    
  else if (prec%iprcparm(mld_smoother_sweeps_) > 1) then 
    !
    !
    ! Apply prec%iprcparm(mld_smoother_sweeps_) sweeps of a block-Jacobi solver
    ! to compute an approximate solution of a linear system.
    !
    !

    if (size(prec%av) < mld_ap_nd_) then 
      info = 4011
      goto 9999
    endif

    allocate(tx(n_col),ty(n_col),stat=info)
    if (info /= 0) then 
      info=4025
      call psb_errpush(info,name,i_err=(/2*n_col,0,0,0,0/),&
           & a_err='real(psb_spk_)')
      goto 9999      
    end if

    tx = szero
    ty = szero
    do i=1, prec%iprcparm(mld_smoother_sweeps_) 
      !
      ! Compute Y(j+1) = D^(-1)*(X-ND*Y(j)), where D and ND are the
      ! block diagonal part and the remaining part of the local matrix
      ! and Y(j) is the approximate solution at sweep j.
      !
      ty(1:n_row) = x(1:n_row)
      call psb_spmm(-sone,prec%av(mld_ap_nd_),tx,sone,ty,&
           &   prec%desc_data,info,work=aux,trans=trans_)

      if (info /=0) exit
      
      call mld_sub_solve(sone,prec,ty,szero,tx,desc_data,trans_,aux,info) 

      if (info /=0) exit
    end do

    if (info == 0) call psb_geaxpby(alpha,tx,beta,y,desc_data,info)

    if (info /= 0) then 
      info=4001
      call psb_errpush(info,name,a_err='subsolve with Jacobi sweeps > 1')
      goto 9999      
    end if

    deallocate(tx,ty,stat=info)
    if (info /= 0) then 
      info=4001
      call psb_errpush(info,name,a_err='final cleanup with Jacobi sweeps > 1')
      goto 9999      
    end if

  else

    info = 10
    call psb_errpush(info,name,&
         & i_err=(/2,prec%iprcparm(mld_smoother_sweeps_),0,0,0/))
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

end subroutine mld_ssub_aply

