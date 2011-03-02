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
! File: mld_zbaseprec_aply.f90
!
! Subroutine: mld_zbaseprec_aply
! Version:    complex
!
!  This routine applies a base preconditioner by computing
!
!                          Y = beta*Y + alpha*op(K^(-1))*X,
!  where
!  - K is the base preconditioner, stored in prec,
!  - op(K^(-1)) is K^(-1) or its transpose, according to the value of trans,
!  - X and Y are vectors,
!  - alpha and beta are scalars.
!
!  The routine is used by mld_dmlprec_aply, to apply the multilevel preconditioners,
!  or directly by mld_dprec_aply, to apply the basic one-level preconditioners (diagonal,
!  block-Jacobi or additive Schwarz). It also manages the case of no preconditioning.
!
!
! Arguments:
!   alpha      -  complex(psb_dpk_), input.
!                 The scalar alpha.
!   prec       -  type(mld_zbaseprec_type), input.
!                 The base preconditioner data structure containing the local part
!                 of the preconditioner K.
!   x          -  complex(psb_dpk_), dimension(:), input.
!                 The local part of the vector X.
!   beta       -  complex(psb_dpk_), input.
!                 The scalar beta.
!   y          -  complex(psb_dpk_), dimension(:), input/output.
!                 The local part of the vector Y.
!   desc_data  -  type(psb_desc_type), input.
!                 The communication descriptor associated to the matrix to be
!                 preconditioned.
!   trans      -  character, optional.
!                 If trans='N','n' then op(K^(-1)) = K^(-1);
!                 if trans='T','t' then op(K^(-1)) = K^(-T) (transpose of K^(-1)).
!   work       -  real(psb_dpk_), dimension (:), optional, target.
!                 Workspace. Its size must be at least 4*psb_cd_get_local_cols(desc_data).
!   info       -  integer, output.
!                 Error code.
!  
subroutine mld_zbaseprec_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)

  use psb_sparse_mod
  use mld_z_inner_mod, mld_protect_name => mld_zbaseprec_aply

  implicit none 

  ! Arguments
  type(psb_desc_type),intent(in)      :: desc_data
  type(mld_zbaseprec_type), intent(in) :: prec
  complex(psb_dpk_),intent(in)      :: x(:)
  complex(psb_dpk_),intent(inout)   :: y(:)
  complex(psb_dpk_),intent(in)      :: alpha,beta
  character(len=1)                    :: trans
  complex(psb_dpk_),target          :: work(:)
  integer, intent(out)                :: info

  ! Local variables
  complex(psb_dpk_), pointer :: ww(:)
  integer           :: ictxt, np, me, err_act
  integer           :: n_row, int_err(5)
  character(len=20) :: name, ch_err
  character         :: trans_

  name='mld_zbaseprec_aply'
  info = psb_success_
  call psb_erractionsave(err_act)

  ictxt = psb_cd_get_context(desc_data)

  call psb_info(ictxt, me, np)

  trans_= psb_toupper(trans)
  select case(trans_)
  case('N','T','C')
    ! Ok
  case default
    info=psb_err_iarg_invalid_i_
    int_err(1)=6
    ch_err(2:2)=trans
    goto 9999
  end select

  select case(prec%iprcparm(mld_smoother_type_))

  case(mld_noprec_)
    !
    ! No preconditioner
    !

    call psb_geaxpby(alpha,x,beta,y,desc_data,info)

  case(mld_diag_)
    !
    ! Diagonal preconditioner
    !

    if (size(work) >= size(x)) then 
      ww => work
    else
      allocate(ww(size(x)),stat=info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_alloc_request_,name,i_err=(/size(x),0,0,0,0/),a_err='complex(psb_dpk_)')
        goto 9999      
      end if
    end if

    n_row = psb_cd_get_local_rows(desc_data)
    if (trans_ == 'C') then 
      ww(1:n_row) = x(1:n_row)*conjg(prec%d(1:n_row))
    else
      ww(1:n_row) = x(1:n_row)*prec%d(1:n_row)
    end if
    call psb_geaxpby(alpha,ww,beta,y,desc_data,info)

    if (size(work) < size(x)) then 
      deallocate(ww,stat=info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_from_subroutine_,name,a_err='Deallocate')
        goto 9999      
      end if
    end if

  case(mld_bjac_,mld_as_)
    !
    ! Additive Schwarz preconditioner (including block-Jacobi as special case)
    !
    call mld_as_aply(alpha,prec,x,beta,y,desc_data,trans_,work,info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='mld_as_aply'
      goto 9999
    end if

  case default
    call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_smoother_type_')
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name,i_err=int_err,a_err=ch_err)
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_zbaseprec_aply

