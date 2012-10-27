!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010,2012
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
subroutine mld_c_jac_smoother_apply_vect(alpha,sm,x,beta,y,desc_data,trans,sweeps,work,info)
  
  use psb_base_mod
  use mld_c_jac_smoother, mld_protect_name => mld_c_jac_smoother_apply_vect
  implicit none 
  type(psb_desc_type), intent(in)               :: desc_data
  class(mld_c_jac_smoother_type), intent(inout) :: sm
  type(psb_c_vect_type),intent(inout)           :: x
  type(psb_c_vect_type),intent(inout)           :: y
  complex(psb_spk_),intent(in)                     :: alpha,beta
  character(len=1),intent(in)                   :: trans
  integer, intent(in)                           :: sweeps
  complex(psb_spk_),target, intent(inout)          :: work(:)
  integer, intent(out)                          :: info

  integer    :: n_row,n_col
  type(psb_c_vect_type)  :: tx, ty
  complex(psb_spk_), pointer :: ww(:), aux(:)
  integer    :: ictxt,np,me,i, err_act
  character          :: trans_
  character(len=20)  :: name='c_jac_smoother_apply'

  call psb_erractionsave(err_act)

  info = psb_success_

  trans_ = psb_toupper(trans)
  select case(trans_)
  case('N')
  case('T')
  case('C')
  case default
    call psb_errpush(psb_err_iarg_invalid_i_,name)
    goto 9999
  end select

  if (.not.allocated(sm%sv)) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  end if

  n_row = desc_data%get_local_rows()
  n_col = desc_data%get_local_cols()

  if (n_col <= size(work)) then 
    ww => work(1:n_col)
    if ((4*n_col+n_col) <= size(work)) then 
      aux => work(n_col+1:)
    else
      allocate(aux(4*n_col),stat=info)
      if (info /= psb_success_) then 
        info=psb_err_alloc_request_
        call psb_errpush(info,name,i_err=(/4*n_col,0,0,0,0/),&
             & a_err='complex(psb_spk_)')
        goto 9999      
      end if
    endif
  else
    allocate(ww(n_col),aux(4*n_col),stat=info)
    if (info /= psb_success_) then 
      info=psb_err_alloc_request_
      call psb_errpush(info,name,i_err=(/5*n_col,0,0,0,0/),&
           & a_err='complex(psb_spk_)')
      goto 9999      
    end if
  endif

  if ((sweeps == 1).or.(sm%nnz_nd_tot==0)) then 

    call sm%sv%apply(alpha,x,beta,y,desc_data,trans_,aux,info) 

    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,&
           & name,a_err='Error in sub_aply Jacobi Sweeps = 1')
      goto 9999
    endif

  else if (sweeps  > 1) then 

    !
    !
    ! Apply multiple sweeps of a block-Jacobi solver
    ! to compute an approximate solution of a linear system.
    !
    !
    call tx%bld(x%get_nrows(),mold=x%v)
    call tx%set(czero)
    call ty%bld(x%get_nrows(),mold=x%v)

    do i=1, sweeps
      !
      ! Compute Y(j+1) = D^(-1)*(X-ND*Y(j)), where D and ND are the
      ! block diagonal part and the remaining part of the local matrix
      ! and Y(j) is the approximate solution at sweep j.
      !
      call psb_geaxpby(cone,x,czero,ty,desc_data,info)
      call psb_spmm(-cone,sm%nd,tx,cone,ty,desc_data,info,work=aux,trans=trans_)

      if (info /= psb_success_) exit

      call sm%sv%apply(cone,ty,czero,tx,desc_data,trans_,aux,info) 

      if (info /= psb_success_) exit
    end do

    if (info == psb_success_) call psb_geaxpby(alpha,tx,beta,y,desc_data,info)

    if (info /= psb_success_) then 
      info=psb_err_internal_error_
      call psb_errpush(info,name,a_err='subsolve with Jacobi sweeps > 1')
      goto 9999      
    end if

    call tx%free(info) 
    if (info == psb_success_) call ty%free(info) 
    if (info /= psb_success_) then 
      info=psb_err_internal_error_
      call psb_errpush(info,name,a_err='final cleanup with Jacobi sweeps > 1')
      goto 9999      
    end if

  else

    info = psb_err_iarg_neg_
    call psb_errpush(info,name,&
         & i_err=(/2,sweeps,0,0,0/))
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
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_c_jac_smoother_apply_vect
