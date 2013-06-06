!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010,2012,2013
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
!
!
!
!
!
!
subroutine mld_d_diag_solver_apply_vect(alpha,sv,x,beta,y,desc_data,trans,work,info)
  use psb_base_mod
  use mld_d_diag_solver, mld_protect_name => mld_d_diag_solver_apply_vect
  implicit none 
  type(psb_desc_type), intent(in)              :: desc_data
  class(mld_d_diag_solver_type), intent(inout) :: sv
  type(psb_d_vect_type), intent(inout)         :: x
  type(psb_d_vect_type), intent(inout)         :: y
  real(psb_dpk_),intent(in)                    :: alpha,beta
  character(len=1),intent(in)                  :: trans
  real(psb_dpk_),target, intent(inout)         :: work(:)
  integer, intent(out)                         :: info

  integer    :: n_row,n_col
  real(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer    :: ictxt,np,me,i, err_act
  character          :: trans_
  character(len=20)  :: name='d_diag_solver_apply'

  call psb_erractionsave(err_act)

  info = psb_success_

  trans_ = psb_toupper(trans)
  select case(trans_)
  case('N')
  case('T','C')
  case default
    call psb_errpush(psb_err_iarg_invalid_i_,name)
    goto 9999
  end select

  n_row = desc_data%get_local_rows()
  n_col = desc_data%get_local_cols()
  if (x%get_nrows() < n_row) then 
    info = 36
    call psb_errpush(info,name,i_err=(/2,n_row,0,0,0/))
    goto 9999
  end if
  if (y%get_nrows() < n_row) then 
    info = 36
    call psb_errpush(info,name,i_err=(/3,n_row,0,0,0/))
    goto 9999
  end if
  if (.not.allocated(sv%dv)) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner: D")
    goto 9999
  end if
  if (sv%dv%get_nrows() < n_row) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner: D")
    goto 9999
  end if

  call y%mlt(alpha,sv%dv,x,beta,info,conjgx=trans_)

  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='vect%mlt')
    goto 9999      
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_d_diag_solver_apply_vect

subroutine mld_d_diag_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info)
  use psb_base_mod
  use mld_d_diag_solver, mld_protect_name => mld_d_diag_solver_apply
  implicit none 
  type(psb_desc_type), intent(in)           :: desc_data
  class(mld_d_diag_solver_type), intent(in) :: sv
  real(psb_dpk_), intent(inout)             :: x(:)
  real(psb_dpk_), intent(inout)             :: y(:)
  real(psb_dpk_),intent(in)                 :: alpha,beta
  character(len=1),intent(in)               :: trans
  real(psb_dpk_),target, intent(inout)      :: work(:)
  integer, intent(out)                      :: info

  integer    :: n_row,n_col
  real(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer    :: ictxt,np,me,i, err_act
  character          :: trans_
  character(len=20)  :: name='d_diag_solver_apply'

  call psb_erractionsave(err_act)

  info = psb_success_

  trans_ = psb_toupper(trans)
  select case(trans_)
  case('N')
  case('T','C')
  case default
    call psb_errpush(psb_err_iarg_invalid_i_,name)
    goto 9999
  end select

  n_row = desc_data%get_local_rows()
  n_col = desc_data%get_local_cols()

  if (trans_ == 'C') then 
    if (beta == dzero) then 

      if (alpha == dzero) then 
        y(1:n_row) = dzero
      else if (alpha == done) then 
        do i=1, n_row
          y(i) = (sv%d(i)) * x(i)
        end do
      else if (alpha == -done) then 
        do i=1, n_row
          y(i) = -(sv%d(i)) * x(i)
        end do
      else
        do i=1, n_row
          y(i) = alpha * (sv%d(i)) * x(i)
        end do
      end if

    else if (beta == done) then 

      if (alpha == dzero) then 
        !y(1:n_row) = dzero
      else if (alpha == done) then 
        do i=1, n_row
          y(i) = (sv%d(i)) * x(i) + y(i)
        end do
      else if (alpha == -done) then 
        do i=1, n_row
          y(i) = -(sv%d(i)) * x(i)  + y(i)
        end do
      else
        do i=1, n_row
          y(i) = alpha * (sv%d(i)) * x(i) + y(i)
        end do
      end if

    else if (beta == -done) then 

      if (alpha == dzero) then 
        y(1:n_row) = -y(1:n_row)        
      else if (alpha == done) then 
        do i=1, n_row
          y(i) = (sv%d(i)) * x(i) - y(i)
        end do
      else if (alpha == -done) then 
        do i=1, n_row
          y(i) = -(sv%d(i)) * x(i)  - y(i)
        end do
      else
        do i=1, n_row
          y(i) = alpha * (sv%d(i)) * x(i) - y(i)
        end do
      end if

    else

      if (alpha == dzero) then 
        y(1:n_row) = beta *y(1:n_row)        
      else if (alpha == done) then 
        do i=1, n_row
          y(i) = (sv%d(i)) * x(i) + beta*y(i)
        end do
      else if (alpha == -done) then 
        do i=1, n_row
          y(i) = -(sv%d(i)) * x(i)  + beta*y(i)
        end do
      else
        do i=1, n_row
          y(i) = alpha * (sv%d(i)) * x(i) + beta*y(i)
        end do
      end if

    end if

  else if (trans_ /= 'C') then 

    if (beta == dzero) then 

      if (alpha == dzero) then 
        y(1:n_row) = dzero
      else if (alpha == done) then 
        do i=1, n_row
          y(i) = sv%d(i) * x(i)
        end do
      else if (alpha == -done) then 
        do i=1, n_row
          y(i) = -sv%d(i) * x(i)
        end do
      else
        do i=1, n_row
          y(i) = alpha * sv%d(i) * x(i)
        end do
      end if

    else if (beta == done) then 

      if (alpha == dzero) then 
        !y(1:n_row) = dzero
      else if (alpha == done) then 
        do i=1, n_row
          y(i) = sv%d(i) * x(i) + y(i)
        end do
      else if (alpha == -done) then 
        do i=1, n_row
          y(i) = -sv%d(i) * x(i)  + y(i)
        end do
      else
        do i=1, n_row
          y(i) = alpha * sv%d(i) * x(i) + y(i)
        end do
      end if

    else if (beta == -done) then 

      if (alpha == dzero) then 
        y(1:n_row) = -y(1:n_row)        
      else if (alpha == done) then 
        do i=1, n_row
          y(i) = sv%d(i) * x(i) - y(i)
        end do
      else if (alpha == -done) then 
        do i=1, n_row
          y(i) = -sv%d(i) * x(i)  - y(i)
        end do
      else
        do i=1, n_row
          y(i) = alpha * sv%d(i) * x(i) - y(i)
        end do
      end if

    else

      if (alpha == dzero) then 
        y(1:n_row) = beta *y(1:n_row)        
      else if (alpha == done) then 
        do i=1, n_row
          y(i) = sv%d(i) * x(i) + beta*y(i)
        end do
      else if (alpha == -done) then 
        do i=1, n_row
          y(i) = -sv%d(i) * x(i)  + beta*y(i)
        end do
      else
        do i=1, n_row
          y(i) = alpha * sv%d(i) * x(i) + beta*y(i)
        end do
      end if

    end if

  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_d_diag_solver_apply

subroutine mld_d_diag_solver_bld(a,desc_a,sv,upd,info,b,amold,vmold)
  use psb_base_mod
  use mld_d_diag_solver, mld_protect_name => mld_d_diag_solver_bld

  Implicit None

  ! Arguments
  type(psb_dspmat_type), intent(in), target           :: a
  Type(psb_desc_type), Intent(in)                     :: desc_a 
  class(mld_d_diag_solver_type), intent(inout)        :: sv
  character, intent(in)                               :: upd
  integer, intent(out)                                :: info
  type(psb_dspmat_type), intent(in), target, optional :: b
  class(psb_d_base_sparse_mat), intent(in), optional  :: amold
  class(psb_d_base_vect_type), intent(in), optional   :: vmold
  ! Local variables
  integer :: n_row,n_col, nrow_a, nztota
  real(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer :: ictxt,np,me,i, err_act, debug_unit, debug_level
  character(len=20)  :: name='d_diag_solver_bld', ch_err

  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ictxt       = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'


  n_row  = desc_a%get_local_rows()
  nrow_a = a%get_nrows()
  if (allocated(sv%d)) then 
    if (size(sv%d) < n_row) then 
      deallocate(sv%d)
    endif
  endif
  if (.not.allocated(sv%d)) then 
    allocate(sv%d(n_row),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999      
    end if

  endif

  call a%get_diag(sv%d,info)
  if (present(b)) then 
    if (info == psb_success_) call b%get_diag(sv%d(nrow_a+1:), info)
  end if
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='get_diag')
    goto 9999      
  end if

  do i=1,n_row
    if (sv%d(i) == dzero) then 
      sv%d(i) = done
    else
      sv%d(i) = done/sv%d(i)
    end if
  end do
  allocate(sv%dv,stat=info) 
  if (info == psb_success_) then 
    if (present(vmold)) then 
      allocate(sv%dv%v,mold=vmold,stat=info) 
    else
      allocate(psb_d_base_vect_type :: sv%dv%v,stat=info) 
    end if
  end if
  if (info == psb_success_) then 
    call sv%dv%bld(sv%d)
  else
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate sv%dv')
    goto 9999      
  end if

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' end'

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine mld_d_diag_solver_bld

