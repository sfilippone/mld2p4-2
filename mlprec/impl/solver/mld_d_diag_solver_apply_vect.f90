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
