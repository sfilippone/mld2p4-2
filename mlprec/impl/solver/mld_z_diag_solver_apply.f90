subroutine mld_z_diag_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info)
  
  use psb_base_mod
  use mld_z_diag_solver, mld_protect_name => mld_z_diag_solver_apply
  implicit none 
  type(psb_desc_type), intent(in)           :: desc_data
  class(mld_z_diag_solver_type), intent(in) :: sv
  complex(psb_dpk_), intent(inout)             :: x(:)
  complex(psb_dpk_), intent(inout)             :: y(:)
  complex(psb_dpk_),intent(in)                 :: alpha,beta
  character(len=1),intent(in)               :: trans
  complex(psb_dpk_),target, intent(inout)      :: work(:)
  integer, intent(out)                      :: info

  integer    :: n_row,n_col
  complex(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer    :: ictxt,np,me,i, err_act
  character          :: trans_
  character(len=20)  :: name='z_diag_solver_apply'

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
    if (beta == zzero) then 

      if (alpha == zzero) then 
        y(1:n_row) = zzero
      else if (alpha == zone) then 
        do i=1, n_row
          y(i) = conjg(sv%d(i)) * x(i)
        end do
      else if (alpha == -zone) then 
        do i=1, n_row
          y(i) = -conjg(sv%d(i)) * x(i)
        end do
      else
        do i=1, n_row
          y(i) = alpha * conjg(sv%d(i)) * x(i)
        end do
      end if

    else if (beta == zone) then 

      if (alpha == zzero) then 
        !y(1:n_row) = zzero
      else if (alpha == zone) then 
        do i=1, n_row
          y(i) = conjg(sv%d(i)) * x(i) + y(i)
        end do
      else if (alpha == -zone) then 
        do i=1, n_row
          y(i) = -conjg(sv%d(i)) * x(i)  + y(i)
        end do
      else
        do i=1, n_row
          y(i) = alpha * conjg(sv%d(i)) * x(i) + y(i)
        end do
      end if

    else if (beta == -zone) then 

      if (alpha == zzero) then 
        y(1:n_row) = -y(1:n_row)        
      else if (alpha == zone) then 
        do i=1, n_row
          y(i) = conjg(sv%d(i)) * x(i) - y(i)
        end do
      else if (alpha == -zone) then 
        do i=1, n_row
          y(i) = -conjg(sv%d(i)) * x(i)  - y(i)
        end do
      else
        do i=1, n_row
          y(i) = alpha * conjg(sv%d(i)) * x(i) - y(i)
        end do
      end if

    else

      if (alpha == zzero) then 
        y(1:n_row) = beta *y(1:n_row)        
      else if (alpha == zone) then 
        do i=1, n_row
          y(i) = conjg(sv%d(i)) * x(i) + beta*y(i)
        end do
      else if (alpha == -zone) then 
        do i=1, n_row
          y(i) = -conjg(sv%d(i)) * x(i)  + beta*y(i)
        end do
      else
        do i=1, n_row
          y(i) = alpha * conjg(sv%d(i)) * x(i) + beta*y(i)
        end do
      end if

    end if

  else if (trans_ /= 'C') then 

    if (beta == zzero) then 

      if (alpha == zzero) then 
        y(1:n_row) = zzero
      else if (alpha == zone) then 
        do i=1, n_row
          y(i) = sv%d(i) * x(i)
        end do
      else if (alpha == -zone) then 
        do i=1, n_row
          y(i) = -sv%d(i) * x(i)
        end do
      else
        do i=1, n_row
          y(i) = alpha * sv%d(i) * x(i)
        end do
      end if

    else if (beta == zone) then 

      if (alpha == zzero) then 
        !y(1:n_row) = zzero
      else if (alpha == zone) then 
        do i=1, n_row
          y(i) = sv%d(i) * x(i) + y(i)
        end do
      else if (alpha == -zone) then 
        do i=1, n_row
          y(i) = -sv%d(i) * x(i)  + y(i)
        end do
      else
        do i=1, n_row
          y(i) = alpha * sv%d(i) * x(i) + y(i)
        end do
      end if

    else if (beta == -zone) then 

      if (alpha == zzero) then 
        y(1:n_row) = -y(1:n_row)        
      else if (alpha == zone) then 
        do i=1, n_row
          y(i) = sv%d(i) * x(i) - y(i)
        end do
      else if (alpha == -zone) then 
        do i=1, n_row
          y(i) = -sv%d(i) * x(i)  - y(i)
        end do
      else
        do i=1, n_row
          y(i) = alpha * sv%d(i) * x(i) - y(i)
        end do
      end if

    else

      if (alpha == zzero) then 
        y(1:n_row) = beta *y(1:n_row)        
      else if (alpha == zone) then 
        do i=1, n_row
          y(i) = sv%d(i) * x(i) + beta*y(i)
        end do
      else if (alpha == -zone) then 
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

end subroutine mld_z_diag_solver_apply
