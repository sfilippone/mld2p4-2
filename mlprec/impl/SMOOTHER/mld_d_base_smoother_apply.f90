subroutine mld_d_base_smoother_apply(alpha,sm,x,beta,y,desc_data,trans,sweeps,work,info)
  
  use psb_base_mod
  use mld_d_base_smoother_mod, mld_protect_name =>  mld_d_base_smoother_apply
  implicit none 
  type(psb_desc_type), intent(in)             :: desc_data
  class(mld_d_base_smoother_type), intent(in) :: sm
  real(psb_dpk_),intent(inout)                :: x(:)
  real(psb_dpk_),intent(inout)                :: y(:)
  real(psb_dpk_),intent(in)                   :: alpha,beta
  character(len=1),intent(in)                 :: trans
  integer, intent(in)                         :: sweeps
  real(psb_dpk_),target, intent(inout)        :: work(:)
  integer, intent(out)                        :: info

  Integer           :: err_act
  character(len=20) :: name='d_base_smoother_apply'

  call psb_erractionsave(err_act)
  info = psb_success_
  if (allocated(sm%sv)) then 
    call sm%sv%apply(alpha,x,beta,y,desc_data,trans,work,info)
  else
    info = 1121
  endif
  if (info /= psb_success_) then 
    call psb_errpush(info,name)
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

end subroutine mld_d_base_smoother_apply
