subroutine mld_s_base_solver_check(sv,info)
  
  use psb_base_mod
  use mld_s_base_solver_mod, mld_protect_name =>  mld_s_base_solver_check
  Implicit None
  ! Arguments
  class(mld_s_base_solver_type), intent(inout) :: sv
  integer, intent(out)                   :: info
  Integer           :: err_act
  character(len=20) :: name='d_base_solver_check'

  call psb_erractionsave(err_act)
  info = psb_success_


  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine mld_s_base_solver_check
