subroutine mld_z_base_solver_descr(sv,info,iout,coarse)
  
  use psb_base_mod
  use mld_z_base_solver_mod, mld_protect_name =>  mld_z_base_solver_descr
  Implicit None
  ! Arguments
  class(mld_z_base_solver_type), intent(in) :: sv
  integer, intent(out)                      :: info
  integer, intent(in), optional             :: iout
  logical, intent(in), optional             :: coarse

  ! Local variables
  integer      :: err_act
  integer      :: ictxt, me, np
  character(len=20), parameter :: name='mld_z_base_solver_descr'
  integer      :: iout_


  call psb_erractionsave(err_act)

  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine mld_z_base_solver_descr
