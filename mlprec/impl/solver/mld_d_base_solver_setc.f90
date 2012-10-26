subroutine mld_d_base_solver_setc(sv,what,val,info)
  
  use psb_base_mod
  use mld_d_base_solver_mod, mld_protect_name =>  mld_d_base_solver_setc
  Implicit None
  ! Arguments
  class(mld_d_base_solver_type), intent(inout) :: sv
  integer, intent(in)                          :: what 
  character(len=*), intent(in)                 :: val
  integer, intent(out)                         :: info
  Integer           :: err_act, ival 
  character(len=20) :: name='d_base_solver_setc'

  call psb_erractionsave(err_act)

  info = psb_success_

  call mld_stringval(val,ival,info)
  if (info == psb_success_) call sv%set(what,ival,info)

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
end subroutine mld_d_base_solver_setc
