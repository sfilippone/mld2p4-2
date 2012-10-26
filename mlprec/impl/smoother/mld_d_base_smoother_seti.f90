subroutine mld_d_base_smoother_seti(sm,what,val,info)
  
  use psb_base_mod
  use mld_d_base_smoother_mod, mld_protect_name =>  mld_d_base_smoother_seti
  Implicit None
  ! Arguments
  class(mld_d_base_smoother_type), intent(inout) :: sm 
  integer, intent(in)                            :: what 
  integer, intent(in)                            :: val
  integer, intent(out)                           :: info
  Integer           :: err_act
  character(len=20) :: name='d_base_smoother_seti'

  call psb_erractionsave(err_act)
  info = psb_success_

  if (allocated(sm%sv)) then 
    call sm%sv%set(what,val,info)
  end if
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
end subroutine mld_d_base_smoother_seti
