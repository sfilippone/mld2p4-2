subroutine mld_z_as_smoother_setc(sm,what,val,info)
  
  use psb_base_mod
  use mld_z_as_smoother, mld_protect_nam => mld_z_as_smoother_setc
  Implicit None
  ! Arguments
  class(mld_z_as_smoother_type), intent(inout) :: sm
  integer, intent(in)                    :: what 
  character(len=*), intent(in)           :: val
  integer, intent(out)                   :: info
  Integer :: err_act, ival
  character(len=20)  :: name='z_as_smoother_setc'

  info = psb_success_
  call psb_erractionsave(err_act)


  call mld_stringval(val,ival,info)
  if (info == psb_success_) call sm%set(what,ival,info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info, name)
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
end subroutine mld_z_as_smoother_setc
