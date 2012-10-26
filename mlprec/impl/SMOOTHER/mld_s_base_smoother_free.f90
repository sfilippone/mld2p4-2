subroutine mld_s_base_smoother_free(sm,info)
  
  use psb_base_mod
  use mld_s_base_smoother_mod, mld_protect_name => mld_s_base_smoother_free
  Implicit None

  ! Arguments
  class(mld_s_base_smoother_type), intent(inout) :: sm
  integer, intent(out)                           :: info
  Integer           :: err_act
  character(len=20) :: name='s_base_smoother_free'

  call psb_erractionsave(err_act)
  info = psb_success_

  if (allocated(sm%sv)) then 
    call sm%sv%free(info)
    if (info == psb_success_) deallocate(sm%sv,stat=info) 
  end if
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
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
end subroutine mld_s_base_smoother_free
