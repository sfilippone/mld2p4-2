subroutine mld_s_base_smoother_setr(sm,what,val,info)
  
  use psb_base_mod
  use mld_s_base_smoother_mod, mld_protect_name =>  mld_s_base_smoother_setr
  Implicit None

  ! Arguments
  class(mld_s_base_smoother_type), intent(inout) :: sm 
  integer, intent(in)                            :: what 
  real(psb_spk_), intent(in)                     :: val
  integer, intent(out)                           :: info
  Integer :: err_act
  character(len=20)  :: name='s_base_smoother_setr'

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
end subroutine mld_s_base_smoother_setr
