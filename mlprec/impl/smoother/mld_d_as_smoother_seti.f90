subroutine mld_d_as_smoother_seti(sm,what,val,info)
  
  use psb_base_mod
  use mld_d_as_smoother, mld_protect_nam => mld_d_as_smoother_seti
  Implicit None

  ! Arguments
  class(mld_d_as_smoother_type), intent(inout) :: sm 
  integer, intent(in)                    :: what 
  integer, intent(in)                    :: val
  integer, intent(out)                   :: info
  Integer :: err_act
  character(len=20)  :: name='d_as_smoother_seti'

  info = psb_success_
  call psb_erractionsave(err_act)

  select case(what) 
!!$    case(mld_smoother_sweeps_) 
!!$      sm%sweeps = val
  case(mld_sub_ovr_) 
    sm%novr   = val
  case(mld_sub_restr_) 
    sm%restr  = val
  case(mld_sub_prol_) 
    sm%prol   = val
  case default
    if (allocated(sm%sv)) then 
      call sm%sv%set(what,val,info)
    end if
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine mld_d_as_smoother_seti
