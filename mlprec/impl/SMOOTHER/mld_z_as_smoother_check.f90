subroutine mld_z_as_smoother_check(sm,info)
  

  use psb_base_mod
  use mld_z_as_smoother, mld_protect_nam => mld_z_as_smoother_check

  Implicit None

  ! Arguments
  class(mld_z_as_smoother_type), intent(inout) :: sm 
  integer, intent(out)                   :: info
  Integer           :: err_act
  character(len=20) :: name='z_as_smoother_check'

  call psb_erractionsave(err_act)
  info = psb_success_

  call mld_check_def(sm%restr,&
       & 'Restrictor',psb_halo_,is_legal_restrict)
  call mld_check_def(sm%prol,&
       & 'Prolongator',psb_none_,is_legal_prolong)
  call mld_check_def(sm%novr,&
       & 'Overlap layers ',0,is_legal_n_ovr)


  if (allocated(sm%sv)) then 
    call sm%sv%check(info)
  else 
    info=3111
    call psb_errpush(info,name)
    goto 9999
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
end subroutine mld_z_as_smoother_check
