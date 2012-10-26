subroutine mld_z_base_smoother_descr(sm,info,iout,coarse)
  
  use psb_base_mod
  use mld_z_base_smoother_mod, mld_protect_name =>  mld_z_base_smoother_descr
  Implicit None

  ! Arguments
  class(mld_z_base_smoother_type), intent(in) :: sm
  integer, intent(out)                        :: info
  integer, intent(in), optional               :: iout
  logical, intent(in), optional               :: coarse

  ! Local variables
  integer      :: err_act
  integer      :: ictxt, me, np
  character(len=20), parameter :: name='mld_z_base_smoother_descr'
  integer :: iout_
  logical      :: coarse_


  call psb_erractionsave(err_act)
  info = psb_success_

  if (present(coarse)) then 
    coarse_ = coarse
  else
    coarse_ = .false.
  end if
  if (present(iout)) then 
    iout_ = iout
  else 
    iout_ = 6
  end if

  if (.not.coarse_) &
       &  write(iout_,*) 'Base smoother with local solver'
  if (allocated(sm%sv)) then 
    call sm%sv%descr(info,iout,coarse)
    if (info /= psb_success_) then 
      info = psb_err_from_subroutine_ 
      call psb_errpush(info,name,a_err='Local solver')
      goto 9999
    end if
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
end subroutine mld_z_base_smoother_descr
