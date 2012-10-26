subroutine mld_d_base_smoother_bld(a,desc_a,sm,upd,info,amold,vmold)
  
  use psb_base_mod
  use mld_d_base_smoother_mod, mld_protect_name =>  mld_d_base_smoother_bld
  Implicit None

  ! Arguments
  type(psb_dspmat_type), intent(in), target      :: a
  Type(psb_desc_type), Intent(in)                :: desc_a 
  class(mld_d_base_smoother_type), intent(inout) :: sm 
  character, intent(in)                          :: upd
  integer, intent(out)                           :: info
  class(psb_d_base_sparse_mat), intent(in), optional :: amold
  class(psb_d_base_vect_type), intent(in), optional  :: vmold
  Integer           :: err_act
  character(len=20) :: name='d_base_smoother_bld'

  call psb_erractionsave(err_act)

  info = psb_success_
  if (allocated(sm%sv)) then 
    call sm%sv%build(a,desc_a,upd,info,amold=amold,vmold=vmold)
  else
    info = 1121
    call psb_errpush(info,name)
  endif
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
end subroutine mld_d_base_smoother_bld
