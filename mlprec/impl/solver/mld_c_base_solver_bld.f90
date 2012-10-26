subroutine mld_c_base_solver_bld(a,desc_a,sv,upd,info,b,amold,vmold)
  
  use psb_base_mod
  use mld_c_base_solver_mod, mld_protect_name =>  mld_c_base_solver_bld
  Implicit None
  ! Arguments
  type(psb_cspmat_type), intent(in), target           :: a
  Type(psb_desc_type), Intent(in)                     :: desc_a 
  class(mld_c_base_solver_type), intent(inout)        :: sv
  character, intent(in)                               :: upd
  integer, intent(out)                                :: info
  type(psb_cspmat_type), intent(in), target, optional :: b
  class(psb_c_base_sparse_mat), intent(in), optional  :: amold
  class(psb_c_base_vect_type), intent(in), optional   :: vmold

  Integer :: err_act
  character(len=20)  :: name='d_base_solver_bld'

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
end subroutine mld_c_base_solver_bld
