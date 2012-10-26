subroutine mld_c_base_solver_setr(sv,what,val,info)
  
  use psb_base_mod
  use mld_c_base_solver_mod, mld_protect_name =>  mld_c_base_solver_setr
  Implicit None
  ! Arguments
  class(mld_c_base_solver_type), intent(inout) :: sv 
  integer, intent(in)                          :: what 
  real(psb_spk_), intent(in)                   :: val
  integer, intent(out)                         :: info
  Integer           :: err_act
  character(len=20) :: name='d_base_solver_setr'


  ! Correct action here is doing nothing. 
  info = 0

  return
end subroutine mld_c_base_solver_setr
