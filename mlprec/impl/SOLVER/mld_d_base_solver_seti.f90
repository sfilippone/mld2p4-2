subroutine mld_d_base_solver_seti(sv,what,val,info)
  
  use psb_base_mod
  use mld_d_base_solver_mod, mld_protect_name =>  mld_d_base_solver_seti
  Implicit None
  ! Arguments
  class(mld_d_base_solver_type), intent(inout) :: sv 
  integer, intent(in)                          :: what 
  integer, intent(in)                          :: val
  integer, intent(out)                         :: info
  Integer           :: err_act
  character(len=20) :: name='d_base_solver_seti'

  ! Correct action here is doing nothing. 
  info = 0

  return
end subroutine mld_d_base_solver_seti
