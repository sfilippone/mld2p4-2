subroutine mld_z_base_solver_dmp(sv,ictxt,level,info,prefix,head,solver)
  
  use psb_base_mod
  use mld_z_base_solver_mod, mld_protect_name =>  mld_z_base_solver_dmp
  implicit none 
  class(mld_z_base_solver_type), intent(in) :: sv
  integer, intent(in)              :: ictxt,level
  integer, intent(out)             :: info
  character(len=*), intent(in), optional :: prefix, head
  logical, optional, intent(in)    :: solver
  integer :: i, j, il1, iln, lname, lev
  integer :: icontxt,iam, np
  character(len=80)  :: prefix_
  character(len=120) :: fname ! len should be at least 20 more than
  logical :: solver_
  !  len of prefix_ 

  info = 0

  if (present(prefix)) then 
    prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
  else
    prefix_ = "dump_slv_z"
  end if

  call psb_info(ictxt,iam,np)

  if (present(solver)) then 
    solver_ = solver
  else
    solver_ = .false. 
  end if
  lname = len_trim(prefix_)
  fname = trim(prefix_)
  write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
  lname = lname + 5

  ! At base level do nothing for the solver

end subroutine mld_z_base_solver_dmp
