subroutine mld_z_ilu_solver_dmp(sv,ictxt,level,info,prefix,head,solver)
  
  use psb_base_mod
  use mld_z_ilu_solver, mld_protect_name => mld_z_ilu_solver_dmp
  implicit none 
  class(mld_z_ilu_solver_type), intent(in) :: sv
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


  call psb_info(ictxt,iam,np)

  if (present(solver)) then 
    solver_ = solver
  else
    solver_ = .false. 
  end if

  if (solver_) then 
    if (present(prefix)) then 
      prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
    else
      prefix_ = "dump_slv_z"
    end if
    lname = len_trim(prefix_)
    fname = trim(prefix_)
    write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
    lname = lname + 5

    write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_lower.mtx'
    if (sv%l%is_asb()) &
         & call sv%l%print(fname,head=head)
    write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_diag.mtx'
    if (allocated(sv%d)) &
         & call psb_geprt(fname,sv%d,head=head)
    write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_upper.mtx'
    if (sv%u%is_asb()) &
         & call sv%u%print(fname,head=head)

  end if

end subroutine mld_z_ilu_solver_dmp
