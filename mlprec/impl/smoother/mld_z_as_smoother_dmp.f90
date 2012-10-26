subroutine mld_z_as_smoother_dmp(sm,ictxt,level,info,prefix,head,smoother,solver)
  
  use psb_base_mod
  use mld_z_as_smoother, mld_protect_nam => mld_z_as_smoother_dmp
  implicit none 
  class(mld_z_as_smoother_type), intent(in) :: sm
  integer, intent(in)              :: ictxt,level
  integer, intent(out)             :: info
  character(len=*), intent(in), optional :: prefix, head
  logical, optional, intent(in)    :: smoother, solver
  integer :: i, j, il1, iln, lname, lev
  integer :: icontxt,iam, np
  character(len=80)  :: prefix_
  character(len=120) :: fname ! len should be at least 20 more than
  logical :: smoother_
  !  len of prefix_ 

  info = 0

  if (present(prefix)) then 
    prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
  else
    prefix_ = "dump_smth_z"
  end if

  call psb_info(ictxt,iam,np)

  if (present(smoother)) then 
    smoother_ = smoother
  else
    smoother_ = .false. 
  end if
  lname = len_trim(prefix_)
  fname = trim(prefix_)
  write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
  lname = lname + 5

  if (smoother_) then 
    write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_nd.mtx'
    if (sm%nd%is_asb()) &
         & call sm%nd%print(fname,head=head)
  end if
  ! At base level do nothing for the smoother
  if (allocated(sm%sv)) &
       & call sm%sv%dump(ictxt,level,info,solver=solver,prefix=prefix)

end subroutine mld_z_as_smoother_dmp
