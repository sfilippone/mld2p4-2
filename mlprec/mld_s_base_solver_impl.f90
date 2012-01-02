!
! Apply: comes in two versions, on plain arrays or on encapsulated
! vectors.
! The base version throws an error, since it means that no explicit
! choice was made. 
! Question: would it make sense to transform the base version into
! the ID version, i.e. "base_solver" is the identity operator? 
! 

subroutine mld_s_base_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info)
  use psb_base_mod
  use mld_s_base_solver_mod, mld_protect_name =>  mld_s_base_solver_apply
  implicit none 
  type(psb_desc_type), intent(in)           :: desc_data
  class(mld_s_base_solver_type), intent(in) :: sv
  real(psb_spk_),intent(inout)              :: x(:)
  real(psb_spk_),intent(inout)              :: y(:)
  real(psb_spk_),intent(in)                 :: alpha,beta
  character(len=1),intent(in)               :: trans
  real(psb_spk_),target, intent(inout)      :: work(:)
  integer, intent(out)                      :: info

  Integer :: err_act
  character(len=20)  :: name='d_base_solver_apply'

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

end subroutine mld_s_base_solver_apply

subroutine mld_s_base_solver_apply_vect(alpha,sv,x,beta,y,desc_data,trans,work,info)
  use psb_base_mod
  use mld_s_base_solver_mod, mld_protect_name =>  mld_s_base_solver_apply_vect
  implicit none 
  type(psb_desc_type), intent(in)              :: desc_data
  class(mld_s_base_solver_type), intent(inout) :: sv
  type(psb_s_vect_type),intent(inout)          :: x
  type(psb_s_vect_type),intent(inout)          :: y
  real(psb_spk_),intent(in)                    :: alpha,beta
  character(len=1),intent(in)                  :: trans
  real(psb_spk_),target, intent(inout)         :: work(:)
  integer, intent(out)                         :: info

  Integer :: err_act
  character(len=20)  :: name='d_base_solver_apply'

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

end subroutine mld_s_base_solver_apply_vect


!
! Build
! The base version throws an error, since it means that no explicit
! choice was made. 
!
subroutine mld_s_base_solver_bld(a,desc_a,sv,upd,info,b,amold,vmold)
  use psb_base_mod
  use mld_s_base_solver_mod, mld_protect_name =>  mld_s_base_solver_bld
  Implicit None
  ! Arguments
  type(psb_sspmat_type), intent(in), target           :: a
  Type(psb_desc_type), Intent(in)                     :: desc_a 
  class(mld_s_base_solver_type), intent(inout)        :: sv
  character, intent(in)                               :: upd
  integer, intent(out)                                :: info
  type(psb_sspmat_type), intent(in), target, optional :: b
  class(psb_s_base_sparse_mat), intent(in), optional  :: amold
  class(psb_s_base_vect_type), intent(in), optional   :: vmold

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
end subroutine mld_s_base_solver_bld

subroutine mld_s_base_solver_check(sv,info)
  use psb_base_mod
  use mld_s_base_solver_mod, mld_protect_name =>  mld_s_base_solver_check
  Implicit None
  ! Arguments
  class(mld_s_base_solver_type), intent(inout) :: sv
  integer, intent(out)                   :: info
  Integer           :: err_act
  character(len=20) :: name='d_base_solver_check'

  call psb_erractionsave(err_act)
  info = psb_success_


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
end subroutine mld_s_base_solver_check

!
! Set.
! The base version does nothing; the principle is that
! SET acts on what is known, and delegates what is unknown.
! Since we are at the bottom of the hierarchy, there's no one
! to delegate, so we do nothing. 
!
subroutine mld_s_base_solver_seti(sv,what,val,info)
  use psb_base_mod
  use mld_s_base_solver_mod, mld_protect_name =>  mld_s_base_solver_seti
  Implicit None
  ! Arguments
  class(mld_s_base_solver_type), intent(inout) :: sv 
  integer, intent(in)                          :: what 
  integer, intent(in)                          :: val
  integer, intent(out)                         :: info
  Integer           :: err_act
  character(len=20) :: name='d_base_solver_seti'

  ! Correct action here is doing nothing. 
  info = 0

  return
end subroutine mld_s_base_solver_seti

subroutine mld_s_base_solver_setc(sv,what,val,info)
  use psb_base_mod
  use mld_s_base_solver_mod, mld_protect_name =>  mld_s_base_solver_setc
  Implicit None
  ! Arguments
  class(mld_s_base_solver_type), intent(inout) :: sv
  integer, intent(in)                          :: what 
  character(len=*), intent(in)                 :: val
  integer, intent(out)                         :: info
  Integer           :: err_act, ival 
  character(len=20) :: name='d_base_solver_setc'

  call psb_erractionsave(err_act)

  info = psb_success_

  call mld_stringval(val,ival,info)
  if (info == psb_success_) call sv%set(what,ival,info)

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
end subroutine mld_s_base_solver_setc

subroutine mld_s_base_solver_setr(sv,what,val,info)
  use psb_base_mod
  use mld_s_base_solver_mod, mld_protect_name =>  mld_s_base_solver_setr
  Implicit None
  ! Arguments
  class(mld_s_base_solver_type), intent(inout) :: sv 
  integer, intent(in)                          :: what 
  real(psb_spk_), intent(in)                   :: val
  integer, intent(out)                         :: info
  Integer           :: err_act
  character(len=20) :: name='d_base_solver_setr'


  ! Correct action here is doing nothing. 
  info = 0

  return
end subroutine mld_s_base_solver_setr

!
! Free
! The base version throws an error, since it means that no explicit
! choice was made. IS THIS CORRECT? I suspect it would be better
! to be silent here, to cover reallocation. 
!
subroutine mld_s_base_solver_free(sv,info)
  use psb_base_mod
  use mld_s_base_solver_mod, mld_protect_name =>  mld_s_base_solver_free
  Implicit None
  ! Arguments
  class(mld_s_base_solver_type), intent(inout) :: sv
  integer, intent(out)                         :: info
  Integer           :: err_act
  character(len=20) :: name='d_base_solver_free'

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
end subroutine mld_s_base_solver_free

subroutine mld_s_base_solver_descr(sv,info,iout,coarse)
  use psb_base_mod
  use mld_s_base_solver_mod, mld_protect_name =>  mld_s_base_solver_descr
  Implicit None
  ! Arguments
  class(mld_s_base_solver_type), intent(in) :: sv
  integer, intent(out)                      :: info
  integer, intent(in), optional             :: iout
  logical, intent(in), optional             :: coarse

  ! Local variables
  integer      :: err_act
  integer      :: ictxt, me, np
  character(len=20), parameter :: name='mld_s_base_solver_descr'
  integer      :: iout_


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
end subroutine mld_s_base_solver_descr

!
! Dump. For debugging purposes. 
!
subroutine mld_s_base_solver_dmp(sv,ictxt,level,info,prefix,head,solver)
  use psb_base_mod
  use mld_s_base_solver_mod, mld_protect_name =>  mld_s_base_solver_dmp
  implicit none 
  class(mld_s_base_solver_type), intent(in) :: sv
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
    prefix_ = "dump_slv_s"
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

end subroutine mld_s_base_solver_dmp

