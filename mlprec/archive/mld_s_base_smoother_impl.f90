!  
!   
!                             MLD2P4  version 2.0
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.3)
!    
!    (C) Copyright 2008, 2010, 2012, 2015
!  
!                        Salvatore Filippone  University of Rome Tor Vergata
!                        Alfredo Buttari      CNRS-IRIT, Toulouse
!                        Pasqua D'Ambra       ICAR-CNR, Naples
!                        Daniela di Serafino  Second University of Naples
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the MLD2P4 group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!  

!
! Apply: comes in two versions, on plain arrays or on encapsulated
! vectors.
! This basic version just applies the local solver, whatever that
! is. 
!

subroutine mld_s_base_smoother_apply(alpha,sm,x,beta,y,desc_data,trans,sweeps,work,info)
  use psb_base_mod
  use mld_s_base_smoother_mod, mld_protect_name =>  mld_s_base_smoother_apply
  implicit none 
  type(psb_desc_type), intent(in)             :: desc_data
  class(mld_s_base_smoother_type), intent(in) :: sm
  real(psb_spk_),intent(inout)                :: x(:)
  real(psb_spk_),intent(inout)                :: y(:)
  real(psb_spk_),intent(in)                   :: alpha,beta
  character(len=1),intent(in)                 :: trans
  integer, intent(in)                         :: sweeps
  real(psb_spk_),target, intent(inout)        :: work(:)
  integer, intent(out)                        :: info

  Integer           :: err_act
  character(len=20) :: name='s_base_smoother_apply'

  call psb_erractionsave(err_act)
  info = psb_success_
  if (allocated(sm%sv)) then 
    call sm%sv%apply(alpha,x,beta,y,desc_data,trans,work,info)
  else
    info = 1121
  endif
  if (info /= psb_success_) then 
    call psb_errpush(info,name)
    goto 9999 
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

end subroutine mld_s_base_smoother_apply

subroutine mld_s_base_smoother_apply_vect(alpha,sm,x,beta,y,desc_data,&
     &  trans,sweeps,work,info)
  use psb_base_mod
  use mld_s_base_smoother_mod, mld_protect_name =>  mld_s_base_smoother_apply_vect
  implicit none 
  type(psb_desc_type), intent(in)                :: desc_data
  class(mld_s_base_smoother_type), intent(inout) :: sm
  type(psb_s_vect_type),intent(inout)            :: x
  type(psb_s_vect_type),intent(inout)            :: y
  real(psb_spk_),intent(in)                      :: alpha,beta
  character(len=1),intent(in)                    :: trans
  integer, intent(in)                            :: sweeps
  real(psb_spk_),target, intent(inout)           :: work(:)
  integer, intent(out)                           :: info

  Integer           :: err_act
  character(len=20) :: name='s_base_smoother_apply'

  call psb_erractionsave(err_act)
  info = psb_success_
  if (allocated(sm%sv)) then 
    call sm%sv%apply(alpha,x,beta,y,desc_data,trans,work,info)
  else
    info = 1121
  endif
  if (info /= psb_success_) then 
    call psb_errpush(info,name)
    goto 9999 
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

end subroutine mld_s_base_smoother_apply_vect

!
! Check:
! 1. Check that we do have a solver object
! 2. Call its check method
!

subroutine mld_s_base_smoother_check(sm,info)
  use psb_base_mod
  use mld_s_base_smoother_mod, mld_protect_name =>  mld_s_base_smoother_check
  Implicit None

  ! Arguments
  class(mld_s_base_smoother_type), intent(inout) :: sm 
  integer, intent(out)                   :: info
  Integer           :: err_act
  character(len=20) :: name='s_base_smoother_check'

  call psb_erractionsave(err_act)
  info = psb_success_

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
end subroutine mld_s_base_smoother_check

!
! Set methods: the come in multiple versions according
! to whether we are setting with integer, real or character
! input.
! The basic rule is: if the input refers to a parameter
! of the smoother, use it, otherwise pass it to the
! solver object for further processing.
! Since there are no parameters in the base smoother
! we just pass everything to the solver object. 
!
subroutine mld_s_base_smoother_seti(sm,what,val,info)
  use psb_base_mod
  use mld_s_base_smoother_mod, mld_protect_name =>  mld_s_base_smoother_seti
  Implicit None
  ! Arguments
  class(mld_s_base_smoother_type), intent(inout) :: sm 
  integer, intent(in)                            :: what 
  integer, intent(in)                            :: val
  integer, intent(out)                           :: info
  Integer           :: err_act
  character(len=20) :: name='s_base_smoother_seti'

  call psb_erractionsave(err_act)
  info = psb_success_

  if (allocated(sm%sv)) then 
    call sm%sv%set(what,val,info)
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
end subroutine mld_s_base_smoother_seti

subroutine mld_s_base_smoother_setc(sm,what,val,info)
  use psb_base_mod
  use mld_s_base_smoother_mod, mld_protect_name =>  mld_s_base_smoother_setc
  Implicit None

  ! Arguments
  class(mld_s_base_smoother_type), intent(inout) :: sm 
  integer, intent(in)                            :: what 
  character(len=*), intent(in)                   :: val
  integer, intent(out)                           :: info
  Integer           :: err_act
  character(len=20) :: name='s_base_smoother_setc'

  call psb_erractionsave(err_act)

  info = psb_success_

  if (allocated(sm%sv)) then 
    call sm%sv%set(what,val,info)
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
end subroutine mld_s_base_smoother_setc

subroutine mld_s_base_smoother_setr(sm,what,val,info)
  use psb_base_mod
  use mld_s_base_smoother_mod, mld_protect_name =>  mld_s_base_smoother_setr
  Implicit None

  ! Arguments
  class(mld_s_base_smoother_type), intent(inout) :: sm 
  integer, intent(in)                            :: what 
  real(psb_spk_), intent(in)                     :: val
  integer, intent(out)                           :: info
  Integer :: err_act
  character(len=20)  :: name='s_base_smoother_setr'

  call psb_erractionsave(err_act)


  info = psb_success_

  if (allocated(sm%sv)) then 
    call sm%sv%set(what,val,info)
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
end subroutine mld_s_base_smoother_setr


!
! Build method.
! At base level we only have to pass data to the inner solver. 
! AMOLD/VMOLD allow to have any relevant sparse matrix or vector
! to be stored in a given format. This is essential e.g.
! when dealing  with GPUs. 
!
subroutine mld_s_base_smoother_bld(a,desc_a,sm,upd,info,amold,vmold)
  use psb_base_mod
  use mld_s_base_smoother_mod, mld_protect_name =>  mld_s_base_smoother_bld
  Implicit None

  ! Arguments
  type(psb_sspmat_type), intent(in), target      :: a
  Type(psb_desc_type), Intent(in)                :: desc_a 
  class(mld_s_base_smoother_type), intent(inout) :: sm 
  character, intent(in)                          :: upd
  integer, intent(out)                           :: info
  class(psb_s_base_sparse_mat), intent(in), optional :: amold
  class(psb_s_base_vect_type), intent(in), optional  :: vmold
  Integer           :: err_act
  character(len=20) :: name='s_base_smoother_bld'

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
end subroutine mld_s_base_smoother_bld

!
! Free method (aka destructor).
! In most cases we could do without; however
! for cases where there are data objects allocated outside
! of the Fortran RTE we need to free them explicitly.
!
! Even in that case, we could do without this if FINAL
! subroutines were supported, which is not the case

! in GNU Fortran up to 4.7. 
!
subroutine mld_s_base_smoother_free(sm,info)
  use psb_base_mod
  use mld_s_base_smoother_mod, mld_protect_name => mld_s_base_smoother_free
  Implicit None

  ! Arguments
  class(mld_s_base_smoother_type), intent(inout) :: sm
  integer, intent(out)                           :: info
  Integer           :: err_act
  character(len=20) :: name='s_base_smoother_free'

  call psb_erractionsave(err_act)
  info = psb_success_

  if (allocated(sm%sv)) then 
    call sm%sv%free(info)
    if (info == psb_success_) deallocate(sm%sv,stat=info) 
  end if
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
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
end subroutine mld_s_base_smoother_free

!
! Print a description
!

subroutine mld_s_base_smoother_descr(sm,info,iout,coarse)
  use psb_base_mod
  use mld_s_base_smoother_mod, mld_protect_name =>  mld_s_base_smoother_descr
  Implicit None

  ! Arguments
  class(mld_s_base_smoother_type), intent(in) :: sm
  integer, intent(out)                        :: info
  integer, intent(in), optional               :: iout
  logical, intent(in), optional               :: coarse

  ! Local variables
  integer      :: err_act
  integer      :: ictxt, me, np
  character(len=20), parameter :: name='mld_s_base_smoother_descr'
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
end subroutine mld_s_base_smoother_descr

!
! Dump 
! to file, for debugging purposes.
!
subroutine mld_s_base_smoother_dmp(sm,ictxt,level,info,prefix,head,smoother,solver)
  use psb_base_mod
  use mld_s_base_smoother_mod, mld_protect_name =>  mld_s_base_smoother_dmp
  implicit none 
  class(mld_s_base_smoother_type), intent(in) :: sm
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
    prefix_ = "dump_smth_s"
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

  ! At base level do nothing for the smoother
  if (allocated(sm%sv)) &
       & call sm%sv%dump(ictxt,level,info,solver=solver,prefix=prefix)

end subroutine mld_s_base_smoother_dmp

