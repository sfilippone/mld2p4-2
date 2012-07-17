!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010
!!$
!!$                      Salvatore Filippone  University of Rome Tor Vergata
!!$                      Alfredo Buttari      CNRS-IRIT, Toulouse
!!$                      Pasqua D'Ambra       ICAR-CNR, Naples
!!$                      Daniela di Serafino  Second University of Naples
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the MLD2P4 group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$

!
! Subroutine: mld_file_onelev_descr
! Version: real
!
!  This routine prints a description of the preconditioner to the standard 
!  output or to a file. It must be called after the preconditioner has been
!  built by mld_precbld.
!
! Arguments:
!  p       -  type(mld_Tprec_type), input.
!             The preconditioner data structure to be printed out.
!  info    -  integer, output.
!             error code.
!  iout    -  integer, input, optional.
!             The id of the file where the preconditioner description
!             will be printed. If iout is not present, then the standard
!             output is condidered.
!
subroutine mld_s_base_onelev_descr(lv,il,nl,info,iout)
  use psb_base_mod
  use mld_s_onelev_mod, mld_protect_name => mld_s_base_onelev_descr
  Implicit None
  ! Arguments
  class(mld_s_onelev_type), intent(in) :: lv
  integer, intent(in)                 :: il,nl
  integer, intent(out)                :: info
  integer, intent(in), optional       :: iout

  ! Local variables
  integer      :: err_act
  integer      :: ictxt, me, np
  character(len=20), parameter :: name='mld_s_base_onelev_descr'
  integer      :: iout_
  logical      :: coarse


  call psb_erractionsave(err_act)


  coarse = (il==nl)

  if (present(iout)) then 
    iout_ = iout
  else 
    iout_ = 6
  end if

  write(iout_,*) 
  if (il == 2) then 
    call lv%parms%mldescr(iout_,info)
    write(iout_,*) 
  end if

  if (coarse)  then 
    write(iout_,*) ' Level ',il,' (coarsest)'
  else
    write(iout_,*) ' Level ',il
  end if

  call lv%parms%descr(iout_,info,coarse=coarse)

  if (nl > 1) then 
    if (allocated(lv%map%naggr)) then
      write(iout_,*) '  Size of coarse matrix: ', &
           &  sum(lv%map%naggr(:))
      write(iout_,*) '  Sizes of aggregates: ', &
           &  lv%map%naggr(:)
    end if
  end if

  if (coarse.and.allocated(lv%sm)) &
       & call lv%sm%descr(info,iout=iout_,coarse=coarse)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine mld_s_base_onelev_descr


!
! Subroutines: mld_T_onelev_precfree
! Version: real
!
!  These routines deallocate the mld_Tonelev_type 
!
! Arguments:
!  p       -  type(mld_Tonelev_type), input.
!             The data structure to be deallocated.
!  info    -  integer, output.
!             error code.
!
subroutine mld_s_base_onelev_free(lv,info)
  use psb_base_mod
  use mld_s_onelev_mod, mld_protect_name =>  mld_s_base_onelev_free
  implicit none 

  class(mld_s_onelev_type), intent(inout) :: lv
  integer, intent(out)                :: info
  integer :: i

  info = psb_success_

  ! We might just deallocate the top level array, except 
  ! that there may be inner objects containing C pointers,
  ! e.g.  UMFPACK, SLU or CUDA stuff.
  ! We really need FINALs. 
  if (allocated(lv%sm)) &
       & call lv%sm%free(info)

  call lv%ac%free()
  if (lv%desc_ac%is_ok()) &
       & call lv%desc_ac%free(info)
  call lv%map%free(info)

  ! This is a pointer to something else, must not free it here. 
  nullify(lv%base_a) 
  ! This is a pointer to something else, must not free it here. 
  nullify(lv%base_desc) 

  call lv%nullify()

end subroutine mld_s_base_onelev_free

!
! Onelevel checks. 
! The number of Jacobi sweeps to be applied is not
! tied to the Jacobi smoother: logically, you have
! a smoother and you can choose to apply it any number
! of times you like. 
!
subroutine mld_s_base_onelev_check(lv,info)
  use psb_base_mod
  use mld_s_onelev_mod, mld_protect_name => mld_s_base_onelev_check

  Implicit None

  ! Arguments
  class(mld_s_onelev_type), intent(inout) :: lv 
  integer, intent(out)                   :: info
  Integer           :: err_act
  character(len=20) :: name='s_base_onelev_check'

  call psb_erractionsave(err_act)
  info = psb_success_

  call mld_check_def(lv%parms%sweeps,&
       & 'Jacobi sweeps',1,is_legal_jac_sweeps)
  call mld_check_def(lv%parms%sweeps_pre,&
       & 'Jacobi sweeps',1,is_legal_jac_sweeps)
  call mld_check_def(lv%parms%sweeps_post,&
       & 'Jacobi sweeps',1,is_legal_jac_sweeps)


  if (allocated(lv%sm)) then 
    call lv%sm%check(info)
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
end subroutine mld_s_base_onelev_check


! 
! Set routines:
! Parameters belonging here are:
!   Number of smoothing sweeps;
!   Smoother position;
!   Aggregation related parameters
!   Record request on coarse level solver,
!    for checks on solver vs. smoother nomenclature
!    reconciliation. 
! 
subroutine mld_s_base_onelev_seti(lv,what,val,info)
  use psb_base_mod
  use mld_s_onelev_mod, mld_protect_name => mld_s_base_onelev_seti

  Implicit None

  ! Arguments
  class(mld_s_onelev_type), intent(inout) :: lv 
  integer, intent(in)                          :: what 
  integer, intent(in)                          :: val
  integer, intent(out)                         :: info
  Integer           :: err_act
  character(len=20) :: name='s_base_onelev_seti'

  call psb_erractionsave(err_act)
  info = psb_success_

  select case (what) 

  case (mld_smoother_sweeps_)
    lv%parms%sweeps      = val
    lv%parms%sweeps_pre  = val
    lv%parms%sweeps_post = val

  case (mld_smoother_sweeps_pre_)
    lv%parms%sweeps_pre  = val

  case (mld_smoother_sweeps_post_)
    lv%parms%sweeps_post = val

  case (mld_ml_type_)
    lv%parms%ml_type       = val

  case (mld_aggr_alg_)
    lv%parms%aggr_alg      = val

  case (mld_aggr_kind_)
    lv%parms%aggr_kind     = val

  case (mld_coarse_mat_)
    lv%parms%coarse_mat    = val

  case (mld_smoother_pos_)
    lv%parms%smoother_pos  = val

  case (mld_aggr_omega_alg_)
    lv%parms%aggr_omega_alg= val

  case (mld_aggr_eig_)
    lv%parms%aggr_eig      = val

  case (mld_aggr_filter_)
    lv%parms%aggr_filter   = val

  case (mld_coarse_solve_)
    lv%parms%coarse_solve    = val

  case default
    if (allocated(lv%sm)) then 
      call lv%sm%set(what,val,info)
    end if
    if (info /= psb_success_) goto 9999
  end select
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine mld_s_base_onelev_seti

subroutine mld_s_base_onelev_setc(lv,what,val,info)
  use psb_base_mod
  use mld_s_onelev_mod, mld_protect_name => mld_s_base_onelev_setc

  Implicit None

  ! Arguments
  class(mld_s_onelev_type), intent(inout) :: lv 
  integer, intent(in)                            :: what 
  character(len=*), intent(in)                   :: val
  integer, intent(out)                           :: info
  Integer           :: err_act
  character(len=20) :: name='s_base_onelev_setc'
  integer :: ival 

  call psb_erractionsave(err_act)

  info = psb_success_

  call mld_stringval(val,ival,info)
  if (info == psb_success_) call lv%set(what,ival,info)

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
end subroutine mld_s_base_onelev_setc

subroutine mld_s_base_onelev_setr(lv,what,val,info)
  use psb_base_mod
  use mld_s_onelev_mod, mld_protect_name => mld_s_base_onelev_setr

  Implicit None

  ! Arguments
  class(mld_s_onelev_type), intent(inout) :: lv 
  integer, intent(in)                            :: what 
  real(psb_spk_), intent(in)                     :: val
  integer, intent(out)                           :: info
  Integer :: err_act
  character(len=20)  :: name='s_base_onelev_setr'

  call psb_erractionsave(err_act)


  info = psb_success_

  select case (what) 

  case (mld_aggr_omega_val_)
    lv%parms%aggr_omega_val= val

  case (mld_aggr_thresh_)
    lv%parms%aggr_thresh   = val

  case default
    if (allocated(lv%sm)) then 
      call lv%sm%set(what,val,info)
    end if
    if (info /= psb_success_) goto 9999
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine mld_s_base_onelev_setr

!
! Dump on file: can be fine-tuned to include the (aggregated) matrix
! as well as smoother and solver. 
!
subroutine mld_s_base_onelev_dump(lv,level,info,prefix,head,ac,rp,smoother,solver)
  use psb_base_mod
  use mld_s_onelev_mod, mld_protect_name => mld_s_base_onelev_dump
  implicit none 
  class(mld_s_onelev_type), intent(in) :: lv
  integer, intent(in)              :: level
  integer, intent(out)             :: info
  character(len=*), intent(in), optional :: prefix, head
  logical, optional, intent(in)    :: ac, rp, smoother, solver
  integer :: i, j, il1, iln, lname, lev
  integer :: icontxt,iam, np
  character(len=80)  :: prefix_
  character(len=120) :: fname ! len should be at least 20 more than
  logical :: ac_, rp_
  !  len of prefix_ 

  info = 0

  if (present(prefix)) then 
    prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
  else
    prefix_ = "dump_lev_s"
  end if

  if (associated(lv%base_desc)) then 
    icontxt = lv%base_desc%get_context()
    call psb_info(icontxt,iam,np)
  else 
    icontxt = -1 
    iam     = -1
  end if
  if (present(ac)) then 
    ac_ = ac
  else
    ac_ = .false. 
  end if
  if (present(rp)) then 
    rp_ = rp
  else
    rp_ = .false. 
  end if
  lname = len_trim(prefix_)
  fname = trim(prefix_)
  write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
  lname = lname + 5

  if (level >= 2) then 
    if (ac_) then 
      write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_ac.mtx'
      write(0,*) 'Filename ',fname
      call lv%ac%print(fname,head=head)
    end if
    if (rp_) then 
      write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_r.mtx'
      write(0,*) 'Filename ',fname
      call lv%map%map_X2Y%print(fname,head=head)
      write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_p.mtx'
      write(0,*) 'Filename ',fname
      call lv%map%map_Y2X%print(fname,head=head)
    end if
  end if
  if (allocated(lv%sm)) &
       & call lv%sm%dump(icontxt,level,info,smoother=smoother, &
       & solver=solver,prefix=prefix)

end subroutine mld_s_base_onelev_dump


