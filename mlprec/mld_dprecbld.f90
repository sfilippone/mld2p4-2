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
! File: mld_dprecbld.f90
!
! Subroutine: mld_dprecbld
! Version:    real
! Contains:   subroutine init_baseprec_av
!
!  This routine builds the preconditioner according to the requirements made by
!  the user through the subroutines mld_precinit and mld_precset. 
! 
!
! Arguments:
!    a       -  type(psb_dspmat_type).
!               The sparse matrix structure containing the local part of the
!               matrix to be preconditioned.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    p       -  type(mld_dprec_type), input/output.
!               The preconditioner data structure containing the local part
!               of the preconditioner to be built.
!    info    -  integer, output.
!               Error code.              
!  
subroutine mld_dprecbld(a,desc_a,p,info)

  use psb_sparse_mod
  use mld_inner_mod
  use mld_prec_mod, mld_protect_name => mld_dprecbld
  use mld_d_jac_smoother
  use mld_d_as_smoother
  use mld_d_diag_solver
  use mld_d_ilu_solver
  
  Implicit None

  ! Arguments
  type(psb_dspmat_type),intent(in), target :: a
  type(psb_desc_type), intent(in), target   :: desc_a
  type(mld_dprec_type),intent(inout),target :: p
  integer, intent(out)                      :: info
!!$  character, intent(in), optional         :: upd

  ! Local Variables
  type(mld_dprec_type)    :: t_prec
  Integer      :: err,i,k,ictxt, me,np, err_act, iszv, newsz
  integer      :: ipv(mld_ifpsz_), val
  integer      :: int_err(5)
  character    :: upd_
  integer            :: debug_level, debug_unit
  character(len=20)  :: name, ch_err

  if (psb_get_errstatus().ne.0) return 
  info=psb_success_
  err=0
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  name = 'mld_dprecbld'
  info = psb_success_
  int_err(1) = 0
  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  p%ictxt = ictxt

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Entering ',desc_a%matrix_data(:)
  !
  ! For the time being we are commenting out the UPDATE argument;
  ! we plan to resurrect it later. 
!!$  if (present(upd)) then 
!!$    if (debug_level >= psb_debug_outer_) &
!!$         & write(debug_unit,*) me,' ',trim(name),'UPD ', upd
!!$
!!$    if ((psb_toupper(upd).eq.'F').or.(psb_toupper(upd).eq.'T')) then
!!$      upd_=psb_toupper(upd)
!!$    else
!!$      upd_='F'
!!$    endif
!!$  else
!!$    upd_='F'
!!$  endif
  upd_ = 'F'

  if (.not.allocated(p%precv)) then 
    !! Error: should have called mld_dprecinit
    info=3111
    call psb_errpush(info,name)
    goto 9999
  end if

  !
  ! Check to ensure all procs have the same 
  !   
  newsz = -1
  iszv  = size(p%precv)
  call psb_bcast(ictxt,iszv)
  if (iszv /= size(p%precv)) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent size of precv')
    goto 9999
  end if
  
  if (iszv <= 0) then 
    ! Is this really possible? probably not.
    info=psb_err_from_subroutine_
    ch_err='size bpv'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  !
  ! Number of levels = 1
  !
  if (iszv == 1) then 
    !
    ! Check on the iprcparm contents: they should be the same
    ! on all processes.
    !
    if (me == psb_root_) ipv(:) = p%precv(1)%iprcparm(:) 
    call psb_bcast(ictxt,ipv) 
    if (any(ipv(:) /=  p%precv(1)%iprcparm(:) )) then
      write(debug_unit,*) me,name,&
           &': Inconsistent arguments among processes, forcing a default'
      p%precv(1)%iprcparm(:) = ipv(:) 
    end if
    !
    ! Remember to fix base_a and base_desc
    ! 
    call init_baseprec_av(p%precv(1)%prec,info)
    p%precv(1)%base_a    => a
    p%precv(1)%base_desc => desc_a

    if (info /= psb_success_) then 
      call psb_errpush(psb_err_internal_error_,name,a_err='Base level precbuild.')
      goto 9999
    end if
    !
    ! Build the base preconditioner
    !
    select case(p%precv(1)%prec%iprcparm(mld_sub_solve_))
    case(mld_ilu_n_,mld_milu_n_)      
      call mld_check_def(p%precv(1)%prec%iprcparm(mld_sub_fillin_),&
           & 'Level',0,is_legal_ml_lev)
    case(mld_ilu_t_)                 
      call mld_check_def(p%precv(1)%prec%rprcparm(mld_sub_iluthrs_),&
           & 'Eps',dzero,is_legal_fact_thrs)
    end select
    call mld_check_def(p%precv(1)%iprcparm(mld_smoother_sweeps_),&
         & 'Jacobi sweeps',1,is_legal_jac_sweeps)

    !
    !  Test version for beginning of OO stuff. 
    ! 
    if (allocated(p%precv(1)%sm)) then 
      call p%precv(1)%sm%free(info)
      if (info == psb_success_) deallocate(p%precv(1)%sm,stat=info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_alloc_dealloc_,name,a_err='One level preconditioner build.')
        goto 9999
      endif
    end if
    select case (p%precv(1)%prec%iprcparm(mld_smoother_type_)) 
    case(mld_jac_, mld_bjac_) 
      allocate(mld_d_jac_smoother_type :: p%precv(1)%sm, stat=info)
    case(mld_as_)
      allocate(mld_d_as_smoother_type  :: p%precv(1)%sm, stat=info)
    case default
      info = -1 
    end select
    if (info /= psb_success_) then 
      write(0,*) ' Smoother allocation error',info,&
           & p%precv(1)%prec%iprcparm(mld_smoother_type_)
      call psb_errpush(psb_err_internal_error_,name,a_err='One level preconditioner build.')
      goto 9999
    endif
    call p%precv(1)%sm%set(mld_sub_restr_,p%precv(1)%prec%iprcparm(mld_sub_restr_),info)
    call p%precv(1)%sm%set(mld_sub_prol_,p%precv(1)%prec%iprcparm(mld_sub_prol_),info)
    call p%precv(1)%sm%set(mld_sub_ovr_,p%precv(1)%prec%iprcparm(mld_sub_ovr_),info)

    select case (p%precv(1)%prec%iprcparm(mld_sub_solve_)) 
    case(mld_ilu_n_,mld_milu_n_,mld_ilu_t_) 
      allocate(mld_d_ilu_solver_type :: p%precv(1)%sm%sv, stat=info)
      if (info == psb_success_) call  p%precv(1)%sm%sv%set(mld_sub_solve_,&
           & p%precv(1)%prec%iprcparm(mld_sub_solve_),info)
      if (info == psb_success_) call  p%precv(1)%sm%sv%set(mld_sub_fillin_,&
           & p%precv(1)%prec%iprcparm(mld_sub_fillin_),info)
      if (info == psb_success_) call  p%precv(1)%sm%sv%set(mld_sub_iluthrs_,&
           & p%precv(1)%prec%rprcparm(mld_sub_iluthrs_),info)
    case(mld_diag_scale_)
      allocate(mld_d_diag_solver_type :: p%precv(1)%sm%sv, stat=info)
    case default
      info = -1 
    end select

    if (info /= psb_success_) then 
      write(0,*) ' Solver allocation error',info,&
           & p%precv(1)%prec%iprcparm(mld_sub_solve_)
      call psb_errpush(psb_err_internal_error_,name,a_err='One level preconditioner build.')
      goto 9999
    endif

    call p%precv(1)%sm%build(a,desc_a,upd_,info)
    if (info /= psb_success_) then 
      write(0,*) ' Smoother build error',info
      call psb_errpush(psb_err_internal_error_,name,a_err='One level preconditioner build.')
      goto 9999
    endif
      
  !
  ! Number of levels > 1
  !
  else if (iszv > 1) then
    !
    ! Build the multilevel preconditioner
    ! 
    call  mld_mlprec_bld(a,desc_a,p,info)
    
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_internal_error_,name,a_err='Multilevel preconditioner build.')
      goto 9999
    endif
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

contains

  subroutine init_baseprec_av(p,info)
    type(mld_dbaseprec_type), intent(inout) :: p
    integer                                :: info
!!$    if (allocated(p%av)) then
!!$      if (size(p%av) /= mld_max_avsz_) then 
!!$        deallocate(p%av,stat=info)
!!$        if (info /= psb_success_) return 
!!$      endif
!!$    end if
!!$    if (.not.(allocated(p%av))) then 
!!$      allocate(p%av(mld_max_avsz_),stat=info)
!!$      if (info /= psb_success_) return
!!$    end if
!!$    do k=1,size(p%av)
!!$      call psb_nullify_sp(p%av(k))
!!$    end do

  end subroutine init_baseprec_av

end subroutine mld_dprecbld
