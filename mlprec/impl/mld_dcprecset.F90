!!$ 
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.3)
!!$  
!!$  (C) Copyright 2008, 2010, 2012, 2015
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
! File: mld_dprecset.f90
!
! Subroutine: mld_dprecseti
! Version: real
!
!  This routine sets the integer parameters defining the preconditioner. More
!  precisely, the integer parameter identified by 'what' is assigned the value
!  contained in 'val'.
!  For the multilevel preconditioners, the levels are numbered in increasing
!  order starting from the finest one, i.e. level 1 is the finest level. 
!
!  To set character and real parameters, see mld_dprecsetc and mld_dprecsetr,
!  respectively.
!
!
! Arguments:
!    p       -  type(mld_dprec_type), input/output.
!               The preconditioner data structure.
!    what    -  integer, input.
!               The number identifying the parameter to be set.
!               A mnemonic constant has been associated to each of these
!               numbers, as reported in the MLD2P4 User's and Reference Guide.
!    val     -  integer, input.
!               The value of the parameter to be set. The list of allowed
!               values is reported in the MLD2P4 User's and Reference Guide.
!    info    -  integer, output.
!               Error code.
!    ilev    -  integer, optional, input.
!               For the multilevel preconditioner, the level at which the
!               preconditioner parameter has to be set. 
!               If nlev is not present, the parameter identified by 'what'
!               is set at all the appropriate levels.
!
!  NOTE: currently, the use of the argument ilev is not "safe" and is reserved to
!  MLD2P4 developers. Indeed, by using ilev it is possible to set different values
!  of the same parameter at different levels 1,...,nlev-1, even in cases where
!  the parameter must have the same value at all the levels but the coarsest one.
!  For this reason, the interface mld_precset to this routine has been built in
!  such a way that ilev is not visible to the user (see mld_prec_mod.f90).
!   
subroutine mld_dcprecseti(p,what,val,info,ilev)

  use psb_base_mod
  use mld_d_prec_mod, mld_protect_name => mld_dcprecseti
  use mld_d_jac_smoother
  use mld_d_as_smoother
  use mld_d_diag_solver
  use mld_d_ilu_solver
  use mld_d_id_solver
  use mld_d_gs_solver
#if defined(HAVE_UMF_)
  use mld_d_umf_solver
#endif
#if defined(HAVE_SLUDIST_)
  use mld_d_sludist_solver
#endif
#if defined(HAVE_SLU_)
  use mld_d_slu_solver
#endif
#if defined(HAVE_MUMPS_)  
  use mld_d_mumps_solver
#endif

  implicit none

  ! Arguments
  class(mld_dprec_type), intent(inout)    :: p
  character(len=*), intent(in)            :: what 
  integer(psb_ipk_), intent(in)           :: val
  integer(psb_ipk_), intent(out)          :: info
  integer(psb_ipk_), optional, intent(in) :: ilev

  ! Local variables
  integer(psb_ipk_)                      :: ilev_, nlev_
  character(len=*), parameter            :: name='mld_precseti'

  info = psb_success_
  if (.not.allocated(p%precv)) then 
    info = 3111
    write(psb_err_unit,*) name,&
         & ': Error: uninitialized preconditioner,',&
         &' should call MLD_PRECINIT'
    return 
  endif
  nlev_ = size(p%precv)

  if (present(ilev)) then 
    ilev_ = ilev
  else
    ilev_ = 1 
  end if

  if ((ilev_<1).or.(ilev_ > nlev_)) then 
    info = -1
    write(psb_err_unit,*) name,&
         &': Error: invalid ILEV/NLEV combination',ilev_, nlev_
    return
  endif

  if (psb_toupper(what) == 'COARSE_AGGR_SIZE') then 
    p%coarse_aggr_size = max(val,-1)
    return
  end if
  !
  ! Set preconditioner parameters at level ilev.
  !
  if (present(ilev)) then 

    if (ilev_ == 1) then
      ! 
      ! Rules for fine level are slightly different.
      ! 
      select case(psb_toupper(trim(what))) 
      case('SMOOTHER_TYPE')
        call onelev_set_smoother(p%precv(ilev_),val,info)
      case('SUB_SOLVE')
        call onelev_set_solver(p%precv(ilev_),val,info)
      case('SMOOTHER_SWEEPS','ML_TYPE','AGGR_ALG','AGGR_KIND',&
           & 'SMOOTHER_POS','AGGR_OMEGA_ALG','AGGR_EIG',&
           & 'SMOOTHER_SWEEPS_PRE','SMOOTHER_SWEEPS_POST',&
           & 'SUB_RESTR','SUB_PROL', &
           & 'SUB_REN','SUB_OVR','SUB_FILLIN')
        call p%precv(ilev_)%set(what,val,info)

      case default
        call p%precv(ilev_)%set(what,val,info)
      end select

    else if (ilev_ > 1) then 

      select case(psb_toupper(what)) 
      case('SMOOTHER_TYPE')
        call onelev_set_smoother(p%precv(ilev_),val,info)
      case('SUB_SOLVE')
        call onelev_set_solver(p%precv(ilev_),val,info)
      case('SMOOTHER_SWEEPS','ML_TYPE','AGGR_ALG','AGGR_KIND',&
           & 'SMOOTHER_POS','AGGR_OMEGA_ALG','AGGR_EIG',&
           & 'SMOOTHER_SWEEPS_PRE','SMOOTHER_SWEEPS_POST',&
           & 'SUB_RESTR','SUB_PROL', &
           & 'SUB_REN','SUB_OVR','SUB_FILLIN',&
           & 'COARSE_MAT')
        call p%precv(ilev_)%set(what,val,info)

      case('COARSE_SUBSOLVE')
        if (ilev_ /= nlev_) then 
          write(psb_err_unit,*) name,&
               & ': Error: Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if
        call onelev_set_solver(p%precv(ilev_),val,info)
      case('COARSE_SOLVE')
        if (ilev_ /= nlev_) then 
          write(psb_err_unit,*) name,&
               & ': Error: Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if

        if (nlev_ > 1) then 
          call p%precv(nlev_)%set('COARSE_SOLVE',val,info)
          select case (val) 
          case(mld_bjac_)
            call onelev_set_smoother(p%precv(nlev_),val,info)
#if defined(HAVE_UMF_)
            call onelev_set_solver(p%precv(nlev_),mld_umf_,info)
#elif defined(HAVE_SLU_) 
            call onelev_set_solver(p%precv(nlev_),mld_slu_,info)
#elif defined(HAVE_MUMPS_)
            call onelev_set_solver(p%precv(nlev_),mld_mumps_,info)
#else 
            call onelev_set_solver(p%precv(nlev_),mld_ilu_n_,info)
#endif
            call p%precv(nlev_)%set('COARSE_MAT',mld_distr_mat_,info)
          case(mld_umf_, mld_slu_,mld_ilu_n_, mld_ilu_t_,mld_milu_n_)
            call onelev_set_smoother(p%precv(nlev_),mld_bjac_,info)
            call onelev_set_solver(p%precv(nlev_),val,info)
            call p%precv(nlev_)%set('COARSE_MAT',mld_repl_mat_,info)
          case(mld_sludist_)
            call onelev_set_smoother(p%precv(nlev_),mld_bjac_,info)
            call onelev_set_solver(p%precv(nlev_),val,info)
            call p%precv(nlev_)%set('COARSE_MAT',mld_distr_mat_,info)
          case(mld_mumps_)
            call onelev_set_smoother(p%precv(nlev_),mld_bjac_,info)
            call onelev_set_solver(p%precv(nlev_),val,info)
            call p%precv(nlev_)%set('COARSE_MAT',mld_distr_mat_,info)
          case(mld_jac_)
            call onelev_set_smoother(p%precv(nlev_),mld_bjac_,info)
            call onelev_set_solver(p%precv(nlev_),mld_diag_scale_,info)
            call p%precv(nlev_)%set('COARSE_MAT',mld_distr_mat_,info)
          end select

        endif
      case('COARSE_SWEEPS')
        if (ilev_ /= nlev_) then 
          write(psb_err_unit,*) name,&
               & ': Error: Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if
        call p%precv(nlev_)%set('SMOOTHER_SWEEPS',val,info)

      case('COARSE_FILLIN')
        if (ilev_ /= nlev_) then 
          write(psb_err_unit,*) name,&
               & ': Error: Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if
        call p%precv(nlev_)%set('SUB_FILLIN',val,info)
      case default
        call p%precv(ilev_)%set(what,val,info)
      end select

    endif

  else if (.not.present(ilev)) then 
    !
    ! ilev not specified: set preconditioner parameters at all the appropriate
    ! levels
    !
    select case(psb_toupper(trim(what))) 
    case('SUB_SOLVE')
      do ilev_=1,max(1,nlev_-1)
        if (.not.allocated(p%precv(ilev_)%sm)) then 
          write(psb_err_unit,*) name,&
               & ': Error: uninitialized preconditioner component,',&
               & ' should call MLD_PRECINIT' 
          info = -1 
          return 
        endif
        call onelev_set_solver(p%precv(ilev_),val,info)

      end do

    case('SUB_RESTR','SUB_PROL',&
         & 'SUB_REN','SUB_OVR','SUB_FILLIN')
      do ilev_=1,max(1,nlev_-1)
        call p%precv(ilev_)%set(what,val,info)
      end do

    case('SMOOTHER_SWEEPS')
      do ilev_=1,max(1,nlev_-1)
        call p%precv(ilev_)%set(what,val,info)
      end do

    case('SMOOTHER_TYPE')
      do ilev_=1,max(1,nlev_-1)
        call onelev_set_smoother(p%precv(ilev_),val,info)
      end do

    case('ML_TYPE','AGGR_ALG','AGGR_KIND',&
         & 'SMOOTHER_SWEEPS_PRE','SMOOTHER_SWEEPS_POST',&
         & 'SMOOTHER_POS','AGGR_OMEGA_ALG',&
         & 'AGGR_EIG','AGGR_FILTER')
      do ilev_=1,nlev_
        call p%precv(ilev_)%set(what,val,info)
      end do

    case('COARSE_MAT')
      if (nlev_ > 1) then 
        call p%precv(nlev_)%set('COARSE_MAT',val,info)
      end if

    case('COARSE_SOLVE')
      if (nlev_ > 1) then 

        call p%precv(nlev_)%set('COARSE_SOLVE',val,info)
        select case (val) 
        case(mld_bjac_)
          call onelev_set_smoother(p%precv(nlev_),mld_bjac_,info)
#if defined(HAVE_UMF_)
          call onelev_set_solver(p%precv(nlev_),mld_umf_,info)
#elif defined(HAVE_SLU_) 
          call onelev_set_solver(p%precv(nlev_),mld_slu_,info)
#elif defined(HAVE_MUMPS_) 
          call onelev_set_solver(p%precv(nlev_),mld_mumps_,info)
#else 
          call onelev_set_solver(p%precv(nlev_),mld_ilu_n_,info)
#endif
          call p%precv(nlev_)%set('COARSE_MAT',mld_distr_mat_,info)
        case(mld_umf_, mld_slu_,mld_ilu_n_, mld_ilu_t_,mld_milu_n_)
          call onelev_set_smoother(p%precv(nlev_),mld_bjac_,info)
          call onelev_set_solver(p%precv(nlev_),val,info)
          call p%precv(nlev_)%set('COARSE_MAT',mld_repl_mat_,info)
        case(mld_sludist_)
          call onelev_set_smoother(p%precv(nlev_),mld_bjac_,info)
          call onelev_set_solver(p%precv(nlev_),val,info)
          call p%precv(nlev_)%set('COARSE_MAT',mld_distr_mat_,info)
        case(mld_mumps_)
            call onelev_set_smoother(p%precv(nlev_),mld_bjac_,info)
            call onelev_set_solver(p%precv(nlev_),val,info)
            call p%precv(nlev_)%set('COARSE_MAT',mld_distr_mat_,info)
        case(mld_jac_)
          call onelev_set_smoother(p%precv(nlev_),mld_bjac_,info)
          call onelev_set_solver(p%precv(nlev_),mld_diag_scale_,info)
          call p%precv(nlev_)%set('COARSE_MAT',mld_distr_mat_,info)
        end select

      endif

    case('COARSE_SUBSOLVE')
      if (nlev_ > 1) then 
        call onelev_set_solver(p%precv(nlev_),val,info)
      endif

    case('COARSE_SWEEPS')

      if (nlev_ > 1) then
        call p%precv(nlev_)%set('SMOOTHER_SWEEPS',val,info)
      end if

    case('COARSE_FILLIN')
      if (nlev_ > 1) then 
        call p%precv(nlev_)%set('SUB_FILLIN',val,info)
      end if
    case default
      do ilev_=1,nlev_
        call p%precv(ilev_)%set(what,val,info)
      end do
    end select

  endif

contains

  subroutine onelev_set_smoother(level,val,info)
    type(mld_d_onelev_type), intent(inout) :: level
    integer(psb_ipk_), intent(in)          :: val
    integer(psb_ipk_), intent(out)         :: info
    info = psb_success_

    !
    ! This here requires a bit more attention.
    !
    select case (val) 
    case (mld_noprec_)
      if (allocated(level%sm)) then 
        select type (sm => level%sm)
        type is (mld_d_base_smoother_type) 
          ! do nothing
        class default
          call level%sm%free(info)
          if (info == 0) deallocate(level%sm)
          if (info == 0) allocate(mld_d_base_smoother_type ::&
               & level%sm, stat=info)
          if (info == 0) allocate(mld_d_id_solver_type ::&
               & level%sm%sv, stat=info) 
        end select
      else 
        allocate(mld_d_base_smoother_type ::&
             &  level%sm, stat=info)
        if (info ==0) allocate(mld_d_id_solver_type ::&
             & level%sm%sv, stat=info) 
      endif

    case (mld_jac_)
      if (allocated(level%sm)) then 
        select type (sm => level%sm)
        class is (mld_d_jac_smoother_type) 
            ! do nothing
        class default
          call level%sm%free(info)
          if (info == 0) deallocate(level%sm)
          if (info == 0) allocate(mld_d_jac_smoother_type :: &
               & level%sm, stat=info)
          if (info == 0) allocate(mld_d_diag_solver_type :: &
               & level%sm%sv, stat=info)
        end select
      else 
        allocate(mld_d_jac_smoother_type :: level%sm, stat=info)
        if (info == 0) allocate(mld_d_diag_solver_type ::&
             & level%sm%sv, stat=info)
      endif

    case (mld_bjac_)
      if (allocated(level%sm)) then 
        select type (sm => level%sm)
        class is (mld_d_jac_smoother_type) 
            ! do nothing
        class default
          call level%sm%free(info)
          if (info == 0) deallocate(level%sm)
          if (info == 0) allocate(mld_d_jac_smoother_type ::&
               & level%sm, stat=info)
          if (info == 0) allocate(mld_d_ilu_solver_type ::&
               & level%sm%sv, stat=info)
        end select
      else 
        allocate(mld_d_jac_smoother_type :: level%sm, stat=info)
        if (info == 0) allocate(mld_d_ilu_solver_type ::&
             & level%sm%sv, stat=info)
      endif

    case (mld_as_)
      if (allocated(level%sm)) then 
        select type (sm => level%sm)
        class is (mld_d_as_smoother_type) 
            ! do nothing
        class default
          call level%sm%free(info)
          if (info == 0) deallocate(level%sm)
          if (info == 0) allocate(mld_d_as_smoother_type ::&
               & level%sm, stat=info)
          if (info == 0) allocate(mld_d_ilu_solver_type ::&
               & level%sm%sv, stat=info)
        end select
      else 
        allocate(mld_d_as_smoother_type :: level%sm, stat=info)
        if (info == 0) allocate(mld_d_ilu_solver_type ::&
             & level%sm%sv, stat=info)
      endif

    case default
      !
      ! Do nothing and hope for the best :) 
      !
    end select
    if (allocated(level%sm)) &
         & call level%sm%default()

  end subroutine onelev_set_smoother

  subroutine onelev_set_solver(level,val,info)
    type(mld_d_onelev_type), intent(inout) :: level
    integer(psb_ipk_), intent(in)          :: val
    integer(psb_ipk_), intent(out)         :: info
    info = psb_success_

    !
    ! This here requires a bit more attention.
    !
    select case (val) 
    case (mld_f_none_)
      if (allocated(level%sm)) then 
        if (allocated(level%sm%sv)) then 
          select type (sv => level%sm%sv)
          class is (mld_d_id_solver_type) 
            ! do nothing
          class default
            call level%sm%sv%free(info)
            if (info == 0) deallocate(level%sm%sv)
            if (info == 0) allocate(mld_d_id_solver_type ::&
                 & level%sm%sv, stat=info)
          end select
        else 
          allocate(mld_d_id_solver_type :: level%sm%sv, stat=info)
        endif
        if (allocated(level%sm)) then 
          if (allocated(level%sm%sv)) &
               & call level%sm%sv%default()
        end if
      else
        write(0,*) 'Calling set_solver without a smoother?'
        info = -5
      end if
      
    case (mld_diag_scale_)
      if (allocated(level%sm)) then 
        if (allocated(level%sm%sv)) then 
          select type (sv => level%sm%sv)
          class is (mld_d_diag_solver_type) 
            ! do nothing
          class default
            call level%sm%sv%free(info)
            if (info == 0) deallocate(level%sm%sv)
            if (info == 0) allocate(mld_d_diag_solver_type ::&
                 &  level%sm%sv, stat=info)
          end select
        else 
          allocate(mld_d_diag_solver_type :: level%sm%sv, stat=info)
        endif
        if (allocated(level%sm)) then 
          if (allocated(level%sm%sv)) &
               & call level%sm%sv%default()
        end if
      else
        write(0,*) 'Calling set_solver without a smoother?'
        info = -5
      end if

    case (mld_gs_)
      if (allocated(level%sm)) then 
        if (allocated(level%sm%sv)) then 
          select type (sv => level%sm%sv)
            class is (mld_d_gs_solver_type) 
              ! do nothing
            class default
            call level%sm%sv%free(info)
            if (info == 0) deallocate(level%sm%sv)
            if (info == 0) allocate(mld_d_gs_solver_type ::&
                 &  level%sm%sv, stat=info)
          end select
        else 
          allocate(mld_d_gs_solver_type :: level%sm%sv, stat=info)
        endif
        if (allocated(level%sm%sv)) then 
          call level%sm%sv%default()
        else
        endif

      else
        write(0,*) 'Calling set_solver without a smoother?'
        info = -5
      end if

    case (mld_ilu_n_,mld_milu_n_,mld_ilu_t_)
      if (allocated(level%sm)) then 
        if (allocated(level%sm%sv)) then 
          select type (sv => level%sm%sv)
          class is (mld_d_ilu_solver_type) 
            ! do nothing
          class default
            call level%sm%sv%free(info)
            if (info == 0) deallocate(level%sm%sv)
            if (info == 0) allocate(mld_d_ilu_solver_type ::&
                 & level%sm%sv, stat=info)
          end select
        else 
          allocate(mld_d_ilu_solver_type :: level%sm%sv, stat=info)
        endif
        if (allocated(level%sm)) then 
          if (allocated(level%sm%sv)) &
               & call level%sm%sv%default()
        end if
        call level%sm%sv%set('SUB_SOLVE',val,info)
      else
        write(0,*) 'Calling set_solver without a smoother?'
        info = -5
      end if
      
#ifdef HAVE_SLU_
    case (mld_slu_) 
      if (allocated(level%sm)) then 
        if (allocated(level%sm%sv)) then 
          select type (sv => level%sm%sv)
          class is (mld_d_slu_solver_type) 
            ! do nothing
          class default
            call level%sm%sv%free(info)
            if (info == 0) deallocate(level%sm%sv)
            if (info == 0) allocate(mld_d_slu_solver_type ::&
                 & level%sm%sv, stat=info)
          end select
        else 
          allocate(mld_d_slu_solver_type :: level%sm%sv, stat=info)
        endif
        if (allocated(level%sm)) then 
          if (allocated(level%sm%sv)) &
               & call level%sm%sv%default()
        end if
      else
        write(0,*) 'Calling set_solver without a smoother?'
        info = -5
      end if
#endif
#ifdef HAVE_UMF_
    case (mld_umf_) 
      if (allocated(level%sm)) then 
        if (allocated(level%sm%sv)) then 
          select type (sv => level%sm%sv)
          class is (mld_d_umf_solver_type) 
            ! do nothing
          class default
            call level%sm%sv%free(info)
            if (info == 0) deallocate(level%sm%sv)
            if (info == 0) allocate(mld_d_umf_solver_type ::&
                 & level%sm%sv, stat=info)
          end select
        else 
          allocate(mld_d_umf_solver_type :: level%sm%sv, stat=info)
        endif
        if (allocated(level%sm)) then 
          if (allocated(level%sm%sv)) &
               & call level%sm%sv%default()
        end if
      else
        write(0,*) 'Calling set_solver without a smoother?'
        info = -5
      end if
#endif
#ifdef HAVE_SLUDIST_
    case (mld_sludist_) 
      if (allocated(level%sm)) then 
        if (allocated(level%sm%sv)) then 
          select type (sv => level%sm%sv)
          class is (mld_d_sludist_solver_type) 
            ! do nothing
          class default
            call level%sm%sv%free(info)
            if (info == 0) deallocate(level%sm%sv)
            if (info == 0) allocate(mld_d_sludist_solver_type ::&
                 & level%sm%sv, stat=info)
          end select
        else 
          allocate(mld_d_sludist_solver_type :: level%sm%sv, stat=info)
        endif
        if (allocated(level%sm)) then 
          if (allocated(level%sm%sv)) &
               & call level%sm%sv%default()
        end if
      else
        write(0,*) 'Calling set_solver without a smoother?'
        info = -5
      end if
#endif

#ifdef HAVE_MUMPS_
    case (mld_mumps_) 
      if (allocated(level%sm%sv)) then 
        select type (sv => level%sm%sv)
        class is (mld_d_mumps_solver_type) 
            ! do nothing
        class default
          call level%sm%sv%free(info)
          if (info == 0) deallocate(level%sm%sv)
          if (info == 0) allocate(mld_d_mumps_solver_type ::&
               & level%sm%sv, stat=info)
        end select
      else 
        allocate(mld_d_mumps_solver_type :: level%sm%sv, stat=info)
      endif
      if (allocated(level%sm)) then 
        if (allocated(level%sm%sv)) &
             & call level%sm%sv%default()
      end if
#endif
    case default
      !
      ! Do nothing and hope for the best :) 
      !
    end select

  end subroutine onelev_set_solver


end subroutine mld_dcprecseti

!
! Subroutine: mld_dprecsetc
! Version: real
!
!  This routine sets the character parameters defining the preconditioner. More
!  precisely, the character parameter identified by 'what' is assigned the value
!  contained in 'val'.
!  For the multilevel preconditioners, the levels are numbered in increasing
!  order starting from the finest one, i.e. level 1 is the finest level. 
!
!  To set integer and real parameters, see mld_dprecseti and mld_dprecsetr,
!  respectively.
!
!
! Arguments:
!    p       -  type(mld_dprec_type), input/output.
!               The preconditioner data structure.
!    what    -  integer, input.
!               The number identifying the parameter to be set.
!               A mnemonic constant has been associated to each of these
!               numbers, as reported in the MLD2P4 User's and Reference Guide.
!    string  -  character(len=*), input.
!               The value of the parameter to be set. The list of allowed
!               values is reported in the MLD2P4 User's and Reference Guide.
!    info    -  integer, output.
!               Error code.
!    ilev    -  integer, optional, input.
!               For the multilevel preconditioner, the level at which the
!               preconditioner parameter has to be set. 
!               If nlev is not present, the parameter identified by 'what'
!               is set at all the appropriate levels.
!
!  NOTE: currently, the use of the argument ilev is not "safe" and is reserved to
!  MLD2P4 developers. Indeed, by using ilev it is possible to set different values
!  of the same parameter at different levels 1,...,nlev-1, even in cases where
!  the parameter must have the same value at all the levels but the coarsest one.
!  For this reason, the interface mld_precset to this routine has been built in
!  such a way that ilev is not visible to the user (see mld_prec_mod.f90).
!   
subroutine mld_dcprecsetc(p,what,string,info,ilev)

  use psb_base_mod
  use mld_d_prec_mod, mld_protect_name => mld_dcprecsetc

  implicit none

  ! Arguments
  class(mld_dprec_type), intent(inout)    :: p
  character(len=*), intent(in)            :: what 
  character(len=*), intent(in)            :: string
  integer(psb_ipk_), intent(out)          :: info
  integer(psb_ipk_), optional, intent(in) :: ilev

  ! Local variables
  integer(psb_ipk_)                      :: ilev_, nlev_,val
  character(len=*), parameter            :: name='mld_precsetc'

  info = psb_success_

  if (.not.allocated(p%precv)) then 
    info = 3111
    return 
  endif
  nlev_ = size(p%precv)

  if (present(ilev)) then 
    ilev_ = ilev
  else
    ilev_ = 1 
  end if

  if ((ilev_<1).or.(ilev_ > nlev_)) then 
    write(psb_err_unit,*) name,&
         & ': Error: invalid ILEV/NLEV combination',ilev_, nlev_
    info = -1
    return
  endif

  val =  mld_stringval(string)
  if (val >=0)  then 
    call p%set(what,val,info,ilev=ilev)
  else
    call p%precv(ilev_)%set(what,string,info)
  end if

end subroutine mld_dcprecsetc


!
! Subroutine: mld_dprecsetr
! Version: real
!
!  This routine sets the real parameters defining the preconditioner. More
!  precisely, the real parameter identified by 'what' is assigned the value
!  contained in 'val'.
!  For the multilevel preconditioners, the levels are numbered in increasing
!  order starting from the finest one, i.e. level 1 is the finest level. 
!
!  To set integer and character parameters, see mld_dprecseti and mld_dprecsetc,
!  respectively.
!
! Arguments:
!    p       -  type(mld_dprec_type), input/output.
!               The preconditioner data structure.
!    what    -  integer, input.
!               The number identifying the parameter to be set.
!               A mnemonic constant has been associated to each of these
!               numbers, as reported in the MLD2P4 User's and Reference Guide.
!    val     -  real(psb_dpk_), input.
!               The value of the parameter to be set. The list of allowed
!               values is reported in the MLD2P4 User's and Reference Guide.
!    info    -  integer, output.
!               Error code.
!    ilev    -  integer, optional, input.
!               For the multilevel preconditioner, the level at which the
!               preconditioner parameter has to be set. 
!               If nlev is not present, the parameter identified by 'what'
!               is set at all the appropriate levels.
!
!  NOTE: currently, the use of the argument ilev is not "safe" and is reserved to
!  MLD2P4 developers. Indeed, by using ilev it is possible to set different values
!  of the same parameter at different levels 1,...,nlev-1, even in cases where
!  the parameter must have the same value at all the levels but the coarsest one.
!  For this reason, the interface mld_precset to this routine has been built in
!  such a way that ilev is not visible to the user (see mld_prec_mod.f90).
!   
subroutine mld_dcprecsetr(p,what,val,info,ilev)

  use psb_base_mod
  use mld_d_prec_mod, mld_protect_name => mld_dcprecsetr

  implicit none

  ! Arguments
  class(mld_dprec_type), intent(inout)    :: p
  character(len=*), intent(in)            :: what 
  real(psb_dpk_), intent(in)              :: val
  integer(psb_ipk_), intent(out)          :: info
  integer(psb_ipk_), optional, intent(in) :: ilev

! Local variables
  integer(psb_ipk_)                      :: ilev_,nlev_
  real(psb_dpk_)                         :: thr 
  character(len=*), parameter            :: name='mld_precsetr'

  info = psb_success_

  if (present(ilev)) then 
    ilev_ = ilev
  else
    ilev_ = 1 
  end if

  if (.not.allocated(p%precv)) then 
    write(psb_err_unit,*) name,&
         &': Error: uninitialized preconditioner,',&
         &' should call MLD_PRECINIT' 
    info = 3111
    return 
  endif
  nlev_ = size(p%precv)

  if ((ilev_<1).or.(ilev_ > nlev_)) then 
    write(psb_err_unit,*) name,&
         & ': Error: invalid ILEV/NLEV combination',&
         & ilev_, nlev_
    info = -1
    return
  endif

  !
  ! Set preconditioner parameters at level ilev.
  !
  if (present(ilev)) then 
    
    call p%precv(ilev_)%set(what,val,info)

  else if (.not.present(ilev)) then 
      !
      ! ilev not specified: set preconditioner parameters at all the appropriate levels
      !

      select case(psb_toupper(what)) 
      case('COARSE_ILUTHRS')
        ilev_=nlev_
        call p%precv(ilev_)%set('SUB_ILUTHRS',val,info)

      case('AGGR_THRESH')
        thr = val
        do ilev_ = 2, nlev_
          call p%precv(ilev_)%set('AGGR_THRESH',thr,info)
          thr = thr * p%precv(ilev_)%parms%aggr_scale
        end do

      case default

        do ilev_=1,nlev_
          call p%precv(ilev_)%set(what,val,info)
        end do
      end select

  endif

end subroutine mld_dcprecsetr


