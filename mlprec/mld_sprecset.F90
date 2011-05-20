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
! File: mld_sprecset.f90
!
! Subroutine: mld_sprecseti
! Version: real
!
!  This routine sets the integer parameters defining the preconditioner. More
!  precisely, the integer parameter identified by 'what' is assigned the value
!  contained in 'val'.
!  For the multilevel preconditioners, the levels are numbered in increasing
!  order starting from the finest one, i.e. level 1 is the finest level. 
!
!  To set character and real parameters, see mld_sprecsetc and mld_sprecsetr,
!  respectively.
!
!
! Arguments:
!    p       -  type(mld_sprec_type), input/output.
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
subroutine mld_sprecseti(p,what,val,info,ilev)

  use psb_base_mod
  use mld_s_prec_mod, mld_protect_name => mld_sprecseti
  use mld_s_jac_smoother
  use mld_s_as_smoother
  use mld_s_diag_solver
  use mld_s_ilu_solver
  use mld_s_id_solver
#ifdef HAVE_SLU_
  use mld_s_slu_solver
#endif

  implicit none

  ! Arguments
  type(mld_sprec_type), intent(inout)    :: p
  integer, intent(in)                    :: what 
  integer, intent(in)                    :: val
  integer, intent(out)                   :: info
  integer, optional, intent(in)          :: ilev

  ! Local variables
  integer                                :: ilev_, nlev_
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

  !
  ! Set preconditioner parameters at level ilev.
  !
  if (present(ilev)) then 

    if (ilev_ == 1) then
      ! 
      ! Rules for fine level are slightly different.
      ! 
      select case(what) 
      case(mld_smoother_type_)
        call onelev_set_smoother(p%precv(ilev_),val,info)
      case(mld_sub_solve_)
        call onelev_set_solver(p%precv(ilev_),val,info)
      case(mld_smoother_sweeps_,mld_ml_type_,mld_aggr_alg_,mld_aggr_kind_,&
           & mld_smoother_pos_,mld_aggr_omega_alg_,mld_aggr_eig_,&
           & mld_smoother_sweeps_pre_,mld_smoother_sweeps_post_,&
           & mld_sub_restr_,mld_sub_prol_, &
           & mld_sub_ren_,mld_sub_ovr_,mld_sub_fillin_)
        call p%precv(ilev_)%set(what,val,info)

      case default
        call p%precv(ilev_)%set(what,val,info)
      end select

    else if (ilev_ > 1) then 

      select case(what) 
      case(mld_smoother_type_)
        call onelev_set_smoother(p%precv(ilev_),val,info)
      case(mld_sub_solve_)
        call onelev_set_solver(p%precv(ilev_),val,info)
      case(mld_smoother_sweeps_,mld_ml_type_,mld_aggr_alg_,mld_aggr_kind_,&
           & mld_smoother_pos_,mld_aggr_omega_alg_,mld_aggr_eig_,&
           & mld_smoother_sweeps_pre_,mld_smoother_sweeps_post_,&
           & mld_sub_restr_,mld_sub_prol_, &
           & mld_sub_ren_,mld_sub_ovr_,mld_sub_fillin_,&
           & mld_coarse_mat_)
        call p%precv(ilev_)%set(what,val,info)

      case(mld_coarse_subsolve_)
        if (ilev_ /= nlev_) then 
          write(psb_err_unit,*) name,&
               & ': Error: Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if
        call onelev_set_solver(p%precv(ilev_),val,info)
      case(mld_coarse_solve_)
        if (ilev_ /= nlev_) then 
          write(psb_err_unit,*) name,&
               & ': Error: Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if

        if (nlev_ > 1) then 
          call p%precv(nlev_)%set(mld_coarse_solve_,val,info)
          select case (val) 
          case(mld_bjac_)
            call onelev_set_smoother(p%precv(nlev_),val,info)
#if defined(HAVE_SLU_) 
            call onelev_set_solver(p%precv(nlev_),mld_slu_,info)
#else 
            call onelev_set_solver(p%precv(nlev_),mld_ilu_n_,info)
#endif
            call p%precv(nlev_)%set(mld_coarse_mat_,mld_distr_mat_,info)
          case(mld_umf_, mld_slu_,mld_ilu_n_, mld_ilu_t_,mld_milu_n_)
            call onelev_set_smoother(p%precv(nlev_),mld_bjac_,info)
            call onelev_set_solver(p%precv(nlev_),val,info)
            call p%precv(nlev_)%set(mld_coarse_mat_,mld_repl_mat_,info)
          case(mld_sludist_)
            call onelev_set_smoother(p%precv(nlev_),mld_bjac_,info)
            call onelev_set_solver(p%precv(nlev_),val,info)
            call p%precv(nlev_)%set(mld_coarse_mat_,mld_distr_mat_,info)
          case(mld_jac_)
            call onelev_set_smoother(p%precv(nlev_),mld_bjac_,info)
            call onelev_set_solver(p%precv(nlev_),mld_diag_scale_,info)
            call p%precv(nlev_)%set(mld_coarse_mat_,mld_distr_mat_,info)
          end select

        endif
      case(mld_coarse_sweeps_)
        if (ilev_ /= nlev_) then 
          write(psb_err_unit,*) name,&
               & ': Error: Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if
        call p%precv(nlev_)%set(mld_smoother_sweeps_,val,info)

      case(mld_coarse_fillin_)
        if (ilev_ /= nlev_) then 
          write(psb_err_unit,*) name,&
               & ': Error: Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if
        call p%precv(nlev_)%set(mld_sub_fillin_,val,info)
      case default
        call p%precv(ilev_)%set(what,val,info)
      end select

    endif

  else if (.not.present(ilev)) then 
    !
    ! ilev not specified: set preconditioner parameters at all the appropriate
    ! levels
    !
    select case(what) 
    case(mld_sub_solve_)
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

    case(mld_sub_restr_,mld_sub_prol_,&
         & mld_sub_ren_,mld_sub_ovr_,mld_sub_fillin_)
      do ilev_=1,max(1,nlev_-1)
        call p%precv(ilev_)%set(what,val,info)
      end do

    case(mld_smoother_sweeps_)
      do ilev_=1,max(1,nlev_-1)
        call p%precv(ilev_)%set(what,val,info)
      end do

    case(mld_smoother_type_)
      do ilev_=1,max(1,nlev_-1)
        call onelev_set_smoother(p%precv(ilev_),val,info)
      end do

    case(mld_ml_type_,mld_aggr_alg_,mld_aggr_kind_,&
         & mld_smoother_sweeps_pre_,mld_smoother_sweeps_post_,&
         & mld_smoother_pos_,mld_aggr_omega_alg_,&
         & mld_aggr_eig_,mld_aggr_filter_)
      do ilev_=1,nlev_
        call p%precv(ilev_)%set(what,val,info)
      end do

    case(mld_coarse_mat_)
      if (nlev_ > 1) then 
        call p%precv(nlev_)%set(mld_coarse_mat_,val,info)
      end if

    case(mld_coarse_solve_)
      if (nlev_ > 1) then 

        call p%precv(nlev_)%set(mld_coarse_solve_,val,info)
        select case (val) 
        case(mld_bjac_)
          call onelev_set_smoother(p%precv(nlev_),mld_bjac_,info)
#if defined(HAVE_SLU_) 
          call onelev_set_solver(p%precv(nlev_),mld_slu_,info)
#else 
          call onelev_set_solver(p%precv(nlev_),mld_ilu_n_,info)
#endif
          call p%precv(nlev_)%set(mld_coarse_mat_,mld_distr_mat_,info)
        case(mld_umf_, mld_slu_,mld_ilu_n_, mld_ilu_t_,mld_milu_n_)
          call onelev_set_smoother(p%precv(nlev_),mld_bjac_,info)
          call onelev_set_solver(p%precv(nlev_),val,info)
          call p%precv(nlev_)%set(mld_coarse_mat_,mld_repl_mat_,info)
        case(mld_sludist_)
          call onelev_set_smoother(p%precv(nlev_),mld_bjac_,info)
          call onelev_set_solver(p%precv(nlev_),val,info)
          call p%precv(nlev_)%set(mld_coarse_mat_,mld_distr_mat_,info)
        case(mld_jac_)
          call onelev_set_smoother(p%precv(nlev_),mld_bjac_,info)
          call onelev_set_solver(p%precv(nlev_),mld_diag_scale_,info)
          call p%precv(nlev_)%set(mld_coarse_mat_,mld_distr_mat_,info)
        end select

      endif

    case(mld_coarse_subsolve_)
      if (nlev_ > 1) then 
        call onelev_set_solver(p%precv(nlev_),val,info)
      endif

    case(mld_coarse_sweeps_)

      if (nlev_ > 1) then
        call p%precv(nlev_)%set(mld_smoother_sweeps_,val,info)
      end if

    case(mld_coarse_fillin_)
      if (nlev_ > 1) then 
        call p%precv(nlev_)%set(mld_sub_fillin_,val,info)
      end if
    case default
      do ilev_=1,nlev_
        call p%precv(ilev_)%set(what,val,info)
      end do
    end select

  endif

contains

  subroutine onelev_set_smoother(level,val,info)
    type(mld_sonelev_type), intent(inout) :: level
    integer, intent(in)                   :: val
    integer, intent(out)                  :: info
    info = psb_success_

    !
    ! This here requires a bit more attention.
    !
    select case (val) 
    case (mld_noprec_)
      if (allocated(level%sm)) then 
        select type (sm => level%sm)
        type is (mld_s_base_smoother_type) 
          ! do nothing
        class default
          call level%sm%free(info)
          if (info == 0) deallocate(level%sm)
          if (info == 0) allocate(mld_s_base_smoother_type ::&
               & level%sm, stat=info)
          if (info == 0) allocate(mld_s_id_solver_type ::&
               & level%sm%sv, stat=info) 
        end select
      else 
        allocate(mld_s_base_smoother_type ::&
             &  level%sm, stat=info)
        if (info ==0) allocate(mld_s_id_solver_type ::&
             & level%sm%sv, stat=info) 
        call level%sm%default()
      endif

    case (mld_jac_)
      if (allocated(level%sm)) then 
        select type (sm => level%sm)
        class is (mld_s_jac_smoother_type) 
          ! do nothing
        class default
          call level%sm%free(info)
          if (info == 0) deallocate(level%sm)
          if (info == 0) allocate(mld_s_jac_smoother_type :: &
               & level%sm, stat=info)
          if (info == 0) allocate(mld_s_diag_solver_type :: &
               & level%sm%sv, stat=info)
        end select
      else 
        allocate(mld_s_jac_smoother_type :: level%sm, stat=info)
        if (info == 0) allocate(mld_s_diag_solver_type ::&
             & level%sm%sv, stat=info)
      endif

    case (mld_bjac_)
      if (allocated(level%sm)) then 
        select type (sm => level%sm)
        class is (mld_s_jac_smoother_type) 
          ! do nothing
        class default
          call level%sm%free(info)
          if (info == 0) deallocate(level%sm)
          if (info == 0) allocate(mld_s_jac_smoother_type ::&
               & level%sm, stat=info)
          if (info == 0) allocate(mld_s_ilu_solver_type ::&
               & level%sm%sv, stat=info)
        end select
      else 
        allocate(mld_s_jac_smoother_type :: level%sm, stat=info)
        if (info == 0) allocate(mld_s_ilu_solver_type ::&
             & level%sm%sv, stat=info)
      endif

    case (mld_as_)
      if (allocated(level%sm)) then 
        select type (sm => level%sm)
        class is (mld_s_as_smoother_type) 
          ! do nothing
        class default
          call level%sm%free(info)
          if (info == 0) deallocate(level%sm)
          if (info == 0) allocate(mld_s_as_smoother_type ::&
               & level%sm, stat=info)
          if (info == 0) allocate(mld_s_ilu_solver_type ::&
               & level%sm%sv, stat=info)
        end select
      else 
        allocate(mld_s_as_smoother_type :: level%sm, stat=info)
        if (info == 0) allocate(mld_s_ilu_solver_type ::&
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
    type(mld_sonelev_type), intent(inout) :: level
    integer, intent(in)                   :: val
    integer, intent(out)                  :: info
    info = psb_success_

    !
    ! This here requires a bit more attention.
    !
    select case (val) 
    case (mld_f_none_)
      if (allocated(level%sm%sv)) then 
        select type (sv => level%sm%sv)
        class is (mld_s_id_solver_type) 
          ! do nothing
        class default
          call level%sm%sv%free(info)
          if (info == 0) deallocate(level%sm%sv)
          if (info == 0) allocate(mld_s_id_solver_type ::&
               & level%sm%sv, stat=info)
        end select
      else 
        allocate(mld_s_id_solver_type :: level%sm%sv, stat=info)
      endif
      if (allocated(level%sm)) then 
        if (allocated(level%sm%sv)) &
             & call level%sm%sv%default()
      end if
      

    case (mld_diag_scale_)
      if (allocated(level%sm%sv)) then 
        select type (sv => level%sm%sv)
        class is (mld_s_diag_solver_type) 
          ! do nothing
        class default
          call level%sm%sv%free(info)
          if (info == 0) deallocate(level%sm%sv)
          if (info == 0) allocate(mld_s_diag_solver_type ::&
               &  level%sm%sv, stat=info)
        end select
      else 
        allocate(mld_s_diag_solver_type :: level%sm%sv, stat=info)
      endif
      if (allocated(level%sm)) then 
        if (allocated(level%sm%sv)) &
             & call level%sm%sv%default()
      end if
    

    case (mld_ilu_n_,mld_milu_n_,mld_ilu_t_)
      if (allocated(level%sm%sv)) then 
        select type (sv => level%sm%sv)
        class is (mld_s_ilu_solver_type) 
          ! do nothing
        class default
          call level%sm%sv%free(info)
          if (info == 0) deallocate(level%sm%sv)
          if (info == 0) allocate(mld_s_ilu_solver_type ::&
               & level%sm%sv, stat=info)
        end select
      else 
        allocate(mld_s_ilu_solver_type :: level%sm%sv, stat=info)
      endif
      if (allocated(level%sm)) then 
        if (allocated(level%sm%sv)) &
             & call level%sm%sv%default()
      end if
      call level%sm%sv%set(mld_sub_solve_,val,info)

#ifdef HAVE_SLU_
    case (mld_slu_) 
      if (allocated(level%sm%sv)) then 
        select type (sv => level%sm%sv)
        class is (mld_s_slu_solver_type) 
          ! do nothing
        class default
          call level%sm%sv%free(info)
          if (info == 0) deallocate(level%sm%sv)
          if (info == 0) allocate(mld_s_slu_solver_type ::&
               & level%sm%sv, stat=info)
        end select
      else 
        allocate(mld_s_slu_solver_type :: level%sm%sv, stat=info)
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


end subroutine mld_sprecseti

subroutine mld_sprecsetsm(p,val,info,ilev)

  use psb_base_mod
  use mld_s_prec_mod, mld_protect_name => mld_sprecsetsm

  implicit none

  ! Arguments
  type(mld_sprec_type), intent(inout)    :: p
  class(mld_s_base_smoother_type), intent(in) :: val
  integer, intent(out)                   :: info
  integer, optional, intent(in)          :: ilev

  ! Local variables
  integer                                :: ilev_, nlev_, ilmin, ilmax
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
    ilmin = ilev
    ilmax = ilev
  else
    ilev_ = 1 
    ilmin = 1
    ilmax = nlev_
  end if

  if ((ilev_<1).or.(ilev_ > nlev_)) then 
    info = -1
    write(psb_err_unit,*) name,&
         & ': Error: invalid ILEV/NLEV combination',ilev_, nlev_
    return
  endif
  

  do ilev_ = ilmin, ilmax 
    if (allocated(p%precv(ilev_)%sm)) then 
      if (allocated(p%precv(ilev_)%sm%sv)) then 
        deallocate(p%precv(ilev_)%sm%sv)
      endif
      deallocate(p%precv(ilev_)%sm)
    end if
#ifdef HAVE_MOLD 
    allocate(p%precv(ilev_)%sm,mold=val) 
#else
    allocate(p%precv(ilev_)%sm,source=val) 
#endif
    call p%precv(ilev_)%sm%default()
  end do

end subroutine mld_sprecsetsm

subroutine mld_sprecsetsv(p,val,info,ilev)

  use psb_base_mod
  use mld_s_prec_mod, mld_protect_name => mld_sprecsetsv

  implicit none

  ! Arguments
  type(mld_sprec_type), intent(inout)    :: p
  class(mld_s_base_solver_type), intent(in) :: val
  integer, intent(out)                   :: info
  integer, optional, intent(in)          :: ilev

  ! Local variables
  integer                                :: ilev_, nlev_, ilmin, ilmax
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
    ilmin = ilev
    ilmax = ilev
  else
    ilev_ = 1 
    ilmin = 1
    ilmax = nlev_
  end if


  if ((ilev_<1).or.(ilev_ > nlev_)) then 
    info = -1
    write(psb_err_unit,*) name,&
         & ': Error: invalid ILEV/NLEV combination',ilev_, nlev_
    return
  endif


  do ilev_ = ilmin, ilmax 
    if (allocated(p%precv(ilev_)%sm)) then 
      if (allocated(p%precv(ilev_)%sm%sv)) &
           & deallocate(p%precv(ilev_)%sm%sv)
#ifdef HAVE_MOLD 
      allocate(p%precv(ilev_)%sm%sv,mold=val) 
#else
      allocate(p%precv(ilev_)%sm%sv,source=val) 
#endif
      call p%precv(ilev_)%sm%sv%default()
    else
      info = 3111
      write(psb_err_unit,*) name,&
           &': Error: uninitialized preconditioner component,',&
           &' should call MLD_PRECINIT/MLD_PRECSET' 
      return 

    end if

  end do



end subroutine mld_sprecsetsv

!
! Subroutine: mld_sprecsetc
! Version: real
!
!  This routine sets the character parameters defining the preconditioner. More
!  precisely, the character parameter identified by 'what' is assigned the value
!  contained in 'val'.
!  For the multilevel preconditioners, the levels are numbered in increasing
!  order starting from the finest one, i.e. level 1 is the finest level. 
!
!  To set integer and real parameters, see mld_sprecseti and mld_sprecsetr,
!  respectively.
!
!
! Arguments:
!    p       -  type(mld_sprec_type), input/output.
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
!   
subroutine mld_sprecsetc(p,what,string,info,ilev)

  use psb_base_mod
  use mld_s_prec_mod, mld_protect_name => mld_sprecsetc

  implicit none

  ! Arguments
  type(mld_sprec_type), intent(inout)    :: p
  integer, intent(in)                    :: what 
  character(len=*), intent(in)           :: string
  integer, intent(out)                   :: info
  integer, optional, intent(in)          :: ilev

  ! Local variables
  integer                                :: ilev_, nlev_,val
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

  call mld_stringval(string,val,info)
  if (info == psb_success_) call mld_inner_precset(p,what,val,info,ilev=ilev)


end subroutine mld_sprecsetc


!
! Subroutine: mld_sprecsetr
! Version: real
!
!  This routine sets the real parameters defining the preconditioner. More
!  precisely, the real parameter identified by 'what' is assigned the value
!  contained in 'val'.
!  For the multilevel preconditioners, the levels are numbered in increasing
!  order starting from the finest one, i.e. level 1 is the finest level. 
!
!  To set integer and character parameters, see mld_sprecseti and mld_sprecsetc,
!  respectively.
!
! Arguments:
!    p       -  type(mld_sprec_type), input/output.
!               The preconditioner data structure.
!    what    -  integer, input.
!               The number identifying the parameter to be set.
!               A mnemonic constant has been associated to each of these
!               numbers, as reported in the MLD2P4 User's and Reference Guide.
!    val     -  real(psb_spk_), input.
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
!   
subroutine mld_sprecsetr(p,what,val,info,ilev)

  use psb_base_mod
  use mld_s_prec_mod, mld_protect_name => mld_sprecsetr

  implicit none

  ! Arguments
  type(mld_sprec_type), intent(inout)    :: p
  integer, intent(in)                    :: what 
  real(psb_spk_), intent(in)           :: val
  integer, intent(out)                   :: info
  integer, optional, intent(in)          :: ilev

! Local variables
  integer                                :: ilev_,nlev_
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

      select case(what) 
      case(mld_coarse_iluthrs_)
        ilev_=nlev_
        call p%precv(ilev_)%set(mld_sub_iluthrs_,val,info)

      case default

        do ilev_=1,nlev_
          call p%precv(ilev_)%set(what,val,info)
        end do
      end select

  endif

end subroutine mld_sprecsetr
