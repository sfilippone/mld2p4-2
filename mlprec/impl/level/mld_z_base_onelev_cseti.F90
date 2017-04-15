!  
!   
!                             MLD2P4  version 2.1
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.4)
!    
!    (C) Copyright 2008, 2010, 2012, 2015, 2017 
!  
!        Salvatore Filippone    Cranfield University
!        Ambra Abdullahi Hassan University of Rome Tor Vergata
!        Alfredo Buttari        CNRS-IRIT, Toulouse
!        Pasqua D'Ambra         ICAR-CNR, Naples
!        Daniela di Serafino    University of Campania "L. Vanvitelli", Caserta
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
subroutine mld_z_base_onelev_cseti(lv,what,val,info,pos)
  
  use psb_base_mod
  use mld_z_onelev_mod, mld_protect_name => mld_z_base_onelev_cseti
  use mld_z_jac_smoother
  use mld_z_as_smoother
  use mld_z_diag_solver
  use mld_z_ilu_solver
  use mld_z_id_solver
  use mld_z_gs_solver
#if defined(HAVE_UMF_)
  use mld_z_umf_solver
#endif
#if defined(HAVE_SLUDIST_)
  use mld_z_sludist_solver
#endif
#if defined(HAVE_SLU_)
  use mld_z_slu_solver
#endif
#if defined(HAVE_MUMPS_)
  use mld_z_mumps_solver
#endif

  Implicit None

  ! Arguments
  class(mld_z_onelev_type), intent(inout) :: lv 
  character(len=*), intent(in)              :: what 
  integer(psb_ipk_), intent(in)             :: val
  integer(psb_ipk_), intent(out)            :: info
  character(len=*), optional, intent(in)      :: pos
  ! Local 
  integer(psb_ipk_)  :: ipos_, err_act
  character(len=20) :: name='z_base_onelev_cseti'
  type(mld_z_base_smoother_type) :: mld_z_base_smoother_mold
  type(mld_z_jac_smoother_type)  ::  mld_z_jac_smoother_mold
  type(mld_z_as_smoother_type)   ::  mld_z_as_smoother_mold
  type(mld_z_diag_solver_type)   ::  mld_z_diag_solver_mold
  type(mld_z_ilu_solver_type)    ::  mld_z_ilu_solver_mold
  type(mld_z_id_solver_type)     ::  mld_z_id_solver_mold
  type(mld_z_gs_solver_type)     ::  mld_z_gs_solver_mold
  type(mld_z_bwgs_solver_type)   ::  mld_z_bwgs_solver_mold
#if defined(HAVE_UMF_)
  type(mld_z_umf_solver_type)    ::  mld_z_umf_solver_mold
#endif
#if defined(HAVE_SLUDIST_)
  type(mld_z_sludist_solver_type) ::  mld_z_sludist_solver_mold
#endif
#if defined(HAVE_SLU_)
  type(mld_z_slu_solver_type)   ::  mld_z_slu_solver_mold
#endif
#if defined(HAVE_MUMPS_)
  type(mld_z_mumps_solver_type) ::  mld_z_mumps_solver_mold
#endif
  
  call psb_erractionsave(err_act)
  info = psb_success_

  if (present(pos)) then
    select case(psb_toupper(trim(pos)))
    case('PRE')
      ipos_ = mld_pre_smooth_
    case('POST')
      ipos_ = mld_post_smooth_
    case default
      ipos_ = mld_pre_smooth_
    end select
  else
    ipos_ = mld_pre_smooth_
  end if
  
  select case (psb_toupper(what))
  case ('SMOOTHER_TYPE')
    select case (val) 
    case (mld_noprec_)
      call lv%set(mld_z_base_smoother_mold,info,pos=pos)
      if (info == 0) call lv%set(mld_z_id_solver_mold,info,pos=pos)
      
    case (mld_jac_)
      call lv%set(mld_z_jac_smoother_mold,info,pos=pos)
      if (info == 0) call lv%set(mld_z_diag_solver_mold,info,pos=pos)
      
    case (mld_bjac_)
      call lv%set(mld_z_jac_smoother_mold,info,pos=pos)
      if (info == 0) call lv%set(mld_z_ilu_solver_mold,info,pos=pos)

    case (mld_as_)
      call lv%set(mld_z_as_smoother_mold,info,pos=pos)
      if (info == 0) call lv%set(mld_z_ilu_solver_mold,info,pos=pos)

    case (mld_fbgs_)
      call lv%set(mld_z_jac_smoother_mold,info,pos='pre')
      if (info == 0) call lv%set(mld_z_gs_solver_mold,info,pos='pre')
      call lv%set(mld_z_jac_smoother_mold,info,pos='post')
      if (info == 0) call lv%set(mld_z_bwgs_solver_mold,info,pos='post')
      
    case default
      !
      ! Do nothing and hope for the best :) 
      !
    end select
    if (ipos_==mld_pre_smooth_) then 
      if (allocated(lv%sm)) call lv%sm%default()
    else if (ipos_==mld_post_smooth_) then
      if (allocated(lv%sm2a)) call lv%sm2a%default()
    end if
    

  case('SUB_SOLVE')
    select case (val) 
    case (mld_f_none_)
      call lv%set(mld_z_id_solver_mold,info,pos=pos)
      
    case (mld_diag_scale_)
      call lv%set(mld_z_diag_solver_mold,info,pos=pos)
      
    case (mld_gs_)
      call lv%set(mld_z_gs_solver_mold,info,pos=pos)
      
    case (mld_bwgs_)
      call lv%set(mld_z_bwgs_solver_mold,info,pos=pos)
      
    case (mld_ilu_n_,mld_milu_n_,mld_ilu_t_)
      call lv%set(mld_z_ilu_solver_mold,info,pos=pos)
      if (info == 0) then
        select case(ipos_)
        case(mld_pre_smooth_) 
          call lv%sm%sv%set('SUB_SOLVE',val,info)
        case (mld_post_smooth_)
          if (allocated(lv%sm2a)) call lv%sm2a%sv%set('SUB_SOLVE',val,info)
        case default
          ! Impossible!! 
          info = psb_err_internal_error_
        end select
      end if
#ifdef HAVE_SLU_
    case (mld_slu_) 
      call lv%set(mld_z_slu_solver_mold,info,pos=pos)
#endif
#ifdef HAVE_MUMPS_
    case (mld_mumps_) 
      call lv%set(mld_z_mumps_solver_mold,info,pos=pos)
#endif
#ifdef HAVE_SLUDIST_
    case (mld_sludist_)
      call lv%set(mld_z_sludist_solver_mold,info,pos=pos)
#endif
#ifdef HAVE_UMF_
    case (mld_umf_)
      call lv%set(mld_z_umf_solver_mold,info,pos=pos)
#endif
    case default
      !
      ! Do nothing and hope for the best :) 
      !
    end select
    

  case ('SMOOTHER_SWEEPS')
    lv%parms%sweeps      = val
    lv%parms%sweeps_pre  = val
    lv%parms%sweeps_post = val

  case ('SMOOTHER_SWEEPS_PRE')
    lv%parms%sweeps_pre  = val

  case ('SMOOTHER_SWEEPS_POST')
    lv%parms%sweeps_post = val

  case ('ML_CYCLE')
    lv%parms%ml_cycle      = val

  case ('PAR_AGGR_ALG')
    lv%parms%par_aggr_alg  = val

  case ('AGGR_ORD')
    lv%parms%aggr_ord      = val

  case ('AGGR_TYPE')
    lv%parms%aggr_type     = val

  case ('AGGR_PROL')
    lv%parms%aggr_prol     = val

  case ('COARSE_MAT')
    lv%parms%coarse_mat    = val

  case ('AGGR_OMEGA_ALG')
    lv%parms%aggr_omega_alg= val

  case ('AGGR_EIG')
    lv%parms%aggr_eig      = val

  case ('AGGR_FILTER')
    lv%parms%aggr_filter   = val

  case ('COARSE_SOLVE')
    lv%parms%coarse_solve    = val

  case default
    select case(ipos_)
    case(mld_pre_smooth_) 
      if (allocated(lv%sm)) then 
        call lv%sm%set(what,val,info)
      end if
    case (mld_post_smooth_)
      if (allocated(lv%sm2a)) then 
        call lv%sm2a%set(what,val,info)
      end if
    case default
      ! Impossible!! 
      info = psb_err_internal_error_
    end select

  end select
  if (info /= psb_success_) goto 9999
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine mld_z_base_onelev_cseti
