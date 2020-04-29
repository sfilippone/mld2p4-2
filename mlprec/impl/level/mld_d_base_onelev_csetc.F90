!  
!   
!                             MLD2P4  version 2.2
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!    
!    (C) Copyright 2008-2018 
!  
!        Salvatore Filippone  
!        Pasqua D'Ambra   
!        Daniela di Serafino   
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
subroutine mld_d_base_onelev_csetc(lv,what,val,info,pos,idx)
  
  use psb_base_mod
  use mld_d_onelev_mod, mld_protect_name => mld_d_base_onelev_csetc
  use mld_d_base_aggregator_mod
  use mld_d_dec_aggregator_mod
  use mld_d_symdec_aggregator_mod
  use mld_d_jac_smoother
  use mld_d_as_smoother
  use mld_d_diag_solver
  use mld_d_l1_diag_solver
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

  Implicit None

  ! Arguments
  class(mld_d_onelev_type), intent(inout) :: lv 
  character(len=*), intent(in)              :: what 
  character(len=*), intent(in)              :: val
  integer(psb_ipk_), intent(out)            :: info
  character(len=*), optional, intent(in)    :: pos
  integer(psb_ipk_), intent(in), optional   :: idx
  ! Local 
  integer(psb_ipk_)  :: ipos_, err_act
  character(len=20) :: name='d_base_onelev_csetc'
  integer(psb_ipk_) :: ival 
  type(mld_d_base_smoother_type)  :: mld_d_base_smoother_mold
  type(mld_d_jac_smoother_type)   ::  mld_d_jac_smoother_mold
  type(mld_d_as_smoother_type)    ::  mld_d_as_smoother_mold
  type(mld_d_diag_solver_type)    ::  mld_d_diag_solver_mold
  type(mld_d_l1_diag_solver_type) ::  mld_d_l1_diag_solver_mold
  type(mld_d_ilu_solver_type)     ::  mld_d_ilu_solver_mold
  type(mld_d_id_solver_type)      ::  mld_d_id_solver_mold
  type(mld_d_gs_solver_type)      ::  mld_d_gs_solver_mold
  type(mld_d_bwgs_solver_type)    ::  mld_d_bwgs_solver_mold
#if defined(HAVE_UMF_)
  type(mld_d_umf_solver_type)     ::  mld_d_umf_solver_mold
#endif
#if defined(HAVE_SLUDIST_)
  type(mld_d_sludist_solver_type) ::  mld_d_sludist_solver_mold
#endif
#if defined(HAVE_SLU_)
  type(mld_d_slu_solver_type)   ::  mld_d_slu_solver_mold
#endif
#if defined(HAVE_MUMPS_)
  type(mld_d_mumps_solver_type) ::  mld_d_mumps_solver_mold
#endif
  

  call psb_erractionsave(err_act)

  info = psb_success_

  ival = lv%stringval(val)


  if (present(pos)) then
    select case(psb_toupper(trim(pos)))
    case('PRE')
      ipos_ = mld_smooth_pre_
    case('POST')
      ipos_ = mld_smooth_post_
    case default
      ipos_ = mld_smooth_both_
    end select
  else
    ipos_ = mld_smooth_both_
  end if
  
  select case (psb_toupper(trim(what)))
  case ('SMOOTHER_TYPE')
    select case (psb_toupper(trim(val)))
    case ('NOPREC','NONE')
      call lv%set(mld_d_base_smoother_mold,info,pos=pos)
      if (info == 0) call lv%set(mld_d_id_solver_mold,info,pos=pos)
      
    case ('JAC','JACOBI')
      call lv%set(mld_d_jac_smoother_mold,info,pos=pos)
      if (info == 0) call lv%set(mld_d_diag_solver_mold,info,pos=pos)

    case ('L1-JACOBI')
      call lv%set(mld_d_jac_smoother_mold,info,pos=pos)
      if (info == 0) call lv%set(mld_d_l1_diag_solver_mold,info,pos=pos)
      
    case ('BJAC')
      call lv%set(mld_d_jac_smoother_mold,info,pos=pos)
      if (info == 0) call lv%set(mld_d_ilu_solver_mold,info,pos=pos)

    case ('AS')
      call lv%set(mld_d_as_smoother_mold,info,pos=pos)
      if (info == 0) call lv%set(mld_d_ilu_solver_mold,info,pos=pos)

    case ('FBGS')
      call lv%set(mld_d_jac_smoother_mold,info,pos='pre')
      if (info == 0) call lv%set(mld_d_gs_solver_mold,info,pos='pre')
      call lv%set(mld_d_jac_smoother_mold,info,pos='post')
      if (info == 0) call lv%set(mld_d_bwgs_solver_mold,info,pos='post')
      
    case default
      !
      ! Do nothing and hope for the best :) 
      !
    end select
    if ((ipos_==mld_smooth_pre_).or.(ipos_==mld_smooth_both_)) then 
      if (allocated(lv%sm)) call lv%sm%default()
    end if
    if ((ipos_==mld_smooth_post_).or.(ipos_==mld_smooth_both_)) then
      if (allocated(lv%sm2a)) call lv%sm2a%default()
    end if
    

  case('SUB_SOLVE')
    select case (psb_toupper(trim(val)))
    case ('NONE','NOPREC','FACT_NONE')
      call lv%set(mld_d_id_solver_mold,info,pos=pos)
      
    case ('DIAG')
      call lv%set(mld_d_diag_solver_mold,info,pos=pos)
      
    case ('L1-DIAG')
      call lv%set(mld_d_l1_diag_solver_mold,info,pos=pos)
      
    case ('GS','FGS','FWGS')
      call lv%set(mld_d_gs_solver_mold,info,pos=pos)
      
    case ('BGS','BWGS')
      call lv%set(mld_d_bwgs_solver_mold,info,pos=pos)
      
    case ('ILU','ILUT','MILU')
      call lv%set(mld_d_ilu_solver_mold,info,pos=pos)
      if (info == 0) then
        if ((ipos_==mld_smooth_pre_) .or.(ipos_==mld_smooth_both_)) then 
          call lv%sm%sv%set('SUB_SOLVE',val,info)
        end if
        if ((ipos_==mld_smooth_post_).or.(ipos_==mld_smooth_both_))then 
          if (allocated(lv%sm2a)) call lv%sm2a%sv%set('SUB_SOLVE',val,info)
        end if
      end if
#ifdef HAVE_SLU_
    case ('SLU') 
      call lv%set(mld_d_slu_solver_mold,info,pos=pos)
#endif
#ifdef HAVE_MUMPS_
    case ('MUMPS') 
      call lv%set(mld_d_mumps_solver_mold,info,pos=pos)
#endif
#ifdef HAVE_SLUDIST_
    case ('SLUDIST')
      call lv%set(mld_d_sludist_solver_mold,info,pos=pos)
#endif
#ifdef HAVE_UMF_
    case ('UMF')
      call lv%set(mld_d_umf_solver_mold,info,pos=pos)
#endif
    case default
      !
      ! Do nothing and hope for the best :) 
      !
    end select
    
  case ('ML_CYCLE')
    lv%parms%ml_cycle      = mld_stringval(val)

  case ('PAR_AGGR_ALG')
    ival = mld_stringval(val)
    lv%parms%par_aggr_alg  = ival
    if (allocated(lv%aggr)) then
      call lv%aggr%free(info)
      if (info == 0) deallocate(lv%aggr,stat=info)
      if (info /= 0) then
        info = psb_err_internal_error_
        return
      end if
    end if
    
    select case(ival)
    case(mld_dec_aggr_)
      allocate(mld_d_dec_aggregator_type :: lv%aggr, stat=info)
    case(mld_sym_dec_aggr_)
      allocate(mld_d_symdec_aggregator_type :: lv%aggr, stat=info)
    case default
      info =  psb_err_internal_error_
    end select
    if (info == psb_success_) call lv%aggr%default()
    
  case ('AGGR_ORD')
    lv%parms%aggr_ord      = mld_stringval(val)

  case ('AGGR_TYPE')
    lv%parms%aggr_type     = mld_stringval(val)
    if (allocated(lv%aggr)) call lv%aggr%set_aggr_type(lv%parms,info)

  case ('AGGR_PROL')
    lv%parms%aggr_prol     = mld_stringval(val)

  case ('COARSE_MAT')
    lv%parms%coarse_mat    = mld_stringval(val)

  case ('AGGR_OMEGA_ALG')
    lv%parms%aggr_omega_alg= mld_stringval(val)

  case ('AGGR_EIG')
    lv%parms%aggr_eig      = mld_stringval(val)

  case ('AGGR_FILTER')
    lv%parms%aggr_filter   = mld_stringval(val)

  case ('COARSE_SOLVE')
    lv%parms%coarse_solve    = mld_stringval(val)

  case default
    if ((ipos_==mld_smooth_pre_) .or.(ipos_==mld_smooth_both_)) then 
      if (allocated(lv%sm)) then 
        call lv%sm%set(what,val,info,idx=idx)
      end if
    end if
    if ((ipos_==mld_smooth_post_).or.(ipos_==mld_smooth_both_))then 
      if (allocated(lv%sm2a)) then 
        call lv%sm2a%set(what,val,info,idx=idx)
      end if
    end if
    if (allocated(lv%aggr)) call lv%aggr%set(what,val,info,idx=idx)

  end select


  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine mld_d_base_onelev_csetc
