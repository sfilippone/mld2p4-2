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
! File: mld_base_prec_type.f90
!
! Module: mld_base_prec_type
!
!  Constants and utilities in common to all type variants of MLD preconditioners.
!  - integer constants defining the preconditioner;
!  - character constants describing the preconditioner (used by the routines
!    printing out a preconditioner description);
!  - the interfaces to the routines for the management of the preconditioner
!    data structure (see below).
!
!  It contains routines for
!  - converting character constants defining the preconditioner into integer
!    constants; 
!  - checking if the preconditioner is correctly defined;
!  - printing a	description of the preconditioner;
!  - deallocating the preconditioner data structure.  
!

module mld_base_prec_type

  !
  ! This reduces the size of .mod file. Without the ONLY clause compilation 
  ! blows up on some systems.
  !
  use psb_const_mod
  use psb_sparse_mod, only :&
       & psb_d_sparse_mat, psb_z_sparse_mat,&
       & psb_s_sparse_mat, psb_c_sparse_mat,&
       & psb_desc_type,&
       & psb_slinmap_type, psb_dlinmap_type,&
       & psb_clinmap_type, psb_zlinmap_type, &
       & psb_dpk_, psb_spk_, psb_long_int_k_,  &
       & psb_spfree, psb_cdfree, psb_halo_, psb_none_, psb_sum_, psb_avg_, &
       & psb_nohalo_, psb_square_root_, psb_toupper, psb_root_,&
       & psb_sizeof_int, psb_sizeof_long_int, psb_sizeof_sp, psb_sizeof_dp, psb_sizeof,&
       & psb_cd_get_context, psb_info
  use psb_prec_mod, only: psb_sprec_type, psb_dprec_type,&
       & psb_cprec_type, psb_zprec_type, psb_d_base_prec_type


  type mld_aux_onelev_map_type
    integer              :: naggr
    integer, allocatable :: ilaggr(:)
  end type mld_aux_onelev_map_type
  type mld_aux_map_type
    type(mld_aux_onelev_map_type), allocatable :: mapv(:)
  end type mld_aux_map_type
    

  !
  ! Entries in iprcparm
  !
  ! These are in baseprec
  ! 
  integer, parameter :: mld_smoother_type_   =  1         
  integer, parameter :: mld_sub_solve_       =  2
  integer, parameter :: mld_sub_restr_       =  3
  integer, parameter :: mld_sub_prol_        =  4
  integer, parameter :: mld_sub_ren_         =  5
  integer, parameter :: mld_sub_ovr_         =  6
  integer, parameter :: mld_sub_fillin_      =  8
  !! 2 ints for 64 bit versions
  integer, parameter :: mld_slu_ptr_         = 10
  integer, parameter :: mld_umf_symptr_      = 12
  integer, parameter :: mld_umf_numptr_      = 14
  integer, parameter :: mld_slud_ptr_        = 16
  integer, parameter :: mld_prec_status_     = 18 
  !
  ! These are in onelev
  ! 
  integer, parameter :: mld_ml_type_              = 20
  integer, parameter :: mld_smoother_sweeps_pre_  = 21
  integer, parameter :: mld_smoother_sweeps_post_ = 22
  integer, parameter :: mld_smoother_pos_         = 23
  integer, parameter :: mld_aggr_kind_            = 24
  integer, parameter :: mld_aggr_alg_             = 25
  integer, parameter :: mld_aggr_omega_alg_       = 26
  integer, parameter :: mld_aggr_eig_             = 27
  integer, parameter :: mld_aggr_filter_          = 28
  integer, parameter :: mld_coarse_mat_           = 29
  integer, parameter :: mld_coarse_solve_         = 30 
  integer, parameter :: mld_coarse_sweeps_        = 31
  integer, parameter :: mld_coarse_fillin_        = 32
  integer, parameter :: mld_coarse_subsolve_      = 33
  integer, parameter :: mld_smoother_sweeps_      = 34
  integer, parameter :: mld_ifpsz_                = 36

  !
  ! Legal values for entry: mld_smoother_type_
  ! 
  integer, parameter :: mld_min_prec_=0, mld_noprec_=0, mld_jac_=1, mld_bjac_=2,&
       & mld_as_=3, mld_max_prec_=3
  !  VERY IMPORTANT: we are relying on the following to be true: 
  !                  mld_pjac_ == mld_diag_scale_
  !                  mld_bjac_ == mld_milu_n_ (or mld_ilu_n_ would be fine)
  !        mld_diag_scale_ < min(mld_slu_ mld_umf_, mld_sludist_)
  !  This is a quick&dirty fix, but I have nothing better now...
  !
  ! Legal values for entry: mld_sub_solve_
  !
  integer, parameter :: mld_f_none_=0, mld_diag_scale_=1, mld_ilu_n_=2, mld_milu_n_=3
  integer, parameter :: mld_ilu_t_=4, mld_slu_=5, mld_umf_=6, mld_sludist_=7
  integer, parameter :: mld_max_sub_solve_= 7
  !
  ! Legal values for entry: mld_sub_ren_
  !
  integer, parameter :: mld_renum_none_=0, mld_renum_glb_=1, mld_renum_gps_=2
  ! For the time being we are disabling GPS renumbering.
  integer, parameter :: mld_max_renum_=1
  !
  ! Legal values for entry: mld_ml_type_
  !
  integer, parameter :: mld_no_ml_=0, mld_add_ml_=1, mld_mult_ml_=2
  integer, parameter :: mld_new_ml_prec_=3, mld_max_ml_type_=mld_mult_ml_
  !
  ! Legal values for entry: mld_smoother_pos_
  !
  integer, parameter :: mld_pre_smooth_=1, mld_post_smooth_=2,&
       &  mld_twoside_smooth_=3, mld_max_smooth_=mld_twoside_smooth_
  !
  ! Legal values for entry: mld_aggr_kind_
  !
  integer, parameter :: mld_no_smooth_=0, mld_smooth_prol_=1
  integer, parameter :: mld_min_energy_=2, mld_biz_prol_=3
  ! Disabling biz_prol for the time being.
  integer, parameter :: mld_max_aggr_kind_=mld_min_energy_
  !
  ! Legal values for entry: mld_aggr_filter_
  !
  integer, parameter :: mld_no_filter_mat_=0, mld_filter_mat_=1
  integer, parameter :: mld_max_filter_mat_=mld_no_filter_mat_
  !  
  ! Legal values for entry: mld_aggr_alg_
  !
  integer, parameter :: mld_dec_aggr_=0, mld_sym_dec_aggr_=1
  integer, parameter :: mld_glb_aggr_=2, mld_new_dec_aggr_=3
  integer, parameter :: mld_new_glb_aggr_=4
  integer, parameter :: mld_max_aggr_alg_=mld_dec_aggr_

  !
  ! Legal values for entry: mld_aggr_omega_alg_
  !
  integer, parameter :: mld_eig_est_=0, mld_user_choice_=999
  !
  ! Legal values for entry: mld_aggr_eig_
  !
  integer, parameter :: mld_max_norm_=0
  !
  ! Legal values for entry: mld_coarse_mat_
  !
  integer, parameter :: mld_distr_mat_=0, mld_repl_mat_=1
  integer, parameter :: mld_max_coarse_mat_=mld_repl_mat_  
  !
  ! Legal values for entry: mld_prec_status_
  !
  integer, parameter :: mld_prec_built_=98765

  !
  ! Entries in rprcparm: ILU(k,t) threshold, smoothed aggregation omega
  !
  integer, parameter :: mld_sub_iluthrs_    = 1
  integer, parameter :: mld_aggr_omega_val_ = 2
  integer, parameter :: mld_aggr_thresh_    = 3
  integer, parameter :: mld_coarse_iluthrs_ = 4
  integer, parameter :: mld_rfpsz_          = 8

  !
  ! Fields for sparse matrices ensembles stored in av()
  ! 
  integer, parameter :: mld_l_pr_=1, mld_u_pr_=2, mld_bp_ilu_avsz_=2
  integer, parameter :: mld_ap_nd_=3, mld_ac_=4, mld_sm_pr_t_=5, mld_sm_pr_=6
  integer, parameter :: mld_smth_avsz_=6, mld_max_avsz_=mld_smth_avsz_ 

  !
  ! Character constants used by mld_file_prec_descr
  !
  character(len=19), parameter, private :: &
       &  eigen_estimates(0:0)=(/'infinity norm     '/)
  character(len=19), parameter, private :: &
       &  smooth_pos_names(1:3)=(/'pre-smoothing     ','post-smoothing    ',&
       & 'pre/post-smoothing'/)
  character(len=15), parameter, private :: &
       &  aggr_kinds(0:3)=(/'nonsmoothed   ','smoothed      ',&
       &           'min energy    ','bizr. smoothed'/)
  character(len=15), parameter, private :: &
       &  matrix_names(0:1)=(/'distributed   ','replicated    '/)
  character(len=18), parameter, private :: &
       &  aggr_names(0:4)=(/'local aggregation ','sym. local aggr.  ',&
       &     'global aggregation', 'new local aggr.   ','new global aggr.  '/)
  character(len=6), parameter, private :: &
       &  restrict_names(0:4)=(/'none ','halo ','     ','     ','     '/)
  character(len=12), parameter, private :: &
       &  prolong_names(0:3)=(/'none       ','sum        ','average    ','square root'/)
  character(len=15), parameter, private :: &
       &  ml_names(0:3)=(/'none          ','additive      ','multiplicative',&
       & 'new ML        '/)
  character(len=15), parameter, private :: &
       &  fact_names(0:7)=(/'none          ','Point Jacobi  ','ILU(n)        ',&
       &  'MILU(n)       ','ILU(t,n)      ',&
       &  'UMFPACK LU    ',&
       &  'SuperLU_Dist  ','SuperLU       '/)

  interface mld_check_def
    module procedure mld_icheck_def, mld_scheck_def, mld_dcheck_def
  end interface

contains

  !
  ! Subroutine: mld_stringval
  !
  !  This routine converts the string contained into string into the corresponding
  !  integer value.
  !
  ! Arguments:
  !    string  -  character(len=*), input.
  !               The string to be converted.
  !    val     -  integer, output.
  !               The integer value corresponding to the string
  !    info    -  integer, output.
  !               Error code.
  !
  subroutine mld_stringval(string,val,info)
    implicit none 
  ! Arguments
    character(len=*), intent(in) :: string
    integer, intent(out) :: val, info
    character(len=*), parameter :: name='mld_stringval'
    
    info = psb_success_
    select case(psb_toupper(trim(string)))
    case('NONE')
      val = 0
    case('HALO')
      val = psb_halo_ 
    case('SUM')
      val = psb_sum_
    case('AVG')
      val = psb_avg_
    case('FACT_NONE')
      val = mld_f_none_
    case('ILU')
      val = mld_ilu_n_
    case('MILU')
      val = mld_milu_n_
    case('ILUT')
      val = mld_ilu_t_
    case('UMF')
      val = mld_umf_
    case('SLU')
      val = mld_slu_
    case('SLUDIST')
      val = mld_sludist_
    case('DSCALE')
      val = mld_diag_scale_
    case('ADD')
      val = mld_add_ml_
    case('MULT')
      val = mld_mult_ml_
    case('DEC')
      val = mld_dec_aggr_
    case('SYMDEC')
      val = mld_sym_dec_aggr_
    case('GLB')
      val = mld_glb_aggr_
    case('REPL')
      val = mld_repl_mat_
    case('DIST')
      val = mld_distr_mat_
    case('NONSMOOTHED')
      val = mld_no_smooth_
    case('SMOOTHED')
      val = mld_smooth_prol_
    case('MINENERGY')
      val = mld_min_energy_
    case('PRE')
      val = mld_pre_smooth_
    case('POST')
      val = mld_post_smooth_
    case('TWOSIDE')
      val = mld_twoside_smooth_
    case('NOPREC')
      val = mld_noprec_
!!$    case('DIAG')
!!$      val = mld_diag_
    case('BJAC')
      val = mld_bjac_
    case('JAC','JACOBI')
      val = mld_jac_
    case('AS')
      val = mld_as_
    case('RENUM_NONE')
      val = mld_renum_none_
    case('RENUM_GLOBAL')
      val = mld_renum_glb_
    case('RENUM_GPS')
      val = mld_renum_gps_
    case('A_NORMI')
      val = mld_max_norm_
    case('USER_CHOICE')
      val = mld_user_choice_
    case('EIG_EST')
      val = mld_eig_est_
    case('FILTER')
      val = mld_filter_mat_
    case('NO_FILTER')
      val = mld_no_filter_mat_
    case default
      val  = -1
      info = -1
    end select
    if (info /= psb_success_) then 
      write(0,*) name,': Error: unknown request: "',trim(string),'"'
    end if
  end subroutine mld_stringval

    
  !
  ! Routines printing out a description of the preconditioner
  !

  subroutine mld_base_prec_descr(iout,iprcparm, info,rprcparm,dprcparm)
    implicit none 
    integer, intent(in) :: iprcparm(:),iout
    integer, intent(out) :: info
    real(psb_spk_), intent(in), optional :: rprcparm(:)
    real(psb_dpk_), intent(in), optional :: dprcparm(:)
    
    info = psb_success_
    if (count((/ present(rprcparm),present(dprcparm) /)) /= 1) then 
      info=psb_err_no_optional_arg_
!!$      call psb_errpush(info,name,a_err=" rprcparm, dprcparm")
      return
    endif
    
    select case(iprcparm(mld_smoother_type_))
    case(mld_noprec_)
      write(iout,*) '  No preconditioning'
    case(mld_jac_)
      write(iout,*) '  Jacobi '
    case(mld_bjac_)
      write(iout,*) '  Block Jacobi with ',&
           &  fact_names(iprcparm(mld_sub_solve_))
      select case(iprcparm(mld_sub_solve_))
      case(mld_ilu_n_,mld_milu_n_)      
        write(iout,*) '  Fill level:',iprcparm(mld_sub_fillin_)
      case(mld_ilu_t_)         
        write(iout,*) '  Fill level:',iprcparm(mld_sub_fillin_)
        if (present(rprcparm)) then 
          write(iout,*) '  Fill threshold :',rprcparm(mld_sub_iluthrs_)
        else
          write(iout,*) '  Fill threshold :',dprcparm(mld_sub_iluthrs_)
        end if
      case(mld_slu_,mld_umf_,mld_sludist_,mld_diag_scale_) 
      case default
        write(iout,*) '  Should never get here!'
      end select
    case(mld_as_)
      write(iout,*) '  Additive Schwarz with ',&
           &  fact_names(iprcparm(mld_sub_solve_))
      select case(iprcparm(mld_sub_solve_))
      case(mld_ilu_n_,mld_milu_n_)      
        write(iout,*) '  Fill level:',iprcparm(mld_sub_fillin_)
      case(mld_ilu_t_)
        write(iout,*) '  Fill level:',iprcparm(mld_sub_fillin_)
        if (present(rprcparm)) then 
          write(iout,*) '  Fill threshold :',rprcparm(mld_sub_iluthrs_)
        else
          write(iout,*) '  Fill threshold :',dprcparm(mld_sub_iluthrs_)
        end if
      case(mld_slu_,mld_umf_,mld_sludist_,mld_diag_scale_) 
      case default
        write(iout,*) '  Should never get here!'
      end select
      write(iout,*) '  Overlap:',&
           &  iprcparm(mld_sub_ovr_)
      write(iout,*) '  Restriction: ',&
           &  restrict_names(iprcparm(mld_sub_restr_))
      write(iout,*) '  Prolongation: ',&
           &  prolong_names(iprcparm(mld_sub_prol_))
    end select
    return
  end subroutine mld_base_prec_descr

  subroutine mld_ml_alg_descr(iout,ilev,iprcparm, info,rprcparm,dprcparm)
    implicit none 
    integer, intent(in) :: iprcparm(:),iout,ilev
    integer, intent(out) :: info
    real(psb_spk_), intent(in), optional :: rprcparm(:)
    real(psb_dpk_), intent(in), optional :: dprcparm(:)
    integer :: sweeps

    info = psb_success_
    if (count((/ present(rprcparm),present(dprcparm) /)) /= 1) then 
      info=psb_err_no_optional_arg_
!!$      call psb_errpush(info,name,a_err=" rprcparm, dprcparm")
      return
    endif

    if (iprcparm(mld_ml_type_)>mld_no_ml_) then

      write(iout,*) '  Multilevel type: ',&
           &   ml_names(iprcparm(mld_ml_type_))
      write(iout,*) '  Smoother position: ',&
           & smooth_pos_names(iprcparm(mld_smoother_pos_))
      if (iprcparm(mld_ml_type_) == mld_add_ml_) then
        write(iout,*) '  Number of sweeps : ',&
             & iprcparm(mld_smoother_sweeps_) 
      else 
        select case (iprcparm(mld_smoother_pos_))
        case (mld_pre_smooth_)
          write(iout,*) '  Number of sweeps : ',&
               & iprcparm(mld_smoother_sweeps_pre_)
        case (mld_post_smooth_)
          write(iout,*) '  Number of sweeps : ',&
               &  iprcparm(mld_smoother_sweeps_post_)
        case (mld_twoside_smooth_)
          write(iout,*) '  Number of sweeps : pre: ',&
               &  iprcparm(mld_smoother_sweeps_pre_) ,&
               &  '  post: ',&
               &  iprcparm(mld_smoother_sweeps_post_)
        end select
      end if

      write(iout,*) '  Aggregation: ', &
           &   aggr_names(iprcparm(mld_aggr_alg_))
      write(iout,*) '  Aggregation type: ', &
           &  aggr_kinds(iprcparm(mld_aggr_kind_))
      if (present(rprcparm)) then 
        write(iout,*) '  Aggregation threshold: ', &
             &  rprcparm(mld_aggr_thresh_)
      else
        write(iout,*) '  Aggregation threshold: ', &
             &  dprcparm(mld_aggr_thresh_)
      end if
      if (iprcparm(mld_aggr_kind_) /= mld_no_smooth_) then
         if (iprcparm(mld_aggr_omega_alg_) == mld_eig_est_) then 
          write(iout,*) '  Damping omega computation: spectral radius estimate'
          write(iout,*) '  Spectral radius estimate: ', &
               & eigen_estimates(iprcparm(mld_aggr_eig_))
        else if (iprcparm(mld_aggr_omega_alg_) == mld_user_choice_) then 
          write(iout,*) '  Damping omega computation: user defined value.'
        else 
          write(iout,*) '  Damping omega computation: unknown value in iprcparm!!'
        end if
      end if
    end if

    return
  end subroutine mld_ml_alg_descr

  subroutine mld_ml_level_descr(iout,ilev,iprcparm,nlaggr, info,rprcparm,dprcparm)
    implicit none 
    integer, intent(in) :: iprcparm(:),iout,ilev
    integer, intent(in), allocatable :: nlaggr(:)
    integer, intent(out) :: info
    real(psb_spk_), intent(in), optional :: rprcparm(:)
    real(psb_dpk_), intent(in), optional :: dprcparm(:)

    info = psb_success_
    if (count((/ present(rprcparm),present(dprcparm) /)) /= 1) then 
      info=psb_err_no_optional_arg_
!!$      call psb_errpush(info,name,a_err=" rprcparm, dprcparm")
      return
    endif

    if (iprcparm(mld_ml_type_)>mld_no_ml_) then
      write(iout,*) ' Level ',ilev
      if (allocated(nlaggr)) then
        write(iout,*) '  Size of coarse matrix: ', &
             &  sum(nlaggr(:))
        write(iout,*) '  Sizes of aggregates: ', &
             &  nlaggr(:)
      end if
      if (iprcparm(mld_aggr_kind_) /= mld_no_smooth_) then
        if (present(rprcparm)) then 
          write(iout,*) '  Damping omega: ', &
               & rprcparm(mld_aggr_omega_val_)  
        else
          write(iout,*) '  Damping omega: ', &
               & dprcparm(mld_aggr_omega_val_)  
        end if
      end if
    end if
    
    return
  end subroutine mld_ml_level_descr

  subroutine mld_ml_coarse_descr(iout,ilev,iprcparm,iprcparm2,nlaggr,info,&
       & rprcparm,dprcparm, rprcparm2,dprcparm2)
    implicit none 
    integer, intent(in) :: iprcparm(:),iprcparm2(:),iout,ilev
    integer, intent(in), allocatable :: nlaggr(:)
    integer, intent(out) :: info
    real(psb_spk_), intent(in), optional :: rprcparm(:), rprcparm2(:)
    real(psb_dpk_), intent(in), optional :: dprcparm(:), dprcparm2(:)

    info = psb_success_
    if (count((/ present(rprcparm),present(dprcparm) /)) /= 1) then 
      info=psb_err_no_optional_arg_
!!$      call psb_errpush(info,name,a_err=" rprcparm, dprcparm")
      return
    endif
    if (count((/ present(rprcparm2),present(dprcparm2) /)) /= 1) then 
      info=psb_err_no_optional_arg_
!!$      call psb_errpush(info,name,a_err=" rprcparm, dprcparm")
      return
    endif

    if (iprcparm(mld_ml_type_)>mld_no_ml_) then

      write(iout,*) ' Level ',ilev,' (coarsest)'
      write(iout,*) '  Coarsest matrix: ',&
           & matrix_names(iprcparm(mld_coarse_mat_))
      if (allocated(nlaggr)) then 
        write(iout,*) '  Size of coarsest matrix: ', &
             &  sum( nlaggr(:))
        write(iout,*) '  Sizes of aggregates: ', &
             &  nlaggr(:)
      end if
      if (iprcparm(mld_aggr_kind_) /= mld_no_smooth_) then 

        if (present(rprcparm)) then 
          write(iout,*) '  Damping omega: ', &
               & rprcparm(mld_aggr_omega_val_)  
        else
          write(iout,*) '  Damping omega: ', &
               & dprcparm(mld_aggr_omega_val_)  
        end if
      end if
      if (iprcparm(mld_coarse_mat_) == mld_distr_mat_ .and. &
           & iprcparm(mld_sub_solve_) /= mld_sludist_) then
!!$        write(iout,*) '  Coarsest matrix solver: ',&
!!$             & smoother_names(iprcparm2(mld_smoother_type_))
        select case (iprcparm2(mld_smoother_type_))
        case(mld_bjac_,mld_as_)
          write(iout,*) '  subdomain solver: ',&          
               &  fact_names(iprcparm2(mld_sub_solve_))
          write(iout,*) '  Number of smoother sweeps: ', &
               &   (iprcparm2(mld_smoother_sweeps_))
        case(mld_jac_) 
          write(iout,*) '  Number of smoother sweeps: ', &
               &   (iprcparm2(mld_smoother_sweeps_))
        end select
      else
        write(iout,*) '  Coarsest matrix solver: ', &
             &  fact_names(iprcparm2(mld_sub_solve_))
      end if
      select case(iprcparm2(mld_sub_solve_))
      case(mld_ilu_n_,mld_milu_n_)      
        write(iout,*) '  Fill level:',iprcparm2(mld_sub_fillin_)
      case(mld_ilu_t_)
        write(iout,*) '  Fill level:',iprcparm2(mld_sub_fillin_)
        if (present(rprcparm2)) then 
          write(iout,*) '  Fill threshold :',rprcparm2(mld_sub_iluthrs_)
        else if (present(dprcparm2)) then 
          write(iout,*) '  Fill threshold :',dprcparm2(mld_sub_iluthrs_)
        end if
      case(mld_slu_,mld_umf_,mld_sludist_,mld_diag_scale_) 
      case default
        write(iout,*) '  Should never get here!'
      end select
    end if


    return
  end subroutine mld_ml_coarse_descr


  subroutine mld_ml_new_coarse_descr(iout,ilev,iprcparm,nlaggr,info,&
       & rprcparm,dprcparm)
    implicit none 
    integer, intent(in) :: iprcparm(:),iout,ilev
    integer, intent(in), allocatable :: nlaggr(:)
    integer, intent(out) :: info
    real(psb_spk_), intent(in), optional :: rprcparm(:)
    real(psb_dpk_), intent(in), optional :: dprcparm(:)

    info = psb_success_
    if (count((/ present(rprcparm),present(dprcparm) /)) /= 1) then 
      info=psb_err_no_optional_arg_
!!$      call psb_errpush(info,name,a_err=" rprcparm, dprcparm")
      return
    endif
!!$    if (count((/ present(rprcparm2),present(dprcparm2) /)) /= 1) then 
!!$      info=psb_err_no_optional_arg_
! !$      call psb_errpush(info,name,a_err=" rprcparm, dprcparm")
!!$      return
!!$    endif

    if (iprcparm(mld_ml_type_)>mld_no_ml_) then

      write(iout,*) ' Level ',ilev,' (coarsest)'
      write(iout,*) '  Coarsest matrix: ',&
           & matrix_names(iprcparm(mld_coarse_mat_))
      if (allocated(nlaggr)) then 
        write(iout,*) '  Size of coarsest matrix: ', &
             &  sum( nlaggr(:))
        write(iout,*) '  Sizes of aggregates: ', &
             &  nlaggr(:)
      end if
      if (iprcparm(mld_aggr_kind_) /= mld_no_smooth_) then 

        if (present(rprcparm)) then 
          write(iout,*) '  Damping omega: ', &
               & rprcparm(mld_aggr_omega_val_)  
        else
          write(iout,*) '  Damping omega: ', &
               & dprcparm(mld_aggr_omega_val_)  
        end if
      end if
!!$      if (iprcparm(mld_coarse_mat_) == mld_distr_mat_ .and. &
!!$           & iprcparm(mld_sub_solve_) /= mld_sludist_) then
! !$        write(iout,*) '  Coarsest matrix solver: ',&
! !$             & smoother_names(iprcparm2(mld_smoother_type_))
!!$        select case (iprcparm2(mld_smoother_type_))
!!$        case(mld_bjac_,mld_as_)
!!$          write(iout,*) '  subdomain solver: ',&          
!!$               &  fact_names(iprcparm2(mld_sub_solve_))
!!$          write(iout,*) '  Number of smoother sweeps: ', &
!!$               &   (iprcparm2(mld_smoother_sweeps_))
!!$        case(mld_jac_) 
!!$          write(iout,*) '  Number of smoother sweeps: ', &
!!$               &   (iprcparm2(mld_smoother_sweeps_))
!!$        end select
!!$      else
!!$        write(iout,*) '  Coarsest matrix solver: ', &
!!$             &  fact_names(iprcparm2(mld_sub_solve_))
!!$      end if
!!$      select case(iprcparm2(mld_sub_solve_))
!!$      case(mld_ilu_n_,mld_milu_n_)      
!!$        write(iout,*) '  Fill level:',iprcparm2(mld_sub_fillin_)
!!$      case(mld_ilu_t_)
!!$        write(iout,*) '  Fill level:',iprcparm2(mld_sub_fillin_)
!!$        if (present(rprcparm2)) then 
!!$          write(iout,*) '  Fill threshold :',rprcparm2(mld_sub_iluthrs_)
!!$        else if (present(dprcparm2)) then 
!!$          write(iout,*) '  Fill threshold :',dprcparm2(mld_sub_iluthrs_)
!!$        end if
!!$      case(mld_slu_,mld_umf_,mld_sludist_,mld_diag_scale_) 
!!$      case default
!!$        write(iout,*) '  Should never get here!'
!!$      end select
    end if


    return
  end subroutine mld_ml_new_coarse_descr


  !
  ! Functions/subroutines checking if the preconditioner is correctly defined
  !

  function is_legal_base_prec(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_legal_base_prec

    is_legal_base_prec = ((ip>=mld_noprec_).and.(ip<=mld_max_prec_))
    return
  end function is_legal_base_prec
  function is_legal_n_ovr(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_legal_n_ovr

    is_legal_n_ovr = (ip >= 0) 
    return
  end function is_legal_n_ovr
  function is_legal_renum(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_legal_renum
    is_legal_renum = ((ip >= 0).and.(ip <= mld_max_renum_))
    return
  end function is_legal_renum
  function is_legal_jac_sweeps(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_legal_jac_sweeps

    is_legal_jac_sweeps = (ip >= 1) 
    return
  end function is_legal_jac_sweeps
  function is_legal_prolong(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_legal_prolong
    is_legal_prolong = ((ip>=psb_none_).and.(ip<=psb_square_root_))
    return
  end function is_legal_prolong
  function is_legal_restrict(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_legal_restrict
    is_legal_restrict = ((ip == psb_nohalo_).or.(ip==psb_halo_))
    return
  end function is_legal_restrict
  function is_legal_ml_type(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_legal_ml_type

    is_legal_ml_type = ((ip>=mld_no_ml_).and.(ip<=mld_max_ml_type_))
    return
  end function is_legal_ml_type
  function is_legal_ml_aggr_alg(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_legal_ml_aggr_alg

    is_legal_ml_aggr_alg = ((ip>=mld_dec_aggr_).and.(ip<=mld_max_aggr_alg_))
    return
  end function is_legal_ml_aggr_alg
  function is_legal_ml_aggr_omega_alg(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_legal_ml_aggr_omega_alg

    is_legal_ml_aggr_omega_alg = ((ip == mld_eig_est_).or.(ip==mld_user_choice_))
    return
  end function is_legal_ml_aggr_omega_alg
  function is_legal_ml_aggr_eig(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_legal_ml_aggr_eig

    is_legal_ml_aggr_eig = (ip == mld_max_norm_)
    return
  end function is_legal_ml_aggr_eig
  function is_legal_ml_smooth_pos(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_legal_ml_smooth_pos

    is_legal_ml_smooth_pos = ((ip>=mld_pre_smooth_).and.(ip<=mld_max_smooth_))
    return
  end function is_legal_ml_smooth_pos
  function is_legal_ml_aggr_kind(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_legal_ml_aggr_kind

    is_legal_ml_aggr_kind = ((ip>=0).and.(ip<=mld_max_aggr_kind_))
    return
  end function is_legal_ml_aggr_kind
  function is_legal_ml_coarse_mat(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_legal_ml_coarse_mat

    is_legal_ml_coarse_mat = ((ip>=0).and.(ip<=mld_max_coarse_mat_))
    return
  end function is_legal_ml_coarse_mat
  function is_legal_aggr_filter(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_legal_aggr_filter

    is_legal_aggr_filter = ((ip>=0).and.(ip<=mld_max_filter_mat_))
    return
  end function is_legal_aggr_filter
  function is_distr_ml_coarse_mat(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_distr_ml_coarse_mat

    is_distr_ml_coarse_mat = (ip == mld_distr_mat_)
    return
  end function is_distr_ml_coarse_mat
  function is_legal_ml_fact(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_legal_ml_fact
    ! Here the minimum is really 1, mld_fact_none_ is not acceptable.
    is_legal_ml_fact = ((ip>=1).and.(ip<=mld_max_sub_solve_))
    return
  end function is_legal_ml_fact
  function is_legal_ml_lev(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_legal_ml_lev

    is_legal_ml_lev = (ip >= 0)
    return
  end function is_legal_ml_lev
  function is_legal_omega(ip)
    implicit none 
    real(psb_dpk_), intent(in) :: ip
    logical             :: is_legal_omega
    is_legal_omega = ((ip>=0.0d0).and.(ip<=2.0d0))
    return
  end function is_legal_omega
  function is_legal_fact_thrs(ip)
    implicit none 
    real(psb_dpk_), intent(in) :: ip
    logical             :: is_legal_fact_thrs

    is_legal_fact_thrs = (ip>=0.0d0)
    return
  end function is_legal_fact_thrs
  function is_legal_aggr_thrs(ip)
    implicit none 
    real(psb_dpk_), intent(in) :: ip
    logical             :: is_legal_aggr_thrs

    is_legal_aggr_thrs = (ip>=0.0d0)
    return
  end function is_legal_aggr_thrs

  function is_legal_s_omega(ip)
    implicit none 
    real(psb_spk_), intent(in) :: ip
    logical             :: is_legal_s_omega
    is_legal_s_omega = ((ip>=0.0).and.(ip<=2.0))
    return
  end function is_legal_s_omega
  function is_legal_s_fact_thrs(ip)
    implicit none 
    real(psb_spk_), intent(in) :: ip
    logical             :: is_legal_s_fact_thrs

    is_legal_s_fact_thrs = (ip>=0.0)
    return
  end function is_legal_s_fact_thrs
  function is_legal_s_aggr_thrs(ip)
    implicit none 
    real(psb_spk_), intent(in) :: ip
    logical             :: is_legal_s_aggr_thrs

    is_legal_s_aggr_thrs = (ip>=0.0)
    return
  end function is_legal_s_aggr_thrs


  subroutine mld_icheck_def(ip,name,id,is_legal)
    implicit none 
    integer, intent(inout) :: ip
    integer, intent(in)    :: id
    character(len=*), intent(in) :: name
    interface 
      function is_legal(i)
        integer, intent(in) :: i
        logical             :: is_legal
      end function is_legal
    end interface
    character(len=20), parameter :: rname='mld_check_def'

    if (.not.is_legal(ip)) then     
      write(0,*)trim(rname),': Error: Illegal value for ',&
           & name,' :',ip, '. defaulting to ',id
      ip = id
    end if
  end subroutine mld_icheck_def

  subroutine mld_scheck_def(ip,name,id,is_legal)
    implicit none 
    real(psb_spk_), intent(inout) :: ip
    real(psb_spk_), intent(in)    :: id
    character(len=*), intent(in) :: name
    interface 
      function is_legal(i)
        use psb_sparse_mod, only : psb_spk_
        real(psb_spk_), intent(in) :: i
        logical             :: is_legal
      end function is_legal
    end interface
    character(len=20), parameter :: rname='mld_check_def'
    
    if (.not.is_legal(ip)) then     
      write(0,*)trim(rname),': Error: Illegal value for ',&
           & name,' :',ip, '. defaulting to ',id
      ip = id
    end if
  end subroutine mld_scheck_def

  subroutine mld_dcheck_def(ip,name,id,is_legal)
    implicit none 
    real(psb_dpk_), intent(inout) :: ip
    real(psb_dpk_), intent(in)    :: id
    character(len=*), intent(in) :: name
    interface 
      function is_legal(i)
        use psb_sparse_mod, only : psb_dpk_
        real(psb_dpk_), intent(in) :: i
        logical             :: is_legal
      end function is_legal
    end interface
    character(len=20), parameter :: rname='mld_check_def'
    
    if (.not.is_legal(ip)) then     
      write(0,*)trim(rname),': Error: Illegal value for ',&
           & name,' :',ip, '. defaulting to ',id
      ip = id
    end if
  end subroutine mld_dcheck_def


  function pr_to_str(iprec)
    implicit none 

    integer, intent(in)  :: iprec
    character(len=10)     :: pr_to_str

    select case(iprec)
    case(mld_noprec_)
      pr_to_str='NOPREC'
    case(mld_jac_)         
      pr_to_str='JAC'
    case(mld_bjac_)         
      pr_to_str='BJAC'
    case(mld_as_)      
      pr_to_str='AS'
    end select

  end function pr_to_str


end module mld_base_prec_type
