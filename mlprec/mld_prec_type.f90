!!$
!!$ 
!!$                           MLD2P4  version 1.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 2.2)
!!$  
!!$  (C) Copyright 2008
!!$
!!$                      Salvatore Filippone  University of Rome Tor Vergata       
!!$                      Alfredo Buttari      University of Rome Tor Vergata
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
! File: mld_prec_type.f90
!
! Package: mld_prec_type
!          Data structure(s) for sparse matrices
!
!  This module defines: 
!  - the mld_prec_type data structure containing the preconditioner;
!  - character constants describing the preconditioner, used by the routine
!    printing out a preconditioner description;
!  - the interfaces to the routines for the management of the preconditioner
!    data structure (see below).
!
!  It contains:
!  - mld_dprec_sizeof, mld_dbaseprc_sizeof, mld_out_prec_descr, mld_file_prec_descr,
!    mld_icheck_def, mld_dcheck_def, mld_dbase_precfree, mld_nullify_dbaseprec,
!    is_legal_..._..., and their complex versions (if applicable).
!  These routines check if the preconditioner is correctly defined,      print a
!  description of the preconditioner, and deallocate its data structure.  
!

module mld_prec_type

  !
  ! This reduces the size of .mod file. Without the ONLY clause compilation 
  ! blows up on some systems.
  !
  use psb_base_mod, only :&
       & psb_dspmat_type, psb_zspmat_type,&
       & psb_sspmat_type, psb_cspmat_type,&
       & psb_desc_type, psb_inter_desc_type,&
       & psb_dpk_, psb_spk_, psb_long_int_k_,  &
       & psb_sp_free, psb_cdfree, psb_halo_, psb_none_, psb_sum_, psb_avg_, &
       & psb_nohalo_, psb_square_root_, psb_toupper, psb_root_,&
       & psb_sizeof_int, psb_sizeof_long_int, psb_sizeof_sp, psb_sizeof_dp, psb_sizeof,&
       & psb_cd_get_context, psb_info

  !
  ! Type: mld_dprec_type, mld_zprec_type, mld_sprec_type, mld_cprec_type
  !
  !  mld_dprec_type and friends are  the real and complex preconditioner
  !  data structures. In the following description 'd', 's', 'c'  and 'z'
  !  are mostly  omitted.
  !
  ! The multilevel preconditioner data structure, mld_Xprec_type, consists
  ! of an array of 'one-level preconditioner' data structures, mld_X_onelev_type,
  ! each containing the local part of the preconditioner associated to a
  ! certain level. For each level ilev, the base preconditioner K(ilev) is 
  ! built from a matrix A(ilev), which is obtained by 'tranferring' the 
  ! original matrix A (i.e. the matrix to be preconditioned) to level ilev,
  ! through smoothed aggregation.
  !
  ! The levels are numbered in increasing order starting from the finest
  ! one, i.e. level 1 is the finest level and A(1) is the matrix A.
  !
  !|  type mld_Xprec_type
  !|    type(mld_X_onelev_prec_type), allocatable :: precv(:) 
  !|  end type mld_Xprec_type
  !|
  ! 
  !   precv(ilev) is the preconditioner at level ilev.
  !   The number of levels is given by size(precv(:)).
  !
  ! Type: mld_X_onelev_prec_type.
  !       The data type containing necessary items for the current level.
  !
  !   type(mld_Xbaseprc_type) -  prec
  !                   The current level preconditioner (aka smoother).
  !   ac           -  The local part of the matrix A(ilev).
  !   desc_ac      -  type(psb_desc_type).
  !                   The communication descriptor associated to the sparse matrix
  !                   A(ilev), stored in ac.
  !   iprcparm     -  integer, dimension(:), allocatable.
  !                   The integer parameters defining the multilevel strategy
  !   rprcparm     -  real(psb_Ypk_), dimension(:), allocatable.
  !                   The real parameters defining the multilevel strategy
  !   mlia         -  integer, dimension(:), allocatable.
  !                   The aggregation map (ilev-1) --> (ilev).
  !                   In case of non-smoothed aggregation, it is used instead of
  !                   mld_sm_pr_.
  !   nlaggr       -  integer, dimension(:), allocatable.
  !                   The number of aggregates (rows of A(ilev)) on the
  !                   various processes.
  !   map_desc     -  Stores the mapping between indices from level(ilev-1) to (ilev).
  !                   Unused at level 1 (finest).
  !   base_a       -  type(psb_zspmat_type), pointer.
  !                   Pointer (really a pointer!) to the local part of the base matrix 
  !                   of the current level, i.e. A(ilev); so we have a unified treatment
  !                   of residuals. We need this to avoid passing explicitly the matrix
  !                   A(ilev) to the routine which applies the preconditioner.
  !   base_desc    -  type(psb_desc_type), pointer.
  !                   Pointer to the communication descriptor associated to the sparse
  !                   matrix pointed by base_a. 
  ! Type: mld_Xbaseprc_type  
  !       The smoother. 
  !
  !    av         -  type(psb_Xspmat_type), dimension(:), allocatable(:).
  !                  The sparse matrices needed to apply the preconditioner at
  !                  the current level ilev. 
  !      av(mld_l_pr_)     -  The L factor of the ILU factorization of the local
  !                           diagonal block of A(ilev).
  !      av(mld_u_pr_)     -  The U factor of the ILU factorization of the local
  !                           diagonal block of A(ilev), except its diagonal entries
  !                           (stored in d).
  !      av(mld_ap_nd_)    -  The entries of the local part of A(ilev) outside
  !                           the diagonal block, for block-Jacobi sweeps.
  !   d            -  real/complex(psb_Ypk_), dimension(:), allocatable.
  !                   The diagonal entries of the U factor in the ILU factorization
  !                   of A(ilev).
  !   desc_data    -  type(psb_desc_type).
  !                   The communication descriptor associated to the base preconditioner,
  !                   i.e. to the sparse matrices needed to apply the base preconditioner
  !                   at the current level.
  !   iprcparm     -  integer, dimension(:), allocatable.
  !                   The integer parameters defining the base preconditioner K(ilev)
  !                   (the iprcparm entries and values are specified below).
  !   rprcparm     -  real(psb_Ypk_), dimension(:), allocatable.
  !                   The real parameters defining the base preconditioner K(ilev)
  !                   (the rprcparm entries and values are specified below).
  !   perm         -  integer, dimension(:), allocatable.
  !                   The row and column permutations applied to the local part of
  !                   A(ilev) (defined only if iprcparm(mld_sub_ren_)>0). 
  !   invperm      -  integer, dimension(:), allocatable.
  !                   The inverse of the permutation stored in perm.
  !
  !   Note that when the LU factorization of the matrix A(ilev) is computed instead of
  !   the ILU one, by using UMFPACK or SuperLU_dist, the corresponding L and U factors
  !   are stored in data structures provided by UMFPACK or SuperLU_dist and pointed by
  !   prec%iprcparm(mld_umf_ptr) or prec%iprcparm(mld_slu_ptr), respectively.
  !

  type mld_sbaseprc_type
    type(psb_sspmat_type), allocatable :: av(:) 
    real(psb_spk_), allocatable        :: d(:)  
    type(psb_desc_type)                :: desc_data
    integer, allocatable               :: iprcparm(:) 
    real(psb_spk_), allocatable        :: rprcparm(:) 
    integer, allocatable               :: perm(:),  invperm(:) 
  end type mld_sbaseprc_type

  type mld_s_onelev_prec_type
    type(mld_sbaseprc_type)            :: prec
    integer, allocatable               :: iprcparm(:) 
    real(psb_spk_), allocatable        :: rprcparm(:) 
    type(psb_sspmat_type)              :: ac
    type(psb_desc_type)                :: desc_ac
    integer, allocatable               :: mlia(:), nlaggr(:) 
    type(psb_sspmat_type), pointer     :: base_a    => null() 
    type(psb_desc_type), pointer       :: base_desc => null() 
    type(psb_inter_desc_type)          :: map_desc
  end type mld_s_onelev_prec_type

  type mld_sprec_type
    type(mld_s_onelev_prec_type), allocatable :: precv(:) 
  end type mld_sprec_type

  type mld_dbaseprc_type
    type(psb_dspmat_type), allocatable :: av(:) 
    real(psb_dpk_), allocatable        :: d(:)  
    type(psb_desc_type)                :: desc_data
    integer, allocatable               :: iprcparm(:) 
    real(psb_dpk_), allocatable        :: rprcparm(:) 
    integer, allocatable               :: perm(:),  invperm(:) 
  end type mld_dbaseprc_type

  type mld_d_onelev_prec_type
    type(mld_dbaseprc_type)            :: prec
    integer, allocatable               :: iprcparm(:) 
    real(psb_dpk_), allocatable        :: rprcparm(:) 
    type(psb_dspmat_type)              :: ac
    type(psb_desc_type)                :: desc_ac
    integer, allocatable               :: mlia(:), nlaggr(:) 
    type(psb_dspmat_type), pointer     :: base_a    => null() 
    type(psb_desc_type), pointer       :: base_desc => null() 
    type(psb_inter_desc_type)          :: map_desc
  end type mld_d_onelev_prec_type

  type mld_dprec_type
    type(mld_d_onelev_prec_type), allocatable :: precv(:) 
  end type mld_dprec_type


  type mld_cbaseprc_type
    type(psb_cspmat_type), allocatable :: av(:) 
    complex(psb_spk_), allocatable     :: d(:)  
    type(psb_desc_type)                :: desc_data
    integer, allocatable               :: iprcparm(:) 
    real(psb_spk_), allocatable        :: rprcparm(:) 
    integer, allocatable               :: perm(:),  invperm(:) 
  end type mld_cbaseprc_type

  type mld_c_onelev_prec_type
    type(mld_cbaseprc_type)            :: prec
    integer, allocatable               :: iprcparm(:) 
    real(psb_spk_), allocatable        :: rprcparm(:) 
    type(psb_cspmat_type)              :: ac
    type(psb_desc_type)                :: desc_ac
    integer, allocatable               :: mlia(:), nlaggr(:) 
    type(psb_cspmat_type), pointer     :: base_a    => null() 
    type(psb_desc_type), pointer       :: base_desc => null() 
    type(psb_inter_desc_type)          :: map_desc
  end type mld_c_onelev_prec_type

  type mld_cprec_type
    type(mld_c_onelev_prec_type), allocatable :: precv(:) 
  end type mld_cprec_type

  type mld_zbaseprc_type
    type(psb_zspmat_type), allocatable :: av(:) 
    complex(psb_dpk_), allocatable     :: d(:)  
    type(psb_desc_type)                :: desc_data
    integer, allocatable               :: iprcparm(:) 
    real(psb_dpk_), allocatable        :: rprcparm(:) 
    integer, allocatable               :: perm(:),  invperm(:) 
  end type mld_zbaseprc_type

  type mld_z_onelev_prec_type
    type(mld_zbaseprc_type)            :: prec
    integer, allocatable               :: iprcparm(:) 
    real(psb_dpk_), allocatable        :: rprcparm(:) 
    type(psb_zspmat_type)              :: ac
    type(psb_desc_type)                :: desc_ac
    integer, allocatable               :: mlia(:), nlaggr(:) 
    type(psb_zspmat_type), pointer     :: base_a    => null() 
    type(psb_desc_type), pointer       :: base_desc => null() 
    type(psb_inter_desc_type)          :: map_desc
  end type mld_z_onelev_prec_type

  type mld_zprec_type
    type(mld_z_onelev_prec_type), allocatable :: precv(:) 
  end type mld_zprec_type


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
  integer, parameter :: mld_smoother_sweeps_ =  9
  !! 2 ints for 64 bit versions
  integer, parameter :: mld_slu_ptr_         = 10
  integer, parameter :: mld_umf_symptr_      = 12
  integer, parameter :: mld_umf_numptr_      = 14
  integer, parameter :: mld_slud_ptr_        = 16
  integer, parameter :: mld_prec_status_     = 18 
  !
  ! These are in onelev_prec
  ! 
  integer, parameter :: mld_ml_type_         = 20
  integer, parameter :: mld_smoother_pos_    = 21
  integer, parameter :: mld_aggr_kind_       = 22
  integer, parameter :: mld_aggr_alg_        = 23
  integer, parameter :: mld_aggr_omega_alg_  = 24
  integer, parameter :: mld_aggr_eig_        = 25
  integer, parameter :: mld_coarse_mat_      = 26
  integer, parameter :: mld_coarse_solve_    = 27 
  integer, parameter :: mld_coarse_sweeps_   = 28
  integer, parameter :: mld_coarse_fillin_   = 29
  integer, parameter :: mld_coarse_subsolve_ = 30
  integer, parameter :: mld_ifpsz_           = 32

  !
  ! Legal values for entry: mld_smoother_type_
  ! 
  integer, parameter :: mld_min_prec_=0, mld_noprec_=0, mld_diag_=1, mld_bjac_=2,&
       & mld_as_=3, mld_max_prec_=3
  !
  ! Legal values for entry: mld_sub_solve_
  !
  integer, parameter :: mld_f_none_=0,mld_ilu_n_=1,mld_milu_n_=2, mld_ilu_t_=3
  integer, parameter :: mld_slu_=4, mld_umf_=5, mld_sludist_=6, mld_max_sub_solve_=mld_sludist_
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
  integer, parameter :: mld_no_smooth_=0, mld_smooth_prol_=1, mld_biz_prol_=2
  ! Disabling biz_prol for the time being.
  integer, parameter :: mld_max_aggr_kind_=mld_smooth_prol_
  !  
  ! Legal values for entry: mld_aggr_alg_
  !
  integer, parameter :: mld_dec_aggr_=0, mld_sym_dec_aggr_=1
  integer, parameter :: mld_glb_aggr_=2, mld_new_dec_aggr_=3
  integer, parameter :: mld_new_glb_aggr_=4
  integer, parameter :: mld_max_aggr_alg_=mld_new_glb_aggr_

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
       &  smooth_names(1:3)=(/'pre-smoothing     ','post-smoothing    ',&
       & 'pre/post-smoothing'/)
  character(len=15), parameter, private :: &
       &  aggr_kinds(0:2)=(/'no  smoother  ','omega smoother',&
       &           'bizr. smoother'/)
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
       &  fact_names(0:6)=(/'none          ','ILU(n)        ',&
       &  'MILU(n)       ','ILU(t,n)      ',&
       &  'SuperLU       ','UMFPACK LU    ',&
       &  'SuperLU_Dist  '/)

  !
  ! Interfaces to routines for checking the definition of the preconditioner,
  ! for printing its description and for deallocating its data structure
  !

  interface mld_base_precfree
    module procedure mld_sbase_precfree, mld_cbase_precfree,&
         &  mld_dbase_precfree, mld_zbase_precfree
  end interface

  interface mld_onelev_precfree
    module procedure mld_s_onelev_precfree, mld_d_onelev_precfree, &
         & mld_c_onelev_precfree, mld_z_onelev_precfree
  end interface

  interface mld_nullify_baseprec
    module procedure mld_nullify_sbaseprec, mld_nullify_cbaseprec,&
         &  mld_nullify_dbaseprec, mld_nullify_zbaseprec
  end interface

  interface mld_nullify_onelevprec
    module procedure  mld_nullify_s_onelevprec, mld_nullify_d_onelevprec,&
         & mld_nullify_c_onelevprec, mld_nullify_z_onelevprec
  end interface

  interface mld_check_def
    module procedure mld_icheck_def, mld_scheck_def, mld_dcheck_def
  end interface

  interface mld_precdescr
    module procedure mld_file_prec_descr, &
         &  mld_zfile_prec_descr,&
         &  mld_sfile_prec_descr,&
         &  mld_cfile_prec_descr
  end interface

  interface mld_prec_short_descr
    module procedure mld_prec_short_descr, mld_zprec_short_descr
  end interface

  interface mld_sizeof
    module procedure mld_sprec_sizeof, mld_cprec_sizeof, &
         & mld_dprec_sizeof, mld_zprec_sizeof, &
         & mld_sbaseprc_sizeof, mld_cbaseprc_sizeof,&
         & mld_dbaseprc_sizeof, mld_zbaseprc_sizeof, &
         & mld_s_onelev_prec_sizeof, mld_d_onelev_prec_sizeof,&
         & mld_c_onelev_prec_sizeof, mld_z_onelev_prec_sizeof
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
    
    info = 0
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
    case('RAW')
      val = mld_no_smooth_
    case('SMOOTH')
      val = mld_smooth_prol_
    case('PRE')
      val = mld_pre_smooth_
    case('POST')
      val = mld_post_smooth_
    case('TWOSIDE')
      val = mld_twoside_smooth_
    case('NOPREC')
      val = mld_noprec_
    case('DIAG')
      val = mld_diag_
    case('BJAC')
      val = mld_bjac_
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
    case default
      val  = -1
      info = -1
    end select
    if (info /= 0) then 
      write(0,*) name,': Error: unknown request: "',trim(string),'"'
    end if
  end subroutine mld_stringval

  !
  ! Function returning the size of the mld_prec_type data structure
  !

  function mld_sprec_sizeof(prec)  result(val)
    implicit none 
    type(mld_sprec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    val = 0
    if (allocated(prec%precv)) then 
      do i=1, size(prec%precv)
        val = val + mld_sizeof(prec%precv(i))
      end do
    end if
  end function mld_sprec_sizeof

  function mld_dprec_sizeof(prec) result(val)
    implicit none 
    type(mld_dprec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    val = 0
    if (allocated(prec%precv)) then 
      do i=1, size(prec%precv)
        val = val + mld_sizeof(prec%precv(i))
      end do
    end if
  end function mld_dprec_sizeof

  function mld_cprec_sizeof(prec) result(val)
    implicit none 
    type(mld_cprec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    val = 0
    if (allocated(prec%precv)) then 
      do i=1, size(prec%precv)
        val = val + mld_sizeof(prec%precv(i))
      end do
    end if
  end function mld_cprec_sizeof

  function mld_zprec_sizeof(prec) result(val)
    implicit none 
    type(mld_zprec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    val = 0
    if (allocated(prec%precv)) then 
      do i=1, size(prec%precv)
        val = val + mld_sizeof(prec%precv(i))
      end do
    end if
  end function mld_zprec_sizeof

  !
  ! Function returning the size of the mld_baseprc_type data structure
  !

  function mld_sbaseprc_sizeof(prec) result(val)
    implicit none 
    type(mld_sbaseprc_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    
    val = 0
    if (allocated(prec%iprcparm)) then 
      val = val + psb_sizeof_int * size(prec%iprcparm)
      if (prec%iprcparm(mld_prec_status_) == mld_prec_built_) then 
        select case(prec%iprcparm(mld_sub_solve_)) 
        case(mld_ilu_n_,mld_ilu_t_)
          ! do nothing
        case(mld_slu_)
        case(mld_umf_)
        case(mld_sludist_)
        case default
        end select
        
      end if
    end if
    if (allocated(prec%rprcparm)) val = val + psb_sizeof_sp * size(prec%rprcparm)
    if (allocated(prec%d))        val = val + psb_sizeof_sp * size(prec%d)
    if (allocated(prec%perm))     val = val + psb_sizeof_int * size(prec%perm)
    if (allocated(prec%invperm))  val = val + psb_sizeof_int * size(prec%invperm)
                                  val = val + psb_sizeof(prec%desc_data)
    if (allocated(prec%av))  then 
      do i=1,size(prec%av)
        val = val + psb_sizeof(prec%av(i))
      end do
    end if
    
  end function mld_sbaseprc_sizeof

  function mld_dbaseprc_sizeof(prec) result(val)
    implicit none 
    type(mld_dbaseprc_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    
    val = 0
    if (allocated(prec%iprcparm)) then 
      val = val + psb_sizeof_int * size(prec%iprcparm)
      if (prec%iprcparm(mld_prec_status_) == mld_prec_built_) then 
        select case(prec%iprcparm(mld_sub_solve_)) 
        case(mld_ilu_n_,mld_ilu_t_)
          ! do nothing
        case(mld_slu_)
        case(mld_umf_)
        case(mld_sludist_)
        case default
        end select
        
      end if
    end if
    if (allocated(prec%rprcparm)) val = val + psb_sizeof_dp * size(prec%rprcparm)
    if (allocated(prec%d))        val = val + psb_sizeof_dp * size(prec%d)
    if (allocated(prec%perm))     val = val + psb_sizeof_int * size(prec%perm)
    if (allocated(prec%invperm))  val = val + psb_sizeof_int * size(prec%invperm)
                                  val = val + psb_sizeof(prec%desc_data)
    if (allocated(prec%av))  then 
      do i=1,size(prec%av)
        val = val + psb_sizeof(prec%av(i))
      end do
    end if


  end function mld_dbaseprc_sizeof

  function mld_cbaseprc_sizeof(prec) result(val)
    implicit none 
    type(mld_cbaseprc_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    
    val = 0
    if (allocated(prec%iprcparm)) then 
      val = val + psb_sizeof_int * size(prec%iprcparm)
      if (prec%iprcparm(mld_prec_status_) == mld_prec_built_) then 
        select case(prec%iprcparm(mld_sub_solve_)) 
        case(mld_ilu_n_,mld_ilu_t_)
          ! do nothing
        case(mld_slu_)
        case(mld_umf_)
        case(mld_sludist_)
        case default
        end select
        
      end if
    end if
    if (allocated(prec%rprcparm)) val = val + psb_sizeof_sp * size(prec%rprcparm)
    if (allocated(prec%d))        val = val + 2 * psb_sizeof_sp * size(prec%d)
    if (allocated(prec%perm))     val = val + psb_sizeof_int * size(prec%perm)
    if (allocated(prec%invperm))  val = val + psb_sizeof_int * size(prec%invperm)
                                  val = val + psb_sizeof(prec%desc_data)
    if (allocated(prec%av))  then 
      do i=1,size(prec%av)
        val = val + psb_sizeof(prec%av(i))
      end do
    end if
    
  end function mld_cbaseprc_sizeof

  function mld_zbaseprc_sizeof(prec) result(val)
    implicit none 
    type(mld_zbaseprc_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    
    val = 0
    if (allocated(prec%iprcparm)) then 
      val = val + psb_sizeof_int * size(prec%iprcparm)
      if (prec%iprcparm(mld_prec_status_) == mld_prec_built_) then 
        select case(prec%iprcparm(mld_sub_solve_)) 
        case(mld_ilu_n_,mld_ilu_t_)
          ! do nothing
        case(mld_slu_)
        case(mld_umf_)
        case(mld_sludist_)
        case default
        end select
        
      end if
    end if
    if (allocated(prec%rprcparm)) val = val + psb_sizeof_dp * size(prec%rprcparm)
    if (allocated(prec%d))        val = val + 2 * psb_sizeof_dp * size(prec%d)
    if (allocated(prec%perm))     val = val + psb_sizeof_int * size(prec%perm)
    if (allocated(prec%invperm))  val = val + psb_sizeof_int * size(prec%invperm)
                                  val = val + psb_sizeof(prec%desc_data)
    if (allocated(prec%av))  then 
      do i=1,size(prec%av)
        val = val + psb_sizeof(prec%av(i))
      end do
    end if
    
  end function mld_zbaseprc_sizeof

  function mld_s_onelev_prec_sizeof(prec) result(val)
    implicit none 
    type(mld_s_onelev_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    
    val = mld_sizeof(prec%prec)
    if (allocated(prec%iprcparm)) then 
      val = val + psb_sizeof_int * size(prec%iprcparm)
    end if
    if (allocated(prec%rprcparm)) val = val + psb_sizeof_sp * size(prec%rprcparm)
    val = val + psb_sizeof(prec%desc_ac)
    val = val + psb_sizeof(prec%ac)
    val = val + psb_sizeof(prec%map_desc) 

  end function mld_s_onelev_prec_sizeof

  function mld_d_onelev_prec_sizeof(prec) result(val)
    implicit none 
    type(mld_d_onelev_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    
    val = mld_sizeof(prec%prec)
    if (allocated(prec%iprcparm)) then 
      val = val + psb_sizeof_int * size(prec%iprcparm)
    end if
    if (allocated(prec%rprcparm)) val = val + psb_sizeof_dp * size(prec%rprcparm)
    val = val + psb_sizeof(prec%desc_ac)
    val = val + psb_sizeof(prec%ac)
    val = val + psb_sizeof(prec%map_desc) 

  end function mld_d_onelev_prec_sizeof

  function mld_c_onelev_prec_sizeof(prec) result(val)
    implicit none 
    type(mld_c_onelev_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    
    val = mld_sizeof(prec%prec)
    if (allocated(prec%iprcparm)) then 
      val = val + psb_sizeof_int * size(prec%iprcparm)
    end if
    if (allocated(prec%rprcparm)) val = val + psb_sizeof_sp * size(prec%rprcparm)
    val = val + psb_sizeof(prec%desc_ac)
    val = val + psb_sizeof(prec%ac)
    val = val + psb_sizeof(prec%map_desc) 

  end function mld_c_onelev_prec_sizeof

  function mld_z_onelev_prec_sizeof(prec) result(val)
    implicit none 
    type(mld_z_onelev_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    
    val = mld_sizeof(prec%prec)
    if (allocated(prec%iprcparm)) then 
      val = val + psb_sizeof_int * size(prec%iprcparm)
    end if
    if (allocated(prec%rprcparm)) val = val + psb_sizeof_dp * size(prec%rprcparm)
    val = val + psb_sizeof(prec%desc_ac)
    val = val + psb_sizeof(prec%ac)
    val = val + psb_sizeof(prec%map_desc) 

  end function mld_z_onelev_prec_sizeof
    
  !
  ! Routines printing out a description of the preconditioner
  !

  subroutine mld_base_prec_descr(iout,iprcparm, info,rprcparm,dprcparm)
    implicit none 
    integer, intent(in) :: iprcparm(:),iout
    integer, intent(out) :: info
    real(psb_spk_), intent(in), optional :: rprcparm(:)
    real(psb_dpk_), intent(in), optional :: dprcparm(:)
    
    info = 0
    if (count((/ present(rprcparm),present(dprcparm) /)) /= 1) then 
      info=581
!!$      call psb_errpush(info,name,a_err=" rprcparm, dprcparm")
      return
    endif
    
    select case(iprcparm(mld_smoother_type_))
    case(mld_noprec_)
      write(iout,*) '  No preconditioning'
    case(mld_diag_)
      write(iout,*) '  Diagonal scaling'
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
      case(mld_slu_,mld_umf_,mld_sludist_) 
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
      case(mld_slu_,mld_umf_,mld_sludist_) 
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

    info = 0
    if (count((/ present(rprcparm),present(dprcparm) /)) /= 1) then 
      info=581
!!$      call psb_errpush(info,name,a_err=" rprcparm, dprcparm")
      return
    endif

    if (iprcparm(mld_ml_type_)>mld_no_ml_) then

      write(iout,*) '  Multilevel type: ',&
           &   ml_names(iprcparm(mld_ml_type_))
      write(iout,*) '  Smoother position: ',&
           & smooth_names(iprcparm(mld_smoother_pos_))
      write(iout,*) '  Aggregation: ', &
           &   aggr_names(iprcparm(mld_aggr_alg_))
      write(iout,*) '  Aggregation smoothing: ', &
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

    info = 0
    if (count((/ present(rprcparm),present(dprcparm) /)) /= 1) then 
      info=581
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

    info = 0
    if (count((/ present(rprcparm),present(dprcparm) /)) /= 1) then 
      info=581
!!$      call psb_errpush(info,name,a_err=" rprcparm, dprcparm")
      return
    endif
    if (count((/ present(rprcparm2),present(dprcparm2) /)) /= 1) then 
      info=581
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
        write(iout,*) '  Coarsest matrix solver: block Jacobi with ', &
             &  fact_names(iprcparm2(mld_sub_solve_))
        write(iout,*) '  Number of Jacobi sweeps: ', &
             &   (iprcparm2(mld_smoother_sweeps_))
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
      case(mld_slu_,mld_umf_,mld_sludist_) 
      case default
        write(iout,*) '  Should never get here!'
      end select
    end if


    return
  end subroutine mld_ml_coarse_descr


  !
  ! Subroutine: mld_file_prec_descr
  ! Version: real
  !
  !  This routine prints a description of the preconditioner to the standard 
  !  output or to a file. It must be called after the preconditioner has been
  !  built by mld_precbld.
  !
  ! Arguments:
  !  p       -  type(mld_dprec_type), input.
  !             The preconditioner data structure to be printed out.
  !  info    -  integer, output.
  !             error code.
  !  iout    -  integer, input, optional.
  !             The id of the file where the preconditioner description
  !             will be printed. If iout is not present, then the standard
  !             output is condidered.
  !
  subroutine mld_file_prec_descr(p,info,iout)
    implicit none 
    ! Arguments
    type(mld_dprec_type), intent(in) :: p
    integer, intent(out)             :: info
    integer, intent(in), optional    :: iout

    ! Local variables
    integer      :: ilev, nlev
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_file_prec_descr'
    integer :: iout_

    info = 0
    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if
    if (iout_ < 0) iout_ = 6 

    if (allocated(p%precv)) then
      ictxt = psb_cd_get_context(p%precv(1)%prec%desc_data)
      
      call psb_info(ictxt,me,np)
      
      !
      ! The preconditioner description is printed by processor psb_root_.
      ! This agrees with the fact that all the parameters defining the
      ! preconditioner have the same values on all the procs (this is
      ! ensured by mld_precbld).
      !
      if (me==psb_root_) then
        
        write(iout_,*) 
        write(iout_,'(a)') 'Preconditioner description'
        nlev = size(p%precv)
        if (nlev >= 1) then
          !
          ! Print description of base preconditioner
          !

          write(iout_,*) ' '

          if (nlev > 1) then
            write(iout_,*) 'Multilevel Schwarz'
            write(iout_,*) 
            write(iout_,*) 'Base preconditioner (smoother) details'
          endif

          ilev = 1 
          call mld_base_prec_descr(iout_,p%precv(ilev)%prec%iprcparm,info,&
               & dprcparm=p%precv(ilev)%prec%rprcparm)

        end if

        if (nlev > 1) then

          !
          ! Print multilevel details
          !
          write(iout_,*) 
          write(iout_,*) 'Multilevel details'

          do ilev = 2, nlev 
            if (.not.allocated(p%precv(ilev)%iprcparm)) then 
              info = 3111
              write(iout_,*) ' ',name,': error: inconsistent MLPREC part, should call MLD_PRECINIT'
              return
            endif
          end do

          write(iout_,*) ' Number of levels: ',nlev

          !
          ! Currently, all the preconditioner parameters must have the same value at levels
          ! 2,...,nlev-1, hence only the values at level 2 are printed
          !

          ilev=2
          call mld_ml_alg_descr(iout_,ilev,p%precv(ilev)%iprcparm, info,&
               & dprcparm=p%precv(ilev)%rprcparm)

          !
          ! Coarse matrices are different at levels 2,...,nlev-1, hence related
          ! info is printed separately
          !
          write(iout_,*) 
          do ilev = 2, nlev-1
            call mld_ml_level_descr(iout_,ilev,p%precv(ilev)%iprcparm,&
                 & p%precv(ilev)%nlaggr,info,&
                 & dprcparm=p%precv(ilev)%rprcparm)
          end do

          !
          ! Print coarsest level details
          !

          ilev = nlev
          write(iout_,*) 
          call mld_ml_coarse_descr(iout_,ilev,&
               & p%precv(ilev)%iprcparm,p%precv(ilev)%prec%iprcparm,&
               & p%precv(ilev)%nlaggr,info,&
               & dprcparm=p%precv(ilev)%rprcparm,&
               & dprcparm2=p%precv(ilev)%prec%rprcparm)
        end if
        
      endif
      write(iout_,*) 
    else
      write(iout_,*) trim(name), &
           & ': Error: no base preconditioner available, something is wrong!'
      info = -2
      return
    endif


  end subroutine mld_file_prec_descr

  subroutine mld_sfile_prec_descr(p,info,iout)
    implicit none 

    ! Arguments
    type(mld_sprec_type), intent(in) :: p
    integer, intent(out)             :: info
    integer, intent(in), optional    :: iout

    ! Local variables
    integer      :: ilev, nlev
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_file_prec_descr'
    integer :: iout_

    info = 0
    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if
    if (iout_ < 0) iout_ = 6 

    if (allocated(p%precv)) then
      ictxt = psb_cd_get_context(p%precv(1)%prec%desc_data)
      
      call psb_info(ictxt,me,np)
      
      !
      ! The preconditioner description is printed by processor psb_root_.
      ! This agrees with the fact that all the parameters defining the
      ! preconditioner have the same values on all the procs (this is
      ! ensured by mld_precbld).
      !
      if (me==psb_root_) then
        
        write(iout_,*) 
        write(iout_,*) 'Preconditioner description'
        nlev = size(p%precv)
        if (nlev >= 1) then
          !
          ! Print description of base preconditioner
          !

          write(iout_,*) ' '

          if (nlev > 1) then
            write(iout_,*) 'Multilevel Schwarz'
            write(iout_,*) 
            write(iout_,*) 'Base preconditioner (smoother) details'
          endif

          ilev = 1 
          call mld_base_prec_descr(iout_,p%precv(ilev)%prec%iprcparm,info,&
               & rprcparm=p%precv(ilev)%prec%rprcparm)

        end if

        if (nlev > 1) then

          !
          ! Print multilevel details
          !
          write(iout_,*) 
          write(iout_,*) 'Multilevel details'

          do ilev = 2, nlev 
            if (.not.allocated(p%precv(ilev)%iprcparm)) then 
              info = 3111
              write(iout_,*) ' ',name,': error: inconsistent MLPREC part, should call MLD_PRECINIT'
              return
            endif
          end do

          write(iout_,*) ' Number of levels: ',nlev

          !
          ! Currently, all the preconditioner parameters must have the same value at levels
          ! 2,...,nlev-1, hence only the values at level 2 are printed
          !

          ilev=2
          call mld_ml_alg_descr(iout_,ilev,p%precv(ilev)%iprcparm, info,&
               & rprcparm=p%precv(ilev)%rprcparm)

          !
          ! Coarse matrices are different at levels 2,...,nlev-1, hence related
          ! info is printed separately
          !
          write(iout_,*)                       
          do ilev = 2, nlev-1
            call mld_ml_level_descr(iout_,ilev,p%precv(ilev)%iprcparm,&
                 & p%precv(ilev)%nlaggr,info,&
                 & rprcparm=p%precv(ilev)%rprcparm)
          end do

          !
          ! Print coarsest level details
          !

          ilev = nlev
          write(iout_,*) 
          call mld_ml_coarse_descr(iout_,ilev,&
               & p%precv(ilev)%iprcparm,p%precv(ilev)%prec%iprcparm,&
               & p%precv(ilev)%nlaggr,info,&
               & rprcparm=p%precv(ilev)%rprcparm,  &
               & rprcparm2=p%precv(ilev)%prec%rprcparm)

        end if
        
      endif
      write(iout_,*) 
    else

      write(iout_,*) trim(name), &
           & ': Error: no base preconditioner available, something is wrong!'
      info = -2
      return
    endif


  end subroutine mld_sfile_prec_descr


  function  mld_prec_short_descr(p)
    implicit none 
    type(mld_dprec_type), intent(in) :: p
    character(len=20) :: mld_prec_short_descr
    mld_prec_short_descr = ' '
  end function mld_prec_short_descr


  !
  ! Subroutine: mld_zfile_prec_descr
  ! Version: complex
  !
  !  This routine prints to a file a description of the preconditioner.
  !
  ! Arguments:
  !  p       -  type(mld_zprec_type), input.
  !             The preconditioner data structure to be printed out.
  !  iout    -  integer, input.
  !             The id of the file where the preconditioner description
  !             will be printed.
  !
  subroutine mld_zfile_prec_descr(p,info,iout)
    implicit none 

    ! Arguments
    type(mld_zprec_type), intent(in) :: p
    integer, intent(out)             :: info
    integer, intent(in), optional    :: iout

    ! Local variables
    integer      :: ilev, nlev
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_file_prec_descr'
    integer :: iout_

    info = 0
    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if
    if (iout_ < 0) iout_ = 6 

    if (allocated(p%precv)) then
      ictxt = psb_cd_get_context(p%precv(1)%prec%desc_data)
      
      call psb_info(ictxt,me,np)
      
      !
      ! The preconditioner description is printed by processor psb_root_.
      ! This agrees with the fact that all the parameters defining the
      ! preconditioner have the same values on all the procs (this is
      ! ensured by mld_precbld).
      !
      if (me==psb_root_) then
        
        write(iout_,*) 
        write(iout_,*) 'Preconditioner description'
        nlev = size(p%precv)
        if (nlev >= 1) then
          !
          ! Print description of base preconditioner
          !

          write(iout_,*) ' '

          if (nlev > 1) then
            write(iout_,*) 'Multilevel Schwarz'
            write(iout_,*) 
            write(iout_,*) 'Base preconditioner (smoother) details'
          endif

          ilev = 1 
          call mld_base_prec_descr(iout_,p%precv(ilev)%prec%iprcparm,info,&
               & dprcparm=p%precv(ilev)%prec%rprcparm)

        end if

        if (nlev > 1) then

          !
          ! Print multilevel details
          !
          write(iout_,*) 
          write(iout_,*) 'Multilevel details'

          do ilev = 2, nlev 
            if (.not.allocated(p%precv(ilev)%iprcparm)) then 
              info = 3111
              write(iout_,*) ' ',name,': error: inconsistent MLPREC part, should call MLD_PRECINIT'
              return
            endif
          end do

          write(iout_,*) ' Number of levels: ',nlev

          !
          ! Currently, all the preconditioner parameters must have the same value at levels
          ! 2,...,nlev-1, hence only the values at level 2 are printed
          !

          ilev=2
          call mld_ml_alg_descr(iout_,ilev,p%precv(ilev)%iprcparm, info,&
               & dprcparm=p%precv(ilev)%rprcparm)

          !
          ! Coarse matrices are different at levels 2,...,nlev-1, hence related
          ! info is printed separately
          !
          write(iout_,*) 
          do ilev = 2, nlev-1
            call mld_ml_level_descr(iout_,ilev,p%precv(ilev)%iprcparm,&
                 & p%precv(ilev)%nlaggr,info,&
                 & dprcparm=p%precv(ilev)%rprcparm)
          end do

          !
          ! Print coarsest level details
          !

          ilev = nlev
          write(iout_,*) 
          call mld_ml_coarse_descr(iout_,ilev,&
               & p%precv(ilev)%iprcparm,p%precv(ilev)%prec%iprcparm,&
               & p%precv(ilev)%nlaggr,info,&
               & dprcparm=p%precv(ilev)%rprcparm,&
               & dprcparm2=p%precv(ilev)%prec%rprcparm)
        end if
        
      endif
      write(iout_,*) 
    else
      write(iout_,*) trim(name), &
           & ': Error: no base preconditioner available, something is wrong!'
      info = -2
      return
    endif

  end subroutine mld_zfile_prec_descr

  subroutine mld_cfile_prec_descr(p,info,iout)
    implicit none 

    ! Arguments
    type(mld_cprec_type), intent(in) :: p
    integer, intent(out)             :: info
    integer, intent(in), optional    :: iout

    ! Local variables
    integer      :: ilev, nlev
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_file_prec_descr'
    integer :: iout_

    info = 0
    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if
    if (iout_ < 0) iout_ = 6 

    if (allocated(p%precv)) then
      ictxt = psb_cd_get_context(p%precv(1)%prec%desc_data)
      
      call psb_info(ictxt,me,np)
      
      !
      ! The preconditioner description is printed by processor psb_root_.
      ! This agrees with the fact that all the parameters defining the
      ! preconditioner have the same values on all the procs (this is
      ! ensured by mld_precbld).
      !
      if (me==psb_root_) then
        write(iout_,*)            
        write(iout_,*) 'Preconditioner description'
        nlev = size(p%precv)
        if (nlev >= 1) then
          !
          ! Print description of base preconditioner
          !

          write(iout_,*) ' '

          if (nlev > 1) then
            write(iout_,*) 'Multilevel Schwarz'
            write(iout_,*) 
            write(iout_,*) 'Base preconditioner (smoother) details'
          endif

          ilev = 1 
          call mld_base_prec_descr(iout_,p%precv(ilev)%prec%iprcparm,info,&
               & rprcparm=p%precv(ilev)%prec%rprcparm)

        end if

        if (nlev > 1) then

          !
          ! Print multilevel details
          !
          write(iout_,*) 
          write(iout_,*) 'Multilevel details'

          do ilev = 2, nlev 
            if (.not.allocated(p%precv(ilev)%iprcparm)) then 
              info = 3111
              write(iout_,*) ' ',name,': error: inconsistent MLPREC part, should call MLD_PRECINIT'
              return
            endif
          end do

          write(iout_,*) ' Number of levels: ',nlev

          !
          ! Currently, all the preconditioner parameters must have the same value at levels
          ! 2,...,nlev-1, hence only the values at level 2 are printed
          !

          ilev=2
          call mld_ml_alg_descr(iout_,ilev,p%precv(ilev)%iprcparm, info,&
               & rprcparm=p%precv(ilev)%rprcparm)

          !
          ! Coarse matrices are different at levels 2,...,nlev-1, hence related
          ! info is printed separately
          !
          write(iout_,*) 

          do ilev = 2, nlev-1
            call mld_ml_level_descr(iout_,ilev,p%precv(ilev)%iprcparm,&
                 & p%precv(ilev)%nlaggr,info,&
                 & rprcparm=p%precv(ilev)%rprcparm)
          end do

          !
          ! Print coarsest level details
          !

          ilev = nlev
          write(iout_,*) 
          call mld_ml_coarse_descr(iout_,ilev,&
               & p%precv(ilev)%iprcparm,p%precv(ilev)%prec%iprcparm,&
               & p%precv(ilev)%nlaggr,info,&
               & rprcparm=p%precv(ilev)%rprcparm,&
               & rprcparm2=p%precv(ilev)%prec%rprcparm)
        end if
        
      endif
      write(iout_,*) 
    else
      write(iout_,*) trim(name), &
           & ': Error: no base preconditioner available, something is wrong!'
      info = -2
      return
    endif

  end subroutine mld_cfile_prec_descr


  function  mld_zprec_short_descr(p)
    implicit none 
    type(mld_zprec_type), intent(in) :: p
    character(len=20) :: mld_zprec_short_descr
    mld_zprec_short_descr = ' '

  end function mld_zprec_short_descr


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
    is_legal_restrict = ((ip==psb_nohalo_).or.(ip==psb_halo_))
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

    is_legal_ml_aggr_omega_alg = ((ip==mld_eig_est_).or.(ip==mld_user_choice_))
    return
  end function is_legal_ml_aggr_omega_alg
  function is_legal_ml_aggr_eig(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_legal_ml_aggr_eig

    is_legal_ml_aggr_eig = (ip==mld_max_norm_)
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
  function is_distr_ml_coarse_mat(ip)
    implicit none 
    integer, intent(in) :: ip
    logical             :: is_distr_ml_coarse_mat

    is_distr_ml_coarse_mat = (ip==mld_distr_mat_)
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
        use psb_base_mod, only : psb_spk_
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
        use psb_base_mod, only : psb_dpk_
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

  subroutine mld_sbase_precfree(p,info)
    implicit none 

    type(mld_sbaseprc_type), intent(inout) :: p
    integer, intent(out)                :: info
    integer :: i

    info = 0

    ! Actually we might just deallocate the top level array, except 
    ! for the inner UMFPACK or SLU stuff

    if (allocated(p%d)) then 
      deallocate(p%d,stat=info)
    end if

    if (allocated(p%av))  then 
      do i=1,size(p%av) 
        call psb_sp_free(p%av(i),info)
        if (info /= 0) then 
          ! Actually, we don't care here about this.
          ! Just let it go.
          ! return
        end if
      enddo
      deallocate(p%av,stat=info)
    end if

    if (allocated(p%desc_data%matrix_data)) &
         & call psb_cdfree(p%desc_data,info)
    
    if (allocated(p%rprcparm)) then 
      deallocate(p%rprcparm,stat=info)
    end if


    if (allocated(p%perm)) then 
      deallocate(p%perm,stat=info)
    endif

    if (allocated(p%invperm)) then 
      deallocate(p%invperm,stat=info)
    endif

    if (allocated(p%iprcparm)) then 
      if (p%iprcparm(mld_sub_solve_)==mld_slu_) then 
!!$        call mld_sslu_free(p%iprcparm(mld_slu_ptr_),info)
      end if
!!$      if (p%iprcparm(mld_sub_solve_)==mld_sludist_) then 
!!$        call mld_ssludist_free(p%iprcparm(mld_slud_ptr_),info)
!!$      end if
!!$      if (p%iprcparm(mld_sub_solve_)==mld_umf_) then 
!!$        call mld_dumf_free(p%iprcparm(mld_umf_symptr_),&
!!$             & p%iprcparm(mld_umf_numptr_),info)
!!$      end if
      deallocate(p%iprcparm,stat=info)
    end if
    call mld_nullify_baseprec(p)
  end subroutine mld_sbase_precfree


  subroutine mld_s_onelev_precfree(p,info)
    implicit none 

    type(mld_s_onelev_prec_type), intent(inout) :: p
    integer, intent(out)                :: info
    integer :: i

    info = 0

    ! Actually we might just deallocate the top level array, except 
    ! for the inner UMFPACK or SLU stuff
    call mld_base_precfree(p%prec,info)
    
    call psb_sp_free(p%ac,info)
    if (allocated(p%desc_ac%matrix_data)) &
         & call psb_cdfree(p%desc_ac,info)
    
    if (allocated(p%rprcparm)) then 
      deallocate(p%rprcparm,stat=info)
    end if
    ! This is a pointer to something else, must not free it here. 
    nullify(p%base_a) 
    ! This is a pointer to something else, must not free it here. 
    nullify(p%base_desc) 

    if (allocated(p%mlia)) then 
      deallocate(p%mlia,stat=info)
    endif

    if (allocated(p%nlaggr)) then 
      deallocate(p%nlaggr,stat=info)
    endif

    !
    ! free explicitly map_desc???
    ! For now thanks to allocatable semantics
    ! works anyway. 
    !

    call mld_nullify_onelevprec(p)
  end subroutine mld_s_onelev_precfree


  subroutine mld_nullify_s_onelevprec(p)
    implicit none 

    type(mld_s_onelev_prec_type), intent(inout) :: p

    nullify(p%base_a) 
    nullify(p%base_desc) 

  end subroutine mld_nullify_s_onelevprec

  subroutine mld_nullify_sbaseprec(p)
    implicit none 

    type(mld_sbaseprc_type), intent(inout) :: p

!!$    nullify(p%base_a) 
!!$    nullify(p%base_desc) 

  end subroutine mld_nullify_sbaseprec


  subroutine mld_dbase_precfree(p,info)
    implicit none 

    type(mld_dbaseprc_type), intent(inout) :: p
    integer, intent(out)                :: info
    integer :: i

    info = 0

    ! Actually we might just deallocate the top level array, except 
    ! for the inner UMFPACK or SLU stuff

    if (allocated(p%d)) then 
      deallocate(p%d,stat=info)
    end if

    if (allocated(p%av))  then 
      do i=1,size(p%av) 
        call psb_sp_free(p%av(i),info)
        if (info /= 0) then 
          ! Actually, we don't care here about this.
          ! Just let it go.
          ! return
        end if
      enddo
      deallocate(p%av,stat=info)
    end if

    if (allocated(p%desc_data%matrix_data)) &
         & call psb_cdfree(p%desc_data,info)
    
    if (allocated(p%rprcparm)) then 
      deallocate(p%rprcparm,stat=info)
    end if

    if (allocated(p%perm)) then 
      deallocate(p%perm,stat=info)
    endif

    if (allocated(p%invperm)) then 
      deallocate(p%invperm,stat=info)
    endif

    if (allocated(p%iprcparm)) then 
      if (p%iprcparm(mld_sub_solve_)==mld_slu_) then 
        call mld_dslu_free(p%iprcparm(mld_slu_ptr_),info)
      end if
      if (p%iprcparm(mld_sub_solve_)==mld_sludist_) then 
        call mld_dsludist_free(p%iprcparm(mld_slud_ptr_),info)
      end if
      if (p%iprcparm(mld_sub_solve_)==mld_umf_) then 
        call mld_dumf_free(p%iprcparm(mld_umf_symptr_),&
             & p%iprcparm(mld_umf_numptr_),info)
      end if
      deallocate(p%iprcparm,stat=info)
    end if
    call mld_nullify_baseprec(p)
  end subroutine mld_dbase_precfree

  subroutine mld_d_onelev_precfree(p,info)
    implicit none 

    type(mld_d_onelev_prec_type), intent(inout) :: p
    integer, intent(out)                :: info
    integer :: i

    info = 0

    ! Actually we might just deallocate the top level array, except 
    ! for the inner UMFPACK or SLU stuff
    call mld_base_precfree(p%prec,info)
    
    call psb_sp_free(p%ac,info)
    if (allocated(p%desc_ac%matrix_data)) &
         & call psb_cdfree(p%desc_ac,info)
    
    if (allocated(p%rprcparm)) then 
      deallocate(p%rprcparm,stat=info)
    end if
    ! This is a pointer to something else, must not free it here. 
    nullify(p%base_a) 
    ! This is a pointer to something else, must not free it here. 
    nullify(p%base_desc) 

    if (allocated(p%mlia)) then 
      deallocate(p%mlia,stat=info)
    endif

    if (allocated(p%nlaggr)) then 
      deallocate(p%nlaggr,stat=info)
    endif

    !
    ! free explicitly map_desc???
    ! For now thanks to allocatable semantics
    ! works anyway. 
    !

    call mld_nullify_onelevprec(p)
  end subroutine mld_d_onelev_precfree

  subroutine mld_nullify_dbaseprec(p)
    implicit none 

    type(mld_dbaseprc_type), intent(inout) :: p
!!$
!!$    nullify(p%base_a) 
!!$    nullify(p%base_desc) 

  end subroutine mld_nullify_dbaseprec

  subroutine mld_nullify_d_onelevprec(p)
    implicit none 

    type(mld_d_onelev_prec_type), intent(inout) :: p

    nullify(p%base_a) 
    nullify(p%base_desc) 

  end subroutine mld_nullify_d_onelevprec

  subroutine mld_cbase_precfree(p,info)
    implicit none 
    type(mld_cbaseprc_type), intent(inout) :: p
    integer, intent(out)                :: info
    integer :: i

    info = 0

    if (allocated(p%d)) then 
      deallocate(p%d,stat=info)
    end if

    if (allocated(p%av))  then 
      do i=1,size(p%av) 
        call psb_sp_free(p%av(i),info)
        if (info /= 0) then 
          ! Actually, we don't care here about this.
          ! Just let it go.
          ! return
        end if
      enddo
      deallocate(p%av,stat=info)

    end if
    if (allocated(p%desc_data%matrix_data)) &
         & call psb_cdfree(p%desc_data,info)
    
    if (allocated(p%rprcparm)) then 
      deallocate(p%rprcparm,stat=info)
    end if

    if (allocated(p%perm)) then 
      deallocate(p%perm,stat=info)
    endif

    if (allocated(p%invperm)) then 
      deallocate(p%invperm,stat=info)
    endif

    if (allocated(p%iprcparm)) then 
      if (p%iprcparm(mld_sub_solve_)==mld_slu_) then 
!!$        call mld_cslu_free(p%iprcparm(mld_slu_ptr_),info)
      end if
!!$      if (p%iprcparm(mld_sub_solve_)==mld_umf_) then 
!!$        call mld_zumf_free(p%iprcparm(mld_umf_symptr_),&
!!$             & p%iprcparm(mld_umf_numptr_),info)
!!$      end if
      deallocate(p%iprcparm,stat=info)
    end if
    call mld_nullify_baseprec(p)
  end subroutine mld_cbase_precfree

  subroutine mld_c_onelev_precfree(p,info)
    implicit none 

    type(mld_c_onelev_prec_type), intent(inout) :: p
    integer, intent(out)                :: info
    integer :: i

    info = 0

    ! Actually we might just deallocate the top level array, except 
    ! for the inner UMFPACK or SLU stuff
    call mld_base_precfree(p%prec,info)
    
    call psb_sp_free(p%ac,info)
    if (allocated(p%desc_ac%matrix_data)) &
         & call psb_cdfree(p%desc_ac,info)
    
    if (allocated(p%rprcparm)) then 
      deallocate(p%rprcparm,stat=info)
    end if
    ! This is a pointer to something else, must not free it here. 
    nullify(p%base_a) 
    ! This is a pointer to something else, must not free it here. 
    nullify(p%base_desc) 

    if (allocated(p%mlia)) then 
      deallocate(p%mlia,stat=info)
    endif

    if (allocated(p%nlaggr)) then 
      deallocate(p%nlaggr,stat=info)
    endif

    !
    ! free explicitly map_desc???
    ! For now thanks to allocatable semantics
    ! works anyway. 
    !

    call mld_nullify_onelevprec(p)
  end subroutine mld_c_onelev_precfree

  subroutine mld_nullify_c_onelevprec(p)
    implicit none 

    type(mld_c_onelev_prec_type), intent(inout) :: p

    nullify(p%base_a) 
    nullify(p%base_desc) 

  end subroutine mld_nullify_c_onelevprec

  subroutine mld_nullify_cbaseprec(p)
    implicit none 

    type(mld_cbaseprc_type), intent(inout) :: p

!!$    nullify(p%base_a) 
!!$    nullify(p%base_desc) 

  end subroutine mld_nullify_cbaseprec

  subroutine mld_zbase_precfree(p,info)
    implicit none 
    type(mld_zbaseprc_type), intent(inout) :: p
    integer, intent(out)                :: info
    integer :: i

    info = 0

    if (allocated(p%d)) then 
      deallocate(p%d,stat=info)
    end if

    if (allocated(p%av))  then 
      do i=1,size(p%av) 
        call psb_sp_free(p%av(i),info)
        if (info /= 0) then 
          ! Actually, we don't care here about this.
          ! Just let it go.
          ! return
        end if
      enddo
      deallocate(p%av,stat=info)

    end if
    if (allocated(p%desc_data%matrix_data)) &
         & call psb_cdfree(p%desc_data,info)
    
    if (allocated(p%rprcparm)) then 
      deallocate(p%rprcparm,stat=info)
    end if

    if (allocated(p%perm)) then 
      deallocate(p%perm,stat=info)
    endif

    if (allocated(p%invperm)) then 
      deallocate(p%invperm,stat=info)
    endif

    if (allocated(p%iprcparm)) then 
      if (p%iprcparm(mld_sub_solve_)==mld_slu_) then 
        call mld_zslu_free(p%iprcparm(mld_slu_ptr_),info)
      end if
      if (p%iprcparm(mld_sub_solve_)==mld_umf_) then 
        call mld_zumf_free(p%iprcparm(mld_umf_symptr_),&
             & p%iprcparm(mld_umf_numptr_),info)
      end if
      deallocate(p%iprcparm,stat=info)
    end if
    call mld_nullify_baseprec(p)
  end subroutine mld_zbase_precfree

  subroutine mld_z_onelev_precfree(p,info)
    implicit none 

    type(mld_z_onelev_prec_type), intent(inout) :: p
    integer, intent(out)                :: info
    integer :: i

    info = 0

    ! Actually we might just deallocate the top level array, except 
    ! for the inner UMFPACK or SLU stuff
    call mld_base_precfree(p%prec,info)
    
    call psb_sp_free(p%ac,info)
    if (allocated(p%desc_ac%matrix_data)) &
         & call psb_cdfree(p%desc_ac,info)
    
    if (allocated(p%rprcparm)) then 
      deallocate(p%rprcparm,stat=info)
    end if
    ! This is a pointer to something else, must not free it here. 
    nullify(p%base_a) 
    ! This is a pointer to something else, must not free it here. 
    nullify(p%base_desc) 

    if (allocated(p%mlia)) then 
      deallocate(p%mlia,stat=info)
    endif

    if (allocated(p%nlaggr)) then 
      deallocate(p%nlaggr,stat=info)
    endif

    !
    ! free explicitly map_desc???
    ! For now thanks to allocatable semantics
    ! works anyway. 
    !

    call mld_nullify_onelevprec(p)
  end subroutine mld_z_onelev_precfree

  subroutine mld_nullify_z_onelevprec(p)
    implicit none 

    type(mld_z_onelev_prec_type), intent(inout) :: p

    nullify(p%base_a) 
    nullify(p%base_desc) 

  end subroutine mld_nullify_z_onelevprec


  subroutine mld_nullify_zbaseprec(p)
    implicit none 

    type(mld_zbaseprc_type), intent(inout) :: p

!!$    nullify(p%base_a) 
!!$    nullify(p%base_desc) 

  end subroutine mld_nullify_zbaseprec


  function pr_to_str(iprec)
    implicit none 

    integer, intent(in)  :: iprec
    character(len=10)     :: pr_to_str

    select case(iprec)
    case(mld_noprec_)
      pr_to_str='NOPREC'
    case(mld_diag_)         
      pr_to_str='DIAG'
    case(mld_bjac_)         
      pr_to_str='BJAC'
    case(mld_as_)      
      pr_to_str='AS'
    end select

  end function pr_to_str

end module mld_prec_type
