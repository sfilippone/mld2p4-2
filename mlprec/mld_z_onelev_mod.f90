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
! File: mld_z_onelev_mod.f90
!
! Module: mld_z_onelev_mod
!
!  This module defines: 
!  - the mld_z_onelev_type data structure containing one level
!    of a multilevel  preconditioner and related
!    data structures;
!
!  It contains routines for
!  - Building and applying; 
!  - checking if the preconditioner is correctly defined;
!  - printing a	description of the preconditioner;
!  - deallocating the preconditioner data structure.  
!

module mld_z_onelev_mod

  use mld_base_prec_type
  use mld_z_base_smoother_mod
  use mld_z_dec_aggregator_mod
  use psb_base_mod, only : psb_zspmat_type, psb_z_vect_type, &
       & psb_z_base_vect_type, psb_zlinmap_type, psb_dpk_, &
       & psb_ipk_, psb_long_int_k_, psb_desc_type, psb_i_base_vect_type, &
       & psb_erractionsave, psb_error_handler
  !
  !
  ! Type: mld_zonelev_type.
  !
  !  It is the data type containing the necessary items for the	current
  !  level (essentially, the smoother, the current-level matrix
  !  and the restriction and prolongation operators).
  !
  !  type mld_zonelev_type
  !    class(mld_z_base_smoother_type), allocatable   :: sm, sm2a
  !    class(mld_z_base_smoother_type), pointer       :: sm2 => null()
  !    class(mld_zmlprec_wrk_type), allocatable       :: wrk
  !    class(mld_z_base_aggregator_type), allocatable :: aggr
  !    type(mld_dml_parms)             :: parms 
  !    type(psb_zspmat_type)           :: ac
  !    type(psb_zesc_type)             :: desc_ac
  !    type(psb_zspmat_type), pointer  :: base_a    => null() 
  !    type(psb_desc_type), pointer    :: base_desc => null() 
  !    type(psb_zlinmap_type)          :: map
  !  end type mld_zonelev_type
  !
  !  Note that d denotes the kind of the real data type to be chosen
  !  according to single/double precision version of MLD2P4.
  !
  !   sm,sm2a      -  class(mld_z_base_smoother_type), allocatable
  !                   The current level pre- and post-smooother.
  !   sm2          -  class(mld_z_base_smoother_type), pointer
  !                   The current level post-smooother; if sm2a is allocated
  !                   explicitly, then sm2 => sm2a, otherwise sm2 => sm.
  !   wrk          -  class(mld_zmlprec_wrk_type), allocatable
  !                   Workspace for application of preconditioner; may be
  !                   pre-allocated to save time in the application within a
  !                   Krylov solver.
  !   aggr         -  class(mld_z_base_aggregator_type), allocatable 
  !                   The aggregator object: holds the algorithmic choices and
  !                   (possibly) additional data for building the aggregation.
  !   parms        -  type(mld_dml_parms)
  !                   The parameters defining the multilevel strategy.
  !   ac           -  The local part of the current-level matrix, built by
  !                   coarsening the previous-level matrix.
  !   desc_ac      -  type(psb_desc_type).
  !                   The communication descriptor associated to the matrix
  !                   stored in ac.
  !   base_a       -  type(psb_zspmat_type), pointer.
  !                   Pointer (really a pointer!) to the local part of the current 
  !                   matrix (so we have a unified treatment of residuals).
  !                   We need this to avoid passing explicitly the current matrix
  !                   to the routine which applies the preconditioner.
  !   base_desc    -  type(psb_desc_type), pointer.
  !                   Pointer to the communication descriptor associated to the
  !                   matrix pointed by base_a.
  !   map          -  Stores the maps (restriction and prolongation) between the
  !                   vector spaces associated to the index spaces of the previous
  !                   and current levels.
  !
  !   Methods:  
  !     Most methods follow the encapsulation hierarchy: they take whatever action
  !     is appropriate for the current object, then call the corresponding method for
  !     the contained object.
  !     As an example: the descr() method prints out a description of the
  !     level. It starts by invoking the descr() method of the parms object,
  !     then calls the descr() method of the smoother object. 
  !
  !    descr      -   Prints a description of the object.
  !    default    -   Set default values
  !    dump       -   Dump to file object contents
  !    set        -   Sets various parameters; when a request is unknown
  !                   it is passed to the smoother object for further processing.
  !    check        -  Sanity checks.
  !    sizeof       -  Total memory occupation in bytes
  !    get_nzeros   -  Number of nonzeros 
  !    get_wrksz    -  How many workspace vector does apply_vect need
  !    allocate_wrk -  Allocate auxiliary workspace
  !    free_wrk     -  Free     auxiliary workspace
  !    bld_tprol    -  Invoke the aggr method to build the tentative prolongator
  !    mat_asb      -  Build the final (possibly smoothed) prolongator and coarse matrix. 
  !
  !  
  type mld_zmlprec_wrk_type
    complex(psb_dpk_), allocatable  :: tx(:), ty(:), x2l(:), y2l(:)
    type(psb_z_vect_type)  :: vtx, vty, vx2l, vy2l
    integer(psb_ipk_)        :: wvsz = 0
    type(psb_z_vect_type), allocatable :: wv(:)
  contains
    procedure, pass(wk) :: alloc      => z_wrk_alloc
    procedure, pass(wk) :: free       => z_wrk_free
    procedure, pass(wk) :: clone      => z_wrk_clone
    procedure, pass(wk) :: move_alloc => z_wrk_move_alloc
    procedure, pass(wk) :: cnv        => z_wrk_cnv
    procedure, pass(wk) :: sizeof     => z_wrk_sizeof    
  end type mld_zmlprec_wrk_type
  private :: z_wrk_alloc, z_wrk_free, &
       & z_wrk_clone, z_wrk_move_alloc, z_wrk_cnv, z_wrk_sizeof    
  
  type mld_z_onelev_type
    class(mld_z_base_smoother_type), allocatable   :: sm, sm2a
    class(mld_z_base_smoother_type), pointer       :: sm2 => null()
    class(mld_zmlprec_wrk_type), allocatable       :: wrk
    class(mld_z_base_aggregator_type), allocatable :: aggr
    type(mld_dml_parms)              :: parms 
    type(psb_zspmat_type)            :: ac
    integer(psb_ipk_)                :: ac_nz_loc, ac_nz_tot
    type(psb_desc_type)              :: desc_ac
    type(psb_zspmat_type), pointer   :: base_a    => null() 
    type(psb_desc_type), pointer     :: base_desc => null() 
    type(psb_zspmat_type)            :: tprol
    type(psb_zlinmap_type)           :: map
    real(psb_dpk_)                     :: szratio
  contains
    procedure, pass(lv) :: bld_tprol   => z_base_onelev_bld_tprol
    procedure, pass(lv) :: mat_asb     => mld_z_base_onelev_mat_asb
    procedure, pass(lv) :: backfix     => z_base_onelev_backfix
    procedure, pass(lv) :: update_aggr => z_base_onelev_update_aggr
    procedure, pass(lv) :: bld     => mld_z_base_onelev_build
    procedure, pass(lv) :: clone   => z_base_onelev_clone
    procedure, pass(lv) :: cnv     => mld_z_base_onelev_cnv
    procedure, pass(lv) :: descr   => mld_z_base_onelev_descr
    procedure, pass(lv) :: default => z_base_onelev_default
    procedure, pass(lv) :: free    => mld_z_base_onelev_free
    procedure, pass(lv) :: nullify => z_base_onelev_nullify
    procedure, pass(lv) :: check => mld_z_base_onelev_check
    procedure, pass(lv) :: dump  => mld_z_base_onelev_dump
    procedure, pass(lv) :: cseti => mld_z_base_onelev_cseti
    procedure, pass(lv) :: csetr => mld_z_base_onelev_csetr
    procedure, pass(lv) :: csetc => mld_z_base_onelev_csetc
    procedure, pass(lv) :: setsm => mld_z_base_onelev_setsm
    procedure, pass(lv) :: setsv => mld_z_base_onelev_setsv
    procedure, pass(lv) :: setag => mld_z_base_onelev_setag
    generic, public     :: set   => cseti, csetr, csetc, setsm, setsv, setag 
    procedure, pass(lv) :: sizeof => z_base_onelev_sizeof
    procedure, pass(lv) :: get_nzeros => z_base_onelev_get_nzeros
    procedure, pass(lv) :: get_wrksz => z_base_onelev_get_wrksize
    procedure, pass(lv) :: allocate_wrk   => z_base_onelev_allocate_wrk
    procedure, pass(lv) :: free_wrk       => z_base_onelev_free_wrk
    procedure, nopass   :: stringval => mld_stringval
    procedure, pass(lv) :: move_alloc => z_base_onelev_move_alloc
    
  end type mld_z_onelev_type

  type mld_z_onelev_node
    type(mld_z_onelev_type) :: item
    type(mld_z_onelev_node), pointer :: prev=>null(), next=>null()
  end type mld_z_onelev_node

  private :: z_base_onelev_default, z_base_onelev_sizeof, &
       & z_base_onelev_nullify, z_base_onelev_get_nzeros, &
       & z_base_onelev_clone, z_base_onelev_move_alloc, &
       & z_base_onelev_get_wrksize, z_base_onelev_allocate_wrk, &
       & z_base_onelev_free_wrk

  interface 
    subroutine mld_z_base_onelev_mat_asb(lv,a,desc_a,ilaggr,nlaggr,op_prol,info)
      import :: psb_zspmat_type, psb_desc_type, psb_dpk_, psb_ipk_
      import :: mld_z_onelev_type
      implicit none 
      class(mld_z_onelev_type), intent(inout), target :: lv
      type(psb_zspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)     :: desc_a
      integer(psb_ipk_), intent(inout) :: ilaggr(:),nlaggr(:)
      type(psb_zspmat_type), intent(inout)  :: op_prol
      integer(psb_ipk_), intent(out)      :: info
    end subroutine mld_z_base_onelev_mat_asb
  end interface

  interface
    subroutine mld_z_base_onelev_build(lv,info,amold,vmold,imold)
      import :: psb_z_base_sparse_mat, psb_z_base_vect_type, &
           & psb_i_base_vect_type, psb_dpk_, mld_z_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      implicit none
      class(mld_z_onelev_type), target, intent(inout) :: lv
      integer(psb_ipk_), intent(out) :: info
      class(psb_z_base_sparse_mat), intent(in), optional :: amold
      class(psb_z_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine mld_z_base_onelev_build
  end interface

  interface 
    subroutine mld_z_base_onelev_descr(lv,il,nl,ilmin,info,iout)
      import :: psb_zspmat_type, psb_z_vect_type, psb_z_base_vect_type, &
           & psb_zlinmap_type, psb_dpk_, mld_z_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      ! Arguments
      class(mld_z_onelev_type), intent(in) :: lv
      integer(psb_ipk_), intent(in)                 :: il,nl,ilmin
      integer(psb_ipk_), intent(out)                :: info
      integer(psb_ipk_), intent(in), optional       :: iout
    end subroutine mld_z_base_onelev_descr
  end interface

  interface 
    subroutine mld_z_base_onelev_cnv(lv,info,amold,vmold,imold)
      import :: mld_z_onelev_type, psb_z_base_vect_type, psb_dpk_, &
           & psb_z_base_sparse_mat, psb_ipk_, psb_i_base_vect_type
      ! Arguments
      class(mld_z_onelev_type), intent(inout)            :: lv 
      integer(psb_ipk_), intent(out)                     :: info
      class(psb_z_base_sparse_mat), intent(in), optional :: amold
      class(psb_z_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine mld_z_base_onelev_cnv
  end interface
   
  interface 
    subroutine mld_z_base_onelev_free(lv,info)
      import :: psb_zspmat_type, psb_z_vect_type, psb_z_base_vect_type, &
           & psb_zlinmap_type, psb_dpk_, mld_z_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      implicit none 
      
      class(mld_z_onelev_type), intent(inout) :: lv
      integer(psb_ipk_), intent(out)                :: info
    end subroutine mld_z_base_onelev_free
  end interface
  
  interface 
    subroutine mld_z_base_onelev_check(lv,info)
      import :: psb_zspmat_type, psb_z_vect_type, psb_z_base_vect_type, &
           & psb_zlinmap_type, psb_dpk_, mld_z_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      ! Arguments
      class(mld_z_onelev_type), intent(inout) :: lv 
      integer(psb_ipk_), intent(out)            :: info
    end subroutine mld_z_base_onelev_check
  end interface
  
  interface 
    subroutine mld_z_base_onelev_setsm(lv,val,info,pos)
      import :: psb_dpk_, mld_z_onelev_type, mld_z_base_smoother_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      
      ! Arguments
      class(mld_z_onelev_type), target, intent(inout) :: lv 
      class(mld_z_base_smoother_type), intent(in)     :: val
      integer(psb_ipk_), intent(out)                  :: info
      character(len=*), optional, intent(in)          :: pos
    end subroutine mld_z_base_onelev_setsm
  end interface
  
  interface 
    subroutine mld_z_base_onelev_setsv(lv,val,info,pos)
      import :: psb_dpk_, mld_z_onelev_type, mld_z_base_solver_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      
      ! Arguments
      class(mld_z_onelev_type), target, intent(inout) :: lv 
      class(mld_z_base_solver_type), intent(in)       :: val
      integer(psb_ipk_), intent(out)                  :: info
      character(len=*), optional, intent(in)          :: pos
    end subroutine mld_z_base_onelev_setsv
  end interface
  
  interface 
    subroutine mld_z_base_onelev_setag(lv,val,info,pos)
      import :: psb_dpk_, mld_z_onelev_type, mld_z_base_aggregator_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      
      ! Arguments
      class(mld_z_onelev_type), target, intent(inout) :: lv 
      class(mld_z_base_aggregator_type), intent(in)       :: val
      integer(psb_ipk_), intent(out)                  :: info
      character(len=*), optional, intent(in)          :: pos
    end subroutine mld_z_base_onelev_setag
  end interface
  
  interface 
    subroutine mld_z_base_onelev_cseti(lv,what,val,info,pos,idx)
      import :: psb_zspmat_type, psb_z_vect_type, psb_z_base_vect_type, &
           & psb_zlinmap_type, psb_dpk_, mld_z_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      
      ! Arguments
      class(mld_z_onelev_type), intent(inout) :: lv 
      character(len=*), intent(in)              :: what 
      integer(psb_ipk_), intent(in)             :: val
      integer(psb_ipk_), intent(out)            :: info
      character(len=*), optional, intent(in)      :: pos
    integer(psb_ipk_), intent(in), optional       :: idx
    end subroutine mld_z_base_onelev_cseti
  end interface
  
  interface 
    subroutine mld_z_base_onelev_csetc(lv,what,val,info,pos,idx)
      import :: psb_zspmat_type, psb_z_vect_type, psb_z_base_vect_type, &
           & psb_zlinmap_type, psb_dpk_, mld_z_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      ! Arguments
      class(mld_z_onelev_type), intent(inout) :: lv 
      character(len=*), intent(in)              :: what 
      character(len=*), intent(in)              :: val
      integer(psb_ipk_), intent(out)            :: info
      character(len=*), optional, intent(in)      :: pos
    integer(psb_ipk_), intent(in), optional       :: idx
    end subroutine mld_z_base_onelev_csetc
  end interface
  
  interface 
    subroutine mld_z_base_onelev_csetr(lv,what,val,info,pos,idx)
      import :: psb_zspmat_type, psb_z_vect_type, psb_z_base_vect_type, &
           & psb_zlinmap_type, psb_dpk_, mld_z_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      
      class(mld_z_onelev_type), intent(inout) :: lv 
      character(len=*), intent(in)              :: what 
      real(psb_dpk_), intent(in)                 :: val
      integer(psb_ipk_), intent(out)            :: info
      character(len=*), optional, intent(in)      :: pos
      integer(psb_ipk_), intent(in), optional     :: idx
    end subroutine mld_z_base_onelev_csetr
  end interface

  interface 
    subroutine mld_z_base_onelev_dump(lv,level,info,prefix,head,ac,rp,smoother,&
         & solver,tprol,global_num)
      import :: psb_zspmat_type, psb_z_vect_type, psb_z_base_vect_type, &
           & psb_zlinmap_type, psb_dpk_, mld_z_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      implicit none 
      class(mld_z_onelev_type), intent(in) :: lv
      integer(psb_ipk_), intent(in)          :: level
      integer(psb_ipk_), intent(out)         :: info
      character(len=*), intent(in), optional :: prefix, head
      logical, optional, intent(in)    :: ac, rp, smoother, solver, tprol, global_num
    end subroutine mld_z_base_onelev_dump
  end interface
  
contains
  !
  ! Function returning the size of the mld_prec_type data structure
  ! in bytes or in number of nonzeros of the operator(s) involved. 
  !

  function z_base_onelev_get_nzeros(lv) result(val)
    implicit none 
    class(mld_z_onelev_type), intent(in) :: lv
    integer(psb_long_int_k_) :: val
    integer(psb_ipk_)        :: i
    val = 0
    if (allocated(lv%sm)) &
         &  val =  lv%sm%get_nzeros()
    if (allocated(lv%sm2a)) &
         &  val =  val + lv%sm2a%get_nzeros()
  end function z_base_onelev_get_nzeros

  function z_base_onelev_sizeof(lv) result(val)
    implicit none 
    class(mld_z_onelev_type), intent(in) :: lv
    integer(psb_long_int_k_) :: val
    integer(psb_ipk_)        :: i
    
    val = 0
    val = val + lv%desc_ac%sizeof()
    val = val + lv%ac%sizeof()
    val = val + lv%tprol%sizeof()
    val = val + lv%map%sizeof() 
    if (allocated(lv%sm))   val = val + lv%sm%sizeof()
    if (allocated(lv%sm2a)) val = val + lv%sm2a%sizeof()
    if (allocated(lv%aggr)) val = val + lv%aggr%sizeof()
    if (allocated(lv%wrk))  val = val + lv%wrk%sizeof()
  end function z_base_onelev_sizeof


  subroutine z_base_onelev_nullify(lv)
    implicit none 

    class(mld_z_onelev_type), intent(inout) :: lv

    nullify(lv%base_a) 
    nullify(lv%base_desc) 
    nullify(lv%sm2)
  end subroutine z_base_onelev_nullify

  !
  ! Multilevel defaults: 
  !  multiplicative vs. additive ML framework;
  !  Smoothed decoupled aggregation with zero threshold; 
  !  distributed coarse matrix;
  !  damping omega  computed with the max-norm estimate of the
  !  dominant eigenvalue;
  !  two-sided smoothing (i.e. V-cycle) with 1 smoothing sweep;
  !

  subroutine z_base_onelev_default(lv)

    Implicit None
 
    ! Arguments
    class(mld_z_onelev_type), target, intent(inout) :: lv
    integer(psb_ipk_) :: info 

    lv%parms%sweeps_pre      = 1
    lv%parms%sweeps_post     = 1
    lv%parms%ml_cycle        = mld_vcycle_ml_
    lv%parms%aggr_type       = mld_soc1_
    lv%parms%par_aggr_alg    = mld_dec_aggr_
    lv%parms%aggr_ord        = mld_aggr_ord_nat_
    lv%parms%aggr_prol       = mld_smooth_prol_
    lv%parms%coarse_mat      = mld_distr_mat_
    lv%parms%aggr_omega_alg  = mld_eig_est_
    lv%parms%aggr_eig        = mld_max_norm_
    lv%parms%aggr_filter     = mld_no_filter_mat_
    lv%parms%aggr_omega_val  = dzero
    lv%parms%aggr_thresh     = 0.01_psb_dpk_
    
    if (allocated(lv%sm)) call lv%sm%default()
    if (allocated(lv%sm2a)) then
      call lv%sm2a%default()
      lv%sm2 => lv%sm2a
    else
      lv%sm2 => lv%sm
    end if
    if (.not.allocated(lv%aggr)) allocate(mld_z_dec_aggregator_type :: lv%aggr,stat=info)
    if (allocated(lv%aggr)) call lv%aggr%default()
    
    return

  end subroutine z_base_onelev_default

  subroutine  z_base_onelev_bld_tprol(lv,a,desc_a,&
       & ilaggr,nlaggr,op_prol,ag_data,info)
    implicit none
    class(mld_z_onelev_type), intent(inout), target :: lv
    type(psb_zspmat_type), intent(inout)   :: a
    type(psb_desc_type), intent(inout)       :: desc_a
    integer(psb_ipk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
    type(psb_zspmat_type), intent(out)  :: op_prol
    type(mld_daggr_data), intent(in)    :: ag_data
    integer(psb_ipk_), intent(out)      :: info
    
    call lv%aggr%bld_tprol(lv%parms,ag_data,a,desc_a,ilaggr,nlaggr,op_prol,info)
    
  end subroutine z_base_onelev_bld_tprol


  subroutine  z_base_onelev_update_aggr(lv,lvnext,info)
    implicit none
    class(mld_z_onelev_type), intent(inout), target :: lv, lvnext
    integer(psb_ipk_), intent(out)      :: info

    call lv%aggr%update_next(lvnext%aggr,info)
    
  end subroutine z_base_onelev_update_aggr


  subroutine  z_base_onelev_backfix(lv,lvprev,info)
    implicit none
    class(mld_z_onelev_type), intent(inout), target :: lv, lvprev
    integer(psb_ipk_), intent(out)      :: info

    info = psb_success_
    if (lv%aggr%xt_desc()) then
      call lv%aggr%backfix(lvprev%base_a,lvprev%ac,&
           & lvprev%base_desc,lvprev%desc_ac,info)
    end if
    
  end subroutine z_base_onelev_backfix


  subroutine z_base_onelev_clone(lv,lvout,info)

    Implicit None

    ! Arguments
    class(mld_z_onelev_type), target, intent(inout) :: lv 
    class(mld_z_onelev_type), target, intent(inout) :: lvout
    integer(psb_ipk_), intent(out)                    :: info 

    info = psb_success_
    if (allocated(lv%sm)) then 
      call lv%sm%clone(lvout%sm,info)
    else 
      if (allocated(lvout%sm)) then 
        call lvout%sm%free(info)
        if (info==psb_success_) deallocate(lvout%sm,stat=info)
      end if
    end if
    if (allocated(lv%sm2a)) then 
      call lv%sm%clone(lvout%sm2a,info)
      lvout%sm2 => lvout%sm2a
    else 
      if (allocated(lvout%sm2a)) then 
        call lvout%sm2a%free(info)
        if (info==psb_success_) deallocate(lvout%sm2a,stat=info)
      end if
      lvout%sm2 => lvout%sm
    end if
    if (allocated(lv%aggr)) then 
      call  lv%aggr%clone(lvout%aggr,info)
    else
      if (allocated(lvout%aggr)) then 
        call lvout%aggr%free(info)
        if (info==psb_success_) deallocate(lvout%aggr,stat=info)
      end if
    end if
    if (info == psb_success_) call lv%parms%clone(lvout%parms,info)
    if (info == psb_success_) call lv%ac%clone(lvout%ac,info)
    if (info == psb_success_) call lv%tprol%clone(lvout%tprol,info)
    if (info == psb_success_) call lv%desc_ac%clone(lvout%desc_ac,info)
    if (info == psb_success_) call lv%map%clone(lvout%map,info)
    lvout%base_a    => lv%base_a
    lvout%base_desc => lv%base_desc
    
    return

  end subroutine z_base_onelev_clone

  subroutine z_base_onelev_move_alloc(lv, b,info)
    use psb_base_mod
    implicit none
    class(mld_z_onelev_type), target, intent(inout) :: lv, b
    integer(psb_ipk_), intent(out) :: info 
    
    call b%free(info)
    b%parms  = lv%parms
    b%szratio = lv%szratio
    if (associated(lv%sm2,lv%sm2a)) then 
      call move_alloc(lv%sm,b%sm)
      call move_alloc(lv%sm2a,b%sm2a)
      b%sm2 =>b%sm2a
    else
      call move_alloc(lv%sm,b%sm)
      call move_alloc(lv%sm2a,b%sm2a)
      b%sm2 =>b%sm
    end if

    call move_alloc(lv%aggr,b%aggr)
    if (info == psb_success_) call psb_move_alloc(lv%ac,b%ac,info) 
    if (info == psb_success_) call psb_move_alloc(lv%tprol,b%tprol,info) 
    if (info == psb_success_) call psb_move_alloc(lv%desc_ac,b%desc_ac,info) 
    if (info == psb_success_) call psb_move_alloc(lv%map,b%map,info) 
    b%base_a    => lv%base_a
    b%base_desc => lv%base_desc
    
  end subroutine z_base_onelev_move_alloc

  
  function z_base_onelev_get_wrksize(lv) result(val)
    implicit none 
    class(mld_z_onelev_type), intent(inout) :: lv
    integer(psb_ipk_)  :: val

    val = 0
    ! SM and SM2A can share work vectors
    if (allocated(lv%sm))   val = val + lv%sm%get_wrksz()
    if (allocated(lv%sm2a)) val = max(val,lv%sm2a%get_wrksz())
    !
    ! Now for the ML application itself
    !

    !  VTX/VTY/VX2L/VY2L are stored explicitly
    !

    !
    ! additions for specific ML/cycles
    !
    select case(lv%parms%ml_cycle)
    case(mld_add_ml_,mld_mult_ml_,mld_vcycle_ml_, mld_wcycle_ml_)
      ! We're good
      
    case(mld_kcycle_ml_, mld_kcyclesym_ml_)
      !
      ! We need 7 in inneritkcycle.
      !  Can we reuse vtx? 
      ! 
      val = val + 7
      
    case default
      ! Need a better error signaling ?
      val = -1
    end select
    
  end function z_base_onelev_get_wrksize

  subroutine z_base_onelev_allocate_wrk(lv,info,vmold)
    use psb_base_mod
    implicit none
    class(mld_z_onelev_type), target, intent(inout) :: lv
    integer(psb_ipk_), intent(out) :: info 
    class(psb_z_base_vect_type), intent(in), optional  :: vmold
    !
    integer(psb_ipk_) :: nwv, i
    info = psb_success_
    nwv = lv%get_wrksz()
    if (.not.allocated(lv%wrk)) allocate(lv%wrk,stat=info)
    if (info == 0) call lv%wrk%alloc(nwv,lv%base_desc,info,vmold=vmold)
    
  end subroutine z_base_onelev_allocate_wrk

  
  subroutine z_base_onelev_free_wrk(lv,info)
    use psb_base_mod
    implicit none
    class(mld_z_onelev_type), target, intent(inout) :: lv
    integer(psb_ipk_), intent(out) :: info 
    !
    integer(psb_ipk_) :: nwv,i 
    info = psb_success_

    if (allocated(lv%wrk)) then
      call lv%wrk%free(info)
      if (info == 0) deallocate(lv%wrk,stat=info)
    end if
  end subroutine z_base_onelev_free_wrk
  
  subroutine z_wrk_alloc(wk,nwv,desc,info,vmold)
    use psb_base_mod
    
    Implicit None
    
    ! Arguments
    class(mld_zmlprec_wrk_type), target, intent(inout) :: wk
    integer(psb_ipk_), intent(in)                     :: nwv
    type(psb_desc_type), intent(in)                   :: desc
    integer(psb_ipk_), intent(out)                    :: info 
    class(psb_z_base_vect_type), intent(in), optional  :: vmold
    !
    integer(psb_ipk_) :: i

    info = psb_success_
    call wk%free(info)
    call psb_geasb(wk%vx2l,desc,info,&
         & scratch=.true.,mold=vmold)
    call psb_geasb(wk%vy2l,desc,info,&
         & scratch=.true.,mold=vmold)
    call psb_geasb(wk%vtx,desc,info,&
         & scratch=.true.,mold=vmold)
    call psb_geasb(wk%vty,desc,info,&
         & scratch=.true.,mold=vmold)
    allocate(wk%wv(nwv),stat=info)
    do i=1,nwv
      call psb_geasb(wk%wv(i),desc,info,&
           & scratch=.true.,mold=vmold)
    end do
    
  end subroutine z_wrk_alloc
    
  subroutine z_wrk_free(wk,info)

    Implicit None

    ! Arguments
    class(mld_zmlprec_wrk_type), target, intent(inout) :: wk
    integer(psb_ipk_), intent(out)                    :: info 
    !
    integer(psb_ipk_) :: i
    info = psb_success_

    if (allocated(wk%tx)) deallocate(wk%tx, stat=info)
    if (allocated(wk%ty)) deallocate(wk%ty, stat=info)
    if (allocated(wk%x2l)) deallocate(wk%x2l, stat=info)
    if (allocated(wk%y2l)) deallocate(wk%y2l, stat=info)
    call wk%vtx%free(info)
    call wk%vty%free(info)
    call wk%vx2l%free(info)
    call wk%vy2l%free(info)
    if (allocated(wk%wv)) then
      do i=1,size(wk%wv)
        call wk%wv(i)%free(info)
      end do
      deallocate(wk%wv, stat=info)
    end if

  end subroutine z_wrk_free
    
  subroutine z_wrk_clone(wk,wkout,info)
    use psb_base_mod
    Implicit None

    ! Arguments
    class(mld_zmlprec_wrk_type), target, intent(inout) :: wk
    class(mld_zmlprec_wrk_type), target, intent(inout) :: wkout
    integer(psb_ipk_), intent(out)                    :: info 
    !
    integer(psb_ipk_) :: i
    info = psb_success_
    
    call psb_safe_ab_cpy(wk%tx,wkout%tx,info)
    call psb_safe_ab_cpy(wk%ty,wkout%ty,info)
    call psb_safe_ab_cpy(wk%x2l,wkout%x2l,info)
    call psb_safe_ab_cpy(wk%y2l,wkout%y2l,info)
    call wk%vtx%clone(wkout%vtx,info)
    call wk%vty%clone(wkout%vty,info)
    call wk%vx2l%clone(wkout%vx2l,info)
    call wk%vy2l%clone(wkout%vy2l,info)
    if (allocated(wkout%wv)) then
      do i=1,size(wkout%wv)
        call wkout%wv(i)%free(info)
      end do
      deallocate( wkout%wv)
    end if
    allocate(wkout%wv(size(wk%wv)),stat=info)
    do i=1,size(wk%wv)
      call wk%wv(i)%clone(wkout%wv(i),info)
    end do
    return

  end subroutine z_wrk_clone
  
  subroutine z_wrk_move_alloc(wk, b,info)
    implicit none
    class(mld_zmlprec_wrk_type), target, intent(inout) :: wk, b
    integer(psb_ipk_), intent(out) :: info 
    
    call b%free(info)
    call move_alloc(wk%tx,b%tx)
    call move_alloc(wk%ty,b%ty)
    call move_alloc(wk%x2l,b%x2l)
    call move_alloc(wk%y2l,b%y2l)
    !
    ! Should define V%move_alloc....
    call move_alloc(wk%vtx%v,b%vtx%v)
    call move_alloc(wk%vty%v,b%vty%v)
    call move_alloc(wk%vx2l%v,b%vx2l%v)
    call move_alloc(wk%vy2l%v,b%vy2l%v)
    call move_alloc(wk%wv,b%wv)
        
  end subroutine z_wrk_move_alloc

  subroutine z_wrk_cnv(wk,info,vmold)
    use psb_base_mod
    
    Implicit None
    
    ! Arguments
    class(mld_zmlprec_wrk_type), target, intent(inout) :: wk
    integer(psb_ipk_), intent(out)                    :: info 
    class(psb_z_base_vect_type), intent(in), optional  :: vmold
    !
    integer(psb_ipk_) :: i

    info = psb_success_
    if (present(vmold)) then
      call wk%vtx%cnv(vmold)
      call wk%vty%cnv(vmold)
      call wk%vx2l%cnv(vmold)
      call wk%vy2l%cnv(vmold)
      if (allocated(wk%wv)) then
        do i=1,size(wk%wv)
          call wk%wv(i)%cnv(vmold)
        end do
      end if
    end if
  end subroutine z_wrk_cnv

  function z_wrk_sizeof(wk) result(val)
    use psb_realloc_mod
    implicit none 
    class(mld_zmlprec_wrk_type), intent(in) :: wk
    integer(psb_long_int_k_) :: val
    integer :: i
    val = 0
    val = val + psb_size(wk%tx)
    val = val + psb_size(wk%ty)
    val = val + psb_size(wk%x2l)
    val = val + psb_size(wk%y2l)
    val = val + wk%vtx%sizeof()
    val = val + wk%vty%sizeof()
    val = val + wk%vx2l%sizeof()
    val = val + wk%vy2l%sizeof()
    if (allocated(wk%wv)) then
      do i=1, size(wk%wv)
        val = val + wk%wv(i)%sizeof()
      end do
    end if
  end function z_wrk_sizeof
 
end module mld_z_onelev_mod
