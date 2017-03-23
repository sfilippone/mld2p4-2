!  
!   
!                             MLD2P4  version 2.1
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.3)
!    
!    (C) Copyright 2008, 2010, 2012, 2015, 2017 
!  
!                        Salvatore Filippone  Cranfield University
!  		      Ambra Abdullahi Hassan University of Rome Tor Vergata
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
! File: mld_c_onelev_mod.f90
!
! Module: mld_c_onelev_mod
!
!  This module defines: 
!  - the mld_c_onelev_type data structure containing one level
!    of a multilevel  preconditioner and related
!    data structures;
!
!  It contains routines for
!  - Building and applying; 
!  - checking if the preconditioner is correctly defined;
!  - printing a	description of the preconditioner;
!  - deallocating the preconditioner data structure.  
!

module mld_c_onelev_mod

  use mld_base_prec_type
  use mld_c_base_smoother_mod
  use mld_c_base_aggregator_mod
  use psb_base_mod, only : psb_cspmat_type, psb_c_vect_type, &
       & psb_c_base_vect_type, psb_clinmap_type, psb_spk_, &
       & psb_ipk_, psb_long_int_k_, psb_desc_type, psb_i_base_vect_type, &
       & psb_erractionsave, psb_error_handler
  !
  !
  ! Type: mld_Tonelev_type.
  !
  !  It is the data type containing the necessary items for the	current
  !  level (essentially, the smoother, the current-level matrix
  !  and the restriction and prolongation operators).
  !
  !  type mld_Tonelev_type
  !    class(mld_T_base_smoother_type), allocatable :: sm
  !    type(mld_RTml_parms)            :: parms 
  !    type(psb_Tspmat_type)           :: ac
  !    type(psb_Tesc_type)             :: desc_ac
  !    type(psb_Tspmat_type), pointer  :: base_a    => null() 
  !    type(psb_Tesc_type), pointer    :: base_desc => null() 
  !    type(psb_Tlinmap_type)          :: map
  !  end type mld_Tonelev_type
  !
  !  Note that psb_Tpk denotes the kind of the real data type to be chosen
  !  according to single/double precision version of MLD2P4.
  !
  !   sm           -  class(mld_T_base_smoother_type), allocatable
  !                   The current level preconditioner (aka smoother).
  !   parms        -  type(mld_RTml_parms)
  !                   The parameters defining the multilevel strategy.
  !   ac           -  The local part of the current-level matrix, built by
  !                   coarsening the previous-level matrix.
  !   desc_ac      -  type(psb_desc_type).
  !                   The communication descriptor associated to the matrix
  !                   stored in ac.
  !   base_a       -  type(psb_Tspmat_type), pointer.
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
  !    check      -   Sanity checks.
  !    sizeof     -   Total memory occupation in bytes
  !    get_nzeros -   Number of nonzeros 
  !
  !
  type mld_c_onelev_type
    class(mld_c_base_smoother_type), allocatable :: sm, sm2a
    class(mld_c_base_smoother_type), pointer :: sm2 => null()
    class(mld_c_base_aggregator_type), allocatable :: aggr
    type(mld_sml_parms)              :: parms 
    type(psb_cspmat_type)            :: ac
    integer(psb_ipk_)                :: ac_nz_loc, ac_nz_tot
    type(psb_desc_type)              :: desc_ac
    type(psb_cspmat_type), pointer   :: base_a    => null() 
    type(psb_desc_type), pointer     :: base_desc => null() 
    type(psb_cspmat_type)            :: tprol
    type(psb_clinmap_type)           :: map
    real(psb_spk_)                     :: szratio
  contains
    procedure, pass(lv) :: bld     => mld_c_base_onelev_build
    procedure, pass(lv) :: clone   => c_base_onelev_clone
    procedure, pass(lv) :: cnv     => mld_c_base_onelev_cnv
    procedure, pass(lv) :: descr   => mld_c_base_onelev_descr
    procedure, pass(lv) :: default => c_base_onelev_default
    procedure, pass(lv) :: free    => mld_c_base_onelev_free
    procedure, pass(lv) :: nullify => c_base_onelev_nullify
    procedure, pass(lv) :: check => mld_c_base_onelev_check
    procedure, pass(lv) :: dump  => mld_c_base_onelev_dump
    procedure, pass(lv) :: seti  => mld_c_base_onelev_seti
    procedure, pass(lv) :: setr  => mld_c_base_onelev_setr
    procedure, pass(lv) :: setc  => mld_c_base_onelev_setc
    procedure, pass(lv) :: cseti => mld_c_base_onelev_cseti
    procedure, pass(lv) :: csetr => mld_c_base_onelev_csetr
    procedure, pass(lv) :: csetc => mld_c_base_onelev_csetc
    procedure, pass(lv) :: setsm => mld_c_base_onelev_setsm
    procedure, pass(lv) :: setsv => mld_c_base_onelev_setsv
    procedure, pass(lv) :: setag => mld_c_base_onelev_setag
    generic, public     :: set   => seti, setr, setc, &
         & cseti, csetr, csetc, setsm, setsv, setag 
    procedure, pass(lv) :: sizeof => c_base_onelev_sizeof
    procedure, pass(lv) :: get_nzeros => c_base_onelev_get_nzeros
    procedure, nopass   :: stringval => mld_stringval
    procedure, pass(lv) :: move_alloc => c_base_onelev_move_alloc
  end type mld_c_onelev_type

  type mld_c_onelev_node
    type(mld_c_onelev_type) :: item
    type(mld_c_onelev_node), pointer :: prev=>null(), next=>null()
  end type mld_c_onelev_node

  private :: c_base_onelev_default, c_base_onelev_sizeof, &
       & c_base_onelev_nullify, c_base_onelev_get_nzeros, &
       & c_base_onelev_clone, c_base_onelev_move_alloc



  interface
    subroutine mld_c_base_onelev_build(lv,info,amold,vmold,imold)
      import :: psb_c_base_sparse_mat, psb_c_base_vect_type, &
           & psb_i_base_vect_type, psb_spk_, mld_c_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      implicit none
      class(mld_c_onelev_type), target, intent(inout) :: lv
      integer(psb_ipk_), intent(out) :: info
      class(psb_c_base_sparse_mat), intent(in), optional :: amold
      class(psb_c_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine mld_c_base_onelev_build
  end interface

  interface 
    subroutine mld_c_base_onelev_descr(lv,il,nl,ilmin,info,iout)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_clinmap_type, psb_spk_, mld_c_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      ! Arguments
      class(mld_c_onelev_type), intent(in) :: lv
      integer(psb_ipk_), intent(in)                 :: il,nl,ilmin
      integer(psb_ipk_), intent(out)                :: info
      integer(psb_ipk_), intent(in), optional       :: iout
    end subroutine mld_c_base_onelev_descr
  end interface

  interface 
    subroutine mld_c_base_onelev_cnv(lv,info,amold,vmold,imold)
      import :: mld_c_onelev_type, psb_c_base_vect_type, psb_spk_, &
           & psb_c_base_sparse_mat, psb_ipk_, psb_i_base_vect_type
      ! Arguments
      class(mld_c_onelev_type), intent(inout)            :: lv 
      integer(psb_ipk_), intent(out)                     :: info
      class(psb_c_base_sparse_mat), intent(in), optional :: amold
      class(psb_c_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine mld_c_base_onelev_cnv
  end interface
   
  interface 
    subroutine mld_c_base_onelev_free(lv,info)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_clinmap_type, psb_spk_, mld_c_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      implicit none 
      
      class(mld_c_onelev_type), intent(inout) :: lv
      integer(psb_ipk_), intent(out)                :: info
    end subroutine mld_c_base_onelev_free
  end interface
  
  interface 
    subroutine mld_c_base_onelev_check(lv,info)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_clinmap_type, psb_spk_, mld_c_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      ! Arguments
      class(mld_c_onelev_type), intent(inout) :: lv 
      integer(psb_ipk_), intent(out)            :: info
    end subroutine mld_c_base_onelev_check
  end interface
  
  interface 
    subroutine mld_c_base_onelev_seti(lv,what,val,info,pos)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_clinmap_type, psb_spk_, mld_c_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      
      ! Arguments
      class(mld_c_onelev_type), intent(inout) :: lv 
      integer(psb_ipk_), intent(in)             :: what 
      integer(psb_ipk_), intent(in)             :: val
      integer(psb_ipk_), intent(out)            :: info
      character(len=*), optional, intent(in)      :: pos
    end subroutine mld_c_base_onelev_seti
  end interface

  interface 
    subroutine mld_c_base_onelev_setsm(lv,val,info,pos)
      import :: psb_spk_, mld_c_onelev_type, mld_c_base_smoother_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      
      ! Arguments
      class(mld_c_onelev_type), target, intent(inout) :: lv 
      class(mld_c_base_smoother_type), intent(in)     :: val
      integer(psb_ipk_), intent(out)                  :: info
      character(len=*), optional, intent(in)          :: pos
    end subroutine mld_c_base_onelev_setsm
  end interface
  
  interface 
    subroutine mld_c_base_onelev_setsv(lv,val,info,pos)
      import :: psb_spk_, mld_c_onelev_type, mld_c_base_solver_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      
      ! Arguments
      class(mld_c_onelev_type), target, intent(inout) :: lv 
      class(mld_c_base_solver_type), intent(in)       :: val
      integer(psb_ipk_), intent(out)                  :: info
      character(len=*), optional, intent(in)          :: pos
    end subroutine mld_c_base_onelev_setsv
  end interface
  
  interface 
    subroutine mld_c_base_onelev_setag(lv,val,info,pos)
      import :: psb_spk_, mld_c_onelev_type, mld_c_base_aggregator_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      
      ! Arguments
      class(mld_c_onelev_type), target, intent(inout) :: lv 
      class(mld_c_base_aggregator_type), intent(in)       :: val
      integer(psb_ipk_), intent(out)                  :: info
      character(len=*), optional, intent(in)          :: pos
    end subroutine mld_c_base_onelev_setag
  end interface
  
  interface 
    subroutine mld_c_base_onelev_setc(lv,what,val,info,pos)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_clinmap_type, psb_spk_, mld_c_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      ! Arguments
      class(mld_c_onelev_type), intent(inout) :: lv 
      integer(psb_ipk_), intent(in)             :: what 
      character(len=*), intent(in)              :: val
      integer(psb_ipk_), intent(out)            :: info
      character(len=*), optional, intent(in)      :: pos
    end subroutine mld_c_base_onelev_setc
  end interface
  
  interface 
    subroutine mld_c_base_onelev_setr(lv,what,val,info,pos)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_clinmap_type, psb_spk_, mld_c_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      
      class(mld_c_onelev_type), intent(inout) :: lv 
      integer(psb_ipk_), intent(in)             :: what 
      real(psb_spk_), intent(in)                 :: val
      integer(psb_ipk_), intent(out)            :: info
      character(len=*), optional, intent(in)      :: pos
    end subroutine mld_c_base_onelev_setr
  end interface

  
  interface 
    subroutine mld_c_base_onelev_cseti(lv,what,val,info,pos)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_clinmap_type, psb_spk_, mld_c_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      
      ! Arguments
      class(mld_c_onelev_type), intent(inout) :: lv 
      character(len=*), intent(in)              :: what 
      integer(psb_ipk_), intent(in)             :: val
      integer(psb_ipk_), intent(out)            :: info
      character(len=*), optional, intent(in)      :: pos
    end subroutine mld_c_base_onelev_cseti
  end interface
  
  interface 
    subroutine mld_c_base_onelev_csetc(lv,what,val,info,pos)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_clinmap_type, psb_spk_, mld_c_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      ! Arguments
      class(mld_c_onelev_type), intent(inout) :: lv 
      character(len=*), intent(in)              :: what 
      character(len=*), intent(in)              :: val
      integer(psb_ipk_), intent(out)            :: info
      character(len=*), optional, intent(in)      :: pos
    end subroutine mld_c_base_onelev_csetc
  end interface
  
  interface 
    subroutine mld_c_base_onelev_csetr(lv,what,val,info,pos)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_clinmap_type, psb_spk_, mld_c_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      Implicit None
      
      class(mld_c_onelev_type), intent(inout) :: lv 
      character(len=*), intent(in)              :: what 
      real(psb_spk_), intent(in)                 :: val
      integer(psb_ipk_), intent(out)            :: info
      character(len=*), optional, intent(in)      :: pos
    end subroutine mld_c_base_onelev_csetr
  end interface

  interface 
    subroutine mld_c_base_onelev_dump(lv,level,info,prefix,head,ac,rp,smoother,&
         & solver,global_num)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_clinmap_type, psb_spk_, mld_c_onelev_type, &
           & psb_ipk_, psb_long_int_k_, psb_desc_type
      implicit none 
      class(mld_c_onelev_type), intent(in) :: lv
      integer(psb_ipk_), intent(in)          :: level
      integer(psb_ipk_), intent(out)         :: info
      character(len=*), intent(in), optional :: prefix, head
      logical, optional, intent(in)    :: ac, rp, smoother, solver, global_num
    end subroutine mld_c_base_onelev_dump
  end interface
  
contains
  !
  ! Function returning the size of the mld_prec_type data structure
  ! in bytes or in number of nonzeros of the operator(s) involved. 
  !

  function c_base_onelev_get_nzeros(lv) result(val)
    implicit none 
    class(mld_c_onelev_type), intent(in) :: lv
    integer(psb_long_int_k_) :: val
    integer(psb_ipk_)        :: i
    val = 0
    if (allocated(lv%sm)) &
         &  val =  lv%sm%get_nzeros()
    if (allocated(lv%sm2a)) &
         &  val =  val + lv%sm2a%get_nzeros()
  end function c_base_onelev_get_nzeros

  function c_base_onelev_sizeof(lv) result(val)
    implicit none 
    class(mld_c_onelev_type), intent(in) :: lv
    integer(psb_long_int_k_) :: val
    integer(psb_ipk_)        :: i
    
    val = 0
    val = val + lv%desc_ac%sizeof()
    val = val + lv%ac%sizeof()
    val = val + lv%tprol%sizeof()
    val = val + lv%map%sizeof() 
    if (allocated(lv%sm))  val = val + lv%sm%sizeof()
    if (allocated(lv%sm2a))  val = val + lv%sm2a%sizeof()
  end function c_base_onelev_sizeof


  subroutine c_base_onelev_nullify(lv)
    implicit none 

    class(mld_c_onelev_type), intent(inout) :: lv

    nullify(lv%base_a) 
    nullify(lv%base_desc) 
    nullify(lv%sm2)
  end subroutine c_base_onelev_nullify

  !
  ! Multilevel defaults: 
  !  multiplicative vs. additive ML framework;
  !  Smoothed decoupled aggregation with zero threshold; 
  !  distributed coarse matrix;
  !  damping omega  computed with the max-norm estimate of the
  !  dominant eigenvalue;
  !  two-sided smoothing (i.e. V-cycle) with 1 smoothing sweep;
  !

  subroutine c_base_onelev_default(lv)

    Implicit None
 
    ! Arguments
    class(mld_c_onelev_type), target, intent(inout) :: lv
    integer(psb_ipk_) :: info 

    lv%parms%sweeps          = 1
    lv%parms%sweeps_pre      = 1
    lv%parms%sweeps_post     = 1
    lv%parms%ml_type         = mld_mult_ml_
    lv%parms%aggr_alg        = mld_dec_aggr_
    lv%parms%aggr_ord        = mld_aggr_ord_nat_
    lv%parms%aggr_kind       = mld_smooth_prol_
    lv%parms%coarse_mat      = mld_distr_mat_
    lv%parms%smoother_pos    = mld_twoside_smooth_
    lv%parms%aggr_omega_alg  = mld_eig_est_
    lv%parms%aggr_eig        = mld_max_norm_
    lv%parms%aggr_filter     = mld_no_filter_mat_
    lv%parms%aggr_omega_val  = czero
    lv%parms%aggr_thresh     = czero
    
    if (allocated(lv%sm)) call lv%sm%default()
    if (allocated(lv%sm2a)) then
      call lv%sm2a%default()
      lv%sm2 => lv%sm2a
    else
      lv%sm2 => lv%sm
    end if
    if (.not.allocated(lv%aggr)) allocate(mld_c_base_aggregator_type :: lv%aggr,stat=info)
    if (allocated(lv%aggr)) call lv%aggr%default()
    
    return

  end subroutine c_base_onelev_default



  subroutine c_base_onelev_clone(lv,lvout,info)

    Implicit None

    ! Arguments
    class(mld_c_onelev_type), target, intent(inout) :: lv 
    class(mld_c_onelev_type), target, intent(inout) :: lvout
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

  end subroutine c_base_onelev_clone


  subroutine c_base_onelev_move_alloc(lv, b,info)
    use psb_base_mod
    implicit none
    class(mld_c_onelev_type), target, intent(inout) :: lv, b
    integer(psb_ipk_), intent(out) :: info 
    
    call b%free(info)
    b%parms  = lv%parms
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
    
  end subroutine c_base_onelev_move_alloc

end module mld_c_onelev_mod
