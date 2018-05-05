!  
!   
!                             MLD2P4  version 2.2
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!    
!    (C) Copyright 2008, 2010, 2012, 2015, 2017 , 2017 
!  
!                        Salvatore Filippone  Cranfield University
!  		      Ambra Abdullahi Hassan University of Rome Tor Vergata
!        Pasqua D'Ambra         IAC-CNR, Naples, IT
!        Daniela di Serafino    University of Campania "L. Vanvitelli", Caserta, IT
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
!
!
!  The aggregator object hosts the aggregation method for building
!  the multilevel hierarchy. This variant is based on the hybrid method
!  presented in 
!
!    S. Gratton, P. Henon, P. Jiranek and X. Vasseur:
!    Reducing complexity of algebraic multigrid by aggregation
!    Numerical Lin. Algebra with Applications, 2016, 23:501-518
!    
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

module bcm_csr_type_mod
 use iso_c_binding
 type, bind(c)::  bcm_Vector
   type(c_ptr) :: data
   integer(c_int) :: size
   integer(c_int) :: owns_data
 end type 

 type, bind(c)::  bcm_CSRMatrix
   type(c_ptr) :: i
   type(c_ptr) :: j
   integer(c_int) :: num_rows
   integer(c_int) :: num_cols
   integer(c_int) :: num_nonzeros
   integer(c_int) :: owns_data
   type(c_ptr) :: data
 end type 
end module bcm_csr_type_mod

module mld_d_bcmatch_aggregator_mod
  use mld_d_base_aggregator_mod
  use bcm_csr_type_mod

  type, extends(mld_d_base_aggregator_type) :: mld_d_bcmatch_aggregator_type
    integer(psb_ipk_) :: matching_alg
    integer(psb_ipk_) :: n_sweeps
    real(psb_dpk_), allocatable :: w_tmp(:) 
    type(bcm_Vector) :: w_par
    integer(psb_ipk_) :: max_csize
    integer(psb_ipk_) :: max_nlevels
    !type(psb_d_vect_type) :: w
  contains
    procedure, pass(ag) :: bld_tprol => mld_d_bcmatch_aggregator_build_tprol
    procedure, pass(ag) :: set    => d_bcmatch_aggr_cseti
    procedure, pass(ag) :: default    =>d_bcmatch_aggr_set_default
    procedure, pass(ag) :: mat_asb   => mld_d_bcmatch_aggregator_mat_asb
    procedure, pass(ag) :: update_level => d_bcmatch_aggregator_update_level
!!$    procedure, pass(ag) :: clone        => mld_d_base_aggregator_clone
!!$    procedure, pass(ag) :: free         => mld_d_bcmatch_aggregator_free
!!$    procedure, pass(ag) :: default      => mld_d_base_aggregator_default
    procedure, nopass   :: fmt          => mld_d_bcmatch_aggregator_fmt
  end type mld_d_bcmatch_aggregator_type


  interface
    subroutine  mld_d_bcmatch_aggregator_build_tprol(ag,parms,a,desc_a,ilaggr,nlaggr,op_prol,info)
      import :: mld_d_bcmatch_aggregator_type, psb_desc_type, psb_dspmat_type, psb_dpk_,  &
           & psb_ipk_, psb_long_int_k_, mld_dml_parms
      implicit none
      class(mld_d_bcmatch_aggregator_type), target, intent(inout) :: ag
      type(mld_dml_parms), intent(inout)  :: parms 
      type(psb_dspmat_type), intent(in)   :: a
      type(psb_desc_type), intent(in)     :: desc_a
      integer(psb_ipk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
      type(psb_dspmat_type), intent(out)  :: op_prol
      integer(psb_ipk_), intent(out)      :: info
    end subroutine mld_d_bcmatch_aggregator_build_tprol
  end interface

  interface
    subroutine  mld_d_bcmatch_aggregator_mat_asb(ag,parms,a,desc_a,ilaggr,nlaggr,ac,&
         & op_prol,op_restr,info)
      import :: mld_d_bcmatch_aggregator_type, psb_desc_type, psb_dspmat_type, psb_dpk_,  &
           & psb_ipk_, psb_long_int_k_, mld_dml_parms
      implicit none
      class(mld_d_bcmatch_aggregator_type), target, intent(inout) :: ag
      type(mld_dml_parms), intent(inout)   :: parms 
      type(psb_dspmat_type), intent(in)    :: a
      type(psb_desc_type), intent(in)      :: desc_a
      integer(psb_ipk_), intent(inout)     :: ilaggr(:), nlaggr(:)
      type(psb_dspmat_type), intent(inout)   :: op_prol
      type(psb_dspmat_type), intent(out)   :: ac,op_restr
      integer(psb_ipk_), intent(out)       :: info
    end subroutine mld_d_bcmatch_aggregator_mat_asb
  end interface  
  
  
contains


  function mld_d_bcmatch_aggregator_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "BootCMatch aggregation"
  end function mld_d_bcmatch_aggregator_fmt

  subroutine  d_bcmatch_aggregator_update_level(ag,agnext,info)
    implicit none 
    class(mld_d_bcmatch_aggregator_type), target, intent(inout) :: ag
    class(mld_d_base_aggregator_type), target, intent(inout) :: agnext
    integer(psb_ipk_), intent(out)       :: info

    !
    !
    select type(agnext)
    type is (mld_d_bcmatch_aggregator_type)
      agnext%matching_alg=ag%matching_alg
      agnext%n_sweeps=ag%n_sweeps 
      agnext%max_csize=ag%max_csize 
      agnext%max_nlevels=ag%max_nlevels 
      agnext%w_par=ag%w_par 
    end select
    info = 0 
  end subroutine d_bcmatch_aggregator_update_level

  subroutine d_bcmatch_aggr_cseti(ag,what,val,info)

    Implicit None

    ! Arguments
    class(mld_d_bcmatch_aggregator_type), intent(inout) :: ag
    character(len=*), intent(in)                  :: what
    integer(psb_ipk_), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_)  :: err_act, iwhat
    character(len=20)  :: name='d_bcmatch_aggr_cseti'
    info = psb_success_

    select case(what)
      case('BCM_MATCH_ALG')
        ag%matching_alg=val
      case('BCM_SWEEPS')
        ag%n_sweeps=val
      case('BCM_MAX_CSIZE')
        ag%max_csize=val
      case('BCM_MAX_NLEVELS')
        ag%max_nlevels=val
      case('BCM_W_SIZE')
        ag%w_par%size=val
        ag%w_par%owns_data=0
        allocate(ag%w_tmp(val))
        ag%w_tmp = 1.0_psb_dpk_
        call set_cloc(ag%w_tmp, ag%w_par)
      case default
    end select
    return
  contains
    subroutine set_cloc(vect,w_par)
      real(psb_dpk_), target :: vect(:)
      type(bcm_Vector) :: w_par

      w_par%data = c_loc(vect)
    end subroutine set_cloc

  end subroutine d_bcmatch_aggr_cseti

  subroutine d_bcmatch_aggr_set_default(ag)

    Implicit None

    ! Arguments
    class(mld_d_bcmatch_aggregator_type), intent(inout) :: ag
    character(len=20)  :: name='d_bcmatch_aggr_set_default'
    ag%matching_alg=0
    ag%n_sweeps=1
    ag%max_nlevels=36
    ag%max_csize=10

    return

  end subroutine d_bcmatch_aggr_set_default

!!$  subroutine  d_bcmatch_aggregator_free(ag,info)
!!$    implicit none 
!!$    class(mld_d_bcmatch_aggregator_type), target, intent(inout) :: ag
!!$    integer(psb_ipk_), intent(out)       :: info
!!$
!!$    info = 0 
!!$  end subroutine d_bcmatch_aggregator_free

  
end module mld_d_bcmatch_aggregator_mod
