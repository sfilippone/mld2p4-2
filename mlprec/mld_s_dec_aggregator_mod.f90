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
!
!
!  The aggregator object hosts the aggregation method for building
!  the multilevel hierarchy. 
!
!    M. Brezina and P. Vanek, A black-box iterative solver based on a 
!    two-level Schwarz method, Computing,  63 (1999), 233-263.
!    P. D'Ambra, D. di Serafino and S. Filippone, On the development of
!    PSBLAS-based parallel two-level Schwarz preconditioners, Appl. Num. Math.
!    57 (2007), 1181-1196.
!    
module mld_s_dec_aggregator_mod

  use mld_s_base_aggregator_mod
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
  type, extends(mld_s_base_aggregator_type) :: mld_s_dec_aggregator_type
    procedure(mld_s_map_bld), nopass, pointer :: map_bld => null()
    
  contains
    procedure, pass(ag) :: bld_tprol     => mld_s_dec_aggregator_build_tprol
    procedure, pass(ag) :: mat_asb       => mld_s_dec_aggregator_mat_asb
    procedure, pass(ag) :: default       => mld_s_dec_aggregator_default
    procedure, pass(ag) :: set_aggr_type => mld_s_dec_aggregator_set_aggr_type
    procedure, nopass   :: fmt           => mld_s_dec_aggregator_fmt
  end type mld_s_dec_aggregator_type


  abstract interface  
    subroutine mld_s_map_bld(iorder,theta,a,desc_a,nlaggr,ilaggr,info)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, psb_ipk_, psb_lpk_
      implicit none 
      integer(psb_ipk_), intent(in)     :: iorder
      type(psb_sspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)    :: desc_a
      real(psb_spk_), intent(in)         :: theta
      integer(psb_lpk_), allocatable, intent(out)  :: ilaggr(:),nlaggr(:)
      integer(psb_ipk_), intent(out)               :: info
    end subroutine mld_s_map_bld
  end interface

  procedure(mld_s_map_bld) ::  mld_s_soc1_map_bld, mld_s_soc2_map_bld

  interface
    subroutine  mld_s_dec_aggregator_build_tprol(ag,parms,a,desc_a,ilaggr,nlaggr,op_prol,info)
      import :: mld_s_dec_aggregator_type, psb_desc_type, psb_sspmat_type, psb_spk_,  &
           & psb_ipk_, psb_lpk_,  psb_lsspmat_type, mld_sml_parms 
      implicit none
      class(mld_s_dec_aggregator_type), target, intent(inout) :: ag
      type(mld_sml_parms), intent(inout)  :: parms 
      type(psb_sspmat_type), intent(in)   :: a
      type(psb_desc_type), intent(in)     :: desc_a
      integer(psb_lpk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
      type(psb_lsspmat_type), intent(out)  :: op_prol
      integer(psb_ipk_), intent(out)      :: info
    end subroutine mld_s_dec_aggregator_build_tprol
  end interface

  interface
    subroutine  mld_s_dec_aggregator_mat_asb(ag,parms,a,desc_a,ilaggr,nlaggr,ac,&
         & op_prol,op_restr,info)
      import :: mld_s_dec_aggregator_type, psb_desc_type, psb_sspmat_type, psb_spk_,  &
           & psb_ipk_, psb_lpk_, psb_lsspmat_type, mld_sml_parms
      implicit none
      class(mld_s_dec_aggregator_type), target, intent(inout) :: ag
      type(mld_sml_parms), intent(inout)   :: parms 
      type(psb_sspmat_type), intent(in)    :: a
      type(psb_desc_type), intent(in)      :: desc_a
      integer(psb_lpk_), intent(inout)     :: ilaggr(:), nlaggr(:)
      type(psb_lsspmat_type), intent(inout)   :: op_prol
      type(psb_sspmat_type), intent(out)   :: ac,op_restr
      integer(psb_ipk_), intent(out)       :: info
    end subroutine mld_s_dec_aggregator_mat_asb
  end interface  


contains

  subroutine  mld_s_dec_aggregator_set_aggr_type(ag,parms,info)
    use mld_base_prec_type
    implicit none 
    class(mld_s_dec_aggregator_type), intent(inout) :: ag
    type(mld_sml_parms), intent(in)   :: parms
    integer(psb_ipk_), intent(out) :: info

    select case(parms%aggr_type)
    case (mld_noalg_)
      ag%map_bld => null()
    case (mld_soc1_)
      ag%map_bld => mld_s_soc1_map_bld
    case (mld_soc2_)
      ag%map_bld => mld_s_soc2_map_bld
    case default
      write(0,*) 'Unknown aggregation type, defaulting to SOC1'
      ag%map_bld => mld_s_soc1_map_bld
    end select
    
    return
  end subroutine mld_s_dec_aggregator_set_aggr_type

  
  subroutine  mld_s_dec_aggregator_default(ag)
    implicit none 
    class(mld_s_dec_aggregator_type), intent(inout) :: ag

    ag%map_bld => mld_s_soc1_map_bld
    
    return
  end subroutine mld_s_dec_aggregator_default

  function mld_s_dec_aggregator_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "Decoupled aggregation"
  end function mld_s_dec_aggregator_fmt
  
end module mld_s_dec_aggregator_mod
