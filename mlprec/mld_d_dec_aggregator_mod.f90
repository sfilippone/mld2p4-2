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
!  Basic (decoupled) aggregation algorithm. Based on the ideas in
!    M. Brezina and P. Vanek, A black-box iterative solver based on a 
!    two-level Schwarz method, Computing,  63 (1999), 233-263.
!    P. Vanek, J. Mandel and M. Brezina, Algebraic Multigrid by Smoothed
!    Aggregation for Second and Fourth Order Elliptic Problems, Computing, 56
!    (1996), 179-196.
!    
module mld_d_dec_aggregator_mod

  use mld_d_base_aggregator_mod
  !> \namespace  mld_d_dec_aggregator_mod  \class  mld_d_dec_aggregator_type
  !! \extends mld_d_base_aggregator_mod::mld_d_base_aggregator_type
  !!
  !!   type, extends(mld_d_base_aggregator_type) :: mld_d_dec_aggregator_type
  !!       procedure(mld_d_soc_map_bld), nopass, pointer :: soc_map_bld => null()
  !!   end type
  !!  
  !!   This is the simplest aggregation method: starting from the
  !!   strength-of-connection measure for defining the aggregation
  !!   presented in 
  !!
  !!    M. Brezina and P. Vanek, A black-box iterative solver based on a 
  !!    two-level Schwarz method, Computing,  63 (1999), 233-263.
  !!    P. Vanek, J. Mandel and M. Brezina, Algebraic Multigrid by Smoothed
  !!    Aggregation for Second and Fourth Order Elliptic Problems, Computing, 56
  !!    (1996), 179-196.
  !!
  !!   it achieves parallelization by simply acting on the local matrix,
  !!   i.e. by "decoupling" the subdomains.
  !!   The data structure hosts a "map_bld" function pointer which allows to
  !!   choose other ways to measure "strength-of-connection", of which the
  !!   Vanek-Brezina-Mandel is the default. More details are available in 
  !!
  !!    P. D'Ambra, D. di Serafino and S. Filippone, On the development of
  !!    PSBLAS-based parallel two-level Schwarz preconditioners, Appl. Num. Math.
  !!    57 (2007), 1181-1196.
  !!
  !!    The soc_map_bld method is used inside the implementation of build_tprol
  !!
  !
  !
  type, extends(mld_d_base_aggregator_type) :: mld_d_dec_aggregator_type
    procedure(mld_d_soc_map_bld), nopass, pointer :: soc_map_bld => null()
    
  contains
    procedure, pass(ag) :: bld_tprol     => mld_d_dec_aggregator_build_tprol
    procedure, pass(ag) :: mat_bld       => mld_d_dec_aggregator_mat_bld
    procedure, pass(ag) :: mat_asb       => mld_d_dec_aggregator_mat_asb
    procedure, pass(ag) :: default       => mld_d_dec_aggregator_default
    procedure, pass(ag) :: set_aggr_type => mld_d_dec_aggregator_set_aggr_type
    procedure, pass(ag) :: descr         => mld_d_dec_aggregator_descr
    procedure, nopass   :: fmt           => mld_d_dec_aggregator_fmt
  end type mld_d_dec_aggregator_type


  procedure(mld_d_soc_map_bld) ::  mld_d_soc1_map_bld, mld_d_soc2_map_bld

  interface
    subroutine  mld_d_dec_aggregator_build_tprol(ag,parms,ag_data,&
         & a,desc_a,ilaggr,nlaggr,op_prol,info)
      import :: mld_d_dec_aggregator_type, psb_desc_type, psb_dspmat_type, psb_dpk_,  &
           & psb_ipk_, psb_long_int_k_, mld_dml_parms, mld_daggr_data
      implicit none
      class(mld_d_dec_aggregator_type), target, intent(inout) :: ag
      type(mld_dml_parms), intent(inout)  :: parms 
      type(mld_daggr_data), intent(in)    :: ag_data
      type(psb_dspmat_type), intent(inout) :: a
      type(psb_desc_type), intent(inout)     :: desc_a
      integer(psb_ipk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
      type(psb_dspmat_type), intent(out)  :: op_prol
      integer(psb_ipk_), intent(out)      :: info
    end subroutine mld_d_dec_aggregator_build_tprol
  end interface

  interface
    subroutine  mld_d_dec_aggregator_mat_bld(ag,parms,a,desc_a,ilaggr,nlaggr,ac,&
         & op_prol,op_restr,info)
      import :: mld_d_dec_aggregator_type, psb_desc_type, psb_dspmat_type, psb_dpk_,  &
           & psb_ipk_, psb_long_int_k_, mld_dml_parms
      implicit none
      class(mld_d_dec_aggregator_type), target, intent(inout) :: ag
      type(mld_dml_parms), intent(inout)   :: parms 
      type(psb_dspmat_type), intent(in)    :: a
      type(psb_desc_type), intent(in)      :: desc_a
      integer(psb_ipk_), intent(inout)     :: ilaggr(:), nlaggr(:)
      type(psb_dspmat_type), intent(inout)   :: op_prol
      type(psb_dspmat_type), intent(out)   :: ac,op_restr
      integer(psb_ipk_), intent(out)       :: info
    end subroutine mld_d_dec_aggregator_mat_bld
  end interface  

  interface
    subroutine  mld_d_dec_aggregator_mat_asb(ag,parms,a,desc_a,ilaggr,nlaggr,&
         & ac,desc_ac,op_prol,op_restr,info)
      import :: mld_d_dec_aggregator_type, psb_desc_type, psb_dspmat_type, psb_dpk_,  &
           & psb_ipk_, psb_long_int_k_, mld_dml_parms
      implicit none
      class(mld_d_dec_aggregator_type), target, intent(inout) :: ag
      type(mld_dml_parms), intent(inout)  :: parms 
      type(psb_dspmat_type), intent(in)    :: a
      type(psb_desc_type), intent(in)        :: desc_a
      integer(psb_ipk_), intent(inout)       :: ilaggr(:), nlaggr(:)
      type(psb_dspmat_type), intent(inout) :: op_prol, ac, op_restr
      type(psb_desc_type), intent(inout)     :: desc_ac
      integer(psb_ipk_), intent(out)         :: info
    end subroutine mld_d_dec_aggregator_mat_asb
  end interface  


contains

  subroutine  mld_d_dec_aggregator_set_aggr_type(ag,parms,info)
    use mld_base_prec_type
    implicit none 
    class(mld_d_dec_aggregator_type), intent(inout) :: ag
    type(mld_dml_parms), intent(in)   :: parms
    integer(psb_ipk_), intent(out) :: info

    select case(parms%aggr_type)
    case (mld_noalg_)
      ag%soc_map_bld => null()
    case (mld_soc1_)
      ag%soc_map_bld => mld_d_soc1_map_bld
    case (mld_soc2_)
      ag%soc_map_bld => mld_d_soc2_map_bld
    case default
      write(0,*) 'Unknown aggregation type, defaulting to SOC1'
      ag%soc_map_bld => mld_d_soc1_map_bld
    end select
    
    return
  end subroutine mld_d_dec_aggregator_set_aggr_type

  
  subroutine  mld_d_dec_aggregator_default(ag)
    implicit none 
    class(mld_d_dec_aggregator_type), intent(inout) :: ag

    call ag%mld_d_base_aggregator_type%default()
    ag%soc_map_bld => mld_d_soc1_map_bld
    
    return
  end subroutine mld_d_dec_aggregator_default

  function mld_d_dec_aggregator_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "Decoupled aggregation"
  end function mld_d_dec_aggregator_fmt
  
  subroutine  mld_d_dec_aggregator_descr(ag,parms,iout,info)
    implicit none 
    class(mld_d_dec_aggregator_type), intent(in) :: ag
    type(mld_dml_parms), intent(in)   :: parms
    integer(psb_ipk_), intent(in)  :: iout
    integer(psb_ipk_), intent(out) :: info

    write(iout,*) 'Decoupled Aggregator'
    write(iout,*) 'Aggregator object type: ',ag%fmt()
    call parms%mldescr(iout,info)
    
    return
  end subroutine mld_d_dec_aggregator_descr

end module mld_d_dec_aggregator_mod
