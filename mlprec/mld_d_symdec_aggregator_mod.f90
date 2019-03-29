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
!  Locally symmetrized (decoupled) aggregation algorithm.
!  This version differs from the basic decoupled aggregation algorithm
!  only because it works on (the pattern of)  A+A^T  instead of A. 
!
!    
module mld_d_symdec_aggregator_mod

  use mld_d_dec_aggregator_mod
  !> \namespace  mld_d_symdec_aggregator_mod  \class  mld_d_symdec_aggregator_type
  !! \extends mld_d_dec_aggregator_mod::mld_d_dec_aggregator_type
  !!
  !!  This version differs from the basic decoupled aggregation algorithm
  !!   only because it works on (the pattern of)  A+A^T  instead of A.
  !!   
  !
  type, extends(mld_d_dec_aggregator_type) :: mld_d_symdec_aggregator_type
    
  contains
    procedure, pass(ag) :: bld_tprol    => mld_d_symdec_aggregator_build_tprol
    procedure, pass(ag) :: descr        => mld_d_symdec_aggregator_descr
    procedure, nopass   :: fmt          => mld_d_symdec_aggregator_fmt
  end type mld_d_symdec_aggregator_type


  interface
    subroutine  mld_d_symdec_aggregator_build_tprol(ag,parms,ag_data,&
         & a,desc_a,ilaggr,nlaggr,op_prol,info)
      import :: mld_d_symdec_aggregator_type, psb_desc_type, psb_dspmat_type, psb_dpk_,  &
           & psb_ipk_, psb_long_int_k_, mld_dml_parms, mld_daggr_data
      implicit none
      class(mld_d_symdec_aggregator_type), target, intent(inout) :: ag
      type(mld_dml_parms), intent(inout)  :: parms 
      type(mld_daggr_data), intent(in)    :: ag_data
      type(psb_dspmat_type), intent(in)   :: a
      type(psb_desc_type), intent(in)     :: desc_a
      integer(psb_ipk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
      type(psb_dspmat_type), intent(out)  :: op_prol
      integer(psb_ipk_), intent(out)      :: info
    end subroutine mld_d_symdec_aggregator_build_tprol
  end interface


contains

  function mld_d_symdec_aggregator_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "Symmetric Decoupled aggregation"
  end function mld_d_symdec_aggregator_fmt
  
  subroutine  mld_d_symdec_aggregator_descr(ag,parms,iout,info)
    implicit none 
    class(mld_d_symdec_aggregator_type), intent(in) :: ag
    type(mld_dml_parms), intent(in)   :: parms
    integer(psb_ipk_), intent(in)  :: iout
    integer(psb_ipk_), intent(out) :: info

    write(iout,*) 'Decoupled Aggregator locally-symmetrized'
    write(iout,*) 'Aggregator object type: ',ag%fmt()
    call parms%mldescr(iout,info)
    
    return
  end subroutine mld_d_symdec_aggregator_descr

end module mld_d_symdec_aggregator_mod
