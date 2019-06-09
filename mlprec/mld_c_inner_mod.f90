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
! File: mld_inner_mod.f90
!
! Module: mld_inner_mod
!
!  This module defines the interfaces to the real/complex, single/double
!  precision versions of inner MLD2P4 routines.
!  The interfaces  of the user level routines are defined in mld_prec_mod.f90.
!
module mld_c_inner_mod

  use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_i_base_vect_type, &
       & psb_spk_, psb_c_base_sparse_mat, psb_c_csr_sparse_mat, &
       & psb_c_base_vect_type, psb_ipk_, psb_c_vect_type
  use mld_c_prec_type, only : mld_cprec_type, mld_sml_parms, &
       & mld_c_onelev_type, mld_cmlprec_wrk_type

  interface mld_mlprec_bld
    subroutine mld_cmlprec_bld(a,desc_a,prec,info, amold, vmold,imold)
      import :: psb_cspmat_type, psb_desc_type, psb_i_base_vect_type, &
           & psb_spk_, psb_c_base_sparse_mat, psb_c_base_vect_type, psb_ipk_
      import :: mld_cprec_type
      implicit none
      type(psb_cspmat_type), intent(in), target          :: a
      type(psb_desc_type), intent(inout), target           :: desc_a
      type(mld_cprec_type), intent(inout), target        :: prec
      integer(psb_ipk_), intent(out)                       :: info
      class(psb_c_base_sparse_mat), intent(in), optional :: amold
      class(psb_c_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine mld_cmlprec_bld
  end interface mld_mlprec_bld

  interface mld_mlprec_aply
    subroutine mld_cmlprec_aply(alpha,p,x,beta,y,desc_data,trans,work,info)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      import :: mld_cprec_type
      implicit none 
      type(psb_desc_type),intent(in)        :: desc_data
      type(mld_cprec_type), intent(inout) :: p
      complex(psb_spk_),intent(in)         :: alpha,beta
      complex(psb_spk_),intent(inout)      :: x(:)
      complex(psb_spk_),intent(inout)      :: y(:)
      character,intent(in)               :: trans
      complex(psb_spk_),target             :: work(:)
      integer(psb_ipk_), intent(out)     :: info
    end subroutine mld_cmlprec_aply
    subroutine mld_cmlprec_aply_vect(alpha,p,x,beta,y,desc_data,trans,work,info)
      import :: psb_cspmat_type, psb_desc_type, &
           & psb_spk_, psb_c_vect_type, psb_ipk_
      import :: mld_cprec_type
      implicit none 
      type(psb_desc_type),intent(in)        :: desc_data
      type(mld_cprec_type), intent(inout) :: p
      complex(psb_spk_),intent(in)            :: alpha,beta
      type(psb_c_vect_type),intent(inout) :: x
      type(psb_c_vect_type),intent(inout) :: y
      character,intent(in)                  :: trans
      complex(psb_spk_),target                :: work(:)
      integer(psb_ipk_), intent(out)        :: info
    end subroutine mld_cmlprec_aply_vect
  end interface mld_mlprec_aply
 
  interface mld_map_to_tprol
    subroutine mld_c_map_to_tprol(desc_a,ilaggr,nlaggr,op_prol,info)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      use mld_c_prec_type, only : mld_c_onelev_type
      implicit none 
      type(psb_desc_type), intent(in)     :: desc_a
      integer(psb_ipk_), allocatable, intent(inout) :: ilaggr(:),nlaggr(:)
      type(psb_cspmat_type), intent(out)  :: op_prol
      integer(psb_ipk_), intent(out)      :: info
    end subroutine mld_c_map_to_tprol
  end interface mld_map_to_tprol

  abstract interface
    subroutine mld_caggrmat_var_bld(a,desc_a,ilaggr,nlaggr,parms,ac,op_prol,op_restr,info)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      import ::  mld_c_onelev_type, mld_sml_parms
      implicit none 
      type(psb_cspmat_type), intent(in)           :: a
      type(psb_desc_type), intent(in)               :: desc_a
      integer(psb_ipk_), intent(inout)              :: ilaggr(:), nlaggr(:)
      type(mld_sml_parms), intent(inout)         :: parms 
      type(psb_cspmat_type), intent(out)          :: ac,op_prol,op_restr
      integer(psb_ipk_), intent(out)                :: info
    end subroutine mld_caggrmat_var_bld
  end interface

  procedure(mld_caggrmat_var_bld) ::  mld_caggrmat_nosmth_bld, &
       & mld_caggrmat_smth_bld, mld_caggrmat_minnrg_bld, &
       & mld_caggrmat_biz_bld


  interface mld_par_spspmm
    subroutine mld_c_par_csr_spspmm(acsr,desc_a,bcsr,ccsr,desc_c,info,data)
      import :: psb_c_csr_sparse_mat, psb_desc_type, psb_ipk_
      Implicit None
      type(psb_c_csr_sparse_mat),intent(in)    :: acsr
      type(psb_c_csr_sparse_mat),intent(inout) :: bcsr
      type(psb_c_csr_sparse_mat),intent(out)   :: ccsr      
      type(psb_desc_type),intent(in)           :: desc_a
      type(psb_desc_type),intent(inout)        :: desc_c
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), intent(in), optional  :: data
    End Subroutine mld_c_par_csr_spspmm
  end interface mld_par_spspmm
  
end module mld_c_inner_mod
