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
module mld_s_inner_mod

  use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_i_base_vect_type, &
       & psb_spk_, psb_s_base_sparse_mat, psb_s_base_vect_type, psb_ipk_, &
       & psb_s_vect_type
  use mld_s_prec_type, only : mld_sprec_type, mld_sml_parms, &
       & mld_s_onelev_type, mld_smlprec_wrk_type

  interface mld_mlprec_bld
    subroutine mld_smlprec_bld(a,desc_a,prec,info, amold, vmold,imold)
      import :: psb_sspmat_type, psb_desc_type, psb_i_base_vect_type, &
           & psb_spk_, psb_s_base_sparse_mat, psb_s_base_vect_type, psb_ipk_
      import :: mld_sprec_type
      implicit none
      type(psb_sspmat_type), intent(in), target          :: a
      type(psb_desc_type), intent(inout), target           :: desc_a
      type(mld_sprec_type), intent(inout), target        :: prec
      integer(psb_ipk_), intent(out)                       :: info
      class(psb_s_base_sparse_mat), intent(in), optional :: amold
      class(psb_s_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine mld_smlprec_bld
  end interface mld_mlprec_bld

  interface mld_mlprec_aply
    subroutine mld_smlprec_aply(alpha,p,x,beta,y,desc_data,trans,work,info)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      import :: mld_sprec_type
      implicit none 
      type(psb_desc_type),intent(in)        :: desc_data
      type(mld_sprec_type), intent(inout) :: p
      real(psb_spk_),intent(in)         :: alpha,beta
      real(psb_spk_),intent(inout)      :: x(:)
      real(psb_spk_),intent(inout)      :: y(:)
      character,intent(in)               :: trans
      real(psb_spk_),target             :: work(:)
      integer(psb_ipk_), intent(out)     :: info
    end subroutine mld_smlprec_aply
    subroutine mld_smlprec_aply_vect(alpha,p,x,beta,y,desc_data,trans,work,info)
      import :: psb_sspmat_type, psb_desc_type, &
           & psb_spk_, psb_s_vect_type, psb_ipk_
      import :: mld_sprec_type
      implicit none 
      type(psb_desc_type),intent(in)        :: desc_data
      type(mld_sprec_type), intent(inout) :: p
      real(psb_spk_),intent(in)            :: alpha,beta
      type(psb_s_vect_type),intent(inout) :: x
      type(psb_s_vect_type),intent(inout) :: y
      character,intent(in)                  :: trans
      real(psb_spk_),target                :: work(:)
      integer(psb_ipk_), intent(out)        :: info
    end subroutine mld_smlprec_aply_vect
  end interface mld_mlprec_aply
 
  interface mld_aggrmap_bld
    subroutine mld_s_lev_aggrmap_bld(p,a,desc_a,ilaggr,nlaggr,op_prol,info)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      import :: mld_s_onelev_type
      implicit none 
      type(mld_s_onelev_type), intent(inout), target :: p
      type(psb_sspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)     :: desc_a
      integer(psb_ipk_), intent(out) :: ilaggr(:),nlaggr(:)
      type(psb_sspmat_type), intent(out)  :: op_prol
      integer(psb_ipk_), intent(out)      :: info
    end subroutine mld_s_lev_aggrmap_bld
    subroutine mld_saggrmap_bld(aggr_type,iorder,theta,a,desc_a,ilaggr,nlaggr,op_prol,info)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      implicit none 
      integer(psb_ipk_), intent(in)     :: iorder
      integer(psb_ipk_), intent(in)       :: aggr_type
      real(psb_spk_), intent(in)           :: theta
      type(psb_sspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)     :: desc_a
      integer(psb_ipk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
      type(psb_sspmat_type), intent(out)  :: op_prol        
      integer(psb_ipk_), intent(out)      :: info
    end subroutine mld_saggrmap_bld
  end interface mld_aggrmap_bld


  interface  mld_dec_map_bld
    subroutine mld_s_dec_map_bld(iorder,theta,a,desc_a,nlaggr,ilaggr,info)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      implicit none 
      integer(psb_ipk_), intent(in)     :: iorder
      type(psb_sspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)    :: desc_a
      real(psb_spk_), intent(in)         :: theta
      integer(psb_ipk_), allocatable, intent(out)  :: ilaggr(:),nlaggr(:)
      integer(psb_ipk_), intent(out)               :: info
    end subroutine mld_s_dec_map_bld
  end interface mld_dec_map_bld

  interface  mld_hyb_map_bld
    subroutine mld_s_hyb_map_bld(iorder,theta,a,desc_a,nlaggr,ilaggr,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      implicit none 
      integer(psb_ipk_), intent(in)     :: iorder
      type(psb_sspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)    :: desc_a
      real(psb_spk_), intent(in)         :: theta
      integer(psb_ipk_), allocatable, intent(out)  :: ilaggr(:),nlaggr(:)
      integer(psb_ipk_), intent(out)               :: info
    end subroutine mld_s_hyb_map_bld
  end interface mld_hyb_map_bld

  interface mld_map_to_tprol
    subroutine mld_s_map_to_tprol(desc_a,ilaggr,nlaggr,op_prol,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      use mld_s_prec_type, only : mld_s_onelev_type
      implicit none 
      type(psb_desc_type), intent(in)     :: desc_a
      integer(psb_ipk_), allocatable, intent(inout) :: ilaggr(:),nlaggr(:)
      type(psb_sspmat_type), intent(out)  :: op_prol
      integer(psb_ipk_), intent(out)      :: info
    end subroutine mld_s_map_to_tprol
  end interface mld_map_to_tprol
  

  interface mld_lev_mat_asb
    subroutine mld_s_lev_aggrmat_asb(p,a,desc_a,ilaggr,nlaggr,op_prol,info)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      import :: mld_s_onelev_type
      implicit none 
      type(mld_s_onelev_type), intent(inout), target :: p
      type(psb_sspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)     :: desc_a
      integer(psb_ipk_), intent(inout) :: ilaggr(:),nlaggr(:)
      type(psb_sspmat_type), intent(inout)  :: op_prol
      integer(psb_ipk_), intent(out)      :: info
    end subroutine mld_s_lev_aggrmat_asb
  end interface mld_lev_mat_asb

  interface mld_aggrmat_asb
    subroutine mld_saggrmat_asb(a,desc_a,ilaggr,nlaggr,parms,ac,op_prol,op_restr,info)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      import :: mld_sml_parms
      implicit none 
      type(psb_sspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                  :: desc_a
      integer(psb_ipk_), intent(inout)                 :: ilaggr(:), nlaggr(:)
      type(mld_sml_parms), intent(inout)         :: parms 
      type(psb_sspmat_type), intent(out)          :: ac,op_prol,op_restr
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine mld_saggrmat_asb
  end interface mld_aggrmat_asb

  abstract interface
    subroutine mld_saggrmat_var_asb(a,desc_a,ilaggr,nlaggr,parms,ac,op_prol,op_restr,info)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      import ::  mld_s_onelev_type, mld_sml_parms
      implicit none 
      type(psb_sspmat_type), intent(in)           :: a
      type(psb_desc_type), intent(in)               :: desc_a
      integer(psb_ipk_), intent(inout)              :: ilaggr(:), nlaggr(:)
      type(mld_sml_parms), intent(inout)         :: parms 
      type(psb_sspmat_type), intent(out)          :: ac,op_prol,op_restr
      integer(psb_ipk_), intent(out)                :: info
    end subroutine mld_saggrmat_var_asb
  end interface

  procedure(mld_saggrmat_var_asb) ::  mld_saggrmat_nosmth_asb, &
       & mld_saggrmat_smth_asb, mld_saggrmat_minnrg_asb, &
       & mld_saggrmat_biz_asb

end module mld_s_inner_mod
