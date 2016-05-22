!!$ 
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.3)
!!$  
!!$  (C) Copyright 2008, 2010, 2012, 2015
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
! File: mld_inner_mod.f90
!
! Module: mld_inner_mod
!
!  This module defines the interfaces to the real/complex, single/double
!  precision versions of the MLD2P4 routines, except those of the user level,
!  whose interfaces are defined in mld_prec_mod.f90.
!
module mld_z_inner_mod
  use mld_z_prec_type

  interface mld_mlprec_bld
    subroutine mld_zmlprec_bld(a,desc_a,prec,info, amold, vmold,imold)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_i_base_vect_type, &
           & psb_dpk_, psb_z_base_sparse_mat, psb_z_base_vect_type, psb_ipk_
      use mld_z_prec_type, only : mld_zprec_type
      implicit none
      type(psb_zspmat_type), intent(in), target          :: a
      type(psb_desc_type), intent(inout), target           :: desc_a
      type(mld_zprec_type), intent(inout), target        :: prec
      integer(psb_ipk_), intent(out)                       :: info
      class(psb_z_base_sparse_mat), intent(in), optional :: amold
      class(psb_z_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
!!$      character, intent(in),optional             :: upd
    end subroutine mld_zmlprec_bld
  end interface mld_mlprec_bld


  interface mld_mlprec_aply
    subroutine mld_zmlprec_aply(alpha,p,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_, psb_ipk_
      use mld_z_prec_type, only : mld_zprec_type
      implicit none 
      type(psb_desc_type),intent(in)        :: desc_data
      type(mld_zprec_type), intent(inout) :: p
      complex(psb_dpk_),intent(in)         :: alpha,beta
      complex(psb_dpk_),intent(inout)      :: x(:)
      complex(psb_dpk_),intent(inout)      :: y(:)
      character,intent(in)               :: trans
      complex(psb_dpk_),target             :: work(:)
      integer(psb_ipk_), intent(out)     :: info
    end subroutine mld_zmlprec_aply
    subroutine mld_zmlprec_aply_vect(alpha,p,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, &
           & psb_dpk_, psb_z_vect_type, psb_ipk_
      use mld_z_prec_type, only : mld_zprec_type
      implicit none 
      type(psb_desc_type),intent(in)        :: desc_data
      type(mld_zprec_type), intent(inout) :: p
      complex(psb_dpk_),intent(in)            :: alpha,beta
      type(psb_z_vect_type),intent(inout) :: x
      type(psb_z_vect_type),intent(inout) :: y
      character,intent(in)                  :: trans
      complex(psb_dpk_),target                :: work(:)
      integer(psb_ipk_), intent(out)        :: info
    end subroutine mld_zmlprec_aply_vect
  end interface mld_mlprec_aply


  interface mld_coarse_bld
    subroutine mld_zcoarse_bld(a,desc_a,p,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_, psb_ipk_
      use mld_z_prec_type, only : mld_z_onelev_type
      implicit none 
      type(psb_zspmat_type), intent(in), target      :: a
      type(psb_desc_type), intent(in), target          :: desc_a
      type(mld_z_onelev_type), intent(inout), target :: p
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine mld_zcoarse_bld
  end interface mld_coarse_bld


  interface mld_bld_mlhier_aggsize
    subroutine mld_z_bld_mlhier_aggsize(casize,mxplevs,mnaggratio,a,desc_a,precv,info)
      use psb_base_mod, only : psb_ipk_, psb_zspmat_type, psb_desc_type,psb_dpk_
      use mld_z_prec_type, only : mld_z_onelev_type
      implicit none 
      integer(psb_ipk_), intent(in)    :: casize,mxplevs
      real(psb_dpk_)                   :: mnaggratio
      type(psb_zspmat_type),intent(in), target           :: a
      type(psb_desc_type), intent(inout), target         :: desc_a
      type(mld_z_onelev_type), allocatable, target, intent(inout)  :: precv(:)
      integer(psb_ipk_), intent(out)   :: info
    end subroutine mld_z_bld_mlhier_aggsize
  end interface mld_bld_mlhier_aggsize

  interface mld_bld_mlhier_array
    subroutine mld_z_bld_mlhier_array(nplevs,a,desc_a,precv,info)
      use psb_base_mod, only : psb_ipk_, psb_zspmat_type, psb_desc_type
      use mld_z_prec_type, only : mld_z_onelev_type
      implicit none 
      integer(psb_ipk_), intent(inout)            :: nplevs
      type(psb_zspmat_type),intent(in), target :: a
      type(psb_desc_type), intent(inout), target  :: desc_a
      type(mld_z_onelev_type), allocatable, target, intent(inout)  :: precv(:)
      integer(psb_ipk_), intent(out)   :: info
    end subroutine mld_z_bld_mlhier_array
  end interface mld_bld_mlhier_array
  
  interface mld_aggrmap_bld
    subroutine mld_zaggrmap_bld(aggr_type,iorder,theta,a,desc_a,ilaggr,nlaggr,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_, psb_ipk_
      implicit none 
      integer(psb_ipk_), intent(in)     :: iorder
      integer(psb_ipk_), intent(in)       :: aggr_type
      real(psb_dpk_), intent(in)           :: theta
      type(psb_zspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)     :: desc_a
      integer(psb_ipk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
      integer(psb_ipk_), intent(out)      :: info
    end subroutine mld_zaggrmap_bld
  end interface mld_aggrmap_bld


  interface  mld_dec_map_bld
    subroutine mld_z_dec_map_bld(iorder,theta,a,desc_a,nlaggr,ilaggr,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_, psb_ipk_
      implicit none 
      integer(psb_ipk_), intent(in)     :: iorder
      type(psb_zspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)    :: desc_a
      real(psb_dpk_), intent(in)         :: theta
      integer(psb_ipk_), allocatable, intent(out)  :: ilaggr(:),nlaggr(:)
      integer(psb_ipk_), intent(out)               :: info
    end subroutine mld_z_dec_map_bld
  end interface mld_dec_map_bld


  interface mld_aggrmat_asb
    subroutine mld_zaggrmat_asb(a,desc_a,ilaggr,nlaggr,p,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_, psb_ipk_
      use mld_z_prec_type, only : mld_z_onelev_type
      implicit none 
      type(psb_zspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                  :: desc_a
      integer(psb_ipk_), intent(inout)                 :: ilaggr(:), nlaggr(:)
      type(mld_z_onelev_type), intent(inout), target :: p
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine mld_zaggrmat_asb
  end interface mld_aggrmat_asb

  

  abstract interface
    subroutine mld_zaggrmat_var_asb(a,desc_a,ilaggr,nlaggr,parms,ac,op_prol,op_restr,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_, psb_ipk_
      use mld_z_prec_type, only :  mld_z_onelev_type, mld_dml_parms
      implicit none 
      type(psb_zspmat_type), intent(in)           :: a
      type(psb_desc_type), intent(in)               :: desc_a
      integer(psb_ipk_), intent(inout)              :: ilaggr(:), nlaggr(:)
      type(mld_dml_parms), intent(inout)         :: parms 
      type(psb_zspmat_type), intent(out)          :: ac,op_prol,op_restr
      integer(psb_ipk_), intent(out)                :: info
    end subroutine mld_zaggrmat_var_asb
  end interface


  procedure(mld_zaggrmat_var_asb) ::  mld_zaggrmat_nosmth_asb, &
       & mld_zaggrmat_smth_asb, mld_zaggrmat_minnrg_asb, &
       & mld_zaggrmat_biz_asb


end module mld_z_inner_mod
