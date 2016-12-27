!!$ 
!!$ 
!!$                           MLD2P4  version 2.1
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.3)
!!$  
!!$  (C) Copyright 2008, 2010, 2012, 2015, 2017 
!!$
!!$                      Salvatore Filippone  Cranfield University
!!$		      Ambra Abdullahi Hassan University of Rome Tor Vergata
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
! File: mld_d_prec_mod.f90
!
! Module: mld_d_prec_mod
!
!  This module defines the interfaces to the real/complex, single/double
!  precision versions of the user-level MLD2P4 routines.
!
module mld_d_prec_mod

  use mld_d_prec_type
  use mld_d_jac_smoother
  use mld_d_as_smoother
  use mld_d_id_solver
  use mld_d_diag_solver
  use mld_d_ilu_solver
  use mld_d_gs_solver

  interface mld_precinit
    subroutine mld_dprecinit(p,ptype,info,nlev)
      import :: psb_dspmat_type, psb_desc_type, psb_dpk_, &
           & mld_dprec_type, psb_ipk_
      type(mld_dprec_type), intent(inout)    :: p
      character(len=*), intent(in)             :: ptype
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), optional, intent(in)  :: nlev
    end subroutine mld_dprecinit
  end interface mld_precinit

  interface mld_precset
    module procedure mld_d_iprecsetsm, mld_d_iprecsetsv, &
         & mld_d_iprecseti, mld_d_iprecsetc, mld_d_iprecsetr, &
         & mld_d_cprecseti, mld_d_cprecsetc, mld_d_cprecsetr, &
         & mld_d_iprecsetag
  end interface mld_precset

  interface mld_precbld
    subroutine mld_dprecbld(a,desc_a,prec,info,amold,vmold,imold)
      import :: psb_dspmat_type, psb_desc_type, psb_dpk_, &
           & psb_d_base_sparse_mat, psb_d_base_vect_type, &
           & psb_i_base_vect_type, mld_dprec_type, psb_ipk_
      implicit none
      type(psb_dspmat_type), intent(in), target          :: a
      type(psb_desc_type), intent(inout), target           :: desc_a
      type(mld_dprec_type), intent(inout), target        :: prec
      integer(psb_ipk_), intent(out)                       :: info
      class(psb_d_base_sparse_mat), intent(in), optional :: amold
      class(psb_d_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
      !      character, intent(in),optional             :: upd
    end subroutine mld_dprecbld
  end interface mld_precbld

  interface mld_hierarchy_bld
    subroutine mld_d_hierarchy_bld(a,desc_a,prec,info)
      import :: psb_dspmat_type, psb_desc_type, psb_dpk_, &
           & mld_dprec_type, psb_ipk_
      implicit none
      type(psb_dspmat_type), intent(in), target          :: a
      type(psb_desc_type), intent(inout), target           :: desc_a
      type(mld_dprec_type), intent(inout), target        :: prec
      integer(psb_ipk_), intent(out)                       :: info
      !      character, intent(in),optional             :: upd
    end subroutine mld_d_hierarchy_bld
  end interface mld_hierarchy_bld

  interface mld_extprol_bld
    subroutine mld_d_extprol_bld(a,desc_a,p,prolv,restrv,info,amold,vmold,imold)
      import :: psb_dspmat_type, psb_desc_type, psb_dpk_, &
           & psb_d_base_sparse_mat, psb_d_base_vect_type, &
           & psb_i_base_vect_type, mld_dprec_type, psb_ipk_

      ! Arguments
      type(psb_dspmat_type),intent(in), target           :: a
      type(psb_dspmat_type),intent(inout), target        :: prolv(:)
      type(psb_dspmat_type),intent(inout), target        :: restrv(:)
      type(psb_desc_type), intent(inout), target         :: desc_a
      type(mld_dprec_type),intent(inout),target          :: p
      integer(psb_ipk_), intent(out)                       :: info
      class(psb_d_base_sparse_mat), intent(in), optional :: amold
      class(psb_d_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
      ! !$  character, intent(in), optional         :: upd
    end subroutine mld_d_extprol_bld
  end interface mld_extprol_bld
  
  interface mld_smoothers_bld
    subroutine mld_d_smoothers_bld(a,desc_a,prec,info,amold,vmold,imold)
      import :: psb_dspmat_type, psb_desc_type, psb_dpk_, &
           & psb_d_base_sparse_mat, psb_d_base_vect_type, &
           & psb_i_base_vect_type, mld_dprec_type, psb_ipk_
      implicit none
      type(psb_dspmat_type), intent(in), target          :: a
      type(psb_desc_type), intent(inout), target           :: desc_a
      type(mld_dprec_type), intent(inout), target        :: prec
      integer(psb_ipk_), intent(out)                       :: info
      class(psb_d_base_sparse_mat), intent(in), optional :: amold
      class(psb_d_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
      !      character, intent(in),optional             :: upd
    end subroutine mld_d_smoothers_bld
  end interface mld_smoothers_bld

contains

  subroutine mld_d_iprecsetsm(p,val,info,pos)
    type(mld_dprec_type), intent(inout)    :: p
    class(mld_d_base_smoother_type), intent(in)   :: val
    integer(psb_ipk_), intent(out)           :: info
    character(len=*), optional, intent(in)      :: pos

    call p%set(val,info,pos=pos)
  end subroutine mld_d_iprecsetsm

  subroutine mld_d_iprecsetsv(p,val,info,pos)
    type(mld_dprec_type), intent(inout)    :: p
    class(mld_d_base_solver_type), intent(in)   :: val
    integer(psb_ipk_), intent(out)                :: info
    character(len=*), optional, intent(in)      :: pos
    call p%set(val,info, pos=pos)
  end subroutine mld_d_iprecsetsv

  subroutine mld_d_iprecsetag(p,val,info,pos)
    type(mld_dprec_type), intent(inout)    :: p
    class(mld_d_base_aggregator_type), intent(in)   :: val
    integer(psb_ipk_), intent(out)                :: info
    character(len=*), optional, intent(in)      :: pos
    call p%set(val,info, pos=pos)
  end subroutine mld_d_iprecsetag

  subroutine mld_d_iprecseti(p,what,val,info,pos)
    type(mld_dprec_type), intent(inout)    :: p
    integer(psb_ipk_), intent(in)            :: what 
    integer(psb_ipk_), intent(in)            :: val
    integer(psb_ipk_), intent(out)           :: info
    character(len=*), optional, intent(in)   :: pos

    call p%set(what,val,info,pos=pos)
  end subroutine mld_d_iprecseti

  subroutine mld_d_iprecsetr(p,what,val,info,pos)
    type(mld_dprec_type), intent(inout)    :: p
    integer(psb_ipk_), intent(in)            :: what 
    real(psb_dpk_), intent(in)             :: val
    integer(psb_ipk_), intent(out)           :: info
    character(len=*), optional, intent(in)      :: pos

    call p%set(what,val,info,pos=pos)
  end subroutine mld_d_iprecsetr

  subroutine mld_d_iprecsetc(p,what,val,info,pos)
    type(mld_dprec_type), intent(inout)   :: p
    integer(psb_ipk_), intent(in)           :: what 
    character(len=*), intent(in)            :: val
    integer(psb_ipk_), intent(out)          :: info
    character(len=*), optional, intent(in)      :: pos

    call p%set(what,val,info,pos=pos)
  end subroutine mld_d_iprecsetc

  subroutine mld_d_cprecseti(p,what,val,info,pos)
    type(mld_dprec_type), intent(inout)   :: p
    character(len=*), intent(in)            :: what 
    integer(psb_ipk_), intent(in)           :: val
    integer(psb_ipk_), intent(out)          :: info
    character(len=*), optional, intent(in)      :: pos

    call p%set(what,val,info,pos=pos)
  end subroutine mld_d_cprecseti

  subroutine mld_d_cprecsetr(p,what,val,info,pos)
    type(mld_dprec_type), intent(inout)   :: p
    character(len=*), intent(in)            :: what 
    real(psb_dpk_), intent(in)             :: val
    integer(psb_ipk_), intent(out)          :: info
    character(len=*), optional, intent(in)      :: pos

    call p%set(what,val,info,pos=pos)
  end subroutine mld_d_cprecsetr

  subroutine mld_d_cprecsetc(p,what,val,info,pos)
    type(mld_dprec_type), intent(inout)   :: p
    character(len=*), intent(in)            :: what 
    character(len=*), intent(in)            :: val
    integer(psb_ipk_), intent(out)          :: info
    character(len=*), optional, intent(in)      :: pos

    call p%set(what,val,info,pos=pos)
  end subroutine mld_d_cprecsetc

end module mld_d_prec_mod
