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
! File: mld_c_prec_mod.f90
!
! Module: mld_c_prec_mod
!
!  This module defines the interfaces to the real/complex, single/double
!  precision versions of the user-level MLD2P4 routines.
!
module mld_c_prec_mod

  use mld_c_prec_type
  use mld_c_jac_smoother
  use mld_c_as_smoother
  use mld_c_id_solver
  use mld_c_diag_solver
  use mld_c_l1_diag_solver
  use mld_c_ilu_solver
  use mld_c_gs_solver

  interface mld_precset
    module procedure mld_c_iprecsetsm, mld_c_iprecsetsv, &
         & mld_c_cprecseti, mld_c_cprecsetc, mld_c_cprecsetr, &
         & mld_c_iprecsetag
  end interface mld_precset

  interface mld_extprol_bld
    subroutine mld_c_extprol_bld(a,desc_a,p,prolv,restrv,info,amold,vmold,imold)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, &
           & psb_c_base_sparse_mat, psb_c_base_vect_type, &
           & psb_i_base_vect_type, mld_cprec_type, psb_ipk_

      ! Arguments
      type(psb_cspmat_type),intent(in), target           :: a
      type(psb_cspmat_type),intent(inout), target        :: prolv(:)
      type(psb_cspmat_type),intent(inout), target        :: restrv(:)
      type(psb_desc_type), intent(inout), target         :: desc_a
      type(mld_cprec_type),intent(inout),target          :: p
      integer(psb_ipk_), intent(out)                       :: info
      class(psb_c_base_sparse_mat), intent(in), optional :: amold
      class(psb_c_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
      ! !$  character, intent(in), optional         :: upd
    end subroutine mld_c_extprol_bld
  end interface mld_extprol_bld
  
contains

  subroutine mld_c_iprecsetsm(p,val,info,pos)
    type(mld_cprec_type), intent(inout)    :: p
    class(mld_c_base_smoother_type), intent(in)   :: val
    integer(psb_ipk_), intent(out)           :: info
    character(len=*), optional, intent(in)      :: pos

    call p%set(val,info,pos=pos)
  end subroutine mld_c_iprecsetsm

  subroutine mld_c_iprecsetsv(p,val,info,pos)
    type(mld_cprec_type), intent(inout)    :: p
    class(mld_c_base_solver_type), intent(in)   :: val
    integer(psb_ipk_), intent(out)                :: info
    character(len=*), optional, intent(in)      :: pos
    call p%set(val,info, pos=pos)
  end subroutine mld_c_iprecsetsv

  subroutine mld_c_iprecsetag(p,val,info,pos)
    type(mld_cprec_type), intent(inout)    :: p
    class(mld_c_base_aggregator_type), intent(in)   :: val
    integer(psb_ipk_), intent(out)                :: info
    character(len=*), optional, intent(in)      :: pos
    call p%set(val,info, pos=pos)
  end subroutine mld_c_iprecsetag

  subroutine mld_c_cprecseti(p,what,val,info,pos)
    type(mld_cprec_type), intent(inout)   :: p
    character(len=*), intent(in)            :: what 
    integer(psb_ipk_), intent(in)           :: val
    integer(psb_ipk_), intent(out)          :: info
    character(len=*), optional, intent(in)      :: pos

    call p%set(what,val,info,pos=pos)
  end subroutine mld_c_cprecseti

  subroutine mld_c_cprecsetr(p,what,val,info,pos)
    type(mld_cprec_type), intent(inout)   :: p
    character(len=*), intent(in)            :: what 
    real(psb_spk_), intent(in)             :: val
    integer(psb_ipk_), intent(out)          :: info
    character(len=*), optional, intent(in)      :: pos

    call p%set(what,val,info,pos=pos)
  end subroutine mld_c_cprecsetr

  subroutine mld_c_cprecsetc(p,what,val,info,pos)
    type(mld_cprec_type), intent(inout)   :: p
    character(len=*), intent(in)            :: what 
    character(len=*), intent(in)            :: val
    integer(psb_ipk_), intent(out)          :: info
    character(len=*), optional, intent(in)      :: pos

    call p%set(what,val,info,pos=pos)
  end subroutine mld_c_cprecsetc

end module mld_c_prec_mod
