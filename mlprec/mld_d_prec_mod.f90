!!$ 
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010
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
! File: mld_prec_mod.f90
!
! Module: mld_prec_mod
!
!  This module defines the interfaces to the real/complex, single/double
!  precision versions of the user-level MLD2P4 routines.
!
module mld_d_prec_mod

  use mld_d_prec_type
  use mld_d_move_alloc_mod

  interface mld_precinit
    subroutine mld_dprecinit(p,ptype,info,nlev)
      use psb_sparse_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_d_prec_type, only : mld_dprec_type
      type(mld_dprec_type), intent(inout)    :: p
      character(len=*), intent(in)           :: ptype
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: nlev
    end subroutine mld_dprecinit
  end interface

  interface mld_precset
    module procedure mld_i_dprecsetsm, mld_i_dprecsetsv, &
         & mld_i_dprecseti, mld_i_dprecsetc, mld_i_dprecsetr
  end interface

  interface mld_inner_precset
    subroutine mld_dprecsetsm(p,val,info,ilev)
      use psb_sparse_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_d_prec_type, only : mld_dprec_type, mld_d_base_smoother_type
      type(mld_dprec_type), intent(inout)    :: p
      class(mld_d_base_smoother_type), intent(in) :: val
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: ilev
    end subroutine mld_dprecsetsm
    subroutine mld_dprecsetsv(p,val,info,ilev)
      use psb_sparse_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_d_prec_type, only : mld_dprec_type, mld_d_base_solver_type
      type(mld_dprec_type), intent(inout)    :: p
      class(mld_d_base_solver_type), intent(in) :: val
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: ilev
    end subroutine mld_dprecsetsv
    subroutine mld_dprecseti(p,what,val,info,ilev)
      use psb_sparse_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_d_prec_type, only : mld_dprec_type
      type(mld_dprec_type), intent(inout)    :: p
      integer, intent(in)                    :: what 
      integer, intent(in)                    :: val
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: ilev
    end subroutine mld_dprecseti
    subroutine mld_dprecsetr(p,what,val,info,ilev)
      use psb_sparse_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_d_prec_type, only : mld_dprec_type
      type(mld_dprec_type), intent(inout)    :: p
      integer, intent(in)                    :: what 
      real(psb_dpk_), intent(in)             :: val
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: ilev
    end subroutine mld_dprecsetr
    subroutine mld_dprecsetc(p,what,string,info,ilev)
      use psb_sparse_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_d_prec_type, only : mld_dprec_type
      type(mld_dprec_type), intent(inout)    :: p
      integer, intent(in)                    :: what 
      character(len=*), intent(in)           :: string
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: ilev
    end subroutine mld_dprecsetc
  end interface

  interface mld_precbld
    subroutine mld_dprecbld(a,desc_a,prec,info)
      use psb_sparse_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_d_prec_type, only : mld_dprec_type
      implicit none
      type(psb_dspmat_type), intent(in), target   :: a
      type(psb_desc_type), intent(in), target     :: desc_a
      type(mld_dprec_type), intent(inout), target :: prec
      integer, intent(out)                        :: info
!!$      character, intent(in),optional             :: upd
    end subroutine mld_dprecbld
  end interface

contains

  subroutine mld_i_dprecsetsm(p,val,info)
    use psb_sparse_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
    use mld_d_prec_type, only : mld_dprec_type, mld_d_base_smoother_type
    type(mld_dprec_type), intent(inout)    :: p
    class(mld_d_base_smoother_type), intent(in)   :: val
    integer, intent(out)                   :: info

    call mld_inner_precset(p,val,info)
  end subroutine mld_i_dprecsetsm

  subroutine mld_i_dprecsetsv(p,val,info)
    use psb_sparse_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
    use mld_d_prec_type, only : mld_dprec_type, mld_d_base_solver_type
    type(mld_dprec_type), intent(inout)    :: p
    class(mld_d_base_solver_type), intent(in)   :: val
    integer, intent(out)                   :: info

    call mld_inner_precset(p,val,info)
  end subroutine mld_i_dprecsetsv

  subroutine mld_i_dprecseti(p,what,val,info)
    use psb_sparse_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
    use mld_d_prec_type, only : mld_dprec_type
    type(mld_dprec_type), intent(inout)    :: p
    integer, intent(in)                    :: what 
    integer, intent(in)                    :: val
    integer, intent(out)                   :: info

    call mld_inner_precset(p,what,val,info)
  end subroutine mld_i_dprecseti

  subroutine mld_i_dprecsetr(p,what,val,info)
    use psb_sparse_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
    use mld_d_prec_type, only : mld_dprec_type
    type(mld_dprec_type), intent(inout)    :: p
    integer, intent(in)                    :: what 
    real(psb_dpk_), intent(in)             :: val
    integer, intent(out)                   :: info

    call mld_inner_precset(p,what,val,info)
  end subroutine mld_i_dprecsetr

  subroutine mld_i_dprecsetc(p,what,val,info)
    use psb_sparse_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
    use mld_d_prec_type, only : mld_dprec_type
    type(mld_dprec_type), intent(inout)    :: p
    integer, intent(in)                    :: what 
    character(len=*), intent(in)           :: val
    integer, intent(out)                   :: info

    call mld_inner_precset(p,what,val,info)
  end subroutine mld_i_dprecsetc

end module mld_d_prec_mod
