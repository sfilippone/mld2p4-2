!!$ 
!!$ 
!!$                           MLD2P4  version 1.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 2.2)
!!$  
!!$  (C) Copyright 2008
!!$
!!$                      Salvatore Filippone  University of Rome Tor Vergata       
!!$                      Alfredo Buttari      University of Rome Tor Vergata
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
!  This module defines the interfaces to the real and complex versions of the
!  MLD2P4 routines.
!
module mld_prec_mod

  use mld_prec_type

  interface mld_precinit
    subroutine mld_dprecinit(p,ptype,info,nlev)
      use psb_base_mod
      use mld_prec_type
      type(mld_dprec_type), intent(inout)    :: p
      character(len=*), intent(in)           :: ptype
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: nlev
    end subroutine mld_dprecinit
    subroutine mld_zprecinit(p,ptype,info,nlev)
      use psb_base_mod
      use mld_prec_type
      type(mld_zprec_type), intent(inout)    :: p
      character(len=*), intent(in)           :: ptype
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: nlev
    end subroutine mld_zprecinit
  end interface

  interface mld_precset
    subroutine mld_dprecseti(p,what,val,info,ilev)
      use psb_base_mod
      use mld_prec_type
      type(mld_dprec_type), intent(inout)    :: p
      integer, intent(in)                    :: what 
      integer, intent(in)                    :: val
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: ilev
    end subroutine mld_dprecseti
    subroutine mld_dprecsetd(p,what,val,info,ilev)
      use psb_base_mod
      use mld_prec_type
      type(mld_dprec_type), intent(inout)    :: p
      integer, intent(in)                    :: what 
      real(kind(1.d0)), intent(in)           :: val
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: ilev
    end subroutine mld_dprecsetd
    subroutine mld_dprecsetc(p,what,string,info,ilev)
      use psb_base_mod
      use mld_prec_type
      type(mld_dprec_type), intent(inout)    :: p
      integer, intent(in)                    :: what 
      character(len=*), intent(in)           :: string
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: ilev
    end subroutine mld_dprecsetc
    subroutine mld_zprecseti(p,what,val,info,ilev)
      use psb_base_mod
      use mld_prec_type
      type(mld_zprec_type), intent(inout)    :: p
      integer, intent(in)                    :: what 
      integer, intent(in)                    :: val
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: ilev
    end subroutine mld_zprecseti
    subroutine mld_zprecsetd(p,what,val,info,ilev)
      use psb_base_mod
      use mld_prec_type
      type(mld_zprec_type), intent(inout)    :: p
      integer, intent(in)                    :: what 
      real(kind(1.d0)), intent(in)           :: val
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: ilev
    end subroutine mld_zprecsetd
    subroutine mld_zprecsetc(p,what,string,info,ilev)
      use psb_base_mod
      use mld_prec_type
      type(mld_zprec_type), intent(inout)    :: p
      integer, intent(in)                    :: what 
      character(len=*), intent(in)           :: string
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: ilev
    end subroutine mld_zprecsetc
  end interface

  interface mld_precfree
    subroutine mld_dprecfree(p,info)
      use psb_base_mod
      use mld_prec_type
      type(mld_dprec_type), intent(inout) :: p
      integer, intent(out)                :: info
    end subroutine mld_dprecfree
    subroutine mld_zprecfree(p,info)
      use psb_base_mod
      use mld_prec_type
      type(mld_zprec_type), intent(inout) :: p
      integer, intent(out)                :: info
    end subroutine mld_zprecfree
  end interface

  interface mld_precaply
    subroutine mld_dprec_aply(prec,x,y,desc_data,info,trans,work)
      use psb_base_mod
      use mld_prec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(mld_dprec_type), intent(in)  :: prec
      real(kind(0.d0)),intent(in)       :: x(:)
      real(kind(0.d0)),intent(inout)    :: y(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
      real(kind(0.d0)),intent(inout), optional, target :: work(:)
    end subroutine mld_dprec_aply
    subroutine mld_dprec_aply1(prec,x,desc_data,info,trans)
      use psb_base_mod
      use mld_prec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(mld_dprec_type), intent(in)  :: prec
      real(kind(0.d0)),intent(inout)    :: x(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
    end subroutine mld_dprec_aply1
    subroutine mld_zprec_aply(prec,x,y,desc_data,info,trans,work)
      use psb_base_mod
      use mld_prec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(mld_zprec_type), intent(in)  :: prec
      complex(kind(0.d0)),intent(in)    :: x(:)
      complex(kind(0.d0)),intent(inout) :: y(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
      complex(kind(0.d0)),intent(inout), optional, target :: work(:)
    end subroutine mld_zprec_aply
    subroutine mld_zprec_aply1(prec,x,desc_data,info,trans)
      use psb_base_mod
      use mld_prec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(mld_zprec_type), intent(in)  :: prec
      complex(kind(0.d0)),intent(inout) :: x(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
    end subroutine mld_zprec_aply1
  end interface

  interface mld_precbld
    subroutine mld_dprecbld(a,desc_a,prec,info)
      use psb_base_mod
      use mld_prec_type
      implicit none
      type(psb_dspmat_type), intent(in), target  :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(mld_dprec_type), intent(inout)        :: prec
      integer, intent(out)                       :: info
!!$      character, intent(in),optional             :: upd
    end subroutine mld_dprecbld
    subroutine mld_zprecbld(a,desc_a,prec,info)
      use psb_base_mod
      use mld_prec_type
      implicit none
      type(psb_zspmat_type), intent(in), target  :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(mld_zprec_type), intent(inout)        :: prec
      integer, intent(out)                       :: info
!!$      character, intent(in),optional             :: upd
    end subroutine mld_zprecbld
  end interface

end module mld_prec_mod
