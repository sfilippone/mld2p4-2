!!$ 
!!$ 
!!$                           MLD2P4  version 1.1
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 2.3.1)
!!$  
!!$  (C) Copyright 2008,2009
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
!  This module defines the interfaces to the real/complex, single/double
!  precision versions of the user-level MLD2P4 routines.
!
module mld_prec_mod

  use mld_prec_type

  interface mld_precinit
!!$    subroutine mld_sprecinit(p,ptype,info,nlev)
!!$      use psb_sparse_mod, only : psb_s_sparse_mat, psb_desc_type, psb_spk_
!!$      use mld_prec_type, only : mld_sprec_type
!!$      type(mld_sprec_type), intent(inout)    :: p
!!$      character(len=*), intent(in)           :: ptype
!!$      integer, intent(out)                   :: info
!!$      integer, optional, intent(in)          :: nlev
!!$    end subroutine mld_sprecinit
    subroutine mld_dprecinit(p,ptype,info,nlev)
      use psb_sparse_mod, only : psb_d_sparse_mat, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dprec_type
      type(mld_dprec_type), intent(inout)    :: p
      character(len=*), intent(in)           :: ptype
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: nlev
    end subroutine mld_dprecinit
!!$    subroutine mld_cprecinit(p,ptype,info,nlev)
!!$      use psb_sparse_mod, only : psb_c_sparse_mat, psb_desc_type, psb_spk_
!!$      use mld_prec_type, only : mld_cprec_type
!!$      type(mld_cprec_type), intent(inout)    :: p
!!$      character(len=*), intent(in)           :: ptype
!!$      integer, intent(out)                   :: info
!!$      integer, optional, intent(in)          :: nlev
!!$    end subroutine mld_cprecinit
!!$    subroutine mld_zprecinit(p,ptype,info,nlev)
!!$      use psb_sparse_mod, only : psb_z_sparse_mat, psb_desc_type, psb_dpk_
!!$      use mld_prec_type, only : mld_zprec_type
!!$      type(mld_zprec_type), intent(inout)    :: p
!!$      character(len=*), intent(in)           :: ptype
!!$      integer, intent(out)                   :: info
!!$      integer, optional, intent(in)          :: nlev
!!$    end subroutine mld_zprecinit
  end interface

  interface mld_precset
    module procedure  mld_i_dprecseti, mld_i_dprecsetc, mld_i_dprecsetr
!!$    module procedure mld_i_sprecseti, mld_i_sprecsetc, mld_i_sprecsetr,&
!!$         & mld_i_dprecseti, mld_i_dprecsetc, mld_i_dprecsetr,&
!!$         & mld_i_cprecseti, mld_i_cprecsetc, mld_i_cprecsetr,&
!!$         & mld_i_zprecseti, mld_i_zprecsetc, mld_i_zprecsetr
  end interface

  interface mld_inner_precset
!!$    subroutine mld_sprecseti(p,what,val,info,ilev)
!!$      use psb_sparse_mod, only : psb_s_sparse_mat, psb_desc_type, psb_spk_
!!$      use mld_prec_type, only : mld_sprec_type
!!$      type(mld_sprec_type), intent(inout)    :: p
!!$      integer, intent(in)                    :: what 
!!$      integer, intent(in)                    :: val
!!$      integer, intent(out)                   :: info
!!$      integer, optional, intent(in)          :: ilev
!!$    end subroutine mld_sprecseti
!!$    subroutine mld_sprecsetr(p,what,val,info,ilev)
!!$      use psb_sparse_mod, only : psb_s_sparse_mat, psb_desc_type, psb_spk_
!!$      use mld_prec_type, only : mld_sprec_type
!!$      type(mld_sprec_type), intent(inout)    :: p
!!$      integer, intent(in)                    :: what 
!!$      real(psb_spk_), intent(in)           :: val
!!$      integer, intent(out)                   :: info
!!$      integer, optional, intent(in)          :: ilev
!!$    end subroutine mld_sprecsetr
!!$    subroutine mld_sprecsetc(p,what,string,info,ilev)
!!$      use psb_sparse_mod, only : psb_s_sparse_mat, psb_desc_type, psb_spk_
!!$      use mld_prec_type, only : mld_sprec_type
!!$      type(mld_sprec_type), intent(inout)    :: p
!!$      integer, intent(in)                    :: what 
!!$      character(len=*), intent(in)           :: string
!!$      integer, intent(out)                   :: info
!!$      integer, optional, intent(in)          :: ilev
!!$    end subroutine mld_sprecsetc
    subroutine mld_dprecsetsm(p,what,val,info,ilev)
      use psb_sparse_mod, only : psb_d_sparse_mat, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dprec_type, mld_d_base_smoother_type
      type(mld_dprec_type), intent(inout)    :: p
      integer, intent(in)                    :: what 
      class(mld_d_base_smoother_type), intent(in) :: val
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: ilev
    end subroutine mld_dprecsetsm
    subroutine mld_dprecsetsv(p,what,val,info,ilev)
      use psb_sparse_mod, only : psb_d_sparse_mat, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dprec_type, mld_d_base_solver_type
      type(mld_dprec_type), intent(inout)    :: p
      integer, intent(in)                    :: what 
      class(mld_d_base_solver_type), intent(in) :: val
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: ilev
    end subroutine mld_dprecsetsv
    subroutine mld_dprecseti(p,what,val,info,ilev)
      use psb_sparse_mod, only : psb_d_sparse_mat, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dprec_type
      type(mld_dprec_type), intent(inout)    :: p
      integer, intent(in)                    :: what 
      integer, intent(in)                    :: val
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: ilev
    end subroutine mld_dprecseti
    subroutine mld_dprecsetr(p,what,val,info,ilev)
      use psb_sparse_mod, only : psb_d_sparse_mat, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dprec_type
      type(mld_dprec_type), intent(inout)    :: p
      integer, intent(in)                    :: what 
      real(psb_dpk_), intent(in)             :: val
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: ilev
    end subroutine mld_dprecsetr
    subroutine mld_dprecsetc(p,what,string,info,ilev)
      use psb_sparse_mod, only : psb_d_sparse_mat, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dprec_type
      type(mld_dprec_type), intent(inout)    :: p
      integer, intent(in)                    :: what 
      character(len=*), intent(in)           :: string
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: ilev
    end subroutine mld_dprecsetc
!!$    subroutine mld_cprecseti(p,what,val,info,ilev)
!!$      use psb_sparse_mod, only : psb_c_sparse_mat, psb_desc_type, psb_spk_
!!$      use mld_prec_type, only : mld_cprec_type
!!$      type(mld_cprec_type), intent(inout)    :: p
!!$      integer, intent(in)                    :: what 
!!$      integer, intent(in)                    :: val
!!$      integer, intent(out)                   :: info
!!$      integer, optional, intent(in)          :: ilev
!!$    end subroutine mld_cprecseti
!!$    subroutine mld_cprecsetr(p,what,val,info,ilev)
!!$      use psb_sparse_mod, only : psb_c_sparse_mat, psb_desc_type, psb_spk_
!!$      use mld_prec_type, only : mld_cprec_type
!!$      type(mld_cprec_type), intent(inout)    :: p
!!$      integer, intent(in)                    :: what 
!!$      real(psb_spk_), intent(in)             :: val
!!$      integer, intent(out)                   :: info
!!$      integer, optional, intent(in)          :: ilev
!!$    end subroutine mld_cprecsetr
!!$    subroutine mld_cprecsetc(p,what,string,info,ilev)
!!$      use psb_sparse_mod, only : psb_c_sparse_mat, psb_desc_type, psb_spk_
!!$      use mld_prec_type, only : mld_cprec_type
!!$      type(mld_cprec_type), intent(inout)    :: p
!!$      integer, intent(in)                    :: what 
!!$      character(len=*), intent(in)           :: string
!!$      integer, intent(out)                   :: info
!!$      integer, optional, intent(in)          :: ilev
!!$    end subroutine mld_cprecsetc
!!$    subroutine mld_zprecseti(p,what,val,info,ilev)
!!$      use psb_sparse_mod, only : psb_z_sparse_mat, psb_desc_type, psb_dpk_
!!$      use mld_prec_type, only : mld_zprec_type
!!$      type(mld_zprec_type), intent(inout)    :: p
!!$      integer, intent(in)                    :: what 
!!$      integer, intent(in)                    :: val
!!$      integer, intent(out)                   :: info
!!$      integer, optional, intent(in)          :: ilev
!!$    end subroutine mld_zprecseti
!!$    subroutine mld_zprecsetr(p,what,val,info,ilev)
!!$      use psb_sparse_mod, only : psb_z_sparse_mat, psb_desc_type, psb_dpk_
!!$      use mld_prec_type, only : mld_zprec_type
!!$      type(mld_zprec_type), intent(inout)    :: p
!!$      integer, intent(in)                    :: what 
!!$      real(psb_dpk_), intent(in)             :: val
!!$      integer, intent(out)                   :: info
!!$      integer, optional, intent(in)          :: ilev
!!$    end subroutine mld_zprecsetr
!!$    subroutine mld_zprecsetc(p,what,string,info,ilev)
!!$      use psb_sparse_mod, only : psb_z_sparse_mat, psb_desc_type, psb_dpk_
!!$      use mld_prec_type, only : mld_zprec_type
!!$      type(mld_zprec_type), intent(inout)    :: p
!!$      integer, intent(in)                    :: what 
!!$      character(len=*), intent(in)           :: string
!!$      integer, intent(out)                   :: info
!!$      integer, optional, intent(in)          :: ilev
!!$    end subroutine mld_zprecsetc
  end interface
!!$
!!$  interface mld_precaply
!!$    subroutine mld_sprecaply(prec,x,y,desc_data,info,trans,work)
!!$      use psb_sparse_mod, only : psb_s_sparse_mat, psb_desc_type, psb_spk_
!!$      use mld_prec_type, only : mld_sprec_type
!!$      type(psb_desc_type),intent(in)    :: desc_data
!!$      type(mld_sprec_type), intent(in)  :: prec
!!$      real(psb_spk_),intent(in)       :: x(:)
!!$      real(psb_spk_),intent(inout)    :: y(:)
!!$      integer, intent(out)              :: info
!!$      character(len=1), optional        :: trans
!!$      real(psb_spk_),intent(inout), optional, target :: work(:)
!!$    end subroutine mld_sprecaply
!!$    subroutine mld_sprecaply1(prec,x,desc_data,info,trans)
!!$      use psb_sparse_mod, only : psb_s_sparse_mat, psb_desc_type, psb_spk_
!!$      use mld_prec_type, only : mld_sprec_type
!!$      type(psb_desc_type),intent(in)    :: desc_data
!!$      type(mld_sprec_type), intent(in)  :: prec
!!$      real(psb_spk_),intent(inout)    :: x(:)
!!$      integer, intent(out)              :: info
!!$      character(len=1), optional        :: trans
!!$    end subroutine mld_sprecaply1
!!$    subroutine mld_dprecaply(prec,x,y,desc_data,info,trans,work)
!!$      use psb_sparse_mod, only : psb_d_sparse_mat, psb_desc_type, psb_dpk_
!!$      use mld_prec_type, only : mld_dprec_type
!!$      type(psb_desc_type),intent(in)    :: desc_data
!!$      type(mld_dprec_type), intent(in)  :: prec
!!$      real(psb_dpk_),intent(in)       :: x(:)
!!$      real(psb_dpk_),intent(inout)    :: y(:)
!!$      integer, intent(out)              :: info
!!$      character(len=1), optional        :: trans
!!$      real(psb_dpk_),intent(inout), optional, target :: work(:)
!!$    end subroutine mld_dprecaply
!!$    subroutine mld_dprecaply1(prec,x,desc_data,info,trans)
!!$      use psb_sparse_mod, only : psb_d_sparse_mat, psb_desc_type, psb_dpk_
!!$      use mld_prec_type, only : mld_dprec_type
!!$      type(psb_desc_type),intent(in)    :: desc_data
!!$      type(mld_dprec_type), intent(in)  :: prec
!!$      real(psb_dpk_),intent(inout)    :: x(:)
!!$      integer, intent(out)              :: info
!!$      character(len=1), optional        :: trans
!!$    end subroutine mld_dprecaply1
!!$    subroutine mld_cprecaply(prec,x,y,desc_data,info,trans,work)
!!$      use psb_sparse_mod, only : psb_c_sparse_mat, psb_desc_type, psb_spk_
!!$      use mld_prec_type, only : mld_cprec_type
!!$      type(psb_desc_type),intent(in)    :: desc_data
!!$      type(mld_cprec_type), intent(in)  :: prec
!!$      complex(psb_spk_),intent(in)    :: x(:)
!!$      complex(psb_spk_),intent(inout) :: y(:)
!!$      integer, intent(out)              :: info
!!$      character(len=1), optional        :: trans
!!$      complex(psb_spk_),intent(inout), optional, target :: work(:)
!!$    end subroutine mld_cprecaply
!!$    subroutine mld_cprecaply1(prec,x,desc_data,info,trans)
!!$      use psb_sparse_mod, only : psb_c_sparse_mat, psb_desc_type, psb_spk_
!!$      use mld_prec_type, only : mld_cprec_type
!!$      type(psb_desc_type),intent(in)    :: desc_data
!!$      type(mld_cprec_type), intent(in)  :: prec
!!$      complex(psb_spk_),intent(inout) :: x(:)
!!$      integer, intent(out)              :: info
!!$      character(len=1), optional        :: trans
!!$    end subroutine mld_cprecaply1
!!$    subroutine mld_zprecaply(prec,x,y,desc_data,info,trans,work)
!!$      use psb_sparse_mod, only : psb_z_sparse_mat, psb_desc_type, psb_dpk_
!!$      use mld_prec_type, only : mld_zprec_type
!!$      type(psb_desc_type),intent(in)    :: desc_data
!!$      type(mld_zprec_type), intent(in)  :: prec
!!$      complex(psb_dpk_),intent(in)    :: x(:)
!!$      complex(psb_dpk_),intent(inout) :: y(:)
!!$      integer, intent(out)              :: info
!!$      character(len=1), optional        :: trans
!!$      complex(psb_dpk_),intent(inout), optional, target :: work(:)
!!$    end subroutine mld_zprecaply
!!$    subroutine mld_zprecaply1(prec,x,desc_data,info,trans)
!!$      use psb_sparse_mod, only : psb_z_sparse_mat, psb_desc_type, psb_dpk_
!!$      use mld_prec_type, only : mld_zprec_type
!!$      type(psb_desc_type),intent(in)    :: desc_data
!!$      type(mld_zprec_type), intent(in)  :: prec
!!$      complex(psb_dpk_),intent(inout) :: x(:)
!!$      integer, intent(out)              :: info
!!$      character(len=1), optional        :: trans
!!$    end subroutine mld_zprecaply1
!!$  end interface
!!$
  interface mld_precbld
    subroutine mld_sprecbld(a,desc_a,prec,info)
      use psb_sparse_mod, only : psb_s_sparse_mat, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sprec_type
      implicit none
      type(psb_s_sparse_mat), intent(in), target   :: a
      type(psb_desc_type), intent(in), target     :: desc_a
      type(mld_sprec_type), intent(inout), target :: prec
      integer, intent(out)                        :: info
!!$      character, intent(in),optional             :: upd
    end subroutine mld_sprecbld
    subroutine mld_dprecbld(a,desc_a,prec,info)
      use psb_sparse_mod, only : psb_d_sparse_mat, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dprec_type
      implicit none
      type(psb_d_sparse_mat), intent(in), target   :: a
      type(psb_desc_type), intent(in), target     :: desc_a
      type(mld_dprec_type), intent(inout), target :: prec
      integer, intent(out)                        :: info
!!$      character, intent(in),optional             :: upd
    end subroutine mld_dprecbld
    subroutine mld_cprecbld(a,desc_a,prec,info)
      use psb_sparse_mod, only : psb_c_sparse_mat, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cprec_type
      implicit none
      type(psb_c_sparse_mat), intent(in), target   :: a
      type(psb_desc_type), intent(in), target     :: desc_a
      type(mld_cprec_type), intent(inout), target :: prec
      integer, intent(out)                        :: info
!!$      character, intent(in),optional             :: upd
    end subroutine mld_cprecbld
    subroutine mld_zprecbld(a,desc_a,prec,info)
      use psb_sparse_mod, only : psb_z_sparse_mat, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zprec_type
      implicit none
      type(psb_z_sparse_mat), intent(in), target  :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(mld_zprec_type), intent(inout)        :: prec
      integer, intent(out)                       :: info
!!$      character, intent(in),optional             :: upd
    end subroutine mld_zprecbld
  end interface

contains

!!$  subroutine mld_i_sprecseti(p,what,val,info)
!!$    use psb_sparse_mod, only : psb_s_sparse_mat, psb_desc_type, psb_spk_
!!$    use mld_prec_type, only : mld_sprec_type
!!$    type(mld_sprec_type), intent(inout)    :: p
!!$    integer, intent(in)                    :: what 
!!$    integer, intent(in)                    :: val
!!$    integer, intent(out)                   :: info
!!$
!!$    call mld_inner_precset(p,what,val,info)
!!$  end subroutine mld_i_sprecseti
!!$
!!$  subroutine mld_i_sprecsetr(p,what,val,info)
!!$    use psb_sparse_mod, only : psb_s_sparse_mat, psb_desc_type, psb_spk_
!!$    use mld_prec_type, only : mld_sprec_type
!!$    type(mld_sprec_type), intent(inout)    :: p
!!$    integer, intent(in)                    :: what 
!!$    real(psb_spk_), intent(in)           :: val
!!$    integer, intent(out)                   :: info
!!$
!!$    call mld_inner_precset(p,what,val,info)
!!$  end subroutine mld_i_sprecsetr
!!$
!!$  subroutine mld_i_sprecsetc(p,what,val,info)
!!$    use psb_sparse_mod, only : psb_s_sparse_mat, psb_desc_type, psb_spk_
!!$    use mld_prec_type, only : mld_sprec_type
!!$    type(mld_sprec_type), intent(inout)    :: p
!!$    integer, intent(in)                    :: what 
!!$    character(len=*), intent(in)           :: val
!!$    integer, intent(out)                   :: info
!!$
!!$    call mld_inner_precset(p,what,val,info)
!!$  end subroutine mld_i_sprecsetc

  subroutine mld_i_dprecseti(p,what,val,info)
    use psb_sparse_mod, only : psb_d_sparse_mat, psb_desc_type, psb_dpk_
    use mld_prec_type, only : mld_dprec_type
    type(mld_dprec_type), intent(inout)    :: p
    integer, intent(in)                    :: what 
    integer, intent(in)                    :: val
    integer, intent(out)                   :: info

    call mld_inner_precset(p,what,val,info)
  end subroutine mld_i_dprecseti

  subroutine mld_i_dprecsetr(p,what,val,info)
    use psb_sparse_mod, only : psb_d_sparse_mat, psb_desc_type, psb_dpk_
    use mld_prec_type, only : mld_dprec_type
    type(mld_dprec_type), intent(inout)    :: p
    integer, intent(in)                    :: what 
    real(psb_dpk_), intent(in)             :: val
    integer, intent(out)                   :: info

    call mld_inner_precset(p,what,val,info)
  end subroutine mld_i_dprecsetr

  subroutine mld_i_dprecsetc(p,what,val,info)
    use psb_sparse_mod, only : psb_d_sparse_mat, psb_desc_type, psb_dpk_
    use mld_prec_type, only : mld_dprec_type
    type(mld_dprec_type), intent(inout)    :: p
    integer, intent(in)                    :: what 
    character(len=*), intent(in)           :: val
    integer, intent(out)                   :: info

    call mld_inner_precset(p,what,val,info)
  end subroutine mld_i_dprecsetc

!!$  subroutine mld_i_cprecseti(p,what,val,info)
!!$    use psb_sparse_mod, only : psb_c_sparse_mat, psb_desc_type, psb_spk_
!!$    use mld_prec_type, only : mld_cprec_type
!!$    type(mld_cprec_type), intent(inout)    :: p
!!$    integer, intent(in)                    :: what 
!!$    integer, intent(in)                    :: val
!!$    integer, intent(out)                   :: info
!!$
!!$    call mld_inner_precset(p,what,val,info)
!!$  end subroutine mld_i_cprecseti
!!$
!!$  subroutine mld_i_cprecsetr(p,what,val,info)
!!$    use psb_sparse_mod, only : psb_c_sparse_mat, psb_desc_type, psb_spk_
!!$    use mld_prec_type, only : mld_cprec_type
!!$    type(mld_cprec_type), intent(inout)    :: p
!!$    integer, intent(in)                    :: what 
!!$    real(psb_spk_), intent(in)             :: val
!!$    integer, intent(out)                   :: info
!!$
!!$    call mld_inner_precset(p,what,val,info)
!!$  end subroutine mld_i_cprecsetr
!!$
!!$  subroutine mld_i_cprecsetc(p,what,val,info)
!!$    use psb_sparse_mod, only : psb_c_sparse_mat, psb_desc_type, psb_spk_
!!$    use mld_prec_type, only : mld_cprec_type
!!$    type(mld_cprec_type), intent(inout)    :: p
!!$    integer, intent(in)                    :: what 
!!$    character(len=*), intent(in)           :: val
!!$    integer, intent(out)                   :: info
!!$
!!$    call mld_inner_precset(p,what,val,info)
!!$  end subroutine mld_i_cprecsetc
!!$
!!$  subroutine mld_i_zprecseti(p,what,val,info)
!!$    use psb_sparse_mod, only : psb_z_sparse_mat, psb_desc_type, psb_dpk_
!!$    use mld_prec_type, only : mld_zprec_type
!!$    type(mld_zprec_type), intent(inout)    :: p
!!$    integer, intent(in)                    :: what 
!!$    integer, intent(in)                    :: val
!!$    integer, intent(out)                   :: info
!!$
!!$    call mld_inner_precset(p,what,val,info)
!!$  end subroutine mld_i_zprecseti
!!$
!!$  subroutine mld_i_zprecsetr(p,what,val,info)
!!$    use psb_sparse_mod, only : psb_z_sparse_mat, psb_desc_type, psb_dpk_
!!$    use mld_prec_type, only : mld_zprec_type
!!$    type(mld_zprec_type), intent(inout)    :: p
!!$    integer, intent(in)                    :: what 
!!$    real(psb_dpk_), intent(in)             :: val
!!$    integer, intent(out)                   :: info
!!$
!!$    call mld_inner_precset(p,what,val,info)
!!$  end subroutine mld_i_zprecsetr
!!$
!!$  subroutine mld_i_zprecsetc(p,what,val,info)
!!$    use psb_sparse_mod, only : psb_z_sparse_mat, psb_desc_type, psb_dpk_
!!$    use mld_prec_type, only : mld_zprec_type
!!$    type(mld_zprec_type), intent(inout)    :: p
!!$    integer, intent(in)                    :: what 
!!$    character(len=*), intent(in)           :: val
!!$    integer, intent(out)                   :: info
!!$
!!$    call mld_inner_precset(p,what,val,info)
!!$  end subroutine mld_i_zprecsetc
!!$


end module mld_prec_mod
