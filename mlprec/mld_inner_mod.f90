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
! File: mld_inner_mod.f90
!
! Module: mld_inner_mod
!
!  This module defines the interfaces to the real/complex, single/double
!  precision versions of the MLD2P4 routines, except those of the user level,
!  whose interfaces are defined in mld_prec_mod.f90.
!
module mld_inner_mod
  use mld_prec_type
  use mld_move_alloc_mod

  interface mld_baseprec_aply
    subroutine mld_sbaseprec_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_sbaseprec_type), intent(in) :: prec
      real(psb_spk_),intent(in)         :: x(:)
      real(psb_spk_),intent(inout)      :: y(:)
      real(psb_spk_),intent(in)         :: alpha,beta
      character(len=1)                    :: trans
      real(psb_spk_),target             :: work(:)
      integer, intent(out)                :: info
    end subroutine mld_sbaseprec_aply
    subroutine mld_dbaseprec_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_dbaseprec_type), intent(in) :: prec
      real(psb_dpk_),intent(in)         :: x(:)
      real(psb_dpk_),intent(inout)      :: y(:)
      real(psb_dpk_),intent(in)         :: alpha,beta
      character(len=1)                    :: trans
      real(psb_dpk_),target             :: work(:)
      integer, intent(out)                :: info
    end subroutine mld_dbaseprec_aply
    subroutine mld_cbaseprec_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_cbaseprec_type), intent(in) :: prec
      complex(psb_spk_),intent(in)      :: x(:)
      complex(psb_spk_),intent(inout)   :: y(:)
      complex(psb_spk_),intent(in)      :: alpha,beta
      character(len=1)                    :: trans
      complex(psb_spk_),target          :: work(:)
      integer, intent(out)                :: info
    end subroutine mld_cbaseprec_aply
    subroutine mld_zbaseprec_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_zbaseprec_type), intent(in) :: prec
      complex(psb_dpk_),intent(in)      :: x(:)
      complex(psb_dpk_),intent(inout)   :: y(:)
      complex(psb_dpk_),intent(in)      :: alpha,beta
      character(len=1)                    :: trans
      complex(psb_dpk_),target          :: work(:)
      integer, intent(out)                :: info
    end subroutine mld_zbaseprec_aply
  end interface

  interface mld_mlprec_bld
    subroutine mld_smlprec_bld(a,desc_a,prec,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sprec_type
      implicit none
      type(psb_sspmat_type), intent(in), target   :: a
      type(psb_desc_type), intent(in), target     :: desc_a
      type(mld_sprec_type), intent(inout), target :: prec
      integer, intent(out)                        :: info
!!$      character, intent(in),optional             :: upd
    end subroutine mld_smlprec_bld
    subroutine mld_dmlprec_bld(a,desc_a,prec,info)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dprec_type
      implicit none
      type(psb_dspmat_type), intent(in), target   :: a
      type(psb_desc_type), intent(in), target     :: desc_a
      type(mld_dprec_type), intent(inout), target :: prec
      integer, intent(out)                        :: info
!!$      character, intent(in),optional             :: upd
    end subroutine mld_dmlprec_bld
    subroutine mld_cmlprec_bld(a,desc_a,prec,info)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cprec_type
      implicit none
      type(psb_cspmat_type), intent(in), target   :: a
      type(psb_desc_type), intent(in), target     :: desc_a
      type(mld_cprec_type), intent(inout), target :: prec
      integer, intent(out)                        :: info
!!$      character, intent(in),optional             :: upd
    end subroutine mld_cmlprec_bld
    subroutine mld_zmlprec_bld(a,desc_a,prec,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zprec_type
      implicit none
      type(psb_zspmat_type), intent(in), target  :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(mld_zprec_type), intent(inout)        :: prec
      integer, intent(out)                       :: info
!!$      character, intent(in),optional             :: upd
    end subroutine mld_zmlprec_bld
  end interface

  interface mld_as_aply
    subroutine mld_sas_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_sbaseprec_type), intent(in) :: prec
      real(psb_spk_),intent(in)         :: x(:)
      real(psb_spk_),intent(inout)      :: y(:)
      real(psb_spk_),intent(in)         :: alpha,beta
      character(len=1)                    :: trans
      real(psb_spk_),target             :: work(:)
      integer, intent(out)                :: info
    end subroutine mld_sas_aply
    subroutine mld_das_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_dbaseprec_type), intent(in) :: prec
      real(psb_dpk_),intent(in)         :: x(:)
      real(psb_dpk_),intent(inout)      :: y(:)
      real(psb_dpk_),intent(in)         :: alpha,beta
      character(len=1)                    :: trans
      real(psb_dpk_),target             :: work(:)
      integer, intent(out)                :: info
    end subroutine mld_das_aply
    subroutine mld_cas_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_cbaseprec_type), intent(in) :: prec
      complex(psb_spk_),intent(in)      :: x(:)
      complex(psb_spk_),intent(inout)   :: y(:)
      complex(psb_spk_),intent(in)      :: alpha,beta
      character(len=1)                    :: trans
      complex(psb_spk_),target          :: work(:)
      integer, intent(out)                :: info
    end subroutine mld_cas_aply
    subroutine mld_zas_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_zbaseprec_type), intent(in) :: prec
      complex(psb_dpk_),intent(in)      :: x(:)
      complex(psb_dpk_),intent(inout)   :: y(:)
      complex(psb_dpk_),intent(in)      :: alpha,beta
      character(len=1)                    :: trans
      complex(psb_dpk_),target          :: work(:)
      integer, intent(out)                :: info
    end subroutine mld_zas_aply
  end interface

  interface mld_mlprec_aply
    subroutine mld_smlprec_aply(alpha,p,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type, mld_sprec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_sprec_type), intent(in)  :: p
      real(psb_spk_),intent(in)         :: alpha,beta
      real(psb_spk_),intent(in)         :: x(:)
      real(psb_spk_),intent(inout)      :: y(:)
      character,intent(in)              :: trans
      real(psb_spk_),target             :: work(:)
      integer, intent(out)              :: info
    end subroutine mld_smlprec_aply
    subroutine mld_dmlprec_aply(alpha,p,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type, mld_dprec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(mld_dprec_type), intent(in)  :: p
      real(psb_dpk_),intent(in)         :: alpha,beta
      real(psb_dpk_),intent(in)         :: x(:)
      real(psb_dpk_),intent(inout)      :: y(:)
      character,intent(in)              :: trans
      real(psb_dpk_),target             :: work(:)
      integer, intent(out)              :: info
    end subroutine mld_dmlprec_aply
    subroutine mld_cmlprec_aply(alpha,p,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type, mld_cprec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(mld_cprec_type), intent(in)  :: p
      complex(psb_spk_),intent(in)      :: alpha,beta
      complex(psb_spk_),intent(in)      :: x(:)
      complex(psb_spk_),intent(inout)   :: y(:)
      character,intent(in)              :: trans
      complex(psb_spk_),target          :: work(:)
      integer, intent(out)              :: info
    end subroutine mld_cmlprec_aply
    subroutine mld_zmlprec_aply(alpha,p,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type, mld_zprec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(mld_zprec_type), intent(in)  :: p
      complex(psb_dpk_),intent(in)      :: alpha,beta
      complex(psb_dpk_),intent(in)      :: x(:)
      complex(psb_dpk_),intent(inout)   :: y(:)
      character,intent(in)              :: trans
      complex(psb_dpk_),target          :: work(:)
      integer, intent(out)              :: info
    end subroutine mld_zmlprec_aply
  end interface


  interface mld_sub_aply
    subroutine mld_ssub_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type
      type(psb_desc_type), intent(in)       :: desc_data
      type(mld_sbaseprec_type), intent(in)   :: prec
      real(psb_spk_),intent(in)           :: x(:)
      real(psb_spk_),intent(inout)        :: y(:)
      real(psb_spk_),intent(in)           :: alpha,beta
      character(len=1),intent(in)           :: trans
      real(psb_spk_),target,intent(inout) :: work(:)
      integer, intent(out)                  :: info
    end subroutine mld_ssub_aply
    subroutine mld_dsub_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type
      type(psb_desc_type), intent(in)       :: desc_data
      type(mld_dbaseprec_type), intent(in)   :: prec
      real(psb_dpk_),intent(in)           :: x(:)
      real(psb_dpk_),intent(inout)        :: y(:)
      real(psb_dpk_),intent(in)           :: alpha,beta
      character(len=1),intent(in)           :: trans
      real(psb_dpk_),target,intent(inout) :: work(:)
      integer, intent(out)                  :: info
    end subroutine mld_dsub_aply
    subroutine mld_csub_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type
      type(psb_desc_type), intent(in)          :: desc_data
      type(mld_cbaseprec_type), intent(in)      :: prec
      complex(psb_spk_),intent(in)           :: x(:)
      complex(psb_spk_),intent(inout)        :: y(:)
      complex(psb_spk_),intent(in)           :: alpha,beta
      character(len=1),intent(in)              :: trans
      complex(psb_spk_),target,intent(inout) :: work(:)
      integer, intent(out)                     :: info
    end subroutine mld_csub_aply
    subroutine mld_zsub_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type
      type(psb_desc_type), intent(in)          :: desc_data
      type(mld_zbaseprec_type), intent(in)      :: prec
      complex(psb_dpk_),intent(in)           :: x(:)
      complex(psb_dpk_),intent(inout)        :: y(:)
      complex(psb_dpk_),intent(in)           :: alpha,beta
      character(len=1),intent(in)              :: trans
      complex(psb_dpk_),target,intent(inout) :: work(:)
      integer, intent(out)                     :: info
    end subroutine mld_zsub_aply
  end interface


  interface mld_sub_solve
    subroutine mld_ssub_solve(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type
      type(psb_desc_type), intent(in)       :: desc_data
      type(mld_sbaseprec_type), intent(in)   :: prec
      real(psb_spk_),intent(in)           :: x(:)
      real(psb_spk_),intent(inout)        :: y(:)
      real(psb_spk_),intent(in)           :: alpha,beta
      character(len=1),intent(in)           :: trans
      real(psb_spk_),target,intent(inout) :: work(:)
      integer, intent(out)                  :: info
    end subroutine mld_ssub_solve
    subroutine mld_dsub_solve(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type
      type(psb_desc_type), intent(in)       :: desc_data
      type(mld_dbaseprec_type), intent(in)   :: prec
      real(psb_dpk_),intent(in)           :: x(:)
      real(psb_dpk_),intent(inout)        :: y(:)
      real(psb_dpk_),intent(in)           :: alpha,beta
      character(len=1),intent(in)           :: trans
      real(psb_dpk_),target,intent(inout) :: work(:)
      integer, intent(out)                  :: info
    end subroutine mld_dsub_solve
    subroutine mld_csub_solve(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type
      type(psb_desc_type), intent(in)          :: desc_data
      type(mld_cbaseprec_type), intent(in)      :: prec
      complex(psb_spk_),intent(in)           :: x(:)
      complex(psb_spk_),intent(inout)        :: y(:)
      complex(psb_spk_),intent(in)           :: alpha,beta
      character(len=1),intent(in)              :: trans
      complex(psb_spk_),target,intent(inout) :: work(:)
      integer, intent(out)                     :: info
    end subroutine mld_csub_solve
    subroutine mld_zsub_solve(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type
      type(psb_desc_type), intent(in)          :: desc_data
      type(mld_zbaseprec_type), intent(in)      :: prec
      complex(psb_dpk_),intent(in)           :: x(:)
      complex(psb_dpk_),intent(inout)        :: y(:)
      complex(psb_dpk_),intent(in)           :: alpha,beta
      character(len=1),intent(in)              :: trans
      complex(psb_dpk_),target,intent(inout) :: work(:)
      integer, intent(out)                     :: info
    end subroutine mld_zsub_solve
  end interface


  interface mld_asmat_bld
    Subroutine mld_sasmat_bld(ptype,novr,a,blk,desc_data,upd,desc_p,info,outfmt)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type
      integer, intent(in)                 :: ptype,novr
      Type(psb_sspmat_type), Intent(in)   ::  a
      Type(psb_sspmat_type), Intent(out)  ::  blk
      Type(psb_desc_type), Intent(inout)  :: desc_p
      Type(psb_desc_type), Intent(in)     :: desc_data 
      Character, Intent(in)               :: upd
      integer, intent(out)                :: info
      character(len=5), optional          :: outfmt
    end Subroutine mld_sasmat_bld
    Subroutine mld_dasmat_bld(ptype,novr,a,blk,desc_data,upd,desc_p,info,outfmt)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type
      integer, intent(in)                 :: ptype,novr
      Type(psb_dspmat_type), Intent(in)   ::  a
      Type(psb_dspmat_type), Intent(out)  ::  blk
      Type(psb_desc_type), Intent(inout)  :: desc_p
      Type(psb_desc_type), Intent(in)     :: desc_data 
      Character, Intent(in)               :: upd
      integer, intent(out)                :: info
      character(len=5), optional          :: outfmt
    end Subroutine mld_dasmat_bld
    Subroutine mld_casmat_bld(ptype,novr,a,blk,desc_data,upd,desc_p,info,outfmt)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type
      integer, intent(in)                 :: ptype,novr
      Type(psb_cspmat_type), Intent(in)   ::  a
      Type(psb_cspmat_type), Intent(out)  ::  blk
      Type(psb_desc_type), Intent(inout)  :: desc_p
      Type(psb_desc_type), Intent(in)     :: desc_data 
      Character, Intent(in)               :: upd
      integer, intent(out)                :: info
      character(len=5), optional          :: outfmt
    end Subroutine mld_casmat_bld
    Subroutine mld_zasmat_bld(ptype,novr,a,blk,desc_data,upd,desc_p,info,outfmt)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type
      integer, intent(in)                 :: ptype,novr
      Type(psb_zspmat_type), Intent(in)   ::  a
      Type(psb_zspmat_type), Intent(out)  ::  blk
      Type(psb_desc_type), Intent(inout)  :: desc_p
      Type(psb_desc_type), Intent(in)     :: desc_data 
      Character, Intent(in)               :: upd
      integer, intent(out)                :: info
      character(len=5), optional          :: outfmt
    end Subroutine mld_zasmat_bld
  end interface

  interface mld_sp_renum
    subroutine mld_ssp_renum(a,blck,p,atmp,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type
      type(psb_sspmat_type), intent(in)      :: a,blck
      type(psb_sspmat_type), intent(out)     :: atmp
      type(mld_sbaseprec_type), intent(inout) :: p
      integer, intent(out)   :: info
    end subroutine mld_ssp_renum
    subroutine mld_dsp_renum(a,blck,p,atmp,info)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type
      type(psb_dspmat_type), intent(in)      :: a,blck
      type(psb_dspmat_type), intent(out)     :: atmp
      type(mld_dbaseprec_type), intent(inout) :: p
      integer, intent(out)   :: info
    end subroutine mld_dsp_renum
    subroutine mld_csp_renum(a,blck,p,atmp,info)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type
      type(psb_cspmat_type), intent(in)      :: a,blck
      type(psb_cspmat_type), intent(out)     :: atmp
      type(mld_cbaseprec_type), intent(inout) :: p
      integer, intent(out)   :: info
    end subroutine mld_csp_renum
    subroutine mld_zsp_renum(a,blck,p,atmp,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type
      type(psb_zspmat_type), intent(in)      :: a,blck
      type(psb_zspmat_type), intent(out)     :: atmp
      type(mld_zbaseprec_type), intent(inout) :: p
      integer, intent(out)   :: info
    end subroutine mld_zsp_renum
  end interface

  interface mld_coarse_bld
    subroutine mld_scoarse_bld(a,desc_a,p,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type, mld_sonelev_type
      type(psb_sspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      type(mld_sonelev_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_scoarse_bld
    subroutine mld_dcoarse_bld(a,desc_a,p,info)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type, mld_donelev_type
      type(psb_dspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      type(mld_donelev_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_dcoarse_bld
    subroutine mld_ccoarse_bld(a,desc_a,p,info)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type, mld_conelev_type
      type(psb_cspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      type(mld_conelev_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_ccoarse_bld
    subroutine mld_zcoarse_bld(a,desc_a,p,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type, mld_zonelev_type
      type(psb_zspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      type(mld_zonelev_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_zcoarse_bld
  end interface

  interface mld_aggrmap_bld
    subroutine mld_saggrmap_bld(aggr_type,theta,a,desc_a,ilaggr,nlaggr,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type
      integer, intent(in)               :: aggr_type
      real(psb_spk_), intent(in)        :: theta
      type(psb_sspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer, allocatable, intent(out) :: ilaggr(:),nlaggr(:)
      integer, intent(out)              :: info
    end subroutine mld_saggrmap_bld
    subroutine mld_daggrmap_bld(aggr_type,theta,a,desc_a,ilaggr,nlaggr,info)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type
      integer, intent(in)               :: aggr_type
      real(psb_dpk_), intent(in)        :: theta
      type(psb_dspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer, allocatable, intent(out) :: ilaggr(:),nlaggr(:)
      integer, intent(out)              :: info
    end subroutine mld_daggrmap_bld
    subroutine mld_caggrmap_bld(aggr_type,theta,a,desc_a,ilaggr,nlaggr,info)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type
      integer, intent(in)               :: aggr_type
      real(psb_spk_), intent(in)        :: theta
      type(psb_cspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer, allocatable, intent(out) :: ilaggr(:),nlaggr(:)
      integer, intent(out)              :: info
    end subroutine mld_caggrmap_bld
    subroutine mld_zaggrmap_bld(aggr_type,theta,a,desc_a,ilaggr,nlaggr,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type
      integer, intent(in)               :: aggr_type
      real(psb_dpk_), intent(in)        :: theta
      type(psb_zspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer, allocatable, intent(out) :: ilaggr(:),nlaggr(:)
      integer, intent(out)              :: info
    end subroutine mld_zaggrmap_bld
  end interface

  interface mld_aggrmat_asb
    subroutine mld_saggrmat_asb(a,desc_a,ilaggr,nlaggr,p,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type, mld_sonelev_type
      type(psb_sspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      integer, intent(inout)                         :: ilaggr(:), nlaggr(:)
      type(mld_sonelev_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_saggrmat_asb
    subroutine mld_daggrmat_asb(a,desc_a,ilaggr,nlaggr,p,info)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type, mld_donelev_type
      type(psb_dspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      integer, intent(inout)                         :: ilaggr(:), nlaggr(:)
      type(mld_donelev_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_daggrmat_asb
    subroutine mld_caggrmat_asb(a,desc_a,ilaggr,nlaggr,p,info)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type, mld_conelev_type
      type(psb_cspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      integer, intent(inout)                         :: ilaggr(:), nlaggr(:)
      type(mld_conelev_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_caggrmat_asb
    subroutine mld_zaggrmat_asb(a,desc_a,ilaggr,nlaggr,p,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type, mld_zonelev_type
      type(psb_zspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      integer, intent(inout)                         :: ilaggr(:), nlaggr(:)
      type(mld_zonelev_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_zaggrmat_asb
  end interface

  interface mld_aggrmat_nosmth_asb
    subroutine mld_saggrmat_nosmth_asb(a,desc_a,ilaggr,nlaggr,p,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type, mld_sonelev_type
      type(psb_sspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      integer, intent(inout)                         :: ilaggr(:), nlaggr(:)
      type(mld_sonelev_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_saggrmat_nosmth_asb
    subroutine mld_daggrmat_nosmth_asb(a,desc_a,ilaggr,nlaggr,p,info)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type, mld_donelev_type
      type(psb_dspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      integer, intent(inout)                         :: ilaggr(:), nlaggr(:)
      type(mld_donelev_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_daggrmat_nosmth_asb
    subroutine mld_caggrmat_nosmth_asb(a,desc_a,ilaggr,nlaggr,p,info)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type, mld_conelev_type
      type(psb_cspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      integer, intent(inout)                         :: ilaggr(:), nlaggr(:)
      type(mld_conelev_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_caggrmat_nosmth_asb
    subroutine mld_zaggrmat_nosmth_asb(a,desc_a,ilaggr,nlaggr,p,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type, mld_zonelev_type
      type(psb_zspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      integer, intent(inout)                         :: ilaggr(:), nlaggr(:)
      type(mld_zonelev_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_zaggrmat_nosmth_asb
  end interface

  interface mld_aggrmat_smth_asb
    subroutine mld_saggrmat_smth_asb(a,desc_a,ilaggr,nlaggr,p,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type, mld_sonelev_type
      type(psb_sspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      integer, intent(inout)                         :: ilaggr(:), nlaggr(:)
      type(mld_sonelev_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_saggrmat_smth_asb
    subroutine mld_daggrmat_smth_asb(a,desc_a,ilaggr,nlaggr,p,info)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type, mld_donelev_type
      type(psb_dspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      integer, intent(inout)                         :: ilaggr(:), nlaggr(:)
      type(mld_donelev_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_daggrmat_smth_asb
    subroutine mld_caggrmat_smth_asb(a,desc_a,ilaggr,nlaggr,p,info)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type, mld_conelev_type
      type(psb_cspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      integer, intent(inout)                         :: ilaggr(:), nlaggr(:)
      type(mld_conelev_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_caggrmat_smth_asb
    subroutine mld_zaggrmat_smth_asb(a,desc_a,ilaggr,nlaggr,p,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type, mld_zonelev_type
      type(psb_zspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      integer, intent(inout)                         :: ilaggr(:), nlaggr(:)
      type(mld_zonelev_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_zaggrmat_smth_asb
  end interface

  interface mld_baseprec_bld
    subroutine mld_sbaseprec_bld(a,desc_a,p,info,upd)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type
      type(psb_sspmat_type), target              :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(mld_sbaseprec_type),intent(inout)      :: p
      integer, intent(out)                       :: info
      character, intent(in), optional            :: upd
    end subroutine mld_sbaseprec_bld
    subroutine mld_dbaseprec_bld(a,desc_a,p,info,upd)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type
      type(psb_dspmat_type), target              :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(mld_dbaseprec_type),intent(inout)      :: p
      integer, intent(out)                       :: info
      character, intent(in), optional            :: upd
    end subroutine mld_dbaseprec_bld
    subroutine mld_cbaseprec_bld(a,desc_a,p,info,upd)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type
      type(psb_cspmat_type), target              :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(mld_cbaseprec_type),intent(inout)      :: p
      integer, intent(out)                       :: info
      character, intent(in), optional            :: upd
    end subroutine mld_cbaseprec_bld
    subroutine mld_zbaseprec_bld(a,desc_a,p,info,upd)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type
      type(psb_zspmat_type), target              :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(mld_zbaseprec_type),intent(inout)      :: p
      integer, intent(out)                       :: info
      character, intent(in), optional            :: upd
    end subroutine mld_zbaseprec_bld
  end interface

  interface mld_as_bld
    subroutine mld_sas_bld(a,desc_a,p,upd,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type
      type(psb_sspmat_type),intent(in), target :: a
      type(psb_desc_type), intent(in), target  :: desc_a
      type(mld_sbaseprec_type),intent(inout)   :: p
      character, intent(in)                    :: upd
      integer, intent(out)                     :: info
    end subroutine mld_sas_bld
    subroutine mld_das_bld(a,desc_a,p,upd,info)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type
      type(psb_dspmat_type),intent(in), target :: a
      type(psb_desc_type), intent(in), target  :: desc_a
      type(mld_dbaseprec_type),intent(inout)   :: p
      character, intent(in)                    :: upd
      integer, intent(out)                     :: info
    end subroutine mld_das_bld
    subroutine mld_cas_bld(a,desc_a,p,upd,info)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type
      type(psb_cspmat_type),intent(in), target :: a
      type(psb_desc_type), intent(in), target  :: desc_a
      type(mld_cbaseprec_type),intent(inout)   :: p
      character, intent(in)                    :: upd
      integer, intent(out)                     :: info
    end subroutine mld_cas_bld
    subroutine mld_zas_bld(a,desc_a,p,upd,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type
      type(psb_zspmat_type),intent(in), target :: a
      type(psb_desc_type), intent(in), target  :: desc_a
      type(mld_zbaseprec_type),intent(inout)   :: p
      character, intent(in)                    :: upd
      integer, intent(out)                     :: info
    end subroutine mld_zas_bld
  end interface

  interface mld_diag_bld
    subroutine mld_sdiag_bld(a,desc_data,p,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type
      integer, intent(out) :: info
      type(psb_sspmat_type), intent(in), target :: a
      type(psb_desc_type),intent(in)            :: desc_data
      type(mld_sbaseprec_type), intent(inout)    :: p
    end subroutine mld_sdiag_bld
    subroutine mld_ddiag_bld(a,desc_data,p,info)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type
      integer, intent(out) :: info
      type(psb_dspmat_type), intent(in), target :: a
      type(psb_desc_type),intent(in)            :: desc_data
      type(mld_dbaseprec_type), intent(inout)    :: p
    end subroutine mld_ddiag_bld
    subroutine mld_cdiag_bld(a,desc_data,p,info)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type
      integer, intent(out) :: info
      type(psb_cspmat_type), intent(in), target :: a
      type(psb_desc_type),intent(in)            :: desc_data
      type(mld_cbaseprec_type), intent(inout)    :: p
    end subroutine mld_cdiag_bld
    subroutine mld_zdiag_bld(a,desc_data,p,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type
      integer, intent(out) :: info
      type(psb_zspmat_type), intent(in), target :: a
      type(psb_desc_type),intent(in)            :: desc_data
      type(mld_zbaseprec_type), intent(inout)    :: p
    end subroutine mld_zdiag_bld
  end interface

  interface mld_fact_bld
    subroutine mld_sfact_bld(a,p,upd,info,blck)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type
      type(psb_sspmat_type), intent(in), target :: a
      type(mld_sbaseprec_type), intent(inout)    :: p
      integer, intent(out)                      :: info
      character, intent(in)                     :: upd
      type(psb_sspmat_type), intent(in), target, optional  :: blck
    end subroutine mld_sfact_bld
    subroutine mld_dfact_bld(a,p,upd,info,blck)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type
      type(psb_dspmat_type), intent(in), target :: a
      type(mld_dbaseprec_type), intent(inout)    :: p
      integer, intent(out)                      :: info
      character, intent(in)                     :: upd
      type(psb_dspmat_type), intent(in), target, optional  :: blck
    end subroutine mld_dfact_bld
    subroutine mld_cfact_bld(a,p,upd,info,blck)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type
      type(psb_cspmat_type), intent(in), target :: a
      type(mld_cbaseprec_type), intent(inout)    :: p
      integer, intent(out)                      :: info
      character, intent(in)                     :: upd
      type(psb_cspmat_type), intent(in), target, optional  :: blck
    end subroutine mld_cfact_bld
    subroutine mld_zfact_bld(a,p,upd,info,blck)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type
      type(psb_zspmat_type), intent(in), target :: a
      type(mld_zbaseprec_type), intent(inout)    :: p
      integer, intent(out)                      :: info
      character, intent(in)                     :: upd
      type(psb_zspmat_type), intent(in), target, optional  :: blck
    end subroutine mld_zfact_bld
  end interface

  interface mld_ilu_bld
    subroutine mld_silu_bld(a,p,upd,info,blck)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type
      integer, intent(out) :: info
      type(psb_sspmat_type), intent(in), target :: a
      type(mld_sbaseprec_type), intent(inout)    :: p
      character, intent(in)                     :: upd
      type(psb_sspmat_type), intent(in), optional :: blck
    end subroutine mld_silu_bld
    subroutine mld_dilu_bld(a,p,upd,info,blck)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type
      integer, intent(out) :: info
      type(psb_dspmat_type), intent(in), target :: a
      type(mld_dbaseprec_type), intent(inout)    :: p
      character, intent(in)                     :: upd
      type(psb_dspmat_type), intent(in), optional :: blck
    end subroutine mld_dilu_bld
    subroutine mld_cilu_bld(a,p,upd,info,blck)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type
      integer, intent(out) :: info
      type(psb_cspmat_type), intent(in), target :: a
      type(mld_cbaseprec_type), intent(inout)    :: p
      character, intent(in)                     :: upd
      type(psb_cspmat_type), intent(in), optional :: blck
    end subroutine mld_cilu_bld
    subroutine mld_zilu_bld(a,p,upd,info,blck)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type
      integer, intent(out) :: info
      type(psb_zspmat_type), intent(in), target :: a
      type(mld_zbaseprec_type), intent(inout)    :: p
      character, intent(in)                     :: upd
      type(psb_zspmat_type), intent(in), optional :: blck
    end subroutine mld_zilu_bld
  end interface

  interface mld_sludist_bld
    subroutine mld_ssludist_bld(a,desc_a,p,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type
      type(psb_sspmat_type), intent(inout)   :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_sbaseprec_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_ssludist_bld
    subroutine mld_dsludist_bld(a,desc_a,p,info)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type
      type(psb_dspmat_type), intent(inout)   :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_dbaseprec_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_dsludist_bld
    subroutine mld_csludist_bld(a,desc_a,p,info)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type
      type(psb_cspmat_type), intent(inout)   :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_cbaseprec_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_csludist_bld
    subroutine mld_zsludist_bld(a,desc_a,p,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type
      type(psb_zspmat_type), intent(inout)   :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_zbaseprec_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_zsludist_bld
  end interface

  interface mld_slu_bld
    subroutine mld_sslu_bld(a,desc_a,p,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type
      type(psb_sspmat_type), intent(inout)      :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_sbaseprec_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_sslu_bld
    subroutine mld_dslu_bld(a,desc_a,p,info)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type
      type(psb_dspmat_type), intent(inout)      :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_dbaseprec_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_dslu_bld
    subroutine mld_cslu_bld(a,desc_a,p,info)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type
      type(psb_cspmat_type), intent(inout)   :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_cbaseprec_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_cslu_bld
    subroutine mld_zslu_bld(a,desc_a,p,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type
      type(psb_zspmat_type), intent(inout)   :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_zbaseprec_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_zslu_bld
  end interface

  interface mld_umf_bld
    subroutine mld_sumf_bld(a,desc_a,p,info)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_sbaseprec_type
      type(psb_sspmat_type), intent(inout)    :: a
      type(psb_desc_type), intent(in)         :: desc_a
      type(mld_sbaseprec_type), intent(inout) :: p
      integer, intent(out)                    :: info
    end subroutine mld_sumf_bld
    subroutine mld_dumf_bld(a,desc_a,p,info)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_dbaseprec_type
      type(psb_dspmat_type), intent(inout)    :: a
      type(psb_desc_type), intent(in)         :: desc_a
      type(mld_dbaseprec_type), intent(inout) :: p
      integer, intent(out)                    :: info
    end subroutine mld_dumf_bld
    subroutine mld_cumf_bld(a,desc_a,p,info)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      use mld_prec_type, only : mld_cbaseprec_type
      type(psb_cspmat_type), intent(inout)    :: a
      type(psb_desc_type), intent(in)         :: desc_a
      type(mld_cbaseprec_type), intent(inout) :: p
      integer, intent(out)                    :: info
    end subroutine mld_cumf_bld
    subroutine mld_zumf_bld(a,desc_a,p,info)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      use mld_prec_type, only : mld_zbaseprec_type
      type(psb_zspmat_type), intent(inout)    :: a
      type(psb_desc_type), intent(in)         :: desc_a
      type(mld_zbaseprec_type), intent(inout) :: p
      integer, intent(out)                    :: info
    end subroutine mld_zumf_bld
  end interface

  interface mld_ilu0_fact
    subroutine mld_silu0_fact(ialg,a,l,u,d,info,blck)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      integer, intent(in)                 :: ialg
      integer, intent(out)                :: info
      type(psb_sspmat_type),intent(in)    :: a
      type(psb_sspmat_type),intent(inout) :: l,u
      type(psb_sspmat_type),intent(in), optional, target :: blck
      real(psb_spk_), intent(inout)     ::  d(:)
    end subroutine mld_silu0_fact
    subroutine mld_dilu0_fact(ialg,a,l,u,d,info,blck)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      integer, intent(in)                 :: ialg
      integer, intent(out)                :: info
      type(psb_dspmat_type),intent(in)    :: a
      type(psb_dspmat_type),intent(inout) :: l,u
      type(psb_dspmat_type),intent(in), optional, target :: blck
      real(psb_dpk_), intent(inout)     ::  d(:)
    end subroutine mld_dilu0_fact
    subroutine mld_cilu0_fact(ialg,a,l,u,d,info,blck)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      integer, intent(in)                 :: ialg
      integer, intent(out)                :: info
      type(psb_cspmat_type),intent(in)    :: a
      type(psb_cspmat_type),intent(inout) :: l,u
      type(psb_cspmat_type),intent(in), optional, target :: blck
      complex(psb_spk_), intent(inout)     ::  d(:)
    end subroutine mld_cilu0_fact
    subroutine mld_zilu0_fact(ialg,a,l,u,d,info,blck)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      integer, intent(in)                 :: ialg
      integer, intent(out)                :: info
      type(psb_zspmat_type),intent(in)    :: a
      type(psb_zspmat_type),intent(inout) :: l,u
      type(psb_zspmat_type),intent(in), optional, target :: blck
      complex(psb_dpk_), intent(inout)     ::  d(:)
    end subroutine mld_zilu0_fact
  end interface

  interface mld_iluk_fact
    subroutine mld_siluk_fact(fill_in,ialg,a,l,u,d,info,blck)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      integer, intent(in)                 :: fill_in,ialg
      integer, intent(out)                :: info
      type(psb_sspmat_type),intent(in)    :: a
      type(psb_sspmat_type),intent(inout) :: l,u
      type(psb_sspmat_type),intent(in), optional, target :: blck
      real(psb_spk_), intent(inout)     ::  d(:)
    end subroutine mld_siluk_fact
    subroutine mld_diluk_fact(fill_in,ialg,a,l,u,d,info,blck)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      integer, intent(in)                 :: fill_in,ialg
      integer, intent(out)                :: info
      type(psb_dspmat_type),intent(in)    :: a
      type(psb_dspmat_type),intent(inout) :: l,u
      type(psb_dspmat_type),intent(in), optional, target :: blck
      real(psb_dpk_), intent(inout)     ::  d(:)
    end subroutine mld_diluk_fact
    subroutine mld_ciluk_fact(fill_in,ialg,a,l,u,d,info,blck)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      integer, intent(in)                 :: fill_in,ialg
      integer, intent(out)                :: info
      type(psb_cspmat_type),intent(in)    :: a
      type(psb_cspmat_type),intent(inout) :: l,u
      type(psb_cspmat_type),intent(in), optional, target :: blck
      complex(psb_spk_), intent(inout)     ::  d(:)
    end subroutine mld_ciluk_fact
    subroutine mld_ziluk_fact(fill_in,ialg,a,l,u,d,info,blck)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      integer, intent(in)                 :: fill_in,ialg
      integer, intent(out)                :: info
      type(psb_zspmat_type),intent(in)    :: a
      type(psb_zspmat_type),intent(inout) :: l,u
      type(psb_zspmat_type),intent(in), optional, target :: blck
      complex(psb_dpk_), intent(inout)     ::  d(:)
    end subroutine mld_ziluk_fact
  end interface

  interface mld_ilut_fact
    subroutine mld_silut_fact(fill_in,thres,a,l,u,d,info,blck)
      use psb_base_mod, only : psb_sspmat_type, psb_desc_type, psb_spk_
      integer, intent(in)                 :: fill_in
      real(psb_spk_), intent(in)        :: thres
      integer, intent(out)                :: info
      type(psb_sspmat_type),intent(in)    :: a
      type(psb_sspmat_type),intent(inout) :: l,u
      type(psb_sspmat_type),intent(in), optional, target :: blck
      real(psb_spk_), intent(inout)     ::  d(:)
    end subroutine mld_silut_fact
    subroutine mld_dilut_fact(fill_in,thres,a,l,u,d,info,blck)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      integer, intent(in)                 :: fill_in
      real(psb_dpk_), intent(in)        :: thres
      integer, intent(out)                :: info
      type(psb_dspmat_type),intent(in)    :: a
      type(psb_dspmat_type),intent(inout) :: l,u
      type(psb_dspmat_type),intent(in), optional, target :: blck
      real(psb_dpk_), intent(inout)     ::  d(:)
    end subroutine mld_dilut_fact
    subroutine mld_cilut_fact(fill_in,thres,a,l,u,d,info,blck)
      use psb_base_mod, only : psb_cspmat_type, psb_desc_type, psb_spk_
      integer, intent(in)                 :: fill_in
      real(psb_spk_), intent(in)        :: thres
      integer, intent(out)                :: info
      type(psb_cspmat_type),intent(in)    :: a
      type(psb_cspmat_type),intent(inout) :: l,u
      type(psb_cspmat_type),intent(in), optional, target :: blck
      complex(psb_spk_), intent(inout)  ::  d(:)
    end subroutine mld_cilut_fact
    subroutine mld_zilut_fact(fill_in,thres,a,l,u,d,info,blck)
      use psb_base_mod, only : psb_zspmat_type, psb_desc_type, psb_dpk_
      integer, intent(in)                 :: fill_in
      real(psb_dpk_), intent(in)        :: thres
      integer, intent(out)                :: info
      type(psb_zspmat_type),intent(in)    :: a
      type(psb_zspmat_type),intent(inout) :: l,u
      type(psb_zspmat_type),intent(in), optional, target :: blck
      complex(psb_dpk_), intent(inout)  ::  d(:)
    end subroutine mld_zilut_fact
  end interface

end module mld_inner_mod
