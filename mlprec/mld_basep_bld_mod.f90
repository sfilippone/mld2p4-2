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
!
! mld_inner_mod is split in two because otherwise compilation blows
! up on XL Fortran 10.0. 
!
!
module mld_basep_bld_mod
  interface mld_baseprc_bld
    subroutine mld_sbaseprc_bld(a,desc_a,p,info,upd)
      use psb_base_mod
      use mld_prec_type
      type(psb_sspmat_type), target              :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(mld_sbaseprc_type),intent(inout)      :: p
      integer, intent(out)                       :: info
      character, intent(in), optional            :: upd
    end subroutine mld_sbaseprc_bld
    subroutine mld_dbaseprc_bld(a,desc_a,p,info,upd)
      use psb_base_mod
      use mld_prec_type
      type(psb_dspmat_type), target              :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(mld_dbaseprc_type),intent(inout)      :: p
      integer, intent(out)                       :: info
      character, intent(in), optional            :: upd
    end subroutine mld_dbaseprc_bld
    subroutine mld_cbaseprc_bld(a,desc_a,p,info,upd)
      use psb_base_mod
      use mld_prec_type
      type(psb_cspmat_type), target              :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(mld_cbaseprc_type),intent(inout)      :: p
      integer, intent(out)                       :: info
      character, intent(in), optional            :: upd
    end subroutine mld_cbaseprc_bld
    subroutine mld_zbaseprc_bld(a,desc_a,p,info,upd)
      use psb_base_mod
      use mld_prec_type
      type(psb_zspmat_type), target              :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(mld_zbaseprc_type),intent(inout)      :: p
      integer, intent(out)                       :: info
      character, intent(in), optional            :: upd
    end subroutine mld_zbaseprc_bld
  end interface

  interface mld_as_bld
    subroutine mld_sas_bld(a,desc_a,p,upd,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_sspmat_type), target           :: a
      type(psb_desc_type), intent(in), target :: desc_a
      type(mld_sbaseprc_type),intent(inout)   :: p
      character, intent(in)                   :: upd
      integer, intent(out)                    :: info
    end subroutine mld_sas_bld
    subroutine mld_das_bld(a,desc_a,p,upd,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_dspmat_type), target           :: a
      type(psb_desc_type), intent(in), target :: desc_a
      type(mld_dbaseprc_type),intent(inout)   :: p
      character, intent(in)                   :: upd
      integer, intent(out)                    :: info
    end subroutine mld_das_bld
    subroutine mld_cas_bld(a,desc_a,p,upd,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_cspmat_type), target           :: a
      type(psb_desc_type), intent(in), target :: desc_a
      type(mld_cbaseprc_type),intent(inout)   :: p
      character, intent(in)                   :: upd
      integer, intent(out)                    :: info
    end subroutine mld_cas_bld
    subroutine mld_zas_bld(a,desc_a,p,upd,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_zspmat_type), target           :: a
      type(psb_desc_type), intent(in), target :: desc_a
      type(mld_zbaseprc_type),intent(inout)   :: p
      character, intent(in)                   :: upd
      integer, intent(out)                    :: info
    end subroutine mld_zas_bld
  end interface

  interface mld_mlprec_bld
    subroutine mld_smlprec_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_sspmat_type), intent(inout), target :: a
      type(psb_desc_type), intent(in), target   :: desc_a
      type(mld_sbaseprc_type), intent(inout), target :: p
      integer, intent(out)                      :: info
    end subroutine mld_smlprec_bld
    subroutine mld_dmlprec_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_dspmat_type), intent(inout), target :: a
      type(psb_desc_type), intent(in), target   :: desc_a
      type(mld_dbaseprc_type), intent(inout), target :: p
      integer, intent(out)                      :: info
    end subroutine mld_dmlprec_bld
    subroutine mld_cmlprec_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_cspmat_type), intent(inout), target :: a
      type(psb_desc_type), intent(in), target   :: desc_a
      type(mld_cbaseprc_type), intent(inout),target :: p
      integer, intent(out)                      :: info
    end subroutine mld_cmlprec_bld
    subroutine mld_zmlprec_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_zspmat_type), intent(inout), target :: a
      type(psb_desc_type), intent(in), target   :: desc_a
      type(mld_zbaseprc_type), intent(inout),target :: p
      integer, intent(out)                      :: info
    end subroutine mld_zmlprec_bld
  end interface

  interface mld_diag_bld
    subroutine mld_sdiag_bld(a,desc_data,p,info)
      use psb_base_mod
      use mld_prec_type
      integer, intent(out) :: info
      type(psb_sspmat_type), intent(in), target :: a
      type(psb_desc_type),intent(in)            :: desc_data
      type(mld_sbaseprc_type), intent(inout)    :: p
    end subroutine mld_sdiag_bld
    subroutine mld_ddiag_bld(a,desc_data,p,info)
      use psb_base_mod
      use mld_prec_type
      integer, intent(out) :: info
      type(psb_dspmat_type), intent(in), target :: a
      type(psb_desc_type),intent(in)            :: desc_data
      type(mld_dbaseprc_type), intent(inout)    :: p
    end subroutine mld_ddiag_bld
    subroutine mld_cdiag_bld(a,desc_data,p,info)
      use psb_base_mod
      use mld_prec_type
      integer, intent(out) :: info
      type(psb_cspmat_type), intent(in), target :: a
      type(psb_desc_type),intent(in)            :: desc_data
      type(mld_cbaseprc_type), intent(inout)    :: p
    end subroutine mld_cdiag_bld
    subroutine mld_zdiag_bld(a,desc_data,p,info)
      use psb_base_mod
      use mld_prec_type
      integer, intent(out) :: info
      type(psb_zspmat_type), intent(in), target :: a
      type(psb_desc_type),intent(in)            :: desc_data
      type(mld_zbaseprc_type), intent(inout)    :: p
    end subroutine mld_zdiag_bld
  end interface

  interface mld_fact_bld
    subroutine mld_sfact_bld(a,p,upd,info,blck)
      use psb_base_mod
      use mld_prec_type
      type(psb_sspmat_type), intent(in), target :: a
      type(mld_sbaseprc_type), intent(inout)    :: p
      integer, intent(out)                      :: info
      character, intent(in)                     :: upd
      type(psb_sspmat_type), intent(in), target, optional  :: blck
    end subroutine mld_sfact_bld
    subroutine mld_dfact_bld(a,p,upd,info,blck)
      use psb_base_mod
      use mld_prec_type
      type(psb_dspmat_type), intent(in), target :: a
      type(mld_dbaseprc_type), intent(inout)    :: p
      integer, intent(out)                      :: info
      character, intent(in)                     :: upd
      type(psb_dspmat_type), intent(in), target, optional  :: blck
    end subroutine mld_dfact_bld
    subroutine mld_cfact_bld(a,p,upd,info,blck)
      use psb_base_mod
      use mld_prec_type
      type(psb_cspmat_type), intent(in), target :: a
      type(mld_cbaseprc_type), intent(inout)    :: p
      integer, intent(out)                      :: info
      character, intent(in)                     :: upd
      type(psb_cspmat_type), intent(in), target, optional  :: blck
    end subroutine mld_cfact_bld
    subroutine mld_zfact_bld(a,p,upd,info,blck)
      use psb_base_mod
      use mld_prec_type
      type(psb_zspmat_type), intent(in), target :: a
      type(mld_zbaseprc_type), intent(inout)    :: p
      integer, intent(out)                      :: info
      character, intent(in)                     :: upd
      type(psb_zspmat_type), intent(in), target, optional  :: blck
    end subroutine mld_zfact_bld
  end interface

  interface mld_ilu_bld
    subroutine mld_silu_bld(a,p,upd,info,blck)
      use psb_base_mod
      use mld_prec_type
      integer, intent(out) :: info
      type(psb_sspmat_type), intent(in), target :: a
      type(mld_sbaseprc_type), intent(inout)    :: p
      character, intent(in)                     :: upd
      type(psb_sspmat_type), intent(in), optional :: blck
    end subroutine mld_silu_bld
    subroutine mld_dilu_bld(a,p,upd,info,blck)
      use psb_base_mod
      use mld_prec_type
      integer, intent(out) :: info
      type(psb_dspmat_type), intent(in), target :: a
      type(mld_dbaseprc_type), intent(inout)    :: p
      character, intent(in)                     :: upd
      type(psb_dspmat_type), intent(in), optional :: blck
    end subroutine mld_dilu_bld
    subroutine mld_cilu_bld(a,p,upd,info,blck)
      use psb_base_mod
      use mld_prec_type
      integer, intent(out) :: info
      type(psb_cspmat_type), intent(in), target :: a
      type(mld_cbaseprc_type), intent(inout)    :: p
      character, intent(in)                     :: upd
      type(psb_cspmat_type), intent(in), optional :: blck
    end subroutine mld_cilu_bld
    subroutine mld_zilu_bld(a,p,upd,info,blck)
      use psb_base_mod
      use mld_prec_type
      integer, intent(out) :: info
      type(psb_zspmat_type), intent(in), target :: a
      type(mld_zbaseprc_type), intent(inout)    :: p
      character, intent(in)                     :: upd
      type(psb_zspmat_type), intent(in), optional :: blck
    end subroutine mld_zilu_bld
  end interface

  interface mld_sludist_bld
    subroutine mld_ssludist_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_sspmat_type), intent(inout)   :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_sbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_ssludist_bld
    subroutine mld_dsludist_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_dspmat_type), intent(inout)   :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_dbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_dsludist_bld
    subroutine mld_csludist_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_cspmat_type), intent(inout)   :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_cbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_csludist_bld
    subroutine mld_zsludist_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_zspmat_type), intent(inout)   :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_zbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_zsludist_bld
  end interface

  interface mld_slu_bld
    subroutine mld_sslu_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_sspmat_type), intent(inout)      :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_sbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_sslu_bld
    subroutine mld_dslu_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_dspmat_type), intent(inout)      :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_dbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_dslu_bld
    subroutine mld_cslu_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_cspmat_type), intent(inout)   :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_cbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_cslu_bld
    subroutine mld_zslu_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_zspmat_type), intent(inout)   :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_zbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_zslu_bld
  end interface

  interface mld_umf_bld
    subroutine mld_sumf_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_sspmat_type), intent(inout)      :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_sbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_sumf_bld
    subroutine mld_dumf_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_dspmat_type), intent(inout)      :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_dbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_dumf_bld
    subroutine mld_zumf_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_zspmat_type), intent(in)      :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_zbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_zumf_bld
    subroutine mld_cumf_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_cspmat_type), intent(in)      :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_cbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_cumf_bld
  end interface

  interface mld_ilu0_fact
    subroutine mld_silu0_fact(ialg,a,l,u,d,info,blck)
      use psb_base_mod
      integer, intent(in)                 :: ialg
      integer, intent(out)                :: info
      type(psb_sspmat_type),intent(in)    :: a
      type(psb_sspmat_type),intent(inout) :: l,u
      type(psb_sspmat_type),intent(in), optional, target :: blck
      real(psb_spk_), intent(inout)     ::  d(:)
    end subroutine mld_silu0_fact
    subroutine mld_dilu0_fact(ialg,a,l,u,d,info,blck)
      use psb_base_mod
      integer, intent(in)                 :: ialg
      integer, intent(out)                :: info
      type(psb_dspmat_type),intent(in)    :: a
      type(psb_dspmat_type),intent(inout) :: l,u
      type(psb_dspmat_type),intent(in), optional, target :: blck
      real(psb_dpk_), intent(inout)     ::  d(:)
    end subroutine mld_dilu0_fact
    subroutine mld_cilu0_fact(ialg,a,l,u,d,info,blck)
      use psb_base_mod
      integer, intent(in)                 :: ialg
      integer, intent(out)                :: info
      type(psb_cspmat_type),intent(in)    :: a
      type(psb_cspmat_type),intent(inout) :: l,u
      type(psb_cspmat_type),intent(in), optional, target :: blck
      complex(psb_spk_), intent(inout)     ::  d(:)
    end subroutine mld_cilu0_fact
    subroutine mld_zilu0_fact(ialg,a,l,u,d,info,blck)
      use psb_base_mod
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
      use psb_base_mod
      integer, intent(in)                 :: fill_in,ialg
      integer, intent(out)                :: info
      type(psb_sspmat_type),intent(in)    :: a
      type(psb_sspmat_type),intent(inout) :: l,u
      type(psb_sspmat_type),intent(in), optional, target :: blck
      real(psb_spk_), intent(inout)     ::  d(:)
    end subroutine mld_siluk_fact
    subroutine mld_diluk_fact(fill_in,ialg,a,l,u,d,info,blck)
      use psb_base_mod
      integer, intent(in)                 :: fill_in,ialg
      integer, intent(out)                :: info
      type(psb_dspmat_type),intent(in)    :: a
      type(psb_dspmat_type),intent(inout) :: l,u
      type(psb_dspmat_type),intent(in), optional, target :: blck
      real(psb_dpk_), intent(inout)     ::  d(:)
    end subroutine mld_diluk_fact
    subroutine mld_ciluk_fact(fill_in,ialg,a,l,u,d,info,blck)
      use psb_base_mod
      integer, intent(in)                 :: fill_in,ialg
      integer, intent(out)                :: info
      type(psb_cspmat_type),intent(in)    :: a
      type(psb_cspmat_type),intent(inout) :: l,u
      type(psb_cspmat_type),intent(in), optional, target :: blck
      complex(psb_spk_), intent(inout)     ::  d(:)
    end subroutine mld_ciluk_fact
    subroutine mld_ziluk_fact(fill_in,ialg,a,l,u,d,info,blck)
      use psb_base_mod
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
      use psb_base_mod
      integer, intent(in)                 :: fill_in
      real(psb_spk_), intent(in)        :: thres
      integer, intent(out)                :: info
      type(psb_sspmat_type),intent(in)    :: a
      type(psb_sspmat_type),intent(inout) :: l,u
      type(psb_sspmat_type),intent(in), optional, target :: blck
      real(psb_spk_), intent(inout)     ::  d(:)
    end subroutine mld_silut_fact
    subroutine mld_dilut_fact(fill_in,thres,a,l,u,d,info,blck)
      use psb_base_mod
      integer, intent(in)                 :: fill_in
      real(psb_dpk_), intent(in)        :: thres
      integer, intent(out)                :: info
      type(psb_dspmat_type),intent(in)    :: a
      type(psb_dspmat_type),intent(inout) :: l,u
      type(psb_dspmat_type),intent(in), optional, target :: blck
      real(psb_dpk_), intent(inout)     ::  d(:)
    end subroutine mld_dilut_fact
    subroutine mld_cilut_fact(fill_in,thres,a,l,u,d,info,blck)
      use psb_base_mod
      integer, intent(in)                 :: fill_in
      real(psb_spk_), intent(in)        :: thres
      integer, intent(out)                :: info
      type(psb_cspmat_type),intent(in)    :: a
      type(psb_cspmat_type),intent(inout) :: l,u
      type(psb_cspmat_type),intent(in), optional, target :: blck
      complex(psb_spk_), intent(inout)  ::  d(:)
    end subroutine mld_cilut_fact
    subroutine mld_zilut_fact(fill_in,thres,a,l,u,d,info,blck)
      use psb_base_mod
      integer, intent(in)                 :: fill_in
      real(psb_dpk_), intent(in)        :: thres
      integer, intent(out)                :: info
      type(psb_zspmat_type),intent(in)    :: a
      type(psb_zspmat_type),intent(inout) :: l,u
      type(psb_zspmat_type),intent(in), optional, target :: blck
      complex(psb_dpk_), intent(inout)  ::  d(:)
    end subroutine mld_zilut_fact
  end interface
end module mld_basep_bld_mod

