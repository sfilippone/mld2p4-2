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
module mld_basep_bld
  interface mld_baseprc_bld
    subroutine mld_dbaseprc_bld(a,desc_a,p,info,upd)
      use psb_base_mod
      use mld_prec_type
      type(psb_dspmat_type), target              :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(mld_dbaseprc_type),intent(inout)      :: p
      integer, intent(out)                       :: info
      character, intent(in), optional            :: upd
    end subroutine mld_dbaseprc_bld
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
    subroutine mld_das_bld(a,desc_a,p,upd,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_dspmat_type), target           :: a
      type(psb_desc_type), intent(in), target :: desc_a
      type(mld_dbaseprc_type),intent(inout)   :: p
      character, intent(in)                   :: upd
      integer, intent(out)                    :: info
    end subroutine mld_das_bld
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
    subroutine mld_dmlprec_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_dspmat_type), intent(inout), target :: a
      type(psb_desc_type), intent(in), target   :: desc_a
      type(mld_dbaseprc_type), intent(inout), target :: p
      integer, intent(out)                      :: info
    end subroutine mld_dmlprec_bld
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
    subroutine mld_ddiag_bld(a,desc_data,p,info)
      use psb_base_mod
      use mld_prec_type
      integer, intent(out) :: info
      type(psb_dspmat_type), intent(in), target :: a
      type(psb_desc_type),intent(in)            :: desc_data
      type(mld_dbaseprc_type), intent(inout)    :: p
    end subroutine mld_ddiag_bld
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
    subroutine mld_dfact_bld(a,p,upd,info,blck)
      use psb_base_mod
      use mld_prec_type
      type(psb_dspmat_type), intent(in), target :: a
      type(mld_dbaseprc_type), intent(inout)    :: p
      integer, intent(out)                      :: info
      character, intent(in)                     :: upd
      type(psb_dspmat_type), intent(in), target, optional  :: blck
    end subroutine mld_dfact_bld
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
    subroutine mld_dilu_bld(a,p,upd,info,blck)
      use psb_base_mod
      use mld_prec_type
      integer, intent(out) :: info
      type(psb_dspmat_type), intent(in), target :: a
      type(mld_dbaseprc_type), intent(inout)    :: p
      character, intent(in)                     :: upd
      type(psb_dspmat_type), intent(in), optional :: blck
    end subroutine mld_dilu_bld
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
    subroutine mld_dsludist_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_dspmat_type), intent(inout)   :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_dbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_dsludist_bld
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
    subroutine mld_dslu_bld(a,desc_a,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_dspmat_type), intent(inout)      :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(mld_dbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine mld_dslu_bld
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
  end interface

  interface mld_ilu0_fact
    subroutine mld_dilu0_fact(ialg,a,l,u,d,info,blck)
      use psb_base_mod
      integer, intent(in)                 :: ialg
      integer, intent(out)                :: info
      type(psb_dspmat_type),intent(in)    :: a
      type(psb_dspmat_type),intent(inout) :: l,u
      type(psb_dspmat_type),intent(in), optional, target :: blck
      real(kind(1.d0)), intent(inout)     ::  d(:)
    end subroutine mld_dilu0_fact
    subroutine mld_zilu0_fact(ialg,a,l,u,d,info,blck)
      use psb_base_mod
      integer, intent(in)                 :: ialg
      integer, intent(out)                :: info
      type(psb_zspmat_type),intent(in)    :: a
      type(psb_zspmat_type),intent(inout) :: l,u
      type(psb_zspmat_type),intent(in), optional, target :: blck
      complex(kind(1.d0)), intent(inout)     ::  d(:)
    end subroutine mld_zilu0_fact
  end interface

  interface mld_iluk_fact
    subroutine mld_diluk_fact(fill_in,ialg,a,l,u,d,info,blck)
      use psb_base_mod
      integer, intent(in)                 :: fill_in,ialg
      integer, intent(out)                :: info
      type(psb_dspmat_type),intent(in)    :: a
      type(psb_dspmat_type),intent(inout) :: l,u
      type(psb_dspmat_type),intent(in), optional, target :: blck
      real(kind(1.d0)), intent(inout)     ::  d(:)
    end subroutine mld_diluk_fact
    subroutine mld_ziluk_fact(fill_in,ialg,a,l,u,d,info,blck)
      use psb_base_mod
      integer, intent(in)                 :: fill_in,ialg
      integer, intent(out)                :: info
      type(psb_zspmat_type),intent(in)    :: a
      type(psb_zspmat_type),intent(inout) :: l,u
      type(psb_zspmat_type),intent(in), optional, target :: blck
      complex(kind(1.d0)), intent(inout)     ::  d(:)
    end subroutine mld_ziluk_fact
  end interface

  interface mld_ilut_fact
    subroutine mld_dilut_fact(fill_in,thres,a,l,u,d,info,blck)
      use psb_base_mod
      integer, intent(in)                 :: fill_in
      real(kind(1.d0)), intent(in)        :: thres
      integer, intent(out)                :: info
      type(psb_dspmat_type),intent(in)    :: a
      type(psb_dspmat_type),intent(inout) :: l,u
      type(psb_dspmat_type),intent(in), optional, target :: blck
      real(kind(1.d0)), intent(inout)     ::  d(:)
    end subroutine mld_dilut_fact
    subroutine mld_zilut_fact(fill_in,thres,a,l,u,d,info,blck)
      use psb_base_mod
      integer, intent(in)                 :: fill_in
      real(kind(1.d0)), intent(in)        :: thres
      integer, intent(out)                :: info
      type(psb_zspmat_type),intent(in)    :: a
      type(psb_zspmat_type),intent(inout) :: l,u
      type(psb_zspmat_type),intent(in), optional, target :: blck
      complex(kind(1.d0)), intent(inout)  ::  d(:)
    end subroutine mld_zilut_fact
  end interface
end module mld_basep_bld

module mld_inner_mod
  use mld_prec_type
  use mld_basep_bld
  interface mld_baseprec_aply
    subroutine mld_dbaseprec_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_dbaseprc_type), intent(in) :: prec
      real(kind(0.d0)),intent(in)         :: x(:)
      real(kind(0.d0)),intent(inout)      :: y(:)
      real(kind(0.d0)),intent(in)         :: alpha,beta
      character(len=1)                    :: trans
      real(kind(0.d0)),target             :: work(:)
      integer, intent(out)                :: info
    end subroutine mld_dbaseprec_aply
    subroutine mld_zbaseprec_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_zbaseprc_type), intent(in) :: prec
      complex(kind(1.d0)),intent(in)      :: x(:)
      complex(kind(1.d0)),intent(inout)   :: y(:)
      complex(kind(1.d0)),intent(in)      :: alpha,beta
      character(len=1)                    :: trans
      complex(kind(1.d0)),target          :: work(:)
      integer, intent(out)                :: info
    end subroutine mld_zbaseprec_aply
  end interface

  interface mld_as_aply
    subroutine mld_das_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_dbaseprc_type), intent(in) :: prec
      real(kind(0.d0)),intent(in)         :: x(:)
      real(kind(0.d0)),intent(inout)      :: y(:)
      real(kind(0.d0)),intent(in)         :: alpha,beta
      character(len=1)                    :: trans
      real(kind(0.d0)),target             :: work(:)
      integer, intent(out)                :: info
    end subroutine mld_das_aply
    subroutine mld_zas_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_zbaseprc_type), intent(in) :: prec
      complex(kind(1.d0)),intent(in)      :: x(:)
      complex(kind(1.d0)),intent(inout)   :: y(:)
      complex(kind(1.d0)),intent(in)      :: alpha,beta
      character(len=1)                    :: trans
      complex(kind(1.d0)),target          :: work(:)
      integer, intent(out)                :: info
    end subroutine mld_zas_aply
  end interface

  interface mld_mlprec_aply
    subroutine mld_dmlprec_aply(alpha,baseprecv,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_dbaseprc_type), intent(in) :: baseprecv(:)
      real(kind(0.d0)),intent(in)         :: alpha,beta
      real(kind(0.d0)),intent(in)         :: x(:)
      real(kind(0.d0)),intent(inout)      :: y(:)
      character                           :: trans
      real(kind(0.d0)),target             :: work(:)
      integer, intent(out)                :: info
    end subroutine mld_dmlprec_aply
    subroutine mld_zmlprec_aply(alpha,baseprecv,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_zbaseprc_type), intent(in) :: baseprecv(:)
      complex(kind(0.d0)),intent(in)      :: alpha,beta
      complex(kind(0.d0)),intent(in)      :: x(:)
      complex(kind(0.d0)),intent(inout)   :: y(:)
      character                           :: trans
      complex(kind(0.d0)),target          :: work(:)
      integer, intent(out)                :: info
    end subroutine mld_zmlprec_aply
  end interface


  interface mld_sub_aply
    subroutine mld_dsub_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_desc_type), intent(in)       :: desc_data
      type(mld_dbaseprc_type), intent(in)   :: prec
      real(kind(0.d0)),intent(in)           :: x(:)
      real(kind(0.d0)),intent(inout)        :: y(:)
      real(kind(0.d0)),intent(in)           :: alpha,beta
      character(len=1),intent(in)           :: trans
      real(kind(0.d0)),target,intent(inout) :: work(:)
      integer, intent(out)                  :: info
    end subroutine mld_dsub_aply
    subroutine mld_zsub_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_desc_type), intent(in)          :: desc_data
      type(mld_zbaseprc_type), intent(in)      :: prec
      complex(kind(0.d0)),intent(in)           :: x(:)
      complex(kind(0.d0)),intent(inout)        :: y(:)
      complex(kind(0.d0)),intent(in)           :: alpha,beta
      character(len=1),intent(in)              :: trans
      complex(kind(0.d0)),target,intent(inout) :: work(:)
      integer, intent(out)                     :: info
    end subroutine mld_zsub_aply
  end interface


  interface mld_sub_solve
    subroutine mld_dsub_solve(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_desc_type), intent(in)       :: desc_data
      type(mld_dbaseprc_type), intent(in)   :: prec
      real(kind(0.d0)),intent(in)           :: x(:)
      real(kind(0.d0)),intent(inout)        :: y(:)
      real(kind(0.d0)),intent(in)           :: alpha,beta
      character(len=1),intent(in)           :: trans
      real(kind(0.d0)),target,intent(inout) :: work(:)
      integer, intent(out)                  :: info
    end subroutine mld_dsub_solve
    subroutine mld_zsub_solve(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_desc_type), intent(in)          :: desc_data
      type(mld_zbaseprc_type), intent(in)      :: prec
      complex(kind(0.d0)),intent(in)           :: x(:)
      complex(kind(0.d0)),intent(inout)        :: y(:)
      complex(kind(0.d0)),intent(in)           :: alpha,beta
      character(len=1),intent(in)              :: trans
      complex(kind(0.d0)),target,intent(inout) :: work(:)
      integer, intent(out)                     :: info
    end subroutine mld_zsub_solve
  end interface


  interface mld_asmat_bld
    Subroutine mld_dasmat_bld(ptype,novr,a,blk,desc_data,upd,desc_p,info,outfmt)
      use psb_base_mod
      use mld_prec_type
      integer, intent(in)                 :: ptype,novr
      Type(psb_dspmat_type), Intent(in)   ::  a
      Type(psb_dspmat_type), Intent(out)  ::  blk
      Type(psb_desc_type), Intent(inout)  :: desc_p
      Type(psb_desc_type), Intent(in)     :: desc_data 
      Character, Intent(in)               :: upd
      integer, intent(out)                :: info
      character(len=5), optional          :: outfmt
    end Subroutine mld_dasmat_bld
    Subroutine mld_zasmat_bld(ptype,novr,a,blk,desc_data,upd,desc_p,info,outfmt)
      use psb_base_mod
      use mld_prec_type
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
    subroutine mld_dsp_renum(a,blck,p,atmp,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_dspmat_type), intent(in)      :: a,blck
      type(psb_dspmat_type), intent(out)     :: atmp
      type(mld_dbaseprc_type), intent(inout) :: p
      integer, intent(out)   :: info
    end subroutine mld_dsp_renum
    subroutine mld_zsp_renum(a,blck,p,atmp,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_zspmat_type), intent(in)      :: a,blck
      type(psb_zspmat_type), intent(out)     :: atmp
      type(mld_zbaseprc_type), intent(inout) :: p
      integer, intent(out)   :: info
    end subroutine mld_zsp_renum
  end interface

  interface mld_aggrmap_bld
    subroutine mld_daggrmap_bld(aggr_type,a,desc_a,nlaggr,ilaggr,info)
      use psb_base_mod
      use mld_prec_type
      integer, intent(in)               :: aggr_type
      type(psb_dspmat_type), intent(in), target :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer, allocatable, intent(out) :: ilaggr(:),nlaggr(:)
      integer, intent(out)              :: info
    end subroutine mld_daggrmap_bld
    subroutine mld_zaggrmap_bld(aggr_type,a,desc_a,nlaggr,ilaggr,info)
      use psb_base_mod
      use mld_prec_type
      integer, intent(in)               :: aggr_type
      type(psb_zspmat_type), intent(in), target :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer, allocatable, intent(out) :: ilaggr(:),nlaggr(:)
      integer, intent(out)              :: info
    end subroutine mld_zaggrmap_bld
  end interface

  interface mld_aggrmat_asb
    subroutine mld_daggrmat_asb(a,desc_a,ac,desc_ac,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_dspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      type(psb_dspmat_type), intent(out)             :: ac
      type(psb_desc_type), intent(out)               :: desc_ac
      type(mld_dbaseprc_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_daggrmat_asb
    subroutine mld_zaggrmat_asb(a,desc_a,ac,desc_ac,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_zspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      type(psb_zspmat_type), intent(out)             :: ac
      type(psb_desc_type), intent(out)               :: desc_ac
      type(mld_zbaseprc_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_zaggrmat_asb
  end interface

  interface mld_aggrmat_raw_asb
    subroutine mld_daggrmat_raw_asb(a,desc_a,ac,desc_ac,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_dspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      type(psb_dspmat_type), intent(out)             :: ac
      type(psb_desc_type), intent(out)               :: desc_ac
      type(mld_dbaseprc_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_daggrmat_raw_asb
    subroutine mld_zaggrmat_raw_asb(a,desc_a,ac,desc_ac,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_zspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      type(psb_zspmat_type), intent(out)             :: ac
      type(psb_desc_type), intent(out)               :: desc_ac
      type(mld_zbaseprc_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_zaggrmat_raw_asb
  end interface

  interface mld_aggrmat_smth_asb
    subroutine mld_daggrmat_smth_asb(a,desc_a,ac,desc_ac,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_dspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      type(psb_dspmat_type), intent(out)             :: ac
      type(psb_desc_type), intent(out)               :: desc_ac
      type(mld_dbaseprc_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_daggrmat_smth_asb
    subroutine mld_zaggrmat_smth_asb(a,desc_a,ac,desc_ac,p,info)
      use psb_base_mod
      use mld_prec_type
      type(psb_zspmat_type), intent(in)              :: a
      type(psb_desc_type), intent(in)                :: desc_a
      type(psb_zspmat_type), intent(out)             :: ac
      type(psb_desc_type), intent(out)               :: desc_ac
      type(mld_zbaseprc_type), intent(inout), target :: p
      integer, intent(out)                           :: info
    end subroutine mld_zaggrmat_smth_asb
  end interface

end module mld_inner_mod
