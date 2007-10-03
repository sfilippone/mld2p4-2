!!$ 
!!$ 
!!$                    MD2P4
!!$    Multilevel Domain Decomposition Parallel Preconditioner Package for PSBLAS
!!$                      for 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$                       Daniela di Serafino    Second University of Naples
!!$                       Pasqua D'Ambra         ICAR-CNR                      
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the MD2P4 group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MD2P4 GROUP OR ITS CONTRIBUTORS
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
! This version of the preconditioner module PSB_PREC_MODE renames 
! "on the fly" the MLD preconditioner routines and data types so that 
! the Krylov iterations in PSBLAS can be tricked into using the MLD versions
! instead of the original ones. Since there is no native runtime polymorphism 
! this implies a recompilation, but thanks to the renaming feature of 
! Fortran 95 (and to the compatibility of the calling sequences) there is 
! no need to change the source code for the Krylov methods.
!
module psb_prec_mod

#if (__GNUC__==4) && (__GNUC_MINOR__<=2)
  ! GNU Fortran 4.2.
  ! Workaround for PR 32634, it is fixed in GNU Fortran 4.3, will 
  ! it be fixed in 4.2??? 
  use mld_prec_type, &
       & psb_dbaseprc_type    => mld_dbaseprc_type,&
       & psb_zbaseprc_type    => mld_zbaseprc_type,&
       & psb_dprec_type       => mld_dprec_type,&
       & psb_zprec_type       => mld_zprec_type,&
       & psb_base_precfree    => mld_base_precfree,&
       & psb_nullify_baseprec => mld_nullify_baseprec,&
       & psb_prec_descr       => mld_prec_descr,&
       & psb_prec_short_descr => mld_prec_short_descr

  use mld_prec_mod


  interface psb_precbld
    module procedure mld_dprecbld, mld_zprecbld
  end interface

  interface psb_precinit
    module procedure  mld_dprecinit, mld_zprecinit
  end interface

  interface psb_precset
    module procedure  mld_dprecseti, mld_dprecsetd,&
         &  mld_zprecseti,  mld_zprecsetd
  end interface

  interface psb_precfree
    module procedure  mld_dprecfree,  mld_zprecfree
  end interface

  interface psb_precaply
    module procedure  mld_dprec_aply,  mld_dprec_aply1, &
         &  mld_zprec_aply,  mld_zprec_aply1
  end interface

#else 

  use mld_prec_mod, &
       & psb_dbaseprc_type    => mld_dbaseprc_type,&
       & psb_zbaseprc_type    => mld_zbaseprc_type,&
       & psb_dprec_type       => mld_dprec_type,&
       & psb_zprec_type       => mld_zprec_type,&
       & psb_base_precfree    => mld_base_precfree,&
       & psb_nullify_baseprec => mld_nullify_baseprec,&
       & psb_prec_descr       => mld_prec_descr,&
       & psb_prec_short_descr => mld_prec_short_descr,&
       & psb_precbld          => mld_precbld,   &
       & psb_precinit         => mld_precinit,  &
       & psb_precfree         => mld_precfree, &
       & psb_precset          => mld_precset,  &
       & psb_precaply         => mld_precaply

#endif

  interface psb_sizeof
    module procedure mld_dprec_sizeof, mld_zprec_sizeof, &
         & mld_dbaseprc_sizeof, mld_zbaseprc_sizeof
  end interface


end module psb_prec_mod
