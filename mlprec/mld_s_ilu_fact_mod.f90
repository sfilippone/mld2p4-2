!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010,2012
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
module mld_s_ilu_fact_mod

  use mld_s_base_solver_mod 

  interface mld_ilu0_fact
    subroutine mld_silu0_fact(ialg,a,l,u,d,info,blck,upd)
      import psb_sspmat_type, psb_spk_, psb_ipk_
      integer(psb_ipk_), intent(in)         :: ialg
      integer(psb_ipk_), intent(out)        :: info
      type(psb_sspmat_type),intent(in)      :: a
      type(psb_sspmat_type),intent(inout)   :: l,u
      type(psb_sspmat_type),intent(in), optional, target :: blck
      character, intent(in), optional       :: upd
      real(psb_spk_), intent(inout)      :: d(:)
    end subroutine mld_silu0_fact
  end interface

  interface mld_iluk_fact
    subroutine mld_siluk_fact(fill_in,ialg,a,l,u,d,info,blck)
      import psb_sspmat_type, psb_spk_, psb_ipk_
      integer(psb_ipk_), intent(in)        :: fill_in,ialg
      integer(psb_ipk_), intent(out)       :: info
      type(psb_sspmat_type),intent(in)     :: a
      type(psb_sspmat_type),intent(inout)  :: l,u
      type(psb_sspmat_type),intent(in), optional, target :: blck
      real(psb_spk_), intent(inout)     ::  d(:)
    end subroutine mld_siluk_fact
  end interface

  interface mld_ilut_fact
    subroutine mld_silut_fact(fill_in,thres,a,l,u,d,info,blck,iscale)
      import  psb_sspmat_type, psb_spk_, psb_ipk_
      integer(psb_ipk_), intent(in)        :: fill_in
      real(psb_spk_), intent(in)           :: thres
      integer(psb_ipk_), intent(out)       :: info
      type(psb_sspmat_type),intent(in)     :: a
      type(psb_sspmat_type),intent(inout)  :: l,u
      real(psb_spk_), intent(inout)     :: d(:)
      type(psb_sspmat_type),intent(in), optional, target :: blck
      integer(psb_ipk_), intent(in), optional  :: iscale
    end subroutine mld_silut_fact
  end interface

end module mld_s_ilu_fact_mod
