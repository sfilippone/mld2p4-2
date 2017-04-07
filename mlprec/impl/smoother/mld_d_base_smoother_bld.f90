!  
!   
!                             MLD2P4  version 2.1
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.4)
!    
!    (C) Copyright 2008, 2010, 2012, 2015, 2017 
!  
!        Salvatore Filippone    Cranfield University
!        Ambra Abdullahi Hassan University of Rome Tor Vergata
!        Alfredo Buttari        CNRS-IRIT, Toulouse
!        Pasqua D'Ambra         ICAR-CNR, Naples
!        Daniela di Serafino    University of Campania "L. Vanvitelli", Caserta
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
subroutine mld_d_base_smoother_bld(a,desc_a,sm,upd,info,amold,vmold,imold)
  
  use psb_base_mod
  use mld_d_base_smoother_mod, mld_protect_name =>  mld_d_base_smoother_bld
  Implicit None

  ! Arguments
  type(psb_dspmat_type), intent(in), target      :: a
  Type(psb_desc_type), Intent(inout)               :: desc_a 
  class(mld_d_base_smoother_type), intent(inout) :: sm 
  character, intent(in)                            :: upd
  integer(psb_ipk_), intent(out)                   :: info
  class(psb_d_base_sparse_mat), intent(in), optional :: amold
  class(psb_d_base_vect_type), intent(in), optional  :: vmold
  class(psb_i_base_vect_type), intent(in), optional  :: imold
  integer(psb_ipk_) :: err_act
  character(len=20) :: name='d_base_smoother_bld'

  call psb_erractionsave(err_act)

  info = psb_success_
  if (allocated(sm%sv)) then 
    call sm%sv%build(a,desc_a,upd,info,amold=amold,vmold=vmold)
  else
    info = 1121
    call psb_errpush(info,name)
  endif
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end subroutine mld_d_base_smoother_bld
