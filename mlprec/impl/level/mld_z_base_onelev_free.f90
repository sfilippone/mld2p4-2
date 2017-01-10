!!$
!!$ 
!!$                           MLD2P4  version 2.1
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.4)
!!$  
!!$  (C) Copyright 2008, 2010, 2012, 2015, 2017 
!!$
!!$                      Salvatore Filippone    Cranfield University
!!$		         Ambra Abdullahi Hassan University of Rome Tor Vergata
!!$                      Alfredo Buttari        CNRS-IRIT, Toulouse
!!$                      Pasqua D'Ambra         ICAR-CNR, Naples
!!$                      Daniela di Serafino    Second University of Naples
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
subroutine mld_z_base_onelev_free(lv,info)
  
  use psb_base_mod
  use mld_z_onelev_mod, mld_protect_name =>  mld_z_base_onelev_free
  implicit none 

  class(mld_z_onelev_type), intent(inout) :: lv
  integer(psb_ipk_), intent(out)          :: info
  integer(psb_ipk_) :: i

  info = psb_success_

  ! We might just deallocate the top level array, except 
  ! that there may be inner objects containing C pointers,
  ! e.g.  UMFPACK, SLU or CUDA stuff.
  ! We really need FINALs. 
  if (allocated(lv%sm)) &
       & call lv%sm%free(info)

  call lv%ac%free()
  if (lv%desc_ac%is_ok()) &
       & call lv%desc_ac%free(info)
  call lv%map%free(info)

  ! This is a pointer to something else, must not free it here. 
  nullify(lv%base_a) 
  ! This is a pointer to something else, must not free it here. 
  nullify(lv%base_desc) 

  call lv%nullify()

end subroutine mld_z_base_onelev_free
