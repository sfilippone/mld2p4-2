!  
!   
!                             MLD2P4  version 2.1
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!    
!    (C) Copyright 2008, 2010, 2012, 2015, 2017 , 2017 
!  
!                        Salvatore Filippone  Cranfield University
!  		      Ambra Abdullahi Hassan University of Rome Tor Vergata
!        Pasqua D'Ambra         IAC-CNR, Naples, IT
!        Daniela di Serafino    University of Campania "L. Vanvitelli", Caserta, IT
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
subroutine mld_d_base_onelev_setag(lev,val,info,pos)

  use psb_base_mod
  use mld_d_onelev_mod, mld_protect_name => mld_d_base_onelev_setag

  implicit none

  ! Arguments
  class(mld_d_onelev_type), target, intent(inout) :: lev
  class(mld_d_base_aggregator_type), intent(in)   :: val
  integer(psb_ipk_), intent(out)                  :: info
  character(len=*), optional, intent(in)          :: pos
  
  ! Local variables
  integer(psb_ipk_)                :: ipos_
  character(len=*), parameter      :: name='mld_base_onelev_setag'

  info = psb_success_

  ! Ignore pos for aggregator
  
  if (allocated(lev%aggr)) then 
    if (.not.same_type_as(lev%aggr,val))  then
      call lev%aggr%free(info)
      deallocate(lev%aggr,stat=info)
      if (info /= 0) then
        info = 3111
        return
      end if
    end if
  end if
      
  if (.not.allocated(lev%aggr)) then 
    allocate(lev%aggr,mold=val,stat=info) 
    if (info /= 0) then
      info = 3111
      return
    end if
  end if
  
end subroutine mld_d_base_onelev_setag

