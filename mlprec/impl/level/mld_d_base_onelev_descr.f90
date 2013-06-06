!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010,2012,2013
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
subroutine mld_d_base_onelev_descr(lv,il,nl,info,iout)
  
  use psb_base_mod
  use mld_d_onelev_mod, mld_protect_name => mld_d_base_onelev_descr
  Implicit None
  ! Arguments
  class(mld_d_onelev_type), intent(in)  :: lv
  integer(psb_ipk_), intent(in)           :: il,nl
  integer(psb_ipk_), intent(out)          :: info
  integer(psb_ipk_), intent(in), optional :: iout

  ! Local variables
  integer(psb_ipk_)      :: err_act
  character(len=20), parameter :: name='mld_d_base_onelev_descr'
  integer(psb_ipk_)      :: iout_
  logical      :: coarse


  call psb_erractionsave(err_act)


  coarse = (il==nl)

  if (present(iout)) then 
    iout_ = iout
  else 
    iout_ = 6
  end if

  write(iout_,*) 
  if (il == 2) then 
    call lv%parms%mldescr(iout_,info)
    write(iout_,*) 
  end if

  if (coarse)  then 
    write(iout_,*) ' Level ',il,' (coarsest)'
  else
    write(iout_,*) ' Level ',il
  end if

  call lv%parms%descr(iout_,info,coarse=coarse)

  if (nl > 1) then 
    if (allocated(lv%map%naggr)) then
      write(iout_,*) '  Size of coarse matrix: ', &
           &  sum(lv%map%naggr(:))
      write(iout_,*) '  Sizes of aggregates: ', &
           &  lv%map%naggr(:)
    end if
  end if

  if (coarse.and.allocated(lv%sm)) &
       & call lv%sm%descr(info,iout=iout_,coarse=coarse)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine mld_d_base_onelev_descr
