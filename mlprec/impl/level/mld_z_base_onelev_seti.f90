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
subroutine mld_z_base_onelev_seti(lv,what,val,info)
  
  use psb_base_mod
  use mld_z_onelev_mod, mld_protect_name => mld_z_base_onelev_seti

  Implicit None

  ! Arguments
  class(mld_z_onelev_type), intent(inout) :: lv 
  integer, intent(in)                          :: what 
  integer, intent(in)                          :: val
  integer, intent(out)                         :: info
  Integer           :: err_act
  character(len=20) :: name='z_base_onelev_seti'

  call psb_erractionsave(err_act)
  info = psb_success_

  select case (what) 

  case (mld_smoother_sweeps_)
    lv%parms%sweeps      = val
    lv%parms%sweeps_pre  = val
    lv%parms%sweeps_post = val

  case (mld_smoother_sweeps_pre_)
    lv%parms%sweeps_pre  = val

  case (mld_smoother_sweeps_post_)
    lv%parms%sweeps_post = val

  case (mld_ml_type_)
    lv%parms%ml_type       = val

  case (mld_aggr_alg_)
    lv%parms%aggr_alg      = val

  case (mld_aggr_kind_)
    lv%parms%aggr_kind     = val

  case (mld_coarse_mat_)
    lv%parms%coarse_mat    = val

  case (mld_smoother_pos_)
    lv%parms%smoother_pos  = val

  case (mld_aggr_omega_alg_)
    lv%parms%aggr_omega_alg= val

  case (mld_aggr_eig_)
    lv%parms%aggr_eig      = val

  case (mld_aggr_filter_)
    lv%parms%aggr_filter   = val

  case (mld_coarse_solve_)
    lv%parms%coarse_solve    = val

  case default
    if (allocated(lv%sm)) then 
      call lv%sm%set(what,val,info)
    end if
    if (info /= psb_success_) goto 9999
  end select
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine mld_z_base_onelev_seti
