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
subroutine mld_z_jac_smoother_bld(a,desc_a,sm,upd,info,amold,vmold)

  use psb_base_mod
  use mld_z_diag_solver
  use mld_z_jac_smoother, mld_protect_name => mld_z_jac_smoother_bld
  Implicit None

  ! Arguments
  type(psb_zspmat_type), intent(in), target          :: a
  Type(psb_desc_type), Intent(inout)                   :: desc_a 
  class(mld_z_jac_smoother_type), intent(inout)      :: sm
  character, intent(in)                                :: upd
  integer(psb_ipk_), intent(out)                       :: info
  class(psb_z_base_sparse_mat), intent(in), optional :: amold
  class(psb_z_base_vect_type), intent(in), optional  :: vmold
  ! Local variables
  integer(psb_ipk_) :: n_row,n_col, nrow_a, nztota, nzeros
  complex(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer(psb_ipk_) :: ictxt,np,me,i, err_act, debug_unit, debug_level
  character(len=20) :: name='z_jac_smoother_bld', ch_err

  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ictxt       = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'


  n_row  = desc_a%get_local_rows()

  nrow_a = a%get_nrows()
  nztota = a%get_nzeros()
  select type (smsv => sm%sv)
  type is (mld_z_diag_solver_type)
    call a%clip_diag(sm%nd,info)
    class default
    call a%csclip(sm%nd,info,&
         & jmin=nrow_a+1,rscale=.false.,cscale=.false.)
  end select
  if (info == psb_success_) then 
    if (present(amold)) then 
      call sm%nd%cscnv(info,&
           & mold=amold,dupl=psb_dupl_add_)
    else
      call sm%nd%cscnv(info,&
           & type='csr',dupl=psb_dupl_add_)
    endif
  end if
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,&
         & a_err='clip & psb_spcnv csr 4')
    goto 9999
  end if

  call sm%sv%build(a,desc_a,upd,info,amold=amold,vmold=vmold)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,&
         & a_err='solver build')
    goto 9999
  end if
  nzeros = sm%nd%get_nzeros()
  call psb_sum(ictxt,nzeros)
  sm%nnz_nd_tot = nzeros
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' end'

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_z_jac_smoother_bld
