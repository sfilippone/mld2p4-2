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
subroutine mld_zmlprec_bld(a,desc_a,p,info)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_zmlprec_bld

  implicit none 

  type(psb_zspmat_type), intent(in), target :: a
  type(psb_desc_type), intent(in), target   :: desc_a
  type(mld_zbaseprc_type), intent(inout),target    :: p
  integer, intent(out)                      :: info

  type(psb_desc_type)                       :: desc_ac

  integer :: i, nrg, nzg, err_act,k
  character(len=20) :: name, ch_err
  logical, parameter :: debug=.false.
  type(psb_zspmat_type)                     :: ac
  integer :: ictxt, np, me

  name='psb_zmlprec_bld'
  if (psb_get_errstatus().ne.0) return 
  info = 0
  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt,me,np)
  call psb_erractionsave(err_act)
  call psb_nullify_sp(ac)


  if (.not.allocated(p%iprcparm)) then 
    info = 2222
    call psb_errpush(info,name)
    goto 9999
  endif
  call mld_check_def(p%iprcparm(mld_ml_type_),'Multilevel type',&
       &   mld_mult_ml_,is_legal_ml_type)
  call mld_check_def(p%iprcparm(mld_aggr_alg_),'aggregation',&
       &   mld_dec_aggr_,is_legal_ml_aggr_kind)
  call mld_check_def(p%iprcparm(mld_aggr_kind_),'Smoother kind',&
       &   mld_smooth_prol_,is_legal_ml_smth_kind)
  call mld_check_def(p%iprcparm(mld_coarse_mat_),'Coarse matrix',&
       &   mld_distr_mat_,is_legal_ml_coarse_mat)
  call mld_check_def(p%iprcparm(mld_smooth_pos_),'smooth_pos',&
       &   mld_pre_smooth_,is_legal_ml_smooth_pos)


!!$  nullify(p%desc_data)
  select case(p%iprcparm(mld_sub_solve_))
  case(mld_ilu_n_)      
    call mld_check_def(p%iprcparm(mld_sub_fill_in_),'Level',0,is_legal_ml_lev)
  case(mld_ilu_t_)                 
    call mld_check_def(p%dprcparm(mld_fact_thrs_),'Eps',dzero,is_legal_fact_thrs)
  end select
  call mld_check_def(p%dprcparm(mld_aggr_damp_),'omega',dzero,is_legal_omega)
  call mld_check_def(p%iprcparm(mld_smooth_sweeps_),'Jacobi sweeps',&
       & 1,is_legal_jac_sweeps)


  ! Currently this is ignored by gen_aggrmap, but it could be 
  ! changed in the future. Need to package nlaggr & mlia in a 
  ! private data structure? 
  call mld_aggrmap_bld(p%iprcparm(mld_aggr_alg_),a,desc_a,p%nlaggr,p%mlia,info)
    
  if(info /= 0) then
    info=4010
    ch_err='mld_aggrmap_bld'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (debug) write(0,*) 'Out from genaggrmap',p%nlaggr

  call psb_nullify_desc(desc_ac)
  call mld_aggrmat_asb(a,desc_a,ac,desc_ac,p,info)
  if(info /= 0) then
    info=4010
    ch_err='mld_aggrmat_asb'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  if (debug) write(0,*) 'Out from bldaggrmat',desc_ac%matrix_data(:)



  call mld_baseprc_bld(ac,desc_ac,p,info)
  if (debug) write(0,*) 'Out from baseprcbld',info
  if(info /= 0) then
    info=4010
    ch_err='mld_baseprc_bld'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  

  !
  ! We have used a separate ac because:
  ! 1. We want to reuse the same routines mld_ilu_bld etc.
  ! 2. We do NOT want to pass an argument twice to them 
  !    p%av(mld_ac_) and p, as this would violate the Fortran standard
  ! Hence a separate AC and a TRANSFER function at the end. 
  !
  call psb_sp_transfer(ac,p%av(mld_ac_),info)
  p%base_a => p%av(mld_ac_)
  call psb_cdtransfer(desc_ac,p%desc_ac,info)

  if (info /= 0) then 
    info=4010
    ch_err='psb_cdtransfer'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  p%base_desc => p%desc_ac

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  Return

end subroutine mld_zmlprec_bld
