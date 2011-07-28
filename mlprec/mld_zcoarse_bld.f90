!!$ 
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010
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
! File: mld_zcoarse_bld.f90
!
! Subroutine: mld_zcoarse_bld
! Version:    complex
!
!  This routine builds the matrix associated to the current level of the
!  multilevel preconditioner from the matrix associated to the previous level,
!  by using a smoothed aggregation technique (therefore, it also builds the
!  prolongation and restriction operators mapping the current level to the
!  previous one and vice versa). Then the routine builds the base preconditioner
!  at the current level.
!  The current level is regarded as the coarse one, while the previous as
!  the fine one. This is in agreement with the fact that the routine is called,
!  by mld_mlprec_bld, only on levels >=2.
!
! 
! Arguments:
!    a       -  type(psb_dspmat_type).
!               The sparse matrix structure containing the local part of the
!               matrix to be preconditioned.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    p       -  type(mld_zonelev_type), input/output.
!               The 'one-level' data structure containing the local part
!               of the base preconditioner to be built as well as
!               information concerning the prolongator and its transpose.
!    info    -  integer, output.
!               Error code.         
!  
subroutine mld_zcoarse_bld(a,desc_a,p,info)

  use psb_base_mod
  use mld_z_inner_mod, mld_protect_name => mld_zcoarse_bld

  implicit none

  ! Arguments
  type(psb_zspmat_type), intent(in), target :: a
  type(psb_desc_type), intent(in), target   :: desc_a
  type(mld_zonelev_type), intent(inout),target    :: p
  integer, intent(out)                      :: info

  ! Local variables
  character(len=20)      :: name
  integer                :: ictxt, np, me, err_act
  integer, allocatable   :: ilaggr(:), nlaggr(:)

  name='mld_zcoarse_bld'
  if (psb_get_errstatus().ne.0) return 
  call psb_erractionsave(err_act)
  info = psb_success_
  ictxt = desc_a%get_context()
  call psb_info(ictxt,me,np)

  call mld_check_def(p%parms%ml_type,'Multilevel type',&
       &   mld_mult_ml_,is_legal_ml_type)
  call mld_check_def(p%parms%aggr_alg,'Aggregation',&
       &   mld_dec_aggr_,is_legal_ml_aggr_alg)
  call mld_check_def(p%parms%aggr_kind,'Smoother',&
       &   mld_smooth_prol_,is_legal_ml_aggr_kind)
  call mld_check_def(p%parms%coarse_mat,'Coarse matrix',&
       &   mld_distr_mat_,is_legal_ml_coarse_mat)
  call mld_check_def(p%parms%aggr_filter,'Use filtered matrix',&
       &   mld_no_filter_mat_,is_legal_aggr_filter)
  call mld_check_def(p%parms%smoother_pos,'smooth_pos',&
       &   mld_pre_smooth_,is_legal_ml_smooth_pos)
  call mld_check_def(p%parms%aggr_omega_alg,'Omega Alg.',&
       &   mld_eig_est_,is_legal_ml_aggr_omega_alg)
  call mld_check_def(p%parms%aggr_eig,'Eigenvalue estimate',&
       &   mld_max_norm_,is_legal_ml_aggr_eig)
  call mld_check_def(p%parms%aggr_omega_val,'Omega',dzero,is_legal_omega)
  call mld_check_def(p%parms%aggr_thresh,'Aggr_Thresh',dzero,is_legal_aggr_thrs)


  !
  !  Build a mapping between the row indices of the fine-level matrix 
  !  and the row indices of the coarse-level matrix, according to a decoupled 
  !  aggregation algorithm. This also defines a tentative prolongator from
  !  the coarse to the fine level.
  ! 
  call mld_aggrmap_bld(p%parms%aggr_alg,p%parms%aggr_thresh,&
       & a,desc_a,ilaggr,nlaggr,info)

  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_aggrmap_bld')
    goto 9999
  end if

  !
  ! Build the coarse-level matrix from the fine-level one, starting from 
  ! the mapping defined by mld_aggrmap_bld and applying the aggregation
  ! algorithm specified by p%iprcparm(mld_aggr_kind_)
  !
  call mld_aggrmat_asb(a,desc_a,ilaggr,nlaggr,p,info)

  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_aggrmat_asb')
    goto 9999
  end if
  !
  ! Fix the base_a and base_desc pointers for handling of residuals.
  ! This is correct because this routine is only called at levels >=2.
  !
  p%base_a    => p%ac
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

end subroutine mld_zcoarse_bld
