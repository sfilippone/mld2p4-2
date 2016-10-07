!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.3)
!!$  
!!$  (C) Copyright 2008, 2010, 2012, 2015
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
! File: mld_d_base_aggregator_mat_asb.f90
!
! Subroutine: mld_d_base_aggregator_mat_asb
! Version:    real
!
!  This routine builds the matrix associated to the current level of the
!  multilevel preconditioner from the matrix associated to the previous level,
!  by using the user-specified aggregation technique (therefore, it also builds the
!  prolongation and restriction operators mapping the current level to the
!  previous one and vice versa). 
!  The current level is regarded as the coarse one, while the previous as
!  the fine one. This is in agreement with the fact that the routine is called,
!  by mld_mlprec_bld, only on levels >=2.
!  The main structure is:
!  1. Perform sanity checks;
!  2. Call mld_Xaggrmat_asb to compute prolongator/restrictor/AC
!  3. According to the choice of DIST/REPL for AC, build a descriptor DESC_AC,
!     and adjust the column numbering of AC/OP_PROL/OP_RESTR
!  4. Pack restrictor and prolongator into p%map
!  5. Fix base_a and base_desc pointers.
!
! 
! Arguments:
!    ag       -  type(mld_d_base_aggregator_type), input/output.
!               The aggregator object
!    parms   -  type(mld_dml_parms), input 
!               The aggregation parameters
!    a          -  type(psb_dspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!                  The 'one-level' data structure that will contain the local
!                  part of the matrix to be built as well as the information
!                  concerning the prolongator and its transpose.
!    ilaggr     -  integer, dimension(:), input
!                  The mapping between the row indices of the coarse-level
!                  matrix and the row indices of the fine-level matrix.
!                  ilaggr(i)=j means that node i in the adjacency graph
!                  of the fine-level matrix is mapped onto node j in the
!                  adjacency graph of the coarse-level matrix. Note that the indices
!                  are assumed to be shifted so as to make sure the ranges on
!                  the various processes do not   overlap.
!    nlaggr     -  integer, dimension(:) input
!                  nlaggr(i) contains the aggregates held by process i.
!    ac         -  type(psb_dspmat_type), output
!                  The coarse matrix on output 
!                  
!    op_prol    -  type(psb_dspmat_type), input/output
!                  The tentative prolongator on input, the computed prolongator on output
!               
!    op_restr    -  type(psb_dspmat_type), output
!                  The restrictor operator; normally, it is the transpose of the prolongator. 
!               
!    info       -  integer, output.
!                  Error code.
!  
subroutine  mld_d_base_aggregator_mat_asb(ag,parms,a,desc_a,ilaggr,nlaggr,ac,op_prol,op_restr,info)
  use psb_base_mod
  use mld_d_inner_mod, mld_protect_name => mld_d_base_aggregator_mat_asb
  implicit none
  
  class(mld_d_base_aggregator_type), target, intent(inout) :: ag
  type(mld_dml_parms), intent(inout)      :: parms 
  type(psb_dspmat_type), intent(in)    :: a
  type(psb_desc_type), intent(in)      :: desc_a
  integer(psb_ipk_), intent(inout)     :: ilaggr(:), nlaggr(:)
  type(psb_dspmat_type), intent(out)   :: ac,op_prol,op_restr
  integer(psb_ipk_), intent(out)       :: info

  ! Local variables
  character(len=20)             :: name
  integer(psb_mpik_)            :: ictxt, np, me
  integer(psb_ipk_)             :: err_act
  integer(psb_ipk_)            :: debug_level, debug_unit

  name='mld_d_base_aggregator_mat_asb'
  if (psb_get_errstatus().ne.0) return 
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  info  = psb_success_
  ictxt = desc_a%get_context()
  call psb_info(ictxt,me,np)

  call mld_check_def(parms%aggr_kind,'Smoother',&
       &   mld_smooth_prol_,is_legal_ml_aggr_kind)
  call mld_check_def(parms%coarse_mat,'Coarse matrix',&
       &   mld_distr_mat_,is_legal_ml_coarse_mat)
  call mld_check_def(parms%aggr_filter,'Use filtered matrix',&
       &   mld_no_filter_mat_,is_legal_aggr_filter)
  call mld_check_def(parms%smoother_pos,'smooth_pos',&
       &   mld_pre_smooth_,is_legal_ml_smooth_pos)
  call mld_check_def(parms%aggr_omega_alg,'Omega Alg.',&
       &   mld_eig_est_,is_legal_ml_aggr_omega_alg)
  call mld_check_def(parms%aggr_eig,'Eigenvalue estimate',&
       &   mld_max_norm_,is_legal_ml_aggr_eig)
  call mld_check_def(parms%aggr_omega_val,'Omega',dzero,is_legal_d_omega)

  !
  ! Build the coarse-level matrix from the fine-level one, starting from 
  ! the mapping defined by mld_aggrmap_bld and applying the aggregation
  ! algorithm specified by p%iprcparm(mld_aggr_kind_)
  !
  call mld_daggrmat_asb(a,desc_a,ilaggr,nlaggr,parms,ac,op_prol,op_restr,info)
  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

  
end subroutine mld_d_base_aggregator_mat_asb
