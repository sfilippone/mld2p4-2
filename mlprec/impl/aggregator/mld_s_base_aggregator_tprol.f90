
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
! File: mld_s_base_aggregator_tprol.f90
!
! Subroutine: mld_s_base_aggregator_tprol
! Version:    real
!
!  This routine is mainly an interface to dec_map_bld where the real work is performed. 
!  It takes care of some consistency checking, and calls map_to_tprol, which is
!  refactored and shared among all the aggregation methods that produce a simple
!  integer mapping.
!
! 
! Arguments:
!    p       -  type(mld_s_onelev_type), input/output.
!               The 'one-level' data structure containing the control
!               parameters and (eventually) coarse matrix and prolongator/restrictors. 
!               
!    a       -  type(psb_sspmat_type).
!               The sparse matrix structure containing the local part of the
!               fine-level matrix.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    ilaggr     -  integer, dimension(:), allocatable, output
!                  The mapping between the row indices of the coarse-level
!                  matrix and the row indices of the fine-level matrix.
!                  ilaggr(i)=j means that node i in the adjacency graph
!                  of the fine-level matrix is mapped onto node j in the
!                  adjacency graph of the coarse-level matrix. Note that on exit the indices
!                  will be shifted so as to make sure the ranges on the various processes do not
!                  overlap.
!    nlaggr     -  integer, dimension(:), allocatable, output
!                  nlaggr(i) contains the aggregates held by process i.
!    op_prol    -  type(psb_sspmat_type), output
!               The tentative prolongator, based on ilaggr.
!               
!    info    -  integer, output.
!               Error code.         
!  
subroutine  mld_s_base_aggregator_build_tprol(ag,parms,a,desc_a,ilaggr,nlaggr,op_prol,info)
  use psb_base_mod
  use mld_s_prec_type, mld_protect_name => mld_s_base_aggregator_build_tprol
  use mld_s_inner_mod
  implicit none
  class(mld_s_base_aggregator_type), target, intent(inout) :: ag
  type(mld_sml_parms), intent(inout)  :: parms 
  type(psb_sspmat_type), intent(in)   :: a
  type(psb_desc_type), intent(in)     :: desc_a
  integer(psb_ipk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
  type(psb_sspmat_type), intent(out)  :: op_prol
  integer(psb_ipk_), intent(out)      :: info

  ! Local variables
  character(len=20)            :: name
  integer(psb_mpik_)           :: ictxt, np, me
  integer(psb_ipk_)            :: err_act
  integer(psb_ipk_)            :: ntaggr
  integer(psb_ipk_)            :: debug_level, debug_unit

  name='mld_s_base_aggregator_tprol'
  if (psb_get_errstatus().ne.0) return 
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  info  = psb_success_
  ictxt = desc_a%get_context()
  call psb_info(ictxt,me,np)

  call mld_check_def(parms%ml_cycle,'Multilevel cycle',&
       &   mld_mult_ml_,is_legal_ml_cycle)
  call mld_check_def(parms%par_aggr_alg,'Aggregation',&
       &   mld_dec_aggr_,is_legal_ml_par_aggr_alg)
  call mld_check_def(parms%aggr_ord,'Ordering',&
       &   mld_aggr_ord_nat_,is_legal_ml_aggr_ord)
  call mld_check_def(parms%aggr_thresh,'Aggr_Thresh',szero,is_legal_s_aggr_thrs)


  call mld_dec_map_bld(parms%aggr_ord,parms%aggr_thresh,a,desc_a,nlaggr,ilaggr,info)

  if (info==psb_success_) call mld_map_to_tprol(desc_a,ilaggr,nlaggr,op_prol,info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='dec_map_bld/map_to_tprol')
    goto 9999
  endif
  
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
  
end subroutine mld_s_base_aggregator_build_tprol