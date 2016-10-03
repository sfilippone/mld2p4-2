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
! File: mld_dcoarse_bld.f90
!
! Subroutine: mld_dcoarse_bld
! Version:    real
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
!               fine-level matrix.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    p       -  type(mld_d_onelev_type), input/output.
!               The 'one-level' data structure containing the local part
!               of the base preconditioner to be built as well as
!               information concerning the prolongator and its transpose.
!    info    -  integer, output.
!               Error code.         
!  
subroutine mld_d_lev_aggrmap_bld(p,a,desc_a,ilaggr,nlaggr,op_prol,info)

  use psb_base_mod
  use mld_d_inner_mod, mld_protect_name => mld_d_lev_aggrmap_bld 

  implicit none

  ! Arguments
  type(mld_d_onelev_type), intent(inout), target :: p
  type(psb_dspmat_type), intent(in) :: a
  type(psb_desc_type), intent(in)     :: desc_a
  integer(psb_ipk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
  type(psb_dspmat_type), intent(out)  :: op_prol
  integer(psb_ipk_), intent(out)      :: info
  

  ! Local variables
  character(len=20)                :: name
  integer(psb_mpik_)               :: ictxt, np, me
  integer(psb_ipk_)                :: err_act
  integer(psb_ipk_)                :: nzl, ntaggr
  integer(psb_ipk_)            :: debug_level, debug_unit

  name='mld_d_lev_aggrmap_bld'
  if (psb_get_errstatus().ne.0) return 
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  info  = psb_success_
  ictxt = desc_a%get_context()
  call psb_info(ictxt,me,np)

  call mld_check_def(p%parms%ml_type,'Multilevel type',&
       &   mld_mult_ml_,is_legal_ml_type)
  call mld_check_def(p%parms%aggr_alg,'Aggregation',&
       &   mld_dec_aggr_,is_legal_ml_aggr_alg)
  call mld_check_def(p%parms%aggr_ord,'Ordering',&
       &   mld_aggr_ord_nat_,is_legal_ml_aggr_ord)
  call mld_check_def(p%parms%aggr_thresh,'Aggr_Thresh',dzero,is_legal_d_aggr_thrs)

  select case(p%parms%aggr_alg)
  case (mld_dec_aggr_, mld_sym_dec_aggr_)  
    
    !
    !  Build a mapping between the row indices of the fine-level matrix 
    !  and the row indices of the coarse-level matrix, according to a decoupled 
    !  aggregation algorithm. This also defines a tentative prolongator from
    !  the coarse to the fine level.
    ! 
    call mld_aggrmap_bld(p%parms%aggr_alg,p%parms%aggr_ord,p%parms%aggr_thresh,&
         & a,desc_a,ilaggr,nlaggr,op_prol,info)
    
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_aggrmap_bld')
      goto 9999
    end if

  case (mld_bcmatch_aggr_)
    write(0,*) 'Matching is not implemented yet '
    info = -1111
    call psb_errpush(psb_err_input_value_invalid_i_,name,&
         & i_err=(/ione,p%parms%aggr_alg,izero,izero,izero/))
    goto 9999
    
  case default

    info = -1
    call psb_errpush(psb_err_input_value_invalid_i_,name,&
         & i_err=(/ione,p%parms%aggr_alg,izero,izero,izero/))
    goto 9999

  end select

  
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine mld_d_lev_aggrmap_bld
