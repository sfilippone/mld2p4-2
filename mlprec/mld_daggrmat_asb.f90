!!$ 
!!$ 
!!$                                MLD2P4
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS v.2.0)
!!$  
!!$  (C) Copyright 2007  Alfredo Buttari      University of Rome Tor Vergata
!!$                      Pasqua D'Ambra       ICAR-CNR, Naples
!!$                      Daniela di Serafino  Second University of Naples
!!$                      Salvatore Filippone  University of Rome Tor Vergata       
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
! File: mld_daggrmat_asb.f90
!
! Subroutine: mld_daggrmat_asb
! Version:    real
!
!  This routine builds a coarse-level matrix A_C from a fine-level matrix A
!  by using a Galerkin approach, i.e.
!
!                               A_C = P_C^T A P_C,
!
!  where P_C is a prolongator from the coarse level to the fine one.
! 
!  A mapping from the nodes of the adjacency graph of A to the nodes of the
!  adjacency graph of A_C has been computed by the mld_aggrmap_bld subroutine.
!  The prolongator P_C is built here from this mapping,      according to the
!  value of p%iprcparm(mld_aggr_kind_), specified by the user through
!  mld_dprecinit and mld_dprecset.
!
!  Currently three different prolongators are implemented, corresponding to
!  three aggregation algorithms:
!  1. raw aggregation,
!  2. smoothed aggregation,
!  3. "bizarre" aggregation.
!  1. The raw aggregation uses as prolongator the piecewise constant interpolation
!     operator corresponding to the fine-to-coarse level mapping built by
!     mld_aggrmap_bld. This is called tentative prolongator.
!  2. The smoothed aggregation uses as prolongator the operator obtained by applying
!     a damped Jacobi smoother to the tentative prolongator.
!  3. The "bizarre" aggregation uses a prolongator proposed by the authors of MLD2P4.
!     This prolongator still requires a deep analysis and testing and its use is
!     not recommended.
!
!  For more details see
!    M. Brezina and P. Vanek, A black-box iterative solver based on a two-level
!    Schwarz method, Computing,  63 (1999), 233-263.
!    P. D'Ambra, D. di Serafino and S. Filippone, On the development of PSBLAS-based
!    parallel two-level Schwarz preconditioners, Appl. Num. Math., 57 (2007),
!    1181-1196.
!
!
! Arguments:
!    a          -  type(psb_dspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!    ac         -  type(psb_dspmat_type), output.
!                  The sparse matrix structure containing the local part of
!                  the coarse-level matrix.
!    desc_ac    -  type(psb_desc_type), output.
!                  The communication descriptor of the coarse-level matrix.
!    p          -  type(mld_dbaseprc_type), input/output.
!                  The base preconditioner data structure containing the local
!                  part of the base preconditioner to be built.
!    info       -  integer, output.
!                  Error code.
!
subroutine mld_daggrmat_asb(a,desc_a,ac,desc_ac,p,info)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_daggrmat_asb

  implicit none

! Arguments
  type(psb_dspmat_type), intent(in), target  :: a
  type(psb_desc_type), intent(in)            :: desc_a
  type(psb_dspmat_type), intent(inout), target :: ac    
  type(psb_desc_type), intent(inout)         :: desc_ac 
  type(mld_dbaseprc_type), intent(inout), target  :: p
  integer, intent(out)                       :: info

! Local variables
  logical, parameter :: aggr_dump=.false.
  integer ::ictxt,np,me, err_act, icomm
  character(len=20) :: name

  name='mld_aggrmat_asb'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  ictxt = psb_cd_get_context(desc_a)
  icomm = psb_cd_get_mpic(desc_a)

  call psb_info(ictxt, me, np)

  select case (p%iprcparm(mld_aggr_kind_))
  case (mld_no_smooth_) 

    call mld_aggrmat_raw_asb(a,desc_a,ac,desc_ac,p,info)

    if(info /= 0) then
      call psb_errpush(4010,name,a_err='raw_aggregate')
      goto 9999
    end if
    if (aggr_dump) call psb_csprt(90+me,ac,head='% Raw aggregate.')

  case(mld_smooth_prol_,mld_biz_prol_) 
    if (aggr_dump) call psb_csprt(70+me,a,head='% Input matrix')
    call mld_aggrmat_smth_asb(a,desc_a,ac,desc_ac,p,info)

    if(info /= 0) then
      call psb_errpush(4010,name,a_err='smooth_aggregate')
      goto 9999
    end if
    if (aggr_dump) call psb_csprt(90+me,ac,head='% Smooth aggregate.')
  case default
    call psb_errpush(4010,name,a_err=name)
    goto 9999

  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_daggrmat_asb
