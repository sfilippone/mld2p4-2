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
! File: mld_zaggrmap_bld.f90
!
! Subroutine: mld_zaggrmap_bld
! Version:    complex
!
!  This routine builds a mapping from the row indices of the fine-level matrix
!  to the row indices of the coarse-level matrix, according to a decoupled 
!  aggregation algorithm. This mapping will be used by mld_aggrmat_asb to
!  build the coarse-level matrix.  
!
!  The aggregation algorithm is a parallel version of that described in
!  * M. Brezina and P. Vanek, A black-box iterative solver based on a 
!    two-level Schwarz method, Computing,  63 (1999), 233-263.
!  * P. Vanek, J. Mandel and M. Brezina, Algebraic Multigrid by Smoothed
!    Aggregation for Second and Fourth Order Elliptic Problems, Computing, 56
!    (1996), 179-196.
!  For more details see
!    P. D'Ambra, D. di Serafino and S. Filippone, On the development of
!    PSBLAS-based parallel two-level Schwarz preconditioners, Appl. Num. Math.
!    57 (2007), 1181-1196.
!
!
! Arguments:
!    aggr_type  -  integer, input.
!                  The scalar used to identify the aggregation algorithm.
!    theta      -  real, input.
!                  The aggregation threshold used in the aggregation algorithm.
!    a          -  type(psb_zspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!    ilaggr     -  integer, dimension(:), allocatable.
!                  The mapping between the row indices of the coarse-level
!                  matrix and the row indices of the fine-level matrix.
!                  ilaggr(i)=j means that node i in the adjacency graph
!                  of the fine-level matrix is mapped onto node j in the
!                  adjacency graph of the coarse-level matrix.
!    nlaggr     -  integer, dimension(:), allocatable.
!                  nlaggr(i) contains the aggregates held by process i.
!    info       -  integer, output.
!                  Error code.
!
subroutine mld_zaggrmap_bld(aggr_type,theta,a,desc_a,ilaggr,nlaggr,info)

  use psb_base_mod
  use mld_z_inner_mod, mld_protect_name => mld_zaggrmap_bld
  
  implicit none

  ! Arguments
  integer, intent(in)               :: aggr_type
  real(psb_dpk_), intent(in)        :: theta
  type(psb_zspmat_type), intent(in) :: a
  type(psb_desc_type), intent(in)   :: desc_a
  integer, allocatable, intent(out) :: ilaggr(:),nlaggr(:)
  integer, intent(out)              :: info

  ! Local variables
  integer, allocatable  :: ils(:), neigh(:)
  integer :: icnt,nlp,k,n,ia,isz,nr, naggr,i,j,m
  type(psb_zspmat_type) :: atmp, atrans
  logical :: recovery
  integer :: debug_level, debug_unit
  integer :: ictxt,np,me,err_act
  integer :: nrow, ncol, n_ne
  character(len=20)  :: name, ch_err

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  name = 'mld_aggrmap_bld'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  !
  ictxt=desc_a%get_context()
  call psb_info(ictxt,me,np)
  nrow  = desc_a%get_local_rows()
  ncol  = desc_a%get_local_cols()

  select case (aggr_type)
  case (mld_dec_aggr_)  
    call mld_dec_map_bld(theta,a,desc_a,nlaggr,ilaggr,info)    

  case (mld_sym_dec_aggr_)  
    nr = a%get_nrows()
    call a%csclip(atmp,info,imax=nr,jmax=nr,&
         & rscale=.false.,cscale=.false.)
    call atmp%set_nrows(nr)
    call atmp%set_ncols(nr)
    if (info == psb_success_) call atrans%transp(atmp)
    if (info == psb_success_) call atrans%cscnv(info,type='COO')
    if (info == psb_success_) call psb_rwextd(nr,atmp,info,b=atrans,rowscale=.false.) 
    call atmp%set_nrows(nr)
    call atmp%set_ncols(nr)
    if (info == psb_success_) call atrans%free()
    if (info == psb_success_) call atmp%cscnv(info,type='CSR')
    if (info == psb_success_) call mld_dec_map_bld(theta,atmp,desc_a,nlaggr,ilaggr,info)    
    if (info == psb_success_) call atmp%free()

  case default

    info = -1
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=(/1,aggr_type,0,0,0/))
    goto 9999
  end select

  if (info /= psb_success_) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='dec_map_bld')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_zaggrmap_bld
