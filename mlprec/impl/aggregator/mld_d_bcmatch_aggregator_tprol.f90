!!$
!!$ 
!!$                           MLD2P4  version 2.1
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.3)
!!$  
!!$  (C) Copyright 2008, 2010, 2012, 2015, 2017 
!!$
!!$                      Salvatore Filippone  Cranfield University
!!$		      Ambra Abdullahi Hassan University of Rome Tor Vergata
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
! File: mld_d_bcmatch_aggregator_tprol.f90
!
! Subroutine: mld_d_bcmatch_aggregator_tprol
! Version:    real
!
!
!  This routine is mainly an interface to hyb_map_bld where the real work is performed. 
!  It takes care of some consistency checking, and calls map_to_tprol, which is
!  refactored and shared among all the aggregation methods that produce a simple
!  integer mapping.
!
! 
! Arguments:
!    p       -  type(mld_d_onelev_type), input/output.
!               The 'one-level' data structure containing the control
!               parameters and (eventually) coarse matrix and prolongator/restrictors. 
!               
!    a       -  type(psb_dspmat_type).
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
!    op_prol    -  type(psb_dspmat_type), output
!               The tentative prolongator, based on ilaggr.
!               
!    info    -  integer, output.
!               Error code.         
!

module bcm_CSRMatrix_mod
use psb_base_mod
use psb_util_mod
use iso_c_binding
use bcm_csr_type_mod
implicit none

contains
 subroutine MLD_to_CSR(a,csr_ia, csr_ja, csr_val, C, info)
 type(psb_dspmat_type), intent(in) :: a
 type(bcm_CSRMatrix), intent(out) :: C
 real(c_double), allocatable, target, intent(out) ::  csr_val(:)
 integer(c_int), allocatable, target, intent(out) :: csr_ia(:), csr_ja(:)

 !Local variable
 character(len=20)            :: name
 integer(psb_ipk_) :: info 
 real(psb_dpk_), allocatable :: coo_val(:)
 integer(psb_ipk_), allocatable :: coo_ia(:), coo_ja(:)
 integer(psb_ipk_)  :: x, nz, num_rows, num_cols , i , j, iad, k , k0
 type(psb_d_csr_sparse_mat) :: acsr

 name="MLD_to_CSR"
 num_rows= a%get_nrows()
 num_cols= a%get_ncols()
 nz= a%get_nzeros()
 call a%csgetrow(1,num_rows,nz,coo_ia,coo_ja,coo_val,info)
 !allocate(csr_ia(0:nz-1), csr_ja(0:nz-1), csr_val(0:nz-1), STAT=info)
 if (info /= psb_success_) then 
   info=psb_err_alloc_request_
   call psb_errpush(info,name,i_err=(/nz,izero,izero,izero,izero/),&
        & a_err='integer')
   return
 end if

 call a%cp_to(acsr)
 allocate(csr_ia(0:nz-1), csr_ja(0:nz-1), csr_val(0:nz-1), STAT=info)


 csr_ia(0:min(nz,size(acsr%irp,1))-1)=acsr%irp(1:min(nz,size(acsr%irp,1)))-1
 csr_ja(0:nz-1)=acsr%ja(1:nz)-1
 csr_val(0:nz-1)=acsr%val(1:nz)

 call acsr%free()

  C%num_rows=num_rows
  C%num_cols=num_cols
  C%num_nonzeros=nz
  C%owns_data=0
  C%i=c_loc(csr_ia)
  C%j=c_loc(csr_ja)
  C%data=c_loc(csr_val)
 end subroutine MLD_to_CSR

 subroutine bcm_to_op_prol(P, ilaggr, valaggr, info)
 type(bcm_CSRMatrix), intent(in) :: P
 integer(psb_ipk_), allocatable, intent(out) :: ilaggr(:)
 integer(psb_ipk_), intent(out) :: info
 real(psb_dpk_), allocatable, intent(out) :: valaggr(:)

 ! Local variables
 integer(psb_ipk_), pointer :: point_ia(:), point_ja(:)
 real(psb_dpk_), pointer :: point_val(:)
 integer(psb_ipk_) :: i, j, k
 integer(psb_ipk_) :: n, num_rows,  num_cols, num_nz
 character(len=20)            :: name
 integer(c_int), allocatable, target :: coo_ia(:),coo_ja(:)
 real(c_double), allocatable, target ::  coo_val(:)

  name="bcm_to_op_prol"
  call c_f_pointer(P%i,point_ia,(/1/))
  call c_f_pointer(P%j,point_ja,(/1/))
  call c_f_pointer(P%data,point_val,(/1/))

  !These are I, J, VAL.
  !These are I, J, VAL. 
  num_nz = P%num_nonzeros
  num_rows = P%num_rows
  num_cols = P%num_cols
  if (allocated(coo_ia)) deallocate(coo_ia)
  if (allocated(coo_ja)) deallocate(coo_ja)
  if (allocated(coo_val)) deallocate(coo_val)
  allocate(coo_ia(num_nz),coo_ja(num_nz),coo_val(num_nz), STAT=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/num_nz,izero,izero,izero,izero/),&
         & a_err='integer')
    return
  end if

  n = 1
  !coo_ia=-123
  !coo_ja=-123
  coo_val=point_val(1:num_nz)
  do i=1, num_rows
    do j=point_ia(i)+1, point_ia(i+1)
	coo_ia(n)=i
	coo_ja(n)=point_ja(j) + 1
	n = n + 1
    enddo
  enddo
  if (allocated(ilaggr)) deallocate(ilaggr)
  if (allocated(valaggr)) deallocate(valaggr)
  allocate(ilaggr(num_rows),valaggr(num_rows), STAT=info)
  ilaggr=0
  if (info /= psb_success_) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/num_rows,izero,izero,izero,izero/),&
         & a_err='integer')
    return
  end if

  do k=1,num_nz
    i=coo_ia(k)
    j=coo_ja(k)
    ilaggr(i)=j
    valaggr(i)=coo_val(i)
  enddo
  if (allocated(coo_ia)) deallocate(coo_ia)
  if (allocated(coo_ja)) deallocate(coo_ja)
  if (allocated(coo_val)) deallocate(coo_val)
  nullify(point_ia)
  nullify(point_ja)
  nullify(point_val)
 end subroutine bcm_to_op_prol
end module bcm_CSRMatrix_mod

  
subroutine  mld_d_bcmatch_aggregator_build_tprol(ag,parms,a,desc_a,ilaggr,nlaggr,op_prol,info)
  use psb_base_mod
  use mld_d_bcmatch_aggregator_mod, mld_protect_name => mld_d_bcmatch_aggregator_build_tprol
  use mld_d_inner_mod
  use bcm_CSRMatrix_mod 
  use iso_c_binding      
  implicit none
  class(mld_d_bcmatch_aggregator_type), target, intent(inout) :: ag
  type(mld_dml_parms), intent(inout)  :: parms 
  type(psb_dspmat_type), intent(in)   :: a
  type(psb_desc_type), intent(in)     :: desc_a
  integer(psb_ipk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
  real(psb_dpk_), allocatable:: valaggr(:)
  type(psb_dspmat_type), intent(out)  :: op_prol
  integer(psb_ipk_), intent(out)      :: info


  ! Local variables
  type(psb_dspmat_type)   :: a_tmp
  type(bcm_CSRMatrix) :: C, P
  integer(c_int) :: match_algorithm, n_sweeps, max_csize, max_nlevels
  character(len=20)            :: name, ch_err
  integer(psb_mpik_)           :: ictxt, np, me
  integer(psb_ipk_)            :: err_act, ierr
  integer(psb_ipk_)            :: debug_level, debug_unit
  integer(psb_ipk_)            :: i, j, k
  integer(psb_ipk_), allocatable, target ::  csr_ia(:), csr_ja(:)
  integer(psb_ipk_), allocatable :: aux(:)
  real(psb_dpk_), allocatable, target::  csr_val(:)
  interface
     function bootCMatch(C,match_alg,n_sweeps,max_nlevels,max_csize,w)bind(c,name='bootCMatch') result(P)
       use iso_c_binding  
       use bcm_CSRMatrix_mod         
       implicit none
       type(bcm_CSRMatrix) :: C, P
       type(bcm_Vector) :: w
       integer(c_int) :: match_alg
       integer(c_int) :: n_sweeps
       integer(c_int) :: max_nlevels
       integer(c_int) :: max_csize
     end function bootCMatch
  end interface

  name='mld_d_bcmatch_aggregator_tprol'
  ictxt = desc_a%get_context()
  call psb_info(ictxt,me,np)
  if (psb_get_errstatus().ne.0) return 
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  info  = psb_success_


  call mld_check_def(parms%ml_type,'Multilevel type',&
       &   mld_mult_ml_,is_legal_ml_type)
  call mld_check_def(parms%aggr_alg,'Aggregation',&
       &   mld_dec_aggr_,is_legal_ml_aggr_alg)
  call mld_check_def(parms%aggr_ord,'Ordering',&
       &   mld_aggr_ord_nat_,is_legal_ml_aggr_ord)
  call mld_check_def(parms%aggr_thresh,'Aggr_Thresh',dzero,is_legal_d_aggr_thrs)

  call a%csclip(b=a_tmp, info=info, jmax=a%get_nrows(), imax=a%get_nrows())

  call MLD_to_CSR(a_tmp,csr_ia, csr_ja, csr_val, C, info)

  call a_tmp%free()

  match_algorithm=ag%matching_alg
  n_sweeps=ag%n_sweeps
  max_csize=ag%max_csize
  max_nlevels=ag%max_nlevels
  P = bootCMatch(C, match_algorithm, n_sweeps, max_nlevels, max_csize, ag%w_par)

  call bcm_to_op_prol(P, ilaggr, valaggr, info)

  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='bcm_to_op_prol')
    goto 9999
  end if
  call psb_realloc(a%get_ncols(),ilaggr,info)
  if (info /= psb_success_) then 
    info=psb_err_from_subroutine_
    ch_err='psb_realloc'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call psb_realloc(a%get_ncols(),valaggr,info)
  if (info /= psb_success_) then 
    info=psb_err_from_subroutine_
    ch_err='psb_realloc'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (allocated(nlaggr)) deallocate(nlaggr)
  allocate(nlaggr(np), STAT=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/np,izero,izero,izero,izero/),&
         & a_err='integer')
    goto 9999
  end if

  nlaggr(:)=0
  nlaggr(me+1) = P%num_cols
  call psb_sum(ictxt,nlaggr(1:np))

  call mld_bcmatch_map_to_tprol(desc_a,ilaggr,nlaggr,valaggr,op_prol,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_bcmatch_map_to_tprol')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
  
end subroutine mld_d_bcmatch_aggregator_build_tprol
