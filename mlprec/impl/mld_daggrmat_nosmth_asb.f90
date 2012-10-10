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
! File: mld_daggrmat_nosmth_asb.F90
!
! Subroutine: mld_daggrmat_nosmth_asb
! Version:    real
!
!  This routine builds a coarse-level matrix A_C from a fine-level matrix A
!  by using the Galerkin approach, i.e.
!
!                               A_C = P_C^T A P_C,
!
!  where P_C is the piecewise constant interpolation operator corresponding
!  the fine-to-coarse level mapping built by mld_aggrmap_bld.
! 
!  The coarse-level matrix A_C is distributed among the parallel processes or
!  replicated on each of them, according to the value of p%parms%coarse_mat
!  specified by the user through mld_dprecinit and mld_zprecset.
!
!  For details see
!    P. D'Ambra, D. di Serafino and  S. Filippone, On the development of
!    PSBLAS-based parallel two-level Schwarz preconditioners, Appl. Num. Math.,
!    57 (2007), 1181-1196.
!
!
!
! Arguments:
!    a          -  type(psb_dspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!    p          -  type(mld_d_onelev_type), input/output.
!                  The 'one-level' data structure that will contain the local
!                  part of the matrix to be built as well as the information 
!                  concerning the prolongator and its transpose.
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
subroutine mld_daggrmat_nosmth_asb(a,desc_a,ilaggr,nlaggr,parms,ac,op_prol,op_restr,info)
  use psb_base_mod
  use mld_d_inner_mod, mld_protect_name => mld_daggrmat_nosmth_asb

  implicit none

  ! Arguments
  type(psb_dspmat_type), intent(in)          :: a
  type(psb_desc_type), intent(in)            :: desc_a
  integer, intent(inout)                     :: ilaggr(:), nlaggr(:)
  type(mld_dml_parms), intent(inout)             :: parms 
  type(psb_dspmat_type), intent(out)             :: ac,op_prol,op_restr
  integer, intent(out)                       :: info

  ! Local variables
  integer :: ictxt,np,me, err_act
  integer(psb_mpik_) :: icomm, ndx, minfo
  character(len=20) :: name
  integer(psb_ipk_) :: ierr(5) 
  type(psb_d_coo_sparse_mat) :: ac_coo, acoo
  type(psb_d_csr_sparse_mat) :: acsr1, acsr2
  integer            :: debug_level, debug_unit
  integer :: nrow, nglob, ncol, ntaggr, nzl, ip, &
       & naggr, nzt, naggrm1, i

  name='mld_aggrmat_nosmth_asb'
  if(psb_get_errstatus().ne.0) return 
  info=psb_success_
  call psb_erractionsave(err_act)


  ictxt = desc_a%get_context()
  icomm = desc_a%get_mpic()
  call psb_info(ictxt, me, np)
  nglob = desc_a%get_global_rows()
  nrow  = desc_a%get_local_rows()
  ncol  = desc_a%get_local_cols()


  naggr  = nlaggr(me+1)
  ntaggr = sum(nlaggr)
  naggrm1=sum(nlaggr(1:me))

  do i=1, nrow
    ilaggr(i) = ilaggr(i) + naggrm1
  end do
  call psb_halo(ilaggr,desc_a,info)

  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_halo')
    goto 9999
  end if

  call acoo%allocate(ncol,ntaggr,ncol)

  do i=1,nrow
    acoo%val(i) = done
    acoo%ia(i)  = i
    acoo%ja(i)  = ilaggr(i)  
  end do

  call acoo%set_dupl(psb_dupl_add_)
  call acoo%set_nzeros(nrow)
  call acoo%set_asb()
  call acoo%fix(info)


  call op_prol%mv_from(acoo)
  call op_prol%cscnv(info,type='csr',dupl=psb_dupl_add_)
  if (info == psb_success_)   call op_prol%transp(op_restr)

  call a%csclip(ac_coo,info,jmax=nrow)

  nzt = ac_coo%get_nzeros()
  do i=1, nzt 
    ac_coo%ia(i) = ilaggr(ac_coo%ia(i))
    ac_coo%ja(i) = ilaggr(ac_coo%ja(i))
  enddo
  call ac_coo%set_nrows(naggr)
  call ac_coo%set_ncols(naggr)
  call ac_coo%set_dupl(psb_dupl_add_)
  call ac_coo%fix(info)
  call ac%mv_from(ac_coo) 


  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_daggrmat_nosmth_asb
