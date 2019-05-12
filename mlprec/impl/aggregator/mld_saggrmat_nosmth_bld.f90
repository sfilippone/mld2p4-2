!   
!   
!                             MLD2P4  version 2.2
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!    
!    (C) Copyright 2008-2018 
!  
!        Salvatore Filippone  
!        Pasqua D'Ambra   
!        Daniela di Serafino   
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
! File: mld_saggrmat_nosmth_bld.F90
!
! Subroutine: mld_saggrmat_nosmth_bld
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
!  specified by the user through mld_sprecinit and mld_zprecset.
!  On output from this routine the entries of AC, op_prol, op_restr
!  are still in "global numbering" mode; this is fixed in the calling routine
!  aggregator%mat_bld.
!
!  For details see
!    P. D'Ambra, D. di Serafino and  S. Filippone, On the development of
!    PSBLAS-based parallel two-level Schwarz preconditioners, Appl. Num. Math.,
!    57 (2007), 1181-1196.
!
!
! Arguments:
!    a          -  type(psb_sspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!    p          -  type(mld_s_onelev_type), input/output.
!                  The 'one-level' data structure that will contain the local
!                  part of the matrix to be built as well as the information
!                  concerning the prolongator and its transpose.
!    parms      -   type(mld_sml_parms), input
!                  Parameters controlling the choice of algorithm
!    ac         -  type(psb_sspmat_type), output
!                  The coarse matrix on output 
!                  
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
!    op_prol    -  type(psb_sspmat_type), input/output
!                  The tentative prolongator on input, the computed prolongator on output
!               
!    op_restr    -  type(psb_sspmat_type), output
!                  The restrictor operator; normally, it is the transpose of the prolongator. 
!               
!    info       -  integer, output.
!                  Error code.
!
!
subroutine mld_saggrmat_nosmth_bld(a,desc_a,ilaggr,nlaggr,parms,ac,op_prol,op_restr,info)
  use psb_base_mod
  use mld_base_prec_type
  use mld_s_inner_mod, mld_protect_name => mld_saggrmat_nosmth_bld

  implicit none

  ! Arguments
  type(psb_sspmat_type), intent(in)        :: a
  type(psb_desc_type), intent(in)            :: desc_a
  integer(psb_ipk_), intent(inout)           :: ilaggr(:), nlaggr(:)
  type(mld_sml_parms), intent(inout)      :: parms 
  type(psb_sspmat_type), intent(inout)     :: op_prol
  type(psb_sspmat_type), intent(out)       :: ac,op_restr
  integer(psb_ipk_), intent(out)             :: info

  ! Local variables
  integer(psb_ipk_)  :: err_act
  integer(psb_ipk_)  :: ictxt,np,me, icomm, ndx, minfo
  character(len=20)  :: name
  integer(psb_ipk_)  :: ierr(5) 
  type(psb_s_coo_sparse_mat) :: ac_coo, tmpcoo
  type(psb_s_csr_sparse_mat) :: acsr1, acsr2
  integer(psb_ipk_) :: debug_level, debug_unit
  integer(psb_ipk_) :: nrow, nglob, ncol, ntaggr, nzl, ip, &
       & naggr, nzt, naggrm1, naggrp1, i, k

  name='mld_aggrmat_nosmth_bld'
  if(psb_get_errstatus().ne.0) return 
  info=psb_success_
  call psb_erractionsave(err_act)


  ictxt = desc_a%get_context()
  icomm = desc_a%get_mpic()
  call psb_info(ictxt, me, np)
  nglob = desc_a%get_global_rows()
  nrow  = desc_a%get_local_rows()
  ncol  = desc_a%get_local_cols()


  naggr   = nlaggr(me+1)
  ntaggr  = sum(nlaggr)
  naggrm1 = sum(nlaggr(1:me))
  naggrp1 = sum(nlaggr(1:me+1))

  
  call op_prol%cp_to(tmpcoo)
  call op_prol%cscnv(info,type='csr',dupl=psb_dupl_add_)
  
  call tmpcoo%transp()
  nzl = tmpcoo%get_nzeros()
  i=0
  !
  ! Now we have to fix this.  The only rows of tmpcoo/op_restr that are correct 
  ! are those corresponding to "local" aggregates, i.e. indices in ilaggr(:)
  !
  do k=1, nzl
    if ((naggrm1 < tmpcoo%ia(k)) .and.(tmpcoo%ia(k) <= naggrp1)) then
      i = i+1
      tmpcoo%val(i) = tmpcoo%val(k)
      tmpcoo%ia(i)  = tmpcoo%ia(k)
      tmpcoo%ja(i)  = tmpcoo%ja(k)
    end if
  end do
  call tmpcoo%set_nzeros(i)
  !  call tmpcoo%trim()
  call op_restr%mv_from(tmpcoo)
  call op_restr%cscnv(info,type='csr',dupl=psb_dupl_add_)

  if (info /= psb_success_) goto 9999

  call a%cp_to(ac_coo)

  nzt = ac_coo%get_nzeros()
  k = 0
  do i=1, nzt 
    k = k + 1 
    ac_coo%ia(k)  = ilaggr(ac_coo%ia(i))
    ac_coo%ja(k)  = ilaggr(ac_coo%ja(i))
    ac_coo%val(k) = ac_coo%val(i)
    ! At this point, there may be negative entries,
    ! because that's how ILAGGR marks singletons
    ! If this is the case, roll back K
    if ((ac_coo%ia(k)<=0).or.(ac_coo%ja(k)<=0)) k = k-1      
  enddo
  call ac_coo%set_nrows(naggr)
  call ac_coo%set_ncols(naggr)
  call ac_coo%set_nzeros(k)
  call ac_coo%set_dupl(psb_dupl_add_)
  call ac_coo%fix(info)
  call ac%mv_from(ac_coo) 


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine mld_saggrmat_nosmth_bld
