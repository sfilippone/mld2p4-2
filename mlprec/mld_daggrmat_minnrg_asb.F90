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
! File: mld_daggrmat_minnrg_asb.F90
!
! Subroutine: mld_daggrmat_minnrg_asb
! Version:    real
!
!  This routine builds a coarse-level matrix A_C from a fine-level matrix A
!  by using the Galerkin approach, i.e.
!
!                               A_C = P_C^T A P_C,
!
!  where P_C is a prolongator from the coarse level to the fine one.
! 
!  The prolongator P_C is built according to a smoothed aggregation algorithm,
!  i.e. it is obtained by applying a damped Jacobi smoother to the piecewise
!  constant interpolation operator P corresponding to the fine-to-coarse level 
!  mapping built by the mld_aggrmap_bld subroutine:
!
!                            P_C = (I - omega*D^(-1)A) * P,
!
!  where D is the diagonal matrix with main diagonal equal to the main diagonal
!  of A, and omega is a suitable smoothing parameter. An estimate of the spectral
!  radius of D^(-1)A, to be used in the computation of omega, is provided, 
!  according to the value of p%parms%aggr_omega_alg, specified by the user
!  through mld_dprecinit and mld_dprecset.
!
!  This routine can also build A_C according to a "bizarre" aggregation algorithm,
!  using a "naive" prolongator proposed by the authors of MLD2P4. However, this
!  prolongator still requires a deep analysis and testing and its use is not
!  recommended.
!
!  The coarse-level matrix A_C is distributed among the parallel processes or
!  replicated on each of them, according to the value of p%parms%coarse_mat,
!  specified by the user through mld_dprecinit and mld_dprecset.
!
!  For more details see
!    M. Brezina and P. Vanek, A black-box iterative solver based on a 
!    two-level Schwarz method, Computing,  63 (1999), 233-263.
!    P. D'Ambra, D. di Serafino and S. Filippone, On the development of
!    PSBLAS-based parallel two-level Schwarz preconditioners, Appl. Num. Math.
!    57 (2007), 1181-1196.
!
! Arguments:
!    a          -  type(psb_dspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!    p          -  type(mld_donelev_type), input/output.
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
subroutine mld_daggrmat_minnrg_asb(a,desc_a,ilaggr,nlaggr,p,info)
  use psb_sparse_mod
  use mld_d_inner_mod, mld_protect_name => mld_daggrmat_minnrg_asb

#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

  ! Arguments
  type(psb_dspmat_type), intent(in)            :: a
  type(psb_desc_type), intent(in)               :: desc_a
  integer, intent(inout)                        :: ilaggr(:), nlaggr(:)
  type(mld_donelev_type), intent(inout), target :: p
  integer, intent(out)                          :: info

  ! Local variables
  type(psb_dspmat_type)  :: b
  integer, allocatable :: nzbr(:), idisp(:)
  integer :: nrow, nglob, ncol, ntaggr, nzac, ip, ndx,&
       & naggr, nzl,naggrm1,naggrp1, i, j, k, jd, icolF, nrt
  integer ::ictxt,np,me, err_act, icomm
  character(len=20) :: name
!!$  type(psb_dspmat_type) :: am1,am2, af, ptilde, rtilde, atran, atp, atdatp
!!$  type(psb_dspmat_type) :: am3,am4, ap, adap,atmp,rada, ra, atmp2
  real(psb_dpk_), allocatable :: adiag(:), pj(:), xj(:), yj(:), omf(:),omp(:),omi(:),&
       & oden(:), adinv(:)
  logical            :: filter_mat
  integer            :: debug_level, debug_unit
  integer, parameter :: ncmax=16
  real(psb_dpk_)   :: omega, anorm, tmp, dg, theta, alpha,beta, ommx

  name='mld_aggrmat_minnrg'
  if(psb_get_errstatus().ne.0) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = psb_cd_get_context(desc_a)
  icomm = psb_cd_get_mpic(desc_a)
  ictxt = psb_cd_get_context(desc_a)

  call psb_info(ictxt, me, np)
!!$
!!$
!!$  call psb_nullify_sp(b)
!!$  call psb_nullify_sp(am3)
!!$  call psb_nullify_sp(am4)
!!$  call psb_nullify_sp(am1)
!!$  call psb_nullify_sp(am2)
!!$  call psb_nullify_sp(Ap)
!!$  call psb_nullify_sp(Adap)
!!$  call psb_nullify_sp(Atmp)
!!$  call psb_nullify_sp(Atmp2)
!!$  call psb_nullify_sp(Atran)
!!$  call psb_nullify_sp(Atp)
!!$  call psb_nullify_sp(atdatp)
!!$  call psb_nullify_sp(AF)
!!$  call psb_nullify_sp(ra)
!!$  call psb_nullify_sp(rada)
!!$  call psb_nullify_sp(ptilde)
!!$  call psb_nullify_sp(rtilde)
!!$
!!$  nglob = psb_cd_get_global_rows(desc_a)
!!$  nrow  = psb_cd_get_local_rows(desc_a)
!!$  ncol  = psb_cd_get_local_cols(desc_a)
!!$
!!$  theta = p%rprcparm(mld_aggr_thresh_)
!!$
!!$  naggr  = nlaggr(me+1)
!!$  ntaggr = sum(nlaggr)
!!$
!!$  allocate(nzbr(np), idisp(np),stat=info)
!!$  if (info /= psb_success_) then 
!!$    info=psb_err_alloc_request_
!!$    call psb_errpush(info,name,i_err=(/2*np,0,0,0,0/),&
!!$         & a_err='integer')
!!$    goto 9999      
!!$  end if
!!$
!!$  naggrm1 = sum(nlaggr(1:me))
!!$  naggrp1 = sum(nlaggr(1:me+1))
!!$
!!$  filter_mat = (p%parms%aggr_filter == mld_filter_mat_)
!!$
!!$  ilaggr(1:nrow) = ilaggr(1:nrow) + naggrm1
!!$  call psb_halo(ilaggr,desc_a,info)
!!$
!!$  if (info /= psb_success_) then
!!$    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_halo')
!!$    goto 9999
!!$  end if
!!$
!!$  ! naggr: number of local aggregates
!!$  ! nrow: local rows. 
!!$  ! 
!!$  allocate(adiag(ncol),adinv(ncol),xj(ncol),&
!!$       & yj(ncol),omf(ncol),omp(ntaggr),oden(ntaggr),omi(ncol),stat=info)
!!$
!!$  if (info /= psb_success_) then 
!!$    info=psb_err_alloc_request_
!!$    call psb_errpush(info,name,i_err=(/6*ncol+ntaggr,0,0,0,0/),&
!!$         & a_err='real(psb_dpk_)')
!!$    goto 9999      
!!$  end if
!!$
!!$  ! Get the diagonal D
!!$  call psb_sp_getdiag(a,adiag,info)
!!$  if (info == psb_success_) &
!!$       & call psb_halo(adiag,desc_a,info)
!!$
!!$  if(info /= psb_success_) then
!!$    call psb_errpush(psb_err_from_subroutine_,name,a_err='sp_getdiag')
!!$    goto 9999
!!$  end if
!!$
!!$  ! 1. Allocate Ptilde in sparse matrix form 
!!$  ptilde%fida='COO'
!!$  ptilde%m=ncol
!!$  ptilde%k=ntaggr
!!$  call psb_sp_all(ncol,ntaggr,ptilde,ncol,info)
!!$
!!$
!!$  if (info /= psb_success_) then
!!$    call psb_errpush(psb_err_from_subroutine_,name,a_err='spall')
!!$    goto 9999
!!$  end if
!!$
!!$  do i=1,ncol
!!$    ptilde%aspk(i) = done
!!$    ptilde%ia1(i)  = i
!!$    ptilde%ia2(i)  = ilaggr(i)  
!!$  end do
!!$  ptilde%infoa(psb_nnz_) = ncol
!!$
!!$  call psb_spcnv(ptilde,info,afmt='csr',dupl=psb_dupl_add_)
!!$  if (info == psb_success_) call psb_spcnv(a,am3,info,afmt='csr',dupl=psb_dupl_add_)
!!$  if (info /= psb_success_) then
!!$    call psb_errpush(psb_err_from_subroutine_,name,a_err='spcnv')
!!$    goto 9999
!!$  end if
!!$  if (debug_level >= psb_debug_outer_) &
!!$       & write(debug_unit,*) me,' ',trim(name),&
!!$       & ' Initial copies done.'
!!$
!!$  call psb_symbmm(am3,ptilde,ap,info)
!!$  if (info == psb_success_) call psb_numbmm(am3,ptilde,ap)
!!$
!!$  if(info /= psb_success_) then
!!$    call psb_errpush(psb_err_from_subroutine_,name,a_err='symbmm 1')
!!$    goto 9999
!!$  end if
!!$
!!$  call psb_sp_clone(ap,atmp,info)
!!$
!!$
!!$  do i=1,size(adiag)
!!$    if (adiag(i) /= dzero) then
!!$      adinv(i) = done / adiag(i)
!!$    else
!!$      adinv(i) = done
!!$    end if
!!$  end do
!!$  call psb_sp_scal(adinv,atmp,info)
!!$  call psb_sphalo(atmp,desc_a,am4,info,&
!!$       & colcnv=.false.,rowscale=.true.,outfmt='CSR  ')
!!$  if (info == psb_success_) call psb_rwextd(ncol,atmp,info,b=am4)      
!!$  if (info == psb_success_) call psb_sp_free(am4,info)
!!$
!!$  call psb_symbmm(am3,atmp,adap,info)
!!$  call psb_numbmm(am3,atmp,adap)
!!$  call psb_sp_free(atmp,info)
!!$
!!$! !$  write(0,*) 'Columns of AP',psb_sp_get_ncols(ap)
!!$! !$  write(0,*) 'Columns of ADAP',psb_sp_get_ncols(adap)
!!$  call psb_spcnv(ap,info,afmt='coo')
!!$  if (info == psb_success_) call psb_spcnv(ap,info,afmt='csc')
!!$  if (info == psb_success_) call psb_spcnv(adap,info,afmt='coo')
!!$  if (info == psb_success_) call psb_spcnv(adap,info,afmt='csc')
!!$  if (info /= psb_success_) then 
!!$    write(0,*) 'Failed conversion to CSC'
!!$  end if
!!$
!!$  call csc_mat_col_prod(ap,adap,omp,info)
!!$  call csc_mat_col_prod(adap,adap,oden,info)
!!$  call psb_sum(ictxt,omp)
!!$  call psb_sum(ictxt,oden)
!!$! !$  write(debug_unit,*) trim(name),' OMP :',omp
!!$! !$  write(debug_unit,*) trim(name),' ODEN:',oden
!!$  omp = omp/oden
!!$! !$  write(0,*) 'Check on output prolongator ',omp(1:min(size(omp),10))
!!$  if (debug_level >= psb_debug_outer_) &
!!$       & write(debug_unit,*) me,' ',trim(name),&
!!$       & 'Done NUMBMM 1'
!!$
!!$  ! Compute omega_int
!!$  ommx = -1d300
!!$  do i=1, ncol
!!$    omi(i) = omp(ilaggr(i))
!!$    ommx = max(ommx,omi(i))
!!$  end do
!!$  ! Compute omega_fine
!!$  do i=1, nrow
!!$    omf(i) = ommx
!!$    do j=am3%ia2(i),am3%ia2(i+1)-1
!!$      omf(i) = min(omf(i),omi(am3%ia1(j)))
!!$    end do
!!$    omf(i) = max(dzero,omf(i))
!!$  end do
!!$
!!$
!!$  if (filter_mat) then
!!$    !
!!$    ! Build the filtered matrix Af from A
!!$    ! 
!!$    call psb_spcnv(a,af,info,afmt='csr',dupl=psb_dupl_add_)
!!$
!!$    do i=1,nrow
!!$      tmp = dzero
!!$      jd  = -1 
!!$      do j=af%ia2(i),af%ia2(i+1)-1
!!$        if (af%ia1(j) == i) jd = j 
!!$        if (abs(af%aspk(j)) < theta*sqrt(abs(adiag(i)*adiag(af%ia1(j))))) then
!!$          tmp=tmp+af%aspk(j)
!!$          af%aspk(j)=dzero
!!$        endif
!!$      enddo
!!$      if (jd == -1) then 
!!$        write(0,*) 'Wrong input: we need the diagonal!!!!', i
!!$      else
!!$        af%aspk(jd)=af%aspk(jd)-tmp
!!$      end if
!!$    enddo
!!$    ! Take out zeroed terms 
!!$    call psb_spcnv(af,info,afmt='coo')
!!$    k = 0
!!$    do j=1,psb_sp_get_nnzeros(af)
!!$      if ((af%aspk(j) /= dzero) .or. (af%ia1(j) == af%ia2(j))) then 
!!$        k = k + 1
!!$        af%aspk(k) = af%aspk(j)
!!$        af%ia1(k)  = af%ia1(j)
!!$        af%ia2(k)  = af%ia2(j)
!!$      end if
!!$    end do
!!$! !$  write(debug_unit,*) me,' ',trim(name),' Non zeros from filtered matrix:',k,af%m,af%k
!!$    call psb_sp_setifld(k,psb_nnz_,af,info)
!!$    call psb_spcnv(af,info,afmt='csr')
!!$  end if
!!$
!!$  omf(1:nrow) = omf(1:nrow) * adinv(1:nrow)
!!$
!!$  if (filter_mat) then
!!$    !
!!$    ! Build the smoothed prolongator using the filtered matrix
!!$    ! 
!!$    if (psb_toupper(af%fida) == 'CSR') then 
!!$      do i=1,af%m
!!$        do j=af%ia2(i),af%ia2(i+1)-1
!!$          if (af%ia1(j) == i) then 
!!$            af%aspk(j) = done - omf(i)*af%aspk(j) 
!!$          else
!!$            af%aspk(j) = - omf(i)*af%aspk(j) 
!!$          end if
!!$        end do
!!$      end do
!!$    else 
!!$      call psb_errpush(psb_err_internal_error_,name,a_err='Invalid AF storage format')
!!$      goto 9999
!!$    end if
!!$
!!$    if (debug_level >= psb_debug_outer_) &
!!$         & write(debug_unit,*) me,' ',trim(name),&
!!$         & 'Done gather, going for SYMBMM 1'
!!$    !
!!$    ! Symbmm90 does the allocation for its result.
!!$    ! 
!!$    ! am1 = (I-w*D*Af)Ptilde
!!$    ! Doing it this way means to consider diag(Af_i)
!!$    ! 
!!$    !
!!$    call psb_symbmm(af,ptilde,am1,info)
!!$    if(info /= psb_success_) then
!!$      call psb_errpush(psb_err_from_subroutine_,name,a_err='symbmm 1')
!!$      goto 9999
!!$    end if
!!$
!!$    call psb_numbmm(af,ptilde,am1)
!!$
!!$    if (debug_level >= psb_debug_outer_) &
!!$         & write(debug_unit,*) me,' ',trim(name),&
!!$         & 'Done NUMBMM 1'
!!$  else
!!$    !
!!$    ! Build the smoothed prolongator using the original matrix
!!$    !
!!$    if (psb_toupper(am3%fida) == 'CSR') then 
!!$      do i=1,am3%m
!!$        do j=am3%ia2(i),am3%ia2(i+1)-1
!!$          if (am3%ia1(j) == i) then 
!!$            am3%aspk(j) = done - omf(i)*am3%aspk(j) 
!!$          else
!!$            am3%aspk(j) = - omf(i)*am3%aspk(j) 
!!$          end if
!!$        end do
!!$      end do
!!$    else 
!!$      call psb_errpush(psb_err_internal_error_,name,a_err='Invalid AM3 storage format')
!!$      goto 9999
!!$    end if
!!$
!!$    if (debug_level >= psb_debug_outer_) &
!!$         & write(debug_unit,*) me,' ',trim(name),&
!!$         & 'Done gather, going for SYMBMM 1'
!!$    !
!!$    ! Symbmm90 does the allocation for its result.
!!$    ! 
!!$    ! am1 = (I-w*D*A)Ptilde
!!$    ! 
!!$    !
!!$    call psb_symbmm(am3,ptilde,am1,info)
!!$    if(info /= psb_success_) then
!!$      call psb_errpush(psb_err_from_subroutine_,name,a_err='symbmm 1')
!!$      goto 9999
!!$    end if
!!$
!!$    call psb_numbmm(am3,ptilde,am1)
!!$
!!$    if (debug_level >= psb_debug_outer_) &
!!$         & write(debug_unit,*) me,' ',trim(name),&
!!$         & 'Done NUMBMM 1'
!!$
!!$  end if
!!$
!!$  !
!!$  ! Ok, let's start over with the restrictor
!!$  ! 
!!$  if (.false.) then 
!!$    i = 4
!!$    select case (i)
!!$    case(1)
!!$
!!$      call psb_transp(ptilde,rtilde,fmt='CSR')
!!$      call psb_spcnv(a,atmp,info,afmt='CSR')
!!$
!!$      am4%fida='COO'
!!$      am4%m=ncol-nrow
!!$      am4%k=ncol
!!$      call psb_sp_all(ncol,ntaggr,am4,ncol,info)
!!$
!!$      do i=1,ncol-nrow
!!$        am4%aspk(i) = dzero
!!$        am4%ia1(i)  = i
!!$        am4%ia2(i)  = nrow+i
!!$      end do
!!$      call psb_sp_setifld(nrow-ncol,psb_nnz_,am4,info)
!!$      call psb_spcnv(am4,info,afmt='CSR')
!!$      if (info == psb_success_) call psb_rwextd(ncol,atmp,info,b=am4)      
!!$      if (info == psb_success_) call psb_sp_free(am4,info)
!!$
!!$    case(2)
!!$
!!$      call psb_transp(ptilde,rtilde,fmt='CSR')
!!$      call psb_spcnv(a,atmp,info,afmt='CSR') 
!!$      call psb_sphalo(atmp,desc_a,am4,info,&
!!$           & colcnv=.true.,rowscale=.true.)
!!$      nrt  = psb_sp_get_nrows(am4) 
!!$      call psb_sp_clip(am4,atmp2,info,1,nrt,1,ncol)
!!$      call psb_spcnv(atmp2,info,afmt='CSR')
!!$      atmp2%aspk(:) = dzero
!!$      if (info == psb_success_) call psb_rwextd(ncol,atmp,info,b=atmp2)      
!!$      if (info == psb_success_) call psb_sp_free(am4,info)
!!$      if (info == psb_success_) call psb_sp_free(atmp2,info)
!!$
!!$    case (3)
!!$
!!$      ! We are doing the product only on the local
!!$      ! rows, the non-local contributions will be handled
!!$      ! through the global sum. 
!!$      call psb_transp(ptilde,am4,fmt='CSR')
!!$      nrt  = psb_sp_get_nrows(am4) 
!!$      call psb_sp_clip(am4,rtilde,info,1,nrt,1,nrow)
!!$      call psb_spcnv(a,atmp,info,afmt='CSR') 
!!$
!!$    case(4)
!!$
!!$      call psb_transp(ptilde,rtilde,fmt='COO')
!!$      do i=1, psb_sp_get_nnzeros(rtilde)
!!$        if (rtilde%ia2(i) > nrow) then 
!!$          rtilde%aspk(i) = dzero
!!$        end if
!!$      end do
!!$      call psb_spcnv(rtilde,info,afmt='CSR') 
!!$      call psb_spcnv(a,atmp,info,afmt='CSR') 
!!$      call psb_sphalo(atmp,desc_a,am4,info,&
!!$           & colcnv=.true.,rowscale=.true.)
!!$      nrt  = psb_sp_get_nrows(am4) 
!!$      call psb_sp_clip(am4,atmp2,info,1,nrt,1,ncol)
!!$      call psb_spcnv(atmp2,info,afmt='CSR')
!!$! !$    atmp2%aspk(:) = dzero
!!$      if (info == psb_success_) call psb_rwextd(ncol,atmp,info,b=atmp2)      
!!$      if (info == psb_success_) call psb_sp_free(am4,info)
!!$      if (info == psb_success_) call psb_sp_free(atmp2,info)
!!$
!!$    case default
!!$      write(0,*) 'Not building rtilde/atmp, this will blow up'
!!$      info = psb_err_from_subroutine_
!!$      goto 9999 
!!$    end select
!!$
!!$    if (info == psb_success_) call psb_symbmm(rtilde,atmp,ra,info)
!!$    if (info == psb_success_) call psb_numbmm(rtilde,atmp,ra)
!!$    if (info /= psb_success_) then 
!!$      write(0,*) 'From symbmm 1:',info
!!$      goto 9999
!!$    end if
!!$    call psb_sp_scal(adinv,atmp,info)
!!$    if (info == psb_success_) call psb_symbmm(ra,atmp,rada,info)
!!$    if (info == psb_success_) call psb_numbmm(ra,atmp,rada)
!!$    if (info /= psb_success_) then 
!!$      write(0,*) 'From symbmm 2:',info
!!$      goto 9999
!!$    end if
!!$
!!$    call csr_mat_row_prod(ra,rada,omp,info)
!!$    call csr_mat_row_prod(rada,rada,oden,info)
!!$    call psb_sum(ictxt,omp)
!!$    call psb_sum(ictxt,oden)
!!$  else
!!$
!!$    call psb_transp(ptilde,rtilde,fmt='CSR')
!!$    call psb_spcnv(a,atmp,info,afmt='CSR') 
!!$    call psb_sphalo(atmp,desc_a,am4,info,&
!!$         & colcnv=.true.,rowscale=.true.)
!!$    nrt  = psb_sp_get_nrows(am4) 
!!$    call psb_sp_clip(am4,atmp2,info,1,nrt,1,ncol)
!!$    call psb_spcnv(atmp2,info,afmt='CSR')
!!$    if (info == psb_success_) call psb_rwextd(ncol,atmp,info,b=atmp2)      
!!$    if (info == psb_success_) call psb_sp_free(am4,info)
!!$    if (info == psb_success_) call psb_sp_free(atmp2,info)
!!$    ! This is to compute the transpose. It ONLY works if the
!!$    ! original A has a symmetric pattern.
!!$    call psb_transp(atmp,atmp2) 
!!$    call psb_sp_clip(atmp2,atran,info,1,nrow,1,ncol)
!!$    call psb_sp_free(atmp2,info) 
!!$    ! Now for the product. 
!!$    call psb_symbmm(atran,ptilde,atp,info)
!!$    if (info == psb_success_) call psb_numbmm(atran,ptilde,atp)
!!$    call psb_sp_clone(atp,atmp2,info)
!!$    call psb_sp_scal(adinv,atmp2,info)
!!$    call psb_sphalo(atmp2,desc_a,am4,info,&
!!$         & colcnv=.false.,rowscale=.true.,outfmt='CSR  ')
!!$    if (info == psb_success_) call psb_rwextd(ncol,atmp2,info,b=am4)      
!!$    if (info == psb_success_) call psb_sp_free(am4,info)
!!$    
!!$    call psb_symbmm(atran,atmp2,atdatp,info)
!!$    call psb_numbmm(atran,atmp2,atdatp)
!!$    call psb_sp_free(atmp2,info)
!!$    
!!$    call psb_spcnv(atp,info,afmt='coo')
!!$    if (info == psb_success_) call psb_spcnv(atp,info,afmt='csc')
!!$    if (info == psb_success_) call psb_spcnv(atdatp,info,afmt='coo')
!!$    if (info == psb_success_) call psb_spcnv(atdatp,info,afmt='csc')
!!$    if (info /= psb_success_) then 
!!$      write(0,*) 'Failed conversion to CSC'
!!$    end if
!!$    
!!$    call csc_mat_col_prod(atp,atdatp,omp,info)
!!$    call csc_mat_col_prod(atdatp,atdatp,oden,info)
!!$    call psb_sum(ictxt,omp)
!!$    call psb_sum(ictxt,oden)
!!$    
!!$ 
!!$  end if
!!$! !$  write(debug_unit,*) trim(name),' OMP_R :',omp
!!$! ! $  write(debug_unit,*) trim(name),' ODEN_R:',oden
!!$  omp = omp/oden
!!$! !$  write(0,*) 'Check on output restrictor',omp(1:min(size(omp),10))
!!$  ! Compute omega_int
!!$  ommx = -1d300
!!$  do i=1, ncol
!!$    omi(i) = omp(ilaggr(i))
!!$    ommx = max(ommx,omi(i))
!!$  end do
!!$  ! Compute omega_fine
!!$  ! Going over the columns of atmp means going over the rows
!!$  ! of A^T. Hopefully ;-) 
!!$  call psb_spcnv(atmp,atmp2,info,afmt='coo')
!!$  if (info == psb_success_) call psb_spcnv(atmp2,info,afmt='csc')
!!$
!!$  do i=1, nrow
!!$    omf(i) = ommx
!!$    do j=atmp2%ia2(i),atmp2%ia2(i+1)-1
!!$      omf(i) = min(omf(i),omi(atmp2%ia1(j)))
!!$    end do
!!$    omf(i) = max(dzero,omf(i))
!!$  end do
!!$  omf(1:nrow) = omf(1:nrow)*adinv(1:nrow)
!!$  call psb_halo(omf,desc_a,info)
!!$  call psb_sp_free(atmp2,info) 
!!$
!!$
!!$  if (psb_toupper(atmp%fida) == 'CSR') then 
!!$    do i=1,atmp%m
!!$      do j=atmp%ia2(i),atmp%ia2(i+1)-1
!!$        if (atmp%ia1(j) == i) then 
!!$          atmp%aspk(j) = done - atmp%aspk(j)*omf(atmp%ia1(j))
!!$        else
!!$          atmp%aspk(j) =      - atmp%aspk(j)*omf(atmp%ia1(j))
!!$        end if
!!$      end do
!!$    end do
!!$  else
!!$    call psb_errpush(psb_err_internal_error_,name,a_err='Invalid ATMP storage format')
!!$    goto 9999
!!$  end if
!!$  call psb_symbmm(rtilde,atmp,am2,info)
!!$  call psb_numbmm(rtilde,atmp,am2)
!!$
!!$  !
!!$  ! Now we have to gather the halo of am1, and add it to itself
!!$  ! to multiply it by A,
!!$  !
!!$  call psb_sphalo(am1,desc_a,am4,info,&
!!$       & colcnv=.false.,rowscale=.true.)
!!$  if (info == psb_success_) call psb_rwextd(ncol,am1,info,b=am4)      
!!$  if (info == psb_success_) call psb_sp_free(am4,info)
!!$
!!$  if(info /= psb_success_) then
!!$    call psb_errpush(psb_err_internal_error_,name,a_err='Halo of am1')
!!$    goto 9999
!!$  end if
!!$
!!$
!!$
!!$  call psb_symbmm(a,am1,am3,info)
!!$  if(info /= psb_success_) then
!!$    call psb_errpush(psb_err_from_subroutine_,name,a_err='symbmm 2')
!!$    goto 9999
!!$  end if
!!$
!!$  call psb_numbmm(a,am1,am3)
!!$  if (debug_level >= psb_debug_outer_) &
!!$       & write(debug_unit,*) me,' ',trim(name),&
!!$       & 'Done NUMBMM 2'
!!$
!!$  !
!!$  ! Now we have to fix this.  The only rows of B that are correct 
!!$  ! are those corresponding to "local" aggregates, i.e. indices in ilaggr(:)
!!$  !
!!$  call psb_spcnv(am2,info,afmt='COO')
!!$  nzl = psb_sp_get_nnzeros(am2)
!!$  i=0
!!$  do k=1, nzl
!!$    if ((naggrm1 < am2%ia1(k)) .and. (am2%ia1(k) <= naggrp1)) then
!!$      i = i+1
!!$      am2%aspk(i) = am2%aspk(k)
!!$      am2%ia1(i)  = am2%ia1(k)
!!$      am2%ia2(i)  = am2%ia2(k)
!!$    end if
!!$  end do
!!$  am2%infoa(psb_nnz_) = i
!!$  call psb_spcnv(am2,info,afmt='csr',dupl=psb_dupl_add_)
!!$  if (info /= psb_success_) then 
!!$    call psb_errpush(psb_err_from_subroutine_,name,a_err='spcnv am2')
!!$    goto 9999
!!$  end if
!!$
!!$  if (debug_level >= psb_debug_outer_) &
!!$       & write(debug_unit,*) me,' ',trim(name),&
!!$       & 'starting sphalo/ rwxtd'
!!$
!!$  ! am2 = ((i-wDA)Ptilde)^T
!!$  call psb_sphalo(am3,desc_a,am4,info,&
!!$       & colcnv=.false.,rowscale=.true.)
!!$  if (info == psb_success_) call psb_rwextd(ncol,am3,info,b=am4)      
!!$  if (info == psb_success_) call psb_sp_free(am4,info)
!!$
!!$  if(info /= psb_success_) then
!!$    call psb_errpush(psb_err_internal_error_,name,a_err='Extend am3')
!!$    goto 9999
!!$  end if
!!$
!!$
!!$  if (debug_level >= psb_debug_outer_) &
!!$       & write(debug_unit,*) me,' ',trim(name),&
!!$       & 'starting symbmm 3'
!!$  call psb_symbmm(am2,am3,b,info)
!!$  if (info == psb_success_) call psb_numbmm(am2,am3,b)
!!$  if (info == psb_success_) call psb_sp_free(am3,info)
!!$  if (info == psb_success_) call psb_spcnv(b,info,afmt='coo',dupl=psb_dupl_add_)
!!$  if (info /= psb_success_) then
!!$    call psb_errpush(psb_err_internal_error_,name,a_err='Build b = am2 x am3')
!!$    goto 9999
!!$  end if
!!$
!!$
!!$
!!$  select case(p%parms%coarse_mat)
!!$
!!$  case(mld_distr_mat_) 
!!$
!!$    call psb_sp_clone(b,p%ac,info)
!!$    nzac = p%ac%infoa(psb_nnz_) 
!!$    nzl =  p%ac%infoa(psb_nnz_) 
!!$    if (info == psb_success_) call psb_cdall(ictxt,p%desc_ac,info,nl=nlaggr(me+1))
!!$    if (info == psb_success_) call psb_cdins(nzl,p%ac%ia1,p%ac%ia2,p%desc_ac,info)
!!$    if (info == psb_success_) call psb_cdasb(p%desc_ac,info)
!!$    if (info == psb_success_) call psb_glob_to_loc(p%ac%ia1(1:nzl),p%desc_ac,info,iact='I')
!!$    if (info == psb_success_) call psb_glob_to_loc(p%ac%ia2(1:nzl),p%desc_ac,info,iact='I')
!!$    if (info /= psb_success_) then
!!$      call psb_errpush(psb_err_internal_error_,name,a_err='Creating p%desc_ac and converting ac')
!!$      goto 9999
!!$    end if
!!$    if (debug_level >= psb_debug_outer_) &
!!$         & write(debug_unit,*) me,' ',trim(name),&
!!$         & 'Assembld aux descr. distr.'
!!$
!!$
!!$    p%ac%m=psb_cd_get_local_rows(p%desc_ac)
!!$    p%ac%k=psb_cd_get_local_cols(p%desc_ac)
!!$    p%ac%fida='COO'
!!$    p%ac%descra='GUN'
!!$
!!$    call psb_sp_free(b,info)
!!$    if (info == psb_success_) deallocate(nzbr,idisp,stat=info)
!!$    if (info /= psb_success_) then
!!$      call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_sp_free')
!!$      goto 9999
!!$    end if
!!$
!!$    if (np>1) then 
!!$      nzl = psb_sp_get_nnzeros(am1)
!!$      call psb_glob_to_loc(am1%ia1(1:nzl),p%desc_ac,info,'I')
!!$      if(info /= psb_success_) then
!!$        call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_glob_to_loc')
!!$        goto 9999
!!$      end if
!!$    endif
!!$    am1%k=psb_cd_get_local_cols(p%desc_ac)
!!$
!!$    if (np>1) then 
!!$      call psb_spcnv(am2,info,afmt='coo',dupl=psb_dupl_add_)
!!$      nzl = am2%infoa(psb_nnz_) 
!!$      if (info == psb_success_) call psb_glob_to_loc(am2%ia1(1:nzl),p%desc_ac,info,'I')
!!$      if (info == psb_success_) call psb_spcnv(am2,info,afmt='csr',dupl=psb_dupl_add_)        
!!$      if(info /= psb_success_) then
!!$        call psb_errpush(psb_err_internal_error_,name,a_err='Converting am2 to local')
!!$        goto 9999
!!$      end if
!!$    end if
!!$    am2%m=psb_cd_get_local_cols(p%desc_ac)
!!$
!!$    if (debug_level >= psb_debug_outer_) &
!!$         & write(debug_unit,*) me,' ',trim(name),&
!!$         & 'Done ac '
!!$
!!$  case(mld_repl_mat_) 
!!$    !
!!$    !
!!$    call psb_cdall(ictxt,p%desc_ac,info,mg=ntaggr,repl=.true.)
!!$    nzbr(:) = 0
!!$    nzbr(me+1) = b%infoa(psb_nnz_)
!!$
!!$    call psb_sum(ictxt,nzbr(1:np))
!!$    nzac = sum(nzbr)
!!$    if (info == psb_success_) call psb_sp_all(ntaggr,ntaggr,p%ac,nzac,info)
!!$    if (info /= psb_success_) goto 9999
!!$
!!$    do ip=1,np
!!$      idisp(ip) = sum(nzbr(1:ip-1))
!!$    enddo
!!$    ndx = nzbr(me+1) 
!!$
!!$    call mpi_allgatherv(b%aspk,ndx,mpi_double_precision,p%ac%aspk,nzbr,idisp,&
!!$         & mpi_double_precision,icomm,info)
!!$    if (info == psb_success_) call mpi_allgatherv(b%ia1,ndx,mpi_integer,p%ac%ia1,nzbr,idisp,&
!!$         & mpi_integer,icomm,info)
!!$    if (info == psb_success_) call mpi_allgatherv(b%ia2,ndx,mpi_integer,p%ac%ia2,nzbr,idisp,&
!!$         & mpi_integer,icomm,info)
!!$
!!$    if (info /= psb_success_) then 
!!$      call psb_errpush(psb_err_internal_error_,name,a_err=' from mpi_allgatherv')
!!$      goto 9999
!!$    end if
!!$
!!$    p%ac%m = ntaggr
!!$    p%ac%k = ntaggr
!!$    p%ac%infoa(psb_nnz_) = nzac
!!$    p%ac%fida='COO'
!!$    p%ac%descra='GUN'
!!$    call psb_spcnv(p%ac,info,afmt='coo',dupl=psb_dupl_add_)
!!$    if(info /= psb_success_) goto 9999
!!$    call psb_sp_free(b,info)
!!$    if(info /= psb_success_) goto 9999
!!$
!!$    deallocate(nzbr,idisp,stat=info)
!!$    if (info /= psb_success_) then 
!!$      info = psb_err_alloc_dealloc_
!!$      call psb_errpush(info,name)
!!$      goto 9999
!!$    end if
!!$  case default 
!!$    info = psb_err_internal_error_
!!$    call psb_errpush(info,name,a_err='invalid mld_coarse_mat_')
!!$    goto 9999
!!$  end select
!!$
!!$
!!$
!!$  call psb_spcnv(p%ac,info,afmt='csr',dupl=psb_dupl_add_)
!!$  if(info /= psb_success_) then
!!$    call psb_errpush(psb_err_from_subroutine_,name,a_err='spcnv')
!!$    goto 9999
!!$  end if
!!$
!!$  !
!!$  ! Copy the prolongation/restriction matrices into the descriptor map.
!!$  !  am2 => R    i.e. restriction  operator
!!$  !  am1 => P    i.e. prolongation operator
!!$  !  
!!$  p%map = psb_linmap(psb_map_aggr_,desc_a,&
!!$       & p%desc_ac,am2,am1,ilaggr,nlaggr)
!!$  if (info == psb_success_) call psb_sp_free(am1,info)
!!$  if (info == psb_success_) call psb_sp_free(am2,info)
!!$  if(info /= psb_success_) then
!!$    call psb_errpush(psb_err_from_subroutine_,name,a_err='sp_Free')
!!$    goto 9999
!!$  end if


  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done smooth_aggregate '
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name)
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return


contains

!!$  subroutine csc_mat_col_prod(a,b,v,info)
!!$    type(psb_dspmat_type), intent(in) :: a, b 
!!$    real(psb_dpk_), intent(out)       :: v(:)
!!$    integer, intent(out)              :: info
!!$
!!$    integer                           :: i,j,k, nr, nc,iap,nra,ibp,nrb
!!$    logical                           :: csca, cscb
!!$
!!$    info = psb_success_
!!$    nc   = psb_sp_get_ncols(a)
!!$    if (nc /= psb_sp_get_ncols(b)) then 
!!$      write(0,*) 'Matrices A and B should have same columns'
!!$      info = -1
!!$      return
!!$    end if
!!$    csca = (psb_toupper(a%fida(1:3)) == 'CSC')
!!$    cscb = (psb_toupper(b%fida(1:3)) == 'CSC')
!!$
!!$    if (.not.(csca.and.cscb)) then 
!!$      write(0,*) 'Matrices A and B should be in CSC'
!!$      info = -2
!!$      return
!!$    end if
!!$
!!$    do j=1, nc
!!$      iap  = a%ia2(j)
!!$      nra  = a%ia2(j+1)-iap
!!$      ibp  = b%ia2(j)
!!$      nrb  = b%ia2(j+1)-ibp
!!$      v(j) = sparse_srtd_dot(nra,a%ia1(iap:iap+nra-1),a%aspk(iap:iap+nra-1),&
!!$           & nrb,b%ia1(ibp:ibp+nrb-1),b%aspk(ibp:ibp+nrb-1))
!!$    end do
!!$
!!$  end subroutine csc_mat_col_prod
!!$
!!$
!!$  subroutine csr_mat_row_prod(a,b,v,info)
!!$    type(psb_dspmat_type), intent(in) :: a, b 
!!$    real(psb_dpk_), intent(out)       :: v(:)
!!$    integer, intent(out)              :: info
!!$
!!$    integer                           :: i,j,k, nr, nc,iap,nca,ibp,ncb
!!$    logical                           :: csra, csrb
!!$
!!$    info = psb_success_
!!$    nr   = psb_sp_get_nrows(a)
!!$    if (nr /= psb_sp_get_nrows(b)) then 
!!$      write(0,*) 'Matrices A and B should have same rows'
!!$      info = -1
!!$      return
!!$    end if
!!$    csra = (psb_toupper(a%fida(1:3)) == 'CSR')
!!$    csrb = (psb_toupper(b%fida(1:3)) == 'CSR')
!!$
!!$    if (.not.(csra.and.csrb)) then 
!!$      write(0,*) 'Matrices A and B should be in CSR'
!!$      info = -2
!!$      return
!!$    end if
!!$
!!$    do j=1, nr
!!$      iap  = a%ia2(j)
!!$      nca  = a%ia2(j+1)-iap
!!$      ibp  = b%ia2(j)
!!$      ncb  = b%ia2(j+1)-ibp
!!$      v(j) = sparse_srtd_dot(nca,a%ia1(iap:iap+nca-1),a%aspk(iap:iap+nca-1),&
!!$           & ncb,b%ia1(ibp:ibp+ncb-1),b%aspk(ibp:ibp+ncb-1))
!!$    end do
!!$
!!$  end subroutine csr_mat_row_prod

  function sparse_srtd_dot(nv1,iv1,v1,nv2,iv2,v2) result(dot) 
    integer, intent(in) :: nv1,nv2
    integer, intent(in) :: iv1(:), iv2(:)
    real(psb_dpk_), intent(in) :: v1(:),v2(:)
    real(psb_dpk_)      :: dot

    integer :: i,j,k, ip1, ip2

    dot = dzero 
    ip1 = 1
    ip2 = 1

    do 
      if (ip1 > nv1) exit
      if (ip2 > nv2) exit
      if (iv1(ip1) == iv2(ip2)) then 
        dot = dot + v1(ip1)*v2(ip2)
        ip1 = ip1 + 1
        ip2 = ip2 + 1
      else if (iv1(ip1) < iv2(ip2)) then 
        ip1 = ip1 + 1 
      else
        ip2 = ip2 + 1 
      end if
    end do

  end function sparse_srtd_dot

end subroutine mld_daggrmat_minnrg_asb
