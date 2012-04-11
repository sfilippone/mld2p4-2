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
! File: mld_saggrmat_minnrg_asb.F90
!
! Subroutine: mld_saggrmat_minnrg_asb
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
!  through mld_sprecinit and mld_sprecset.
!
!  This routine can also build A_C according to a "bizarre" aggregation algorithm,
!  using a "naive" prolongator proposed by the authors of MLD2P4. However, this
!  prolongator still requires a deep analysis and testing and its use is not
!  recommended.
!
!  The coarse-level matrix A_C is distributed among the parallel processes or
!  replicated on each of them, according to the value of p%parms%coarse_mat,
!  specified by the user through mld_sprecinit and mld_zprecset.
!
!  For more details see
!    M. Brezina and P. Vanek, A black-box iterative solver based on a 
!    two-level Schwarz method, Computing,  63 (1999), 233-263.
!    P. D'Ambra, D. di Serafino and S. Filippone, On the development of
!    PSBLAS-based parallel two-level Schwarz preconditioners, Appl. Num. Math.
!    57 (2007), 1181-1196.
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
subroutine mld_saggrmat_minnrg_asb(a,desc_a,ilaggr,nlaggr,parms,ac,op_prol,op_restr,info)
  use psb_base_mod
  use mld_s_inner_mod, mld_protect_name => mld_saggrmat_minnrg_asb

  implicit none

  ! Arguments
  type(psb_sspmat_type), intent(in)             :: a
  type(psb_desc_type), intent(in)               :: desc_a
  integer, intent(inout)                        :: ilaggr(:), nlaggr(:)
  type(mld_sml_parms), intent(inout)             :: parms 
  type(psb_sspmat_type), intent(out)             :: ac,op_prol,op_restr
  integer, intent(out)                          :: info

  ! Local variables
  integer(psb_mpik_), allocatable       :: nzbr(:), idisp(:)
  integer :: nrow, nglob, ncol, ntaggr, nzac, ip, ndx,&
       & naggr, nzl,naggrm1,naggrp1, i, j, k, jd, icolF, nrt
  integer                    :: ictxt,np,me, err_act, icomm
  character(len=20)          :: name
  type(psb_sspmat_type)      :: af, ptilde, rtilde, atran, atp, atdatp
  type(psb_sspmat_type)      :: am3,am4, ap, adap,atmp,rada, ra, atmp2, dap, dadap, da
  type(psb_sspmat_type)      :: dat, datp, datdatp, atmp3
  type(psb_s_coo_sparse_mat) :: tmpcoo
  type(psb_s_csr_sparse_mat) :: acsr1, acsr2, acsr3, acsr, acsrf
  type(psb_s_csc_sparse_mat) :: csc_dap, csc_dadap, csc_datp, csc_datdatp, acsc
  real(psb_spk_), allocatable :: adiag(:), adinv(:)
  real(psb_spk_), allocatable :: omf(:), omp(:), omi(:), oden(:)
  logical                    :: filter_mat
  integer(psb_ipk_)          :: ierr(5)
  integer                    :: debug_level, debug_unit
  integer, parameter         :: ncmax=16
  real(psb_spk_)             :: anorm, theta
  real(psb_spk_)          :: tmp, alpha, beta, ommx

  name='mld_aggrmat_minnrg'
  if(psb_get_errstatus().ne.0) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = desc_a%get_context()
  icomm = desc_a%get_mpic()
  ictxt = desc_a%get_context()

  call psb_info(ictxt, me, np)

  nglob = desc_a%get_global_rows()
  nrow  = desc_a%get_local_rows()
  ncol  = desc_a%get_local_cols()

  theta = parms%aggr_thresh

  naggr  = nlaggr(me+1)
  ntaggr = sum(nlaggr)

  allocate(nzbr(np), idisp(np),stat=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_request_; ierr(1)=2*np;
    call psb_errpush(info,name,i_err=ierr,a_err='integer')
    goto 9999      
  end if

  naggrm1 = sum(nlaggr(1:me))
  naggrp1 = sum(nlaggr(1:me+1))

  filter_mat = (parms%aggr_filter == mld_filter_mat_)

  ilaggr(1:nrow) = ilaggr(1:nrow) + naggrm1
  call psb_halo(ilaggr,desc_a,info)

  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_halo')
    goto 9999
  end if

  ! naggr: number of local aggregates
  ! nrow: local rows. 
  ! 
  allocate(adiag(ncol),adinv(ncol),&
       & omf(ncol),omp(ntaggr),oden(ntaggr),omi(ncol),stat=info)

  if (info /= psb_success_) then 
    info=psb_err_alloc_request_; ierr(1)=6*ncol+ntaggr;
    call psb_errpush(info,name,i_err=ierr,a_err='real(psb_spk_)')
    goto 9999      
  end if

  ! Get the diagonal D
  call a%get_diag(adiag,info)
  if (info == psb_success_) &
       & call psb_halo(adiag,desc_a,info)

  do i=1,size(adiag)
    if (adiag(i) /= szero) then
      adinv(i) = sone / adiag(i)
    else
      adinv(i) = sone
    end if
  end do

  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='sp_getdiag')
    goto 9999
  end if


  ! 1. Allocate Ptilde in sparse matrix form 
  call tmpcoo%allocate(ncol,ntaggr,ncol)
  do i=1,ncol
    tmpcoo%val(i) = sone
    tmpcoo%ia(i)  = i
    tmpcoo%ja(i)  = ilaggr(i)  
  end do
  call tmpcoo%set_nzeros(ncol)
  call tmpcoo%set_dupl(psb_dupl_add_)
  call tmpcoo%set_asb()
  call ptilde%mv_from(tmpcoo)
  call ptilde%cscnv(info,type='csr')

!!$  call local_dump(me,ptilde,'csr-ptilde','Ptilde-1')

  if (info == psb_success_) call a%cscnv(am3,info,type='csr',dupl=psb_dupl_add_)
  if (info == psb_success_) call a%cscnv(da,info,type='csr',dupl=psb_dupl_add_)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='spcnv')
    goto 9999
  end if
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ' Initial copies done.'

  call da%scal(adinv,info)

  call psb_symbmm(da,ptilde,dap,info)
  if (info == psb_success_) call psb_numbmm(da,ptilde,dap)

  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='symbmm 1')
    goto 9999
  end if

  call dap%clone(atmp,info)

  call psb_sphalo(atmp,desc_a,am4,info,&
       & colcnv=.false.,rowscale=.true.,outfmt='CSR  ')
  if (info == psb_success_) call psb_rwextd(ncol,atmp,info,b=am4)      
  if (info == psb_success_) call am4%free()

  call psb_symbmm(da,atmp,dadap,info)
  call psb_numbmm(da,atmp,dadap)
  call atmp%free()

  !  !$  write(0,*) 'Columns of AP',psb_sp_get_ncols(ap)
  !  !$  write(0,*) 'Columns of ADAP',psb_sp_get_ncols(adap)
  call dap%mv_to(csc_dap)
  call dadap%mv_to(csc_dadap)


  call csc_mat_col_prod(csc_dap,csc_dadap,omp,info)
  call csc_mat_col_prod(csc_dadap,csc_dadap,oden,info)
  call psb_sum(ictxt,omp)
  call psb_sum(ictxt,oden)
  ! !$  write(0,*) trim(name),' OMP :',omp
  ! !$  write(0,*) trim(name),' ODEN:',oden

  omp = omp/oden

  ! !$  write(0,*) 'Check on output prolongator ',omp(1:min(size(omp),10))
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done NUMBMM 1'

  call am3%mv_to(acsr3)
  ! Compute omega_int
  ommx = cmplx(szero,szero)
  do i=1, ncol
    omi(i) = omp(ilaggr(i))
    if(abs(omi(i)) .gt. abs(ommx)) ommx = omi(i)
  end do
  ! Compute omega_fine
  do i=1, nrow
    omf(i) = ommx
    do j=acsr3%irp(i),acsr3%irp(i+1)-1
      if(abs(omi(acsr3%ja(j))) .lt. abs(omf(i))) omf(i)=omi(acsr3%ja(j))
    end do
!!$    if(min(real(omf(i)),aimag(omf(i))) < szero) omf(i) = szero
    if(psb_minreal(omf(i)) < szero) omf(i) = szero
  end do

  omf(1:nrow) = omf(1:nrow) * adinv(1:nrow)

  if (filter_mat) then
    !
    ! Build the filtered matrix Af from A
    ! 
    call a%cscnv(acsrf,info,dupl=psb_dupl_add_)

    do i=1,nrow
      tmp = szero
      jd  = -1 
      do j=acsrf%irp(i),acsrf%irp(i+1)-1
        if (acsrf%ja(j) == i) jd = j 
        if (abs(acsrf%val(j)) < theta*sqrt(abs(adiag(i)*adiag(acsrf%ja(j))))) then
          tmp=tmp+acsrf%val(j)
          acsrf%val(j)=szero
        endif
      enddo
      if (jd == -1) then 
        write(0,*) 'Wrong input: we need the diagonal!!!!', i
      else
        acsrf%val(jd)=acsrf%val(jd)-tmp
      end if
    enddo
    ! Take out zeroed terms 
    call acsrf%mv_to_coo(tmpcoo,info)
    k = 0
    do j=1,tmpcoo%get_nzeros() 
      if ((tmpcoo%val(j) /= szero) .or. (tmpcoo%ia(j) == tmpcoo%ja(j))) then 
        k = k + 1
        tmpcoo%val(k) = tmpcoo%val(j)
        tmpcoo%ia(k)  = tmpcoo%ia(j)
        tmpcoo%ja(k)  = tmpcoo%ja(j)
      end if
    end do
    call tmpcoo%set_nzeros(k)
    call acsrf%mv_from_coo(tmpcoo,info)

    !
    ! Build the smoothed prolongator using the filtered matrix
    ! 
    do i=1,acsrf%get_nrows()
      do j=acsrf%irp(i),acsrf%irp(i+1)-1
        if (acsrf%ja(j) == i) then 
          acsrf%val(j) = sone - omf(i)*acsrf%val(j) 
        else
          acsrf%val(j) = - omf(i)*acsrf%val(j) 
        end if
      end do
    end do

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done gather, going for SYMBMM 1'

    call af%mv_from(acsrf)
    !
    ! Symbmm90 does the allocation for its result.
    ! 
    ! op_prol = (I-w*D*Af)Ptilde
    ! Doing it this way means to consider diag(Af_i)
    ! 
    !
    call psb_symbmm(af,ptilde,op_prol,info)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='symbmm 1')
      goto 9999
    end if

    call psb_numbmm(af,ptilde,op_prol)

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done NUMBMM 1'
  else
    !
    ! Build the smoothed prolongator using the original matrix
    !
    do i=1,acsr3%get_nrows()
      do j=acsr3%irp(i),acsr3%irp(i+1)-1
        if (acsr3%ja(j) == i) then 
          acsr3%val(j) = sone - omf(i)*acsr3%val(j) 
        else
          acsr3%val(j) = - omf(i)*acsr3%val(j) 
        end if
      end do
    end do

    call am3%mv_from(acsr3)
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done gather, going for SYMBMM 1'
    !
    ! Symbmm90 does the allocation for its result.
    ! 
    ! op_prol = (I-w*D*A)Ptilde
    ! 
    !
    call psb_symbmm(am3,ptilde,op_prol,info)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='symbmm 1')
      goto 9999
    end if

    call psb_numbmm(am3,ptilde,op_prol)

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done NUMBMM 1'

  end if


  !
  ! Ok, let's start over with the restrictor
  ! 
  call ptilde%transc(rtilde)
  call a%cscnv(atmp,info,type='csr')
  call psb_sphalo(atmp,desc_a,am4,info,&
       & colcnv=.true.,rowscale=.true.)
  nrt  = am4%get_nrows() 
  call am4%csclip(atmp2,info,1,nrt,1,ncol)
  call atmp2%cscnv(info,type='CSR')
  if (info == psb_success_) call psb_rwextd(ncol,atmp,info,b=atmp2)      
  call am4%free()
  call atmp2%free()

  ! This is to compute the transpose. It ONLY works if the
  ! original A has a symmetric pattern.
  call atmp%transc(atmp2) 
  call atmp2%csclip(dat,info,1,nrow,1,ncol)
  call dat%cscnv(info,type='csr')
  call dat%scal(adinv,info)

  ! Now for the product. 
  call psb_symbmm(dat,ptilde,datp,info)
  if (info == psb_success_) call psb_numbmm(dat,ptilde,datp)

  call datp%clone(atmp2,info)
  call psb_sphalo(atmp2,desc_a,am4,info,&
       & colcnv=.false.,rowscale=.true.,outfmt='CSR  ')
  if (info == psb_success_) call psb_rwextd(ncol,atmp2,info,b=am4)      
  if (info == psb_success_) call am4%free()


  call psb_symbmm(dat,atmp2,datdatp,info)
  call psb_numbmm(dat,atmp2,datdatp)
  call atmp2%free()

  call datp%mv_to(csc_datp)    
  call datdatp%mv_to(csc_datdatp)    

  call csc_mat_col_prod(csc_datp,csc_datdatp,omp,info)
  call csc_mat_col_prod(csc_datdatp,csc_datdatp,oden,info)
  call psb_sum(ictxt,omp)
  call psb_sum(ictxt,oden)


  ! !$  write(debug_unit,*) trim(name),' OMP_R :',omp
  ! ! $  write(debug_unit,*) trim(name),' ODEN_R:',oden
  omp = omp/oden
  ! !$  write(0,*) 'Check on output restrictor',omp(1:min(size(omp),10))
  ! Compute omega_int
  ommx = cmplx(szero,szero)
  do i=1, ncol
    omi(i) = omp(ilaggr(i))
    if(abs(omi(i)) .gt. abs(ommx)) ommx = omi(i)
  end do
  ! Compute omega_fine
  ! Going over the columns of atmp means going over the rows
  ! of A^T. Hopefully ;-) 
  call atmp%cp_to(acsc)

  do i=1, nrow
    omf(i) = ommx
    do j= acsc%icp(i),acsc%icp(i+1)-1
      if(abs(omi(acsc%ia(j))) .lt. abs(omf(i))) omf(i)=omi(acsc%ia(j))
    end do
!!$    if(min(real(omf(i)),aimag(omf(i))) < szero) omf(i) = szero
    if(psb_minreal(omf(i)) < szero) omf(i) = szero
  end do
  omf(1:nrow) = omf(1:nrow)*adinv(1:nrow)
  call psb_halo(omf,desc_a,info)
  call acsc%free() 


  call atmp%mv_to(acsr1)

  do i=1,acsr1%get_nrows()
    do j=acsr1%irp(i),acsr1%irp(i+1)-1
      if (acsr1%ja(j) == i) then 
        acsr1%val(j) = sone - acsr1%val(j)*omf(acsr1%ja(j))
      else
        acsr1%val(j) =      - acsr1%val(j)*omf(acsr1%ja(j))
      end if
    end do
  end do
  call atmp%mv_from(acsr1)

  call rtilde%mv_to(tmpcoo)
  nzl = tmpcoo%get_nzeros()
  i=0
  do k=1, nzl
    if ((naggrm1 < tmpcoo%ia(k)) .and. (tmpcoo%ia(k) <= naggrp1)) then
      i = i+1
      tmpcoo%val(i) = tmpcoo%val(k)
      tmpcoo%ia(i)  = tmpcoo%ia(k)
      tmpcoo%ja(i)  = tmpcoo%ja(k)
    end if
  end do
  call tmpcoo%set_nzeros(i)
  call rtilde%mv_from(tmpcoo)
  call rtilde%cscnv(info,type='csr')

  call psb_symbmm(rtilde,atmp,op_restr,info)
  call psb_numbmm(rtilde,atmp,op_restr)

  !
  ! Now we have to gather the halo of op_prol, and add it to itself
  ! to multiply it by A,
  !
  call psb_sphalo(op_prol,desc_a,am4,info,&
       & colcnv=.false.,rowscale=.true.)
  if (info == psb_success_) call psb_rwextd(ncol,op_prol,info,b=am4)      
  if (info == psb_success_) call am4%free()

  if(info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='Halo of op_prol')
    goto 9999
  end if

  !
  ! Now we have to fix this.  The only rows of B that are correct 
  ! are those corresponding to "local" aggregates, i.e. indices in ilaggr(:)
  !
  call op_restr%mv_to(tmpcoo)

  nzl = tmpcoo%get_nzeros()
  i=0
  do k=1, nzl
    if ((naggrm1 < tmpcoo%ia(k)) .and. (tmpcoo%ia(k) <= naggrp1)) then
      i = i+1
      tmpcoo%val(i) = tmpcoo%val(k)
      tmpcoo%ia(i)  = tmpcoo%ia(k)
      tmpcoo%ja(i)  = tmpcoo%ja(k)
    end if
  end do
  call tmpcoo%set_nzeros(i)
  call op_restr%mv_from(tmpcoo)
  call op_restr%cscnv(info,type='csr')


  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'starting sphalo/ rwxtd'

  call psb_symbmm(a,op_prol,am3,info)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,&
         & a_err='symbmm 2')
    goto 9999
  end if
  call psb_numbmm(a,op_prol,am3)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done NUMBMM 2'

  call psb_sphalo(am3,desc_a,am4,info,&
       & colcnv=.false.,rowscale=.true.)
  if (info == psb_success_) call psb_rwextd(ncol,am3,info,b=am4)      
  if (info == psb_success_) call am4%free()

  if(info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,&
         & a_err='Extend am3')
    goto 9999
  end if
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done sphalo/ rwxtd'

  call psb_symbmm(op_restr,am3,ac,info)
  if (info == psb_success_) call psb_numbmm(op_restr,am3,ac)
  if (info == psb_success_) call am3%free()
  if (info == psb_success_) call ac%cscnv(info,type='coo',dupl=psb_dupl_add_)

  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,&
         &a_err='Build ac = op_restr x am3')
    goto 9999
  end if



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

  subroutine csc_mat_col_prod(a,b,v,info)
    type(psb_s_csc_sparse_mat), intent(in) :: a, b 
    real(psb_spk_), intent(out)       :: v(:)
    integer, intent(out)              :: info

    integer                           :: i,j,k, nr, nc,iap,nra,ibp,nrb

    info = psb_success_
    nc   = a%get_ncols()
    if (nc /= b%get_ncols()) then 
      write(0,*) 'Matrices A and B should have same columns'
      info = -1
      return
    end if

    do j=1, nc
      iap  = a%icp(j)
      nra  = a%icp(j+1)-iap
      ibp  = b%icp(j)
      nrb  = b%icp(j+1)-ibp
      v(j) = sparse_srtd_dot(nra,a%ia(iap:iap+nra-1),a%val(iap:iap+nra-1),&
           & nrb,b%ia(ibp:ibp+nrb-1),b%val(ibp:ibp+nrb-1))
    end do

  end subroutine csc_mat_col_prod


  subroutine csr_mat_row_prod(a,b,v,info)
    type(psb_s_csr_sparse_mat), intent(in) :: a, b 
    real(psb_spk_), intent(out)       :: v(:)
    integer, intent(out)              :: info

    integer                           :: i,j,k, nr, nc,iap,nca,ibp,ncb

    info = psb_success_
    nr   = a%get_nrows()
    if (nr /= b%get_nrows()) then 
      write(0,*) 'Matrices A and B should have same rows'
      info = -1
      return
    end if

    do j=1, nr
      iap  = a%irp(j)
      nca  = a%irp(j+1)-iap
      ibp  = b%irp(j)
      ncb  = b%irp(j+1)-ibp
      v(j) = sparse_srtd_dot(nca,a%ja(iap:iap+nca-1),a%val(iap:iap+nca-1),&
           & ncb,b%ja(ibp:ibp+ncb-1),b%val(ibp:ibp+ncb-1))
    end do

  end subroutine csr_mat_row_prod


  function sparse_srtd_dot(nv1,iv1,v1,nv2,iv2,v2) result(dot) 
    integer, intent(in) :: nv1,nv2
    integer, intent(in) :: iv1(:), iv2(:)
    real(psb_spk_), intent(in) :: v1(:),v2(:)
    real(psb_spk_)      :: dot

    integer :: i,j,k, ip1, ip2

    dot = szero 
    ip1 = 1
    ip2 = 1

    do 
      if (ip1 > nv1) exit
      if (ip2 > nv2) exit
      if (iv1(ip1) == iv2(ip2)) then 
        dot = dot + (v1(ip1))*v2(ip2)
        ip1 = ip1 + 1
        ip2 = ip2 + 1
      else if (iv1(ip1) < iv2(ip2)) then 
        ip1 = ip1 + 1 
      else
        ip2 = ip2 + 1 
      end if
    end do

  end function sparse_srtd_dot

  subroutine local_dump(me,mat,name,header)
    type(psb_sspmat_type), intent(in) :: mat
    integer, intent(in)               :: me
    character(len=*), intent(in)      :: name
    character(len=*), intent(in)      :: header
    character(len=80) :: filename

    write(filename,'(a,a,i0,a,i0,a)') trim(name),'.p',me
    open(20+me,file=filename)
    call mat%print(20+me,head=trim(header))
    close(20+me)
  end subroutine local_dump

end subroutine mld_saggrmat_minnrg_asb
