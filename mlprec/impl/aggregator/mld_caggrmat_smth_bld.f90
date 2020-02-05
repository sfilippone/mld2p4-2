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
! File: mld_caggrmat_smth_bld.F90
!
! Subroutine: mld_caggrmat_smth_bld
! Version:    complex
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
!  through mld_cprecinit and mld_zprecset.
!
!  The coarse-level matrix A_C is distributed among the parallel processes or
!  replicated on each of them, according to the value of p%parms%coarse_mat,
!  specified by the user through mld_cprecinit and mld_zprecset.
!  On output from this routine the entries of AC, op_prol, op_restr
!  are still in "global numbering" mode; this is fixed in the calling routine
!  aggregator%mat_bld.
!
!
! Arguments:
!    a          -  type(psb_cspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!    p          -  type(mld_c_onelev_type), input/output.
!                  The 'one-level' data structure that will contain the local
!                  part of the matrix to be built as well as the information
!                  concerning the prolongator and its transpose.
!    parms      -   type(mld_sml_parms), input
!                  Parameters controlling the choice of algorithm
!    ac         -  type(psb_cspmat_type), output
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
!    op_prol    -  type(psb_cspmat_type), input/output
!                  The tentative prolongator on input, the computed prolongator on output
!               
!    op_restr    -  type(psb_cspmat_type), output
!                  The restrictor operator; normally, it is the transpose of the prolongator. 
!               
!    info       -  integer, output.
!                  Error code.
!
subroutine mld_caggrmat_smth_bld(a,desc_a,ilaggr,nlaggr,parms,ac,op_prol,op_restr,info)
  use psb_base_mod
  use mld_base_prec_type
  use mld_c_inner_mod, mld_protect_name => mld_caggrmat_smth_bld
  use mld_c_base_aggregator_mod

  implicit none

  ! Arguments
  type(psb_cspmat_type), intent(in)      :: a
  type(psb_desc_type), intent(in)          :: desc_a
  integer(psb_lpk_), intent(inout)         :: ilaggr(:), nlaggr(:)
  type(mld_sml_parms), intent(inout)    :: parms 
  type(psb_lcspmat_type), intent(inout)  :: op_prol
  type(psb_lcspmat_type), intent(out)    :: ac,op_restr
  integer(psb_ipk_), intent(out)           :: info

  ! Local variables
  integer(psb_lpk_) :: nrow, nglob, ncol, ntaggr, ip, &
       & naggr, nzl,naggrm1,naggrp1, i, j, k, jd, icolF, nrw
  integer(psb_ipk_) :: inaggr, nzlp
  integer(psb_ipk_) :: ictxt, np, me
  character(len=20) :: name
  type(psb_desc_type) :: tmp_desc
  type(psb_lc_coo_sparse_mat) :: tmpcoo
  type(psb_lc_csr_sparse_mat) :: acsr1,  acsrf, csr_prol, acsr
  complex(psb_spk_), allocatable :: adiag(:)
  real(psb_spk_), allocatable :: arwsum(:)
  integer(psb_ipk_)  :: ierr(5)
  logical            :: filter_mat
  integer(psb_ipk_)            :: debug_level, debug_unit, err_act
  integer(psb_ipk_), parameter :: ncmax=16
  real(psb_spk_)     :: anorm, omega, tmp, dg, theta
  logical, parameter :: debug_new=.false.
  character(len=80) :: filename

  name='mld_aggrmat_smth_bld'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = desc_a%get_context()

  call psb_info(ictxt, me, np)

  nglob = desc_a%get_global_rows()
  nrow  = desc_a%get_local_rows()
  ncol  = desc_a%get_local_cols()

  theta = parms%aggr_thresh

  naggr  = nlaggr(me+1)
  ntaggr = sum(nlaggr)

  naggrm1 = sum(nlaggr(1:me))
  naggrp1 = sum(nlaggr(1:me+1))
  filter_mat = (parms%aggr_filter == mld_filter_mat_)

  !
  ! naggr: number of local aggregates
  ! nrow: local rows. 
  ! 

  ! Get the diagonal D
  adiag = a%get_diag(info)
  if (info == psb_success_) &
       & call psb_realloc(ncol,adiag,info)
  if (info == psb_success_) &
       & call psb_halo(adiag,desc_a,info)
  if (info == psb_success_) call a%cp_to(acsr)
  call op_prol%cp_to(tmpcoo)

  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='sp_getdiag')
    goto 9999
  end if

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ' Initial copies done.'
  
  call acsr%cp_to_fmt(acsrf,info)

  
  if (filter_mat) then
    !
    ! Build the filtered matrix Af from A
    ! 

    do i=1, nrow
      tmp = czero
      jd  = -1 
      do j=acsrf%irp(i),acsrf%irp(i+1)-1
        if (acsrf%ja(j) == i) jd = j 
        if (abs(acsrf%val(j)) < theta*sqrt(abs(adiag(i)*adiag(acsrf%ja(j))))) then
          tmp=tmp+acsrf%val(j)
          acsrf%val(j)=czero
        endif

      enddo
      if (jd == -1) then 
        write(0,*) 'Wrong input: we need the diagonal!!!!', i
      else
        acsrf%val(jd)=acsrf%val(jd)-tmp
      end if
    enddo
    ! Take out zeroed terms 
    call acsrf%clean_zeros(info)
  end if


  do i=1,size(adiag)
    if (adiag(i) /= czero) then
      adiag(i) = cone / adiag(i)
    else
      adiag(i) = cone
    end if
  end do

  if (parms%aggr_omega_alg == mld_eig_est_) then 

    if (parms%aggr_eig == mld_max_norm_) then 
      allocate(arwsum(nrow))
      call acsr%arwsum(arwsum)      
      anorm = maxval(abs(adiag(1:nrow)*arwsum(1:nrow)))
      call psb_amx(ictxt,anorm)
      omega = 4.d0/(3.d0*anorm)
      parms%aggr_omega_val = omega 

    else 
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='invalid mld_aggr_eig_')
      goto 9999
    end if

  else if (parms%aggr_omega_alg == mld_user_choice_) then 

    omega = parms%aggr_omega_val 

  else if (parms%aggr_omega_alg /= mld_user_choice_) then 
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='invalid mld_aggr_omega_alg_')
    goto 9999
  end if

  
  call acsrf%scal(adiag,info)
  if (info /= psb_success_) goto 9999

  inaggr = naggr
  call psb_cdall(ictxt,tmp_desc,info,nl=inaggr)
  nzlp = tmpcoo%get_nzeros()
  call tmp_desc%indxmap%g2lip_ins(tmpcoo%ja(1:nzlp),info) 
  call tmpcoo%set_ncols(tmp_desc%get_local_cols())
  call tmpcoo%mv_to_fmt(csr_prol,info)  
  !
  ! Build the smoothed prolongator using either A or Af
  !    acsr1 = (I-w*D*A) Prol      acsr1 = (I-w*D*Af) Prol 
  ! This is always done through the variable acsrf which
  ! is a bit less readable, butsaves space and one extra matrix copy
  ! 
  call omega_smooth(omega,acsrf)
  call psb_par_spspmm(acsrf,desc_a,csr_prol,acsr1,tmp_desc,info)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='spspmm 1')
    goto 9999
  end if
  
  
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done SPSPMM 1'

  nzl = acsr1%get_nzeros()
  call tmp_desc%l2gip(acsr1%ja(1:nzl),info)
  call acsr1%set_dupl(psb_dupl_add_)
  call acsr1%set_ncols(ntaggr)
  call op_prol%cp_from(acsr1)

  call mld_spmm_bld_inner(acsr,desc_a,ilaggr,nlaggr,parms,ac,&
       & op_prol,op_restr,info)


  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done smooth_aggregate '
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name)
  call psb_error_handler(err_act)
  return

contains
  
  subroutine omega_smooth(omega,acsr)
    implicit none 
    real(psb_spk_),intent(in) :: omega
    type(psb_lc_csr_sparse_mat), intent(inout) :: acsr
    !
    integer(psb_lpk_) :: i,j
    do i=1,acsr%get_nrows()
      do j=acsr%irp(i),acsr%irp(i+1)-1
        if (acsr%ja(j) == i) then 
          acsr%val(j) = cone - omega*acsr%val(j) 
        else
          acsr%val(j) = - omega*acsr%val(j) 
        end if
      end do
    end do
  end subroutine omega_smooth

end subroutine mld_caggrmat_smth_bld
