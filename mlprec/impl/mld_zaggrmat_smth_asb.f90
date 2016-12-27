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
! File: mld_zaggrmat_smth_asb.F90
!
! Subroutine: mld_zaggrmat_smth_asb
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
!  through mld_zprecinit and mld_zprecset.
!
!  The coarse-level matrix A_C is distributed among the parallel processes or
!  replicated on each of them, according to the value of p%parms%coarse_mat,
!  specified by the user through mld_zprecinit and mld_zprecset.
!  On output from this routine the entries of AC, op_prol, op_restr
!  are still in "global numbering" mode; this is fixed in the calling routine
!  mld_z_lev_aggrmat_asb.
!
!  For more details see
!    M. Brezina and P. Vanek, A black-box iterative solver based on a 
!    two-level Schwarz method, Computing,  63 (1999), 233-263.
!    P. D'Ambra, D. di Serafino and S. Filippone, On the development of
!    PSBLAS-based parallel two-level Schwarz preconditioners, Appl. Num. Math.
!    57 (2007), 1181-1196.
!
!
! Arguments:
!    a          -  type(psb_zspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!    p          -  type(mld_z_onelev_type), input/output.
!                  The 'one-level' data structure that will contain the local
!                  part of the matrix to be built as well as the information
!                  concerning the prolongator and its transpose.
!    parms      -   type(mld_dml_parms), input
!                  Parameters controlling the choice of algorithm
!    ac         -  type(psb_zspmat_type), output
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
!    op_prol    -  type(psb_zspmat_type), input/output
!                  The tentative prolongator on input, the computed prolongator on output
!               
!    op_restr    -  type(psb_zspmat_type), output
!                  The restrictor operator; normally, it is the transpose of the prolongator. 
!               
!    info       -  integer, output.
!                  Error code.
!
subroutine mld_zaggrmat_smth_asb(a,desc_a,ilaggr,nlaggr,parms,ac,op_prol,op_restr,info)
  use psb_base_mod
  use mld_z_inner_mod, mld_protect_name => mld_zaggrmat_smth_asb

  implicit none

  ! Arguments
  type(psb_zspmat_type), intent(in)        :: a
  type(psb_desc_type), intent(in)            :: desc_a
  integer(psb_ipk_), intent(inout)           :: ilaggr(:), nlaggr(:)
  type(mld_dml_parms), intent(inout)      :: parms 
  type(psb_zspmat_type), intent(inout)     :: op_prol
  type(psb_zspmat_type), intent(out)       :: ac,op_restr
  integer(psb_ipk_), intent(out)             :: info

  ! Local variables
  integer(psb_ipk_) :: nrow, nglob, ncol, ntaggr, ip, ndx,&
       & naggr, nzl,naggrm1,naggrp1, i, j, k, jd, icolF, nrw, err_act
  integer(psb_ipk_) ::ictxt, np, me
  character(len=20) :: name
  type(psb_zspmat_type) :: am3, am4, tmp_prol
  type(psb_z_coo_sparse_mat) :: tmpcoo
  type(psb_z_csr_sparse_mat) :: acsr1, acsr2, acsr3, acsrf, ptilde
  complex(psb_dpk_), allocatable :: adiag(:)
  integer(psb_ipk_)  :: ierr(5)
  logical            :: filter_mat
  integer(psb_ipk_)            :: debug_level, debug_unit
  integer(psb_ipk_), parameter :: ncmax=16
  real(psb_dpk_)     :: anorm, omega, tmp, dg, theta

  name='mld_aggrmat_smth_asb'
  if(psb_get_errstatus().ne.0) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = desc_a%get_context()
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

  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='sp_getdiag')
    goto 9999
  end if

  ! 1. Allocate Ptilde in sparse matrix form 
  call op_prol%mv_to(tmpcoo)
  call ptilde%mv_from_coo(tmpcoo,info)
  if (info == psb_success_) call a%cscnv(acsr3,info,dupl=psb_dupl_add_)
  if (info /= psb_success_) goto 9999

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ' Initial copies done.'

  if (filter_mat) then
    !
    ! Build the filtered matrix Af from A
    ! 
    if (info == psb_success_) call acsr3%cp_to_fmt(acsrf,info)

    do i=1,nrow
      tmp = zzero
      jd  = -1 
      do j=acsrf%irp(i),acsrf%irp(i+1)-1
        if (acsrf%ja(j) == i) jd = j 
        if (abs(acsrf%val(j)) < theta*sqrt(abs(adiag(i)*adiag(acsrf%ja(j))))) then
          tmp=tmp+acsrf%val(j)
          acsrf%val(j)=zzero
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
    if (adiag(i) /= zzero) then
      adiag(i) = zone / adiag(i)
    else
      adiag(i) = zone
    end if
  end do

  if (filter_mat) call acsrf%scal(adiag,info)
  if (info == psb_success_) call acsr3%scal(adiag,info)
  if (info /= psb_success_) goto 9999


  if (parms%aggr_omega_alg == mld_eig_est_) then 

    if (parms%aggr_eig == mld_max_norm_) then 

      anorm = acsr3%spnmi()
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

  if (filter_mat) then
    !
    ! Build the smoothed prolongator using the filtered matrix
    ! 
    do i=1,acsrf%get_nrows()
      do j=acsrf%irp(i),acsrf%irp(i+1)-1
        if (acsrf%ja(j) == i) then 
          acsrf%val(j) = zone - omega*acsrf%val(j) 
        else
          acsrf%val(j) = - omega*acsrf%val(j) 
        end if
      end do
    end do

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done gather, going for SPSPMM 1'
    !
    ! 
    ! acsrm1 = (I-w*D*Af)Ptilde
    ! Doing it this way means to consider diag(Af_i)
    ! 
    !
    call psb_spspmm(acsrf,ptilde,acsr1,info)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='spspmm 1')
      goto 9999
    end if

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done SPSPMM 1'

  else

    !
    ! Build the smoothed prolongator using the original matrix
    !
    do i=1,acsr3%get_nrows()
      do j=acsr3%irp(i),acsr3%irp(i+1)-1
        if (acsr3%ja(j) == i) then 
          acsr3%val(j) = zone - omega*acsr3%val(j) 
        else
          acsr3%val(j) = - omega*acsr3%val(j) 
        end if
      end do
    end do

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done gather, going for SPSPMM 1'
    !
    ! acsrm1 = (I-w*D*A)Ptilde
    ! Doing it this way means to consider diag(A_i)
    ! 
    !
    call psb_spspmm(acsr3,ptilde,acsr1,info)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='spspmm 1')
      goto 9999
    end if

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done SPSPMM 1'

  end if
  call ptilde%free()
  call acsr1%set_dupl(psb_dupl_add_)

  call op_prol%cp_from(acsr1)
  call tmp_prol%mv_from(acsr1)
  !
  ! Now we have to gather the halo of tmp_prol, and add it to itself
  ! to multiply it by A,
  !
  call psb_sphalo(tmp_prol,desc_a,am4,info,&
       & colcnv=.false.,rowscale=.true.)
  if (info == psb_success_) call psb_rwextd(ncol,tmp_prol,info,b=am4)      
  if (info == psb_success_) call am4%free()
  if(info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='Halo of tmp_prol')
    goto 9999
  end if

  call psb_spspmm(a,tmp_prol,am3,info)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='spspmm 2')
    goto 9999
  end if

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done SPSPMM 2',parms%aggr_kind, mld_smooth_prol_

  call tmp_prol%cp_to(tmpcoo)
  call tmpcoo%transp()

  nzl = tmpcoo%get_nzeros()
  i=0
  !
  ! Now we have to fix this.  The only rows of B that are correct 
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

  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='spcnv op_restr')
    goto 9999
  end if
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'starting sphalo/ rwxtd'

  ! op_restr = ((i-wDA)Ptilde)^T
  call psb_sphalo(am3,desc_a,am4,info,&
       & colcnv=.false.,rowscale=.true.)
  if (info == psb_success_) call psb_rwextd(ncol,am3,info,b=am4)      
  if (info == psb_success_) call am4%free()
  if(info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='Extend am3')
    goto 9999
  end if


  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'starting spspmm 3'
  call psb_spspmm(op_restr,am3,ac,info)
  if (info == psb_success_) call am3%free()
  if (info == psb_success_) call ac%cscnv(info,type='csr',dupl=psb_dupl_add_)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='Build ac = op_restr x am3')
    goto 9999
  end if

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done smooth_aggregate '
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name)
  call psb_error_handler(err_act)
  return

end subroutine mld_zaggrmat_smth_asb
