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
! File: mld_daggrmat_biz_asb.F90
!
! Subroutine: mld_daggrmat_biz_asb
! Version:    real
!
!  This routine builds a coarse-level matrix A_C from a fine-level matrix A
!  by using the Galerkin approach, i.e.
!
!                               A_C = P_C^T A P_C,
!
!  where P_C is a prolongator from the coarse level to the fine one.
! 
!  This routine builds A_C according to a "bizarre" aggregation algorithm,
!  using a "naive" prolongator proposed by the authors of MLD2P4. However, this
!  prolongator still requires a deep analysis and testing and its use is not
!  recommended.
!
!  The coarse-level matrix A_C is distributed among the parallel processes or
!  replicated on each of them, according to the value of p%parms%coarse_mat,
!  specified by the user through mld_dprecinit and mld_zprecset.
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
subroutine mld_daggrmat_biz_asb(a,desc_a,ilaggr,nlaggr,parms,ac,op_prol,op_restr,info)
  use psb_base_mod
  use mld_d_inner_mod, mld_protect_name => mld_daggrmat_biz_asb

  implicit none

  ! Arguments
  type(psb_dspmat_type), intent(in)             :: a
  type(psb_desc_type), intent(in)               :: desc_a
  integer(psb_ipk_), intent(inout)                        :: ilaggr(:), nlaggr(:)
  type(mld_dml_parms), intent(inout)             :: parms 
  type(psb_dspmat_type), intent(inout)             :: op_prol
  type(psb_dspmat_type), intent(out)             :: ac,op_restr
  integer(psb_ipk_), intent(out)                          :: info

  ! Local variables
  integer(psb_ipk_) :: nrow, nglob, ncol, ntaggr, ip, ndx,&
       & naggr, nzl,naggrm1,naggrp1, i, j, k, jd, icolF, nrw, err_act
  integer(psb_ipk_) ::ictxt, np, me
  character(len=20) :: name
  type(psb_dspmat_type) :: am3, am4
  type(psb_d_coo_sparse_mat) :: tmpcoo
  type(psb_d_csr_sparse_mat) :: acsr1, acsr2, acsr3, acsrf, ptilde
  real(psb_dpk_), allocatable :: adiag(:)
  integer(psb_ipk_)  :: ierr(5)
  logical            :: filter_mat
  integer(psb_ipk_)            :: debug_level, debug_unit
  integer(psb_ipk_), parameter :: ncmax=16
  real(psb_dpk_)     :: anorm, omega, tmp, dg, theta

  name='mld_aggrmat_biz_asb'
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

  filter_mat = (parms%aggr_filter == mld_filter_mat_)

  ! naggr: number of local aggregates
  ! nrow: local rows. 
  ! 
  ! Get the diagonal D
  adiag =  a%get_diag(info)
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
       & ' Initial copies sone.'

  if (filter_mat) then
    !
    ! Build the filtered matrix Af from A
    ! 
    if (info == psb_success_) call a%cscnv(acsrf,info,dupl=psb_dupl_add_)

    do i=1,nrow
      tmp = dzero
      jd  = -1 
      do j=acsrf%irp(i),acsrf%irp(i+1)-1
        if (acsrf%ja(j) == i) jd = j 
        if (abs(acsrf%val(j)) < theta*sqrt(abs(adiag(i)*adiag(acsrf%ja(j))))) then
          tmp=tmp+acsrf%val(j)
          acsrf%val(j)=dzero
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
      if ((tmpcoo%val(j) /= dzero) .or. (tmpcoo%ia(j) == tmpcoo%ja(j))) then 
        k = k + 1
        tmpcoo%val(k) = tmpcoo%val(j)
        tmpcoo%ia(k)  = tmpcoo%ia(j)
        tmpcoo%ja(k)  = tmpcoo%ja(j)
      end if
    end do
    call tmpcoo%set_nzeros(k)
    call tmpcoo%set_dupl(psb_dupl_add_)
    call acsrf%mv_from_coo(tmpcoo,info)
  end if


  do i=1,size(adiag)
    if (adiag(i) /= dzero) then
      adiag(i) = done / adiag(i)
    else
      adiag(i) = done
    end if
  end do

  if (filter_mat) call acsrf%scal(adiag,info)
  if (info == psb_success_) call acsr3%scal(adiag,info)
  if (info /= psb_success_) goto 9999


  if (parms%aggr_omega_alg == mld_eig_est_) then 

    if (parms%aggr_eig == mld_max_norm_) then 

      ! 
      ! This only works with CSR
      !
      anorm = dzero
      dg    = done
      nrw = acsr3%get_nrows()
      do i=1, nrw
        tmp = dzero
        do j=acsr3%irp(i),acsr3%irp(i+1)-1
          if (acsr3%ja(j) <= nrw) then 
            tmp = tmp + abs(acsr3%val(j))
          endif
          if (acsr3%ja(j) == i ) then 
            dg = abs(acsr3%val(j))
          end if
        end do
        anorm = max(anorm,tmp/dg) 
      enddo

      call psb_amx(ictxt,anorm)     
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_internal_error_,name,a_err='Invalid AM3 storage format')
        goto 9999
      end if
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
          acsrf%val(j) = done - omega*acsrf%val(j) 
        else
          acsrf%val(j) = - omega*acsrf%val(j) 
        end if
      end do
    end do

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done gather, going for SYMBMM 1'
    !
    ! Symbmm90 does the allocation for its result.
    ! 
    ! acsrm1 = (I-w*D*Af)Ptilde
    ! Doing it this way means to consider diag(Af_i)
    ! 
    !
    call psb_symbmm(acsrf,ptilde,acsr1,info)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='symbmm 1')
      goto 9999
    end if

    call psb_numbmm(acsrf,ptilde,acsr1)

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
          acsr3%val(j) = done - omega*acsr3%val(j) 
        else
          acsr3%val(j) = - omega*acsr3%val(j) 
        end if
      end do
    end do

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done gather, going for SYMBMM 1'
    !
    ! Symbmm90 does the allocation for its result.
    ! 
    ! acsrm1 = (I-w*D*A)Ptilde
    ! Doing it this way means to consider diag(A_i)
    ! 
    !
    call psb_symbmm(acsr3,ptilde,acsr1,info)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='symbmm 1')
      goto 9999
    end if

    call psb_numbmm(acsr3,ptilde,acsr1)

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done NUMBMM 1'

  end if
  call ptilde%free()
  call acsr1%set_dupl(psb_dupl_add_)

  call op_prol%mv_from(acsr1)

  call psb_rwextd(ncol,op_prol,info)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='Halo of op_prol')
    goto 9999
  end if

  call psb_symbmm(a,op_prol,am3,info)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='symbmm 2')
    goto 9999
  end if

  call psb_numbmm(a,op_prol,am3)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done NUMBMM 2',parms%aggr_kind, mld_smooth_prol_

  call op_prol%transp(op_restr)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'starting sphalo/ rwxtd'

  call psb_rwextd(ncol,am3,info)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='Extend am3')
    goto 9999
  end if


  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'starting symbmm 3'
  call psb_symbmm(op_restr,am3,ac,info)
  if (info == psb_success_) call psb_numbmm(op_restr,am3,ac)
  if (info == psb_success_) call am3%free()
  if (info == psb_success_) call ac%cscnv(info,type='coo',dupl=psb_dupl_add_)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='Build b = op_restr x am3')
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

end subroutine mld_daggrmat_biz_asb
