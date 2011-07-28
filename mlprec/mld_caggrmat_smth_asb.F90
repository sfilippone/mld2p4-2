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
! File: mld_caggrmat_smth_asb.F90
!
! Subroutine: mld_caggrmat_smth_asb
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
!  through mld_cprecinit and mld_cprecset.
!
!  This routine can also build A_C according to a "bizarre" aggregation algorithm,
!  using a "naive" prolongator proposed by the authors of MLD2P4. However, this
!  prolongator still requires a deep analysis and testing and its use is not
!  recommended.
!
!  The coarse-level matrix A_C is distributed among the parallel processes or
!  replicated on each of them, according to the value of p%parms%coarse_mat,
!  specified by the user through mld_cprecinit and mld_cprecset.
!
!  For more details see
!    M. Brezina and P. Vanek, A black-box iterative solver based on a 
!    two-level Schwarz method, Computing,  63 (1999), 233-263.
!    P. D'Ambra, D. di Serafino and S. Filippone, On the development of
!    PSBLAS-based parallel two-level Schwarz preconditioners, Appl. Num. Math.
!    57 (2007), 1181-1196.
!
! Arguments:
!    a          -  type(psb_cspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!    p          -  type(mld_conelev_type), input/output.
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
subroutine mld_caggrmat_smth_asb(a,desc_a,ilaggr,nlaggr,p,info)
  use psb_base_mod
  use mld_c_inner_mod, mld_protect_name => mld_caggrmat_smth_asb

#ifdef MPI_MOD
  use mpi
#endif
  implicit none 
#ifdef MPI_H
  include 'mpif.h'
#endif

  ! Arguments
  type(psb_cspmat_type), intent(in)              :: a
  type(psb_desc_type), intent(in)                :: desc_a
  integer, intent(inout)                          :: ilaggr(:), nlaggr(:)
  type(mld_conelev_type), intent(inout), target :: p
  integer, intent(out)                           :: info

  ! Local variables
  type(psb_cspmat_type) :: b
  integer, allocatable  :: nzbr(:), idisp(:)
  integer :: nrow, nglob, ncol, ntaggr, nzac, ip, ndx,&
       & naggr, nzl,naggrm1,naggrp1, i, j, k, jd, icolF, nrw
  integer ::ictxt,np,me, err_act, icomm
  character(len=20) :: name
  type(psb_cspmat_type) :: am1,am2, am3, am4
  type(psb_c_coo_sparse_mat) :: acoo1, acoo2, acoof, acoo3,acoo4, bcoo, cootmp
  type(psb_c_csr_sparse_mat) :: acsr1, acsr2, acsrf, acsr3,acsr4, bcsr
  complex(psb_spk_), allocatable :: adiag(:)
  logical            :: ml_global_nmb, filter_mat
  integer            :: debug_level, debug_unit
  integer, parameter :: ncmax=16
  real(psb_spk_)   :: omega, anorm, tmp, dg, theta

  name='mld_aggrmat_smth_asb'
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

  theta = p%parms%aggr_thresh

  naggr  = nlaggr(me+1)
  ntaggr = sum(nlaggr)

  allocate(nzbr(np), idisp(np),stat=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/2*np,0,0,0,0/),&
         & a_err='integer')
    goto 9999      
  end if

  naggrm1 = sum(nlaggr(1:me))
  naggrp1 = sum(nlaggr(1:me+1))
  ml_global_nmb = ( (p%parms%aggr_kind == mld_smooth_prol_).or.&
       & ( (p%parms%aggr_kind == mld_biz_prol_).and.&
       &    (p%parms%coarse_mat == mld_repl_mat_)) ) 

  filter_mat = (p%parms%aggr_filter == mld_filter_mat_)

  if (ml_global_nmb) then 
    ilaggr(1:nrow) = ilaggr(1:nrow) + naggrm1
    call psb_halo(ilaggr,desc_a,info)

    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_halo')
      goto 9999
    end if
  end if

  ! naggr: number of local aggregates
  ! nrow: local rows. 
  ! 
  allocate(adiag(ncol),stat=info)

  if (info /= psb_success_) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/nrow,0,0,0,0/),&
         & a_err='complex(psb_spk_)')
    goto 9999      
  end if

  ! Get the diagonal D
  call a%get_diag(adiag,info)
  if (info == psb_success_) &
       & call psb_halo(adiag,desc_a,info)

  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='sp_getdiag')
    goto 9999
  end if

  ! 1. Allocate Ptilde in sparse matrix form 
  if (ml_global_nmb) then 
    call acoo4%allocate(ncol,ntaggr,ncol)
    do i=1,ncol
      acoo4%val(i) = cone
      acoo4%ia(i)  = i
      acoo4%ja(i)  = ilaggr(i)  
    end do
    call acoo4%set_nzeros(ncol)
  else 
    call acoo4%allocate(ncol,naggr,ncol)
    do i=1,nrow
      acoo4%val(i) = cone
      acoo4%ia(i)  = i
      acoo4%ja(i)  = ilaggr(i)  
    end do
    call acoo4%set_nzeros(nrow)
  endif
  call acoo4%set_dupl(psb_dupl_add_)
  
  call acsr4%mv_from_coo(acoo4,info)
  if (info == psb_success_) call a%cscnv(acsr3,info,dupl=psb_dupl_add_)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ' Initial copies done.'
  
  if (filter_mat) then
    !
    ! Build the filtered matrix Af from A
    ! 
    if (info == psb_success_) call a%cscnv(acsrf,info,dupl=psb_dupl_add_)

    do i=1,nrow
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
    call acsrf%mv_to_coo(acoof,info)
    k = 0
    do j=1,acoof%get_nzeros()
      if ((acoof%val(j) /= czero) .or. (acoof%ia(j) == acoof%ja(j))) then 
        k = k + 1
        acoof%val(k) = acoof%val(j)
        acoof%ia(k)  = acoof%ia(j)
        acoof%ja(k)  = acoof%ja(j)
      end if
    end do
    call acoof%set_nzeros(k)
    call acoof%set_dupl(psb_dupl_add_)
    call acsrf%mv_from_coo(acoof,info)
  end if


  do i=1,size(adiag)
    if (adiag(i) /= czero) then
      adiag(i) = cone / adiag(i)
    else
      adiag(i) = cone
    end if
  end do

  if (filter_mat) call acsrf%scal(adiag,info)
  if (info == psb_success_) call acsr3%scal(adiag,info)
  if (info /= psb_success_) goto 9999


  if (p%parms%aggr_omega_alg == mld_eig_est_) then 

    if (p%parms%aggr_eig == mld_max_norm_) then 

      if (p%parms%aggr_kind == mld_biz_prol_) then 

        ! 
        ! This only works with CSR
        !
        anorm = szero
        dg    = sone
        nrw = acsr3%get_nrows()
        do i=1, nrw
          tmp = szero
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
      else
        anorm = acsr3%csnmi()
      endif
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_internal_error_,name,a_err='Invalid AM3 storage format')
        goto 9999
      end if
      omega = 4.d0/(3.d0*anorm)
      p%parms%aggr_omega_val = omega 

    else 
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='invalid mld_aggr_eig_')
      goto 9999
    end if

  else if (p%parms%aggr_omega_alg == mld_user_choice_) then 

    omega = p%parms%aggr_omega_val 

  else if (p%parms%aggr_omega_alg /= mld_user_choice_) then 
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
          acsrf%val(j) = cone - omega*acsrf%val(j) 
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
    call psb_symbmm(acsrf,acsr4,acsr1,info)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='symbmm 1')
      goto 9999
    end if

    call psb_numbmm(acsrf,acsr4,acsr1)

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
          acsr3%val(j) = cone - omega*acsr3%val(j) 
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
    call psb_symbmm(acsr3,acsr4,acsr1,info)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='symbmm 1')
      goto 9999
    end if

    call psb_numbmm(acsr3,acsr4,acsr1)

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done NUMBMM 1'

  end if
  call acsr4%free()
  call acsr1%set_dupl(psb_dupl_add_)

  call am1%mv_from(acsr1)
  if (ml_global_nmb) then 
    !
    ! Now we have to gather the halo of am1, and add it to itself
    ! to multiply it by A,
    !
    call psb_sphalo(am1,desc_a,am4,info,&
         & colcnv=.false.,rowscale=.true.)
    if (info == psb_success_) call psb_rwextd(ncol,am1,info,b=am4)      
    if (info == psb_success_) call am4%free()
  else 
    call psb_rwextd(ncol,am1,info)
  endif
  if(info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='Halo of am1')
    goto 9999
  end if

  call psb_symbmm(a,am1,am3,info)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='symbmm 2')
    goto 9999
  end if

  call psb_numbmm(a,am1,am3)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done NUMBMM 2',p%parms%aggr_kind, mld_smooth_prol_

  if  (p%parms%aggr_kind == mld_smooth_prol_) then 
    call am2%transp(am1)
    call am2%mv_to(acoo2)
    nzl = acoo2%get_nzeros()
    i=0
    !
    ! Now we have to fix this.  The only rows of B that are correct 
    ! are those corresponding to "local" aggregates, i.e. indices in ilaggr(:)
    !
    do k=1, nzl
      if ((naggrm1 < acoo2%ia(k)) .and.(acoo2%ia(k) <= naggrp1)) then
        i = i+1
        acoo2%val(i) = acoo2%val(k)
        acoo2%ia(i)  = acoo2%ia(k)
        acoo2%ja(i)  = acoo2%ja(k)
      end if
    end do
    call acoo2%set_nzeros(i)
    call acoo2%trim()
    call am2%mv_from(acoo2)
    call am2%cscnv(info,type='csr',dupl=psb_dupl_add_)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='spcnv am2')
      goto 9999
    end if
  else
    call am2%transp(am1)
  endif
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'starting sphalo/ rwxtd'

  if (p%parms%aggr_kind == mld_smooth_prol_) then 
    ! am2 = ((i-wDA)Ptilde)^T
    call psb_sphalo(am3,desc_a,am4,info,&
         & colcnv=.false.,rowscale=.true.)
    if (info == psb_success_) call psb_rwextd(ncol,am3,info,b=am4)      
    if (info == psb_success_) call am4%free()
  else if  (p%parms%aggr_kind == mld_biz_prol_) then 
    call psb_rwextd(ncol,am3,info)
  endif
  if(info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='Extend am3')
    goto 9999
  end if


  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'starting symbmm 3'
  call psb_symbmm(am2,am3,b,info)
  if (info == psb_success_) call psb_numbmm(am2,am3,b)
  if (info == psb_success_) call am3%free()
  if (info == psb_success_) call b%cscnv(info,type='coo',dupl=psb_dupl_add_)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='Build b = am2 x am3')
    goto 9999
  end if



  select case(p%parms%aggr_kind)

  case(mld_smooth_prol_) 

    select case(p%parms%coarse_mat)

    case(mld_distr_mat_) 

      nzac = b%get_nzeros()
      nzl =  nzac
      call b%mv_to(bcoo)

      if (info == psb_success_) call psb_cdall(ictxt,p%desc_ac,info,nl=nlaggr(me+1))
      if (info == psb_success_) call psb_cdins(nzl,bcoo%ia,bcoo%ja,p%desc_ac,info)
      if (info == psb_success_) call psb_cdasb(p%desc_ac,info)
      if (info == psb_success_) call psb_glob_to_loc(bcoo%ia(1:nzl),p%desc_ac,info,iact='I')
      if (info == psb_success_) call psb_glob_to_loc(bcoo%ja(1:nzl),p%desc_ac,info,iact='I')
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,a_err='Creating p%desc_ac and converting ac')
        goto 9999
      end if
      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & 'Assembld aux descr. distr.'
      call p%ac%mv_from(bcoo)

      call p%ac%set_nrows(p%desc_ac%get_local_rows())
      call p%ac%set_ncols(p%desc_ac%get_local_cols())
      call p%ac%set_asb()

      if (info == psb_success_) deallocate(nzbr,idisp,stat=info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_sp_free')
        goto 9999
      end if

      if (np>1) then 
        call am1%mv_to(acsr1)
        nzl = acsr1%get_nzeros()
        call psb_glob_to_loc(acsr1%ja(1:nzl),p%desc_ac,info,'I')
        if(info /= psb_success_) then
          call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_glob_to_loc')
          goto 9999
        end if
        call am1%mv_from(acsr1)
      endif
      call am1%set_ncols(p%desc_ac%get_local_cols())

      if (np>1) then 
        call am2%cscnv(info,type='coo',dupl=psb_dupl_add_)
        call am2%mv_to(acoo2)
        nzl = acoo2%get_nzeros()
        if (info == psb_success_) call psb_glob_to_loc(acoo2%ia(1:nzl),p%desc_ac,info,'I')
        call acoo2%set_dupl(psb_dupl_add_)
        if (info == psb_success_) call am2%mv_from(acoo2)
        if (info == psb_success_) call am2%cscnv(info,type='csr')        
        if(info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_,name,a_err='Converting am2 to local')
          goto 9999
        end if
      end if
      call am2%set_nrows(p%desc_ac%get_local_cols())

      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & 'Done ac '

    case(mld_repl_mat_) 
      !
      !
      call psb_cdall(ictxt,p%desc_ac,info,mg=ntaggr,repl=.true.)
      if (info == psb_success_) call psb_cdasb(p%desc_ac,info)
      if (info == psb_success_) &
           & call psb_gather(p%ac,b,p%desc_ac,info,dupl=psb_dupl_add_,keeploc=.false.)

      if (info /= psb_success_) goto 9999

      deallocate(nzbr,idisp,stat=info)
      if (info /= psb_success_) then 
        info = psb_err_alloc_dealloc_
        call psb_errpush(info,name)
        goto 9999
      end if
    case default 
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='invalid mld_coarse_mat_')
      goto 9999
    end select


  case(mld_biz_prol_) 

    select case(p%parms%coarse_mat)

    case(mld_distr_mat_) 

      call psb_move_alloc(b,p%ac,info)
      if (info == psb_success_) call psb_cdall(ictxt,p%desc_ac,info,nl=naggr)
      if (info == psb_success_) call psb_cdasb(p%desc_ac,info)
      if (info /=  psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='Build desc_ac, ac')
        goto 9999
      end if


    case(mld_repl_mat_) 
      !
      !
      call psb_cdall(ictxt,p%desc_ac,info,mg=ntaggr,repl=.true.)
      if (info == psb_success_) call psb_cdasb(p%desc_ac,info)
      if(info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_cdall')
        goto 9999
      end if
      call psb_gather(p%ac,b,p%desc_ac,info,dupl=psb_dupl_add_,keeploc=.false.)
      if(info /= psb_success_) goto 9999        

      deallocate(nzbr,idisp,stat=info)
      if (info /= psb_success_) then 
        info = psb_err_alloc_dealloc_
        call psb_errpush(info,name)
        goto 9999
      end if

    case default 
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='invalid mld_coarse_mat_')
      goto 9999
    end select

    deallocate(nzbr,idisp,stat=info)
    if (info /= psb_success_) then 
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    end if

  case default 
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='invalid mld_smooth_prol_')
    goto 9999

  end select

  call p%ac%cscnv(info,type='csr',dupl=psb_dupl_add_)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='spcnv')
    goto 9999
  end if

  !
  ! Copy the prolongation/restriction matrices into the descriptor map.
  !  am2 => PR^T   i.e. restriction  operator
  !  am1 => PR     i.e. prolongation operator
  !  
  p%map = psb_linmap(psb_map_aggr_,desc_a,&
       & p%desc_ac,am2,am1,ilaggr,nlaggr)
  if (info == psb_success_) call am1%free()
  if (info == psb_success_) call am2%free()
  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='sp_Free')
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



end subroutine mld_caggrmat_smth_asb
