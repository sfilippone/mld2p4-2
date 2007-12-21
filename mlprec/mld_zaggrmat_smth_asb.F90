!!$ 
!!$ 
!!$                                MLD2P4
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS v.2.0)
!!$  
!!$  (C) Copyright 2007  Alfredo Buttari      University of Rome Tor Vergata
!!$                      Pasqua D'Ambra       ICAR-CNR, Naples
!!$                      Daniela di Serafino  Second University of Naples
!!$                      Salvatore Filippone  University of Rome Tor Vergata       
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
!  by using a Galerkin approach, i.e.
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
!  according to the value of p%iprcparm(mld_aggr_eig_), specified by the user
!  through mld_dprecinit and mld_dprecset.
!
!  This routine can also build A_C according to a "bizarre" aggregation algorithm,
!  using a "naive" prolongator proposed by the authors of MLD2P4. However, this
!  prolongator still requires a deep analysis and testing and its use is not
!  recommended.
!
!  The coarse-level matrix A_C is distributed among the parallel processes or
!  replicated on each of them, according to the value of p%iprcparm(mld_coarse_mat_),
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
!    a          -  type(psb_zspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!    ac         -  type(psb_zspmat_type), output.
!                  The sparse matrix structure containing the local part of
!                  the coarse-level matrix.
!    desc_ac    -  type(psb_desc_type), output.
!                  The communication descriptor of the coarse-level matrix.
!    p          -  type(mld_zbaseprc_type), input/output.
!                  The base preconditioner data structure containing the local
!                  part of the base preconditioner to be built.
!    info       -  integer, output.
!                  Error code.
!
subroutine mld_zaggrmat_smth_asb(a,desc_a,ac,desc_ac,p,info)
  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_zaggrmat_smth_asb

#ifdef MPI_MOD
  use mpi
#endif
  implicit none 
#ifdef MPI_H
  include 'mpif.h'
#endif

! Arguments
  type(psb_zspmat_type), intent(in), target  :: a
  type(psb_desc_type), intent(in)            :: desc_a
  type(psb_zspmat_type), intent(inout), target :: ac    
  type(psb_desc_type), intent(inout)         :: desc_ac 
  type(mld_zbaseprc_type), intent(inout), target  :: p
  integer, intent(out)                       :: info

! Local variables
  type(psb_zspmat_type)  :: b
  integer, pointer :: nzbr(:), idisp(:)
  integer :: nrow, nglob, ncol, ntaggr, nzac, ip, ndx,&
       & naggr, nzl,naggrm1,naggrp1, i, j, k
  integer ::ictxt,np,me, err_act, icomm
  character(len=20) :: name
  type(psb_zspmat_type), pointer  :: am1,am2
  type(psb_zspmat_type) :: am3,am4
  logical       :: ml_global_nmb
  integer                       :: nz
  integer, allocatable          :: ia(:), ja(:)
  complex(kind(1.d0)), allocatable :: val(:)
  integer            :: debug_level, debug_unit
  integer, parameter :: ncmax=16
  real(kind(1.d0))   :: omega, anorm, tmp, dg

  name='mld_aggrmat_smth_asb'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = psb_cd_get_context(desc_a)
  icomm = psb_cd_get_mpic(desc_a)
  ictxt = psb_cd_get_context(desc_a)

  call psb_info(ictxt, me, np)


  call psb_nullify_sp(b)
  call psb_nullify_sp(am3)
  call psb_nullify_sp(am4)

  am2 => p%av(mld_sm_pr_t_)
  am1 => p%av(mld_sm_pr_)
  call psb_nullify_sp(am1)
  call psb_nullify_sp(am2)

  nglob = psb_cd_get_global_rows(desc_a)
  nrow  = psb_cd_get_local_rows(desc_a)
  ncol  = psb_cd_get_local_cols(desc_a)

  naggr  = p%nlaggr(me+1)
  ntaggr = sum(p%nlaggr)

  allocate(nzbr(np), idisp(np),stat=info)
  if (info /= 0) then 
    info=4025
    call psb_errpush(info,name,i_err=(/2*np,0,0,0,0/),&
         & a_err='integer')
    goto 9999      
  end if

  naggrm1 = sum(p%nlaggr(1:me))
  naggrp1 = sum(p%nlaggr(1:me+1))
  ml_global_nmb = ( (p%iprcparm(mld_aggr_kind_) == mld_smooth_prol_).or.&
       & ( (p%iprcparm(mld_aggr_kind_) == mld_biz_prol_).and.&
       &    (p%iprcparm(mld_coarse_mat_) == mld_repl_mat_)) ) 

  if (ml_global_nmb) then 
    p%mlia(1:nrow) = p%mlia(1:nrow) + naggrm1
    call psb_halo(p%mlia,desc_a,info)

    if (info /= 0) then
      call psb_errpush(4010,name,a_err='psb_halo')
      goto 9999
    end if
  end if

  ! naggr: number of local aggregates
  ! nrow: local rows. 
  ! 
  allocate(p%dorig(nrow),stat=info)

  if (info /= 0) then 
    info=4025
    call psb_errpush(info,name,i_err=(/nrow,0,0,0,0/),&
         & a_err='real(kind(1.d0))')
    goto 9999      
  end if

  ! Get diagonal D
  call psb_sp_getdiag(a,p%dorig,info)
  if(info /= 0) then
    call psb_errpush(4010,name,a_err='sp_getdiag')
    goto 9999
  end if

  do i=1,size(p%dorig)
    if (p%dorig(i) /= zzero) then
      p%dorig(i) = zone / p%dorig(i)
    else
      p%dorig(i) = zone
    end if
  end do

  ! 1. Allocate Ptilde in sparse matrix form 
  am4%fida='COO'
  am4%m=ncol
  if (ml_global_nmb) then 
    am4%k=ntaggr
    call psb_sp_all(ncol,ntaggr,am4,ncol,info)
  else 
    am4%k=naggr
    call psb_sp_all(ncol,naggr,am4,ncol,info)
  endif

  if (info /= 0) then
    call psb_errpush(4010,name,a_err='spall')
    goto 9999
  end if

  if (ml_global_nmb) then 
    do i=1,ncol
      am4%aspk(i) = zone
      am4%ia1(i)  = i
      am4%ia2(i)  = p%mlia(i)  
    end do
    am4%infoa(psb_nnz_) = ncol
  else
    do i=1,nrow
      am4%aspk(i) = zone
      am4%ia1(i)  = i
      am4%ia2(i)  = p%mlia(i)  
    end do
    am4%infoa(psb_nnz_) = nrow
  endif


  call psb_spcnv(am4,info,afmt='csr',dupl=psb_dupl_add_)
  if (info==0) call psb_spcnv(a,am3,info,afmt='csr',dupl=psb_dupl_add_)
  if (info /= 0) then
    call psb_errpush(4010,name,a_err='spcnv')
    goto 9999
  end if

  !
  ! WARNING: the cycles below assume that AM3 does have 
  ! its diagonal elements stored explicitly!!! 
  ! Should we switch to something safer? 
  !
  call psb_sp_scal(am3,p%dorig,info)
  if (info /= 0) goto 9999

  if (p%iprcparm(mld_aggr_eig_) == mld_max_norm_) then 

    if (p%iprcparm(mld_aggr_kind_) == mld_biz_prol_) then 

      ! 
      ! This only works with CSR.
      !
      if (toupper(am3%fida)=='CSR') then 
        anorm = dzero
        dg    = done
        do i=1,am3%m
          tmp = dzero
          do j=am3%ia2(i),am3%ia2(i+1)-1
            if (am3%ia1(j) <= am3%m) then 
              tmp = tmp + abs(am3%aspk(j))
            endif
            if (am3%ia1(j) == i ) then 
              dg = abs(am3%aspk(j))
            end if
          end do
          anorm = max(anorm,tmp/dg) 
        enddo
        
        call psb_amx(ictxt,anorm)     
      else
        info = 4001
      endif
    else
      anorm = psb_spnrmi(am3,desc_a,info)
    endif
    if (info /= 0) then 
      call psb_errpush(4001,name,a_err='Invalid AM3 storage format')
      goto 9999
    end if
    omega = 4.d0/(3.d0*anorm)
    p%dprcparm(mld_aggr_damp_) = omega 

  else if (p%iprcparm(mld_aggr_eig_) == mld_user_choice_) then 

    omega = p%dprcparm(mld_aggr_damp_) 

  else if (p%iprcparm(mld_aggr_eig_) /= mld_user_choice_) then 
    info = 4001
    call psb_errpush(info,name,a_err='invalid mld_aggr_eig_')
    goto 9999
  end if


  if (toupper(am3%fida)=='CSR') then 
    do i=1,am3%m
      do j=am3%ia2(i),am3%ia2(i+1)-1
        if (am3%ia1(j) == i) then 
          am3%aspk(j) = zone - omega*am3%aspk(j) 
        else
          am3%aspk(j) = - omega*am3%aspk(j) 
        end if
      end do
    end do
  else 
    call psb_errpush(4001,name,a_err='Invalid AM3 storage format')
    goto 9999
  end if

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done gather, going for SYMBMM 1'
  !
  ! Symbmm90 does the allocation for its result.
  ! 
  ! am1 = (i-wDA)Ptilde
  ! Doing it this way means to consider diag(Ai)
  ! 
  !
  call psb_symbmm(am3,am4,am1,info)
  if(info /= 0) then
    call psb_errpush(4010,name,a_err='symbmm 1')
    goto 9999
  end if

  call psb_numbmm(am3,am4,am1)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done NUMBMM 1'

  call psb_sp_free(am4,info)
  if(info /= 0) then
    call psb_errpush(4010,name,a_err='sp_free')
    goto 9999
  end if

  if (ml_global_nmb) then 
    !
    ! Now we have to gather the halo of am1, and add it to itself
    ! to multiply it by A,
    !
    call psb_sphalo(am1,desc_a,am4,info,&
         & colcnv=.false.,rowscale=.true.)
    if (info == 0) call psb_rwextd(ncol,am1,info,b=am4)      
    if (info == 0) call psb_sp_free(am4,info)
  else 
    call psb_rwextd(ncol,am1,info)
  endif
  if(info /= 0) then
    call psb_errpush(4001,name,a_err='Halo of am1')
    goto 9999
  end if
  
  call psb_symbmm(a,am1,am3,info)
  if(info /= 0) then
    call psb_errpush(4010,name,a_err='symbmm 2')
    goto 9999
  end if

  call psb_numbmm(a,am1,am3)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done NUMBMM 2'

  if  (p%iprcparm(mld_aggr_kind_) == mld_smooth_prol_) then 
    call psb_transp(am1,am2,fmt='COO')
    nzl = am2%infoa(psb_nnz_)
    i=0
    !
    ! Now we have to fix this.  The only rows of B that are correct 
    ! are those corresponding to "local" aggregates, i.e. indices in p%mlia(:)
    !
    do k=1, nzl
      if ((naggrm1 < am2%ia1(k)) .and.(am2%ia1(k) <= naggrp1)) then
        i = i+1
        am2%aspk(i) = am2%aspk(k)
        am2%ia1(i)  = am2%ia1(k)
        am2%ia2(i)  = am2%ia2(k)
      end if
    end do
    am2%infoa(psb_nnz_) = i
    call psb_spcnv(am2,info,afmt='csr',dupl=psb_dupl_add_)
    if (info /=0) then 
      call psb_errpush(4010,name,a_err='spcnv am2')
      goto 9999
    end if
  else
    call psb_transp(am1,am2)
  endif
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'starting sphalo/ rwxtd'

  if (p%iprcparm(mld_aggr_kind_) == mld_smooth_prol_) then 
    ! am2 = ((i-wDA)Ptilde)^T
    call psb_sphalo(am3,desc_a,am4,info,&
         & colcnv=.false.,rowscale=.true.)
    if (info == 0) call psb_rwextd(ncol,am3,info,b=am4)      
    if (info == 0) call psb_sp_free(am4,info)
  else if  (p%iprcparm(mld_aggr_kind_) == mld_biz_prol_) then 
    call psb_rwextd(ncol,am3,info)
  endif
  if(info /= 0) then
    call psb_errpush(4001,name,a_err='Extend am3')
    goto 9999
  end if


  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'starting symbmm 3'
  call psb_symbmm(am2,am3,b,info)
  if (info == 0) call psb_numbmm(am2,am3,b)
  if (info == 0) call psb_sp_free(am3,info)
  if (info == 0) call psb_spcnv(b,info,afmt='coo',dupl=psb_dupl_add_)
  if (info /= 0) then
    call psb_errpush(4001,name,a_err='Build b = am2 x am3')
    goto 9999
  end if



  select case(p%iprcparm(mld_aggr_kind_))

  case(mld_smooth_prol_) 

    select case(p%iprcparm(mld_coarse_mat_))

    case(mld_distr_mat_) 

      call psb_sp_clone(b,ac,info)
      nzac = ac%infoa(psb_nnz_) 
      nzl =  ac%infoa(psb_nnz_) 
      if (info == 0) call psb_cdall(ictxt,desc_ac,info,nl=p%nlaggr(me+1))
      if (info == 0) call psb_cdins(nzl,ac%ia1,ac%ia2,desc_ac,info)
      if (info == 0) call psb_cdasb(desc_ac,info)
      if (info == 0) call psb_glob_to_loc(ac%ia1(1:nzl),desc_ac,info,iact='I')
      if (info == 0) call psb_glob_to_loc(ac%ia2(1:nzl),desc_ac,info,iact='I')
      if (info /= 0) then
        call psb_errpush(4001,name,a_err='Creating desc_ac and converting ac')
        goto 9999
      end if
      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & 'Assembld aux descr. distr.'


      ac%m=desc_ac%matrix_data(psb_n_row_)
      ac%k=desc_ac%matrix_data(psb_n_col_)
      ac%fida='COO'
      ac%descra='GUN'

      call psb_sp_free(b,info)
      if (info == 0) deallocate(nzbr,idisp,stat=info)
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='psb_sp_free')
        goto 9999
      end if

      if (np>1) then 
        nzl = psb_sp_get_nnzeros(am1)
        call psb_glob_to_loc(am1%ia1(1:nzl),desc_ac,info,'I')
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_glob_to_loc')
          goto 9999
        end if
      endif
      am1%k=desc_ac%matrix_data(psb_n_col_)

      if (np>1) then 
        call psb_spcnv(am2,info,afmt='coo',dupl=psb_dupl_add_)
        nzl = am2%infoa(psb_nnz_) 
        if (info == 0) call psb_glob_to_loc(am2%ia1(1:nzl),desc_ac,info,'I')
        if (info == 0) call psb_spcnv(am2,info,afmt='csr',dupl=psb_dupl_add_)        
        if(info /= 0) then
          call psb_errpush(4001,name,a_err='Converting am2 to local')
          goto 9999
        end if
      end if
      am2%m=desc_ac%matrix_data(psb_n_col_)

      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & 'Done ac '

    case(mld_repl_mat_) 
      !
      !
      call psb_cdall(ictxt,desc_ac,info,mg=ntaggr,repl=.true.)
      nzbr(:) = 0
      nzbr(me+1) = b%infoa(psb_nnz_)

      call psb_sum(ictxt,nzbr(1:np))
      nzac = sum(nzbr)
      if (info == 0) call psb_sp_all(ntaggr,ntaggr,ac,nzac,info)
      if (info /= 0) goto 9999

      do ip=1,np
        idisp(ip) = sum(nzbr(1:ip-1))
      enddo
      ndx = nzbr(me+1) 

      call mpi_allgatherv(b%aspk,ndx,mpi_double_complex,ac%aspk,nzbr,idisp,&
           & mpi_double_complex,icomm,info)
      if (info == 0) call mpi_allgatherv(b%ia1,ndx,mpi_integer,ac%ia1,nzbr,idisp,&
           & mpi_integer,icomm,info)
      if (info == 0) call mpi_allgatherv(b%ia2,ndx,mpi_integer,ac%ia2,nzbr,idisp,&
           & mpi_integer,icomm,info)

      if (info /= 0) then 
        call psb_errpush(4001,name,a_err=' from mpi_allgatherv')
        goto 9999
      end if

      ac%m = ntaggr
      ac%k = ntaggr
      ac%infoa(psb_nnz_) = nzac
      ac%fida='COO'
      ac%descra='GUN'
      call psb_spcnv(ac,info,afmt='coo',dupl=psb_dupl_add_)
      if(info /= 0) goto 9999
      call psb_sp_free(b,info)
      if(info /= 0) goto 9999

      deallocate(nzbr,idisp,stat=info)
      if (info /= 0) then 
        info = 4000
        call psb_errpush(info,name)
        goto 9999
      end if
    case default 
      info = 4001
      call psb_errpush(info,name,a_err='invalid mld_coarse_mat_')
      goto 9999
    end select


  case(mld_biz_prol_) 

    select case(p%iprcparm(mld_coarse_mat_))

    case(mld_distr_mat_) 

      call psb_sp_clone(b,ac,info)
      if (info == 0) call psb_cdall(ictxt,desc_ac,info,nl=naggr)
      if (info == 0) call psb_cdasb(desc_ac,info)
      if (info == 0) call psb_sp_free(b,info)
      if (info /=  0) then
        call psb_errpush(4010,name,a_err='Build desc_ac, ac')
        goto 9999
      end if


    case(mld_repl_mat_) 
      !
      !
      call psb_cdall(ictxt,desc_ac,info,mg=ntaggr,repl=.true.)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='psb_cdall')
        goto 9999
      end if

      nzbr(:) = 0
      nzbr(me+1) = b%infoa(psb_nnz_)
      call psb_sum(ictxt,nzbr(1:np))
      nzac = sum(nzbr)
      call psb_sp_all(ntaggr,ntaggr,ac,nzac,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='psb_sp_all')
        goto 9999
      end if

      do ip=1,np
        idisp(ip) = sum(nzbr(1:ip-1))
      enddo
      ndx = nzbr(me+1) 

      call mpi_allgatherv(b%aspk,ndx,mpi_double_complex,ac%aspk,nzbr,idisp,&
           & mpi_double_complex,icomm,info)
      if (info == 0) call mpi_allgatherv(b%ia1,ndx,mpi_integer,ac%ia1,nzbr,idisp,&
           & mpi_integer,icomm,info)
      if (info == 0) call mpi_allgatherv(b%ia2,ndx,mpi_integer,ac%ia2,nzbr,idisp,&
           & mpi_integer,icomm,info)
      if (info /= 0) then 
        call psb_errpush(4001,name,a_err=' from mpi_allgatherv')
        goto 9999
      end if


      ac%m = ntaggr
      ac%k = ntaggr
      ac%infoa(psb_nnz_) = nzac
      ac%fida='COO'
      ac%descra='GUN'
      call psb_spcnv(ac,info,afmt='coo',dupl=psb_dupl_add_)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='spcnv')
        goto 9999
      end if
      call psb_sp_free(b,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='psb_sp_free')
        goto 9999
      end if

    case default 
      info = 4001
      call psb_errpush(info,name,a_err='invalid mld_coarse_mat_')
      goto 9999
    end select

    deallocate(nzbr,idisp,stat=info)
    if (info /= 0) then 
      info = 4000
      call psb_errpush(info,name)
      goto 9999
    end if

  case default 
    info = 4001
    call psb_errpush(info,name,a_err='invalid mld_smooth_prol_')
    goto 9999

  end select

  call psb_spcnv(ac,info,afmt='csr',dupl=psb_dupl_add_)
  if(info /= 0) then
    call psb_errpush(4010,name,a_err='spcnv')
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



end subroutine mld_zaggrmat_smth_asb
