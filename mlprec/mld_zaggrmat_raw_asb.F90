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
! File: mld_zaggrmat_raw_asb.F90
!
! Subroutine: mld_zaggrmat_raw_asb
! Version:    complex
!
!  This routine builds a coarse-level matrix A_C from a fine-level matrix A
!  by using a Galerkin approach, i.e.
!
!                               A_C = P_C^T A P_C,
!
!  where P_C is the piecewise constant interpolation operator corresponding
!  the fine-to-coarse level mapping built by mld_aggrmap_bld.
! 
!  The coarse-level matrix A_C is distributed among the parallel processes or
!  replicated on each of them, according to the value of p%iprcparm(mld_coarse_mat_),
!  specified by the user through mld_dprecinit and mld_dprecset.
!
!  For details see
!    P. D'Ambra, D. di Serafino and  S. Filippone, On the development of
!    PSBLAS-based parallel two-level Schwarz preconditioners, Appl. Num. Math.,
!    57 (2007), 1181-1196.
!
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
subroutine mld_zaggrmat_raw_asb(a,desc_a,ac,desc_ac,p,info)
  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_zaggrmat_raw_asb

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
  logical, parameter :: aggr_dump=.false.
  integer ::ictxt,np,me, err_act,icomm
  character(len=20) :: name
  type(psb_zspmat_type)  :: b
  integer, pointer :: nzbr(:), idisp(:)
  integer :: nrow, nglob, ncol, ntaggr, nzac, ip, ndx,&
       & naggr, nzt,naggrm1, i

  name='mld_aggrmat_raw_asb'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  call psb_nullify_sp(b)

  ictxt = psb_cd_get_context(desc_a)
  icomm = psb_cd_get_mpic(desc_a)
  call psb_info(ictxt, me, np)
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

  naggrm1=sum(p%nlaggr(1:me))

  if (p%iprcparm(mld_coarse_mat_) == mld_repl_mat_) then
    do i=1, nrow
      p%mlia(i) = p%mlia(i) + naggrm1
    end do
    call psb_halo(p%mlia,desc_a,info)
  end if

  if(info /= 0) then
    call psb_errpush(4010,name,a_err='psb_halo')
    goto 9999
  end if

  call psb_sp_clip(a,b,info,jmax=nrow)
  if(info /= 0) then
    call psb_errpush(4010,name,a_err='spclip')
    goto 9999
  end if
  ! Out from sp_clip is always in COO, but just in case..
  if (tolower(b%fida) /= 'coo') then 
    call psb_errpush(4010,name,a_err='spclip NOT COO')
    goto 9999
  end if
    
  nzt = psb_sp_get_nnzeros(b)
  do i=1, nzt 
    b%ia1(i) = p%mlia(b%ia1(i))
    b%ia2(i) = p%mlia(b%ia2(i))
  enddo
  b%m = naggr
  b%k = naggr
  ! This is to minimize data exchange
  call psb_spcnv(b,info,afmt='coo',dupl=psb_dupl_add_)

  if (p%iprcparm(mld_coarse_mat_) == mld_repl_mat_) then 

    call psb_cdall(ictxt,desc_ac,info,mg=ntaggr,repl=.true.)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psb_cdall')
      goto 9999
    end if

    nzbr(:) = 0
    nzbr(me+1) = nzt
    call psb_sum(ictxt,nzbr(1:np))
    nzac = sum(nzbr)
    
    call psb_sp_all(ntaggr,ntaggr,ac,nzac,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='sp_all')
      goto 9999
    end if

    do ip=1,np
      idisp(ip) = sum(nzbr(1:ip-1))
    enddo
    ndx = nzbr(me+1) 

    call mpi_allgatherv(b%aspk,ndx,mpi_double_complex,ac%aspk,nzbr,idisp,&
         & mpi_double_complex,icomm,info)
    call mpi_allgatherv(b%ia1,ndx,mpi_integer,ac%ia1,nzbr,idisp,&
         & mpi_integer,icomm,info)
    call mpi_allgatherv(b%ia2,ndx,mpi_integer,ac%ia2,nzbr,idisp,&
         & mpi_integer,icomm,info)
    if(info /= 0) then
      info=-1
      call psb_errpush(info,name)
      goto 9999
    end if

    ac%m = ntaggr
    ac%k = ntaggr
    ac%infoa(psb_nnz_) = nzac
    ac%fida='COO'
    ac%descra='G'
    call psb_spcnv(ac,info,afmt='coo',dupl=psb_dupl_add_)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='sp_free')
      goto 9999
    end if

  else if (p%iprcparm(mld_coarse_mat_) == mld_distr_mat_) then 

    call psb_cdall(ictxt,desc_ac,info,nl=naggr)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psb_cdall')
      goto 9999
    end if
    call psb_cdasb(desc_ac,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psb_cdasb')
      goto 9999
    end if

    call psb_sp_clone(b,ac,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='spclone')
      goto 9999
    end if
    call psb_sp_free(b,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='sp_free')
      goto 9999
    end if

  else

    write(0,*) 'Unknown p%iprcparm(coarse_mat) in aggregate_sp',p%iprcparm(mld_coarse_mat_)
  end if

  deallocate(nzbr,idisp)
  
  call psb_spcnv(ac,info,afmt='csr',dupl=psb_dupl_add_)
  if(info /= 0) then
    call psb_errpush(4010,name,a_err='ipcoo2csr')
    goto 9999
  end if


  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_zaggrmat_raw_asb
