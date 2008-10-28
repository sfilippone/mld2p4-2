!!$ 
!!$ 
!!$                           MLD2P4  version 1.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 2.2)
!!$  
!!$  (C) Copyright 2008
!!$
!!$                      Salvatore Filippone  University of Rome Tor Vergata       
!!$                      Alfredo Buttari      University of Rome Tor Vergata
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
! File: mld_caggrmat_raw_asb.F90
!
! Subroutine: mld_caggrmat_raw_asb
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
!  specified by the user through mld_cprecinit and mld_cprecset.
!
!  For details see
!    P. D'Ambra, D. di Serafino and  S. Filippone, On the development of
!    PSBLAS-based parallel two-level Schwarz preconditioners, Appl. Num. Math.,
!    57 (2007), 1181-1196.
!
!
!
! Arguments:
!    a          -  type(psb_cspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!    p          -  type(mld_c_onelev_prec_type), input/output.
!                  The one-level preconditioner data structure containing the local
!                  part of the base preconditioner to be built as well as the
!                  aggregate matrices.
!    info       -  integer, output.
!                  Error code.
!
subroutine mld_caggrmat_raw_asb(a,desc_a,p,info)
  use psb_base_mod
  use mld_inner_mod, mld_protect_name => mld_caggrmat_raw_asb

#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

! Arguments
  type(psb_cspmat_type), intent(in)               :: a
  type(psb_desc_type), intent(in)                 :: desc_a
  type(mld_c_onelev_prec_type), intent(inout), target  :: p
  integer, intent(out)                       :: info

! Local variables
  integer ::ictxt,np,me, err_act, icomm
  character(len=20) :: name
  type(psb_cspmat_type)  :: b
  integer, allocatable   :: nzbr(:), idisp(:)
  type(psb_cspmat_type)  :: am1,am2
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

  call psb_nullify_sp(am1)
  call psb_nullify_sp(am2)


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

  if (p%iprcparm(mld_coarse_mat_) == mld_repl_mat_) then
    call psb_sp_all(ncol,ntaggr,am1,ncol,info)
  else
    call psb_sp_all(ncol,naggr,am1,ncol,info)
  end if

  if (info /= 0) then
    call psb_errpush(4010,name,a_err='spall')
    goto 9999
  end if

  do i=1,nrow
    am1%aspk(i) = cone
    am1%ia1(i)  = i
    am1%ia2(i)  = p%mlia(i)  
  end do
  am1%infoa(psb_nnz_) = nrow

  call psb_spcnv(am1,info,afmt='csr',dupl=psb_dupl_add_)
  call psb_transp(am1,am2)


  call psb_sp_clip(a,b,info,jmax=nrow)
  if(info /= 0) then
    call psb_errpush(4010,name,a_err='spclip')
    goto 9999
  end if
  ! Out from sp_clip is always in COO, but just in case..
  if (psb_tolower(b%fida) /= 'coo') then 
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

    call psb_cdall(ictxt,p%desc_ac,info,mg=ntaggr,repl=.true.)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psb_cdall')
      goto 9999
    end if

    nzbr(:) = 0
    nzbr(me+1) = nzt
    call psb_sum(ictxt,nzbr(1:np))
    nzac = sum(nzbr)
    
    call psb_sp_all(ntaggr,ntaggr,p%ac,nzac,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='sp_all')
      goto 9999
    end if

    do ip=1,np
      idisp(ip) = sum(nzbr(1:ip-1))
    enddo
    ndx = nzbr(me+1) 

    call mpi_allgatherv(b%aspk,ndx,mpi_complex,p%ac%aspk,nzbr,idisp,&
         & mpi_complex,icomm,info)
    call mpi_allgatherv(b%ia1,ndx,mpi_integer,p%ac%ia1,nzbr,idisp,&
         & mpi_integer,icomm,info)
    call mpi_allgatherv(b%ia2,ndx,mpi_integer,p%ac%ia2,nzbr,idisp,&
         & mpi_integer,icomm,info)
    if(info /= 0) then
      info=-1
      call psb_errpush(info,name)
      goto 9999
    end if

    p%ac%m = ntaggr
    p%ac%k = ntaggr
    p%ac%infoa(psb_nnz_) = nzac
    p%ac%fida='COO'
    p%ac%descra='GUN'
    call psb_spcnv(p%ac,info,afmt='coo',dupl=psb_dupl_add_)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='sp_free')
      goto 9999
    end if

  else if (p%iprcparm(mld_coarse_mat_) == mld_distr_mat_) then 

    call psb_cdall(ictxt,p%desc_ac,info,nl=naggr)
    if (info == 0) call psb_cdasb(p%desc_ac,info)
    if (info == 0) call psb_sp_clone(b,p%ac,info)
    if(info /= 0) then
      call psb_errpush(4001,name,a_err='Build ac, desc_ac')
      goto 9999
    end if
    call psb_sp_free(b,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='sp_free')
      goto 9999
    end if

  else
    info = 4001
    call psb_errpush(4001,name,a_err='invalid mld_coarse_mat_')
    goto 9999
  end if

  deallocate(nzbr,idisp)
  
  call psb_spcnv(p%ac,info,afmt='csr',dupl=psb_dupl_add_)
  if(info /= 0) then
    call psb_errpush(4010,name,a_err='spcnv')
    goto 9999
  end if

  !
  ! Copy the prolongation/restriction matrices into the descriptor map.
  !  am2 => PR^T   i.e. restriction  operator
  !  am1 => PR     i.e. prolongation operator
  !  
  p%map_desc = psb_inter_desc(psb_map_aggr_,desc_a,&
       & p%desc_ac,am2,am1)
  if (info == 0) call psb_sp_free(am1,info)
  if (info == 0) call psb_sp_free(am2,info)
  if(info /= 0) then
    call psb_errpush(4010,name,a_err='sp_Free')
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

end subroutine mld_caggrmat_raw_asb
