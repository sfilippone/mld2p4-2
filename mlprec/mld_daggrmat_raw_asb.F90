!!$ 
!!$ 
!!$                    MD2P4
!!$    Multilevel Domain Decomposition Parallel Preconditioner Package for PSBLAS
!!$                      for 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$                       Daniela di Serafino    Second University of Naples
!!$                       Pasqua D'Ambra         ICAR-CNR                      
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the MD2P4 group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MD2P4 GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
subroutine mld_daggrmat_raw_asb(a,desc_a,ac,desc_ac,p,info)
  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_daggrmat_raw_asb

#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif
  type(psb_dspmat_type), intent(in), target  :: a
  type(mld_dbaseprc_type), intent(inout), target  :: p
  type(psb_dspmat_type), intent(inout), target :: ac
  type(psb_desc_type), intent(in)            :: desc_a
  type(psb_desc_type), intent(inout)         :: desc_ac
  integer, intent(out)                       :: info

  logical, parameter :: aggr_dump=.false.
  integer ::ictxt,np,me, err_act, icomm
  character(len=20) :: name, ch_err

  type(psb_dspmat_type)  :: b
  integer, pointer :: nzbr(:), idisp(:)
  integer :: nrow, nglob, ncol, ntaggr, nzac, ip, ndx,&
       & naggr, nzt,jl,nzl,nlr, naggrm1, i, j, k

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

    call psb_cdrep(ntaggr,ictxt,desc_ac,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psb_cdrep')
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

    call mpi_allgatherv(b%aspk,ndx,mpi_double_precision,ac%aspk,nzbr,idisp,&
         & mpi_double_precision,icomm,info)
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

end subroutine mld_daggrmat_raw_asb
