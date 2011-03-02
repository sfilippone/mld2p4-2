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
! File: mld_zaggrmat_nosmth_asb.F90
!
! Subroutine: mld_zaggrmat_nosmth_asb
! Version:    complex
!
!  This routine builds a coarse-level matrix A_C from a fine-level matrix A
!  by using the Galerkin approach, i.e.
!
!                               A_C = P_C^T A P_C,
!
!  where P_C is the piecewise constant interpolation operator corresponding
!  the fine-to-coarse level mapping built by mld_aggrmap_bld.
! 
!  The coarse-level matrix A_C is distributed among the parallel processes or
!  replicated on each of them, according to the value of p%parms%coarse_mat
!  specified by the user through mld_zprecinit and mld_zprecset.
!
!  For details see
!    P. D'Ambra, D. di Serafino and  S. Filippone, On the development of
!    PSBLAS-based parallel two-level Schwarz preconditioners, Appl. Num. Math.,
!    57 (2007), 1181-1196.
!
!
!
! Arguments:
!    a          -  type(psb_zspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!    p          -  type(mld_zonelev_type), input/output.
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
subroutine mld_zaggrmat_nosmth_asb(a,desc_a,ilaggr,nlaggr,p,info)
  use psb_sparse_mod
  use mld_z_inner_mod, mld_protect_name => mld_zaggrmat_nosmth_asb

#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

  ! Arguments
  type(psb_zspmat_type), intent(in)               :: a
  type(psb_desc_type), intent(in)                 :: desc_a
  integer, intent(inout)                          :: ilaggr(:), nlaggr(:)
  type(mld_zonelev_type), intent(inout), target  :: p
  integer, intent(out)                       :: info

  ! Local variables
  integer ::ictxt,np,me, err_act, icomm
  character(len=20) :: name
  type(psb_zspmat_type)  :: b
  integer, allocatable :: nzbr(:), idisp(:)
  type(psb_zspmat_type) :: am1,am2
  type(psb_z_coo_sparse_mat) :: acoo1, acoo2, bcoo, ac_coo
  integer :: nrow, nglob, ncol, ntaggr, nzac, ip, ndx,&
       & naggr, nzt, naggrm1, i

  name='mld_aggrmat_nosmth_asb'
  if(psb_get_errstatus().ne.0) return 
  info=psb_success_
  call psb_erractionsave(err_act)


  ictxt = psb_cd_get_context(desc_a)
  icomm = psb_cd_get_mpic(desc_a)
  call psb_info(ictxt, me, np)
  nglob = psb_cd_get_global_rows(desc_a)
  nrow  = psb_cd_get_local_rows(desc_a)
  ncol  = psb_cd_get_local_cols(desc_a)


  naggr  = nlaggr(me+1)
  ntaggr = sum(nlaggr)
  allocate(nzbr(np), idisp(np),stat=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/2*np,0,0,0,0/),&
         & a_err='integer')
    goto 9999      
  end if

  naggrm1=sum(nlaggr(1:me))

  if (p%parms%coarse_mat == mld_repl_mat_) then
    do i=1, nrow
      ilaggr(i) = ilaggr(i) + naggrm1
    end do
    call psb_halo(ilaggr,desc_a,info)
  end if

  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_halo')
    goto 9999
  end if

  if (p%parms%coarse_mat == mld_repl_mat_) then
    call acoo1%allocate(ncol,ntaggr,ncol)
  else
    call acoo1%allocate(ncol,naggr,ncol)
  end if

  do i=1,nrow
    acoo1%val(i) = zone
    acoo1%ia(i)  = i
    acoo1%ja(i)  = ilaggr(i)  
  end do

  call acoo1%set_dupl(psb_dupl_add_)
  call acoo1%set_nzeros(nrow)
  call acoo1%set_asb()
  call acoo1%fix(info)
  call acoo2%transp(acoo1) 

  call a%csclip(bcoo,info,jmax=nrow)


  nzt = bcoo%get_nzeros()
  do i=1, nzt 
    bcoo%ia(i) = ilaggr(bcoo%ia(i))
    bcoo%ja(i) = ilaggr(bcoo%ja(i))
  enddo
  call bcoo%set_nrows(naggr)
  call bcoo%set_ncols(naggr)
  call bcoo%set_dupl(psb_dupl_add_)
  call bcoo%fix(info)


  if (p%parms%coarse_mat == mld_repl_mat_) then 

    call psb_cdall(ictxt,p%desc_ac,info,mg=ntaggr,repl=.true.)
    if (info == psb_success_) call psb_cdasb(p%desc_ac,info)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_cdall')
      goto 9999
    end if

    nzbr(:) = 0
    nzbr(me+1) = nzt
    call psb_sum(ictxt,nzbr(1:np))
    nzac = sum(nzbr)

    call ac_coo%allocate(ntaggr,ntaggr,nzac)

    do ip=1,np
      idisp(ip) = sum(nzbr(1:ip-1))
    enddo
    ndx = nzbr(me+1) 

    call mpi_allgatherv(bcoo%val,ndx,mpi_double_complex,ac_coo%val,nzbr,idisp,&
         & mpi_double_precision,icomm,info)
    call mpi_allgatherv(bcoo%ia,ndx,mpi_integer,ac_coo%ia,nzbr,idisp,&
         & mpi_integer,icomm,info)
    call mpi_allgatherv(bcoo%ja,ndx,mpi_integer,ac_coo%ja,nzbr,idisp,&
         & mpi_integer,icomm,info)
    if(info /= psb_success_) then
      info=-1
      call psb_errpush(info,name)
      goto 9999
    end if
    call ac_coo%set_nzeros(nzac)
    call ac_coo%set_dupl(psb_dupl_add_)
    call ac_coo%fix(info)
    call p%ac%mv_from(ac_coo)

  else if (p%parms%coarse_mat == mld_distr_mat_) then 

    call psb_cdall(ictxt,p%desc_ac,info,nl=naggr)
    if (info == psb_success_) call psb_cdasb(p%desc_ac,info)
    call p%ac%mv_from(bcoo)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,name,a_err='Build ac, desc_ac')
      goto 9999

    end if
    
  else
    info = psb_err_internal_error_
    call psb_errpush(psb_err_internal_error_,name,a_err='invalid mld_coarse_mat_')
    goto 9999
  end if

  call bcoo%free()

  deallocate(nzbr,idisp)

  call p%ac%cscnv(info,type='csr',dupl=psb_dupl_add_)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='cscnv')
    goto 9999
  end if

  call am1%mv_from(acoo1)
  call am1%cscnv(info,type='csr',dupl=psb_dupl_add_)
  if (info == psb_success_) call am2%mv_from(acoo2)
  if (info == psb_success_) call am2%cscnv(info,type='csr',dupl=psb_dupl_add_)
  !
  ! Copy the prolongation/restriction matrices into the descriptor map.
  !  am2 => PR^T   i.e. restriction  operator
  !  am1 => PR     i.e. prolongation operator
  !  
  if (info == psb_success_) &
       & p%map = psb_linmap(psb_map_aggr_,desc_a,&
       & p%desc_ac,am2,am1,ilaggr,nlaggr)
  if (info == psb_success_) call am1%free()
  if (info == psb_success_) call am2%free()
  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='linmap build')
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

end subroutine mld_zaggrmat_nosmth_asb
