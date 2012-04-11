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
! File: mld_saggrmat_nosmth_asb.F90
!
! Subroutine: mld_saggrmat_nosmth_asb
! Version:    real
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
!  specified by the user through mld_sprecinit and mld_zprecset.
!
!  For details see
!    P. D'Ambra, D. di Serafino and  S. Filippone, On the development of
!    PSBLAS-based parallel two-level Schwarz preconditioners, Appl. Num. Math.,
!    57 (2007), 1181-1196.
!
!
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
subroutine mld_saggrmat_nosmth_asb(a,desc_a,ilaggr,nlaggr,p,info)
  use psb_base_mod
  use mld_s_inner_mod, mld_protect_name => mld_saggrmat_nosmth_asb

#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

  ! Arguments
  type(psb_sspmat_type), intent(in)          :: a
  type(psb_desc_type), intent(in)            :: desc_a
  integer, intent(inout)                     :: ilaggr(:), nlaggr(:)
  type(mld_s_onelev_type), intent(inout), target  :: p
  integer, intent(out)                       :: info
  type(psb_sspmat_type)  :: b, op_prol,op_restr

  ! Local variables
  integer :: ictxt,np,me, err_act
  integer(psb_mpik_) :: icomm, ndx, minfo
  character(len=20) :: name
  integer(psb_ipk_) :: ierr(5) 
  type(psb_s_coo_sparse_mat) :: acoo1, acoo2, bcoo, ac_coo, acoo
  type(psb_s_csr_sparse_mat) :: acsr1, acsr2
  integer            :: debug_level, debug_unit
  integer :: nrow, nglob, ncol, ntaggr, nzl, ip, &
       & naggr, nzt, naggrm1, i

  name='mld_aggrmat_nosmth_asb'
  if(psb_get_errstatus().ne.0) return 
  info=psb_success_
  call psb_erractionsave(err_act)


  ictxt = desc_a%get_context()
  icomm = desc_a%get_mpic()
  call psb_info(ictxt, me, np)
  nglob = desc_a%get_global_rows()
  nrow  = desc_a%get_local_rows()
  ncol  = desc_a%get_local_cols()


  naggr  = nlaggr(me+1)
  ntaggr = sum(nlaggr)
  naggrm1=sum(nlaggr(1:me))

  do i=1, nrow
    ilaggr(i) = ilaggr(i) + naggrm1
  end do
  call psb_halo(ilaggr,desc_a,info)

  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_halo')
    goto 9999
  end if

  call acoo1%allocate(ncol,ntaggr,ncol)

  do i=1,nrow
    acoo1%val(i) = sone
    acoo1%ia(i)  = i
    acoo1%ja(i)  = ilaggr(i)  
  end do

  call acoo1%set_dupl(psb_dupl_add_)
  call acoo1%set_nzeros(nrow)
  call acoo1%set_asb()
  call acoo1%fix(info)


  call op_prol%mv_from(acoo1)
  call op_prol%cscnv(info,type='csr',dupl=psb_dupl_add_)
  if (info == psb_success_)   call op_prol%transp(op_restr)

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


  call b%mv_from(bcoo) 

  if (p%parms%coarse_mat == mld_repl_mat_) then 

    call psb_cdall(ictxt,p%desc_ac,info,mg=ntaggr,repl=.true.)
    if (info == psb_success_) call psb_cdasb(p%desc_ac,info)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_cdall')
      goto 9999
    end if
    call psb_gather(p%ac,b,p%desc_ac,info,dupl=psb_dupl_add_,keeploc=.false.)
    if(info /= psb_success_) goto 9999        

  else if (p%parms%coarse_mat == mld_distr_mat_) then 

    nzl = b%get_nzeros()
    call b%mv_to(bcoo)

    if (info == psb_success_) call psb_cdall(ictxt,p%desc_ac,info,nl=nlaggr(me+1))
    if (info == psb_success_) call psb_cdins(nzl,bcoo%ia,bcoo%ja,p%desc_ac,info)
    if (info == psb_success_) call psb_cdasb(p%desc_ac,info)
    if (info == psb_success_) call psb_glob_to_loc(bcoo%ia(1:nzl),p%desc_ac,info,iact='I')
    if (info == psb_success_) call psb_glob_to_loc(bcoo%ja(1:nzl),p%desc_ac,info,iact='I')
    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Creating p%desc_ac and converting ac')
      goto 9999
    end if
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Assembld aux descr. distr.'
    call p%ac%mv_from(bcoo)

    call p%ac%set_nrows(p%desc_ac%get_local_rows())
    call p%ac%set_ncols(p%desc_ac%get_local_cols())
    call p%ac%set_asb()

    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_sp_free')
      goto 9999
    end if

    if (np>1) then 
      call op_prol%mv_to(acsr1)
      nzl = acsr1%get_nzeros()
      call psb_glob_to_loc(acsr1%ja(1:nzl),p%desc_ac,info,'I')
      if(info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_glob_to_loc')
        goto 9999
      end if
      call op_prol%mv_from(acsr1)
    endif
    call op_prol%set_ncols(p%desc_ac%get_local_cols())

    if (np>1) then 
      call op_restr%cscnv(info,type='coo',dupl=psb_dupl_add_)
      call op_restr%mv_to(acoo)
      nzl = acoo%get_nzeros()
      if (info == psb_success_) call psb_glob_to_loc(acoo%ia(1:nzl),p%desc_ac,info,'I')
      call acoo%set_dupl(psb_dupl_add_)
      if (info == psb_success_) call op_restr%mv_from(acoo)
      if (info == psb_success_) call op_restr%cscnv(info,type='csr')        
      if(info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,a_err='Converting op_restr to local')
        goto 9999
      end if
    end if
    call op_restr%set_nrows(p%desc_ac%get_local_cols())

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done ac '
  else
    info = psb_err_internal_error_
    call psb_errpush(psb_err_internal_error_,name,a_err='invalid mld_coarse_mat_')
    goto 9999
  end if


  call p%ac%cscnv(info,type='csr',dupl=psb_dupl_add_)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='cscnv')
    goto 9999
  end if

  !
  ! Copy the prolongation/restriction matrices into the descriptor map.
  !  op_restr => PR^T   i.e. restriction  operator
  !  op_prol => PR     i.e. prolongation operator
  !  
  if (info == psb_success_) &
       & p%map = psb_linmap(psb_map_aggr_,desc_a,&
       & p%desc_ac,op_restr,op_prol,ilaggr,nlaggr)
  if (info == psb_success_) call op_prol%free()
  if (info == psb_success_) call op_restr%free()
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

end subroutine mld_saggrmat_nosmth_asb
