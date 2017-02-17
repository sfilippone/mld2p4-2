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
! File: mld_daggrmat_asb.f90
!
! Subroutine: mld_daggrmat_asb
! Version:    real
!
!  This routine builds a coarse-level matrix A_C from a fine-level matrix A
!  by using the Galerkin approach, i.e.
!
!                               A_C = P_C^T A P_C,
!
!  where P_C is a prolongator from the coarse level to the fine one.
! 
!  A mapping from the nodes of the adjacency graph of A to the nodes of the
!  adjacency graph of A_C has been computed by the mld_aggrmap_bld subroutine.
!  The prolongator P_C is built here from this mapping, according to the
!  value of p%iprcparm(mld_aggr_kind_), specified by the user through
!  mld_dprecinit and mld_zprecset.
!
!  Currently three different prolongators are implemented, corresponding to
!  three aggregation algorithms:
!  1. non-smoothed aggregation,
!  2. smoothed aggregation,
!  3. "bizarre" aggregation.
!  1. The non-smoothed aggregation uses as prolongator the piecewise constant
!     interpolation operator corresponding to the fine-to-coarse level mapping built
!     by mld_aggrmap_bld. This is called tentative prolongator.
!  2. The smoothed aggregation uses as prolongator the operator obtained by applying
!     a damped Jacobi smoother to the tentative prolongator.
!  3. The "bizarre" aggregation uses a prolongator proposed by the authors of MLD2P4.
!     This prolongator still requires a deep analysis and testing and its use is
!     not recommended.
!
!  For more details see
!    M. Brezina and P. Vanek, A black-box iterative solver based on a two-level
!    Schwarz method, Computing,  63 (1999), 233-263.
!    P. D'Ambra, D. di Serafino and S. Filippone, On the development of PSBLAS-based
!    parallel two-level Schwarz preconditioners, Appl. Num. Math., 57 (2007),
!    1181-1196.
!
!
!
! Arguments:
!    a          -  type(psb_cspmat_type), input.     
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
subroutine mld_daggrmat_asb(a,desc_a,ilaggr,nlaggr,p,info)

  use psb_base_mod
  use mld_d_inner_mod, mld_protect_name => mld_daggrmat_asb

  implicit none

! Arguments
  type(psb_dspmat_type), intent(in)              :: a
  type(psb_desc_type), intent(in)                  :: desc_a
  integer(psb_ipk_), intent(inout)                 :: ilaggr(:), nlaggr(:)
  type(mld_d_onelev_type), intent(inout), target :: p
  integer(psb_ipk_), intent(out)                   :: info

! Local variables
  type(psb_dspmat_type)      :: ac, op_prol,op_restr
  type(psb_d_coo_sparse_mat) :: acoo, bcoo
  type(psb_d_csr_sparse_mat) :: acsr1
  integer(psb_ipk_)            :: nzl,ntaggr, err_act
  integer(psb_ipk_)            :: debug_level, debug_unit
  integer(psb_ipk_)            :: ictxt,np,me
  character(len=20) :: name

  name='mld_aggrmat_asb'
  if(psb_get_errstatus().ne.0) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()


  ictxt = desc_a%get_context()

  call psb_info(ictxt, me, np)

  select case (p%parms%aggr_kind)
  case (mld_no_smooth_) 

    call mld_daggrmat_nosmth_asb(a,desc_a,ilaggr,nlaggr,&
         & p%parms,ac,op_prol,op_restr,info)

  case(mld_smooth_prol_) 

    call mld_daggrmat_smth_asb(a,desc_a,ilaggr,nlaggr, &
         & p%parms,ac,op_prol,op_restr,info)

  case(mld_biz_prol_) 

    call mld_daggrmat_biz_asb(a,desc_a,ilaggr,nlaggr, &
         & p%parms,ac,op_prol,op_restr,info)

  case(mld_min_energy_) 

    call mld_daggrmat_minnrg_asb(a,desc_a,ilaggr,nlaggr, &
         & p%parms,ac,op_prol,op_restr,info)

  case default
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='Invalid aggr kind')
    goto 9999

  end select
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Inner aggrmat asb')
    goto 9999
  end if



  ntaggr = sum(nlaggr)

  select case(p%parms%coarse_mat)

  case(mld_distr_mat_) 

    call ac%mv_to(bcoo)
    nzl = bcoo%get_nzeros()

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
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Converting op_restr to local')
        goto 9999
      end if
    end if
    call op_restr%set_nrows(p%desc_ac%get_local_cols())

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done ac '

  case(mld_repl_mat_) 
    !
    !
    call psb_cdall(ictxt,p%desc_ac,info,mg=ntaggr,repl=.true.)
    if (info == psb_success_) call psb_cdasb(p%desc_ac,info)
    if (info == psb_success_) &
         & call psb_gather(p%ac,ac,p%desc_ac,info,dupl=psb_dupl_add_,keeploc=.false.)

    if (info /= psb_success_) goto 9999

  case default 
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='invalid mld_coarse_mat_')
    goto 9999
  end select

  call p%ac%cscnv(info,type='csr',dupl=psb_dupl_add_)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='spcnv')
    goto 9999
  end if

  !
  ! Copy the prolongation/restriction matrices into the descriptor map.
  !  op_restr => PR^T   i.e. restriction  operator
  !  op_prol => PR     i.e. prolongation operator
  !  

  p%map = psb_linmap(psb_map_aggr_,desc_a,&
       & p%desc_ac,op_restr,op_prol,ilaggr,nlaggr)
  if (info == psb_success_) call op_prol%free()
  if (info == psb_success_) call op_restr%free()
  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='sp_Free')
    goto 9999
  end if


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine mld_daggrmat_asb
