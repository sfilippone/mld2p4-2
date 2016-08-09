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
! File: mld_dcoarse_bld.f90
!
! Subroutine: mld_dcoarse_bld
! Version:    real
!
!  This routine builds the matrix associated to the current level of the
!  multilevel preconditioner from the matrix associated to the previous level,
!  by using a smoothed aggregation technique (therefore, it also builds the
!  prolongation and restriction operators mapping the current level to the
!  previous one and vice versa). Then the routine builds the base preconditioner
!  at the current level.
!  The current level is regarded as the coarse one, while the previous as
!  the fine one. This is in agreement with the fact that the routine is called,
!  by mld_mlprec_bld, only on levels >=2.
!
! 
! Arguments:
!    a       -  type(psb_dspmat_type).
!               The sparse matrix structure containing the local part of the
!               fine-level matrix.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    p       -  type(mld_d_onelev_type), input/output.
!               The 'one-level' data structure containing the local part
!               of the base preconditioner to be built as well as
!               information concerning the prolongator and its transpose.
!    info    -  integer, output.
!               Error code.         
!  
subroutine mld_dcoarse_bld(a,desc_a,p,info)

  use psb_base_mod
  use mld_d_inner_mod, mld_protect_name => mld_dcoarse_bld

  implicit none

  ! Arguments
  type(psb_dspmat_type), intent(in), target     :: a
  type(psb_desc_type), intent(in), target       :: desc_a
  type(mld_d_onelev_type), intent(inout),target :: p
  integer(psb_ipk_), intent(out)                :: info

  ! Local variables
  character(len=20)                :: name
  integer(psb_mpik_)               :: ictxt, np, me
  integer(psb_ipk_)                :: err_act
  integer(psb_ipk_), allocatable   :: ilaggr(:), nlaggr(:)
  type(psb_dspmat_type)            :: ac, op_prol,op_restr
  type(psb_d_coo_sparse_mat)       :: acoo, bcoo
  type(psb_d_csr_sparse_mat)       :: acsr1
  integer(psb_ipk_)                :: nzl, ntaggr
  integer(psb_ipk_)            :: debug_level, debug_unit

  name='mld_dcoarse_bld'
  if (psb_get_errstatus().ne.0) return 
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  info  = psb_success_
  ictxt = desc_a%get_context()
  call psb_info(ictxt,me,np)

  call mld_check_def(p%parms%ml_type,'Multilevel type',&
       &   mld_mult_ml_,is_legal_ml_type)
  call mld_check_def(p%parms%aggr_alg,'Aggregation',&
       &   mld_dec_aggr_,is_legal_ml_aggr_alg)
  call mld_check_def(p%parms%aggr_ord,'Ordering',&
       &   mld_aggr_ord_nat_,is_legal_ml_aggr_ord)
  call mld_check_def(p%parms%aggr_kind,'Smoother',&
       &   mld_smooth_prol_,is_legal_ml_aggr_kind)
  call mld_check_def(p%parms%coarse_mat,'Coarse matrix',&
       &   mld_distr_mat_,is_legal_ml_coarse_mat)
  call mld_check_def(p%parms%aggr_filter,'Use filtered matrix',&
       &   mld_no_filter_mat_,is_legal_aggr_filter)
  call mld_check_def(p%parms%smoother_pos,'smooth_pos',&
       &   mld_pre_smooth_,is_legal_ml_smooth_pos)
  call mld_check_def(p%parms%aggr_omega_alg,'Omega Alg.',&
       &   mld_eig_est_,is_legal_ml_aggr_omega_alg)
  call mld_check_def(p%parms%aggr_eig,'Eigenvalue estimate',&
       &   mld_max_norm_,is_legal_ml_aggr_eig)
  call mld_check_def(p%parms%aggr_omega_val,'Omega',dzero,is_legal_d_omega)
  call mld_check_def(p%parms%aggr_thresh,'Aggr_Thresh',dzero,is_legal_d_aggr_thrs)

  select case(p%parms%aggr_alg)
  case (mld_dec_aggr_, mld_sym_dec_aggr_)  
    
    !
    !  Build a mapping between the row indices of the fine-level matrix 
    !  and the row indices of the coarse-level matrix, according to a decoupled 
    !  aggregation algorithm. This also defines a tentative prolongator from
    !  the coarse to the fine level.
    ! 
    call mld_aggrmap_bld(p%parms%aggr_alg,p%parms%aggr_ord,p%parms%aggr_thresh,&
         & a,desc_a,ilaggr,nlaggr,info)
    
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_aggrmap_bld')
      goto 9999
    end if
    
    !
    ! Build the coarse-level matrix from the fine-level one, starting from 
    ! the mapping defined by mld_aggrmap_bld and applying the aggregation
    ! algorithm specified by p%iprcparm(mld_aggr_kind_)
    !
    call mld_daggrmat_asb(a,desc_a,ilaggr,nlaggr,p%parms,ac,op_prol,op_restr,info)
    
    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_aggrmat_asb')
      goto 9999
    end if

  case (mld_bcmatch_aggr_)
    write(0,*) 'Matching is not implemented yet '
    info = -1111
    call psb_errpush(psb_err_input_value_invalid_i_,name,&
         & i_err=(/ione,p%parms%aggr_alg,izero,izero,izero/))
    goto 9999
    
  case default

    info = -1
    call psb_errpush(psb_err_input_value_invalid_i_,name,&
         & i_err=(/ione,p%parms%aggr_alg,izero,izero,izero/))
    goto 9999

  end select
  

  ! Common code refactored here.
  
  ntaggr = sum(nlaggr)

  select case(p%parms%coarse_mat)

  case(mld_distr_mat_) 

    call ac%mv_to(bcoo)
    if (p%parms%clean_zeros) call bcoo%clean_zeros(info)
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
    !
    ! Clip to local rows.
    !
    call op_restr%set_nrows(p%desc_ac%get_local_rows())

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done ac '

  case(mld_repl_mat_) 
    !
    !
    call psb_cdall(ictxt,p%desc_ac,info,mg=ntaggr,repl=.true.)
    if (info == psb_success_) call psb_cdasb(p%desc_ac,info)
    if ((info == psb_success_).and.p%parms%clean_zeros) call ac%clean_zeros(info)
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
  !
  ! Fix the base_a and base_desc pointers for handling of residuals.
  ! This is correct because this routine is only called at levels >=2.
  !
  p%base_a    => p%ac
  p%base_desc => p%desc_ac
  
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine mld_dcoarse_bld
