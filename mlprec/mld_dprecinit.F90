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
! File: mld_dprecinit.f90
!
! Subroutine: mld_dprecinit
! Version:    real
!
!  This routine allocates and initializes the preconditioner data structure,
!  according to the preconditioner type chosen by the user.
!  
!  A default preconditioner is set for each preconditioner type
!  specified by the user:
!
!    'NOPREC'         - no preconditioner
!
!    'DIAG'           - diagonal preconditioner
!
!    'BJAC'           - block Jacobi preconditioner, with ILU(0)
!                       on the local blocks
!
!    'AS'             - Restricted Additive Schwarz (RAS), with
!                       overlap 1 and ILU(0) on the local submatrices
!
!    'ML'             - Multilevel hybrid preconditioner (additive on the
!                       same level and multiplicative through the levels),
!                       with 2 levels and post-smoothing only. RAS with
!                       overlap 1 and ILU(0) on the local blocks is
!                       applied as post-smoother at each level, but the
!                       coarsest one; four sweeps of the block-Jacobi solver,
!                       with LU from UMFPACK on the blocks, are applied at
!                       the coarsest level, on the distributed coarse matrix. 
!                       The smoothed aggregation algorithm with threshold 0
!                       is used to build the coarse matrix.
!
!  For the multilevel preconditioners, the levels are numbered in increasing
!  order starting from the finest one, i.e. level 1 is the finest level. 
!
!
! Arguments:
!    p       -  type(mld_dprec_type), input/output.
!               The preconditioner data structure.
!    ptype   -  character(len=*), input.
!               The type of preconditioner. Its values are 'NOPREC',
!               'DIAG', 'BJAC', 'AS', 'ML' (and the corresponding
!               lowercase strings).
!    info    -  integer, output.
!               Error code.
!    nlev    -  integer, optional, input.
!               The number of levels of the multilevel preconditioner.
!               If nlev is not present and ptype='ML', then nlev=2
!               is assumed. If ptype/='ML', nlev is ignored.
!  
subroutine mld_dprecinit(p,ptype,info,nlev)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_dprecinit

  implicit none

! Arguments
  type(mld_dprec_type), intent(inout)    :: p
  character(len=*), intent(in)           :: ptype
  integer, intent(out)                   :: info
  integer, optional, intent(in)          :: nlev

! Local variables
  integer                                :: nlev_, ilev_
  character(len=*), parameter            :: name='mld_precinit'
  info = 0
  
  if (allocated(p%precv)) then 
    call mld_precfree(p,info) 
    if (info /=0) then 
      ! Do we want to do something? 
    endif
  endif

  select case(psb_toupper(ptype(1:len_trim(ptype))))
  case ('NOPREC') 
    nlev_ = 1
    ilev_ = 1
    allocate(p%precv(nlev_),stat=info) 
    if (info == 0) call psb_realloc(mld_ifpsz_,p%precv(ilev_)%iprcparm,info)
    if (info == 0) call psb_realloc(mld_rfpsz_,p%precv(ilev_)%rprcparm,info)
    if (info == 0) call psb_realloc(mld_ifpsz_,p%precv(ilev_)%prec%iprcparm,info)
    if (info == 0) call psb_realloc(mld_rfpsz_,p%precv(ilev_)%prec%rprcparm,info)
    if (info /= 0) return
    p%precv(ilev_)%iprcparm(:)      = 0
    p%precv(ilev_)%rprcparm(:)      = dzero
    p%precv(ilev_)%prec%iprcparm(:) = 0
    p%precv(ilev_)%prec%rprcparm(:) = dzero
    p%precv(ilev_)%prec%iprcparm(mld_smoother_type_)   = mld_noprec_
    p%precv(ilev_)%prec%iprcparm(mld_sub_solve_)       = mld_f_none_
    p%precv(ilev_)%prec%iprcparm(mld_sub_restr_)       = psb_none_
    p%precv(ilev_)%prec%iprcparm(mld_sub_prol_)        = psb_none_
    p%precv(ilev_)%prec%iprcparm(mld_sub_ren_)         = 0
    p%precv(ilev_)%prec%iprcparm(mld_sub_ovr_)         = 0
    p%precv(ilev_)%prec%iprcparm(mld_smoother_sweeps_) = 1

  case ('DIAG')
    nlev_ = 1
    ilev_ = 1
    allocate(p%precv(nlev_),stat=info) 
    if (info == 0) call psb_realloc(mld_ifpsz_,p%precv(ilev_)%iprcparm,info)
    if (info == 0) call psb_realloc(mld_rfpsz_,p%precv(ilev_)%rprcparm,info)
    if (info == 0) call psb_realloc(mld_ifpsz_,p%precv(ilev_)%prec%iprcparm,info)
    if (info == 0) call psb_realloc(mld_rfpsz_,p%precv(ilev_)%prec%rprcparm,info)
    if (info /= 0) return
    p%precv(ilev_)%iprcparm(:)      = 0
    p%precv(ilev_)%rprcparm(:)      = dzero
    p%precv(ilev_)%prec%iprcparm(:) = 0
    p%precv(ilev_)%prec%rprcparm(:) = dzero
    p%precv(ilev_)%prec%iprcparm(mld_smoother_type_)   = mld_diag_
    p%precv(ilev_)%prec%iprcparm(mld_sub_solve_)       = mld_f_none_
    p%precv(ilev_)%prec%iprcparm(mld_sub_restr_)       = psb_none_
    p%precv(ilev_)%prec%iprcparm(mld_sub_prol_)        = psb_none_
    p%precv(ilev_)%prec%iprcparm(mld_sub_ren_)         = 0 
    p%precv(ilev_)%prec%iprcparm(mld_sub_ovr_)         = 0
    p%precv(ilev_)%prec%iprcparm(mld_smoother_sweeps_) = 1

  case ('BJAC') 
    nlev_ = 1
    ilev_ = 1
    allocate(p%precv(nlev_),stat=info) 
    if (info == 0) call psb_realloc(mld_ifpsz_,p%precv(ilev_)%iprcparm,info)
    if (info == 0) call psb_realloc(mld_rfpsz_,p%precv(ilev_)%rprcparm,info)
    if (info == 0) call psb_realloc(mld_ifpsz_,p%precv(ilev_)%prec%iprcparm,info)
    if (info == 0) call psb_realloc(mld_rfpsz_,p%precv(ilev_)%prec%rprcparm,info)
    if (info /= 0) return
    p%precv(ilev_)%iprcparm(:)      = 0
    p%precv(ilev_)%rprcparm(:)      = dzero
    p%precv(ilev_)%prec%iprcparm(:) = 0
    p%precv(ilev_)%prec%rprcparm(:) = dzero
    p%precv(ilev_)%prec%iprcparm(mld_smoother_type_)   = mld_bjac_
    p%precv(ilev_)%prec%iprcparm(mld_sub_solve_)       = mld_ilu_n_
    p%precv(ilev_)%prec%iprcparm(mld_sub_restr_)       = psb_none_
    p%precv(ilev_)%prec%iprcparm(mld_sub_prol_)        = psb_none_
    p%precv(ilev_)%prec%iprcparm(mld_sub_ren_)         = 0
    p%precv(ilev_)%prec%iprcparm(mld_sub_ovr_)         = 0
    p%precv(ilev_)%prec%iprcparm(mld_sub_fillin_)      = 0
    p%precv(ilev_)%prec%iprcparm(mld_smoother_sweeps_) = 1

  case ('AS')
    nlev_ = 1
    ilev_ = 1
    allocate(p%precv(nlev_),stat=info) 
    if (info == 0) call psb_realloc(mld_ifpsz_,p%precv(ilev_)%iprcparm,info)
    if (info == 0) call psb_realloc(mld_rfpsz_,p%precv(ilev_)%rprcparm,info)
    if (info == 0) call psb_realloc(mld_ifpsz_,p%precv(ilev_)%prec%iprcparm,info)
    if (info == 0) call psb_realloc(mld_rfpsz_,p%precv(ilev_)%prec%rprcparm,info)
    if (info /= 0) return
    p%precv(ilev_)%iprcparm(:)      = 0
    p%precv(ilev_)%rprcparm(:)      = dzero
    p%precv(ilev_)%prec%iprcparm(:) = 0
    p%precv(ilev_)%prec%rprcparm(:) = dzero
    p%precv(ilev_)%prec%iprcparm(mld_smoother_type_)   = mld_as_ 
    p%precv(ilev_)%prec%iprcparm(mld_sub_solve_)       = mld_ilu_n_
    p%precv(ilev_)%prec%iprcparm(mld_sub_restr_)       = psb_halo_
    p%precv(ilev_)%prec%iprcparm(mld_sub_prol_)        = psb_none_
    p%precv(ilev_)%prec%iprcparm(mld_sub_ren_)         = 0
    p%precv(ilev_)%prec%iprcparm(mld_sub_ovr_)         = 1
    p%precv(ilev_)%prec%iprcparm(mld_sub_fillin_)      = 0
    p%precv(ilev_)%prec%iprcparm(mld_smoother_sweeps_) = 1


  case ('ML')
    
    if (present(nlev)) then 
      nlev_ = max(1,nlev)
    else
      nlev_ = 2
    end if
    ilev_ = 1
    allocate(p%precv(nlev_),stat=info) 
    if (info == 0) call psb_realloc(mld_ifpsz_,p%precv(ilev_)%iprcparm,info)
    if (info == 0) call psb_realloc(mld_rfpsz_,p%precv(ilev_)%rprcparm,info)
    if (info == 0) call psb_realloc(mld_ifpsz_,p%precv(ilev_)%prec%iprcparm,info)
    if (info == 0) call psb_realloc(mld_rfpsz_,p%precv(ilev_)%prec%rprcparm,info)
    if (info /= 0) return
    p%precv(ilev_)%iprcparm(:)      = 0
    p%precv(ilev_)%rprcparm(:)      = dzero
    p%precv(ilev_)%prec%iprcparm(:) = 0
    p%precv(ilev_)%prec%rprcparm(:) = dzero
    p%precv(ilev_)%prec%iprcparm(mld_smoother_type_)   = mld_as_ 
    p%precv(ilev_)%prec%iprcparm(mld_sub_solve_)       = mld_ilu_n_
    p%precv(ilev_)%prec%iprcparm(mld_sub_restr_)       = psb_halo_
    p%precv(ilev_)%prec%iprcparm(mld_sub_prol_)        = psb_none_
    p%precv(ilev_)%prec%iprcparm(mld_sub_ren_)         = 0
    p%precv(ilev_)%prec%iprcparm(mld_sub_ovr_)         = 1
    p%precv(ilev_)%prec%iprcparm(mld_sub_fillin_)      = 0
    p%precv(ilev_)%prec%iprcparm(mld_smoother_sweeps_) = 1
    if (nlev_ == 1) return 

    do ilev_ = 2, nlev_ -1 
      if (info == 0) call psb_realloc(mld_ifpsz_,p%precv(ilev_)%iprcparm,info)
      if (info == 0) call psb_realloc(mld_rfpsz_,p%precv(ilev_)%rprcparm,info)
      if (info == 0) call psb_realloc(mld_ifpsz_,p%precv(ilev_)%prec%iprcparm,info)
      if (info == 0) call psb_realloc(mld_rfpsz_,p%precv(ilev_)%prec%rprcparm,info)
      if (info /= 0) return
      p%precv(ilev_)%iprcparm(:)      = 0
      p%precv(ilev_)%rprcparm(:)      = dzero
      p%precv(ilev_)%prec%iprcparm(:) = 0
      p%precv(ilev_)%prec%rprcparm(:) = dzero
      p%precv(ilev_)%prec%iprcparm(mld_smoother_type_)   = mld_as_
      p%precv(ilev_)%prec%iprcparm(mld_sub_restr_)       = psb_halo_
      p%precv(ilev_)%prec%iprcparm(mld_sub_prol_)        = psb_none_
      p%precv(ilev_)%prec%iprcparm(mld_sub_ren_)         = 0
      p%precv(ilev_)%prec%iprcparm(mld_sub_ovr_)         = 1
      p%precv(ilev_)%prec%iprcparm(mld_sub_solve_)       = mld_ilu_n_
      p%precv(ilev_)%prec%iprcparm(mld_sub_fillin_)      = 0
      p%precv(ilev_)%prec%iprcparm(mld_smoother_sweeps_) = 1
      p%precv(ilev_)%iprcparm(mld_ml_type_)         = mld_mult_ml_
      p%precv(ilev_)%iprcparm(mld_aggr_alg_)        = mld_dec_aggr_
      p%precv(ilev_)%iprcparm(mld_aggr_kind_)       = mld_smooth_prol_
      p%precv(ilev_)%iprcparm(mld_coarse_mat_)      = mld_distr_mat_
      p%precv(ilev_)%iprcparm(mld_smoother_pos_)    = mld_post_smooth_
      p%precv(ilev_)%iprcparm(mld_aggr_omega_alg_)  = mld_eig_est_
      p%precv(ilev_)%iprcparm(mld_aggr_eig_)        = mld_max_norm_
      p%precv(ilev_)%rprcparm(mld_aggr_omega_val_)  = dzero
      p%precv(ilev_)%rprcparm(mld_aggr_thresh_)     = dzero
    end do
    ilev_ = nlev_
    if (info == 0) call psb_realloc(mld_ifpsz_,p%precv(ilev_)%iprcparm,info)
    if (info == 0) call psb_realloc(mld_rfpsz_,p%precv(ilev_)%rprcparm,info)
    if (info == 0) call psb_realloc(mld_ifpsz_,p%precv(ilev_)%prec%iprcparm,info)
    if (info == 0) call psb_realloc(mld_rfpsz_,p%precv(ilev_)%prec%rprcparm,info)
    if (info /= 0) return
    p%precv(ilev_)%iprcparm(:)      = 0
    p%precv(ilev_)%rprcparm(:)      = dzero
    p%precv(ilev_)%prec%iprcparm(:) = 0
    p%precv(ilev_)%prec%rprcparm(:) = dzero
    p%precv(ilev_)%prec%iprcparm(mld_coarse_solve_)    = mld_bjac_
#if defined(HAVE_UMF_)
    p%precv(ilev_)%prec%iprcparm(mld_sub_solve_)       = mld_umf_
#elif defined(HAVE_SLU_)
    p%precv(ilev_)%prec%iprcparm(mld_sub_solve_)       = mld_slu_
#else 
    p%precv(ilev_)%prec%iprcparm(mld_sub_solve_)       = mld_ilu_n_
#endif
    p%precv(ilev_)%prec%iprcparm(mld_smoother_type_)   = mld_bjac_
    p%precv(ilev_)%prec%iprcparm(mld_coarse_mat_)      = mld_distr_mat_
    p%precv(ilev_)%prec%iprcparm(mld_sub_restr_)       = psb_none_
    p%precv(ilev_)%prec%iprcparm(mld_sub_prol_)        = psb_none_
    p%precv(ilev_)%prec%iprcparm(mld_sub_ren_)         = 0
    p%precv(ilev_)%prec%iprcparm(mld_sub_ovr_)         = 0
    p%precv(ilev_)%prec%iprcparm(mld_sub_fillin_)      = 0
    p%precv(ilev_)%prec%iprcparm(mld_smoother_sweeps_) = 4
    p%precv(ilev_)%iprcparm(mld_ml_type_)         = mld_mult_ml_
    p%precv(ilev_)%iprcparm(mld_aggr_alg_)        = mld_dec_aggr_
    p%precv(ilev_)%iprcparm(mld_aggr_kind_)       = mld_smooth_prol_
    p%precv(ilev_)%iprcparm(mld_smoother_pos_)    = mld_post_smooth_
    p%precv(ilev_)%iprcparm(mld_aggr_omega_alg_)  = mld_eig_est_
    p%precv(ilev_)%iprcparm(mld_aggr_eig_)        = mld_max_norm_
    p%precv(ilev_)%rprcparm(mld_aggr_omega_val_)  = dzero
    p%precv(ilev_)%rprcparm(mld_aggr_thresh_)     = dzero
      
  case default
    write(0,*) name,': Warning: Unknown preconditioner type request "',ptype,'"'
    info = 2

  end select


end subroutine mld_dprecinit
