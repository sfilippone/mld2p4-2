!   
!   
!                             MLD2P4  version 2.2
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!    
!    (C) Copyright 2008-2018 
!  
!        Salvatore Filippone  
!        Pasqua D'Ambra   
!        Daniela di Serafino   
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the MLD2P4 group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!  
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
!    'DIAG', 'JACOBI' - diagonal/Jacobi           
!
!    'L1-DIAG', 'L1-JACOBI' - diagonal/Jacobi  with L1 norm correction
!
!    'GS', 'FBGS'     - Hybrid Gauss-Seidel, also symmetrized
!                       
!    'BJAC'           - block Jacobi preconditioner, with ILU(0)
!                       on the local blocks
!                       
!    'L1-BJAC'        - block Jacobi preconditioner, with ILU(0)
!                       on the local blocks and L1 correction for off-diag blocks
!
!    'AS'             - Additive Schwarz (AS), with
!                       overlap 1 and ILU(0) on the local submatrices
!
!    'ML'             - Multilevel hybrid preconditioner (additive on the
!                       same level and multiplicative through the levels),
!                       with 2 levels, pre  and post-smoothing, RAS with
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
!  
subroutine mld_dprecinit(ictxt,prec,ptype,info)

  use psb_base_mod
  use mld_d_prec_mod, mld_protect_name => mld_dprecinit
  use mld_d_jac_smoother
  use mld_d_as_smoother
  use mld_d_id_solver
  use mld_d_diag_solver
  use mld_d_ilu_solver
  use mld_d_gs_solver
#if defined(HAVE_UMF_)
  use mld_d_umf_solver
#endif
#if defined(HAVE_SLU_)
  use mld_d_slu_solver
#endif


  implicit none

  ! Arguments
  integer(psb_ipk_), intent(in)           :: ictxt
  class(mld_dprec_type), intent(inout)  :: prec
  character(len=*), intent(in)            :: ptype
  integer(psb_ipk_), intent(out)          :: info

  ! Local variables
  integer(psb_ipk_)                   :: nlev_, ilev_
  real(psb_dpk_)                      :: thr
  character(len=*), parameter         :: name='mld_precinit'
  info = psb_success_

  if (allocated(prec%precv)) then 
    call prec%free(info) 
    if (info /= psb_success_) then 
      ! Do we want to do something? 
    endif
  endif
  prec%ictxt = ictxt
  prec%ag_data%min_coarse_size = -1

  select case(psb_toupper(trim(ptype)))
  case ('NOPREC','NONE') 
    nlev_ = 1
    ilev_ = 1
    allocate(prec%precv(nlev_),stat=info) 
    allocate(mld_d_base_smoother_type :: prec%precv(ilev_)%sm, stat=info) 
    if (info /= psb_success_) return
    allocate(mld_d_id_solver_type :: prec%precv(ilev_)%sm%sv, stat=info) 
    call prec%precv(ilev_)%default()

  case ('JAC','DIAG','JACOBI') 
    nlev_ = 1
    ilev_ = 1
    allocate(prec%precv(nlev_),stat=info) 
    allocate(mld_d_jac_smoother_type :: prec%precv(ilev_)%sm, stat=info) 
    if (info /= psb_success_) return
    allocate(mld_d_diag_solver_type :: prec%precv(ilev_)%sm%sv, stat=info) 
    call prec%precv(ilev_)%default()

  case ('L1-DIAG','L1-JACOBI','L1_DIAG','L1_JACOBI') 
    nlev_ = 1
    ilev_ = 1
    allocate(prec%precv(nlev_),stat=info) 
    allocate(mld_d_jac_smoother_type :: prec%precv(ilev_)%sm, stat=info) 
    if (info /= psb_success_) return
    allocate(mld_d_l1_diag_solver_type :: prec%precv(ilev_)%sm%sv, stat=info) 
    call prec%precv(ilev_)%default()

  case ('GS','FWGS') 
    nlev_ = 1
    ilev_ = 1
    allocate(prec%precv(nlev_),stat=info) 
    allocate(mld_d_jac_smoother_type :: prec%precv(ilev_)%sm, stat=info) 
    if (info /= psb_success_) return
    allocate(mld_d_gs_solver_type :: prec%precv(ilev_)%sm%sv, stat=info) 
    call prec%precv(ilev_)%default()

  case ('BWGS') 
    nlev_ = 1
    ilev_ = 1
    allocate(prec%precv(nlev_),stat=info) 
    allocate(mld_d_jac_smoother_type :: prec%precv(ilev_)%sm, stat=info) 
    if (info /= psb_success_) return
    allocate(mld_d_bwgs_solver_type :: prec%precv(ilev_)%sm%sv, stat=info) 
    call prec%precv(ilev_)%default()

  case ('FBGS') 
    nlev_ = 1
    ilev_ = 1
    allocate(prec%precv(nlev_),stat=info)
    call prec%set('SMOOTHER_TYPE','FBGS',info)
    call prec%precv(ilev_)%default()

  case ('BJAC') 
    nlev_ = 1
    ilev_ = 1
    allocate(prec%precv(nlev_),stat=info) 
    allocate(mld_d_jac_smoother_type :: prec%precv(ilev_)%sm, stat=info) 
    if (info /= psb_success_) return
    allocate(mld_d_ilu_solver_type :: prec%precv(ilev_)%sm%sv, stat=info) 
    call prec%precv(ilev_)%default()

  case ('L1-BJAC','L1_BJAC') 
    nlev_ = 1
    ilev_ = 1
    allocate(prec%precv(nlev_),stat=info) 
    allocate(mld_d_l1_jac_smoother_type :: prec%precv(ilev_)%sm, stat=info) 
    if (info /= psb_success_) return
    allocate(mld_d_ilu_solver_type :: prec%precv(ilev_)%sm%sv, stat=info) 
    call prec%precv(ilev_)%default()

  case ('AS')
    nlev_ = 1
    ilev_ = 1
    allocate(prec%precv(nlev_),stat=info) 
    allocate(mld_d_as_smoother_type :: prec%precv(ilev_)%sm, stat=info) 
    if (info /= psb_success_) return
    allocate(mld_d_ilu_solver_type :: prec%precv(ilev_)%sm%sv, stat=info) 
    call prec%precv(ilev_)%default()


  case ('ML')

    nlev_ = prec%ag_data%max_levs
    ilev_ = 1
    allocate(prec%precv(nlev_),stat=info)

    do ilev_ = 1, nlev_  
      call prec%precv(ilev_)%default()
    end do
    call prec%set('ML_CYCLE','VCYCLE',info)
    call prec%set('SMOOTHER_TYPE','FBGS',info)
#if defined(HAVE_UMF_) 
    call prec%set('COARSE_SOLVE','UMF',info)
#elif defined(HAVE_MUMPS_)
    call prec%set('COARSE_SOLVE','MUMPS',info)    
#elif defined(HAVE_SLU_)
    call prec%set('COARSE_SOLVE','SLU',info)
#else
    call prec%set('COARSE_SOLVE','ILU',info)
#endif
    
  case default
    write(psb_err_unit,*) name,&
         &': Warning: Unknown preconditioner type request "',ptype,'"'
    info = psb_err_pivot_too_small_

  end select


end subroutine mld_dprecinit
