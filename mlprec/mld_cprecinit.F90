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
! File: mld_cprecinit.f90
!
! Subroutine: mld_cprecinit
! Version:    complex
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
!    'PJAC'           - point  Jacobi preconditioner
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
!                       with LU from SuperLU on the blocks, are applied at
!                       the coarsest level, on the distributed coarse matrix. 
!                       The smoothed aggregation algorithm with threshold 0
!                       is used to build the coarse matrix.
!
!  For the multilevel preconditioners, the levels are numbered in increasing
!  order starting from the finest one, i.e. level 1 is the finest level. 
!
!
! Arguments:
!    p       -  type(mld_cprec_type), input/output.
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
!               is assumed. If ptype /= 'ML', nlev is ignored.
!  
subroutine mld_cprecinit(p,ptype,info,nlev)

  use psb_sparse_mod
  use mld_c_prec_mod, mld_protect_name => mld_cprecinit
  use mld_c_jac_smoother
  use mld_c_as_smoother
  use mld_c_id_solver
  use mld_c_diag_solver
  use mld_c_ilu_solver
#if defined(HAVE_SLU_)
  use mld_c_slu_solver
#endif


  implicit none

! Arguments
  type(mld_cprec_type), intent(inout)    :: p
  character(len=*), intent(in)           :: ptype
  integer, intent(out)                   :: info
  integer, optional, intent(in)          :: nlev

! Local variables
  integer                                :: nlev_, ilev_
  character(len=*), parameter            :: name='mld_precinit'
  info = psb_success_
  
  if (allocated(p%precv)) then 
    call mld_precfree(p,info) 
    if (info /= psb_success_) then 
      ! Do we want to do something? 
    endif
  endif

  select case(psb_toupper(ptype(1:len_trim(ptype))))
  case ('NOPREC','NONE') 
    nlev_ = 1
    ilev_ = 1
    allocate(p%precv(nlev_),stat=info) 
    allocate(mld_c_base_smoother_type :: p%precv(ilev_)%sm, stat=info) 
    if (info /= psb_success_) return
    allocate(mld_c_id_solver_type :: p%precv(ilev_)%sm%sv, stat=info) 
    call p%precv(ilev_)%default()
    
  case ('JAC','DIAG','JACOBI') 
    nlev_ = 1
    ilev_ = 1
    allocate(p%precv(nlev_),stat=info) 
    allocate(mld_c_jac_smoother_type :: p%precv(ilev_)%sm, stat=info) 
    if (info /= psb_success_) return
    allocate(mld_c_diag_solver_type :: p%precv(ilev_)%sm%sv, stat=info) 
    call p%precv(ilev_)%default()

  case ('BJAC') 
    nlev_ = 1
    ilev_ = 1
    allocate(p%precv(nlev_),stat=info) 
    allocate(mld_c_jac_smoother_type :: p%precv(ilev_)%sm, stat=info) 
    if (info /= psb_success_) return
    allocate(mld_c_ilu_solver_type :: p%precv(ilev_)%sm%sv, stat=info) 
    call p%precv(ilev_)%default()
    
  case ('AS')
    nlev_ = 1
    ilev_ = 1
    allocate(p%precv(nlev_),stat=info) 
    allocate(mld_c_as_smoother_type :: p%precv(ilev_)%sm, stat=info) 
    if (info /= psb_success_) return
    allocate(mld_c_ilu_solver_type :: p%precv(ilev_)%sm%sv, stat=info) 
    call p%precv(ilev_)%default()


  case ('ML')
    
    if (present(nlev)) then 
      nlev_ = max(1,nlev)
    else
      nlev_ = 2
    end if
    ilev_ = 1
    allocate(p%precv(nlev_),stat=info) 
    allocate(mld_c_as_smoother_type :: p%precv(ilev_)%sm, stat=info) 
    if (info /= psb_success_) return
    allocate(mld_c_ilu_solver_type :: p%precv(ilev_)%sm%sv, stat=info) 
    call p%precv(ilev_)%default()


    if (nlev_ == 1) return 
    do ilev_ = 2, nlev_ -1 
      allocate(mld_c_as_smoother_type :: p%precv(ilev_)%sm, stat=info) 
      if (info /= psb_success_) return
      allocate(mld_c_ilu_solver_type :: p%precv(ilev_)%sm%sv, stat=info)       
      call p%precv(ilev_)%default()

    end do
    ilev_ = nlev_

    allocate(mld_c_jac_smoother_type :: p%precv(ilev_)%sm, stat=info) 
    if (info /= psb_success_) return
#if defined(HAVE_SLU_) 
    allocate(mld_c_slu_solver_type :: p%precv(ilev_)%sm%sv, stat=info)       
#else 
    allocate(mld_c_ilu_solver_type :: p%precv(ilev_)%sm%sv, stat=info)       
#endif
    call p%precv(ilev_)%default()
    call p%precv(ilev_)%set(mld_smoother_sweeps_,4,info)
    call p%precv(ilev_)%set(mld_sub_restr_,psb_none_,info)
    call p%precv(ilev_)%set(mld_sub_prol_,psb_none_,info)
    call p%precv(ilev_)%set(mld_sub_ovr_,0,info)

!!$    write(0,*) 'Check 5: ',allocated(p%precv(1)%sm)
      
  case default
    write(0,*) name,': Warning: Unknown preconditioner type request "',ptype,'"'
    info = psb_err_pivot_too_small_

  end select


end subroutine mld_cprecinit
