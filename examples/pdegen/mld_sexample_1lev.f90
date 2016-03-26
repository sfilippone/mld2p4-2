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
! File: mld_sexample_1lev.f90
!
! This sample program solves a linear system obtained by discretizing a
! PDE with Dirichlet BCs. The solver is BiCGStab coupled  with one of the
! following multi-level preconditioner, as explained in Section 6.1 of
! the MLD2P4 User's and Reference Guide:
! - choice = 1, default multi-level Schwarz preconditioner (Sec. 6.1, Fig. 2)
! - choice = 2, hybrid three-level Schwarz preconditioner (Sec. 6.1, Fig. 3)
! - choice = 3, additive three-level Schwarz preconditioner (Sec. 6.1, Fig. 4)
!
! The PDE is a general second order equation in 3d
!
!   a1 dd(u)  a2 dd(u)    a3 dd(u)    b1 d(u)   b2 d(u)  b3 d(u)  
! -   ------ -  ------ -  ------ +  -----  +  ------  +  ------ + c u = f
!      dxdx     dydy       dzdz        dx       dy         dz   
!
! with Dirichlet boundary conditions
!   u = g 
!
!  on the unit cube  0<=x,y,z<=1.
!
!
! Note that if b1=b2=b3=c=0., the PDE is the  Laplace equation.
!
! In this sample program the index space of the discretized
! computational domain is first numbered sequentially in a standard way, 
! then the corresponding vector is distributed according to a BLOCK
! data distribution.
!
module spde_mod
contains
  !
  ! functions parametrizing the differential equation 
  !  
  function b1(x,y,z)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) :: b1
    real(psb_spk_), intent(in) :: x,y,z
    b1=1.e0/sqrt(3.e0)
  end function b1
  function b2(x,y,z)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  b2
    real(psb_spk_), intent(in) :: x,y,z
    b2=1.e0/sqrt(3.e0)
  end function b2
  function b3(x,y,z)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  b3
    real(psb_spk_), intent(in) :: x,y,z      
    b3=1.e0/sqrt(3.e0)
  end function b3
  function c(x,y,z)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  c
    real(psb_spk_), intent(in) :: x,y,z      
    c=0.e0
  end function c
  function a1(x,y,z)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  a1   
    real(psb_spk_), intent(in) :: x,y,z
    a1=1.e0/80
  end function a1
  function a2(x,y,z)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  a2
    real(psb_spk_), intent(in) :: x,y,z
    a2=1.e0/80
  end function a2
  function a3(x,y,z)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  a3
    real(psb_spk_), intent(in) :: x,y,z
    a3=1.e0/80
  end function a3
  function g(x,y,z)
    use psb_base_mod, only : psb_spk_, sone, szero
    real(psb_spk_) ::  g
    real(psb_spk_), intent(in) :: x,y,z
    g = szero
    if (x == sone) then
      g = sone
    else if (x == szero) then 
      g = exp(y**2-z**2)
    end if
  end function g
end module spde_mod

program mld_sexample_1lev
  use psb_base_mod
  use mld_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use data_input
  use spde_mod
  use mld_s_mumps_solver
  implicit none


  ! sparse matrices
  type(psb_sspmat_type) :: A

  ! descriptor of sparse matrices
  type(psb_desc_type):: desc_A

  ! preconditioner
  type(mld_sprec_type)  :: P

  ! right-hand side, solution and residual vectors
  type(psb_s_vect_type) :: x, b, r 

  ! solver parameters
  real(psb_spk_)   :: tol, err
  integer :: itmax, iter, itrace, istop
  type(mld_s_mumps_solver_type) :: sv

  ! parallel environment parameters
  integer            :: ictxt, iam, np

  ! other variables
  integer            :: i,info,j
  integer(psb_long_int_k_) :: amatsize, precsize, descsize
  integer            :: idim, nlev, ierr, ircode
  real(psb_dpk_) :: t1, t2, tprec
  real(psb_spk_) :: resmx, resmxp
  character(len=5)   :: afmt='CSR'
  character(len=20)  :: name

  ! initialize the parallel environment
  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif

  name='mld_sexample_ml'
  if(psb_get_errstatus() /= 0) goto 9999
  info=psb_success_
  call psb_set_errverbosity(2)
  !
  ! Hello world
  !
  if (iam == psb_root_) then 
    write(*,*) 'Welcome to MLD2P4 version: ',mld_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
  end if

  ! get parameters

  call get_parms(ictxt,idim,itmax,tol)

  ! allocate and fill in the coefficient matrix, rhs and initial guess

  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call psb_gen_pde3d(ictxt,idim,a,b,x,desc_a,afmt,&
       & a1,a2,a3,b1,b2,b3,c,g,info)  
  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (iam == psb_root_) write(*,'("Overall matrix creation time : ",es12.5)')t2
  if (iam == psb_root_) write(*,'(" ")')



  ! set MUMPS as solver
  call mld_precinit(P,'AS',info)
  call mld_precset(P,mld_sub_ovr_,2,info)
  call sv%default

  call P%set(sv,info)
  call mld_precset(P,mld_mumps_print_err_,10,info)



  ! build the preconditioner

  call psb_barrier(ictxt)
  t1 = psb_wtime()

  call mld_precbld(A,desc_A,P,info)

  tprec = psb_wtime()-t1
  call psb_amx(ictxt, tprec)

  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_precbld')
    goto 9999
  end if

  ! set the initial guess

  call psb_geall(x,desc_A,info)
  call x%set(szero)
  call psb_geasb(x,desc_A,info)

  ! solve Ax=b with preconditioned BiCGSTAB

  call psb_barrier(ictxt)
  t1 = psb_wtime()

  call psb_krylov('BICGSTAB',A,P,b,x,tol,desc_A,info,itmax,iter,err,itrace=1,istop=2)

  t2 = psb_wtime() - t1
  call psb_amx(ictxt,t2)

  call psb_geall(r,desc_A,info)
  call r%set(szero)
  call psb_geasb(r,desc_A,info)
  call psb_geaxpby(sone,b,szero,r,desc_A,info)
  call psb_spmm(-sone,A,x,sone,r,desc_A,info)
  resmx  = psb_genrm2(r,desc_A,info)
  resmxp = psb_geamax(r,desc_A,info)

  amatsize = a%sizeof()
  descsize = desc_a%sizeof()
  precsize = p%sizeof()
  call psb_sum(ictxt,amatsize)
  call psb_sum(ictxt,descsize)
  call psb_sum(ictxt,precsize)

  call mld_precdescr(P,info)

  if (iam == psb_root_) then
    write(*,'(" ")')
    write(*,'("Matrix from PDE example")')
    write(*,'("Computed solution on ",i8," processors")')np
    write(*,'("Iterations to convergence : ",i6)')iter
    write(*,'("Error estimate on exit    : ",es12.5)')err
    write(*,'("Time to build prec.       : ",es12.5)')tprec
    write(*,'("Time to solve system      : ",es12.5)')t2
    write(*,'("Time per iteration        : ",es12.5)')t2/(iter)
    write(*,'("Total time                : ",es12.5)')t2+tprec
    write(*,'("Residual 2-norm           : ",es12.5)')resmx
    write(*,'("Residual inf-norm         : ",es12.5)')resmxp
    write(*,'("Total memory occupation for A      : ",i12)')amatsize
    write(*,'("Total memory occupation for DESC_A : ",i12)')descsize
    write(*,'("Total memory occupation for PREC   : ",i12)')precsize
  end if

  call psb_gefree(b, desc_A,info)
  call psb_gefree(x, desc_A,info)
  call psb_spfree(A, desc_A,info)
  call mld_precfree(P,info)
  call psb_cdfree(desc_A,info)
  call psb_exit(ictxt)
  stop

9999 continue
  call psb_error(ictxt)

contains
  !
  ! get parameters from standard input
  !
  subroutine get_parms(ictxt,idim,itmax,tol)

    use psb_base_mod
    implicit none

    integer             :: idim, ictxt, itmax
    real(psb_spk_)      :: tol
    integer             :: iam, np

    call psb_info(ictxt,iam,np)

    if (iam == psb_root_) then
      ! read input parameters
      call read_data(idim,5)
      call read_data(itmax,5)
      call read_data(tol,5)
    end if

    call psb_bcast(ictxt,idim)
    call psb_bcast(ictxt,itmax)
    call psb_bcast(ictxt,tol)

  end subroutine get_parms
end program mld_sexample_1lev
