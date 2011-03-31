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
! File: mld_sexample_ml.f90
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
!   b1 dd(u)  b2 dd(u)    b3 dd(u)    a1 d(u)   a2 d(u)  a3 d(u)  
! -   ------ -  ------ -  ------ -  -----  -  ------  -  ------ + a4 u  = 0
!      dxdx     dydy       dzdz        dx       dy         dz   
!
! with Dirichlet boundary conditions, on the unit cube  0<=x,y,z<=1.
!
! Example taken from:
!    C.T.Kelley
!    Iterative Methods for Linear and Nonlinear Equations
!    SIAM 1995
!
! In this sample program the index space of the discretized
! computational domain is first numbered sequentially in a standard way, 
! then the corresponding vector is distributed according to a BLOCK
! data distribution.
!
! Boundary conditions are set in a very simple way, by adding 
! equations of the form
!
!   u(x,y) = exp(-x^2-y^2-z^2)
!
! Note that if a1=a2=a3=a4=0., the PDE is the well-known Laplace equation.
!
program mld_sexample_ml
  use psb_base_mod
  use mld_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use data_input

  implicit none

  ! input parameters

  ! sparse matrices
  type(psb_sspmat_type) :: A

  ! sparse matrices descriptor
  type(psb_desc_type):: desc_A

  ! preconditioner
  type(mld_sprec_type)  :: P

  ! right-hand side, solution and residual vectors
  real(psb_spk_), allocatable , save  :: b(:), x(:), r(:)

  ! solver and preconditioner parameters
  real(psb_spk_)   :: tol, err
  integer          :: itmax, iter, istop
  integer          :: nlev

  ! parallel environment parameters
  integer            :: ictxt, iam, np

  ! other variables
  integer            :: choice       
  integer            :: i,info,j
  integer(psb_long_int_k_) :: amatsize, precsize, descsize
  integer            :: idim, ierr, ircode
  real(psb_dpk_)     :: t1, t2, tprec
  real(psb_spk_)     :: resmx, resmxp
  character(len=20)  :: name

  ! initialize the parallel environment

  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif
  !
  ! Hello world
  !
  if (iam == psb_root_) then 
    write(*,*) 'Welcome to MLD2P4 version: ',mld_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
  end if

  name='mld_sexample_ml'
  if(psb_get_errstatus() /= 0) goto 9999
  info=psb_success_
  call psb_set_errverbosity(2)

  ! get parameters

  call get_parms(ictxt,choice,idim,itmax,tol)

  !  allocate and fill in the coefficient matrix, rhs and initial guess

  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call create_matrix(idim,a,b,x,desc_a,ictxt,info)  
  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (iam == psb_root_) write(*,'("Overall matrix creation time : ",es12.5)')t2
  if (iam == psb_root_) write(*,'(" ")')

  select case(choice)

  case(1)

    ! initialize the default multi-level preconditioner, i.e. hybrid
    ! Schwarz, using RAS (with overlap 1 and ILU(0) on the blocks)
    ! as post-smoother and 4 block-Jacobi sweeps (with UMFPACK LU
    ! on the blocks) as distributed coarse-level solver

    call mld_precinit(P,'ML',info)

  case(2)

    ! set a three-level hybrid Schwarz preconditioner, which uses
    ! block Jacobi (with ILU(0) on the blocks) as post-smoother,
    ! a coarsest matrix replicated on the processors, and the
    ! LU factorization from UMFPACK as coarse-level solver

    call mld_precinit(P,'ML',info,nlev=3)
    call mld_precset(P,mld_smoother_type_,'BJAC',info)
    call mld_precset(P,mld_coarse_mat_,'REPL',info)
    call mld_precset(P,mld_coarse_solve_,'UMF',info)

  case(3)

    ! set a three-level additive Schwarz preconditioner, which uses
    ! RAS (with overlap 1 and ILU(0) on the blocks) as pre- and
    ! post-smoother, and 5 block-Jacobi sweeps (with UMFPACK LU
    ! on the blocks) as distributed coarsest-level solver

    call mld_precinit(P,'ML',info,nlev=3)
    call mld_precset(P,mld_ml_type_,'ADD',info)
    call mld_precset(P,mld_smoother_pos_,'TWOSIDE',info)
    call mld_precset(P,mld_coarse_sweeps_,5,info)

  end select

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

  ! set the solver parameters and the initial guess

  call psb_geall(x,desc_A,info)
  x(:) =0.0
  call psb_geasb(x,desc_A,info)

  ! solve Ax=b with preconditioned BiCGSTAB

  call psb_barrier(ictxt)
  t1 = psb_wtime()

  call psb_krylov('BICGSTAB',A,P,b,x,tol,desc_A,info,itmax,iter,err,itrace=1,istop=2)

  t2 = psb_wtime() - t1
  call psb_amx(ictxt,t2)

  call psb_geall(r,desc_A,info)
  r(:) =0.0
  call psb_geasb(r,desc_A,info)
  call psb_geaxpby(sone,b,szero,r,desc_A,info)
  call psb_spmm(-sone,A,x,sone,r,desc_A,info)
  call psb_genrm2s(resmx,r,desc_A,info)
  call psb_geamaxs(resmxp,r,desc_A,info)

  amatsize = psb_sizeof(A)
  descsize = psb_sizeof(desc_A)
  precsize = mld_sizeof(P)
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

9999 continue
  if(info /= psb_success_) then
    call psb_error(ictxt)
  end if
  call psb_exit(ictxt)
  stop

contains
  !
  ! get parameters from standard input
  !
  subroutine get_parms(ictxt,choice,idim,itmax,tol)

    use psb_base_mod
    implicit none

    integer             :: choice, idim, ictxt, itmax
    real(psb_spk_)      :: tol
    integer             :: iam, np

    call psb_info(ictxt,iam,np)

    if (iam == psb_root_) then
      ! read input parameters
      call read_data(choice,5)
      call read_data(idim,5)
      call read_data(itmax,5)
      call read_data(tol,5)
    end if

    call psb_bcast(ictxt,choice)
    call psb_bcast(ictxt,idim)
    call psb_bcast(ictxt,itmax)
    call psb_bcast(ictxt,tol)

  end subroutine get_parms

  !
  !  subroutine to allocate and fill in the coefficient matrix and
  !  the rhs 
  !
  subroutine create_matrix(idim,a,b,xv,desc_a,ictxt,info)
    !
    ! Discretize the partial diferential equation
    ! 
    !   b1 dd(u)  b2 dd(u)    b3 dd(u)    a1 d(u)   a2 d(u)  a3 d(u)  
    ! -   ------ -  ------ -  ------ -  -----  -  ------  -  ------ + a4 u = 0
    !      dxdx     dydy       dzdz        dx       dy         dz   
    !
    ! with Dirichlet boundary conditions, on the unit cube  0<=x,y,z<=1.
    !
    ! Boundary conditions are set in a very simple way, by adding 
    ! equations of the form
    !
    !   u(x,y) = exp(-x^2-y^2-z^2)
    !
    ! Note that if a1=a2=a3=a4=0., the PDE is the well-known Laplace equation.
    !
    use psb_base_mod
    implicit none
    integer                        :: idim
    integer, parameter             :: nb=20
    real(psb_spk_), allocatable    :: b(:),xv(:)
    type(psb_desc_type)            :: desc_a
    integer                        :: ictxt, info
    ! local variables
    type(psb_sspmat_type)    :: a
    real(psb_spk_)           :: zt(nb),glob_x,glob_y,glob_z
    integer                  :: m,n,nnz,glob_row,nlr,i,ii,ib,k,ipoints
    integer                  :: x,y,z,ia,indx_owner
    integer                  :: np, iam, nr, nt
    integer                  :: element
    integer, allocatable     :: irow(:),icol(:),myidx(:)
    real(psb_spk_), allocatable :: val(:)
    ! deltah dimension of each grid cell
    ! deltat discretization time
    real(psb_spk_)         :: deltah
    real(psb_spk_),parameter   :: rhs=0.e0,one=1.e0,zero=0.e0
    real(psb_dpk_)   :: t0, t1, t2, t3, tasb, talc, ttot, tgen 
    real(psb_spk_)   :: a1, a2, a3, a4, b1, b2, b3 
    external           :: a1, a2, a3, a4, b1, b2, b3
    integer            :: err_act

    character(len=20)  :: name

    info = psb_success_
    name = 'create_matrix'
    call psb_erractionsave(err_act)

    call psb_info(ictxt, iam, np)

    deltah = 1.d0/(idim-1)

    ! initialize array descriptor and sparse matrix storage; provide an
    ! estimate of the number of non zeroes 

    ipoints=idim-2
    m   = ipoints*ipoints*ipoints
    n   = m
    nnz = ((n*9)/(np))
    if(iam == psb_root_) write(0,'("Generating Matrix (size=",i0,")...")')n

    !
    ! Using a simple BLOCK distribution.
    !
    nt = (m+np-1)/np
    nr = max(0,min(nt,m-(iam*nt)))

    nt = nr
    call psb_sum(ictxt,nt) 
    if (nt /= m) write(0,*) iam, 'Initialization error ',nr,nt,m
    call psb_barrier(ictxt)
    t0 = psb_wtime()
    call psb_cdall(ictxt,desc_a,info,nl=nr)
    if (info == psb_success_) call psb_spall(a,desc_a,info,nnz=nnz)
    ! define  rhs from boundary conditions; also build initial guess 
    if (info == psb_success_) call psb_geall(b,desc_a,info)
    if (info == psb_success_) call psb_geall(xv,desc_a,info)
    nlr = psb_cd_get_local_rows(desc_a)
    call psb_barrier(ictxt)
    talc = psb_wtime()-t0

    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name)
      goto 9999
    end if

    ! we build an auxiliary matrix consisting of one row at a
    ! time; just a small matrix. might be extended to generate 
    ! a bunch of rows per call. 
    ! 
    allocate(val(20*nb),irow(20*nb),&
         &icol(20*nb),myidx(nlr),stat=info)
    if (info /= psb_success_ ) then 
      info=psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    endif

    do i=1,nlr
      myidx(i) = i
    end do


    call psb_loc_to_glob(myidx,desc_a,info)

    ! loop over rows belonging to current process in a block
    ! distribution.

    call psb_barrier(ictxt)
    t1 = psb_wtime()
    do ii=1, nlr,nb
      ib = min(nb,nlr-ii+1) 
      element = 1
      do k=1,ib
        i=ii+k-1
        ! local matrix pointer 
        glob_row=myidx(i)
        ! compute gridpoint coordinates
        if (mod(glob_row,ipoints*ipoints) == 0) then
          x = glob_row/(ipoints*ipoints)
        else
          x = glob_row/(ipoints*ipoints)+1
        endif
        if (mod((glob_row-(x-1)*ipoints*ipoints),ipoints) == 0) then
          y = (glob_row-(x-1)*ipoints*ipoints)/ipoints
        else
          y = (glob_row-(x-1)*ipoints*ipoints)/ipoints+1
        endif
        z = glob_row-(x-1)*ipoints*ipoints-(y-1)*ipoints
        ! glob_x, glob_y, glob_x coordinates
        glob_x=x*deltah
        glob_y=y*deltah
        glob_z=z*deltah

        ! check on boundary points 
        zt(k) = 0.d0
        ! internal point: build discretization
        !   
        !  term depending on   (x-1,y,z)
        !
        if (x == 1) then 
          val(element)=-b1(glob_x,glob_y,glob_z)&
               & -a1(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          zt(k) = exp(-glob_y**2-glob_z**2)*(-val(element))
        else
          val(element)=-b1(glob_x,glob_y,glob_z)&
               & -a1(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          icol(element) = (x-2)*ipoints*ipoints+(y-1)*ipoints+(z)
          irow(element) = glob_row
          element       = element+1
        endif
        !  term depending on     (x,y-1,z)
        if (y == 1) then 
          val(element)=-b2(glob_x,glob_y,glob_z)&
               & -a2(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          zt(k) = exp(-glob_x**2-glob_z**2)*(-val(element))
        else
          val(element)=-b2(glob_x,glob_y,glob_z)&
               & -a2(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          icol(element) = (x-1)*ipoints*ipoints+(y-2)*ipoints+(z)
          irow(element) = glob_row
          element       = element+1
        endif
        !  term depending on     (x,y,z-1)
        if (z == 1) then 
          val(element)=-b3(glob_x,glob_y,glob_z)&
               & -a3(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          zt(k) = exp(-glob_x**2-glob_y**2)*(-val(element))
        else
          val(element)=-b3(glob_x,glob_y,glob_z)&
               & -a3(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          icol(element) = (x-1)*ipoints*ipoints+(y-1)*ipoints+(z-1)
          irow(element) = glob_row
          element       = element+1
        endif
        !  term depending on     (x,y,z)
        val(element)=2*b1(glob_x,glob_y,glob_z)&
             & +2*b2(glob_x,glob_y,glob_z)&
             & +2*b3(glob_x,glob_y,glob_z)&
             & +a1(glob_x,glob_y,glob_z)&
             & +a2(glob_x,glob_y,glob_z)&
             & +a3(glob_x,glob_y,glob_z)
        val(element) = val(element)/(deltah*&
             & deltah)
        icol(element) = (x-1)*ipoints*ipoints+(y-1)*ipoints+(z)
        irow(element) = glob_row
        element       = element+1                  
        !  term depending on     (x,y,z+1)
        if (z == ipoints) then 
          val(element)=-b1(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          zt(k) = exp(-glob_x**2-glob_y**2)*exp(-glob_z)*(-val(element))  
        else
          val(element)=-b1(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          icol(element) = (x-1)*ipoints*ipoints+(y-1)*ipoints+(z+1)
          irow(element) = glob_row
          element       = element+1
        endif
        !  term depending on     (x,y+1,z)
        if (y == ipoints) then 
          val(element)=-b2(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          zt(k) = exp(-glob_x**2-glob_z**2)*exp(-glob_y)*(-val(element))  
        else
          val(element)=-b2(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          icol(element)=(x-1)*ipoints*ipoints+(y)*ipoints+(z)
          irow(element) = glob_row
          element       = element+1
        endif
        !  term depending on     (x+1,y,z)
        if (x == ipoints) then 
          val(element)=-b3(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          zt(k) = exp(-glob_y**2-glob_z**2)*exp(-glob_x)*(-val(element))  
        else
          val(element)=-b3(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          icol(element) = (x)*ipoints*ipoints+(y-1)*ipoints+(z)
          irow(element) = glob_row
          element       = element+1
        endif

      end do
      call psb_spins(element-1,irow,icol,val,a,desc_a,info)
      if(info /= psb_success_) exit
      call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),b,desc_a,info)
      if(info /= psb_success_) exit
      zt(:)=0.d0
      call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),xv,desc_a,info)
      if(info /= psb_success_) exit
    end do

    tgen = psb_wtime()-t1
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name)
      goto 9999
    end if

    deallocate(val,irow,icol)

    call psb_barrier(ictxt)
    t1 = psb_wtime()
    call psb_cdasb(desc_a,info)
    if (info == psb_success_) &
         & call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_)
    call psb_barrier(ictxt)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name)
      goto 9999
    end if
    call psb_geasb(b,desc_a,info)
    call psb_geasb(xv,desc_a,info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name)
      goto 9999
    end if
    tasb = psb_wtime()-t1
    call psb_barrier(ictxt)
    ttot = psb_wtime() - t0 

    call psb_amx(ictxt,talc)
    call psb_amx(ictxt,tgen)
    call psb_amx(ictxt,tasb)
    call psb_amx(ictxt,ttot)
    if(iam == psb_root_) then
      write(*,'("The matrix has been generated and assembled in ",a3," format.")')&
           &   a%get_fmt()
      write(*,'("-allocation  time : ",es12.5)') talc
      write(*,'("-coeff. gen. time : ",es12.5)') tgen
      write(*,'("-assembly    time : ",es12.5)') tasb
      write(*,'("-total       time : ",es12.5)') ttot

    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine create_matrix
end program mld_sexample_ml
!
! functions parametrizing the differential equation 
!  
function a1(x,y,z)
  use psb_base_mod, only : psb_spk_
  real(psb_spk_) :: a1
  real(psb_spk_) :: x,y,z
!  a1=1.e0
  a1=0.e0
end function a1
function a2(x,y,z)
  use psb_base_mod, only : psb_spk_
  real(psb_spk_) ::  a2
  real(psb_spk_) :: x,y,z
!  a2=2.e1*y
  a2=0.e0
end function a2
function a3(x,y,z)
  use psb_base_mod, only : psb_spk_
  real(psb_spk_) ::  a3
  real(psb_spk_) :: x,y,z      
!  a3=1.e0
  a3=0.e0
end function a3
function a4(x,y,z)
  use psb_base_mod, only : psb_spk_
  real(psb_spk_) ::  a4
  real(psb_spk_) :: x,y,z      
!  a4=1.e0
  a4=0.e0
end function a4
function b1(x,y,z)
  use psb_base_mod, only : psb_spk_
  real(psb_spk_) ::  b1   
  real(psb_spk_) :: x,y,z
  b1=1.e0
end function b1
function b2(x,y,z)
  use psb_base_mod, only : psb_spk_
  real(psb_spk_) ::  b2
  real(psb_spk_) :: x,y,z
  b2=1.e0
end function b2
function b3(x,y,z)
  use psb_base_mod, only : psb_spk_
  real(psb_spk_) ::  b3
  real(psb_spk_) :: x,y,z
  b3=1.e0
end function b3
