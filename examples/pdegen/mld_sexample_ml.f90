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
  integer            :: i,info,j,amatsize,descsize,precsize
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

  name='mld_sexample_ml'
  if(psb_get_errstatus() /= 0) goto 9999
  info=0
  call psb_set_errverbosity(2)

! get parameters

  call get_parms(ictxt,choice,idim,itmax,tol)

!  allocate and fill in the coefficient matrix, rhs and initial guess

  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call create_matrix(idim,A,b,x,desc_A,part_block,ictxt,info)
  t2 = psb_wtime() - t1
  if(info /= 0) then
    info=4010
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_amx(ictxt,t2)
  if (iam == psb_root_) write(*,'("Overall matrix creation time : ",es10.4)')t2
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

  if (info /= 0) then
    call psb_errpush(4010,name,a_err='psb_precbld')
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

  if (iam==psb_root_) then
    write(*,'(" ")')
    write(*,'("Matrix from PDE example")')
    write(*,'("Computed solution on ",i8," processors")')np
    write(*,'("Iterations to convergence : ",i6)')iter
    write(*,'("Error estimate on exit    : ",es10.4)')err
    write(*,'("Time to build prec.       : ",es10.4)')tprec
    write(*,'("Time to solve system      : ",es10.4)')t2
    write(*,'("Time per iteration        : ",es10.4)')t2/(iter)
    write(*,'("Total time                : ",es10.4)')t2+tprec
    write(*,'("Residual 2-norm           : ",es10.4)')resmx
    write(*,'("Residual inf-norm         : ",es10.4)')resmxp
    write(*,'("Total memory occupation for A      : ",i10)')amatsize
    write(*,'("Total memory occupation for DESC_A : ",i10)')descsize
    write(*,'("Total memory occupation for PREC   : ",i10)')precsize
  end if

  call psb_gefree(b, desc_A,info)
  call psb_gefree(x, desc_A,info)
  call psb_spfree(A, desc_A,info)
  call mld_precfree(P,info)
  call psb_cdfree(desc_A,info)

9999 continue
  if(info /= 0) then
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

    if (iam==psb_root_) then
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
  subroutine create_matrix(idim,a,b,xv,desc_a,parts,ictxt,info)
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
    integer, parameter             :: nbmax=10
    real(psb_spk_), allocatable  :: b(:),xv(:)
    type(psb_desc_type)            :: desc_a
    integer                        :: ictxt, info
    interface 
      !   .....user passed subroutine.....
      subroutine parts(global_indx,n,np,pv,nv)
        implicit none
        integer, intent(in)  :: global_indx, n, np
        integer, intent(out) :: nv
        integer, intent(out) :: pv(*) 
      end subroutine parts
    end interface
   ! local variables
    type(psb_sspmat_type)    :: a
    real(psb_spk_)         :: zt(nbmax),glob_x,glob_y,glob_z
    integer                  :: m,n,nnz,glob_row
    integer                  :: x,y,z,ia,indx_owner
    integer                  :: np, iam
    integer                  :: element
    integer                  :: nv, inv
    integer, allocatable     :: irow(:),icol(:)
    real(psb_spk_), allocatable :: val(:)
    integer, allocatable     :: prv(:)
    ! deltah dimension of each grid cell
    ! deltat discretization time
    real(psb_spk_)         :: deltah
    real(psb_spk_),parameter   :: rhs=0.e0,one=1.e0,zero=0.e0
    real(psb_dpk_)   :: t1, t2, t3, tins, tasb
    real(psb_spk_)   :: a1, a2, a3, a4, b1, b2, b3 
    external           :: a1, a2, a3, a4, b1, b2, b3
    integer            :: err_act

    character(len=20)  :: name

    info = 0
    name = 'create_matrix'
    call psb_erractionsave(err_act)

    call psb_info(ictxt, iam, np)

    deltah = 1.e0/(idim-1)

    ! initialize array descriptor and sparse matrix storage; provide an
    ! estimate of the number of non zeroes 

    m   = idim*idim*idim
    n   = m
    nnz = ((n*9)/(np))
    if(iam == psb_root_) write(0,'("Generating Matrix (size=",i0x,")...")')n

    call psb_cdall(ictxt,desc_a,info,mg=n,parts=parts)
    call psb_spall(a,desc_a,info,nnz=nnz)
    ! define  rhs from boundary conditions; also build initial guess 
    call psb_geall(b,desc_a,info)
    call psb_geall(xv,desc_a,info)
    if(info /= 0) then
      info=4010
      call psb_errpush(info,name)
      goto 9999
    end if

    ! we build an auxiliary matrix consisting of one row at a
    ! time; just a small matrix. might be extended to generate 
    ! a bunch of rows per call. 
    ! 
    allocate(val(20*nbmax),irow(20*nbmax),&
         &icol(20*nbmax),prv(np),stat=info)
    if (info /= 0 ) then 
      info=4000
      call psb_errpush(info,name)
      goto 9999
    endif

    tins = 0.e0
    call psb_barrier(ictxt)
    t1 = psb_wtime()

    ! loop over rows belonging to current process in a block
    ! distribution.

    !    icol(1)=1    
    do glob_row = 1, n
      call parts(glob_row,n,np,prv,nv)
      do inv = 1, nv
        indx_owner = prv(inv)
        if (indx_owner == iam) then
          ! local matrix pointer 
          element=1
          ! compute gridpoint coordinates
          if (mod(glob_row,(idim*idim)) == 0) then
            x = glob_row/(idim*idim)
          else
            x = glob_row/(idim*idim)+1
          endif
          if (mod((glob_row-(x-1)*idim*idim),idim) == 0) then
            y = (glob_row-(x-1)*idim*idim)/idim
          else
            y = (glob_row-(x-1)*idim*idim)/idim+1
          endif
          z = glob_row-(x-1)*idim*idim-(y-1)*idim
          ! glob_x, glob_y, glob_x coordinates
          glob_x=x*deltah
          glob_y=y*deltah
          glob_z=z*deltah

          ! check on boundary points 
          zt(1) = 0.e0
          ! internal point: build discretization
          !   
          !  term depending on   (x-1,y,z)
          !
          if (x==1) then 
            val(element)=-b1(glob_x,glob_y,glob_z)&
                 & -a1(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            zt(1) = exp(-glob_y**2-glob_z**2)*(-val(element))
          else
            val(element)=-b1(glob_x,glob_y,glob_z)&
                 & -a1(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            icol(element)=(x-2)*idim*idim+(y-1)*idim+(z)
            element=element+1
          endif
          !  term depending on     (x,y-1,z)
          if (y==1) then 
            val(element)=-b2(glob_x,glob_y,glob_z)&
                 & -a2(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            zt(1) = exp(-glob_y**2-glob_z**2)*exp(-glob_x)*(-val(element))  
          else
            val(element)=-b2(glob_x,glob_y,glob_z)&
                 & -a2(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            icol(element)=(x-1)*idim*idim+(y-2)*idim+(z)
            element=element+1
          endif
          !  term depending on     (x,y,z-1)
          if (z==1) then 
            val(element)=-b3(glob_x,glob_y,glob_z)&
                 & -a3(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            zt(1) = exp(-glob_y**2-glob_z**2)*exp(-glob_x)*(-val(element))  
          else
            val(element)=-b3(glob_x,glob_y,glob_z)&
                 & -a3(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            icol(element)=(x-1)*idim*idim+(y-1)*idim+(z-1)
            element=element+1
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
          icol(element)=(x-1)*idim*idim+(y-1)*idim+(z)
          element=element+1                  
          !  term depending on     (x,y,z+1)
          if (z==idim) then 
            val(element)=-b1(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            zt(1) = exp(-glob_y**2-glob_z**2)*exp(-glob_x)*(-val(element))  
          else
            val(element)=-b1(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            icol(element)=(x-1)*idim*idim+(y-1)*idim+(z+1)
            element=element+1
          endif
          !  term depending on     (x,y+1,z)
          if (y==idim) then 
            val(element)=-b2(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            zt(1) = exp(-glob_y**2-glob_z**2)*exp(-glob_x)*(-val(element))  
          else
            val(element)=-b2(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            icol(element)=(x-1)*idim*idim+(y)*idim+(z)
            element=element+1
          endif
          !  term depending on     (x+1,y,z)
          if (x<idim) then 
            val(element)=-b3(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            icol(element)=(x)*idim*idim+(y-1)*idim+(z)
            element=element+1
          endif
          irow(1:element-1)=glob_row
          ia=glob_row

          t3 = psb_wtime()
          call psb_spins(element-1,irow,icol,val,a,desc_a,info)
          if(info /= 0) exit
          tins = tins + (psb_wtime()-t3)
          call psb_geins(1,(/ia/),zt(1:1),b,desc_a,info)
          if(info /= 0) exit
          zt(1)=0.e0
          call psb_geins(1,(/ia/),zt(1:1),xv,desc_a,info)
          if(info /= 0) exit
        end if
      end do
    end do

    call psb_barrier(ictxt)    
    t2 = psb_wtime()-t1

    if(info /= 0) then
      info=4010
      call psb_errpush(info,name)
      goto 9999
    end if

    deallocate(val,irow,icol)

    t1 = psb_wtime()
    call psb_cdasb(desc_a,info)
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_)
    call psb_barrier(ictxt)
    tasb = psb_wtime()-t1
    if(info /= 0) then
      info=4010
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_amx(ictxt,t2)
    call psb_amx(ictxt,tins)
    call psb_amx(ictxt,tasb)

    if(iam == psb_root_) then
      write(*,'("The matrix has been generated and assembeld in ",a3," format.")')&
           &   a%fida(1:3)
      write(*,'("-pspins time   : ",es10.4)')tins
      write(*,'("-insert time   : ",es10.4)')t2
      write(*,'("-assembly time : ",es10.4)')tasb
    end if

    call psb_geasb(b,desc_a,info)
    call psb_geasb(xv,desc_a,info)
    if(info /= 0) then
      info=4010
      call psb_errpush(info,name)
      goto 9999
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
