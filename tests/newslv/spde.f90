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
! File: ppde.f90
!
! Program: ppde
! This sample program solves a linear system obtained by discretizing a
! PDE with Dirichlet BCs. 
! 
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
program spde
  use psb_sparse_mod
  use mld_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use data_input
  implicit none

  ! input parameters
  character(len=20) :: kmethd, ptype
  character(len=5)  :: afmt
  integer   :: idim

  ! miscellaneous 
  real(psb_spk_), parameter :: one = 1.0
  real(psb_dpk_) :: t1, t2, tprec 

  ! sparse matrix and preconditioner
  type(psb_sspmat_type) :: a
  type(mld_sprec_type)  :: prec
  ! descriptor
  type(psb_desc_type)   :: desc_a
  ! dense matrices
  real(psb_spk_), allocatable :: b(:), x(:)
  ! blacs parameters
  integer            :: ictxt, iam, np

  ! solver parameters
  integer          :: iter, itmax,itrace, istopc, irst, nlv
  integer(psb_long_int_k_) :: amatsize, precsize, descsize
  real(psb_spk_)   :: err, eps

  type precdata
    character(len=20)  :: descr       ! verbose description of the prec
    character(len=10)  :: prec        ! overall prectype
    integer            :: novr        ! number of overlap layers
    integer            :: jsweeps     ! Jacobi/smoother sweeps
    character(len=16)  :: restr       ! restriction  over application of as
    character(len=16)  :: prol        ! prolongation over application of as
    character(len=16)  :: solve       ! Solver  type: ILU, SuperLU, UMFPACK. 
    integer            :: fill1       ! Fill-in for factorization 1
    real(psb_spk_)     :: thr1        ! Threshold for fact. 1 ILU(T)
    character(len=16)  :: smther      ! Smoother                            
    integer            :: nlev        ! Number of levels in multilevel prec. 
    character(len=16)  :: aggrkind    ! smoothed/raw aggregatin
    character(len=16)  :: aggr_alg    ! local or global aggregation
    character(len=16)  :: mltype      ! additive or multiplicative 2nd level prec
    character(len=16)  :: smthpos     ! side: pre, post, both smoothing
    character(len=16)  :: cmat        ! coarse mat
    character(len=16)  :: csolve      ! Coarse solver: bjac, umf, slu, sludist
    character(len=16)  :: csbsolve    ! Coarse subsolver: ILU, ILU(T), SuperLU, UMFPACK. 
    integer            :: cfill       ! Fill-in for factorization 1
    real(psb_spk_)     :: cthres      ! Threshold for fact. 1 ILU(T)
    integer            :: cjswp       ! Jacobi sweeps
    real(psb_spk_)     :: athres      ! smoother aggregation threshold
  end type precdata
  type(precdata)     :: prectype
  ! other variables
  integer            :: info
  character(len=20)  :: name,ch_err

  info=psb_success_


  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif
  if(psb_get_errstatus() /= 0) goto 9999
  name='pde90'
  call psb_set_errverbosity(2)

  !
  !  get parameters
  !
  call get_parms(ictxt,kmethd,prectype,afmt,idim,istopc,itmax,itrace,irst,eps)

  !
  !  allocate and fill in the coefficient matrix, rhs and initial guess 
  !

  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call create_matrix(idim,a,b,x,desc_a,ictxt,afmt,info)  
  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='create_matrix'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (iam == psb_root_) write(*,'("Overall matrix creation time : ",es12.5)')t2
  if (iam == psb_root_) write(*,'(" ")')
  !
  !  prepare the preconditioner.
  !  

  if (psb_toupper(prectype%prec) == 'ML') then 
    nlv = prectype%nlev
    call mld_precinit(prec,prectype%prec,       info,         nlev=nlv)
    call mld_precset(prec,mld_smoother_type_,   prectype%smther,  info)
    call mld_precset(prec,mld_smoother_sweeps_, prectype%jsweeps, info)
    call mld_precset(prec,mld_sub_ovr_,         prectype%novr,    info)
    call mld_precset(prec,mld_sub_restr_,       prectype%restr,   info)
    call mld_precset(prec,mld_sub_prol_,        prectype%prol,    info)
    call mld_precset(prec,mld_sub_solve_,       prectype%solve,   info)
    call mld_precset(prec,mld_sub_fillin_,      prectype%fill1,   info)
    call mld_precset(prec,mld_sub_iluthrs_,     prectype%thr1,    info)
    call mld_precset(prec,mld_aggr_kind_,       prectype%aggrkind,info)
    call mld_precset(prec,mld_aggr_alg_,        prectype%aggr_alg,info)
    call mld_precset(prec,mld_ml_type_,         prectype%mltype,  info)
    call mld_precset(prec,mld_smoother_pos_,    prectype%smthpos, info)
    call mld_precset(prec,mld_aggr_thresh_,     prectype%athres,  info)
    call mld_precset(prec,mld_coarse_solve_,    prectype%csolve,  info)
    call mld_precset(prec,mld_coarse_subsolve_, prectype%csbsolve,info)
    call mld_precset(prec,mld_coarse_mat_,      prectype%cmat,    info)
    call mld_precset(prec,mld_coarse_fillin_,   prectype%cfill,   info)
    call mld_precset(prec,mld_coarse_iluthrs_,  prectype%cthres,  info)
    call mld_precset(prec,mld_coarse_sweeps_,   prectype%cjswp,   info)
  else
    nlv = 1
    call mld_precinit(prec,prectype%prec,       info,         nlev=nlv)
    call mld_precset(prec,mld_smoother_sweeps_, prectype%jsweeps, info)
    call mld_precset(prec,mld_sub_ovr_,         prectype%novr,    info)
    call mld_precset(prec,mld_sub_restr_,       prectype%restr,   info)
    call mld_precset(prec,mld_sub_prol_,        prectype%prol,    info)
    call mld_precset(prec,mld_sub_solve_,       prectype%solve,   info)
    call mld_precset(prec,mld_sub_fillin_,      prectype%fill1,   info)
    call mld_precset(prec,mld_sub_iluthrs_,     prectype%thr1,    info)
  end if  
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call mld_precbld(a,desc_a,prec,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_precbld'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  tprec = psb_wtime()-t1

  call psb_amx(ictxt,tprec)

  if (iam == psb_root_) write(*,'("Preconditioner time : ",es12.5)')tprec
  if (iam == psb_root_) call mld_precdescr(prec,info)
  if (iam == psb_root_) write(*,'(" ")')

  !
  ! iterative method parameters 
  !
  if(iam == psb_root_) write(*,'("Calling iterative method ",a)')kmethd
  call psb_barrier(ictxt)
  t1 = psb_wtime()  
  call psb_krylov(kmethd,a,prec,b,x,eps,desc_a,info,& 
       & itmax=itmax,iter=iter,err=err,itrace=itrace,istop=istopc,irst=irst)     

  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='solver routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1
  call psb_amx(ictxt,t2)

  amatsize = psb_sizeof(a)
  descsize = psb_sizeof(desc_a)
  precsize = mld_sizeof(prec)
  call psb_sum(ictxt,amatsize)
  call psb_sum(ictxt,descsize)
  call psb_sum(ictxt,precsize)
  if (iam == psb_root_) then
    write(*,'(" ")')
    write(*,'("Time to solve matrix          : ",es12.5)')t2
    write(*,'("Time per iteration            : ",es12.5)')t2/iter
    write(*,'("Number of iterations          : ",i0)')iter
    write(*,'("Convergence indicator on exit : ",es12.5)')err
    write(*,'("Info  on exit                 : ",i0)')info
    write(*,'("Total memory occupation for A:      ",i12)')amatsize
    write(*,'("Total memory occupation for DESC_A: ",i12)')descsize
    write(*,'("Total memory occupation for PREC:   ",i12)')precsize
  end if

  !  
  !  cleanup storage and exit
  !
  call psb_gefree(b,desc_a,info)
  call psb_gefree(x,desc_a,info)
  call psb_spfree(a,desc_a,info)
  call mld_precfree(prec,info)
  call psb_cdfree(desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='free routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

9999 continue
  if(info /= psb_success_) then
    call psb_error(ictxt)
  end if
  call psb_exit(ictxt)
  stop

contains
  !
  ! get iteration parameters from standard input
  !
  subroutine  get_parms(ictxt,kmethd,prectype,afmt,idim,istopc,itmax,itrace,irst,eps)
    integer           :: ictxt
    type(precdata)    :: prectype
    character(len=*)  :: kmethd, afmt
    integer           :: idim, istopc,itmax,itrace,irst
    integer           :: np, iam, info
    real(psb_spk_)    :: eps
    character(len=20) :: buffer

    call psb_info(ictxt, iam, np)

    if (iam == psb_root_) then
      call read_data(kmethd,5)
      call read_data(afmt,5)
      call read_data(idim,5)
      call read_data(istopc,5)
      call read_data(itmax,5)
      call read_data(itrace,5)
      call read_data(irst,5)
      call read_data(eps,5)
      call read_data(prectype%descr,5)       ! verbose description of the prec
      call read_data(prectype%prec,5)        ! overall prectype
      call read_data(prectype%novr,5)        ! number of overlap layers
      call read_data(prectype%restr,5)       ! restriction  over application of as
      call read_data(prectype%prol,5)        ! prolongation over application of as
      call read_data(prectype%solve,5)       ! Factorization type: ILU, SuperLU, UMFPACK. 
      call read_data(prectype%fill1,5)       ! Fill-in for factorization 1
      call read_data(prectype%thr1,5)        ! Threshold for fact. 1 ILU(T)
      call read_data(prectype%jsweeps,5)     ! Jacobi sweeps for PJAC
      if (psb_toupper(prectype%prec) == 'ML') then 
        call read_data(prectype%smther,5)      ! Smoother type.
        call read_data(prectype%nlev,5)        ! Number of levels in multilevel prec. 
        call read_data(prectype%aggrkind,5)    ! smoothed/raw aggregatin
        call read_data(prectype%aggr_alg,5)    ! local or global aggregation
        call read_data(prectype%mltype,5)      ! additive or multiplicative 2nd level prec
        call read_data(prectype%smthpos,5)     ! side: pre, post, both smoothing
        call read_data(prectype%cmat,5)        ! coarse mat
        call read_data(prectype%csolve,5)      ! Factorization type: ILU, SuperLU, UMFPACK. 
        call read_data(prectype%csbsolve,5)    ! Factorization type: ILU, SuperLU, UMFPACK. 
        call read_data(prectype%cfill,5)       ! Fill-in for factorization 1
        call read_data(prectype%cthres,5)      ! Threshold for fact. 1 ILU(T)
        call read_data(prectype%cjswp,5)       ! Jacobi sweeps
        call read_data(prectype%athres,5)      ! smoother aggr thresh
      end if
    end if

    ! broadcast parameters to all processors
    call psb_bcast(ictxt,kmethd)
    call psb_bcast(ictxt,afmt)
    call psb_bcast(ictxt,idim)
    call psb_bcast(ictxt,istopc)
    call psb_bcast(ictxt,itmax)
    call psb_bcast(ictxt,itrace)
    call psb_bcast(ictxt,irst)
    call psb_bcast(ictxt,eps)


    call psb_bcast(ictxt,prectype%descr)       ! verbose description of the prec
    call psb_bcast(ictxt,prectype%prec)        ! overall prectype
    call psb_bcast(ictxt,prectype%novr)        ! number of overlap layers
    call psb_bcast(ictxt,prectype%restr)       ! restriction  over application of as
    call psb_bcast(ictxt,prectype%prol)        ! prolongation over application of as
    call psb_bcast(ictxt,prectype%solve)       ! Factorization type: ILU, SuperLU, UMFPACK. 
    call psb_bcast(ictxt,prectype%fill1)       ! Fill-in for factorization 1
    call psb_bcast(ictxt,prectype%thr1)        ! Threshold for fact. 1 ILU(T)
    call psb_bcast(ictxt,prectype%jsweeps)        ! Jacobi sweeps
    if (psb_toupper(prectype%prec) == 'ML') then 
      call psb_bcast(ictxt,prectype%smther)      ! Smoother type.
      call psb_bcast(ictxt,prectype%nlev)        ! Number of levels in multilevel prec. 
      call psb_bcast(ictxt,prectype%aggrkind)    ! smoothed/raw aggregatin
      call psb_bcast(ictxt,prectype%aggr_alg)    ! local or global aggregation
      call psb_bcast(ictxt,prectype%mltype)      ! additive or multiplicative 2nd level prec
      call psb_bcast(ictxt,prectype%smthpos)     ! side: pre, post, both smoothing
      call psb_bcast(ictxt,prectype%cmat)        ! coarse mat
      call psb_bcast(ictxt,prectype%csolve)      ! Factorization type: ILU, SuperLU, UMFPACK. 
      call psb_bcast(ictxt,prectype%csbsolve)    ! Factorization type: ILU, SuperLU, UMFPACK. 
      call psb_bcast(ictxt,prectype%cfill)       ! Fill-in for factorization 1
      call psb_bcast(ictxt,prectype%cthres)      ! Threshold for fact. 1 ILU(T)
      call psb_bcast(ictxt,prectype%cjswp)       ! Jacobi sweeps
      call psb_bcast(ictxt,prectype%athres)      ! smoother aggr thresh
    end if

    if (iam == psb_root_) then 
      write(*,'("Solving matrix       : ell1")')      
      write(*,'("Grid dimensions      : ",i4,"x",i4,"x",i4)')idim,idim,idim
      write(*,'("Number of processors : ",i0)') np
      write(*,'("Data distribution    : BLOCK")')
      write(*,'("Preconditioner       : ",a)') prectype%descr
      write(*,'("Iterative method     : ",a)') kmethd
      write(*,'(" ")')
    endif

    return

  end subroutine get_parms
  !
  !  print an error message 
  !  
  subroutine pr_usage(iout)
    integer :: iout
    write(iout,*)'incorrect parameter(s) found'
    write(iout,*)' usage:  pde90 methd prec dim &
         &[istop itmax itrace]'  
    write(iout,*)' where:'
    write(iout,*)'     methd:    cgstab cgs rgmres bicgstabl' 
    write(iout,*)'     prec :    bjac diag none'
    write(iout,*)'     dim       number of points along each axis'
    write(iout,*)'               the size of the resulting linear '
    write(iout,*)'               system is dim**3'
    write(iout,*)'     istop     stopping criterion  1, 2  '
    write(iout,*)'     itmax     maximum number of iterations [500] '
    write(iout,*)'     itrace    <=0  (no tracing, default) or '  
    write(iout,*)'               >= 1 do tracing every itrace'
    write(iout,*)'               iterations ' 
  end subroutine pr_usage

  !
  !  subroutine to allocate and fill in the coefficient matrix and
  !  the rhs. 
  !
  subroutine create_matrix(idim,a,b,xv,desc_a,ictxt,afmt,info)
    !
    !   discretize the partial diferential equation
    ! 
    !   b1 dd(u)  b2 dd(u)    b3 dd(u)    a1 d(u)   a2 d(u)  a3 d(u)  
    ! -   ------ -  ------ -  ------ -  -----  -  ------  -  ------ + a4 u 
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
    use psb_sparse_mod
    implicit none
    integer                        :: idim
    integer, parameter             :: nb=20
    real(psb_spk_), allocatable    :: b(:),xv(:)
    type(psb_desc_type)            :: desc_a
    integer                        :: ictxt, info
    character                      :: afmt*5
    type(psb_sspmat_type)    :: a
    real(psb_spk_)           :: zt(nb),x,y,z
    integer                  :: m,n,nnz,glob_row,nlr,i,ii,ib,k
    integer                  :: ix,iy,iz,ia,indx_owner
    integer                  :: np, iam, nr, nt
    integer                  :: element
    integer, allocatable     :: irow(:),icol(:),myidx(:)
    real(psb_spk_), allocatable :: val(:)
    ! deltah dimension of each grid cell
    ! deltat discretization time
    real(psb_spk_)         :: deltah
    real(psb_spk_),parameter   :: rhs=0.0,one=1.0,zero=0.0
    real(psb_dpk_)   :: t0, t1, t2, t3, tasb, talc, ttot, tgen 
    real(psb_spk_)   :: a1, a2, a3, a4, b1, b2, b3 
    external           :: a1, a2, a3, a4, b1, b2, b3
    integer            :: err_act

    character(len=20)  :: name, ch_err

    info = psb_success_
    name = 'create_matrix'
    call psb_erractionsave(err_act)

    call psb_info(ictxt, iam, np)

    deltah = 1.0/(idim-1)

    ! initialize array descriptor and sparse matrix storage. provide an
    ! estimate of the number of non zeroes 

    m   = idim*idim*idim
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
      ch_err='allocation rout.'
      call psb_errpush(info,name,a_err=ch_err)
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
        if (mod(glob_row,(idim*idim)) == 0) then
          ix = glob_row/(idim*idim)
        else
          ix = glob_row/(idim*idim)+1
        endif
        if (mod((glob_row-(ix-1)*idim*idim),idim) == 0) then
          iy = (glob_row-(ix-1)*idim*idim)/idim
        else
          iy = (glob_row-(ix-1)*idim*idim)/idim+1
        endif
        iz = glob_row-(ix-1)*idim*idim-(iy-1)*idim
        ! x, y, x coordinates
        x = ix*deltah
        y = iy*deltah
        z = iz*deltah

        ! check on boundary points 
        zt(k) = 0.d0
        ! internal point: build discretization
        !   
        !  term depending on   (x-1,y,z)
        !
        if (ix == 1) then 
          val(element)=-b1(x,y,z)-a1(x,y,z)
          val(element) = val(element)/(deltah*&
               & deltah)
          zt(k) = exp(-y**2-z**2)*(-val(element))
        else
          val(element)=-b1(x,y,z)-a1(x,y,z)
          val(element) = val(element)/(deltah*&
               & deltah)
          icol(element) = (ix-2)*idim*idim+(iy-1)*idim+(iz)
          irow(element) = glob_row
          element       = element+1
        endif
        !  term depending on     (x,y-1,z)
        if (iy == 1) then 
          val(element)=-b2(x,y,z)-a2(x,y,z)
          val(element) = val(element)/(deltah*&
               & deltah)
          zt(k) = exp(-y**2-z**2)*exp(-x)*(-val(element))  
        else
          val(element)=-b2(x,y,z)-a2(x,y,z)
          val(element) = val(element)/(deltah*deltah)
          icol(element) = (ix-1)*idim*idim+(iy-2)*idim+(iz)
          irow(element) = glob_row
          element       = element+1
        endif
        !  term depending on     (x,y,z-1)
        if (iz == 1) then 
          val(element)=-b3(x,y,z)-a3(x,y,z)
          val(element) = val(element)/(deltah*deltah)
          zt(k) = exp(-y**2-z**2)*exp(-x)*(-val(element))  
        else
          val(element)=-b3(x,y,z)-a3(x,y,z)
          val(element) = val(element)/(deltah*deltah)
          icol(element) = (ix-1)*idim*idim+(iy-1)*idim+(iz-1)
          irow(element) = glob_row
          element       = element+1
        endif
        !  term depending on     (x,y,z)
        val(element)=2*b1(x,y,z) + 2*b2(x,y,z)&
             & + 2*b3(x,y,z) + a1(x,y,z)&
             & + a2(x,y,z) + a3(x,y,z)
        val(element) = val(element)/(deltah*deltah)
        icol(element) = (ix-1)*idim*idim+(iy-1)*idim+(iz)
        irow(element) = glob_row
        element       = element+1                  
        !  term depending on     (x,y,z+1)
        if (iz == idim) then 
          val(element)=-b1(x,y,z)
          val(element) = val(element)/(deltah*deltah)
          zt(k) = exp(-y**2-z**2)*exp(-x)*(-val(element))  
        else
          val(element)=-b1(x,y,z)
          val(element) = val(element)/(deltah*deltah)
          icol(element) = (ix-1)*idim*idim+(iy-1)*idim+(iz+1)
          irow(element) = glob_row
          element       = element+1
        endif
        !  term depending on     (x,y+1,z)
        if (iy == idim) then 
          val(element)=-b2(x,y,z)
          val(element) = val(element)/(deltah*deltah)
          zt(k) = exp(-y**2-z**2)*exp(-x)*(-val(element))  
        else
          val(element)=-b2(x,y,z)
          val(element) = val(element)/(deltah*deltah)
          icol(element) = (ix-1)*idim*idim+(iy)*idim+(iz)
          irow(element) = glob_row
          element       = element+1
        endif
        !  term depending on     (x+1,y,z)
        if (ix<idim) then 
          val(element)=-b3(x,y,z)
          val(element) = val(element)/(deltah*deltah)
          icol(element) = (ix)*idim*idim+(iy-1)*idim+(iz)
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
      ch_err='insert rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    deallocate(val,irow,icol)

    call psb_barrier(ictxt)
    t1 = psb_wtime()
    call psb_cdasb(desc_a,info)
    if (info == psb_success_) &
         & call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    call psb_barrier(ictxt)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='asb rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    call psb_geasb(b,desc_a,info)
    call psb_geasb(xv,desc_a,info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='asb rout.'
      call psb_errpush(info,name,a_err=ch_err)
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
      ch_err = a%get_fmt()
      write(*,'("The matrix has been generated and assembled in ",a3," format.")')&
           &   ch_err(1:3)
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
end program spde
!
! functions parametrizing the differential equation 
!  
function a1(x,y,z)
  use psb_sparse_mod, only : psb_spk_
  real(psb_spk_) :: a1
  real(psb_spk_) :: x,y,z
  a1=1.e0
end function a1
function a2(x,y,z)
  use psb_sparse_mod, only : psb_spk_
  real(psb_spk_) ::  a2
  real(psb_spk_) :: x,y,z
  a2=2.e1*y
end function a2
function a3(x,y,z)
  use psb_sparse_mod, only : psb_spk_
  real(psb_spk_) ::  a3
  real(psb_spk_) :: x,y,z      
  a3=1.e0
end function a3
function a4(x,y,z)
  use psb_sparse_mod, only : psb_spk_
  real(psb_spk_) ::  a4
  real(psb_spk_) :: x,y,z      
  a4=1.e0
end function a4
function b1(x,y,z)
  use psb_sparse_mod, only : psb_spk_
  real(psb_spk_) ::  b1   
  real(psb_spk_) :: x,y,z
  b1=1.e0
end function b1
function b2(x,y,z)
  use psb_sparse_mod, only : psb_spk_
  real(psb_spk_) ::  b2
  real(psb_spk_) :: x,y,z
  b2=1.e0
end function b2
function b3(x,y,z)
  use psb_sparse_mod, only : psb_spk_
  real(psb_spk_) ::  b3
  real(psb_spk_) :: x,y,z
  b3=1.e0
end function b3


