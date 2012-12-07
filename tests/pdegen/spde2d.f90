!!$ 
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010,2012
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
! File: spde2d.f90
!
! Program: spde2d
! This sample program solves a linear system obtained by discretizing a
! PDE with Dirichlet BCs. 
! 
!
! The PDE is a general second order equation in 2d
!
!   a1 dd(u)  a2 dd(u)   b1 d(u)   b2 d(u) 
! -   ------ -  ------   -----  +  ------  + c u = f
!      dxdx     dydy        dx       dy    
!
! with Dirichlet boundary conditions
!   u = g 
!
!  on the unit square  0<=x,y<=1.
!
!
! Note that if b1=b2=c=0., the PDE is the  Laplace equation.
!
! In this sample program the index space of the discretized
! computational domain is first numbered sequentially in a standard way, 
! then the corresponding vector is distributed according to a BLOCK
! data distribution.
!
module spde2d_mod
contains

  !
  ! functions parametrizing the differential equation 
  !  
  function b1(x,y)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) :: b1
    real(psb_spk_), intent(in) :: x,y
    b1=1.e0/sqrt(2.e0)
  end function b1
  function b2(x,y)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  b2
    real(psb_spk_), intent(in) :: x,y
    b2=1.e0/sqrt(2.e0)
  end function b2
  function c(x,y)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  c
    real(psb_spk_), intent(in) :: x,y
    c=0.e0
  end function c
  function a1(x,y)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  a1   
    real(psb_spk_), intent(in) :: x,y
    a1=1.e0/80
  end function a1
  function a2(x,y)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  a2
    real(psb_spk_), intent(in) :: x,y
    a2=1.e0/80
  end function a2
  function g(x,y)
    use psb_base_mod, only : psb_spk_, sone, szero
    real(psb_spk_) ::  g
    real(psb_spk_), intent(in) :: x,y
    g = szero
    if (x == sone) then
      g = sone
    else if (x == szero) then 
      g = exp(-y**2)
    end if
  end function g
end module spde2d_mod

program spde2d
  use psb_base_mod
  use mld_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use data_input
  use spde2d_mod
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
  type(psb_s_vect_type)  :: x,b
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
    integer            :: csize       ! aggregation size at which to stop.
    character(len=16)  :: cmat        ! coarse mat
    character(len=16)  :: csolve      ! Coarse solver: bjac, umf, slu, sludist
    character(len=16)  :: csbsolve    ! Coarse subsolver: ILU, ILU(T), SuperLU, UMFPACK. 
    integer            :: cfill       ! Fill-in for factorization 1
    real(psb_spk_)     :: cthres      ! Threshold for fact. 1 ILU(T)
    integer            :: cjswp       ! Jacobi sweeps
    real(psb_spk_)     :: athres      ! smoother aggregation threshold
  end type precdata
  type(precdata)     :: prectype
  type(psb_s_coo_sparse_mat) :: acoo
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
  name='pde2d90'
  call psb_set_errverbosity(2)
  !
  ! Hello world
  !
  if (iam == psb_root_) then 
    write(*,*) 'Welcome to MLD2P4 version: ',mld_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
  end if
  !
  !  get parameters
  !
  call get_parms(ictxt,kmethd,prectype,afmt,idim,istopc,itmax,itrace,irst,eps)

  !
  !  allocate and fill in the coefficient matrix, rhs and initial guess 
  !
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call psb_gen_pde2d(ictxt,idim,a,b,x,desc_a,afmt,a1,a2,b1,b2,c,g,info)  
  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_gen_pde2d'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (iam == psb_root_) &
       & write(psb_out_unit,'("Overall matrix creation time : ",es12.5)')t2
  if (iam == psb_root_) &
       & write(psb_out_unit,'(" ")')
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
    if (prectype%athres >= szero) &
         & call mld_precset(prec,mld_aggr_thresh_,     prectype%athres,  info)
    call mld_precset(prec,mld_coarse_solve_,    prectype%csolve,  info)
    call mld_precset(prec,mld_coarse_subsolve_, prectype%csbsolve,info)
    call mld_precset(prec,mld_coarse_mat_,      prectype%cmat,    info)
    call mld_precset(prec,mld_coarse_fillin_,   prectype%cfill,   info)
    call mld_precset(prec,mld_coarse_iluthrs_,  prectype%cthres,  info)
    call mld_precset(prec,mld_coarse_sweeps_,   prectype%cjswp,   info)
    call mld_precset(prec,mld_coarse_aggr_size_, prectype%csize,  info)
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
!!$  call prec%dump(info,prefix='test-ml',ac=.true.,solver=.true.,smoother=.true.)

  call psb_amx(ictxt,tprec)

  if (iam == psb_root_) &
       & write(psb_out_unit,'("Preconditioner time : ",es12.5)')tprec
  if (iam == psb_root_) call mld_precdescr(prec,info)
  if (iam == psb_root_) &
       & write(psb_out_unit,'(" ")')

  !
  ! iterative method parameters 
  !
  if(iam == psb_root_) &
       & write(psb_out_unit,'("Calling iterative method ",a)')kmethd
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
  amatsize = a%sizeof()
  descsize = desc_a%sizeof()
  precsize = prec%sizeof()
  call psb_sum(ictxt,amatsize)
  call psb_sum(ictxt,descsize)
  call psb_sum(ictxt,precsize)

  if (iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Time to solve matrix          : ",es12.5)')  t2
    write(psb_out_unit,'("Time per iteration            : ",es12.5)')  t2/iter
    write(psb_out_unit,'("Number of iterations          : ",i0)')      iter
    write(psb_out_unit,'("Convergence indicator on exit : ",es12.5)')  err
    write(psb_out_unit,'("Info  on exit                 : ",i0)')      info
    write(psb_out_unit,'("Total memory occupation for A:      ",i12)') amatsize
    write(psb_out_unit,'("Total memory occupation for DESC_A: ",i12)') descsize
    write(psb_out_unit,'("Total memory occupation for PREC:   ",i12)') precsize
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
        call read_data(prectype%csize,5)       ! coarse size
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
      call psb_bcast(ictxt,prectype%csize)       ! coarse size
    end if

    if (iam == psb_root_) then 
      write(psb_out_unit,'("Solving matrix       : ell1")')      
      write(psb_out_unit,'("Grid dimensions      : ",i4,"x",i4)')idim,idim
      write(psb_out_unit,'("Number of processors : ",i0)') np
      write(psb_out_unit,'("Data distribution    : BLOCK")')
      write(psb_out_unit,'("Preconditioner       : ",a)') prectype%descr
      write(psb_out_unit,'("Iterative method     : ",a)') kmethd
      write(psb_out_unit,'(" ")')
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

end program spde2d
