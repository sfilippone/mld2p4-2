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
!
! File: ppde3d.f90
!
! Program: ppde3d
! This sample program solves a linear system obtained by discretizing a
! PDE with Dirichlet BCs. 
! 
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
module ppde3d_mod
contains
  !
  ! functions parametrizing the differential equation 
  !  
  function b1(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) :: b1
    real(psb_dpk_), intent(in) :: x,y,z
    b1=0.d0/sqrt(3.d0)
  end function b1
  function b2(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  b2
    real(psb_dpk_), intent(in) :: x,y,z
    b2=0.d0/sqrt(3.d0)
  end function b2
  function b3(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  b3
    real(psb_dpk_), intent(in) :: x,y,z      
    b3=0.d0/sqrt(3.d0)
  end function b3
  function c(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  c
    real(psb_dpk_), intent(in) :: x,y,z      
    c=0.d0
  end function c
  function a1(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  a1   
    real(psb_dpk_), intent(in) :: x,y,z
    a1=1.d0!/80
  end function a1
  function a2(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  a2
    real(psb_dpk_), intent(in) :: x,y,z
    a2=1.d0!/80
  end function a2
  function a3(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  a3
    real(psb_dpk_), intent(in) :: x,y,z
    a3=1.d0!/80
  end function a3
  function g(x,y,z)
    use psb_base_mod, only : psb_dpk_, done, dzero
    real(psb_dpk_) ::  g
    real(psb_dpk_), intent(in) :: x,y,z
    g = dzero
    if (x == done) then
      g = done
    else if (x == dzero) then 
      g = exp(y**2-z**2)
    end if
  end function g
end module ppde3d_mod

program ppde3d
  use psb_base_mod
  use mld_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use data_input
  use ppde3d_mod
  implicit none

  ! input parameters
  character(len=20) :: kmethd, ptype
  character(len=5)  :: afmt
  integer(psb_ipk_) :: idim

  ! miscellaneous 
  real(psb_dpk_), parameter :: one = 1.d0
  real(psb_dpk_) :: t1, t2, tprec 

  ! sparse matrix and preconditioner
  type(psb_dspmat_type) :: a
  type(mld_dprec_type)  :: prec
  ! descriptor
  type(psb_desc_type)   :: desc_a
  ! dense vectors
  type(psb_d_vect_type) :: x,b
  ! parallel environment
  integer(psb_ipk_) :: ictxt, iam, np

  ! solver parameters
  integer(psb_ipk_)        :: iter, itmax,itrace, istopc, irst, nlv
  integer(psb_long_int_k_) :: amatsize, precsize, descsize
  real(psb_dpk_)   :: err, eps

  type precdata
    character(len=20)  :: descr       ! verbose description of the prec
    character(len=10)  :: prec        ! overall prectype
    integer(psb_ipk_)  :: novr        ! number of overlap layers
    integer(psb_ipk_)  :: jsweeps     ! Jacobi/smoother sweeps
    character(len=16)  :: restr       ! restriction  over application of as
    character(len=16)  :: prol        ! prolongation over application of as
    character(len=16)  :: solve       ! Solver  type: ILU, SuperLU, UMFPACK. 
    integer(psb_ipk_)  :: fill1       ! Fill-in for factorization 1
    integer(psb_ipk_)  :: svsweeps    ! Solver sweeps for GS
    real(psb_dpk_)     :: thr1        ! Threshold for fact. 1 ILU(T)
    character(len=16)  :: smther      ! Smoother                            
    integer(psb_ipk_)  :: nlev        ! Number of levels in multilevel prec. 
    character(len=16)  :: aggrkind    ! smoothed/raw aggregatin
    character(len=16)  :: aggr_alg    ! local or global aggregation
    character(len=16)  :: mltype      ! additive or multiplicative 2nd level prec
    character(len=16)  :: smthpos     ! side: pre, post, both smoothing
    integer(psb_ipk_)  :: csize       ! aggregation size at which to stop.
    character(len=16)  :: cmat        ! coarse mat
    character(len=16)  :: csolve      ! Coarse solver: bjac, umf, slu, sludist
    character(len=16)  :: csbsolve    ! Coarse subsolver: ILU, ILU(T), SuperLU, UMFPACK. 
    integer(psb_ipk_)  :: cfill       ! Fill-in for factorization 1
    real(psb_dpk_)     :: cthres      ! Threshold for fact. 1 ILU(T)
    integer(psb_ipk_)  :: cjswp       ! Jacobi sweeps
    real(psb_dpk_)     :: athres      ! smoother aggregation threshold
  end type precdata
  type(precdata)     :: prectype
  type(psb_d_coo_sparse_mat) :: acoo
  ! other variables
  integer(psb_ipk_)  :: info, i
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
  call psb_set_errverbosity(itwo)
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
  call psb_gen_pde3d(ictxt,idim,a,b,x,desc_a,afmt,&
       & a1,a2,a3,b1,b2,b3,c,g,info)  
  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='create_matrix'
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
    call mld_precset(prec,'smoother_type',   prectype%smther,  info)
    call mld_precset(prec,'smoother_sweeps', prectype%jsweeps, info)
    call mld_precset(prec,'sub_ovr',         prectype%novr,    info)
    call mld_precset(prec,'sub_restr',       prectype%restr,   info)
    call mld_precset(prec,'sub_prol',        prectype%prol,    info)
    call mld_precset(prec,'sub_solve',       prectype%solve,   info)
    call mld_precset(prec,'sub_fillin',      prectype%fill1,   info)
    call mld_precset(prec,'solver_sweeps',   prectype%svsweeps,   info)
    call mld_precset(prec,'sub_iluthrs',     prectype%thr1,    info)
    call mld_precset(prec,'aggr_kind',       prectype%aggrkind,info)
    call mld_precset(prec,'aggr_alg',        prectype%aggr_alg,info)
    call mld_precset(prec,'ml_type',         prectype%mltype,  info)
    call mld_precset(prec,'smoother_pos',    prectype%smthpos, info)
    if (prectype%athres >= dzero) &
         & call mld_precset(prec,'aggr_thresh',     prectype%athres,  info)
    call mld_precset(prec,'coarse_solve',    prectype%csolve,  info)
    call mld_precset(prec,'coarse_subsolve', prectype%csbsolve,info)
    call mld_precset(prec,'coarse_mat',      prectype%cmat,    info)
    call mld_precset(prec,'coarse_fillin',   prectype%cfill,   info)
    call mld_precset(prec,'coarse_iluthrs',  prectype%cthres,  info)
    call mld_precset(prec,'coarse_sweeps',   prectype%cjswp,   info)
    call mld_precset(prec,'coarse_aggr_size', prectype%csize,  info)
  else
    nlv = 1
    call mld_precinit(prec,prectype%prec,       info,       nlev=nlv)
    call mld_precset(prec,'smoother_sweeps', prectype%jsweeps,  info)
    call mld_precset(prec,'sub_ovr',         prectype%novr,     info)
    call mld_precset(prec,'sub_restr',       prectype%restr,    info)
    call mld_precset(prec,'sub_prol',        prectype%prol,     info)
    call mld_precset(prec,'sub_solve',       prectype%solve,    info)
    call mld_precset(prec,'sub_fillin',      prectype%fill1,    info)
    call mld_precset(prec,'solver_sweeps',   prectype%svsweeps, info)
    call mld_precset(prec,'sub_iluthrs',     prectype%thr1,     info)
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
  call psb_exit(ictxt)
  stop

9999 continue
  call psb_error(ictxt)

contains
  !
  ! get iteration parameters from standard input
  !
  subroutine  get_parms(ictxt,kmethd,prectype,afmt,idim,istopc,itmax,itrace,irst,eps)
    integer(psb_ipk_) :: ictxt
    type(precdata)    :: prectype
    character(len=*)  :: kmethd, afmt
    integer(psb_ipk_) :: idim, istopc,itmax,itrace,irst
    integer(psb_ipk_) :: np, iam, info
    real(psb_dpk_)    :: eps
    character(len=20) :: buffer

    call psb_info(ictxt, iam, np)

    if (iam == psb_root_) then
      call read_data(kmethd,psb_inp_unit)
      call read_data(afmt,psb_inp_unit)
      call read_data(idim,psb_inp_unit)
      call read_data(istopc,psb_inp_unit)
      call read_data(itmax,psb_inp_unit)
      call read_data(itrace,psb_inp_unit)
      call read_data(irst,psb_inp_unit)
      call read_data(eps,psb_inp_unit)
      call read_data(prectype%descr,psb_inp_unit)       ! verbose description of the prec
      call read_data(prectype%prec,psb_inp_unit)        ! overall prectype
      call read_data(prectype%novr,psb_inp_unit)        ! number of overlap layers
      call read_data(prectype%restr,psb_inp_unit)       ! restriction  over application of as
      call read_data(prectype%prol,psb_inp_unit)        ! prolongation over application of as
      call read_data(prectype%solve,psb_inp_unit)       ! Factorization type: ILU, SuperLU, UMFPACK. 
      call read_data(prectype%svsweeps,psb_inp_unit)    ! Solver sweeps
      call read_data(prectype%fill1,psb_inp_unit)       ! Fill-in for factorization 1
      call read_data(prectype%thr1,psb_inp_unit)        ! Threshold for fact. 1 ILU(T)
      call read_data(prectype%jsweeps,psb_inp_unit)     ! Jacobi sweeps for PJAC
      if (psb_toupper(prectype%prec) == 'ML') then 
        call read_data(prectype%smther,psb_inp_unit)      ! Smoother type.
        call read_data(prectype%nlev,psb_inp_unit)        ! Number of levels in multilevel prec. 
        call read_data(prectype%aggrkind,psb_inp_unit)    ! smoothed/raw aggregatin
        call read_data(prectype%aggr_alg,psb_inp_unit)    ! local or global aggregation
        call read_data(prectype%mltype,psb_inp_unit)      ! additive or multiplicative 2nd level prec
        call read_data(prectype%smthpos,psb_inp_unit)     ! side: pre, post, both smoothing
        call read_data(prectype%cmat,psb_inp_unit)        ! coarse mat
        call read_data(prectype%csolve,psb_inp_unit)      ! Factorization type: ILU, SuperLU, UMFPACK. 
        call read_data(prectype%csbsolve,psb_inp_unit)    ! Factorization type: ILU, SuperLU, UMFPACK. 
        call read_data(prectype%cfill,psb_inp_unit)       ! Fill-in for factorization 1
        call read_data(prectype%cthres,psb_inp_unit)      ! Threshold for fact. 1 ILU(T)
        call read_data(prectype%cjswp,psb_inp_unit)       ! Jacobi sweeps
        call read_data(prectype%athres,psb_inp_unit)      ! smoother aggr thresh
        call read_data(prectype%csize,psb_inp_unit)       ! coarse size
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
    call psb_bcast(ictxt,prectype%svsweeps)    ! Sweeps for inner GS solver
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
      write(psb_out_unit,'("Grid dimensions      : ",i4,"x",i4,"x",i4)')idim,idim,idim
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
    integer(psb_ipk_) :: iout
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

end program ppde3d

