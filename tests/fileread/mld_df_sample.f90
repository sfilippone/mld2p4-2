!!$ 
!!$ 
!!$                           MLD2P4  version 2.1
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.4)
!!$  
!!$  (C) Copyright 2008, 2010, 2012, 2015, 2017 , 2016
!!$
!!$                      Salvatore Filippone  Cranfield University
!!$		      Ambra Abdullahi Hassan University of Rome Tor Vergata
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
program mld_df_sample
  use psb_base_mod
  use mld_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use data_input
  implicit none


  ! input parameters
  character(len=80) :: kmethd, mtrx_file, rhs_file, sol_file
  character(len=2)  :: filefmt
  type precdata
    character(len=20)  :: descr       ! verbose description of the prec
    character(len=10)  :: prec        ! overall prectype
    integer(psb_ipk_)  :: novr        ! number of overlap layers
    integer(psb_ipk_)  :: jsweeps     ! Jacobi/smoother sweeps
    character(len=16)  :: restr       ! restriction over application of AS
    character(len=16)  :: prol        ! prolongation over application of AS
    character(len=16)  :: solve       ! factorization type: ILU, SuperLU, UMFPACK 
    character(len=16)  :: post_solve  ! Post Solver  type: ILU, SuperLU, UMFPACK. 
    integer(psb_ipk_)  :: fill        ! fillin for factorization 
    integer(psb_ipk_)  :: svsweeps    ! Solver sweeps for GS
    real(psb_dpk_)     :: thr         ! threshold for fact.  ILU(T)
    character(len=16)  :: smther      ! Smoother                            
    integer(psb_ipk_)  :: nlev        ! number of levels in multilevel prec. 
    integer(psb_ipk_)  :: maxlevs     ! Maximum number of levels in multilevel prec. 
    character(len=16)  :: aggrkind    ! smoothed, raw aggregation
    character(len=16)  :: aggr_alg    ! aggregation algorithm (currently only decoupled)
    character(len=16)  :: aggr_ord    ! Ordering for aggregation
    character(len=16)  :: mltype      ! additive or multiplicative multi-level prec
    character(len=16)  :: smthpos     ! side: pre, post, both smoothing
    integer(psb_ipk_)  :: csize       ! aggregation size at which to stop.
    character(len=16)  :: cmat        ! coarse mat: distributed, replicated
    character(len=16)  :: csolve      ! coarse solver: bjac, umf, slu, sludist
    character(len=16)  :: csbsolve    ! coarse subsolver: ILU, ILU(T), SuperLU, UMFPACK 
    integer(psb_ipk_)  :: cfill       ! fillin for coarse factorization 
    real(psb_dpk_)     :: cthres      ! threshold for coarse fact.  ILU(T)
    integer(psb_ipk_)  :: cjswp       ! block-Jacobi sweeps
    real(psb_dpk_)     :: athres      ! smoothed aggregation threshold
    real(psb_dpk_)     :: ascale      ! smoothed aggregation scale factor
    real(psb_dpk_)     :: mnaggratio  ! Minimum aggregation ratio
    integer(psb_ipk_)  :: n_sweeps       
    integer(psb_ipk_)  :: match_algorithm       

  end type precdata
  type(precdata)       :: prec_choice

  ! sparse matrices
  type(psb_dspmat_type) :: a, aux_a

  ! preconditioner data
  type(mld_dprec_type)  :: prec
  ! dense matrices
  real(psb_dpk_), allocatable, target ::  aux_b(:,:), d(:),  aux_x(:,:)
  real(psb_dpk_), allocatable , save  :: x_col_glob(:), r_col_glob(:)
  real(psb_dpk_), pointer  :: b_col_glob(:), ref_col_glob(:)
  type(psb_d_vect_type)    :: b_col, x_col, r_col, ref_col

  ! communications data structure
  type(psb_desc_type):: desc_a

  integer(psb_ipk_)   :: ictxt, iam, np

  ! solver paramters
  integer(psb_ipk_) :: iter, itmax, ierr, itrace, ircode, ipart,&
       & methd, istopc, irst, nlv
  integer(psb_long_int_k_) :: amatsize, precsize, descsize
  real(psb_dpk_)    :: err, eps

  character(len=5)  :: afmt
  character(len=20) :: name, renum
  integer(psb_ipk_), parameter :: iunit=12
  integer(psb_ipk_) :: iparm(20)
  character(len=40) :: fprefix

  ! other variables
  real(psb_dpk_)     :: resmx, resmxp, xdiffn2, xdiffni, xni, xn2
  integer(psb_ipk_)  :: i,info,j,m_problem
  integer(psb_ipk_)  :: internal, m,ii,nnzero
  real(psb_dpk_)     :: t1, t2, tprec, thier, tslv
  real(psb_dpk_)     :: r_amax, b_amax, scale
  integer(psb_ipk_)  :: nrhs, nrow, n_row, dim, nv, ne
  integer(psb_ipk_), allocatable :: ivg(:), ipv(:)
  logical   :: have_guess=.false., have_ref=.false.

  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif


  name='mld_df_sample'
  if(psb_get_errstatus() /= 0) goto 9999
  info=psb_success_
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
  call get_parms(ictxt,mtrx_file,rhs_file,sol_file,filefmt,kmethd,&
       & prec_choice,ipart,afmt,istopc,itmax,itrace,irst,eps)

  call psb_barrier(ictxt)
  t1 = psb_wtime()  
  ! read the input matrix to be processed and (possibly) the rhs 
  nrhs = 1

  if (iam == psb_root_) then
    select case(psb_toupper(filefmt)) 
    case('MM') 
      ! For Matrix Market we have an input file for the matrix
      ! and an (optional) second file for the RHS. 
      call mm_mat_read(aux_a,info,iunit=iunit,filename=mtrx_file)
      if (info == psb_success_) then 
        if (rhs_file /= 'NONE') then
          call mm_array_read(aux_b,info,iunit=iunit,filename=rhs_file)
        end if
      end if
      if ((info == psb_success_).and.(sol_file /= 'NONE')) then
        call mm_array_read(aux_x,info,iunit=iunit,filename=sol_file)
        have_ref = .true.
      end if


    case ('HB')
      ! For Harwell-Boeing we have a single file which may or may not
      ! contain an RHS.
      call hb_read(aux_a,info,iunit=iunit,b=aux_b,filename=mtrx_file)

    case default
      info = -1 
      write(psb_err_unit,*) 'Wrong choice for fileformat ', filefmt
    end select
    if (info /= psb_success_) then
      write(psb_err_unit,*) 'Error while reading input matrix '
      call psb_abort(ictxt)
    end if

    m_problem = aux_a%get_nrows()
    call psb_bcast(ictxt,m_problem)
    
    ! At this point aux_b may still be unallocated
    if (psb_size(aux_b,dim=ione) == m_problem) then
      ! if any rhs were present, broadcast the first one
      write(psb_err_unit,'("Ok, got an rhs ")')
      b_col_glob =>aux_b(:,1)
    else
      write(psb_out_unit,'("Generating an rhs...")')
      write(psb_out_unit,'(" ")')
      call psb_realloc(m_problem,1,aux_b,ircode)
      if (ircode /= 0) then
        call psb_errpush(psb_err_alloc_dealloc_,name)
        goto 9999
      endif

      b_col_glob => aux_b(:,1)
      do i=1, m_problem
        b_col_glob(i) = 1.d0
      enddo
    endif
  else
    call psb_bcast(ictxt,m_problem)
  end if

    if ((have_ref).and.(psb_size(aux_x,dim=ione) == m_problem)) then
      ! if any reference were present, broadcast the first one
      write(psb_err_unit,'("Ok, got a reference solution ")')
      ref_col_glob =>aux_x(:,1)
    else
      write(psb_out_unit,'("No reference solution...")')
      !!! call psb_realloc(m_problem,1,aux_x,ircode)
      !!! if (ircode /= 0) then
      !!!   call psb_errpush(psb_err_alloc_dealloc_,name)
      !!!   goto 9999
      !!! endif
      !!! ref_col_glob => aux_x(:,1)
      !!! do i=1, m_problem
      !!!   ref_col_glob(i) = 0.d0
      !!! enddo
    endif


  ! switch over different partition types
  if (ipart == 0) then 
    call psb_barrier(ictxt)
    if (iam == psb_root_) write(psb_out_unit,'("Partition type: block")')
    allocate(ivg(m_problem),ipv(np))
    do i=1,m_problem
      call part_block(i,m_problem,np,ipv,nv)
      ivg(i) = ipv(1)
    enddo
    call psb_matdist(aux_a, a, ictxt,desc_a,info,fmt=afmt,v=ivg)
  else if (ipart == 2) then 
    if (iam == psb_root_) then 
      write(psb_out_unit,'("Partition type: graph")')
      write(psb_out_unit,'(" ")')
      !      write(psb_err_unit,'("Build type: graph")')
      call build_mtpart(aux_a,np)
    endif
!!$    call psb_barrier(ictxt)
    call distr_mtpart(psb_root_,ictxt)
    call getv_mtpart(ivg)
    call psb_matdist(aux_a, a, ictxt,desc_a,info,fmt=afmt,v=ivg)
  else 
    if (iam == psb_root_) write(psb_out_unit,'("Partition type: block")')
    call psb_matdist(aux_a, a,  ictxt, desc_a,info,fmt=afmt,parts=part_block)
  end if

  call psb_scatter(b_col_glob,b_col,desc_a,info,root=psb_root_)
  call psb_geall(x_col,desc_a,info)
  call x_col%zero()
  call psb_geasb(x_col,desc_a,info)
  call psb_geall(r_col,desc_a,info)
  call r_col%zero()
  call psb_geasb(r_col,desc_a,info)


  !if (have_ref) call psb_geall(ref_col,desc_a,info)
  !if (have_ref) then
  !  if (iam == psb_root_) write(psb_out_unit,'("Scatter reference solution")')
  !  call psb_scatter(ref_col_glob,ref_col,desc_a,info)
  !end if


  t2 = psb_wtime() - t1


  call psb_amx(ictxt, t2)

  if (iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Time to read and partition matrix : ",es12.5)')t2
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,*) 'Preconditioner: ',prec_choice%descr
  end if

  ! 

  if (psb_toupper(prec_choice%prec) == 'ML') then 
    call mld_precinit(prec,prec_choice%prec,       info)
    if (prec_choice%nlev > 0) then
      ! Force number of levels, so disregard the other related arguments.
      call mld_precset(prec,'n_prec_levs', prec_choice%nlev, info)
    else
      if (prec_choice%csize>0)&
           & call mld_precset(prec,'coarse_aggr_size', prec_choice%csize, info)
      if (prec_choice%maxlevs>0)&
           & call mld_precset(prec,'max_prec_levs', prec_choice%maxlevs,  info)
      if (prec_choice%mnaggratio>0)&
           & call mld_precset(prec,'min_aggr_ratio', prec_choice%mnaggratio,  info)
    end if
    if (prec_choice%athres >= dzero) &
         & call mld_precset(prec,'aggr_thresh',     prec_choice%athres,  info)
    call mld_precset(prec,'aggr_kind',       prec_choice%aggrkind,info)
    call mld_precset(prec,'aggr_alg',        prec_choice%aggr_alg,info)
    call mld_precset(prec,'aggr_ord',        prec_choice%aggr_ord,info)
    call mld_precset(prec,'aggr_scale',      prec_choice%ascale,  info)
    call mld_precset(prec,'smoother_type',   prec_choice%smther,  info)
    call mld_precset(prec,'smoother_sweeps', prec_choice%jsweeps, info)
    call mld_precset(prec,'sub_ovr',         prec_choice%novr,    info)
    call mld_precset(prec,'sub_restr',       prec_choice%restr,   info)
    call mld_precset(prec,'sub_prol',        prec_choice%prol,    info)
    call mld_precset(prec,'sub_solve',       prec_choice%solve,   info)
    call mld_precset(prec,'sub_fillin',      prec_choice%fill,    info)
    call mld_precset(prec,'solver_sweeps',   prec_choice%svsweeps,   info)
    call mld_precset(prec,'sub_iluthrs',     prec_choice%thr,     info)
    call mld_precset(prec,'ml_type',         prec_choice%mltype,  info)
    call mld_precset(prec,'smoother_pos',    prec_choice%smthpos, info)
    call mld_precset(prec,'coarse_solve',    prec_choice%csolve,  info)
    call mld_precset(prec,'coarse_subsolve', prec_choice%csbsolve,info)
    call mld_precset(prec,'coarse_mat',      prec_choice%cmat,    info)
    call mld_precset(prec,'coarse_fillin',   prec_choice%cfill,   info)
    call mld_precset(prec,'coarse_iluthrs',  prec_choice%cthres,  info)
    call mld_precset(prec,'coarse_sweeps',   prec_choice%cjswp,   info)

      call prec%set('smoother_type',   prec_choice%smther,   info,pos='post')
      call prec%set('smoother_sweeps', prec_choice%jsweeps,  info,pos='post')
      call prec%set('sub_solve',       prec_choice%post_solve,   info, pos='post')
      call prec%set('sub_ovr',         prec_choice%novr,     info,pos='post')
      call prec%set('sub_restr',       prec_choice%restr,    info,pos='post')
      call prec%set('sub_prol',        prec_choice%prol,     info,pos='post')
      call prec%set('sub_fillin',      prec_choice%fill,     info,pos='post')
      call prec%set('sub_iluthrs',     prec_choice%thr,      info,pos='post')
      call prec%set('solver_sweeps',   prec_choice%svsweeps, info,pos='post')

    call prec%precv(1)%aggr%set('BCM_SWEEPS',prec_choice%n_sweeps, info)
    call prec%precv(1)%aggr%set('BCM_MATCH_ALG',prec_choice%match_algorithm, info)


    ! building the preconditioner
    call psb_barrier(ictxt)
    t1 = psb_wtime()
    call mld_hierarchy_bld(a,desc_a,prec,info)
    thier = psb_wtime()-t1
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_precbld')
      goto 9999
    end if
    call psb_barrier(ictxt)
    t1 = psb_wtime()
    call mld_smoothers_bld(a,desc_a,prec,info)
    tprec = psb_wtime()-t1
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_precbld')
      goto 9999
    end if

  else
    nlv = 1
    call mld_precinit(prec,prec_choice%prec,info)
    if (psb_toupper(prec_choice%prec) /= 'NONE') then 
      call mld_precset(prec,'smoother_sweeps', prec_choice%jsweeps, info)
      call mld_precset(prec,'sub_ovr',         prec_choice%novr,    info)
      call mld_precset(prec,'sub_restr',       prec_choice%restr,   info)
      call mld_precset(prec,'sub_prol',        prec_choice%prol,    info)
      call mld_precset(prec,'sub_solve',       prec_choice%solve,   info)
      call mld_precset(prec,'sub_fillin',      prec_choice%fill,   info)
      call mld_precset(prec,'sub_iluthrs',     prec_choice%thr,    info)
    end if
    ! building the preconditioner
    thier = dzero
    t1 = psb_wtime()
    call mld_precbld(a,desc_a,prec,info)
    tprec = psb_wtime()-t1
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_precbld')
      goto 9999
    end if
  end if

  call psb_amx(ictxt, thier)
  call psb_amx(ictxt, tprec)

  if(iam == psb_root_) then
    write(psb_out_unit,'("Preconditioner time: ",es12.5)')thier+tprec
    write(psb_out_unit,'(" ")')
  end if

  iparm = 0
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call psb_krylov(kmethd,a,prec,b_col,x_col,eps,desc_a,info,& 
       & itmax=itmax,iter=iter,err=err,itrace=itrace,istop=istopc,irst=irst)     
  call psb_barrier(ictxt)
  tslv = psb_wtime() - t1

  call psb_amx(ictxt,tslv)
  call psb_geaxpby(done,b_col,dzero,r_col,desc_a,info)
  call psb_spmm(-done,a,x_col,done,r_col,desc_a,info)
  resmx  = psb_genrm2(r_col,desc_a,info)
  resmxp = psb_geamax(r_col,desc_a,info)
  
  call psb_bcast(ictxt,have_ref)
  if (have_ref) call psb_geall(ref_col,desc_a,info)
  if (have_ref) then
    if (iam == psb_root_) write(psb_out_unit,'("Scatter reference solution")')
    call psb_scatter(ref_col_glob,ref_col,desc_a,info)
  end if

  ! compute error in solution
  if (have_ref) then
    call psb_geaxpby(-done,x_col,done,ref_col,desc_a,info)
    xdiffn2  = psb_genrm2(ref_col,desc_a,info)
    xdiffni  = psb_geamax(ref_col,desc_a,info)
    xn2      = psb_genrm2(ref_col,desc_a,info)
    xni      = psb_geamax(ref_col,desc_a,info)
  end if

  amatsize = a%sizeof()
  descsize = desc_a%sizeof()
  precsize = prec%sizeof()
  call psb_sum(ictxt,amatsize)
  call psb_sum(ictxt,descsize)
  call psb_sum(ictxt,precsize)
  if (iam == psb_root_) then 
    call mld_precdescr(prec,info)
    write(psb_out_unit,'("Matrix: ",a)')mtrx_file
    write(psb_out_unit,'("Computed solution on ",i8," processors")')np
    write(psb_out_unit,'("Iterations to convergence          : ",i6)')iter
    write(psb_out_unit,'("Error estimate on exit             : ",es12.5)')err
    write(psb_out_unit,'("Number of levels in hierarchy      : ",i12)') prec%get_nlevs()
    write(psb_out_unit,'("Time to build hierarchy            : ",es12.5)')thier
    write(psb_out_unit,'("Time to build smoothers            : ",es12.5)')tprec
    write(psb_out_unit,'("Total time for preconditioner      : ",es12.5)')tprec+thier
    write(psb_out_unit,'("Time to solve system               : ",es12.5)')tslv
    write(psb_out_unit,'("Time per iteration                 : ",es12.5)')tslv/(iter)
    write(psb_out_unit,'("Total time                         : ",es12.5)')tslv+tprec
    write(psb_out_unit,'("Residual norm 2                    : ",es12.5)')resmx
    write(psb_out_unit,'("Residual norm inf                  : ",es12.5)')resmxp
    write(psb_out_unit,'("Total memory occupation for A      : ",i12)')amatsize
    write(psb_out_unit,'("Total memory occupation for DESC_A : ",i12)')descsize
    write(psb_out_unit,'("Total memory occupation for PREC   : ",i12)')precsize
    write(psb_out_unit,'("Storage format for A               : ",a  )')a%get_fmt()
    write(psb_out_unit,'("Storage format for DESC_A          : ",a  )')desc_a%get_fmt()
    if (have_ref) then
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'(2x,a10,9x,a8,4x,a20,5x,a8)') &
        & '||X-XREF||','||XREF||','||X-XREF||/||XREF||','(2-norm)'
      write(psb_out_unit,'(1x,3(e12.6,6x))') xdiffn2,xn2,xdiffn2/xn2
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'(2x,a10,9x,a8,4x,a20,4x,a10)') &
        & '||X-XREF||','||XREF||','||X-XREF||/||XREF||','(inf-norm)'
      write(psb_out_unit,'(1x,3(e12.6,6x))') xdiffni,xni,xdiffni/xni
    end if

  end if

  call psb_gather(x_col_glob,x_col,desc_a,info,root=psb_root_)
  if (info == psb_success_) &
       & call psb_gather(r_col_glob,r_col,desc_a,info,root=psb_root_)
  if (info /= psb_success_) goto 9999
  if (iam == psb_root_) then
    write(psb_err_unit,'(" ")')
    write(psb_err_unit,'("Saving x on file")')
    write(20,*) 'matrix: ',mtrx_file
    write(20,*) 'computed solution on ',np,' processors.'
    write(20,*) 'iterations to convergence: ',iter
    write(20,*) 'error estimate (infinity norm) on exit:', &
         & ' ||r||/(||a||||x||+||b||) = ',err
    write(20,'("Residual norm 2          : ",es12.5)')resmx
    write(20,'("Residual norm inf        : ",es12.5)')resmxp
    write(20,'(a8,4(2x,a20))') 'I','X(I)','R(I)','B(I)'
    do i=1,m_problem
      write(20,998) i,x_col_glob(i),r_col_glob(i),b_col_glob(i)
    enddo
  end if
998 format(i8,4(2x,g20.14))
993 format(i6,4(1x,e12.6))


  call psb_gefree(b_col, desc_a,info)
  call psb_gefree(x_col, desc_a,info)
  call psb_spfree(a, desc_a,info)
  call mld_precfree(prec,info)
  call psb_cdfree(desc_a,info)

  call psb_exit(ictxt)
  stop

9999 continue
  call psb_error(ictxt)

contains
  !
  ! get iteration parameters from standard input
  !
  subroutine  get_parms(icontxt,mtrx,rhs,sol,filefmt,kmethd,&
       & prec, ipart,afmt,istopc,itmax,itrace,irst,eps)

    use psb_base_mod
    implicit none

    integer(psb_ipk_)   :: icontxt
    character(len=*)    :: kmethd, mtrx, rhs, afmt,filefmt, sol
    type(precdata)      :: prec
    real(psb_dpk_)      :: eps
    integer(psb_ipk_)   :: iret, istopc,itmax,itrace, ipart, irst
    integer(psb_ipk_)   :: iam, nm, np, i

    call psb_info(icontxt,iam,np)

    if (iam == psb_root_) then
      ! read input parameters
      call read_data(mtrx,psb_inp_unit)
      call read_data(rhs,psb_inp_unit)
      call read_data(sol,psb_inp_unit)             ! solution file (for comparison)
      call read_data(filefmt,psb_inp_unit)
      call read_data(kmethd,psb_inp_unit)
      call read_data(afmt,psb_inp_unit)
      call read_data(ipart,psb_inp_unit)
      call read_data(istopc,psb_inp_unit)
      call read_data(itmax,psb_inp_unit)
      call read_data(itrace,psb_inp_unit)
      call read_data(irst,psb_inp_unit)
      call read_data(eps,psb_inp_unit)
      call read_data(prec%descr,psb_inp_unit)      ! verbose description of the prec
      call read_data(prec%prec,psb_inp_unit)       ! overall prectype
      call read_data(prec%novr,psb_inp_unit)       ! number of overlap layers
      call read_data(prec%restr,psb_inp_unit)      ! restriction  over application of as
      call read_data(prec%prol,psb_inp_unit)       ! prolongation over application of as
      call read_data(prec%solve,psb_inp_unit)      ! Factorization type: ILU, SuperLU, UMFPACK.
      call read_data(prec%post_solve,psb_inp_unit)       ! Subdomain solver:  DSCALE ILU MILU ILUT FWGS BWGS MUMPS UMF SLU
      call read_data(prec%svsweeps,psb_inp_unit)    ! Solver sweeps (GS) 
      call read_data(prec%fill,psb_inp_unit)       ! Fill-in for factorization 
      call read_data(prec%thr,psb_inp_unit)        ! Threshold for fact.  ILU(T)
      call read_data(prec%jsweeps,psb_inp_unit)    ! Jacobi sweeps for PJAC
      if (psb_toupper(prec%prec) == 'ML') then 
        call read_data(prec%nlev,psb_inp_unit)     ! Number of levels in multilevel prec
        call read_data(prec%csize,psb_inp_unit)       ! coarse size
        call read_data(prec%mnaggratio,psb_inp_unit)  ! Minimum aggregation ratio
        call read_data(prec%athres,psb_inp_unit)      ! smoother aggr thresh
        call read_data(prec%maxlevs,psb_inp_unit)     ! Maximum number of levels
        call read_data(prec%smther,psb_inp_unit)   ! Smoother type.
        call read_data(prec%aggrkind,psb_inp_unit) ! smoothed/raw aggregatin
        call read_data(prec%aggr_alg,psb_inp_unit) ! local or global aggregation
        call read_data(prec%aggr_ord,psb_inp_unit) ! Ordering for aggregation
        call read_data(prec%mltype,psb_inp_unit)   ! additive or multiplicative 2nd level prec
        call read_data(prec%smthpos,psb_inp_unit)  ! side: pre, post, both smoothing
        call read_data(prec%cmat,psb_inp_unit)     ! coarse mat
        call read_data(prec%csolve,psb_inp_unit)   ! Factorization type: BJAC, SuperLU, UMFPACK. 
        call read_data(prec%csbsolve,psb_inp_unit) ! Factorization type: ILU, SuperLU, UMFPACK. 
        call read_data(prec%cfill,psb_inp_unit)    ! Fill-in for factorization 
        call read_data(prec%cthres,psb_inp_unit)   ! Threshold for fact.  ILU(T)
        call read_data(prec%cjswp,psb_inp_unit)    ! Jacobi sweeps
        call read_data(prec%ascale,psb_inp_unit)   ! smoother aggr thresh
!BCMATCH parameters
        call read_data(prec%n_sweeps,psb_inp_unit)       
        call read_data(prec%match_algorithm,psb_inp_unit)       

      end if
    end if

    call psb_bcast(icontxt,mtrx)
    call psb_bcast(icontxt,rhs)
    call psb_bcast(icontxt,filefmt)
    call psb_bcast(icontxt,kmethd)
    call psb_bcast(icontxt,afmt)
    call psb_bcast(icontxt,ipart)
    call psb_bcast(icontxt,istopc)
    call psb_bcast(icontxt,itmax)
    call psb_bcast(icontxt,itrace)
    call psb_bcast(icontxt,irst)
    call psb_bcast(icontxt,eps)
    call psb_bcast(icontxt,prec%descr)       ! verbose description of the prec
    call psb_bcast(icontxt,prec%prec)        ! overall prectype
    call psb_bcast(icontxt,prec%novr)        ! number of overlap layers
    call psb_bcast(icontxt,prec%restr)       ! restriction  over application of as
    call psb_bcast(icontxt,prec%prol)        ! prolongation over application of as
    call psb_bcast(icontxt,prec%solve)       ! Factorization type: ILU, SuperLU, UMFPACK. 
    call psb_bcast(icontxt,prec%fill)        ! Fill-in for factorization 
    call psb_bcast(icontxt,prec%thr)         ! Threshold for fact.  ILU(T)
    call psb_bcast(icontxt,prec%jsweeps)       ! Jacobi sweeps
    if (psb_toupper(prec%prec) == 'ML') then 
      call psb_bcast(icontxt,prec%smther)      ! Smoother type.
      call psb_bcast(icontxt,prec%nlev)        ! Number of levels in multilevel prec. 
      call psb_bcast(ictxt,prec%csize)       ! coarse size
      call psb_bcast(ictxt,prec%mnaggratio)  ! Minimum aggregation ratio
      call psb_bcast(ictxt,prec%athres)      ! smoother aggr thresh
      call psb_bcast(ictxt,prec%maxlevs)     ! Maximum number of levels
      call psb_bcast(icontxt,prec%aggrkind)    ! smoothed/raw aggregatin
      call psb_bcast(icontxt,prec%aggr_alg)    ! local or global aggregation
      call psb_bcast(icontxt,prec%aggr_ord)    ! Ordering for aggregation
      call psb_bcast(icontxt,prec%mltype)      ! additive or multiplicative 2nd level prec
      call psb_bcast(icontxt,prec%smthpos)     ! side: pre, post, both smoothing
      call psb_bcast(icontxt,prec%cmat)        ! coarse mat
      call psb_bcast(icontxt,prec%csolve)      ! Factorization type: ILU, SuperLU, UMFPACK. 
      call psb_bcast(icontxt,prec%csbsolve)    ! Factorization type: ILU, SuperLU, UMFPACK. 
      call psb_bcast(icontxt,prec%cfill)       ! Fill-in for factorization 
      call psb_bcast(icontxt,prec%cthres)      ! Threshold for fact.  ILU(T)
      call psb_bcast(icontxt,prec%cjswp)       ! Jacobi sweeps
      call psb_bcast(icontxt,prec%athres)      ! smoother aggr thresh
      call psb_bcast(icontxt,prec%ascale)      ! smoother aggr scale factor
      call psb_bcast(icontxt,prec%n_sweeps)       
      call psb_bcast(icontxt,prec%match_algorithm)       
    end if

  end subroutine get_parms
  subroutine pr_usage(iout)
    integer(psb_ipk_) iout
    write(iout, *) ' number of parameters is incorrect!'
    write(iout, *) ' use: hb_sample mtrx_file methd prec [ptype &
         &itmax istopc itrace]' 
    write(iout, *) ' where:'
    write(iout, *) '     mtrx_file      is stored in hb format'
    write(iout, *) '     methd          may be: cgstab '
    write(iout, *) '     itmax          max iterations [500]        '
    write(iout, *) '     istopc         stopping criterion [1]      '
    write(iout, *) '     itrace         0  (no tracing, default) or '
    write(iout, *) '                    >= 0 do tracing every itrace'
    write(iout, *) '                    iterations ' 
    write(iout, *) '     prec           may be: ilu diagsc none'
    write(iout, *) '     ptype          partition strategy default 0'
    write(iout, *) '                    0: block partition '
  end subroutine pr_usage
end program mld_df_sample
