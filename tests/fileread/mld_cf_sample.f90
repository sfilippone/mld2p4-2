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
program mld_cf_sample
  use psb_base_mod
  use mld_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use data_input
  implicit none


  ! input parameters

  character(len=40) :: kmethd, mtrx_file, rhs_file, guess_file, sol_file, part
  character(len=2)  :: filefmt

  ! Krylov solver data
  type solverdata
    character(len=40)  :: kmethd      ! Krylov solver
    integer(psb_ipk_)  :: istopc      ! stopping criterion
    integer(psb_ipk_)  :: itmax       ! maximum number of iterations
    integer(psb_ipk_)  :: itrace      ! tracing
    integer(psb_ipk_)  :: irst        ! restart
    real(psb_spk_)     :: eps         ! stopping tolerance
  end type solverdata
  type(solverdata)       :: s_choice

  ! preconditioner data
  type precdata

    ! preconditioner type
    character(len=40)  :: descr       ! verbose description of the prec
    character(len=10)  :: ptype       ! preconditioner type

    integer(psb_ipk_)  :: outer_sweeps ! number of outer sweeps: sweeps for 1-level,
                                       ! AMG cycles for ML
    ! general AMG data
    character(len=16)  :: mlcycle      ! AMG cycle type
    integer(psb_ipk_)  :: maxlevs     ! maximum number of levels in AMG preconditioner

    ! AMG aggregation
    character(len=16)  :: aggr_prol    ! aggregation type: SMOOTHED, NONSMOOTHED
    character(len=16)  :: par_aggr_alg    ! parallel aggregation algorithm: DEC, SYMDEC
    character(len=16)  :: aggr_ord    ! ordering for aggregation: NATURAL, DEGREE
    character(len=16)  :: aggr_filter ! filtering: FILTER, NO_FILTER
    real(psb_spk_)     :: mncrratio  ! minimum aggregation ratio
    real(psb_spk_), allocatable :: athresv(:) ! smoothed aggregation threshold vector
    integer(psb_ipk_)  :: thrvsz      ! size of threshold vector
    real(psb_spk_)     :: athres      ! smoothed aggregation threshold
    integer(psb_ipk_)  :: csize       ! minimum size of coarsest matrix

    ! AMG smoother or pre-smoother; also 1-lev preconditioner
    character(len=16)  :: smther      ! (pre-)smoother type: BJAC, AS
    integer(psb_ipk_)  :: jsweeps     ! (pre-)smoother / 1-lev prec. sweeps
    integer(psb_ipk_)  :: novr        ! number of overlap layers
    character(len=16)  :: restr       ! restriction over application of AS
    character(len=16)  :: prol        ! prolongation over application of AS
    character(len=16)  :: solve       ! local subsolver type: ILU, MILU, ILUT,
                                      ! UMF, MUMPS, SLU, FWGS, BWGS, JAC
    integer(psb_ipk_)  :: fill        ! fill-in for incomplete LU factorization
    real(psb_spk_)     :: thr         ! threshold for ILUT factorization

    ! AMG post-smoother; ignored by 1-lev preconditioner
    character(len=16)  :: smther2     ! post-smoother type: BJAC, AS
    integer(psb_ipk_)  :: jsweeps2    ! post-smoother sweeps
    integer(psb_ipk_)  :: novr2       ! number of overlap layers
    character(len=16)  :: restr2      ! restriction  over application of AS
    character(len=16)  :: prol2       ! prolongation over application of AS
    character(len=16)  :: solve2      ! local subsolver type: ILU, MILU, ILUT,
                                      ! UMF, MUMPS, SLU, FWGS, BWGS, JAC
    integer(psb_ipk_)  :: fill2       ! fill-in for incomplete LU factorization
    real(psb_spk_)     :: thr2        ! threshold for ILUT factorization

    ! coarsest-level solver
    character(len=16)  :: cmat        ! coarsest matrix layout: REPL, DIST
    character(len=16)  :: csolve      ! coarsest-lev solver: BJAC, SLUDIST (distr.
                                      ! mat.); UMF, MUMPS, SLU, ILU, ILUT, MILU
                                      ! (repl. mat.)
    character(len=16)  :: csbsolve    ! coarsest-lev local subsolver: ILU, ILUT,
                                      ! MILU, UMF, MUMPS, SLU
    integer(psb_ipk_)  :: cfill       ! fill-in for incomplete LU factorization
    real(psb_spk_)     :: cthres      ! threshold for ILUT factorization
    integer(psb_ipk_)  :: cjswp       ! sweeps for GS or JAC coarsest-lev subsolver

  end type precdata
  type(precdata)       :: p_choice

  ! sparse matrices
  type(psb_cspmat_type) :: a, aux_a

  ! preconditioner data
  type(mld_cprec_type)  :: prec
  ! dense matrices
  complex(psb_spk_), allocatable, target ::  aux_b(:,:), d(:), aux_g(:,:), aux_x(:,:)
  complex(psb_spk_), allocatable , save  :: x_col_glob(:), r_col_glob(:)
  complex(psb_spk_), pointer  :: b_col_glob(:), ref_col_glob(:), guess_col_glob(:)
  type(psb_c_vect_type)    :: b_col, x_col, r_col, ref_col

  ! communications data structure
  type(psb_desc_type):: desc_a

  integer(psb_ipk_)   :: ictxt, iam, np

  ! solver paramters
  integer(psb_ipk_) :: iter, ircode, nlv
  integer(psb_long_int_k_) :: amatsize, precsize, descsize
  real(psb_spk_)    :: err

  character(len=5)  :: afmt
  character(len=20) :: name, renum
  integer(psb_ipk_), parameter :: iunit=12
  integer(psb_ipk_) :: iparm(20)

  ! other variables
  integer(psb_ipk_)  :: i, info, j, k, m_problem
  integer(psb_ipk_)  :: lbw, ubw, prf
  real(psb_dpk_)     :: t1, t2, tprec, thier, tslv
  real(psb_spk_)     :: resmx, resmxp, xdiffn2, xdiffni, xni, xn2
  integer(psb_ipk_)  :: nrhs, nv
  integer(psb_ipk_), allocatable :: ivg(:), ipv(:), perm(:)
  logical   :: have_guess=.false., have_ref=.false.

  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif


  name='mld_cf_sample'
  if(psb_get_errstatus() /= 0) goto 9999
  info=psb_success_
  call psb_set_errverbosity(itwo)
  !
  ! Hello world
  !
  if (iam == psb_root_) then 
    write(psb_out_unit,*) ' '
    write(psb_out_unit,*) 'Welcome to MLD2P4 version: ',mld_version_string_
    write(psb_out_unit,*) 'This is the ',trim(name),' sample test program'
    write(psb_out_unit,*) ' '
  end if
  !
  ! get parameters
  !
  call get_parms(ictxt,mtrx_file,rhs_file,guess_file,sol_file,filefmt, &
       & part,afmt,s_choice,p_choice)

  call psb_barrier(ictxt)
  t1 = psb_wtime()  
  ! read the input matrix to be processed and (possibly) the rhs,
  ! the initial guess and the reference solution
  nrhs = 1

  if (iam == psb_root_) then
    select case(psb_toupper(filefmt))

    case('MM') 
      ! For Matrix Market we have an input file for the matrix
      ! and (optional) separate files for the rhs, the initial guess
      ! and the reference solution
      call mm_mat_read(aux_a,info,iunit=iunit,filename=mtrx_file)
      if ((info == psb_success_).and.(rhs_file /= 'NONE')) &
           & call mm_array_read(aux_b,info,iunit=iunit,filename=rhs_file)
      if ((info == psb_success_).and.(guess_file /= 'NONE')) then
        call mm_array_read(aux_g,info,iunit=iunit,filename=guess_file)
        have_guess = .true.
      end if
      if ((info == psb_success_).and.(sol_file /= 'NONE')) then
        call mm_array_read(aux_x,info,iunit=iunit,filename=sol_file)
        have_ref = .true.
      end if
      if (info == psb_success_) call aux_a%cscnv(info,type='csr')
    case ('HB')
      ! For Harwell-Boeing we have a single file which may or may not
      ! contain rhs, initial guess and reference solution.
      call hb_read(aux_a,info,iunit=iunit,b=aux_b,&
           & g=aux_g,x=aux_x,filename=mtrx_file)
      have_guess = allocated(aux_g)
      have_ref   = allocated(aux_x)

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
    call psb_bcast(ictxt,have_guess)
    call psb_bcast(ictxt,have_ref)
    
    ! At this point aux_b may still be unallocated
    if (psb_size(aux_b,dim=ione) == m_problem) then
      ! if any rhs were present, broadcast the first one
      write(psb_err_unit,'("Ok, got an rhs ")')
      b_col_glob =>aux_b(:,1)
    else
      write(psb_out_unit,'("Generating an rhs...")')
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

    if ((have_guess).and.(psb_size(aux_g,dim=ione) == m_problem)) then
      ! if any initial guess were present, broadcast the first one
      write(psb_err_unit,'("Ok, got an initial guess ")')
      guess_col_glob =>aux_g(:,1)
    else
      write(psb_out_unit,'("Generating an initial guess...")')
      call psb_realloc(m_problem,1,aux_g,ircode)
      if (ircode /= 0) then
        call psb_errpush(psb_err_alloc_dealloc_,name)
        goto 9999
      endif
      guess_col_glob => aux_g(:,1)
      do i=1, m_problem
        guess_col_glob(i) = 0.d0
      enddo
    endif

    if ((have_ref).and.(psb_size(aux_x,dim=ione) == m_problem)) then
      ! if any reference were present, broadcast the first one
      write(psb_err_unit,'("Ok, got a reference solution ")')
      ref_col_glob =>aux_x(:,1)
    else
      write(psb_out_unit,'("No reference solution...")')
    endif

    ! clean zeros in the input matrix
    call aux_a%clean_zeros(info)

  else
    call psb_bcast(ictxt,m_problem)
    call psb_bcast(ictxt,have_guess)
    call psb_bcast(ictxt,have_ref)
  end if
    

  !
  ! Renumbering (NONE for the moment)
  !
  if (iam==psb_root_) then
    renum='NONE'

    call psb_cmp_bwpf(aux_a,lbw,ubw,prf,info)
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,*) 'Bandwidth and profile: ',lbw,ubw,prf
    write(psb_out_unit,*) 'Renumbering algorithm: ',psb_toupper(renum)

    if (trim(psb_toupper(renum))/='NONE') then
      call psb_mat_renum(renum,aux_a,info,perm=perm)
      if (info /= 0) then
        write(psb_err_unit,*) 'Error from RENUM',info
        goto 9999
      end if
      call psb_gelp('N',perm(1:m_problem),b_col_glob(1:m_problem),info)
      call psb_cmp_bwpf(aux_a,lbw,ubw,prf,info)
      write(psb_out_unit,*) 'Bandwidth and profile (renumberd):',lbw,ubw,prf
    end if

    write(psb_out_unit,'(" ")')

  end if

  !
  ! switch over different partition types
  !
  select case (psb_toupper(part))
  case('BLOCK')
    call psb_barrier(ictxt)
    if (iam == psb_root_) write(psb_out_unit,'("Partition type: block")')
    call psb_matdist(aux_a, a,  ictxt, desc_a,info,fmt=afmt,parts=part_block)
  case('GRAPH')
    if (iam == psb_root_) then 
      write(psb_out_unit,'("Partition type: graph")')
      write(psb_out_unit,'(" ")')
      call build_mtpart(aux_a,np)
    endif
    call distr_mtpart(psb_root_,ictxt)
    call getv_mtpart(ivg)
    call psb_matdist(aux_a, a, ictxt,desc_a,info,fmt=afmt,v=ivg)
  case default
    if (iam == psb_root_) write(psb_out_unit,'("Partition type: block")')
    call psb_matdist(aux_a, a,  ictxt, desc_a,info,fmt=afmt,parts=part_block)
  end select

  !
  ! Scatter rhs, initial guess and reference solution
  !
  call psb_geall(b_col,desc_a,info)
  call psb_geall(x_col,desc_a,info)
  if (have_ref) call psb_geall(ref_col,desc_a,info)

  if (iam == psb_root_) write(psb_out_unit,'("Scatter rhs")')
  call psb_scatter(b_col_glob,b_col,desc_a,info)
  if (iam == psb_root_) write(psb_out_unit,'("Scatter initial guess")')
  call psb_scatter(guess_col_glob,x_col,desc_a,info)
  if (have_ref) then
    if (iam == psb_root_) write(psb_out_unit,'("Scatter reference solution")')
    call psb_scatter(ref_col_glob,ref_col,desc_a,info)
  end if

  t2 = psb_wtime() - t1
  call psb_amx(ictxt, t2)

  if (iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Time to read and partition matrix, rhs(, guess, ref sol) : ",es12.5)')t2
    write(psb_out_unit,'(" ")')
  end if

  !
  ! initialize the preconditioner
  !
  call prec%init(ictxt,p_choice%ptype,info)
  select case(trim(psb_toupper(p_choice%ptype)))
  case ('NONE','NOPREC')
    ! Do nothing, keep defaults

  case ('JACOBI','GS','FWGS','FBGS')
    ! 1-level sweeps from "outer_sweeps"
    call prec%set('smoother_sweeps', p_choice%jsweeps, info)
    
  case ('BJAC')
    call prec%set('smoother_sweeps', p_choice%jsweeps, info)
    call prec%set('sub_solve',       p_choice%solve,   info)
    call prec%set('sub_fillin',      p_choice%fill,    info)
    call prec%set('sub_iluthrs',     p_choice%thr,     info)

  case('AS')
    call prec%set('smoother_sweeps', p_choice%jsweeps, info)
    call prec%set('sub_ovr',         p_choice%novr,    info)
    call prec%set('sub_restr',       p_choice%restr,   info)
    call prec%set('sub_prol',        p_choice%prol,    info)
    call prec%set('sub_solve',       p_choice%solve,   info)
    call prec%set('sub_fillin',      p_choice%fill,    info)
    call prec%set('sub_iluthrs',     p_choice%thr,     info)
    
  case ('ML') 
    ! multilevel preconditioner

    call prec%set('ml_cycle',        p_choice%mlcycle,    info)
    call prec%set('outer_sweeps',    p_choice%outer_sweeps,info)
    if (p_choice%csize>0)&
         & call prec%set('min_coarse_size', p_choice%csize,      info)
    if (p_choice%mncrratio>1)&
         & call prec%set('min_cr_ratio',   p_choice%mncrratio, info)
    if (p_choice%maxlevs>0)&
         & call prec%set('max_levs',    p_choice%maxlevs,    info)
    if (p_choice%athres >= dzero) &
         & call prec%set('aggr_thresh',     p_choice%athres,  info)
    if (p_choice%thrvsz>0) then
      do k=1,min(p_choice%thrvsz,size(prec%precv)-1)
        call prec%set('aggr_thresh',     p_choice%athresv(k),  info,ilev=(k+1))
      end do
    end if

    call prec%set('aggr_prol',       p_choice%aggr_prol,   info)
    call prec%set('par_aggr_alg',    p_choice%par_aggr_alg,   info)
    call prec%set('aggr_ord',        p_choice%aggr_ord,   info)
    call prec%set('aggr_filter',     p_choice%aggr_filter,info)


    call prec%set('smoother_type',   p_choice%smther,     info)
    call prec%set('smoother_sweeps', p_choice%jsweeps,    info)

    select case (psb_toupper(p_choice%smther))
    case ('GS','BWGS','FBGS','JACOBI')
      ! do nothing
    case default
      call prec%set('sub_ovr',         p_choice%novr,       info)
      call prec%set('sub_restr',       p_choice%restr,      info)
      call prec%set('sub_prol',        p_choice%prol,       info)
      call prec%set('sub_solve',       p_choice%solve,      info)
      call prec%set('sub_fillin',      p_choice%fill,       info)
      call prec%set('sub_iluthrs',     p_choice%thr,        info)
    end select

    if (psb_toupper(p_choice%smther2) /= 'NONE') then
      call prec%set('smoother_type',   p_choice%smther2,   info,pos='post')
      call prec%set('smoother_sweeps', p_choice%jsweeps2,  info,pos='post')
      select case (psb_toupper(p_choice%smther2))
      case ('GS','BWGS','FBGS','JACOBI')
        ! do nothing
      case default
        call prec%set('sub_ovr',         p_choice%novr2,     info,pos='post')
        call prec%set('sub_restr',       p_choice%restr2,    info,pos='post')
        call prec%set('sub_prol',        p_choice%prol2,     info,pos='post')
        call prec%set('sub_solve',       p_choice%solve2,    info,pos='post')
        call prec%set('sub_fillin',      p_choice%fill2,     info,pos='post')
        call prec%set('sub_iluthrs',     p_choice%thr2,      info,pos='post')
      end select
    end if

    call prec%set('coarse_solve',    p_choice%csolve,    info)
    if (psb_toupper(p_choice%csolve) == 'BJAC') &
         &  call prec%set('coarse_subsolve', p_choice%csbsolve,  info)
    call prec%set('coarse_mat',      p_choice%cmat,      info)
    call prec%set('coarse_fillin',   p_choice%cfill,     info)
    call prec%set('coarse_iluthrs',  p_choice%cthres,    info)
    call prec%set('coarse_sweeps',   p_choice%cjswp,     info)

  end select
  
  ! build the preconditioner
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call prec%hierarchy_build(a,desc_a,info)
  thier = psb_wtime()-t1
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_hierarchy_bld')
    goto 9999
  end if
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call prec%smoothers_build(a,desc_a,info)
  tprec = psb_wtime()-t1
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_smoothers_bld')
    goto 9999
  end if

  call psb_amx(ictxt, thier)
  call psb_amx(ictxt, tprec)

  if(iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Preconditioner: ",a)') trim(p_choice%descr)
    write(psb_out_unit,'("Preconditioner time: ",es12.5)')thier+tprec
    write(psb_out_unit,'(" ")')
  end if

  !
  ! iterative method parameters 
  !
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call psb_krylov(s_choice%kmethd,a,prec,b_col,x_col,s_choice%eps,&
       & desc_a,info,itmax=s_choice%itmax,iter=iter,err=err,itrace=s_choice%itrace,&
       & istop=s_choice%istopc,irst=s_choice%irst)
  call psb_barrier(ictxt)
  tslv = psb_wtime() - t1

  call psb_amx(ictxt,tslv)

  ! compute residual norms
  call psb_geall(r_col,desc_a,info)
  call r_col%zero()
  call psb_geasb(r_col,desc_a,info)
  call psb_geaxpby(cone,b_col,czero,r_col,desc_a,info)
  call psb_spmm(-cone,a,x_col,cone,r_col,desc_a,info)
  resmx  = psb_genrm2(r_col,desc_a,info)
  resmxp = psb_geamax(r_col,desc_a,info)

  ! compute error in solution
  if (have_ref) then
    call psb_geaxpby(-cone,x_col,cone,ref_col,desc_a,info)
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
  call prec%descr(iout=psb_out_unit)
  if (iam == psb_root_) then 
    write(psb_out_unit,'("Matrix: ",a)')mtrx_file
    write(psb_out_unit,'("Computed solution on ",i8," processors")')np
    write(psb_out_unit,'("Krylov method                      : ",a)') trim(s_choice%kmethd)
    write(psb_out_unit,'("Preconditioner                     : ",a)') trim(p_choice%descr)
    write(psb_out_unit,'("Iterations to convergence          : ",i12)')iter
    write(psb_out_unit,'("Relative error estimate on exit    : ",es12.5)') err
    write(psb_out_unit,'("Number of levels in hierarchy      : ",i12)') prec%get_nlevs()
    write(psb_out_unit,'("Time to build hierarchy            : ",es12.5)')thier
    write(psb_out_unit,'("Time to build smoothers            : ",es12.5)')tprec
    write(psb_out_unit,'("Total time for preconditioner      : ",es12.5)')tprec+thier
    write(psb_out_unit,'("Time to solve system               : ",es12.5)')tslv
    write(psb_out_unit,'("Time per iteration                 : ",es12.5)')tslv/iter
    write(psb_out_unit,'("Total time                         : ",es12.5)')tslv+tprec+thier
    write(psb_out_unit,'("Residual 2-norm                    : ",es12.5)')resmx
    write(psb_out_unit,'("Residual inf-norm                  : ",es12.5)')resmxp
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
    write(20,*) 'Matrix: ',mtrx_file
    write(20,*) 'Krylov method:',trim(s_choice%kmethd)
    write(20,*) 'Preconditioner:',trim(p_choice%descr)
    write(20,*) 'Computed solution on ',np,' processors.'
    write(20,*) 'Iterations to convergence: ',iter
    write(20,*) 'Error estimate (infinity norm) on exit:', &
         & ' ||r||/||b|| (inf-norm) = ',err
    write(20,'(" Residual 2-norm 2         : ",es12.5)')resmx
    write(20,'(" Residual inf-norm         : ",es12.5)')resmxp
    write(20,'(a8,4(2x,a20))') 'I','X(I)','R(I)','B(I)'
    do i=1,m_problem
      write(20,998) i,x_col_glob(i),r_col_glob(i),b_col_glob(i)
    enddo
  end if
998 format(i8,4(2x,g20.14))
993 format(i6,4(1x,e12.6))


  call psb_gefree(b_col,desc_a,info)
  call psb_gefree(x_col,desc_a,info)
  call psb_gefree(r_col,desc_a,info)
  call psb_gefree(ref_col,desc_a,info)
  call psb_spfree(a, desc_a,info)
  call prec%free(info)
  call psb_cdfree(desc_a,info)

  call psb_exit(ictxt)
  stop

9999 continue
  call psb_error(ictxt)

contains
  !
  ! get iteration parameters from standard input
  !
  subroutine get_parms(icontxt,mtrx,rhs,guess,sol,filefmt,part,afmt,solve,prec)

    implicit none

    integer(psb_ipk_)   :: icontxt
    character(len=*)    :: mtrx, rhs, guess, sol, filefmt, afmt, part
    type(solverdata)    :: solve
    type(precdata)      :: prec
    integer(psb_ipk_)   :: iam, nm, np, inp_unit
    character(len=1024)   :: filename

    call psb_info(icontxt,iam,np)

    if (iam == psb_root_) then
      ! read input data
      if (command_argument_count()>0) then
        call get_command_argument(1,filename)
        inp_unit = 30
        open(inp_unit,file=filename,action='read',iostat=info)
        if (info /= 0) then
          write(psb_err_unit,*) 'Could not open file ',filename,' for input'
          call psb_abort(icontxt)
          stop
        else
          write(psb_err_unit,*) 'Opened file ',trim(filename),' for input'
        end if
      else
        inp_unit=psb_inp_unit
      end if
      !
      ! input files
      call read_data(mtrx,inp_unit)            ! matrix file
      call read_data(rhs,inp_unit)             ! rhs file
      call read_data(guess,inp_unit)           ! starting guess file
      call read_data(sol,inp_unit)             ! solution file (for comparison)
      call read_data(filefmt,inp_unit)         ! format of files
      call read_data(afmt,inp_unit)            ! matrix storage format
      call read_data(part,inp_unit)            ! partition type
      ! Krylov solver data
      call read_data(solve%kmethd,inp_unit)    ! Krylov solver
      call read_data(solve%istopc,inp_unit)    ! stopping criterion
      call read_data(solve%itmax,inp_unit)     ! max num iterations
      call read_data(solve%itrace,inp_unit)    ! tracing
      call read_data(solve%irst,inp_unit)      ! restart
      call read_data(solve%eps,inp_unit)       ! tolerance
      ! preconditioner type
      call read_data(prec%descr,inp_unit)      ! verbose description of the prec
      call read_data(prec%ptype,inp_unit)      ! preconditioner type
      ! First smoother / 1-lev preconditioner
      call read_data(prec%smther,inp_unit)     ! smoother type
      call read_data(prec%jsweeps,inp_unit)    ! (pre-)smoother / 1-lev prec sweeps
      call read_data(prec%novr,inp_unit)       ! number of overlap layers
      call read_data(prec%restr,inp_unit)      ! restriction  over application of AS
      call read_data(prec%prol,inp_unit)       ! prolongation over application of AS
      call read_data(prec%solve,inp_unit)      ! local subsolver
      call read_data(prec%fill,inp_unit)       ! fill-in for incomplete LU
      call read_data(prec%thr,inp_unit)        ! threshold for ILUT
      ! Second smoother/ AMG post-smoother (if NONE ignored in main)
      call read_data(prec%smther2,inp_unit)     ! smoother type
      call read_data(prec%jsweeps2,inp_unit)    ! (post-)smoother sweeps
      call read_data(prec%novr2,inp_unit)       ! number of overlap layers
      call read_data(prec%restr2,inp_unit)      ! restriction  over application of AS
      call read_data(prec%prol2,inp_unit)       ! prolongation over application of AS
      call read_data(prec%solve2,inp_unit)      ! local subsolver
      call read_data(prec%fill2,inp_unit)       ! fill-in for incomplete LU
      call read_data(prec%thr2,inp_unit)        ! threshold for ILUT
      ! general AMG data
      call read_data(prec%mlcycle,inp_unit)     ! AMG cycle type
      call read_data(prec%outer_sweeps,inp_unit) ! number of 1lev/outer sweeps
      call read_data(prec%maxlevs,inp_unit)    ! max number of levels in AMG prec
      call read_data(prec%csize,inp_unit)       ! min size coarsest mat
      ! aggregation
      call read_data(prec%aggr_prol,inp_unit)    ! aggregation type
      call read_data(prec%par_aggr_alg,inp_unit)    ! parallel aggregation alg
      call read_data(prec%aggr_ord,inp_unit)    ! ordering for aggregation
      call read_data(prec%aggr_filter,inp_unit) ! filtering
      call read_data(prec%mncrratio,inp_unit)  ! minimum aggregation ratio
      call read_data(prec%thrvsz,inp_unit)      ! size of aggr thresh vector
      if (prec%thrvsz > 0) then
        call psb_realloc(prec%thrvsz,prec%athresv,info)
        call read_data(prec%athresv,inp_unit)   ! aggr thresh vector
      else
        read(inp_unit,*)                        ! dummy read to skip a record
      end if
      call read_data(prec%athres,inp_unit)      ! smoothed aggr thresh
      ! coasest-level solver
      call read_data(prec%csolve,inp_unit)      ! coarsest-lev solver
      call read_data(prec%csbsolve,inp_unit)    ! coarsest-lev subsolver
      call read_data(prec%cmat,inp_unit)        ! coarsest mat layout
      call read_data(prec%cfill,inp_unit)       ! fill-in for incompl LU
      call read_data(prec%cthres,inp_unit)      ! Threshold for ILUT
      call read_data(prec%cjswp,inp_unit)       ! sweeps for GS/JAC subsolver
      if (inp_unit /= psb_inp_unit) then
        close(inp_unit)
      end if     
    end if

    call psb_bcast(icontxt,mtrx)
    call psb_bcast(icontxt,rhs)
    call psb_bcast(icontxt,guess)
    call psb_bcast(icontxt,sol)
    call psb_bcast(icontxt,filefmt)
    call psb_bcast(icontxt,afmt)
    call psb_bcast(icontxt,part)

    call psb_bcast(icontxt,solve%kmethd)
    call psb_bcast(icontxt,solve%istopc)
    call psb_bcast(icontxt,solve%itmax)
    call psb_bcast(icontxt,solve%itrace)
    call psb_bcast(icontxt,solve%irst)
    call psb_bcast(icontxt,solve%eps)

    call psb_bcast(icontxt,prec%descr)
    call psb_bcast(icontxt,prec%ptype)

    ! broadcast first (pre-)smoother / 1-lev prec data
    call psb_bcast(icontxt,prec%smther)     
    call psb_bcast(icontxt,prec%jsweeps)
    call psb_bcast(icontxt,prec%novr)
    call psb_bcast(icontxt,prec%restr)
    call psb_bcast(icontxt,prec%prol)
    call psb_bcast(icontxt,prec%solve)
    call psb_bcast(icontxt,prec%fill)
    call psb_bcast(icontxt,prec%thr)
    ! broadcast second (post-)smoother 
    call psb_bcast(icontxt,prec%smther2)
    call psb_bcast(icontxt,prec%jsweeps2)
    call psb_bcast(icontxt,prec%novr2)
    call psb_bcast(icontxt,prec%restr2)
    call psb_bcast(icontxt,prec%prol2)
    call psb_bcast(icontxt,prec%solve2)
    call psb_bcast(icontxt,prec%fill2)
    call psb_bcast(icontxt,prec%thr2)
    
    ! broadcast AMG parameters
    call psb_bcast(icontxt,prec%mlcycle)
    call psb_bcast(icontxt,prec%outer_sweeps)
    call psb_bcast(icontxt,prec%maxlevs)
    
    call psb_bcast(icontxt,prec%aggr_prol)
    call psb_bcast(icontxt,prec%par_aggr_alg)
    call psb_bcast(icontxt,prec%aggr_ord)
    call psb_bcast(icontxt,prec%aggr_filter)
    call psb_bcast(icontxt,prec%mncrratio)
    call psb_bcast(ictxt,prec%thrvsz)
    if (prec%thrvsz > 0) then
      if (iam /= psb_root_) call psb_realloc(prec%thrvsz,prec%athresv,info)
      call psb_bcast(ictxt,prec%athresv)
    end if
    call psb_bcast(ictxt,prec%athres)
    
    call psb_bcast(icontxt,prec%csize)
    call psb_bcast(icontxt,prec%cmat)
    call psb_bcast(icontxt,prec%csolve)
    call psb_bcast(icontxt,prec%csbsolve)
    call psb_bcast(icontxt,prec%cfill)
    call psb_bcast(icontxt,prec%cthres)
    call psb_bcast(icontxt,prec%cjswp)


  end subroutine get_parms
  
end program mld_cf_sample
