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
program df_sample
  use psb_base_mod
  use mld_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use data_input
  implicit none


  ! input parameters
  character(len=40) :: kmethd, mtrx_file, rhs_file
  character(len=2)  :: filefmt
  type precdata
    character(len=20)  :: descr       ! verbose description of the prec
    character(len=10)  :: prec        ! overall prectype
    integer(psb_ipk_)  :: novr        ! number of overlap layers
    integer(psb_ipk_)  :: jsweeps     ! Jacobi/smoother sweeps
    character(len=16)  :: restr       ! restriction over application of AS
    character(len=16)  :: prol        ! prolongation over application of AS
    character(len=16)  :: solve       ! factorization type: ILU, SuperLU, UMFPACK 
    integer(psb_ipk_)  :: fill        ! fillin for factorization 
    real(psb_dpk_)     :: thr         ! threshold for fact.  ILU(T)
    character(len=16)  :: smther      ! Smoother                            
    integer(psb_ipk_)  :: nlev        ! number of levels in multilevel prec. 
    character(len=16)  :: aggrkind    ! smoothed, raw aggregation
    character(len=16)  :: aggr_alg    ! aggregation algorithm (currently only decoupled)
    character(len=16)  :: mltype      ! additive or multiplicative multi-level prec
    character(len=16)  :: smthpos     ! side: pre, post, both smoothing
    character(len=16)  :: cmat        ! coarse mat: distributed, replicated
    character(len=16)  :: csolve      ! coarse solver: bjac, umf, slu, sludist
    character(len=16)  :: csbsolve    ! coarse subsolver: ILU, ILU(T), SuperLU, UMFPACK 
    integer(psb_ipk_)  :: cfill       ! fillin for coarse factorization 
    real(psb_dpk_)     :: cthres      ! threshold for coarse fact.  ILU(T)
    integer(psb_ipk_)  :: cjswp       ! block-Jacobi sweeps
    real(psb_dpk_)     :: athres      ! smoothed aggregation threshold
    real(psb_dpk_)     :: ascale      ! smoothed aggregation scale factor
  end type precdata
  type(precdata)       :: prec_choice

  ! sparse matrices
  type(psb_dspmat_type) :: a, aux_a

  ! preconditioner data
  type(mld_dprec_type)  :: prec
  ! dense matrices
  real(psb_dpk_), allocatable, target ::  aux_b(:,:), d(:)
  real(psb_dpk_), allocatable , save  :: x_col_glob(:), r_col_glob(:)
  real(psb_dpk_), pointer  :: b_col_glob(:)
  type(psb_d_vect_type)    :: b_col, x_col, r_col

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
  integer(psb_ipk_) :: i,info,j,m_problem, lbw,ubw,prf
  integer(psb_ipk_) :: internal, m,ii,nnzero
  real(psb_dpk_)    :: t1, t2, tprec
  real(psb_dpk_)    :: r_amax, b_amax, scale,resmx,resmxp
  integer(psb_ipk_) :: nrhs, nrow, n_row, dim, nv, ne
  integer(psb_ipk_), allocatable :: ivg(:), ipv(:),perm(:)

  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif


  name='df_sample'
  if(psb_get_errstatus() /= 0) goto 9999
  info=psb_success_
  call psb_set_errverbosity(itwo)
  call psb_cd_set_large_threshold(itwo)
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
  call get_parms(ictxt,mtrx_file,rhs_file,filefmt,kmethd,&
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
    call psb_realloc(m_problem,1,aux_b,ircode)
    if (ircode /= 0) then
      call psb_errpush(psb_err_alloc_dealloc_,name)
      goto 9999
    endif
    b_col_glob =>aux_b(:,1)
  end if

  !
  ! Renumbering?   
  !
  if (iam==psb_root_) then 
    renum='NONE'
    call psb_cmp_bwpf(aux_a,lbw,ubw,prf,info)
    write(*,*) 'Bandwidth and profile (original): ',lbw,ubw,prf
    write(*,*) 'Renumbering algorithm           : ',psb_toupper(renum)
    if (trim(psb_toupper(renum))/='NONE') then 
      call psb_mat_renum(renum,aux_a,info,perm=perm)
      if (info /= 0) then 
        write(0,*) 'Error from RENUM',info
        goto 9999
      end if
      
      call psb_gelp('N',perm(1:m_problem),&
           & b_col_glob(1:m_problem),info)
      call psb_cmp_bwpf(aux_a,lbw,ubw,prf,info)
    end if
  
    write(*,*) 'Bandwidth and profile (renumbrd):',lbw,ubw,prf
  end if
  
  call psb_bcast(ictxt,b_col_glob(1:m_problem)) 

  call aux_a%clean_zeros(info)
  
  ! switch over different partition types
  if (ipart == 0) then 
    call psb_barrier(ictxt)
    if (iam == psb_root_) write(psb_out_unit,'("Partition type: block")')
    allocate(ivg(m_problem),ipv(np))
    do i=1,m_problem
      call part_block(i,m_problem,np,ipv,nv)
      ivg(i) = ipv(1)
    enddo
    call psb_matdist(aux_a, a, ictxt, &
         & desc_a,info,b_glob=b_col_glob,b=b_col,fmt=afmt,v=ivg)
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
    call psb_matdist(aux_a, a, ictxt, &
         & desc_a,info,b_glob=b_col_glob,b=b_col,fmt=afmt,v=ivg)
  else 
    if (iam == psb_root_) write(psb_out_unit,'("Partition type: block")')
    call psb_matdist(aux_a, a,  ictxt, &
         & desc_a,info,b_glob=b_col_glob,b=b_col,fmt=afmt,parts=part_block)
  end if

  call psb_geall(x_col,desc_a,info)
  call x_col%set(dzero)
  call psb_geasb(x_col,desc_a,info)
  call psb_geall(r_col,desc_a,info)
  call r_col%set(dzero)
  call psb_geasb(r_col,desc_a,info)
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
    nlv = prec_choice%nlev
    call mld_precinit(prec,prec_choice%prec,info,nlev=nlv)
    
    call mld_precset(prec,mld_smoother_type_,   prec_choice%smther,  info)
    call mld_precset(prec,mld_smoother_sweeps_, prec_choice%jsweeps, info)
    call mld_precset(prec,mld_sub_ovr_,         prec_choice%novr,    info)
    call mld_precset(prec,mld_sub_restr_,       prec_choice%restr,   info)
    call mld_precset(prec,mld_sub_prol_,        prec_choice%prol,    info)
    call mld_precset(prec,mld_sub_solve_,       prec_choice%solve,   info)
    call mld_precset(prec,mld_sub_fillin_,      prec_choice%fill,   info)
    call mld_precset(prec,mld_sub_iluthrs_,     prec_choice%thr,    info)
    call mld_precset(prec,mld_aggr_kind_,       prec_choice%aggrkind,info)
    call mld_precset(prec,mld_aggr_alg_,        prec_choice%aggr_alg,info)
    call mld_precset(prec,mld_ml_type_,         prec_choice%mltype,  info)
    call mld_precset(prec,mld_smoother_pos_,    prec_choice%smthpos, info)
    call mld_precset(prec,mld_aggr_scale_,      prec_choice%ascale,  info)
    call mld_precset(prec,mld_aggr_thresh_,     prec_choice%athres,  info)
    call mld_precset(prec,mld_coarse_solve_,    prec_choice%csolve,  info)
    call mld_precset(prec,mld_ml_type_,'MULT',info)
    call mld_precset(prec,mld_smoother_pos_,'TWOSIDE',info)
    call mld_precset(prec,mld_coarse_sweeps_,4,info)
    call mld_precset(prec,mld_coarse_subsolve_, prec_choice%csbsolve,info)
    call mld_precset(prec,mld_coarse_mat_,      prec_choice%cmat,    info)
    call mld_precset(prec,mld_coarse_fillin_,   prec_choice%cfill,   info)
    call mld_precset(prec,mld_coarse_iluthrs_,  prec_choice%cthres,  info)
    call mld_precset(prec,mld_coarse_sweeps_,   prec_choice%cjswp,   info)
    
  else
    nlv = 1
    call mld_precinit(prec,prec_choice%prec,info)
    if (psb_toupper(prec_choice%prec) /= 'NONE') then 
      call mld_precset(prec,mld_smoother_sweeps_, prec_choice%jsweeps, info)
      call mld_precset(prec,mld_sub_ovr_,         prec_choice%novr,    info)
      call mld_precset(prec,mld_sub_restr_,       prec_choice%restr,   info)
      call mld_precset(prec,mld_sub_prol_,        prec_choice%prol,    info)
      call mld_precset(prec,mld_sub_solve_,       prec_choice%solve,   info)
      call mld_precset(prec,mld_sub_fillin_,      prec_choice%fill,   info)
      call mld_precset(prec,mld_sub_iluthrs_,     prec_choice%thr,    info)
    end if
  end if
  ! building the preconditioner
  t1 = psb_wtime()
  call mld_precbld(a,desc_a,prec,info)
  tprec = psb_wtime()-t1
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_precbld')
    goto 9999
  end if

  call psb_amx(ictxt, tprec)

  if(iam == psb_root_) then
    write(psb_out_unit,'("Preconditioner time: ",es12.5)')tprec
    write(psb_out_unit,'(" ")')
  end if

  write(fprefix,'(a,i3.3,a,i3.3)') 'proc-',iam,'-',np
  call prec%dump(info,head='Pressure test case',ac=.true.)
    

  
  iparm = 0
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call psb_krylov(kmethd,a,prec,b_col,x_col,eps,desc_a,info,& 
       & itmax=itmax,iter=iter,err=err,itrace=itrace,istop=istopc,irst=irst)     
  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1

  call psb_amx(ictxt,t2)
  call psb_geaxpby(done,b_col,dzero,r_col,desc_a,info)
  call psb_spmm(-done,a,x_col,done,r_col,desc_a,info)
  resmx  = psb_genrm2(r_col,desc_a,info)
  resmxp = psb_geamax(r_col,desc_a,info)
  
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
    write(psb_out_unit,'("Iterations to convergence : ",i6)')iter
    write(psb_out_unit,'("Error estimate on exit    : ",es12.5)')err
    write(psb_out_unit,'("Time to buil prec.        : ",es12.5)')tprec
    write(psb_out_unit,'("Time to solve matrix      : ",es12.5)')t2
    write(psb_out_unit,'("Time per iteration        : ",es12.5)')t2/(iter)
    write(psb_out_unit,'("Total time                : ",es12.5)')t2+tprec
    write(psb_out_unit,'("Residual norm 2           : ",es12.5)')resmx
    write(psb_out_unit,'("Residual norm inf         : ",es12.5)')resmxp
    write(psb_out_unit,'("Total memory occupation for A      : ",i12)')amatsize
    write(psb_out_unit,'("Total memory occupation for DESC_A : ",i12)')descsize
    write(psb_out_unit,'("Total memory occupation for PREC   : ",i12)')precsize
    write(psb_out_unit,'("Storage format for A               : ",a  )')a%get_fmt()
    write(psb_out_unit,'("Storage format for DESC_A          : ",a  )')desc_a%get_fmt()
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
  subroutine  get_parms(icontxt,mtrx,rhs,filefmt,kmethd,&
       & prec, ipart,afmt,istopc,itmax,itrace,irst,eps)

    use psb_base_mod
    implicit none

    integer(psb_ipk_)   :: icontxt
    character(len=*)    :: kmethd, mtrx, rhs, afmt,filefmt
    type(precdata)      :: prec
    real(psb_dpk_)      :: eps
    integer(psb_ipk_)   :: iret, istopc,itmax,itrace, ipart, irst
    integer(psb_ipk_)   :: iam, nm, np, i

    call psb_info(icontxt,iam,np)

    if (iam == psb_root_) then
      ! read input parameters
      call read_data(mtrx,psb_inp_unit)
      call read_data(rhs,psb_inp_unit)
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
      call read_data(prec%fill,psb_inp_unit)       ! Fill-in for factorization 
      call read_data(prec%thr,psb_inp_unit)        ! Threshold for fact.  ILU(T)
      call read_data(prec%jsweeps,psb_inp_unit)    ! Jacobi sweeps for PJAC
      if (psb_toupper(prec%prec) == 'ML') then 
        call read_data(prec%nlev,psb_inp_unit)     ! Number of levels in multilevel prec. 
        call read_data(prec%smther,psb_inp_unit)   ! Smoother type.
        call read_data(prec%aggrkind,psb_inp_unit) ! smoothed/raw aggregatin
        call read_data(prec%aggr_alg,psb_inp_unit) ! local or global aggregation
        call read_data(prec%mltype,psb_inp_unit)   ! additive or multiplicative 2nd level prec
        call read_data(prec%smthpos,psb_inp_unit)  ! side: pre, post, both smoothing
        call read_data(prec%cmat,psb_inp_unit)     ! coarse mat
        call read_data(prec%csolve,psb_inp_unit)   ! Factorization type: BJAC, SuperLU, UMFPACK. 
        call read_data(prec%csbsolve,psb_inp_unit) ! Factorization type: ILU, SuperLU, UMFPACK. 
        call read_data(prec%cfill,psb_inp_unit)    ! Fill-in for factorization 
        call read_data(prec%cthres,psb_inp_unit)   ! Threshold for fact.  ILU(T)
        call read_data(prec%cjswp,psb_inp_unit)    ! Jacobi sweeps
        call read_data(prec%athres,psb_inp_unit)   ! smoother aggr thresh
        call read_data(prec%ascale,psb_inp_unit)   ! smoother aggr thresh
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
      call psb_bcast(icontxt,prec%aggrkind)    ! smoothed/raw aggregatin
      call psb_bcast(icontxt,prec%aggr_alg)    ! local or global aggregation
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
end program df_sample
