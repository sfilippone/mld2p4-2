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

module data_input
  
  interface read_data
    module procedure read_char, read_int, read_double
  end interface read_data
  
contains

  subroutine read_char(val,file)
    character(len=*), intent(out) :: val
    integer, intent(in)           :: file
    character(len=1024) :: charbuf
    integer :: idx
    read(file,'(a)')charbuf
    charbuf = adjustl(charbuf)
    idx=index(charbuf,"!")
    read(charbuf(1:idx-1),'(a)') val
!!$    write(0,*) 'read_char got value: "',val,'"'
  end subroutine read_char
  subroutine read_int(val,file)
    integer, intent(out) :: val
    integer, intent(in)  :: file
    character(len=1024) :: charbuf
    integer :: idx
    read(file,'(a)')charbuf
    charbuf = adjustl(charbuf)
    idx=index(charbuf,"!")
    read(charbuf(1:idx-1),*) val
!!$    write(0,*) 'read_int got value: ',val
  end subroutine read_int
  subroutine read_double(val,file)
    use psb_base_mod
    real(psb_dpk_), intent(out) :: val
    integer, intent(in)         :: file
    character(len=1024) :: charbuf
    integer :: idx
    read(file,'(a)')charbuf
    charbuf = adjustl(charbuf)
    idx=index(charbuf,"!")
    read(charbuf(1:idx-1),*) val
!!$    write(0,*) 'read_double got value: ',val
  end subroutine read_double
end module data_input


program df_sample
  use psb_base_mod
  use mld_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use data_input
  implicit none


  ! input parameters
  character(len=40) :: kmethd, mtrx_file, rhs_file
  type precdata
    character(len=20)  :: descr       ! verbose description of the prec
    character(len=10)  :: prec        ! overall prectype
    integer            :: novr        ! number of overlap layers
    character(len=16)  :: restr       ! restriction  over application of as
    character(len=16)  :: prol        ! prolongation over application of as
    character(len=16)  :: solve      ! Factorization type: ILU, SuperLU, UMFPACK. 
    integer            :: fill1       ! Fill-in for factorization 1
    real(psb_dpk_)     :: thr1        ! Threshold for fact. 1 ILU(T)
    integer            :: nlev        ! Number of levels in multilevel prec. 
    character(len=16)  :: aggrkind    ! smoothed/raw aggregatin
    character(len=16)  :: aggr_alg    ! local or global aggregation
    character(len=16)  :: mltype      ! additive or multiplicative 2nd level prec
    character(len=16)  :: smthpos     ! side: pre, post, both smoothing
    character(len=16)  :: cmat        ! coarse mat
    character(len=16)  :: csolve       ! Factorization type: ILU, SuperLU, UMFPACK. 
    integer            :: cfill       ! Fill-in for factorization 1
    real(psb_dpk_)     :: cthres      ! Threshold for fact. 1 ILU(T)
    integer            :: cjswp       ! Jacobi sweeps
    real(psb_dpk_)     :: omega       ! smoother omega
  end type precdata
  type(precdata)        :: prec_choice

  ! sparse matrices
  type(psb_dspmat_type) :: a, aux_a

  ! preconditioner data
  type(mld_dprec_type)  :: prec

  ! dense matrices
  real(kind(1.d0)), allocatable, target ::  aux_b(:,:), d(:)
  real(kind(1.d0)), allocatable , save  :: b_col(:), x_col(:), r_col(:), &
       & x_col_glob(:), r_col_glob(:)
  real(kind(1.d0)), pointer  :: b_col_glob(:)

  ! communications data structure
  type(psb_desc_type):: desc_a

  integer            :: ictxt, iam, np

  ! solver paramters
  integer            :: iter, itmax, ierr, itrace, ircode, ipart,&
       & methd, istopc, irst,amatsize,precsize,descsize, nlv
  real(kind(1.d0))   :: err, eps

  character(len=5)   :: afmt
  character(len=20)  :: name
  integer   :: iparm(20)

  ! other variables
  integer            :: i,info,j,m_problem
  integer            :: internal, m,ii,nnzero
  real(kind(1.d0)) :: t1, t2, tprec, r_amax, b_amax,&
       &scale,resmx,resmxp
  integer :: nrhs, nrow, n_row, dim, nv, ne
  integer, allocatable :: ivg(:), ipv(:)


  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif


  name='df_sample'
  if(psb_get_errstatus() /= 0) goto 9999
  info=0
  call psb_set_errverbosity(2)
  !
  !  get parameters
  !
  call get_parms(ictxt,mtrx_file,rhs_file,kmethd,&
       & prec_choice,ipart,afmt,istopc,itmax,itrace,irst,eps)

  call psb_barrier(ictxt)
  t1 = psb_wtime()  
  ! read the input matrix to be processed and (possibly) the rhs 
  nrhs = 1

  if (iam==psb_root_) then
    call read_mat(mtrx_file, aux_a, ictxt)

    m_problem = aux_a%m
    call psb_bcast(ictxt,m_problem)

    if(rhs_file /= 'NONE') then
      !  reading an rhs
      call read_rhs(rhs_file,aux_b,ictxt)
    end if

    if (psb_size(aux_b,dim=1)==m_problem) then
      ! if any rhs were present, broadcast the first one
      write(0,'("Ok, got an rhs ")')
      b_col_glob =>aux_b(:,1)
    else
      write(*,'("Generating an rhs...")')
      write(*,'(" ")')
      call psb_realloc(m_problem,1,aux_b,ircode)
      if (ircode /= 0) then
        call psb_errpush(4000,name)
        goto 9999
      endif

      b_col_glob => aux_b(:,1)
      do i=1, m_problem
        b_col_glob(i) = 1.d0
      enddo
    endif
    call psb_bcast(ictxt,b_col_glob(1:m_problem))
  else
    call psb_bcast(ictxt,m_problem)
    call psb_realloc(m_problem,1,aux_b,ircode)
    if (ircode /= 0) then
      call psb_errpush(4000,name)
      goto 9999
    endif
    b_col_glob =>aux_b(:,1)
    call psb_bcast(ictxt,b_col_glob(1:m_problem)) 
  end if

  ! switch over different partition types
  if (ipart == 0) then 
    call psb_barrier(ictxt)
    if (iam==psb_root_) write(*,'("Partition type: block")')
    allocate(ivg(m_problem),ipv(np))
    do i=1,m_problem
      call part_block(i,m_problem,np,ipv,nv)
      ivg(i) = ipv(1)
    enddo
    call psb_matdist(aux_a, a, ivg, ictxt, &
         & desc_a,b_col_glob,b_col,info,fmt=afmt)
  else if (ipart == 2) then 
    if (iam==psb_root_) then 
      write(*,'("Partition type: graph")')
      write(*,'(" ")')
      !      write(0,'("Build type: graph")')
      call build_mtpart(aux_a%m,aux_a%fida,aux_a%ia1,aux_a%ia2,np)
    endif
    call psb_barrier(ictxt)
    call distr_mtpart(psb_root_,ictxt)
    call getv_mtpart(ivg)
    call psb_matdist(aux_a, a, ivg, ictxt, &
         & desc_a,b_col_glob,b_col,info,fmt=afmt)
  else 
    if (iam==psb_root_) write(*,'("Partition type: block")')
    call psb_matdist(aux_a, a, part_block, ictxt, &
         & desc_a,b_col_glob,b_col,info,fmt=afmt)
  end if

  call psb_geall(x_col,desc_a,info)
  x_col(:) =0.0
  call psb_geasb(x_col,desc_a,info)
  call psb_geall(r_col,desc_a,info)
  r_col(:) =0.0
  call psb_geasb(r_col,desc_a,info)
  t2 = psb_wtime() - t1


  call psb_amx(ictxt, t2)

  if (iam==psb_root_) then
    write(*,'(" ")')
    write(*,'("Time to read and partition matrix : ",es10.4)')t2
    write(*,'(" ")')
    write(*,*) 'Preconditioner: ',prec_choice%descr
  end if

  ! 

  if (toupper(prec_choice%prec) =='ML') then 
    nlv = prec_choice%nlev
  else
    nlv = 1
  end if
  call mld_precinit(prec,prec_choice%prec,info,nlev=nlv)
  call mld_precset(prec,mld_n_ovr_,prec_choice%novr,info)
  call mld_precset(prec,mld_sub_restr_,prec_choice%restr,info)
  call mld_precset(prec,mld_sub_prol_,prec_choice%prol,info)
  call mld_precset(prec,mld_sub_solve_,prec_choice%solve,info)
  call mld_precset(prec,mld_sub_fill_in_,prec_choice%fill1,info)
  call mld_precset(prec,mld_fact_thrs_,prec_choice%thr1,info)
  if (toupper(prec_choice%prec) =='ML') then 
    call mld_precset(prec,mld_aggr_kind_,prec_choice%aggrkind,info)
    call mld_precset(prec,mld_aggr_alg_,prec_choice%aggr_alg,info)
    call mld_precset(prec,mld_ml_type_,prec_choice%mltype,info)
    call mld_precset(prec,mld_ml_type_,prec_choice%mltype,info)
    call mld_precset(prec,mld_smooth_pos_,prec_choice%smthpos,info)
    call mld_precset(prec,mld_coarse_mat_,prec_choice%cmat,info)
    call mld_precset(prec,mld_coarse_solve_,prec_choice%csolve,info)
    call mld_precset(prec,mld_sub_fill_in_,prec_choice%cfill,info,ilev=nlv)
    call mld_precset(prec,mld_fact_thrs_,prec_choice%cthres,info,ilev=nlv)
    call mld_precset(prec,mld_smooth_sweeps_,prec_choice%cjswp,info,ilev=nlv)
    call mld_precset(prec,mld_smooth_sweeps_,prec_choice%cjswp,info,ilev=nlv)
    if (prec_choice%omega>=0.0) then 
      call mld_precset(prec,mld_aggr_damp_,prec_choice%omega,info,ilev=nlv)
    end if
  end if

  ! building the preconditioner
  t1 = psb_wtime()
  call mld_precbld(a,desc_a,prec,info)
  tprec = psb_wtime()-t1
  if (info /= 0) then
    call psb_errpush(4010,name,a_err='psb_precbld')
    goto 9999
  end if


  call psb_amx(ictxt, tprec)

  if(iam==psb_root_) then
    write(*,'("Preconditioner time: ",es10.4)')tprec
    write(*,'(" ")')
  end if

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
  call psb_genrm2s(resmx,r_col,desc_a,info)
  call psb_geamaxs(resmxp,r_col,desc_a,info)

  amatsize = psb_sizeof(a)
  descsize = psb_sizeof(desc_a)
  precsize = mld_sizeof(prec)
  call psb_sum(ictxt,amatsize)
  call psb_sum(ictxt,descsize)
  call psb_sum(ictxt,precsize)
  if (iam==psb_root_) then 
    call mld_prec_descr(6,prec)
    write(*,'("Matrix: ",a)')mtrx_file
    write(*,'("Computed solution on ",i8," processors")')np
    write(*,'("Iterations to convergence: ",i6)')iter
    write(*,'("Error estimate on exit: ",f7.2)')err
    write(*,'("Time to buil prec.   : ",es10.4)')tprec
    write(*,'("Time to solve matrix : ",es10.4)')t2
    write(*,'("Time per iteration   : ",es10.4)')t2/(iter)
    write(*,'("Total time           : ",es10.4)')t2+tprec
    write(*,'("Residual norm 2   = ",es10.4)')resmx
    write(*,'("Residual norm inf = ",es10.4)')resmxp
    write(*,'("Total memory occupation for A:      ",i10)')amatsize
    write(*,'("Total memory occupation for DESC_A: ",i10)')descsize
    write(*,'("Total memory occupation for PREC:   ",i10)')precsize
  end if

  allocate(x_col_glob(m_problem),r_col_glob(m_problem),stat=ierr)
  if (ierr /= 0) then 
    write(0,*) 'allocation error: no data collection'
  else
    call psb_gather(x_col_glob,x_col,desc_a,info,root=psb_root_)
    call psb_gather(r_col_glob,r_col,desc_a,info,root=psb_root_)
    if (iam==psb_root_) then
      write(0,'(" ")')
      write(0,'("Saving x on file")')
      write(20,*) 'matrix: ',mtrx_file
      write(20,*) 'computed solution on ',np,' processors.'
      write(20,*) 'iterations to convergence: ',iter
      write(20,*) 'error estimate (infinity norm) on exit:', &
           & ' ||r||/(||a||||x||+||b||) = ',err
      write(20,*) 'max residual = ',resmx, resmxp
      write(20,'(a8,4(2x,a20))') 'I','X(I)','R(I)','B(I)'
      do i=1,m_problem
        write(20,998) i,x_col_glob(i),r_col_glob(i),b_col_glob(i)
      enddo
    end if
  end if
998 format(i8,4(2x,g20.14))
993 format(i6,4(1x,e12.6))


  call psb_gefree(b_col, desc_a,info)
  call psb_gefree(x_col, desc_a,info)
  call psb_spfree(a, desc_a,info)
  call mld_precfree(prec,info)
  call psb_cdfree(desc_a,info)

9999 continue
  if(info /= 0) then
    call psb_error(ictxt)
  end if
  call psb_exit(ictxt)
  stop

contains
  !
  ! get iteration parameters from standard input
  !
  subroutine  get_parms(icontxt,mtrx,rhs,kmethd,&
       & prec, ipart,afmt,istopc,itmax,itrace,irst,eps)

    use psb_base_mod
    implicit none

    integer             :: icontxt
    character(len=*)    :: kmethd, mtrx, rhs, afmt
    type(precdata)      :: prec
    integer             :: iret, istopc,itmax,itrace, ipart, irst
    real(psb_dpk_)      :: eps, omega,thr1,thr2
    integer             :: iam, nm, np, i

    call psb_info(icontxt,iam,np)

    if (iam==psb_root_) then
      ! read input parameters
      call read_data(mtrx,5)
      call read_data(rhs,5)
      call read_data(kmethd,5)
      call read_data(afmt,5)
      call read_data(ipart,5)
      call read_data(istopc,5)
      call read_data(itmax,5)
      call read_data(itrace,5)
      call read_data(irst,5)
      call read_data(eps,5)
      call read_data(prec%descr,5)       ! verbose description of the prec
      call read_data(prec%prec,5)        ! overall prectype
      call read_data(prec%novr,5)        ! number of overlap layers
      call read_data(prec%restr,5)       ! restriction  over application of as
      call read_data(prec%prol,5)        ! prolongation over application of as
      call read_data(prec%solve,5)       ! Factorization type: ILU, SuperLU, UMFPACK. 
      call read_data(prec%fill1,5)       ! Fill-in for factorization 1
      call read_data(prec%thr1,5)        ! Threshold for fact. 1 ILU(T)
      if (toupper(prec%prec) == 'ML') then 
        call read_data(prec%nlev,5)        ! Number of levels in multilevel prec. 
        call read_data(prec%aggrkind,5)    ! smoothed/raw aggregatin
        call read_data(prec%aggr_alg,5)    ! local or global aggregation
        call read_data(prec%mltype,5)      ! additive or multiplicative 2nd level prec
        call read_data(prec%smthpos,5)     ! side: pre, post, both smoothing
        call read_data(prec%cmat,5)        ! coarse mat
        call read_data(prec%csolve,5)      ! Factorization type: ILU, SuperLU, UMFPACK. 
        call read_data(prec%cfill,5)       ! Fill-in for factorization 1
        call read_data(prec%cthres,5)      ! Threshold for fact. 1 ILU(T)
        call read_data(prec%cjswp,5)       ! Jacobi sweeps
        call read_data(prec%omega,5)       ! smoother omega
      end if
    end if

    call psb_bcast(icontxt,mtrx)
    call psb_bcast(icontxt,rhs)
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
    call psb_bcast(icontxt,prec%fill1)       ! Fill-in for factorization 1
    call psb_bcast(icontxt,prec%thr1)        ! Threshold for fact. 1 ILU(T)
    if (toupper(prec%prec) == 'ML') then 
      call psb_bcast(icontxt,prec%nlev)        ! Number of levels in multilevel prec. 
      call psb_bcast(icontxt,prec%aggrkind)    ! smoothed/raw aggregatin
      call psb_bcast(icontxt,prec%aggr_alg)    ! local or global aggregation
      call psb_bcast(icontxt,prec%mltype)      ! additive or multiplicative 2nd level prec
      call psb_bcast(icontxt,prec%smthpos)     ! side: pre, post, both smoothing
      call psb_bcast(icontxt,prec%cmat)        ! coarse mat
      call psb_bcast(icontxt,prec%csolve)      ! Factorization type: ILU, SuperLU, UMFPACK. 
      call psb_bcast(icontxt,prec%cfill)       ! Fill-in for factorization 1
      call psb_bcast(icontxt,prec%cthres)      ! Threshold for fact. 1 ILU(T)
      call psb_bcast(icontxt,prec%cjswp)       ! Jacobi sweeps
      call psb_bcast(icontxt,prec%omega)       ! smoother omega
    end if

  end subroutine get_parms
  subroutine pr_usage(iout)
    integer iout
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
  




