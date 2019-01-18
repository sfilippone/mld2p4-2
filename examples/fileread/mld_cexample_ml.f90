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
! File: mld_cexample_ml.f90
!
! This sample program solves a linear system by using CG coupled with
! one of the following multi-level preconditioner, as explained in Section 6.1
! of the MLD2P4 User's and Reference Guide:
!
! - choice = 1, the default multi-level preconditioner solver, i.e., 
! V-cycle with basic smoothed aggregation, 1 hybrid forward/backward
! GS sweep as pre/post-smoother and UMFPACK as coarsest-level
! solver (Sec. 5.1, Fig. 2)
!
! - choice = 2, a V-cycle preconditioner with 1 block-Jacobi sweep
! (with ILU(0) on the blocks) as pre- and post-smoother, and 8 block-Jacobi
! sweeps (with ILU(0) on the blocks) as coarsest-level solver (Sec. 5.1, Fig. 3)
!
! - choice = 3, a W-cycle preconditioner with 2 hybrid forward/backward
!  GS sweeps as pre/post-smoother, a distributed coarsest
! matrix, and MUMPS as coarsest-level solver (Sec. 5.1, Fig. 4)
!
! The matrix and the rhs are read from files (if an rhs is not available, the
! unit rhs is set).
!
program mld_cexample_ml
  use psb_base_mod
  use mld_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use data_input

  implicit none

  ! input file parameters
  character(len=40) :: mtrx_file, rhs_file
  character(len=2)  :: filefmt

  ! sparse matrices
  type(psb_cspmat_type) :: A, aux_A

  ! descriptor of sparse matrices
  type(psb_desc_type):: desc_A

  ! preconditioner
  type(mld_cprec_type)  :: P

  ! right-hand side, solution and residual vectors
  type(psb_c_vect_type)  :: b, x, r
  complex(psb_spk_), allocatable , save  :: x_glob(:), r_glob(:)
  complex(psb_spk_), allocatable, target ::  aux_b(:,:)
  complex(psb_spk_), pointer  :: b_glob(:)

  ! solver and preconditioner parameters
  real(psb_spk_)   :: tol, err
  integer          :: itmax, iter, istop
  integer          :: nlev

  ! parallel environment parameters
  integer            :: ictxt, iam, np

  ! other variables
  integer              :: choice       
  integer              :: i,info,j,m_problem
  integer(psb_long_int_k_) :: amatsize, precsize, descsize
  integer              :: ierr, ircode
  real(psb_spk_)       :: resmx, resmxp
  real(psb_dpk_)       :: t1, t2, tprec
  character(len=20)    :: name, kmethod
  integer, parameter :: iunit=12

  ! initialize the parallel environment

  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif

  name='mld_cexample_ml'
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

  call get_parms(ictxt,mtrx_file,rhs_file,filefmt,choice,itmax,tol)

  call psb_barrier(ictxt)
  t1 = psb_wtime()  

  ! read and assemble the matrix A and the right-hand side b
  ! using PSBLAS routines for sparse matrix / vector management

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
      if (info == psb_success_) call aux_a%cscnv(info,type='csr')
    case ('HB')
      ! For Harwell-Boeing we have a single file which may or may not
      ! contain an RHS.
      call hb_read(aux_a,info,iunit=iunit,b=aux_b,filename=mtrx_file)

    case default
      info = -1 
      write(0,*) 'Wrong choice for fileformat ', filefmt
    end select
    if (info /= psb_success_) then
      write(0,*) 'Error while reading input matrix '
      call psb_abort(ictxt)
    end if

    m_problem = aux_a%get_nrows()
    call psb_bcast(ictxt,m_problem)

    ! At this point aux_b may still be unallocated
    if (psb_size(aux_b,1) == m_problem) then
      ! if any rhs were present, broadcast the first one
      write(0,'("Ok, got an rhs ")')
      b_glob =>aux_b(:,1)
    else
      write(*,'("Generating an rhs...")')
      write(*,'(" ")')
      call psb_realloc(m_problem,1,aux_b,ircode)
      if (ircode /= 0) then
        call psb_errpush(psb_err_alloc_dealloc_,name)
        goto 9999
      endif

      b_glob => aux_b(:,1)
      do i=1, m_problem
        b_glob(i) = 1.d0
      enddo
    endif
  else
    call psb_bcast(ictxt,m_problem)
  end if

  call psb_barrier(ictxt)
  if (iam == psb_root_) write(*,'("Partition type: block")')
  call psb_matdist(aux_A, A, ictxt, desc_A,info,parts=part_block)
  call psb_scatter(b_glob,b,desc_a,info,root=psb_root_)

  t2 = psb_wtime() - t1

  call psb_amx(ictxt, t2)

  if (iam == psb_root_) then
    write(*,'(" ")')
    write(*,'("Time to read and partition matrix : ",es12.5)')t2
    write(*,'(" ")')
  end if

  select case(choice)

  case(1)

    ! initialize the default multi-level preconditioner, i.e. V-cycle
    ! with basic smoothed aggregation, 1 hybrid forward/backward
    ! GS sweep as pre/post-smoother and UMFPACK as coarsest-level
    ! solver

    call P%init(ictxt,'ML',info)
    kmethod = 'CG'

  case(2)

    ! initialize a V-cycle preconditioner with 1 block-Jacobi sweep (with
    ! ILU(0) on the blocks) as pre- and post-smoother, and 8  block-Jacobi
    ! sweeps (with ILU(0) on the blocks) as coarsest-level solver

    call P%init('ML',info)
    call P%set('SMOOTHER_TYPE','BJAC',info)
    call P%set('COARSE_SOLVE','BJAC',info)
    call P%set('COARSE_SWEEPS',8,info)
    kmethod = 'CG'
    
  case(3)

   ! initialize a W-cycle preconditioner with 2 hybrid forward/backward
   ! GS sweeps as pre/post-smoother, a distributed coarsest
   ! matrix, and MUMPS as coarsest-level solver

    call P%init('ML',info)
    call P%set('ML_CYCLE','WCYCLE',info)
    call P%set('SMOOTHER_SWEEPS',2,info)
    call P%set('COARSE_SOLVE','MUMPS',info)
    call P%set('COARSE_MAT','DIST',info)
    kmethod = 'CG'
  end select

  ! build the preconditioner

  call psb_barrier(ictxt)
  t1 = psb_wtime()

  ! build the preconditioner
  call P%hierarchy_build(A,desc_A,info)
  call P%smoothers_build(A,desc_A,info)

  tprec = psb_wtime()-t1
  call psb_amx(ictxt, tprec)

  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_precbld')
    goto 9999
  end if

  ! set the initial guess
  call psb_geall(x,desc_A,info)
  call x%zero()
  call psb_geasb(x,desc_A,info)

  ! solve Ax=b with preconditioned CG

  call psb_barrier(ictxt)
  t1 = psb_wtime()

  call psb_krylov(kmethod,A,P,b,x,tol,desc_A,info,itmax,iter,err,itrace=1,istop=2)

  t2 = psb_wtime() - t1
  call psb_amx(ictxt,t2)

  call psb_geasb(r,desc_A,info,scratch=.true.)
  call psb_geaxpby(cone,b,czero,r,desc_A,info)
  call psb_spmm(-cone,A,x,cone,r,desc_A,info)
  resmx =  psb_genrm2(r,desc_A,info)
  resmxp = psb_geamax(r,desc_A,info)

  amatsize = a%sizeof()
  descsize = desc_a%sizeof()
  precsize = p%sizeof()
  call psb_sum(ictxt,amatsize)
  call psb_sum(ictxt,descsize)
  call psb_sum(ictxt,precsize)

  call P%descr(info)

  if (iam == psb_root_) then
    write(*,'(" ")')
    write(*,'("Matrix: ",A)')mtrx_file
    write(*,'("Computed solution on ",i8," processors")')np
    write(*,'("Krylov method             : ",a)') kmethod
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

  call psb_gather(x_glob,x,desc_a,info,root=psb_root_)
  if (info == psb_success_) &
       & call psb_gather(r_glob,r,desc_a,info,root=psb_root_)
  if (info /= psb_success_) goto 9999
  if (iam == psb_root_) then
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
      write(20,998) i,x_glob(i),r_glob(i),b_glob(i)
    enddo
  end if
998 format(i8,4(2x,g20.14))
993 format(i6,4(1x,e12.6))

  ! deallocate the data structures

  call psb_gefree(b, desc_A,info)
  call psb_gefree(x, desc_A,info)
  call psb_spfree(A, desc_A,info)
  call P%free(info)
  call psb_cdfree(desc_A,info)
  call psb_exit(ictxt)
  stop

9999 continue
  call psb_error(ictxt)

contains
  !
  ! get parameters from standard input
  !
  subroutine get_parms(ictxt,mtrx,rhs,filefmt,choice,itmax,tol)

    implicit none

    integer             :: ictxt, choice, itmax
    real(psb_spk_)      :: tol
    character(len=*)    :: mtrx, rhs,filefmt
    integer             :: iam, np, inp_unit
    character(len=1024)   :: filename

    call psb_info(ictxt,iam,np)

    if (iam == psb_root_) then
      if (command_argument_count()>0) then
        call get_command_argument(1,filename)
        inp_unit = 30
        open(inp_unit,file=filename,action='read',iostat=info)
        if (info /= 0) then
          write(psb_err_unit,*) 'Could not open file ',filename,' for input'
          call psb_abort(ictxt)
          stop
        else
          write(psb_err_unit,*) 'Opened file ',trim(filename),' for input'
        end if
      else
        inp_unit=psb_inp_unit
      end if
      ! read input parameters
      call read_data(mtrx,inp_unit)
      call read_data(rhs,inp_unit)
      call read_data(filefmt,inp_unit)
      call read_data(choice,inp_unit)
      call read_data(itmax,inp_unit)
      call read_data(tol,inp_unit)
      if (inp_unit /= psb_inp_unit) then
        close(inp_unit)
      end if
    end if

    call psb_bcast(ictxt,mtrx)
    call psb_bcast(ictxt,rhs)
    call psb_bcast(ictxt,filefmt)
    call psb_bcast(ictxt,choice)
    call psb_bcast(ictxt,itmax)
    call psb_bcast(ictxt,tol)

  end subroutine get_parms
end program mld_cexample_ml
