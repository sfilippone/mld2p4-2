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
! File: mld_cexample_ml.f90
!
! This sample program solves a linear system by using BiCGStab coupled with
! one of the following multi-level preconditioner, as explained in Section 6.1
! of the MLD2P4 User's and Reference Guide:
! - choice = 1, default multi-level Schwarz preconditioner (Sec. 6.1, Fig. 2)
! - choice = 2, hybrid three-level Schwarz preconditioner (Sec. 6.1, Fig. 3)
! - choice = 3, additive three-level Schwarz preconditioner (Sec. 6.1, Fig. 4)
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

! sparse matrices
  type(psb_cspmat_type) :: A, aux_A

! descriptor of sparse matrices
  type(psb_desc_type):: desc_A

! preconditioner
  type(mld_cprec_type)  :: P

! right-hand side, solution and residual vectors
  complex(psb_spk_), allocatable , save  :: b(:), x(:), r(:), &
       & x_glob(:), r_glob(:)
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
  integer              :: i,info,j,m_problem,amatsize,descsize,precsize
  integer              :: ierr, ircode
  real(psb_dpk_) :: t1, t2, tprec
  real(psb_spk_) :: resmx, resmxp
  character(len=20)    :: name

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
  info=0
  call psb_set_errverbosity(2)

! get parameters

  call get_parms(ictxt,mtrx_file,rhs_file,choice,itmax,tol)

  call psb_barrier(ictxt)
  t1 = psb_wtime()  

! read and assemble the matrix A and the right-hand side b
! using PSBLAS routines for sparse matrix / vector management

  if (iam==psb_root_) then
    call read_mat(mtrx_file, aux_A, ictxt)

    m_problem = aux_A%m
    call psb_bcast(ictxt,m_problem)

    if(rhs_file /= 'NONE') then
      !  reading an rhs
      call read_rhs(rhs_file,aux_b,ictxt)
    end if

    if (psb_size(aux_b,1)==m_problem) then
      ! if any rhs were present, broadcast the first one
      write(0,'("Ok, got an rhs ")')
      b_glob =>aux_b(:,1)
    else
      write(*,'("Generating an rhs...")')
      write(*,'(" ")')
      call psb_realloc(m_problem,1,aux_b,ircode)
      if (ircode /= 0) then
        call psb_errpush(4000,name)
        goto 9999
      endif

      b_glob => aux_b(:,1)
      do i=1, m_problem
        b_glob(i) = 1.d0
      enddo
    endif
    call psb_bcast(ictxt,b_glob(1:m_problem))
  else
    call psb_bcast(ictxt,m_problem)
    call psb_realloc(m_problem,1,aux_b,ircode)
    if (ircode /= 0) then
      call psb_errpush(4000,name)
      goto 9999
    endif
    b_glob =>aux_b(:,1)
    call psb_bcast(ictxt,b_glob(1:m_problem)) 
  end if

  call psb_barrier(ictxt)
  if (iam==psb_root_) write(*,'("Partition type: block")')
    call psb_matdist(aux_A, A, part_block, ictxt, &
         & desc_A,b_glob,b,info)

  t2 = psb_wtime() - t1

  call psb_amx(ictxt, t2)

  if (iam==psb_root_) then
    write(*,'(" ")')
    write(*,'("Time to read and partition matrix : ",es10.4)')t2
    write(*,'(" ")')
  end if

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

! set the initial guess

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
  call psb_geaxpby(cone,b,czero,r,desc_A,info)
  call psb_spmm(-cone,A,x,cone,r,desc_A,info)
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
    write(*,'("Matrix: ",A)')mtrx_file
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

  allocate(x_glob(m_problem),r_glob(m_problem),stat=ierr)
  if (ierr /= 0) then 
    write(0,*) 'allocation error: no data collection'
  else
    call psb_gather(x_glob,x,desc_A,info,root=psb_root_)
    call psb_gather(r_glob,r,desc_A,info,root=psb_root_)
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
        write(20,998) i,x_glob(i),r_glob(i),b_glob(i)
      enddo
    end if
  end if
998 format(i8,4(2x,g20.14))
993 format(i6,4(1x,e12.6))

! deallocate the data structures

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
  subroutine get_parms(ictxt,mtrx,rhs,choice,itmax,tol)

    use psb_base_mod
    implicit none

    integer             :: ictxt, choice, itmax
    real(psb_spk_)      :: tol
    character(len=*)    :: mtrx, rhs
    integer             :: iam, np

    call psb_info(ictxt,iam,np)

    if (iam==psb_root_) then
      ! read input parameters
      call read_data(mtrx,5)
      call read_data(rhs,5)
      call read_data(choice,5)
      call read_data(itmax,5)
      call read_data(tol,5)
    end if

    call psb_bcast(ictxt,mtrx)
    call psb_bcast(ictxt,rhs)
    call psb_bcast(ictxt,choice)
    call psb_bcast(ictxt,itmax)
    call psb_bcast(ictxt,tol)

  end subroutine get_parms
end program mld_cexample_ml
