!!$ 
!!$ 
!!$                           MLD2P4  version 1.1
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 2.3.1)
!!$  
!!$  (C) Copyright 2008,2009
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
! This sample program solves a linear system by using BiCGStab preconditioned by
! RAS with overlap 2 and ILU(0) on the local blocks, as explained in Section 6.1 
! of the MLD2P4 User's and Reference Guide.
!
! The matrix and the rhs are read from files (if an rhs is not available, the
! unit rhs is set). 
!
program mld_sexample_ml
  use psb_sparse_mod
  use mld_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use data_input

  implicit none

! input parameters
  character(len=40) :: mtrx_file, rhs_file
  character(len=2)  :: filefmt

! sparse matrices
  type(psb_sspmat_type) :: A, aux_A

! descriptor of sparse matrices
  type(psb_desc_type):: desc_A

! preconditioner
  type(mld_sprec_type)  :: P

! right-hand side, solution and residual vectors
  real(psb_spk_), allocatable , save  :: b(:), x(:), r(:), &
       & x_glob(:), r_glob(:)
  real(psb_spk_), allocatable, target ::  aux_b(:,:)
  real(psb_spk_), pointer  :: b_glob(:)

! solver and preconditioner parameters
  real(psb_spk_)   :: tol, err
  integer          :: itmax, iter, istop
  integer          :: nlev

! parallel environment parameters
  integer            :: ictxt, iam, np

! other variables
  integer            :: i,info,j,m_problem
  integer(psb_long_int_k_) :: amatsize, precsize, descsize
  integer :: ierr, ircode
  real(psb_dpk_) :: t1, t2, tprec
  real(psb_spk_) :: resmx, resmxp
  character(len=20)  :: name
  integer, parameter :: iunit=12

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

!  get parameters

  call get_parms(ictxt,mtrx_file,rhs_file,filefmt,itmax,tol)

  call psb_barrier(ictxt)
  t1 = psb_wtime()  

! read and assemble the matrix A and the right-hand side b
! using PSBLAS routines for sparse matrix / vector management

  if (iam==psb_root_) then
    select case(psb_toupper(filefmt)) 
    case('MM') 
      ! For Matrix Market we have an input file for the matrix
      ! and an (optional) second file for the RHS. 
      call mm_mat_read(aux_a,info,iunit=iunit,filename=mtrx_file)
      if (info == 0) then 
        if (rhs_file /= 'NONE') then
          call mm_vet_read(aux_b,info,iunit=iunit,filename=rhs_file)
        end if
      end if
      
    case ('HB')
      ! For Harwell-Boeing we have a single file which may or may not
      ! contain an RHS.
      call hb_read(aux_a,info,iunit=iunit,b=aux_b,filename=mtrx_file)
      
    case default
      info = -1 
      write(0,*) 'Wrong choice for fileformat ', filefmt
    end select
    if (info /= 0) then
      write(0,*) 'Error while reading input matrix '
      call psb_abort(ictxt)
    end if
    
    m_problem = aux_a%m
    call psb_bcast(ictxt,m_problem)
    
    ! At this point aux_b may still be unallocated
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
  call psb_matdist(aux_A, A, ictxt, &
         & desc_A,b_glob,b,info, parts=part_block)

  t2 = psb_wtime() - t1

  call psb_amx(ictxt, t2)

  if (iam==psb_root_) then
    write(*,'(" ")')
    write(*,'("Time to read and partition matrix : ",es12.5)')t2
    write(*,'(" ")')
  end if

! set RAS with overlap 2 and ILU(0) on the local blocks

  call mld_precinit(P,'AS',info)
  call mld_precset(P,mld_sub_ovr_,2,info)

! build the preconditioner

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

  call psb_krylov('BICGSTAB',A,P,b,x,tol,desc_A,info,itmax,iter,err,istop=2)

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
    write(*,'("Matrix: ",A)')mtrx_file
    write(*,'("Computed solution on ",i8," processors")')np
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
      write(20,*) 'computed solution on ',np,' processor(s).'
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
  subroutine get_parms(ictxt,mtrx,rhs,filefmt,itmax,tol)

    use psb_sparse_mod
    implicit none

    integer             :: ictxt, itmax
    real(psb_spk_)      :: tol
    character(len=*)    :: mtrx, rhs,filefmt
    integer             :: iam, np

    call psb_info(ictxt,iam,np)

    if (iam==psb_root_) then
      ! read input parameters
      call read_data(mtrx,5)
      call read_data(rhs,5)
      call read_data(filefmt,5)
      call read_data(itmax,5)
      call read_data(tol,5)
    end if

    call psb_bcast(ictxt,mtrx)
    call psb_bcast(ictxt,rhs)
    call psb_bcast(ictxt,filefmt)
    call psb_bcast(ictxt,itmax)
    call psb_bcast(ictxt,tol)

  end subroutine get_parms
end program mld_sexample_ml
