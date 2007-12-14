!!$
!!$ 
!!$                                MLD2P4
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS v.2.0)
!!$  
!!$  (C) Copyright 2007  Alfredo Buttari      University of Rome Tor Vergata
!!$                      Pasqua D'Ambra       ICAR-CNR, Naples
!!$                      Daniela di Serafino  Second University of Naples
!!$                      Salvatore Filippone  University of Rome Tor Vergata  
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
! File: mld_zilut_fct.f90.
!
! Subroutine: mld_zilut_fct.
! Version:    real.
! Contains:   mld_zilut_fctint, ilut_copyin, ilut_fact, ilut_copyout.
!
!  This routine computes the ILU(k,t) factorization of the local part of the
!  matrix stored into a. These factorization is used to build the 'base
!  preconditioner' (block-Jacobi preconditioner/solver, Additive Schwarz
!  preconditioner) corresponding to a certain level of a multilevel preconditioner.
!
!  Details on the above factorization can be found in
!    Y. Saad, Iterative Methods for Sparse Linear Systems, Second Edition,
!    SIAM, 2003, Chapter 10.
!
!  The local matrix to be factorized is stored into a and blck, as specified in
!  the description of the arguments below. The storage format for both the L and
!  U factors is CSR. The diagonal of the U factor is stored separately (actually,
!  the inverse of the diagonal entries is stored; this is then managed in the solve
!  stage associated to the ILU(k,t) factorization).
!  
!
! Arguments:
!    fill_in -  integer, input.
!               The fill-in parameter k in ILU(k,t).
!    thres   -  integer, input.
!               The threshold t, i.e. the drop tolerance, in ILU(k,t).
!    a       -  type(<psb_zspmat_type>), input.
!               The sparse matrix structure containing the local matrix to be
!               factorized. Note that if the 'base' Additive Schwarz preconditioner
!               has overlap greater than 0 and the matrix has not been reordered
!               (see mld_bjac_bld), then a contains only the 'original' local part
!               of the matrix to be factorized, i.e. the rows of the matrix held
!               by the calling process according to the initial data distribution.
!    l       -  type(<psb_zspmat_type>), input/output.
!               The L factor in the incomplete factorization.
!               Note: its allocation is managed by the calling routine mld_ilu_bld,
!               hence it cannot be only intent(out).
!    u       -  type(<psb_zspmat_type>), input/output.
!               The U factor (except its diagonal) in the incomplete factorization.
!               Note: its allocation is managed by the calling routine mld_ilu_bld,
!               hence it cannot be only intent(out).
!    d       -  complex(kind(1.d0)), dimension(:), input/output.
!               The inverse of the diagonal entries of the U factor in the incomplete
!               factorization.
!               Note: its allocation is managed by the calling routine mld_ilu_bld,
!               hence it cannot be only intent(out).
!    info    -  integer, output.                    
!               Error code.
!    blck    -  type(<psb_zspmat_type>), input, optional, target.
!               The sparse matrix structure containing the remote rows of the
!               matrix to be factorized, that have been retrieved by mld_asmat_bld
!               to build an Additive Schwarz base preconditioner with overlap
!               greater than 0. If the overlap is 0 or the matrix has been reordered
!               (see mld_bjac_bld), then blck does not contain any row.
!  
subroutine mld_zilut_fct(fill_in,thres,ialg,a,l,u,d,info,blck)
  
  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_zilut_fct
  implicit none

  ! Arguments
  integer, intent(in)                 :: fill_in, ialg
  real(kind(1.d0)), intent(in)        :: thres
  integer, intent(out)                :: info
  type(psb_zspmat_type),intent(in)    :: a
  type(psb_zspmat_type),intent(inout) :: l,u
  complex(kind(1.d0)), intent(inout)     :: d(:)
  type(psb_zspmat_type),intent(in), optional, target :: blck

  !     Local Variables
  integer   ::  l1, l2, m, err_act
  
  type(psb_zspmat_type), pointer  :: blck_
  character(len=20)   :: name, ch_err
  logical, parameter :: debug=.false.

  name='mld_zilut_fct'
  info = 0
  call psb_erractionsave(err_act)

  if (debug) write(0,*) 'mld_zilut_fct: start'

  ! 
  ! Point to / allocate memory for the incomplete factorization
  !
  if (present(blck)) then 
    blck_ => blck
  else
    allocate(blck_,stat=info) 
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    call psb_sp_all(0,0,blck_,1,info)
    if (info /= 0) then
       info=4010
       ch_err='psb_sp_all'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
    end if

  endif

  !
  ! Compute the ILU(k,t) factorization
  !
  if (debug) write(0,*) 'mld_zilut_fct: calling fctint'
  call mld_zilut_fctint(fill_in,thres,ialg,m,a%m,a,blck_%m,blck_,&
       & d,l%aspk,l%ia1,l%ia2,u%aspk,u%ia1,u%ia2,l1,l2,info)
  if (info /= 0) then
     info=4010
     ch_err='mld_zilut_fctint'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  !
  ! Store information on the L and U sparse matrices
  !
  l%infoa(1) = l1
  l%fida     = 'CSR'
  l%descra   = 'TLU'
  u%infoa(1) = l2
  u%fida     = 'CSR'
  u%descra   = 'TUU'
  l%m = m
  l%k = m
  u%m = m
  u%k = m

  !
  ! Nullify the pointer / deallocate the memory
  !
  if (present(blck)) then 
    blck_ => null() 
  else
    call psb_sp_free(blck_,info)
    if (info /= 0) then
       info=4010
       ch_err='psb_sp_free'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
    end if
    deallocate(blck_) 
  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
     call psb_error()
     return
  end if
  return

contains

  !
  ! Subroutine: mld_zilut_fctint.
  ! Version:    real.
  ! Note: internal subroutine of mld_zilut_fct.
  !
  !  This routine computes the ILU(k,t) factorization of the local part of the
  !  matrix stored into a. These factorization is used to build the 'base
  !  preconditioner' (block-Jacobi preconditioner/solver, Additive Schwarz
  !  preconditioner) corresponding to a certain level of a multilevel preconditioner.
  !
  !  The local matrix to be factorized is stored into a and b, as specified in the
  !  description of the arguments below. The storage format for both the L and U
  !  factors is CSR. The diagonal of the U factor is stored separately (actually,
  !  the inverse of the diagonal entries is stored; this is then managed in the
  !  solve stage associated to the ILU(k,t) factorization).
  !
  !
  ! Arguments:
  !    fill_in -  integer, input.
  !               The fill-in parameter k in ILU(k,t).
  !    thres   -  integer, input.
  !               The threshold t, i.e. the drop tolerance, in ILU(k,t).
  !    m       -  integer, output.
  !               The total number of rows of the local matrix to be factorized,
  !               i.e. ma+mb.
  !    ma      -  integer, input.
  !               The number of rows of the local submatrix stored into a.
  !    a       -  type(<psb_zspmat_type>), input.
  !               The sparse matrix structure containing the local matrix to be
  !               factorized. Note that, if the 'base' Additive Schwarz preconditioner
  !               has overlap greater than 0 and the matrix has not been reordered
  !               (see mld_bjac_bld), then a contains only the 'original' local part
  !               of the matrix to be factorized, i.e. the rows of the matrix held
  !               by the calling process according to the initial data distribution.
  !    mb      -  integer, input.
  !               The number of rows of the local submatrix stored into b.
  !    b       -  type(<psb_zspmat_type>), input.
  !               The sparse matrix structure containing the remote rows of the
  !               matrix to be factorized, that have been retrieved by mld_asmat_bld
  !               to build an Additive Schwarz base preconditioner with overlap
  !               greater than 0. If the overlap is 0 or the matrix has been reordered
  !               (see mld_bjac_bld), then b does not contain   any row.
  !    d       -  complex(kind(1.d0)), dimension(:), output.
  !               The inverse of the diagonal entries of the U factor in the incomplete
  !               factorization.
  !    laspk   -  complex(kind(1.d0)), dimension(:), input/output.
  !               The L factor in the incomplete factorization.
  !    lia1    -  integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the L factor,
  !               according to the CSR storage format.
  !    lia2    -  integer, dimension(:), input/output.
  !               The indices identifying the first nonzero entry of each row
  !               of the L factor in laspk, according to the CSR storage format. 
  !    uaspk   -  complex(kind(1.d0)), dimension(:), input/output.
  !               The U factor in the incomplete factorization.
  !               The entries of U are stored according to the CSR format.
  !    uia1    -  integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the U factor,
  !               according to the CSR storage format.
  !    uia2    -  integer, dimension(:), input/output.
  !               The indices identifying the first nonzero entry of each row
  !               of the U factor in uaspk, according to the CSR storage format. 
  !    l1      -  integer, output
  !               The number of nonzero entries in laspk.
  !    l2      -  integer, output
  !               The number of nonzero entries in uaspk.
  !    info    -  integer, output.           
  !               Error code.
  !
  subroutine mld_zilut_fctint(fill_in,thres,ialg,m,ma,a,mb,b,&
       & d,laspk,lia1,lia2,uaspk,uia1,uia2,l1,l2,info)

    use psb_base_mod

    implicit none 

  ! Arguments
    integer, intent(in)            :: fill_in, ialg
    real(kind(1.d0)), intent(in)   :: thres
    type(psb_zspmat_type)          :: a,b
    integer                        :: m,ma,mb,l1,l2,info
    integer, dimension(:), allocatable          :: lia1,lia2,uia1,uia2
    complex(kind(1.d0)), dimension(:), allocatable :: laspk,uaspk
    complex(kind(1.d0)), dimension(:)              :: d

    ! Local Variables
    integer :: i, ktrw,err_act,nidx,nlw,nup,jmaxup
    real(kind(1.d0)) :: nrmi
    integer, allocatable          :: idxs(:)
    complex(kind(1.d0)), allocatable :: row(:)
    type(psb_int_heap) :: heap
    logical,parameter  :: debug=.false.
    type(psb_zspmat_type) :: trw
    character(len=20), parameter  :: name='mld_zilut_fctint'
    character(len=20)             :: ch_err

    if (psb_get_errstatus() /= 0) return 
    info=0
    call psb_erractionsave(err_act)

    m = ma+mb

    !
    ! Allocate a temporary buffer for the ilut_copyin function 
    !
    if (debug) write(0,*)'LUINT Allocating TRW'
    call psb_sp_all(0,0,trw,1,info)
    if (info==0) call psb_ensure_size(m+1,lia2,info)
    if (info==0) call psb_ensure_size(m+1,uia2,info)

    if (info /= 0) then
      info=4010
      call psb_errpush(info,name,a_err='psb_sp_all')
      goto 9999
    end if
    if (debug) write(0,*)'LUINT Done  Allocating TRW'
    
    l1=0
    l2=0
    lia2(1) = 1
    uia2(1) = 1

    if (debug) write(0,*)'In DCSRLU Begin cycle',m,ma,mb

    !
    ! Allocate memory to hold the entries of a row
    !
    allocate(row(m),stat=info)
    if (info /= 0) then
      info=4010
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999
    end if

    row(:) = zzero

    !
    ! Cycle over the matrix rows
    !
    do i = 1, m

      if (debug) write(0,*)'LUINT: Loop index ',i
      !
      ! At each iteration of the loop we keep in a heap the column indices
      ! affected by the factorization. The heap is initialized and filled
      ! in the ilut_copyin function, and updated during the elimination, in
      ! the ilut_fact routine. The heap is ideal because at each step we need
      ! the lowest index, but we also need to insert new items, and the heap
      ! allows to do both in log time. 
      !
      d(i) = zzero
      if (i<=ma) then 
        call ilut_copyin(i,ma,a,i,1,m,nlw,nup,jmaxup,nrmi,row,heap,ktrw,trw)
      else
        call ilut_copyin(i-ma,mb,b,i,1,m,nlw,nup,jmaxup,nrmi,row,heap,ktrw,trw)
      endif

      if (debug) write(0,*)'LUINT: input Copy done'
      !
      ! Do an elimination step on current row
      !
      call ilut_fact(fill_in,thres,i,m,nrmi,row,heap,&
           & d,uia1,uia2,uaspk,nidx,idxs)
      !
      ! Copy the row into laspk/d(i)/uaspk
      ! 
      call ilut_copyout(fill_in,thres,i,m,nlw,nup,jmaxup,nrmi,row,nidx,idxs,&
           & l1,l2,lia1,lia2,laspk,d,uia1,uia2,uaspk)

    end do

    !
    ! And we're done, so deallocate the memory
    !
    deallocate(row,idxs,stat=info)
    if (info /= 0) then
      info=4010
      call psb_errpush(info,name,a_err='Deallocate')
      goto 9999
    end if
    if (info == 0) call psb_sp_free(trw,info)
    if (info /= 0) then
      info=4010
      ch_err='psb_sp_free'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    if (debug) write(0,*)'Leaving ilu_fct'

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine mld_zilut_fctint

  !
  ! Subroutine: ilut_copyin.
  ! Version:    complex.
  ! Note: internal subroutine of mld_zilut_fct.
  !
  !  This routine performs the following tasks:
  !  - copying a row of a sparse matrix A, stored in the sparse matrix structure a,
  !    into the array row;
  !  - storing into a heap the column indices of the nonzero entries of the copied
  !    row;
  !  - computing the column index of the first entry with maximum absolute value
  !    in the part of the row belonging to the upper triangle;            
  !  - computing the 2-norm of the row.
  !  The output array row is such that it contains a full row of A, i.e. it contains
  !  also the zero entries of the row. This is useful for the elimination step
  !  performed by ilut_fact after the call to ilut_copyin (see mld_ilut_fctint).
  !
  !  If the sparse matrix is in CSR format, a 'straight' copy is performed;
  !  otherwise psb_sp_getblk is used to extract a block of rows, which is then
  !  copied, row by row, into the array row, through successive calls to
  !  ilut_copyin.
  !
  !  This routine is used by mld_zilut_fctint in the computation of the ILU(k,t)
  !  factorization of a local sparse matrix.
  !  
  !
  ! Arguments:
  !    i       -  integer, input.
  !               The local index of the row to be extracted from the 
  !               sparse matrix structure a.
  !    m       -  integer, input.
  !               The number of rows of the local matrix stored into a.
  !    a       -  type(<psb_zspmat_type>), input.
  !               The sparse matrix structure containing the row to be
  !               copied.
  !    jd      -  integer, input.
  !               The column index of the diagonal entry of the row to be
  !               copied.
  !    jmin    -  integer, input.
  !               The minimum valid column index.
  !    jmax    -  integer, input.
  !               The maximum valid column index.
  !               The output matrix will contain a clipped copy taken from
  !               a(1:m,jmin:jmax).
  !    nlw     -  integer, output.
  !               The number of nonzero entries in the part of the row
  !               belonging to the lower triangle of the matrix.
  !    nup     -  integer, output.
  !               The number of nonzero entries in the part of the row
  !               belonging to the upper triangle of the matrix.
  !    jmaxup  -  integer, output.
  !               The column index of the first entry with maximum absolute
  !               value in the part of the row belonging to the upper triangle
  !    nrmi    -  real(kind(1.d0)), output.
  !               The 2-norm of the current row.
  !    row     -  complex(kind(1.d0)), dimension(:), input/output.
  !               In input it is the null vector (see mld_ilut_fctint and
  !               ilut_copyout). In output it contains the row extracted
  !               from the matrix A. It actually contains a full row, i.e.
  !               it contains also the zero entries of the row.
  !    rowlevs -  integer, dimension(:), input/output.
  !               In input rowlevs(k) = -(m+1) for k=1,...,m. In output
  !               rowlevs(k) = 0 for 1 <= k <= jmax and A(i,k) /=0, for
  !               future use in ilut_fact.
  !    heap    -  type(psb_int_heap), input/output.
  !               The heap containing the column indices of the nonzero
  !               entries in the array row.
  !               Note: this argument is intent(inout) and not only intent(out)
  !               to retain its allocation, done by psb_init_heap inside this
  !               routine.
  !    ktrw    -  integer, input/output.
  !               The index identifying the last entry taken from the
  !               staging buffer trw. See below.
  !    trw     -  type(psb_zspmat_type), input/output.
  !               A staging buffer. If the matrix A is not in CSR format, we use
  !               the psb_sp_getblk routine and store its output in trw; when we 
  !               need to call psb_sp_getblk we do it for a block of rows, and then
  !               we consume them from trw in successive calls to this routine,
  !               until we empty the buffer. Thus we will make a call to psb_sp_getblk
  !               every nrb calls to copyin. If A is in CSR format it is unused.
  !
  subroutine ilut_copyin(i,m,a,jd,jmin,jmax,nlw,nup,jmaxup,nrmi,row,heap,ktrw,trw)
    use psb_base_mod
    implicit none 
    type(psb_zspmat_type) :: a,trw
    integer               :: i, m,ktrw,jmin,jmax,jd,nlw,nup,jmaxup
    real(kind(1.d0))      :: nrmi
    complex(kind(1.d0))   :: row(:)
    type(psb_int_heap)    :: heap
    
    integer               :: k,j,info,irb,kin,nz
    integer, parameter    :: nrb=16
    real(kind(1.d0))      :: dmaxup
    real(kind(1.d0)), external    :: dznrm2
    character(len=20), parameter  :: name='mld_zilut_fctint'
    character(len=20)             :: ch_err

    if (psb_get_errstatus() /= 0) return 
    info=0
    call psb_erractionsave(err_act)

    call psb_init_heap(heap,info)

    !
    ! nrmi is the norm of the current sparse row (for the time being,
    ! we use the 2-norm).
    ! NOTE: the 2-norm below includes also elements that are outside
    ! [jmin:jmax] strictly. Is this really important? TO BE CHECKED.
    !

    nlw    = 0
    nup    = 0
    jmaxup = 0
    dmaxup = dzero
    nrmi   = dzero

    if (toupper(a%fida)=='CSR') then

      !
      ! Take a fast shortcut if the matrix is stored in CSR format
      ! 

      do j = a%ia2(i), a%ia2(i+1) - 1
        k          = a%ia1(j)
        if ((jmin<=k).and.(k<=jmax)) then 
          row(k)     = a%aspk(j)
          call psb_insert_heap(k,heap,info)
        end if
        if (k<jd) nlw = nlw + 1 
        if (k>jd) then 
          nup = nup + 1
          if (abs(row(k))>dmaxup) then 
            jmaxup = k
            dmaxup = abs(row(k))
          end if
        end if
      end do
      nz   = a%ia2(i+1) - a%ia2(i)
      nrmi = dznrm2(nz,a%aspk(a%ia2(i)),ione)
    else

      !
      ! Otherwise use psb_sp_getblk, slower but able (in principle) of 
      ! handling any format. In this case, a block of rows is extracted
      ! instead of a single row, for performance reasons, and these
      ! rows are copied one by one into the array row, through successive
      ! calls to ilut_copyin.
      !

      if ((mod(i,nrb) == 1).or.(nrb==1)) then 
        irb = min(m-i+1,nrb)
        call psb_sp_getblk(i,a,trw,info,lrw=i+irb-1)
        if (info /= 0) then
          info=4010
          ch_err='psb_sp_getblk'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        end if
        ktrw=1
      end if
      
      kin = ktrw
      do 
        if (ktrw > trw%infoa(psb_nnz_)) exit
        if (trw%ia1(ktrw) > i) exit
        k          = trw%ia2(ktrw)
        if ((jmin<=k).and.(k<=jmax)) then 
          row(k)     = trw%aspk(ktrw)
          call psb_insert_heap(k,heap,info)
        end if
        if (k<jd) nlw = nlw + 1 
        if (k>jd) then 
          nup = nup + 1
          if (abs(row(k))>dmaxup) then 
            jmaxup = k
            dmaxup = abs(row(k))
          end if
        end if
        ktrw       = ktrw + 1
      enddo
      nz = ktrw - kin
      nrmi = dznrm2(nz,trw%aspk(kin),ione)
    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine ilut_copyin

  !
  ! Subroutine: ilut_fact.
  ! Version:    complex.
  ! Note: internal subroutine of mld_zilut_fct.
  !
  !  This routine does an elimination step of the ILU(k,t) factorization on a single
  !  matrix row (see the calling routine mld_ilut_fctint). Actually, only the dropping
  !  rule based on the threshold is applied here. The dropping rule based on the
  !  fill-in is applied by ilut_copyout.
  !
  !  The routine is used by mld_zilut_fctint in the computation of the ILU(k,t)
  !  factorization of a local sparse matrix.
  !
  !
  ! Arguments
  !    fill_in -  integer, input.
  !               The fill-in parameter k in ILU(k,t).
  !    thres   -  integer, input.
  !               The threshold t, i.e. the drop tolerance, in ILU(k,t).
  !    i       -  integer, input.
  !               The local index of the row to which the factorization is applied.
  !    m       -  integer, input.
  !               The number of rows of the local matrix to which the row belongs.
  !    nrmi    -  real(kind(1.d0)), input.
  !               The 2-norm of the row to which the elimination step has to be
  !               applied.
  !    row     -  complex(kind(1.d0)), dimension(:), input/output.
  !               In input it contains the row to which the elimination step
  !               has to be applied. In output it contains the row after the
  !               elimination step. It actually contains a full row, i.e.
  !               it contains also the zero entries of the row.
  !    heap    -  type(psb_int_heap), input/output.
  !               The heap containing the column indices of the nonzero entries
  !               in the processed row. In input it contains the indices concerning
  !               the row before the elimination step, while in output it contains
  !               the previous indices plus the ones corresponding to transformed
  !               entries in the 'upper part' that have not been dropped.
  !    d       -  complex(kind(1.d0)), input.
  !               The inverse of the diagonal entries of the part of the U factor
  !               above the current row (see ilut_copyout).
  !    uia1    -  integer, dimension(:), input.
  !               The column indices of the nonzero entries of the part of the U
  !               factor above the current row, stored in uaspk row by row (see
  !               ilut_copyout, called by mld_zilut_fctint), according to the CSR
  !               storage format.
  !    uia2    -  integer, dimension(:), input.
  !               The indices identifying the first nonzero entry of each row of
  !               the U factor above the current row, stored in uaspk row by row
  !               (see ilut_copyout, called by mld_zilut_fctint), according to
  !               the CSR storage format.
  !    uaspk   -  complex(kind(1.d0)), dimension(:), input.
  !               The entries of the U factor above the current row (except the
  !               diagonal ones), stored according to the CSR format.
  !    nidx    -  integer, output.
  !               The number of entries of the array row that have been
  !               examined during the elimination step. This will be used
  !               by the routine ilut_copyout.
  !    idxs    -  integer, dimension(:), allocatable, input/output.
  !               The indices of the entries of the array row that have been
  !               examined during the elimination step.This will be used by
  !               by the routine ilut_copyout.
  !               Note: this argument is intent(inout) and not only intent(out)
  !               to retain its allocation, done by this routine.
  !
  subroutine ilut_fact(fill_in,thres,i,m,nrmi,row,heap,&
       & d,uia1,uia2,uaspk,nidx,idxs)

    use psb_base_mod

    implicit none 

  ! Arguments
    type(psb_int_heap)    :: heap 
    integer               :: i,m,fill_in,nidx
    real(kind(1.d0)), intent(in)   :: thres,nrmi
    integer, allocatable  :: idxs(:)
    integer               :: uia1(:),uia2(:)
    complex(kind(1.d0))      :: row(:), uaspk(:),d(:)

    ! Local Variables
    integer               :: k,j,jj,info, lastk
    complex(kind(1.d0))   :: rwk
    logical, parameter    :: debug=.false.

    call psb_ensure_size(200,idxs,info)

    nidx  = 0
    lastk = -1 
    !
    ! Do while there are indices to be processed
    !
    do

      call psb_heap_get_first(k,heap,info) 
      if (info < 0) exit

      ! 
      ! An index may have been put on the heap more than once.
      !
      if (k == lastk) cycle

      lastk = k 
      lowert: if (k<i)  then
 
        !
        ! Dropping rule based on the threshold: compare the absolute
        ! value of each updated entry of row with thres * 2-norm of row.
        !
        rwk    = row(k)
        row(k) = row(k) * d(k)
        if (abs(row(k)) < thres*nrmi) then
          ! 
          ! Drop the entry.
          !
          row(k) = zzero
          cycle
        else
          !
          ! Note: since U is scaled while copying it out (see ilut_copyout),
          ! we can use rwk in the update below.
          !           
          do jj=uia2(k),uia2(k+1)-1
            j = uia1(jj)
            if (j<=k) then 
              write(0,*) 'Error in accessing upper mat???',j,k,jj
            endif
            !
            ! Update row(j) and, if it is not to be discarded, insert
            ! its index into the heap for further processing.
            !
            row(j)     = row(j) - rwk * uaspk(jj)
            if (abs(row(j)) < thres*nrmi) then
              ! 
              ! Drop the entry.
              !
              row(j) = zzero
            else
              !
              ! Do the insertion.
              !
              call psb_insert_heap(j,heap,info)
            endif
          end do
        end if
      end if lowert

      !
      ! If we get here it is an index we need to keep on copyout.
      !
      nidx       = nidx + 1
      call psb_ensure_size(nidx,idxs,info,addsz=psb_heap_resize)      
      idxs(nidx) = k

    end do

    if (debug) then 
      write(0,*) 'At end of factint: ',i,nidx
      write(0,*) idxs(1:nidx)
      write(0,*) row(idxs(1:nidx))
    end if

  end subroutine ilut_fact

  !
  ! Subroutine: ilut_copyout.
  ! Version:    complex.
  ! Note: internal subroutine of mld_zilut_fct.
  !
  !  This routine copies a matrix row, computed by ilut_fact by applying an
  !  elimination step of the ILU(k,t) factorization, into the arrays laspk,
  !  uaspk, d, corresponding to the L factor, the U factor and the diagonal
  !  of U, respectively.
  !
  !  Note that
  !  - the dropping rule based on the fill-in is applied here and not in ilut_fact;
  !    it consists in keeping the nlw+k entries with largest absolute value in
  !    the 'lower part' of the row, and the nup+k ones in the 'upper part';
  !  - the entry in the upper part of the row which has maximum absolute value
  !    in the original matrix is included in the above nup+k entries anyway;
  !  - the part of the row stored into uaspk is scaled by the corresponding
  !    diagonal entry, according to the LDU form of the incomplete factorization;
  !  - the inverse of the diagonal entries of U is actually stored into d; this
  !    is then managed in the solve stage associated to the ILU(k,t) factorization;
  !  - the row entries are stored in laspk and uaspk according to the CSR format;
  !  - the array row is re-initialized for future use in mld_ilut_fct(see also
  !    ilut_copyin and ilut_fact).
  !
  !  This routine is used by mld_zilut_fctint in the computation of the      ILU(k,t)
  !  factorization of a local sparse matrix.
  !  
  !
  ! Arguments:
  !    fill_in -  integer, input.
  !               The fill-in parameter k in ILU(k,t).
  !    thres   -  integer, input.
  !               The threshold t, i.e. the drop tolerance, in ILU(k,t).
  !    i       -  integer, input.
  !               The local index of the row to be copied.
  !    m       -  integer, input.
  !               The number of rows of the local matrix under factorization.
  !    nlw     -  integer, input.
  !               The number of nonzero entries of the 'lower part' of the row
  !               in the initial matrix (i.e. the matrix before the factorization).
  !    nup     -  integer, input.
  !               The number of nonzero entries in the 'upper part' of the row
  !               in the initial matrix.
  !    jmaxup  -  integer, input.
  !               The column index of the first entry with maximum absolute
  !               value in the 'upper part' of the row in the initial matrix.
  !    nrmi    -  real(kind(1.d0)), input.
  !               The 2-norm of the current row in the initial matrix.
  !    row     -  complex(kind(1.d0)), dimension(:), input/output.
  !               It contains, input, the row to be copied, and, in output,
  !               the null vector (the latter is used in the next call to
  !               ilut_copyin in mld_ilut_fact).
  !    nidx    -  integer, input.
  !               The number of entries of the array row that have been examined
  !               during the elimination step carried out by the routine ilut_fact.
  !    idxs    -  integer, dimension(:), allocatable, input.
  !               The indices of the entries of the array row that have been
  !               examined during the elimination step carried out by the routine
  !               ilut_fact.
  !    l1      -  integer, input/output.
  !               Pointer to the last occupied entry of laspk.
  !    l2      -  integer, input/output.
  !               Pointer to the last occupied entry of uaspk.
  !    lia1    -  integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the L factor,
  !               copied in laspk row by row (see mld_zilut_fctint), according
  !               to the CSR storage format.
  !    lia2    -  integer, dimension(:), input/output.
  !               The indices identifying the first nonzero entry of each row
  !               of the L factor, copied in laspk row by row (see 
  !               mld_zilut_fctint), according to the CSR storage format.
  !    laspk   -  complex(kind(1.d0)), dimension(:), input/output.
  !               The array where the entries of the row corresponding to the
  !               L factor are copied.
  !    d       -  complex(kind(1.d0)), dimension(:), input/output.
  !               The array where the inverse of the diagonal entry of the
  !               row is copied (only d(i) is used by the routine). 
  !    uia1    -  integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the U factor
  !               copied in uaspk row by row (see mld_zilut_fctint), according
  !               to the CSR storage format.
  !    uia2    -  integer, dimension(:), input/output.
  !               The indices identifying the first nonzero entry of each row
  !               of the U factor copied in uaspk row by row (see
  !               mld_zilu_fctint), according to the CSR storage format.
  !    uaspk   -  complex(kind(1.d0)), dimension(:), input/output.
  !               The array where the entries of the row corresponding to the
  !               U factor are copied.
  !
  subroutine ilut_copyout(fill_in,thres,i,m,nlw,nup,jmaxup,nrmi,row, &
       & nidx,idxs,l1,l2,lia1,lia2,laspk,d,uia1,uia2,uaspk)

    use psb_base_mod

    implicit none 

    ! Arguments
    integer               :: fill_in,i,l1,l2,m,nidx,idxs(:)
    integer               :: nlw,nup,jmaxup
    integer, allocatable  :: uia1(:),uia2(:), lia1(:),lia2(:)
    real(kind(1.d0)), intent(in) :: thres,nrmi
    complex(kind(1.d0)),allocatable :: uaspk(:), laspk(:)
    complex(kind(1.d0))             :: row(:), d(:)

    ! Local variables
    complex(kind(1.d0)),allocatable :: xw(:)
    integer, allocatable         :: xwid(:), indx(:)
    complex(kind(1.d0))             :: witem
    integer                      :: widx
    integer                      :: k,isz,info,err_act,int_err(5),idxp, nz
    type(psb_dcomplex_idx_heap)  :: heap
    character(len=20), parameter :: name='mld_zilut_fctint'
    character(len=20)            :: ch_err
    logical                      :: fndmaxup
    logical, parameter           :: debug=.false.

    if (psb_get_errstatus() /= 0) return 
    info=0
    call psb_erractionsave(err_act)

    !
    ! Here we need to apply also the dropping rule base on the fill-in. 
    ! We do it by putting into a heap the elements that are not dropped
    ! by using the 2-norm rule, and then copying them out. 
    !
    ! The heap goes down on the entry absolute value, so the first item
    ! is the largest absolute value. 
    !

    call psb_init_heap(heap,info,dir=psb_asort_down_)

    if (info == 0) allocate(xwid(nidx),xw(nidx),indx(nidx),stat=info)
    if (info /= 0) then 
      info=4025
      call psb_errpush(info,name,i_err=(/3*nidx,0,0,0,0/),&
           & a_err='complex(kind(1.d0))')
      goto 9999      
    end if

    !
    ! First the lower part
    !

    nz   = 0
    idxp = 0

    do
 
      idxp = idxp + 1
      if (idxp > nidx) exit
      if (idxs(idxp) >= i) exit
      widx      = idxs(idxp)
      witem     = row(widx)
      if (debug) then 
        write(0,*) 'Lower: Deciding on drop of item ',witem,widx,thres,nrmi,thres*nrmi
      end if

      !
      ! Dropping rule based on the 2-norm
      !
      if (abs(witem) < thres*nrmi) cycle

      nz       = nz + 1 
      xw(nz)   = witem 
      xwid(nz) = widx
      call psb_insert_heap(witem,widx,heap,info)

    end do

    !
    ! Now we have to take out the first nlw+fill_in entries
    ! 
    if (nz <= nlw+fill_in) then
      ! 
      ! Just copy everything from xw, and it is already ordered
      !
    else
      nz = nlw+fill_in
      do k=1,nz
        call psb_heap_get_first(witem,widx,heap,info)
        xw(k)   = witem
        xwid(k) = widx
      end do
    end if

    !
    ! Now put things back into ascending column order
    !
    call psb_msort(xwid(1:nz),indx(1:nz),dir=psb_sort_up_)

    !
    ! Copy out the lower part of the row
    !
    do k=1,nz
      l1     = l1 + 1 
      if (size(laspk) < l1) then
        ! 
        ! Figure out a good reallocation size!
        ! 
        isz  = (max((l1/i)*m,int(1.2*l1),l1+100))
        call psb_realloc(isz,laspk,info) 
        if (info == 0) call psb_realloc(isz,lia1,info) 
        if (info /= 0) then 
          info=4010
          call psb_errpush(info,name,a_err='Allocate')
          goto 9999
        end if
      end if
      lia1(l1)   = xwid(k)
      laspk(l1)  = xw(indx(k))
    end do
    
    !
    ! Make sure idxp points to the diagonal entry
    !
    if (idxp <= size(idxs)) then 
      if (idxs(idxp) < i) then 
        do 
          idxp = idxp + 1
          if (idxp > nidx) exit
          if (idxs(idxp) >= i) exit
        end do
      end if
    end if
    if (idxp > size(idxs)) then 
      write(0,*) 'Warning: missing diagonal element in the row '
    else
      if (idxs(idxp) > i) then 
        write(0,*) 'Warning: missing diagonal element in the row '
      else if (idxs(idxp) /= i) then 
        write(0,*) 'Warning: impossible error: diagonal has vanished'
      else
        !
        ! Copy the diagonal entry
        !
        widx      = idxs(idxp)
        witem     = row(widx)
        d(i)      = witem
        if (abs(d(i)) < epstol) then
          !
          ! Too small pivot: unstable factorization
          !     
          info = 2
          int_err(1) = i
          write(ch_err,'(g20.10)') d(i)
          call psb_errpush(info,name,i_err=int_err,a_err=ch_err)
          goto 9999
        else
          !
          ! Compute 1/pivot
          !
          d(i) = done/d(i)
        end if
      end if
    end if

    !
    ! Now the upper part 
    !

    call psb_init_heap(heap,info,dir=psb_asort_down_)
    nz       = 0
    do
 
      idxp = idxp + 1
      if (idxp > nidx) exit
      widx      = idxs(idxp)
      if (widx <= i) then 
        write(0,*) 'Warning: lower triangle in upper copy',widx,i,idxp,idxs(idxp)
        cycle
      end if
      if (widx > m) then 
        write(0,*) 'Warning: impossible value',widx,i,idxp,idxs(idxp)
        cycle
      end if
      witem     = row(widx)
      if (debug) then 
        write(0,*) 'Upper: Deciding on drop of item ',witem,widx,&
             & jmaxup,thres,nrmi,thres*nrmi
      end if

      !
      ! Dropping rule based on the 2-norm. But keep the jmaxup-th entry anyway.
      !
      if ((widx /= jmaxup) .and. (abs(witem) < thres*nrmi)) then 
        cycle 
      end if

      nz       = nz + 1
      xw(nz)   = witem 
      xwid(nz) = widx
      call psb_insert_heap(witem,widx,heap,info)
      
    end do

    if (debug) then 
      write(0,*) 'Row ',i,' copyout: after first round  at upper:',nz,jmaxup
      write(0,*) xwid(1:nz)
      write(0,*) xw(1:nz)
      write(0,*) 'Dumping heap'
      call psb_dump_heap(0,heap,info)
    end if

    !
    ! Now we have to take out the first nup-fill_in entries. But make sure
    ! we include entry jmaxup.
    !
    if (nz <= nup+fill_in) then
      ! 
      ! Just copy everything from xw
      !
      fndmaxup=.true.
    else
      fndmaxup = .false.
      nz = nup+fill_in
      do k=1,nz
        call psb_heap_get_first(witem,widx,heap,info)
        xw(k)   = witem
        xwid(k) = widx
        if (widx == jmaxup) fndmaxup=.true.
      end do
    end if
    if ((i<jmaxup).and.(jmaxup<=m)) then 
      if (.not.fndmaxup) then 
        ! 
        ! Include entry jmaxup, if it is not already there.
        ! Put it in the place of the smallest coefficient. 
        !
        xw(nz)   = row(jmaxup) 
        xwid(nz) = jmaxup
      endif
    end if

    !
    ! Now we put things back into ascending column order
    !
    call psb_msort(xwid(1:nz),indx(1:nz),dir=psb_sort_up_)
    if (debug) then 
      write(0,*) 'Row ',i,' copyout: after sort at upper:',nz,jmaxup
      write(0,*) xwid(1:nz)
      write(0,*) xw(indx(1:nz))
    end if

    !
    ! Copy out the upper part of the row
    !
    do k=1,nz
      l2     = l2 + 1 
      if (size(uaspk) < l2) then
        ! 
        ! Figure out a good reallocation size!
        ! 
        isz  = max((l2/i)*m,int(1.2*l2),l2+100)
        call psb_realloc(isz,uaspk,info) 
        if (info == 0) call psb_realloc(isz,uia1,info) 
        if (info /= 0) then 
          info=4010
          call psb_errpush(info,name,a_err='Allocate')
          goto 9999
        end if
      end if
      uia1(l2)   = xwid(k)
      uaspk(l2)  = d(i)*xw(indx(k))
    end do

    !
    ! Set row to zero
    !
    do idxp=1,nidx
      row(idxs(idxp)) = zzero
    end do

    !
    ! Store the pointers to the first non occupied entry of in
    ! laspk and uaspk
    !
    lia2(i+1) = l1 + 1
    uia2(i+1) = l2 + 1

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if

  end subroutine ilut_copyout


end subroutine mld_zilut_fct
