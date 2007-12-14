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
! File: mld_ziluk_fct.f90.
!
! Subroutine: mld_ziluk_fct.
! Version:    complex.
! Contains:   mld_ziluk_fctint, iluk_copyin, iluk_fact, iluk_copyout.
!
!  This routine computes either the ILU(k) or the MILU(k) factorization of the
!  local part of the matrix stored into a. These factorizations are used to
!  build the 'base preconditioner' (block-Jacobi preconditioner/solver, Additive
!  Schwarz preconditioner) corresponding to a certain level of a multilevel
!  preconditioner.
!
!  Details on the above factorizations can be found in
!    Y. Saad, Iterative Methods for Sparse Linear Systems, Second Edition,
!    SIAM, 2003, Chapter 10.
!
!  The local matrix to be factorized is stored into a and blck, as specified in
!  the description of the arguments below. The storage format for both the L and
!  U factors is CSR. The diagonal of the U factor is stored separately (actually,
!  the inverse of the diagonal entries is stored; this is then managed in the solve
!  stage associated to the ILU(k)/MILU(k) factorization).
!  
!
! Arguments:
!    fill_in -  integer, input.
!               The fill-in level k in ILU(k)/MILU(k).
!    ialg    -  integer, input.
!               The type of incomplete factorization to be performed.
!               The MILU(k) factorization is computed if ialg = 2 (= mld_milu_n_);
!               the ILU(k) factorization otherwise.
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
subroutine mld_ziluk_fct(fill_in,ialg,a,l,u,d,info,blck)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_ziluk_fct
  implicit none

  ! Arguments
  integer, intent(in)                 :: fill_in, ialg
  integer, intent(out)                :: info
  type(psb_zspmat_type),intent(in)    :: a
  type(psb_zspmat_type),intent(inout) :: l,u
  type(psb_zspmat_type),intent(in), optional, target :: blck
  complex(kind(1.d0)), intent(inout)     ::  d(:)
  !     Local Variables
  integer   :: l1, l2, m, err_act
  
  type(psb_zspmat_type), pointer  :: blck_
  character(len=20)   :: name, ch_err
  logical, parameter :: debug=.false.

  name='mld_ziluk_fct'
  info = 0
  call psb_erractionsave(err_act)

  if (debug) write(0,*) 'mld_diluk_fct: start'

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
  ! Compute the ILU(k) or the MILU(k) factorization, depending on ialg
  !
  if (debug) write(0,*) 'mld_ziluk_fct: calling fctint'
  call mld_ziluk_fctint(fill_in,ialg,m,a%m,a,blck_%m,blck_,&
       & d,l%aspk,l%ia1,l%ia2,u%aspk,u%ia1,u%ia2,l1,l2,info)
  if (info /= 0) then
     info=4010
     ch_err='mld_ziluk_fctint'
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
  ! Subroutine: mld_ziluk_fctint.
  ! Version:    complex.
  ! Note: internal subroutine of mld_ziluk_fct.
  !
  !  This routine computes either the ILU(k) or the MILU(k) factorization of the
  !  local part of the matrix stored into a. These factorizations are used to build
  !  the 'base preconditioner' (block-Jacobi preconditioner/solver, Additive Schwarz
  !  preconditioner) corresponding to a certain level of a multilevel preconditioner.
  !
  !  The local matrix to be factorized is stored into a and b, as specified in the
  !  description of the arguments below. The storage format for both the L and U
  !  factors is CSR. The diagonal of the U factor is stored separately (actually,
  !  the inverse of the diagonal entries is stored; this is then managed in the
  !  solve stage associated to the ILU(k)/MILU(k) factorization).
  !
  !
  ! Arguments:
  !    fill_in -  integer, input.
  !               The fill-in level k in ILU(k)/MILU(k).
  !    ialg    -  integer, input.
  !               The type of incomplete factorization to be performed.
  !               The MILU(k) factorization is computed if ialg = 2 (= mld_milu_n_);
  !               the ILU(k) factorization otherwise.
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
  subroutine mld_ziluk_fctint(fill_in,ialg,m,ma,a,mb,b,&
       & d,laspk,lia1,lia2,uaspk,uia1,uia2,l1,l2,info)

    use psb_base_mod

    implicit none

  ! Arguments 
    integer, intent(in)            :: fill_in, ialg
    type(psb_zspmat_type)          :: a,b
    integer                        :: m,ma,mb,l1,l2,info
    integer, dimension(:), allocatable          :: lia1,lia2,uia1,uia2
    complex(kind(1.d0)), dimension(:), allocatable :: laspk,uaspk
    complex(kind(1.d0)), dimension(:)              :: d

  ! Local variables
    integer :: i, ktrw,err_act, nidx
    integer, allocatable            :: uplevs(:), rowlevs(:),idxs(:)
    complex(kind(1.d0)), allocatable :: row(:)
    type(psb_int_heap) :: heap
    logical,parameter  :: debug=.false.
    type(psb_zspmat_type) :: trw
    character(len=20), parameter  :: name='mld_ziluk_fctint'
    character(len=20)             :: ch_err

    if (psb_get_errstatus() /= 0) return 
    info=0
    call psb_erractionsave(err_act)

    m = ma+mb

    !
    ! Allocate a temporary buffer for the iluk_copyin function 
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
    ! Allocate memory to hold the entries of a row and the corresponding
    ! fill levels
    !
    allocate(uplevs(size(uaspk)),rowlevs(m),row(m),stat=info)
    if (info /= 0) then
      info=4010
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999
    end if

    uplevs(:)  = m+1
    row(:)     = zzero
    rowlevs(:) = -(m+1)
    
    !
    ! Cycle over the matrix rows
    !
    do i = 1, m
      
      if (debug.and.(mod(i,500)==1)) write(0,*)'LUINT: Loop index ',i,ma

      !
      ! At each iteration of the loop we keep in a heap the column indices
      ! affected by the factorization. The heap is initialized and filled
      ! in the iluk_copyin routine, and updated during the elimination, in
      ! the iluk_fact routine. The heap is ideal because at each step we need
      ! the lowest index, but we also need to insert new items, and the heap
      ! allows to do both in log time. 
      !
      d(i) = zzero
      if (i<=ma) then 
        !
        ! Copy into trw the i-th local row of the matrix, stored in a 
        ! 
        call iluk_copyin(i,ma,a,1,m,row,rowlevs,heap,ktrw,trw)
      else
        !
        ! Copy into trw the i-th local row of the matrix, stored in b
        ! (as (i-ma)-th row) 
        ! 
        call iluk_copyin(i-ma,mb,b,1,m,row,rowlevs,heap,ktrw,trw)
      endif

      if (debug) write(0,*)'LUINT: input Copy done'

      ! Do an elimination step on the current row. It turns out we only
      ! need to keep track of fill levels for the upper triangle, hence we
      ! do not have a lowlevs variable.
      !
      call iluk_fact(fill_in,i,m,row,rowlevs,heap,&
           & d,uia1,uia2,uaspk,uplevs,nidx,idxs)
      !
      ! Copy the row into laspk/d(i)/uaspk
      ! 
      call iluk_copyout(fill_in,ialg,i,m,row,rowlevs,nidx,idxs,&
           & l1,l2,lia1,lia2,laspk,d,uia1,uia2,uaspk,uplevs)

    end do

    !
    ! And we're done, so deallocate the memory
    !
    deallocate(uplevs,rowlevs,row,stat=info)
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
  end subroutine mld_ziluk_fctint

  !
  ! Subroutine: iluk_copyin.
  ! Version:    complex.
  ! Note: internal subroutine of mld_ziluk_fct.
  !
  !  This routine copies a row of a sparse matrix A, stored in the sparse matrix
  !  structure a, into the array row and stores into a heap the column indices of
  !  the nonzero entries of the copied row. The output array row is such that it
  !  contains a full row of A, i.e. it contains also the zero entries of the row.
  !  This is useful for the elimination step performed by iluk_fact after the call
  !  to iluk_copyin (see mld_iluk_fctint).
  !  The routine also sets to zero the entries of the array rowlevs corresponding
  !  to the nonzero entries of the copied row (see the description of the arguments
  !  below).
  !
  !  If the sparse matrix is in CSR format, a 'straight' copy is performed;
  !  otherwise psb_sp_getblk is used to extract a block of rows, which is then
  !  copied, row by row, into the array row, through successive calls to
  !  ilu_copyin.
  !
  !  This routine is used by mld_ziluk_fctint in the computation of the
  !  ILU(k)/MILU(k) factorization of a local sparse matrix.
  !  
  !
  ! Arguments:
  !    i       -  integer, input.
  !               The local index of the row to be extracted from the 
  !               sparse matrix structure a.
  !    m       -  integer, input.
  !               The number of rows of the local matrix stored into a.
  !    a       -  type(<psb_zspmat_type>), input.
  !               The sparse matrix structure containing the row to be copied.
  !    jmin    -  integer, input.
  !               The minimum valid column index.
  !    jmax    -  integer, input.
  !               The maximum valid column index.
  !               The output matrix will contain a clipped copy taken from
  !               a(1:m,jmin:jmax).
  !    row     -  complex(kind(1.d0)), dimension(:), input/output.
  !               In input it is the null vector (see mld_iluk_fctint and
  !               iluk_copyout). In output it contains the row extracted
  !               from the matrix A. It actually contains a full row, i.e.
  !               it contains also the zero entries of the row.
  !    rowlevs -  integer, dimension(:), input/output.
  !               In input rowlevs(k) = -(m+1) for k=1,...,m. In output
  !               rowlevs(k) = 0 for 1 <= k <= jmax and A(i,k) /=0, for
  !               future use in iluk_fact.
  !    heap    -  type(psb_int_heap), input/output.
  !               The heap containing the column indices of the nonzero
  !               entries in the array row.
  !               Note: this argument is intent(inout) and not only intent(out)
  !               to retain its allocation, done by psb_init_heap inside this
  !               routine.
  !    ktrw    -  integer, input/output.
  !               The index identifying the last entry taken from the
  !               staging buffer trw. See below.
  !    trw     -  type(psb_dspmat_type), input/output.
  !               A staging buffer. If the matrix A is not in CSR format, we use
  !               the psb_sp_getblk routine and store its output in trw; when we 
  !               need to call psb_sp_getblk we do it for a block of rows, and then
  !               we consume them from trw in successive calls to this routine,
  !               until we empty the buffer. Thus we will make a call to psb_sp_getblk
  !               every nrb calls to copyin. If A is in CSR format it is unused.
  !
  subroutine iluk_copyin(i,m,a,jmin,jmax,row,rowlevs,heap,ktrw,trw)

    use psb_base_mod

    implicit none
  
  ! Arguments 
    type(psb_zspmat_type) :: a,trw
    integer               :: i, rowlevs(:),m,ktrw,jmin,jmax
    complex(kind(1.d0))   :: row(:)
    type(psb_int_heap)    :: heap

  ! Local variables
    integer               :: k,j,info,irb
    integer, parameter    :: nrb=16
    character(len=20), parameter  :: name='mld_ziluk_fctint'
    character(len=20)             :: ch_err

    if (psb_get_errstatus() /= 0) return 
    info=0
    call psb_erractionsave(err_act)
    call psb_init_heap(heap,info) 

    if (toupper(a%fida)=='CSR') then

      !
      ! Take a fast shortcut if the matrix is stored in CSR format
      !
      
      do j = a%ia2(i), a%ia2(i+1) - 1
        k          = a%ia1(j)
        if ((jmin<=k).and.(k<=jmax)) then 
          row(k)     = a%aspk(j)
          rowlevs(k) = 0
          call psb_insert_heap(k,heap,info)
        end if
      end do

    else

      !
      ! Otherwise use psb_sp_getblk, slower but able (in principle) of 
      ! handling any format. In this case, a block of rows is extracted
      ! instead of a single row, for performance reasons, and these
      ! rows are copied one by one into the array row, through successive
      ! calls to iluk_copyin.
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
      
      do 
        if (ktrw > trw%infoa(psb_nnz_)) exit
        if (trw%ia1(ktrw) > i) exit
        k          = trw%ia2(ktrw)
        if ((jmin<=k).and.(k<=jmax)) then 
          row(k)     = trw%aspk(ktrw)
          rowlevs(k) = 0
          call psb_insert_heap(k,heap,info)
        end if
        ktrw       = ktrw + 1
      enddo
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

  end subroutine iluk_copyin

  !
  ! Subroutine: iluk_fact.
  ! Version:    complex.
  ! Note: internal subroutine of mld_ziluk_fct.
  !
  !  This routine does an elimination step of the ILU(k) factorization on a
  !  single matrix row (see the calling routine mld_iluk_fctint).
  !
  !  This step is also the base for a MILU(k) elimination step on the row (see
  !  iluk_copyout). This routine is used by mld_ziluk_fctint in the computation
  !  of the ILU(k)/MILU(k) factorization of a local sparse matrix.
  !
  !  NOTE: it turns out we only need to keep track of the fill levels for
  !  the upper triangle.
  !
  !
  ! Arguments
  !    fill_in -  integer, input.
  !               The fill-in level k in ILU(k).
  !    i       -  integer, input.
  !               The local index of the row to which the factorization is
  !               applied.
  !    m       -  integer, input.
  !               The number of rows of the local matrix to which the row
  !               belongs.
  !    row     -  complex(kind(1.d0)), dimension(:), input/output.
  !               In input it contains the row to which the elimination step
  !               has to be applied. In output it contains the row after the
  !               elimination step. It actually contains a full row, i.e.
  !               it contains also the zero entries of the row.
  !    rowlevs -  integer, dimension(:), input/output.
  !               In input rowlevs(k) = 0 if the k-th entry of the row is
  !               nonzero, and rowlevs(k) = -(m+1) otherwise. In output
  !               rowlevs(k) contains the fill kevel of the k-th entry of
  !               the row after the current elimination step; rowlevs(k) = -(m+1)
  !               means that the k-th row entry is zero throughout the elimination
  !               step.
  !    heap    -  type(psb_int_heap), input/output.
  !               The heap containing the column indices of the nonzero entries
  !               in the processed row. In input it contains the indices concerning
  !               the row before the elimination step, while in output it contains
  !               the indices concerning the transformed row.
  !    d       -  complex(kind(1.d0)), input.
  !               The inverse of the diagonal entries of the part of the U factor
  !               above the current row (see iluk_copyout).
  !    uia1    -  integer, dimension(:), input.
  !               The column indices of the nonzero entries of the part of the U
  !               factor above the current row, stored in uaspk row by row (see
  !               iluk_copyout, called by mld_ziluk_fctint), according to the CSR
  !               storage format.
  !    uia2    -  integer, dimension(:), input.
  !               The indices identifying the first nonzero entry of each row of
  !               the U factor above the current row, stored in uaspk row by row
  !               (see iluk_copyout, called by mld_ziluk_fctint), according to
  !               the CSR storage format.
  !    uaspk   -  complex(kind(1.d0)), dimension(:), input.
  !               The entries of the U factor above the current row (except the
  !               diagonal ones), stored according to the CSR format.
  !    uplevs  -  integer, dimension(:), input.
  !               The fill levels of the nonzero entries in the part of the
  !               U factor above the current row.
  !    nidx    -  integer, output.
  !               The number of entries of the array row that have been
  !               examined during the elimination step. This will be used
  !               by the routine iluk_copyout.
  !    idxs    -  integer, dimension(:), allocatable, input/output.
  !               The indices of the entries of the array row that have been
  !               examined during the elimination step.This will be used by
  !               by the routine iluk_copyout.
  !               Note: this argument is intent(inout) and not only intent(out)
  !               to retain its allocation, done by this routine.
  !
  subroutine iluk_fact(fill_in,i,m,row,rowlevs,heap,d,uia1,uia2,uaspk,uplevs,nidx,idxs)

    use psb_base_mod

    implicit none 

  ! Arguments
    type(psb_int_heap)    :: heap 
    integer               :: i,m, rowlevs(:),fill_in,nidx
    integer, allocatable  :: idxs(:)
    integer               :: uia1(:),uia2(:),uplevs(:)
    complex(kind(1.d0))   :: row(:), uaspk(:),d(:)

    ! Local variables
    integer               :: k,j,lrwk,jj,info, lastk
    complex(kind(1.d0))   :: rwk

    if (.not.allocated(idxs)) then 
      allocate(idxs(200),stat=info)
    endif
    nidx = 0
    lastk = -1

    !
    ! Do while there are indices to be processed
    !
    do

      call psb_heap_get_first(k,heap,info) 
      if (info < 0) exit

      ! 
      ! Just in case an index has been put on the heap more than once.
      !
      if (k == lastk) cycle

      lastk = k 
      nidx = nidx + 1
      if (nidx>size(idxs)) then 
        call psb_realloc(nidx+psb_heap_resize,idxs,info)
      end if
      idxs(nidx) = k
      if ((row(k) /= zzero).and.(rowlevs(k) <= fill_in).and.(k<i)) then 
        !
        ! Note: since U is scaled while copying it out (see iluk_copyout),
        ! we can use rwk in the update below
        ! 
        rwk    = row(k)
        row(k) = row(k) * d(k)    ! d(k) == 1/a(k,k)          
        lrwk   = rowlevs(k)
          
        do jj=uia2(k),uia2(k+1)-1
          j = uia1(jj)
          if (j<=k) then 
            write(0,*) 'Error in accessing upper mat???',j,k,jj
          endif
          !
          ! Insert the index into the heap for further processing.
          ! The fill levels are initialized to a negative value. If we find
          ! one, it means that it is an as yet untouched index, so we need
          ! to insert it; otherwise it is already on the heap, there is no
          ! need to insert it more than once. 
          !
          if (rowlevs(j)<0) then 
            call psb_insert_heap(j,heap,info)
            rowlevs(j) = abs(rowlevs(j))
          end if
          !
          ! Update row(j) and the corresponding fill level
          !
          row(j)     = row(j) - rwk * uaspk(jj)
          rowlevs(j) = min(rowlevs(j),lrwk+uplevs(jj)+1)
        end do

      end if
    end do

  end subroutine iluk_fact

  !
  ! Subroutine: iluk_copyout.
  ! Version:    complex.
  ! Note: internal subroutine of mld_ziluk_fct.
  !
  !  This routine copies a matrix row, computed by iluk_fact by applying an
  !  elimination step of the ILU(k) factorization, into the arrays laspk, uaspk,
  !  d, corresponding to the L factor, the U factor and the diagonal of U,
  !  respectively.
  !
  !  Note that
  !  - the part of the row stored into uaspk is scaled by the corresponding diagonal
  !    entry, according to the LDU form of the incomplete factorization;
  !  - the inverse of the diagonal entries of U is actually stored into d; this is
  !    then managed in the solve stage associated to the ILU(k)/MILU(k) factorization;
  !  - if the MILU(k) factorization has been required (ialg == mld_milu_n_), the
  !    row entries discarded because their fill levels are too high are added to
  !    the diagonal entry of the row;
  !  - the row entries are stored in laspk and uaspk according to the CSR format;
  !  - the arrays row and rowlevs are re-initialized for future use in mld_iluk_fct
  !    (see also iluk_copyin and iluk_fact).
  !
  !  This routine is used by mld_ziluk_fctint in the computation of the
  !  ILU(k)/MILU(k) factorization of a local sparse matrix.
  !  
  !
  ! Arguments:
  !    fill_in -  integer, input.
  !               The fill-in level k in ILU(k)/MILU(k).
  !    ialg    -  integer, input.
  !               The type of incomplete factorization considered. The MILU(k)
  !               factorization is computed if ialg = 2 (= mld_milu_n_); the
  !               ILU(k) factorization otherwise.
  !    i       -  integer, input.
  !               The local index of the row to be copied.
  !    m       -  integer, input.
  !               The number of rows of the local matrix under factorization.
  !    row     -  complex(kind(1.d0)), dimension(:), input/output.
  !               It contains, input, the row to be copied, and, in output,
  !               the null vector (the latter is used in the next call to
  !               iluk_copyin in mld_iluk_fact).
  !    rowlevs -  integer, dimension(:), input/output.
  !               In input rowlevs(k) contains the fill kevel of the k-th entry
  !               of the row to be copied. rowlevs(k) = -(m+1) indicates that
  !               this entry is zero; however, any rowlevs(k) = -(m+1) is not
  !               used by the routine. In output rowlevs(k) = -(m+1) for all k's
  !               (this is an inizialization for the next call to iluk_copyin
  !               in mld_iluk_factint).
  !    nidx    -  integer, input.
  !               The number of entries of the array row that have been examined
  !               during the elimination step carried out by the routine iluk_fact.
  !    idxs    -  integer, dimension(:), allocatable, input.
  !               The indices of the entries of the array row that have been
  !               examined during the elimination step carried out by the routine
  !               iluk_fact.
  !    l1      -  integer, input/output.
  !               Pointer to the last occupied entry of laspk.
  !    l2      -  integer, input/output.
  !               Pointer to the last occupied entry of uaspk.
  !    lia1    -  integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the L factor,
  !               copied in laspk row by row (see mld_ziluk_fctint), according
  !               to the CSR storage format.
  !    lia2    -  integer, dimension(:), input/output.
  !               The indices identifying the first nonzero entry of each row
  !               of the L factor, copied in laspk row by row (see 
  !               mld_ziluk_fctint), according to the CSR storage format.
  !    laspk   -  complex(kind(1.d0)), dimension(:), input/output.
  !               The array where the entries of the row corresponding to the
  !               L factor are copied.
  !    d       -  complex(kind(1.d0)), dimension(:), input/output.
  !               The array where the inverse of the diagonal entry of the
  !               row is copied (only d(i) is used by the routine). 
  !    uia1    -  integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the U factor
  !               copied in uaspk row by row (see mld_ziluk_fctint), according
  !               to the CSR storage format.
  !    uia2    -  integer, dimension(:), input/output.
  !               The indices identifying the first nonzero entry of each row
  !               of the U factor copied in uaspk row by row (see
  !               mld_zilu_fctint), according to the CSR storage format.
  !    uaspk   -  complex(kind(1.d0)), dimension(:), input/output.
  !               The array where the entries of the row corresponding to the
  !               U factor are copied.
  !    uplevs  -  integer, dimension(:), input.
  !               The fill levels of the nonzero entries in the part of the
  !               U factor above the current row.
  !
  subroutine iluk_copyout(fill_in,ialg,i,m,row,rowlevs,nidx,idxs,&
       &  l1,l2,lia1,lia2,laspk,d,uia1,uia2,uaspk,uplevs)

    use psb_base_mod

    implicit none 

    ! Arguments
    integer               :: fill_in,ialg,i, rowlevs(:),l1,l2,m,nidx,idxs(:)
    integer, allocatable  :: uia1(:),uia2(:), lia1(:),lia2(:),uplevs(:)
    complex(kind(1.d0)),allocatable    :: uaspk(:), laspk(:)
    complex(kind(1.d0))   :: row(:), d(:)

    ! Local variables
    integer               :: j,isz,info,err_act,int_err(5),idxp
    character(len=20), parameter  :: name='mld_ziluk_fctint'
    character(len=20)             :: ch_err

    if (psb_get_errstatus() /= 0) return 
    info=0
    call psb_erractionsave(err_act)

    d(i) = dzero

    do idxp=1,nidx

      j = idxs(idxp)

      if (j<i) then
        !
        ! Copy the lower part of the row
        !  
        if (rowlevs(j) <= fill_in) then 
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
          lia1(l1)   = j
          laspk(l1)  = row(j)
        else if (ialg == mld_milu_n_) then
          !
          ! MILU(k): add discarded entries to the diagonal one
          !
          d(i) = d(i) + row(j)
        end if
        !
        ! Re-initialize row(j) and rowlevs(j)
        !
        row(j)     = zzero
        rowlevs(j) = -(m+1)

      else if (j==i) then
        !
        ! Copy the diagonal entry of the row and re-initialize
        ! row(j) and rowlevs(j)
        !  
        d(i)       = d(i) + row(i) 
        row(i)     = zzero
        rowlevs(i) = -(m+1)

      else if (j>i) then 
        !
        ! Copy the upper part of the row
        ! 
        if (rowlevs(j) <= fill_in) then 
          l2     = l2 + 1 
          if (size(uaspk) < l2) then 
            ! 
            ! Figure out a good reallocation size!
            !
            isz  = max((l2/i)*m,int(1.2*l2),l2+100)
            call psb_realloc(isz,uaspk,info) 
            if (info == 0) call psb_realloc(isz,uia1,info) 
            if (info == 0) call psb_realloc(isz,uplevs,info,pad=(m+1))
            if (info /= 0) then 
              info=4010
              call psb_errpush(info,name,a_err='Allocate')
              goto 9999
            end if
          end if
          uia1(l2)   = j
          uaspk(l2)  = row(j)
          uplevs(l2) = rowlevs(j) 
        else if (ialg == mld_milu_n_) then
          !
          ! MILU(k): add discarded entries to the diagonal one
          ! 
          d(i) = d(i) + row(j)
        end if
        !
        ! Re-initialize row(j) and rowlevs(j)
        !
        row(j)     = zzero
        rowlevs(j) = -(m+1)
      end if
    end do

    !
    ! Store the pointers to the first non occupied entry of in
    ! laspk and uaspk
    !
    lia2(i+1) = l1 + 1
    uia2(i+1) = l2 + 1

    !     
    ! Check the pivot size
    !
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
      d(i) = zone/d(i)
    end if

    !
    ! Scale the upper part
    !
    do j=uia2(i), uia2(i+1)-1
      uaspk(j) = d(i)*uaspk(j)
    end do

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if

  end subroutine iluk_copyout


end subroutine mld_ziluk_fct
