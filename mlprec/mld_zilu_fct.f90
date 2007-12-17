!!$
!!$ 
!!$                                MLD2P4
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS v.2.0)
!!$  
!!$  (C) Copyright 2007  Alfredo Buttari      University of Rome Tor Vergata
!!$			 Pasqua D'Ambra       ICAR-CNR, Naples
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
! File: mld_dilu_fct.f90
!
! Subroutine: mld_zilu_fct
! Version:    complex
! Contains:   mld_zilu_fctint, ilu_copyin
!
!  This routine computes either the ILU(0) or the MILU(0) factorization of the
!  local part of the matrix stored into a. These factorizations are used to
!  build the 'base preconditioner' (block-Jacobi preconditioner/solver, Additive
!  Schwarz preconditioner) corresponding to a certain level of a multilevel
!  preconditioner.
!
!  Details on the above factorizations can be found in
!    Y. Saad, Iterative Methods for Sparse Linear Systems, Second Edition,
!    SIAM, 2003, Chapter 10.
!
!  The local matrix to be factorized is stored into a and blck, as specified
!  in the description of the arguments below. The storage format for both the
!  L and U factors is CSR. The diagonal of the U factor is stored separately
!  (actually, the inverse of the diagonal entries is stored; this is then
!  managed in the solve stage associated to the ILU(0)/MILU(0) factorization).
!
!  The routine copies and factors "on the fly" from a and blck into l (L factor),
!  u (U factor, except its diagonal) and d (diagonal of U).
!
!  This implementation of ILU(0)/MILU(0) is faster than the implementation in
!  mld_diluk_fct (the latter routine performs the more general ILU(k)/MILU(k)).
!  
!
! Arguments:
!    ialg    -  integer, input.
!               The type of incomplete factorization to be performed.
!               The MILU(0) factorization is computed if ialg = 2 (= mld_milu_n_);
!               the ILU(0) factorization otherwise.
!    a       -  type(psb_zspmat_type), input.
!               The sparse matrix structure containing the local matrix to be
!               factorized. Note that if the 'base' Additive Schwarz preconditioner
!               has overlap greater than 0 and the matrix has not been reordered
!               (see mld_bjac_bld), then a contains     only the 'original' local part
!               of the matrix to be factorized, i.e. the rows of the matrix held
!               by the calling process according to the initial data distribution.
!    l       -  type(psb_zspmat_type), input/output.
!               The L factor in the incomplete factorization.
!               Note: its allocation is managed by the calling routine mld_ilu_bld,
!               hence it cannot be only intent(out).
!    u       -  type(psb_zspmat_type), input/output.
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
!    blck    -  type(psb_zspmat_type), input, optional, target.
!               The sparse matrix structure containing the remote rows of the
!               matrix to be factorized, that have been retrieved by mld_asmat_bld
!               to build an Additive Schwarz base preconditioner with overlap
!               greater than 0. If the overlap is 0 or the matrix has been reordered
!               (see mld_bjac_bld), then blck is empty.
!  
subroutine mld_zilu_fct(ialg,a,l,u,d,info,blck)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_zilu_fct
  implicit none

  ! Arguments
  integer, intent(in)                 :: ialg
  type(psb_zspmat_type),intent(in)    :: a
  type(psb_zspmat_type),intent(inout) :: l,u
  complex(kind(1.d0)), intent(inout)     :: d(:)
  integer, intent(out)                :: info
  type(psb_zspmat_type),intent(in), optional, target :: blck

  ! Local variables
  integer   :: l1, l2,m,err_act
  type(psb_zspmat_type), pointer  :: blck_
  character(len=20)   :: name, ch_err
  logical, parameter :: debug=.false.

  name='mld_zilu_fct'
  info = 0
  call psb_erractionsave(err_act)

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

    call psb_nullify_sp(blck_)          ! Probably pointless.
    call psb_sp_all(0,0,blck_,1,info)
    if(info.ne.0) then
      info=4010
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    blck_%m=0
  endif

  !
  ! Compute the ILU(0) or the MILU(0) factorization, depending on ialg
  !
  call mld_zilu_fctint(ialg,m,a%m,a,blck_%m,blck_,&
       & d,l%aspk,l%ia1,l%ia2,u%aspk,u%ia1,u%ia2,l1,l2,info)
  if(info.ne.0) then
     info=4010
     ch_err='mld_zilu_fctint'
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
  !     Nullify pointer / deallocate memory
  !
  if (present(blck)) then 
    blck_ => null() 
  else
    call psb_sp_free(blck_,info)
    if(info.ne.0) then
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
  ! Subroutine: mld_zilu_fctint
  ! Version:    complex
  ! Note: internal subroutine of mld_zilu_fct
  !
  !  This routine computes either the ILU(0) or the MILU(0) factorization of the
  !  local part of the matrix stored into a. These factorizations are used to build
  !  the 'base preconditioner' (block-Jacobi preconditioner/solver, Additive Schwarz
  !  preconditioner) corresponding to a certain level of a multilevel preconditioner.
  !
  !  The local matrix to be factorized is stored into a and b, as specified in the
  !  description of the arguments below. The storage format for both the L and U
  !  factors is CSR. The diagonal of the U factor is stored separately (actually,
  !  the inverse of the diagonal entries is stored; this is then managed in the
  !  solve stage associated to the ILU(0)/MILU(0) factorization).
  !
  !  The routine copies and factors "on the fly" from the sparse matrix structures a
  !  and b into the arrays laspk, uaspk, d (L, U without its diagonal, diagonal of U).
  !  
  !
  ! Arguments:
  !    ialg    -  integer, input.
  !               The type of incomplete factorization to be performed.
  !               The MILU(0) factorization is computed if ialg = 2 (= mld_milu_n_);
  !               the ILU(0) factorization otherwise.
  !    m       -  integer, output.
  !               The total number of rows of the local matrix to be factorized,
  !               i.e. ma+mb.
  !    ma      -  integer, input
  !               The number of rows of the local submatrix stored into a.
  !    a       -  type(psb_zspmat_type), input.
  !               The sparse matrix structure containing the local matrix to be
  !               factorized. Note that, if the 'base' Additive Schwarz preconditioner
  !               has overlap greater than 0 and the matrix has not been reordered
  !               (see mld_bjac_bld), then a contains only the 'original' local part
  !               of the matrix to be factorized, i.e. the rows of the matrix held
  !               by the calling process according to the initial data distribution.
  !    mb      -  integer, input.
  !               The number of rows of the local submatrix stored into b.
  !    b       -  type(psb_zspmat_type), input.
  !               The sparse matrix structure containing the remote rows of the
  !               matrix to be factorized, that have been retrieved by mld_asmat_bld
  !               to build an Additive Schwarz base preconditioner with overlap
  !               greater than 0. If the overlap is 0 or the matrix has been
  !               reordered (see mld_bjac_bld), then b does not contain any row.
  !    d       -  complex(kind(1.d0)), dimension(:), output.
  !               The inverse of the diagonal entries of the U factor in the
  !               incomplete factorization.
  !    laspk   -  complex(kind(1.d0)), dimension(:), input/output.
  !               The entries of U are stored according to the CSR format.
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
  !    l1      -  integer, output.
  !               The number of nonzero entries in laspk.
  !    l2      -  integer, output.
  !               The number of nonzero entries in uaspk.
  !    info    -  integer, output.           
  !               Error code.
  !
  subroutine mld_zilu_fctint(ialg,m,ma,a,mb,b,&
       & d,laspk,lia1,lia2,uaspk,uia1,uia2,l1,l2,info)

    implicit none 

    ! Arguments
    integer, intent(in)            :: ialg
    type(psb_zspmat_type)          :: a,b
    integer                        :: m,ma,mb,l1,l2,info
    integer, dimension(:)          :: lia1,lia2,uia1,uia2
    complex(kind(1.d0)), dimension(:) :: laspk,uaspk,d

    ! Local variables
    integer :: i,j,k,l,low1,low2,kk,jj,ll, ktrw,err_act
    complex(kind(1.d0)) :: dia,temp
    integer, parameter :: nrb=16
    logical,parameter  :: debug=.false.
    type(psb_zspmat_type) :: trw
    integer             :: int_err(5) 
    character(len=20)   :: name, ch_err

    name='mld_zilu_fctint'
    if(psb_get_errstatus().ne.0) return 
    info=0
    call psb_erractionsave(err_act)
    call psb_nullify_sp(trw)
    trw%m=0
    trw%k=0
    if(debug) write(0,*)'LUINT Allocating TRW'
    call psb_sp_all(trw,1,info)
    if(info.ne.0) then
      info=4010
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    if(debug) write(0,*)'LUINT Done  Allocating TRW'
    lia2(1) = 1
    uia2(1) = 1
    l1      = 0
    l2      = 0
    m = ma+mb
    if(debug) write(0,*)'In DCSRLU Begin cycle',m,ma,mb

    !
    ! Cycle over the matrix rows
    !
    do i = 1, m

      if(debug) write(0,*)'LUINT: Loop index ',i,ma
      d(i) = zzero

      if (i <= ma) then
        !
        ! Copy the i-th local row of the matrix, stored in a,
        ! into laspk/d(i)/uaspk 
        !
        call ilu_copyin(i,ma,a,i,1,m,l1,lia1,laspk,&
             & d(i),l2,uia1,uaspk,ktrw,trw)
      else
        !
        ! Copy the i-th local row of the matrix, stored in b
        ! (as (i-ma)-th row), into laspk/d(i)/uaspk 
        !
        call ilu_copyin(i-ma,mb,b,i,1,m,l1,lia1,laspk,&
             & d(i),l2,uia1,uaspk,ktrw,trw)
      endif

      lia2(i+1) = l1 + 1
      uia2(i+1) = l2 + 1

      dia = d(i)
      do kk = lia2(i), lia2(i+1) - 1
        !
        ! Compute entry l(i,k) (lower factor L) of the incomplete
        ! factorization
        !
        temp      = laspk(kk)
        k         = lia1(kk)
        laspk(kk) = temp*d(k)
        !
        ! Update the rest of row i (lower and upper factors L and U)
        ! using l(i,k)
        !
        low1 = kk + 1
        low2 = uia2(i)
        !
        updateloop: do  jj = uia2(k), uia2(k+1) - 1
          !
          j = uia1(jj)
          !
          if (j < i) then
            !
            ! search l(i,*) (i-th row of L) for a matching index j
            !
            do  ll = low1, lia2(i+1) - 1
              l = lia1(ll)
              if (l > j) then
                low1 = ll
                exit
              else if (l == j) then
                laspk(ll) = laspk(ll) - temp*uaspk(jj)
                low1 = ll + 1
                cycle updateloop
              end if
            enddo

          else if (j == i) then
            !
            ! j=i: update the diagonal
            !
            dia = dia - temp*uaspk(jj)
            cycle updateloop
            !     
          else if (j > i) then
            !
            ! search u(i,*) (i-th row of U) for a matching index j
            !
            do ll = low2, uia2(i+1) - 1
              l = uia1(ll)
              if (l > j) then
                low2 = ll
                exit
              else if (l == j) then
                uaspk(ll) = uaspk(ll) - temp*uaspk(jj)
                low2 = ll + 1
                cycle updateloop
              end if
            enddo
          end if
          !     
          ! If we get here we missed the cycle updateloop, which means 
          ! that this entry does not match; thus we accumulate on the
          ! diagonal for MILU(0).
          !
          if (ialg == mld_milu_n_) then 
            dia = dia - temp*uaspk(jj)
          end if
        enddo updateloop
      enddo
      !     
      ! Check the pivot size
      !     
      if (abs(dia) < epstol) then
        !
        ! Too small pivot: unstable factorization
        !     
        info = 2
        int_err(1) = i
        write(ch_err,'(g20.10)') abs(dia)
        call psb_errpush(info,name,i_err=int_err,a_err=ch_err)
        goto 9999
      else
        !
        ! Compute 1/pivot
        !
        dia = done/dia
      end if
      d(i) = dia
      !
      ! Scale row i of upper triangle
      !
      do  kk = uia2(i), uia2(i+1) - 1
        uaspk(kk) = uaspk(kk)*dia
      enddo
    enddo

    call psb_sp_free(trw,info)
    if(info.ne.0) then
      info=4010
      ch_err='psb_sp_free'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    if(debug) write(0,*)'Leaving ilu_fct'

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine mld_zilu_fctint

  !
  ! Subroutine: ilu_copyin
  ! Version:    complex
  ! Note: internal subroutine of mld_zilu_fct
  !
  !  This routine copies a row of a sparse matrix A, stored in the psb_dspmat_type 
  !  data structure a, into the arrays laspk and uaspk and into the scalar variable
  !  dia, corresponding to the lower and upper triangles of A and to the diagonal
  !  entry of the row, respectively. The entries in laspk and uaspk are stored
  !  according to the CSR format; the corresponding column indices are stored in
  !  the arrays lia1 and uia1.
  !
  !  If the sparse matrix is in CSR format, a 'straight' copy   is performed;
  !  otherwise psb_sp_getblk is used to extract a block of rows, which is then
  !  copied into laspk, dia, uaspk row by row, through successive calls to
  !  ilu_copyin.
  !
  !  The routine is used by mld_zilu_fctin in the computation of the ILU(0)/MILU(0)
  !  factorization of a local sparse matrix.
  !  
  !  TODO: modify the routine to allow copying into output L and U that are
  !  already filled with indices; this would allow computing an ILU(k) pattern,
  !  then use the ILU(0) internal for subsequent calls with the same pattern. 
  !
  ! Arguments:
  !    i       -  integer, input.
  !               The local index of the row to be extracted from  the 
  !               sparse matrix structure a.
  !    m       -  integer, input.
  !               The number of rows of the local matrix stored into a.
  !    a       -  type(psb_zspmat_type), input.
  !               The sparse matrix structure containing the row to be copied.
  !    jd      -  integer, input.
  !               The column index of the diagonal entry of the row to be
  !               copied.
  !    jmin    -  integer, input.
  !               Minimum valid column index.
  !    jmax    -  integer, input.
  !               Maximum valid column index.
  !               The output matrix will contain a clipped copy taken from
  !               a(1:m,jmin:jmax).
  !    l1      -  integer, input/output.
  !               Pointer to the last occupied entry of laspk.
  !    lia1    -  integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the lower triangle
  !               copied in laspk row by row (see mld_zilu_fctint), according
  !               to the CSR storage format.
  !    laspk   -  complex(kind(1.d0)), dimension(:), input/output.
  !               The array where the entries of the row corresponding to the
  !               lower triangle are copied.
  !    dia     -  complex(kind(1.d0)), output.
  !               The diagonal entry of the copied row.
  !    l2      -  integer, input/output.
  !               Pointer to the last occupied entry of uaspk.
  !    uia1    -  integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the upper triangle
  !               copied in uaspk row by row (see mld_zilu_fctint), according
  !               to the CSR storage format.
  !    uaspk   -  complex(kind(1.d0)), dimension(:), input/output.
  !               The array where the entries of the row corresponding to the
  !               upper triangle are copied. 
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
  subroutine ilu_copyin(i,m,a,jd,jmin,jmax,l1,lia1,laspk,&
       & dia,l2,uia1,uaspk,ktrw,trw)

    use psb_base_mod

    implicit none

    ! Arguments
    type(psb_zspmat_type) :: a,trw
    integer               :: i,m,ktrw,jd,jmin,jmax,l1,l2
    integer               :: lia1(:),uia1(:)
    complex(kind(1.d0))      :: laspk(:), uaspk(:), dia

    ! Local variables
    integer               :: k,j,info,irb
    integer, parameter    :: nrb=16
    character(len=20), parameter  :: name='mld_dilu_fctint'
    character(len=20)             :: ch_err

    if (psb_get_errstatus() /= 0) return 
    info=0
    call psb_erractionsave(err_act)

    if (toupper(a%fida)=='CSR') then

      !
      ! Take a fast shortcut if the matrix is stored in CSR format
      !

      do j = a%ia2(i), a%ia2(i+1) - 1
        k = a%ia1(j)
        !           write(0,*)'KKKKK',k
        if ((k < jd).and.(k >= jmin)) then
          l1 = l1 + 1
          laspk(l1) = a%aspk(j)
          lia1(l1) = k
        else if (k == jd) then
          dia = a%aspk(j)
        else if ((k > jd).and.(k <= jmax)) then
          l2 = l2 + 1
          uaspk(l2) = a%aspk(j)
          uia1(l2) = k
        end if
      enddo

    else

      !
      ! Otherwise use psb_sp_getblk, slower but able (in principle) of 
      ! handling any format. In this case, a block of rows is extracted
      ! instead of a single row, for performance reasons, and these
      ! rows are copied one by one into laspk, dia, uaspk, through
      ! successive calls to ilu_copyin.
      !

      if ((mod(i,nrb) == 1).or.(nrb==1)) then 
        irb = min(m-i+1,nrb)
        call psb_sp_getblk(i,a,trw,info,lrw=i+irb-1)
        if(info.ne.0) then
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
        k = trw%ia2(ktrw)
        if ((k < jd).and.(k >= jmin)) then
          l1 = l1 + 1
          laspk(l1) = trw%aspk(ktrw)
          lia1(l1) = k
        else if (k == jd) then
          dia = trw%aspk(ktrw)
        else if ((k > jd).and.(k <= jmax)) then
          l2 = l2 + 1
          uaspk(l2) = trw%aspk(ktrw)
          uia1(l2) = k
        end if
        ktrw = ktrw + 1
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
  end subroutine ilu_copyin

end subroutine mld_zilu_fct
