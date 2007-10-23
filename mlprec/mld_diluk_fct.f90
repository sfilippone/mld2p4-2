!!$ 
!!$ 
!!$                    MD2P4
!!$    Multilevel Domain Decomposition Parallel Preconditioner Package for PSBLAS
!!$                      for 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$                       Daniela di Serafino    Second University of Naples
!!$                       Pasqua D'Ambra         ICAR-CNR                      
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the MD2P4 group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MD2P4 GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
subroutine mld_diluk_fct(fill_in,ialg,a,l,u,d,info,blck)
  
  !
  ! This routine copies and factors "on the fly" from A and BLCK
  ! into L/D/U. 
  !
  !
  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_diluk_fct
  implicit none
  !     .. Scalar Arguments ..
  integer, intent(in)                 :: fill_in, ialg
  integer, intent(out)                :: info
  !     .. Array Arguments ..
  type(psb_dspmat_type),intent(in)    :: a
  type(psb_dspmat_type),intent(inout) :: l,u
  type(psb_dspmat_type),intent(in), optional, target :: blck
  real(kind(1.d0)), intent(inout)     ::  d(:)
  !     .. Local Scalars ..
  real(kind(1.d0)) ::  dia, temp
  integer   ::  i, j, jj, k, kk, l1, l2, ll, low1, low2,m,ma,err_act
  
  type(psb_dspmat_type), pointer  :: blck_
  character(len=20)   :: name, ch_err
  logical, parameter :: debug=.false.

  name='mld_diluk_fct'
  info = 0
  call psb_erractionsave(err_act)
  !     .. Executable Statements ..
  !
  if (debug) write(0,*) 'mld_diluk_fct: start'
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
  if (debug) write(0,*) 'mld_diluk_fct: calling fctint'
  call mld_diluk_fctint(fill_in,ialg,m,a%m,a,blck_%m,blck_,&
       & d,l%aspk,l%ia1,l%ia2,u%aspk,u%ia1,u%ia2,l1,l2,info)
  if (info /= 0) then
     info=4010
     ch_err='mld_diluk_fctint'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

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
  subroutine mld_diluk_fctint(fill_in,ialg,m,ma,a,mb,b,&
       & d,laspk,lia1,lia2,uaspk,uia1,uia2,l1,l2,info)
    use psb_base_mod
    implicit none 
    integer, intent(in)            :: fill_in, ialg
    type(psb_dspmat_type)          :: a,b
    integer                        :: m,ma,mb,l1,l2,info
    integer, dimension(:), allocatable          :: lia1,lia2,uia1,uia2
    real(kind(1.d0)), dimension(:), allocatable :: laspk,uaspk
    real(kind(1.d0)), dimension(:)              :: d

    integer :: i,j,k,l,low1,low2,kk,jj,ll, ktrw,err_act, &
         & isz,minj,maxj,lrwk,nidx
    real(kind(1.d0)) :: dia,temp,rwk
    integer, allocatable          :: uplevs(:), rowlevs(:),idxs(:)
    real(kind(1.d0)), allocatable :: row(:)
    type(psb_int_heap) :: heap
    
    logical,parameter  :: debug=.false.
    type(psb_dspmat_type) :: trw
    integer             :: int_err(5) 
    character(len=20), parameter  :: name='mld_diluk_fctint'
    character(len=20)             :: ch_err

    if (psb_get_errstatus() /= 0) return 
    info=0
    call psb_erractionsave(err_act)

    m = ma+mb

    !
    ! Temp buffer for copyin function.  
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
    ! Initial allocation
    allocate(uplevs(size(uaspk)),rowlevs(m),row(m),stat=info)
    if (info /= 0) then
      info=4010
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999
    end if
    uplevs(:)  = m+1
    row(:)     = dzero
    rowlevs(:) = -(m+1)

    do i = 1, m
      if (debug.and.(mod(i,500)==1)) write(0,*)'LUINT: Loop index ',i,ma
      !
      ! At each iteration of the loop we keep the indices affected in a heap 
      ! initialized and filled in the copyin function, and updated during 
      ! the elimination. The heap is ideal because at each step we need the 
      ! lowest index, but we also need to insert new items, and the heap allows 
      ! to do both in log time. 
      !
      d(i) = dzero
      if (i<=ma) then 
        call iluk_copyin(i,ma,a,m,row,rowlevs,heap,ktrw,trw)
      else
        call iluk_copyin(i-ma,mb,b,m,row,rowlevs,heap,ktrw,trw)
      endif

      if (debug) write(0,*)'LUINT: input Copy done'
      ! Do an elimination step on current row
      ! Turns out we only need to keep track of levels
      ! for the upper triangle, hence no lowlevs variable.
      !
      call iluk_fact(fill_in,i,m,row,rowlevs,heap,&
           & d,uia1,uia2,uaspk,uplevs,nidx,idxs)
      !
      ! Copy the row into the lower/diag/upper structures.
      ! 
      call iluk_copyout(fill_in,ialg,i,m,row,rowlevs,nidx,idxs,&
           & l1,l2,lia1,lia2,laspk,d,uia1,uia2,uaspk,uplevs)

    end do

    !
    ! And we're done......hopefully :-) 
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
  end subroutine mld_diluk_fctint


  subroutine iluk_copyin(i,m,a,jmax,row,rowlevs,heap,ktrw,trw)
    use psb_base_mod
    implicit none 
    type(psb_dspmat_type) :: a,trw
    integer               :: i, rowlevs(:),m,ktrw,jmax
    real(kind(1.d0))      :: row(:)
    type(psb_int_heap)    :: heap
    
    integer               :: k,j,info,irb
    integer, parameter    :: nrb=16
    character(len=20), parameter  :: name='mld_diluk_fctint'
    character(len=20)             :: ch_err

    if (psb_get_errstatus() /= 0) return 
    info=0
    call psb_erractionsave(err_act)

    !
    ! Here we take a fast shortcut if possible, otherwise 
    ! use spgtblk, slower but able (in principle) to handle 
    ! anything. 
    !
    if (toupper(a%fida)=='CSR') then 
      call psb_init_heap(heap,info) 
      
      do j = a%ia2(i), a%ia2(i+1) - 1
        k          = a%ia1(j)
        if ((1<=k).and.(k<=jmax)) then 
          row(k)     = a%aspk(j)
          rowlevs(k) = 0
          call psb_insert_heap(k,heap,info)
        end if
      end do
    else
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
        if ((1<=k).and.(k<=jmax)) then 
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


  subroutine iluk_fact(fill_in,i,m,row,rowlevs,heap,d,uia1,uia2,uaspk,uplevs,nidx,idxs)
    use psb_base_mod
    implicit none 
    type(psb_dspmat_type) :: a
    type(psb_int_heap)    :: heap 
    integer               :: i,m, rowlevs(:),minj,maxj,fill_in,nidx
    integer, allocatable  :: idxs(:)
    integer               :: uia1(:),uia2(:),uplevs(:)
    real(kind(1.d0))      :: row(:), uaspk(:),d(:)

    integer               :: k,j,lrwk,jj,info, lastk
    real(kind(1.d0))      :: rwk

    ! Do an elimination step on current row
    ! Turns out we only need to keep track of levels
    ! for the upper triangle.
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
      
      if ((row(k) /= dzero).and.(rowlevs(k) <= fill_in).and.(k<i)) then 
        !
        ! Note: since U is scaled while copying out, we can use rwk 
        ! in the update below
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
          ! Insert the index for further processing.
          ! The levels are initialized to a negative value; if we find one,
          ! it means that it's an as yet untouched index, so we need to 
          ! insert it, otherwise it's already on the heap, no need to 
          ! insert more than once. 
          !
          if (rowlevs(j)<0) then 
            call psb_insert_heap(j,heap,info)
            rowlevs(j) = abs(rowlevs(j))
          end if
          row(j)     = row(j) - rwk * uaspk(jj)
          rowlevs(j) = min(rowlevs(j),lrwk+uplevs(jj)+1)
        end do

      end if
    end do

  end subroutine iluk_fact

  subroutine iluk_copyout(fill_in,ialg,i,m,row,rowlevs,nidx,idxs,&
       &  l1,l2,lia1,lia2,laspk,d,uia1,uia2,uaspk,uplevs)
    use psb_base_mod
    implicit none 
    integer               :: fill_in,ialg,i, rowlevs(:),minj,maxj,l1,l2,m,nidx,idxs(:)
    integer, allocatable  :: uia1(:),uia2(:), lia1(:),lia2(:),uplevs(:)
    real(kind(1.d0)),allocatable    :: uaspk(:), laspk(:)
    real(kind(1.d0))      :: row(:), d(:)
    integer               :: k,j,isz,info,err_act,int_err(5),idxp
    character(len=20), parameter  :: name='mld_diluk_fctint'
    character(len=20)             :: ch_err

    if (psb_get_errstatus() /= 0) return 
    info=0
    call psb_erractionsave(err_act)
    !
    ! Copy the row into the lower/diag/upper structures.
    ! 
    ! Lower part 
    d(i) = dzero
    do idxp=1,nidx
      j = idxs(idxp)
      if (j<i) then 
        if (rowlevs(j) <= fill_in) then 
          l1     = l1 + 1 
          if (size(laspk) < l1) then 
            ! Figure out a good reallocation size!! 
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
          d(i) = d(i) + row(j)
        end if
        row(j)     = dzero
        rowlevs(j) = -(m+1)
      else if (j==i) then 
        d(i)       = d(i) + row(i) 
        row(i)     = dzero
        rowlevs(i) = -(m+1)
      else if (j>i) then 
        ! Upper part 
        if (rowlevs(j) <= fill_in) then 
          l2     = l2 + 1 
          if (size(uaspk) < l2) then 
            ! Figure out a good reallocation size!! 
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
          d(i) = d(i) + row(j)
        end if
        row(j)     = dzero
        rowlevs(j) = -(m+1)
      end if
    end do

    lia2(i+1) = l1 + 1
    uia2(i+1) = l2 + 1

    if (abs(d(i)) < epstol) then
      !
      !     Pivot too small: unstable factorization
      !     
      info = 2
      int_err(1) = i
      write(ch_err,'(g20.10)') d(i)
      call psb_errpush(info,name,i_err=int_err,a_err=ch_err)
      goto 9999
    else
      d(i) = done/d(i)
    end if
    ! Now scale upper part
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


end subroutine mld_diluk_fct
