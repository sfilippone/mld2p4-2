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
subroutine mld_zilut_fct(fill_in,thres,ialg,a,l,u,d,info,blck)
  
  !
  ! This routine copies and factors "on the fly" from A and BLCK
  ! into L/D/U. 
  !
  !
  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_zilut_fct
  implicit none
  !     .. Scalar Arguments ..
  integer, intent(in)                 :: fill_in, ialg
  real(kind(1.d0)), intent(in)        :: thres
  integer, intent(out)                :: info
  !     .. Array Arguments ..
  type(psb_zspmat_type),intent(in)    :: a
  type(psb_zspmat_type),intent(inout) :: l,u
  type(psb_zspmat_type),intent(in), optional, target :: blck
  complex(kind(1.d0)), intent(inout)  ::  d(:)
  !     .. Local Scalars ..
  complex(kind(1.d0)) ::  temp
  integer   ::  i, j, jj, k, kk, l1, l2, ll, low1, low2,m,ma,err_act
  
  type(psb_zspmat_type), pointer  :: blck_
  character(len=20)   :: name, ch_err
  logical, parameter :: debug=.false.

  name='mld_zilut_fct'
  info = 0
  call psb_erractionsave(err_act)
  !     .. Executable Statements ..
  !
  if (debug) write(0,*) 'mld_zilut_fct: start'
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
  if (debug) write(0,*) 'mld_zilut_fct: calling fctint'
  call mld_zilut_fctint(fill_in,thres,ialg,m,a%m,a,blck_%m,blck_,&
       & d,l%aspk,l%ia1,l%ia2,u%aspk,u%ia1,u%ia2,l1,l2,info)
  if (info /= 0) then
     info=4010
     ch_err='mld_zilut_fctint'
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
  subroutine mld_zilut_fctint(fill_in,thres,ialg,m,ma,a,mb,b,&
       & d,laspk,lia1,lia2,uaspk,uia1,uia2,l1,l2,info)
    use psb_base_mod
    implicit none 
    integer, intent(in)            :: fill_in, ialg
    real(kind(1.d0)), intent(in)   :: thres
    type(psb_zspmat_type)          :: a,b
    integer                        :: m,ma,mb,l1,l2,info
    integer, dimension(:), allocatable          :: lia1,lia2,uia1,uia2
    complex(kind(1.d0)), dimension(:), allocatable :: laspk,uaspk
    complex(kind(1.d0)), dimension(:)              :: d

    integer :: i,j,k,l,low1,low2,kk,jj,ll, ktrw,err_act, &
         & isz,minj,maxj,lrwk,nidx,nlw,nup,jmaxup
    complex(kind(1.d0)) :: temp,rwk
    real(kind(1.d0)) :: nrmi
    integer, allocatable          :: idxs(:)
    complex(kind(1.d0)), allocatable :: row(:)
    type(psb_int_heap) :: heap
    
    logical,parameter  :: debug=.false.
    type(psb_zspmat_type) :: trw
    integer             :: int_err(5) 
    character(len=20), parameter  :: name='mld_zilut_fctint'
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
    allocate(row(m),stat=info)
    if (info /= 0) then
      info=4010
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999
    end if

    row(:) = zzero

    do i = 1, m
      if (debug) write(0,*)'LUINT: Loop index ',i
      !
      ! At each iteration of the loop we keep the indices affected in a heap 
      ! initialized and filled in the copyin function, and updated during 
      ! the elimination. The heap is ideal because at each step we need the 
      ! lowest index, but we also need to insert new items, and the heap allows 
      ! to do both in log time. 
      !
      d(i) = zzero
      if (i<=ma) then 
        call ilut_copyin(i,ma,a,i,m,nlw,nup,jmaxup,nrmi,row,heap,ktrw,trw)
      else
        call ilut_copyin(i-ma,mb,b,i,m,nlw,nup,jmaxup,nrmi,row,heap,ktrw,trw)
      endif

      if (debug) write(0,*)'LUINT: input Copy done'
      !
      ! Do an elimination step on current row
      !
      call ilut_fact(fill_in,thres,i,m,nrmi,row,heap,&
           & d,uia1,uia2,uaspk,nidx,idxs)
      !
      ! Copy the row into the lower/diag/upper structures.
      ! 
      call ilut_copyout(fill_in,thres,i,m,nlw,nup,jmaxup,nrmi,row,nidx,idxs,&
           & l1,l2,lia1,lia2,laspk,d,uia1,uia2,uaspk)

    end do

    !
    ! And we're done......hopefully :-) 
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


  subroutine ilut_copyin(i,m,a,jd,jmax,nlw,nup,jmaxup,nrmi,row,heap,ktrw,trw)
    use psb_base_mod
    implicit none 
    type(psb_zspmat_type) :: a,trw
    integer               :: i, m,ktrw,jmax,jd,nlw,nup,jmaxup
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

    !
    ! nrmi is the norm of the current sparse row. 
    ! For the time being, use the 2-norm.
    ! Note: below the 2-norm includes also elements 
    ! that are outside [1..JMAX] strictly. 
    ! Is this really important?? TO BE CHECKED.
    !
    nlw    = 0
    nup    = 0
    jmaxup = 0
    dmaxup = dzero
    nrmi   = dzero
    if (toupper(a%fida)=='CSR') then 
      call psb_init_heap(heap,info) 
      do j = a%ia2(i), a%ia2(i+1) - 1
        k          = a%ia1(j)
        if ((1<=k).and.(k<=jmax)) then 
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
        if ((1<=k).and.(k<=jmax)) then 
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


  subroutine ilut_fact(fill_in,thres,i,m,nrmi,row,heap,&
       & d,uia1,uia2,uaspk,nidx,idxs)
    use psb_base_mod
    implicit none 
    type(psb_zspmat_type) :: a
    type(psb_int_heap)    :: heap 
    integer               :: i,m,fill_in,nidx
    real(kind(1.d0)), intent(in)   :: thres,nrmi
    integer, allocatable  :: idxs(:)
    integer               :: uia1(:),uia2(:)
    complex(kind(1.d0))      :: row(:), uaspk(:),d(:)

    integer               :: k,j,lrwk,jj,info, lastk
    complex(kind(1.d0))   :: rwk
    logical, parameter    :: debug=.false.

    ! Do an elimination step on current row
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
        ! Dropping rule 1: compare row(k) with thres*2-norm of row
        !
        rwk    = row(k)
        row(k) = row(k) * d(k)
        
        if (abs(row(k)) < thres*nrmi) then 
          ! drop the element
          row(k) = zzero
          ! We just dropped this index, no need to insert it.
          cycle 
        else
          
          ! Note: since U is scaled while copying out, we can use rwk 
          ! in the update below
          !           
          do jj=uia2(k),uia2(k+1)-1
            j = uia1(jj)
            if (j<=k) then 
              write(0,*) 'Error in accessing upper mat???',j,k,jj
            endif
            !
            ! Insert the index for further processing.
            ! Is there a sensible way to prune the insertion? 
            !
            row(j)     = row(j) - rwk * uaspk(jj)
            if (abs(row(j)) < thres*nrmi) then 
              ! drop the element
              row(j) = zzero
            else
              call psb_insert_heap(j,heap,info)
            endif
          end do
        end if
      end if lowert
      !
      ! If we get here it's an index we need to keep on copyout.
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

  subroutine ilut_copyout(fill_in,thres,i,m,nlw,nup,jmaxup,nrmi,row, &
       & nidx,idxs,l1,l2,lia1,lia2,laspk,d,uia1,uia2,uaspk)
    use psb_base_mod
    implicit none 
    integer               :: fill_in,i, minj,maxj,l1,l2,m,nidx,idxs(:)
    integer               :: nlw,nup,jmaxup
    integer, allocatable  :: uia1(:),uia2(:), lia1(:),lia2(:)
    real(kind(1.d0)), intent(in) :: thres,nrmi
    complex(kind(1.d0)),allocatable :: uaspk(:), laspk(:)
    complex(kind(1.d0))             :: row(:), d(:)
    complex(kind(1.d0)),allocatable :: xw(:),xt(:)
    integer, allocatable         :: xwid(:), indx(:)
    complex(kind(1.d0))             :: witem
    integer                      :: widx
    integer                      :: k,j,isz,info,err_act,int_err(5),idxp, nz
    type(psb_dcomplex_idx_heap)  :: heap
    character(len=20), parameter :: name='mld_zilut_fctint'
    character(len=20)            :: ch_err
    logical                      :: fndmaxup
    logical, parameter           :: debug=.false.

    if (psb_get_errstatus() /= 0) return 
    info=0
    call psb_erractionsave(err_act)
    !
    ! Copy the row into the lower/diag/upper structures.
    ! 

    !
    ! Here we need to apply another dropping rule. We do it 
    ! by putting the nonzero elements in heaps, then copying 
    ! them out. 
    !
    ! The heap goes down on absolute value, so the first item 
    ! is the largest absolute value. 
    !
    call psb_init_heap(heap,info,dir=psb_asort_down_)
    allocate(xwid(nidx),xw(nidx),indx(nidx))


    ! First the lower part. 
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
      ! Dropping rule on 2-norm
      if (abs(witem) < thres*nrmi) cycle 
      nz       = nz + 1 
      xw(nz)   = witem 
      xwid(nz) = widx
      call psb_insert_heap(witem,widx,heap,info)
    end do
    ! Now have to take out the first elements. 
    if (nz <= nlw+fill_in) then 
      ! Just copy everything from xw, and it's already ordered
!!$      write(0,*) 'Lower triang copying everything ',i,nz
    else
      nz = nlw+fill_in
      do k=1,nz
        call psb_heap_get_first(witem,widx,heap,info)
        xw(k)   = witem
        xwid(k) = widx
      end do
    end if
    
    ! Now put things back into ascending column order
    call psb_msort(xwid(1:nz),indx(1:nz),dir=psb_sort_up_)
    do k=1,nz
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
      lia1(l1)   = xwid(k)
      laspk(l1)  = xw(indx(k))
    end do
    
    !
    ! Make sure idxp now points to the diagonal element
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
        ! Copy the diagonal
        widx      = idxs(idxp)
        witem     = row(widx)
        d(i)      = witem
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
          d(i) = zone/d(i)
        end if
      end if
    end if

    !
    ! Now for the upper part 
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
      ! Dropping rule on 2-norm. But keep the entry of jmaxup anyway

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
    
    ! Now have to take out the first elements. But 
    ! make sure we include jmaxup
    if (nz <= nup+fill_in) then 
      ! Just copy everything from xw
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
        ! Include element at jmaxup, if not already there.
        ! Put it in place of smallest coefficient 
        !
        xw(nz)   = row(jmaxup) 
        xwid(nz) = jmaxup
      endif
    end if
    ! Now put things back into ascending column order
    call psb_msort(xwid(1:nz),indx(1:nz),dir=psb_sort_up_)
    if (debug) then 
      write(0,*) 'Row ',i,' copyout: after sort at upper:',nz,jmaxup
      write(0,*) xwid(1:nz)
      write(0,*) xw(indx(1:nz))
    end if
    ! Upper part 
    do k=1,nz
      l2     = l2 + 1 
      if (size(uaspk) < l2) then 
        ! Figure out a good reallocation size!! 
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

    do idxp=1,nidx
      row(idxs(idxp)) = zzero
    end do
    
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
