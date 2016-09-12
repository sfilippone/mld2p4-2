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

subroutine mld_z_dec_map_bld(iorder,theta,a,desc_a,nlaggr,ilaggr,info)

  use psb_base_mod
  use mld_z_inner_mod, mld_protect_name => mld_z_dec_map_bld

  implicit none

  ! Arguments
  integer(psb_ipk_), intent(in)     :: iorder
  type(psb_zspmat_type), intent(in) :: a
  type(psb_desc_type), intent(in)    :: desc_a
  real(psb_dpk_), intent(in)         :: theta
  integer(psb_ipk_), allocatable, intent(out)  :: ilaggr(:),nlaggr(:)
  integer(psb_ipk_), intent(out)               :: info

  ! Local variables
  integer(psb_ipk_), allocatable  :: ils(:), neigh(:), irow(:), icol(:),&
       & ideg(:), idxs(:), tmpaggr(:)
  complex(psb_dpk_), allocatable  :: val(:), diag(:)
  integer(psb_ipk_) :: icnt,nlp,k,n,ia,isz,nr, naggr,i,j,m, nz, ilg, ii, ip
  type(psb_z_csr_sparse_mat) :: acsr
  real(psb_dpk_)  :: cpling, tcl
  logical :: recovery, disjoint
  integer(psb_ipk_) :: debug_level, debug_unit,err_act
  integer(psb_ipk_) :: ictxt,np,me
  integer(psb_ipk_) :: nrow, ncol, n_ne
  character(len=20)  :: name, ch_err

  if (psb_get_errstatus() /= 0) return 
  info=psb_success_
  name = 'mld_dec_map_bld'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  !
  ictxt=desc_a%get_context()
  call psb_info(ictxt,me,np)
  nrow  = desc_a%get_local_rows()
  ncol  = desc_a%get_local_cols()

  nr = a%get_nrows()
  allocate(ilaggr(nr),neigh(nr),ideg(nr),idxs(nr),stat=info)
  if(info /= psb_success_) then
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/2*nr,izero,izero,izero,izero/),&
         & a_err='integer')
    goto 9999
  end if

  diag = a%get_diag(info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_sp_getdiag')
    goto 9999
  end if

  if (iorder == mld_aggr_ord_nat_) then 
    do i=1, nr
      ilaggr(i) = -(nr+1)
      idxs(i)   = i 
    end do
  else 
    call a%cp_to(acsr)
    do i=1, nr
      ilaggr(i) = -(nr+1)
      ideg(i)   = acsr%irp(i+1) - acsr%irp(i)
    end do
    call acsr%free()
    call psb_msort(ideg,ix=idxs,dir=psb_sort_down_)
  end if


  !
  ! Phase one: Start with disjoint groups.
  ! 
  naggr = 0
  icnt = 0
  step1: do ii=1, nr
    i = idxs(ii)

    if (ilaggr(i) == -(nr+1)) then 
      call a%csget(i,i,nz,irow,icol,val,info)
      if (info /= psb_success_) then 
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='csget')
        goto 9999
      end if

      !
      ! Build the set of all strongly coupled nodes 
      !
      ip = 0 
      do k=1, nz
        j   = icol(k)
        if ((1<=j).and.(j<=nr)) then 
          if (abs(val(k)) > theta*sqrt(abs(diag(i)*diag(j)))) then
            ip = ip + 1
            icol(ip) = icol(k)
          end if
        end if
      enddo
      if (ip < 1) then
        write(0,*) "Should at least contain the node itself ! "
        cycle step1
      end if

      !
      ! If the whole strongly coupled neighborhood of I is
      ! as yet unconnected, turn it into the next aggregate
      !
      disjoint = all(ilaggr(icol(1:ip)) == -(nr+1))
      if (disjoint) then 
        icnt      = icnt + 1 
        naggr     = naggr + 1
        do k=1, ip
          ilaggr(icol(k)) = naggr
        end do
        ilaggr(i) = naggr
      end if
    endif
  enddo step1
  
  if (debug_level >= psb_debug_outer_) then 
    write(debug_unit,*) me,' ',trim(name),&
         & ' Check 1:',count(ilaggr == -(nr+1))
  end if

  !
  ! Phase two: join the neighbours
  !
  tmpaggr = ilaggr
  step2: do ii=1,nr
    i = idxs(ii)

    if (ilaggr(i) == -(nr+1)) then         
      call a%csget(i,i,nz,irow,icol,val,info)
      if (info /= psb_success_) then 
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_sp_getrow')
        goto 9999
      end if
      !
      ! Find the most strongly connected neighbour that is
      ! already aggregated, if any, and join its aggregate
      !
      cpling = dzero
      ip = 0
      do k=1, nz
        j   = icol(k)
        if ((1<=j).and.(j<=nr)) then 
          if ((abs(val(k)) > theta*sqrt(abs(diag(i)*diag(j))))&
               & .and. (tmpaggr(j) > 0).and. (abs(val(k)) > cpling)) then
            ip = k
            cpling = abs(val(k))
          end if
        end if
      enddo
      if (ip > 0) then 
        ilaggr(i) = ilaggr(icol(ip))
      end if
    end if
  end do step2


  !
  ! Phase three: sweep over leftovers, if any 
  !
  step3: do ii=1,nr
    i = idxs(ii)

    if (ilaggr(i) < 0) then         
      call a%csget(i,i,nz,irow,icol,val,info)
      if (info /= psb_success_) then 
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_sp_getrow')
        goto 9999
      end if
      !
      ! Find its strongly  connected neighbourhood not 
      ! already aggregated, and make it into a new aggregate.
      !
      cpling = dzero
      ip = 0
      do k=1, nz
        j   = icol(k)
        if ((1<=j).and.(j<=nr)) then 
          if ((abs(val(k)) > theta*sqrt(abs(diag(i)*diag(j))))&
               & .and. (ilaggr(j) < 0))  then
            ip = ip + 1
            icol(ip) = icol(k)
          end if
        end if
      enddo
      if (ip > 0) then
        icnt      = icnt + 1 
        naggr     = naggr + 1
        ilaggr(i) = naggr
        do k=1, ip
          ilaggr(icol(k)) = naggr
        end do
      else
        !
        ! This should not happen: we did not even connect with ourselves.
        ! Create an isolate anyway.
        !
        naggr     = naggr + 1
        ilaggr(i) = naggr
      end if
    end if
  end do step3


  if (count(ilaggr<0) >0) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Fatal error: some leftovers')
    goto 9999
  endif

  if (naggr > ncol) then 
    write(0,*) name,'Error : naggr > ncol',naggr,ncol
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Fatal error: naggr>ncol')
    goto 9999
  end if

  call psb_realloc(ncol,ilaggr,info)
  if (info /= psb_success_) then 
    info=psb_err_from_subroutine_
    ch_err='psb_realloc'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  allocate(nlaggr(np),stat=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/np,izero,izero,izero,izero/),&
         & a_err='integer')
    goto 9999
  end if

  nlaggr(:) = 0
  nlaggr(me+1) = naggr
  call psb_sum(ictxt,nlaggr(1:np))

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine mld_z_dec_map_bld

