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
! File: mld_daggrmap_bld.f90.
!
! Subroutine: mld_daggrmap_bld.
! Version:    real.
!
!  This routine builds a mapping from the row indices of the fine-level matrix
!  to the row indices of the coarse-level matrix, according to a decoupled 
!  aggregation algorithm. This mapping will be used by mld_daggrmat_asb      to
!  build the coarse-level matrix.  
!
!  The aggregation algorithm is a parallel version of that described in
!    M. Brezina and P. Vanek, A black-box iterative solver based on a 
!    two-level Schwarz method, Computing,  63 (1999), 233-263.
!  For more details see
!    P. D'Ambra, D. di Serafino and S. Filippone, On the development of
!    PSBLAS-based parallel two-level Schwarz preconditioners, Appl. Num. Math.
!    57 (2007), 1181-1196.
!
!
! Arguments:
!    aggr_type  -  integer, input.
!                  The scalar used to identify the aggregation algorithm.
!    a          -  type(psb_dspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!    nlaggr     -  integer, dimension(:), allocatable.
!                  nlaggr(i) contains the aggregates held by process i.
!    ilaggr     -  integer, dimension(:), allocatable.
!                  The mapping between the row indices of the coarse-level
!                  matrix and the row indices of the fine-level matrix.
!                  ilaggr(i)=j means that node i in the adjacency graph
!                  of the fine-level matrix is mapped onto node j in the
!                  adjacency graph of the coarse-level matrix.
!    info       -  integer, output.
!                  Error code.
!
subroutine mld_daggrmap_bld(aggr_type,a,desc_a,nlaggr,ilaggr,info)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_daggrmap_bld
  
  implicit none

! Arguments
  integer, intent(in)               :: aggr_type
  type(psb_dspmat_type), intent(in), target  :: a
  type(psb_desc_type), intent(in)   :: desc_a
  integer, allocatable, intent(out) :: ilaggr(:),nlaggr(:)
  integer, intent(out)              :: info

! Local variables
  integer, allocatable  :: ils(:), neigh(:)
  integer :: icnt,nlp,k,n,ia,isz,nr, naggr,i,j,m
  type(psb_dspmat_type), target :: atmp, atrans
  type(psb_dspmat_type), pointer  :: apnt
  logical :: recovery
  logical, parameter :: debug=.false.
  integer :: ictxt,np,me,err_act
  integer :: nrow, ncol, n_ne
  integer, parameter :: one=1, two=2
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  name = 'mld_aggrmap_bld'
  call psb_erractionsave(err_act)
  !
  ! Note. At the time being we are ignoring aggr_type so
  ! that we only have decoupled aggregation. This might 
  ! change in the future. 
  !
  ictxt=psb_cd_get_context(desc_a)
  call psb_info(ictxt,me,np)
  nrow  = psb_cd_get_local_rows(desc_a)
  ncol  = psb_cd_get_local_cols(desc_a)

  select case (aggr_type)
  case (mld_dec_aggr_,mld_sym_dec_aggr_) 
    
    
    nr = a%m
    allocate(ilaggr(nr),neigh(nr),stat=info)
    if(info.ne.0) then
      info=4025
      call psb_errpush(info,name,i_err=(/2*nr,0,0,0,0/),&
           & a_err='integer')
      goto 9999
    end if

    do i=1, nr
      ilaggr(i) = -(nr+1)
    end do
    if (aggr_type == mld_dec_aggr_) then 
      apnt => a
    else
      call psb_sp_clip(a,atmp,info,imax=nr,jmax=nr,&
           & rscale=.false.,cscale=.false.)
      atmp%m=nr
      atmp%k=nr
      call psb_transp(atmp,atrans,fmt='COO')
      call psb_rwextd(nr,atmp,info,b=atrans,rowscale=.false.) 
      atmp%m=nr
      atmp%k=nr
      call psb_sp_free(atrans,info)
      call psb_ipcoo2csr(atmp,info)
      apnt => atmp
    end if


    ! Note: -(nr+1)  Untouched as yet
    !       -i    1<=i<=nr  Adjacent to aggregate i
    !        i    1<=i<=nr  Belonging to aggregate i

    !
    ! Phase one: group nodes together. 
    ! Very simple minded strategy. 
    ! 
    naggr = 0
    nlp   = 0
    do
      icnt = 0
      do i=1, nr 
        if (ilaggr(i) == -(nr+1)) then 
          !
          ! 1. Untouched nodes are marked >0 together 
          !    with their neighbours
          !
          icnt = icnt + 1 
          naggr = naggr + 1 
          ilaggr(i) = naggr

          call psb_neigh(apnt,i,neigh,n_ne,info,lev=one)
          if (info/=0) then 
            info=4010
            ch_err='psb_neigh'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
          do k=1, n_ne
            j = neigh(k)
            if ((1<=j).and.(j<=nr)) then 
              ilaggr(j) = naggr
            endif
          enddo
          !
          ! 2. Untouched neighbours of these nodes are marked <0.
          !
          call psb_neigh(apnt,i,neigh,n_ne,info,lev=two)
          if (info/=0) then 
            info=4010
            ch_err='psb_neigh'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if

          do n = 1, n_ne
            m = neigh(n)
            if ((1<=m).and.(m<=nr)) then
              if (ilaggr(m) == -(nr+1)) ilaggr(m) = -naggr
            endif
          enddo
        endif
      enddo
      nlp = nlp + 1 
      if (icnt == 0) exit 
    enddo
    if (debug) then 
      write(0,*) 'Check 1:',count(ilaggr == -(nr+1)),(a%ia1(i),i=a%ia2(1),a%ia2(2)-1)
    end if

    !
    ! Phase two: sweep over leftovers. 
    !
    allocate(ils(naggr+10),stat=info) 
    if(info.ne.0) then
      info=4025
      call psb_errpush(info,name,i_err=(/naggr+10,0,0,0,0/),&
           & a_err='integer')
      goto 9999
    end if

    do i=1, size(ils)
      ils(i) = 0
    end do
    do i=1, nr 
      n = ilaggr(i)
      if (n>0) then 
        if (n>naggr) then 
          write(0,*) 'loc_Aggregate: n > naggr 1 ? ',n,naggr
        else
          ils(n) = ils(n) + 1 
        end if

      end if
    end do
    if (debug) then 
      write(0,*) 'Phase 1: number of aggregates ',naggr
      write(0,*) 'Phase 1: nodes aggregated     ',sum(ils)
    end if

    recovery=.false.
    do i=1, nr
      if (ilaggr(i) < 0) then 
        !
        ! Now some silly rule to break ties:
        ! Group with smallest adjacent aggregate. 
        !
        isz = nr+1
        ia  = -1

        call psb_neigh(apnt,i,neigh,n_ne,info,lev=one)
        if (info/=0) then 
          info=4010
          ch_err='psb_neigh'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        end if

        do j=1, n_ne
          k = neigh(j)
          if ((1<=k).and.(k<=nr))  then 
            n = ilaggr(k) 
            if (n>0) then 
              if (n>naggr) then 
                write(0,*) 'loc_Aggregate: n > naggr 2? ',n,naggr
              end if

              if (ils(n) < isz) then 
                ia = n
                isz = ils(n)
              endif
            endif
          endif
        enddo
        if (ia == -1) then 
          if (ilaggr(i) > -(nr+1)) then
            ilaggr(i) = abs(ilaggr(i))
            if (ilaggr(I)>naggr) then 
              write(0,*) 'loc_Aggregate: n > naggr 3? ',ilaggr(i),naggr
            end if
            ils(ilaggr(i)) = ils(ilaggr(i)) + 1
            !
            ! This might happen if the pattern is non symmetric. 
            ! Need a better handling. 
            !
            recovery = .true.
          else
            write(0,*) 'Unrecoverable error !!',ilaggr(i), nr
          endif
        else
          ilaggr(i) = ia
          if (ia>naggr) then 
            write(0,*) 'loc_Aggregate: n > naggr 4? ',ia,naggr
          end if

          ils(ia)  = ils(ia) + 1
        endif
      end if
    enddo
    if (recovery) then 
      write(0,*) 'Had to recover from strange situation in loc_aggregate.'
      write(0,*) 'Perhaps an unsymmetric pattern?'
    endif
    if (debug) then 
      write(0,*) 'Phase 2: number of aggregates ',naggr
      write(0,*) 'Phase 2: nodes aggregated     ',sum(ils) 
      do i=1, naggr 
        write(*,*) 'Size of aggregate ',i,' :',count(ilaggr==i), ils(i)
      enddo
      write(*,*) maxval(ils(1:naggr))
      write(*,*) 'Leftovers ',count(ilaggr<0), '  in ',nlp,' loops'
    end if

!!$    write(0,*) 'desc_a loc_aggr 4 : ', desc_a%matrix_data(m_)
    if (count(ilaggr<0) >0) then 
      write(0,*) 'Fatal error: some leftovers!!!'
    endif

    deallocate(ils,neigh,stat=info)
    if (info/=0) then 
      info=4000
      call psb_errpush(info,name)
      goto 9999
    end if

    if (nrow /= size(ilaggr)) then 
      write(0,*) 'SOmething wrong ilaggr ',nrow,size(ilaggr)
    endif
    call psb_realloc(ncol,ilaggr,info)
    if (info/=0) then 
      info=4010
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    allocate(nlaggr(np),stat=info)
    if (info/=0) then 
      info=4025
      call psb_errpush(info,name,i_err=(/np,0,0,0,0/),&
           & a_err='integer')
      goto 9999
    end if

    nlaggr(:) = 0
    nlaggr(me+1) = naggr
    call psb_sum(ictxt,nlaggr(1:np))

    if (aggr_type == mld_sym_dec_aggr_) then 
      call psb_sp_free(atmp,info)
    end if

  case default

    write(0,*) 'Unimplemented aggregation algorithm ',aggr_type
    info = -1
    call psb_errpush(30,name,i_err=(/1,aggr_type,0,0,0/))
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_daggrmap_bld
