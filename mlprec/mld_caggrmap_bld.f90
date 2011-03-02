!!$ 
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010
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
! File: mld_caggrmap_bld.f90
!
! Subroutine: mld_caggrmap_bld
! Version:    complex
!
!  This routine builds a mapping from the row indices of the fine-level matrix
!  to the row indices of the coarse-level matrix, according to a decoupled 
!  aggregation algorithm. This mapping will be used by mld_aggrmat_asb to
!  build the coarse-level matrix.  
!
!  The aggregation algorithm is a parallel version of that described in
!  * M. Brezina and P. Vanek, A black-box iterative solver based on a 
!    two-level Schwarz method, Computing,  63 (1999), 233-263.
!  * P. Vanek, J. Mandel and M. Brezina, Algebraic Multigrid by Smoothed
!    Aggregation for Second and Fourth Order Elliptic Problems, Computing, 56
!    (1996), 179-196.
!  For more details see
!    P. D'Ambra, D. di Serafino and S. Filippone, On the development of
!    PSBLAS-based parallel two-level Schwarz preconditioners, Appl. Num. Math.
!    57 (2007), 1181-1196.
!
!
! Arguments:
!    aggr_type  -  integer, input.
!                  The scalar used to identify the aggregation algorithm.
!    theta      -  real, input.
!                  The aggregation threshold used in the aggregation algorithm.
!    a          -  type(psb_cspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!    ilaggr     -  integer, dimension(:), allocatable.
!                  The mapping between the row indices of the coarse-level
!                  matrix and the row indices of the fine-level matrix.
!                  ilaggr(i)=j means that node i in the adjacency graph
!                  of the fine-level matrix is mapped onto node j in the
!                  adjacency graph of the coarse-level matrix.
!    nlaggr     -  integer, dimension(:), allocatable.
!                  nlaggr(i) contains the aggregates held by process i.
!    info       -  integer, output.
!                  Error code.
!
subroutine mld_caggrmap_bld(aggr_type,theta,a,desc_a,ilaggr,nlaggr,info)

  use psb_sparse_mod
  use mld_c_inner_mod, mld_protect_name => mld_caggrmap_bld
  
  implicit none

  ! Arguments
  integer, intent(in)               :: aggr_type
  real(psb_spk_), intent(in)        :: theta
  type(psb_cspmat_type), intent(in) :: a
  type(psb_desc_type), intent(in)   :: desc_a
  integer, allocatable, intent(out) :: ilaggr(:),nlaggr(:)
  integer, intent(out)              :: info

  ! Local variables
  integer, allocatable  :: ils(:), neigh(:)
  integer :: icnt,nlp,k,n,ia,isz,nr, naggr,i,j,m
  type(psb_cspmat_type) :: atmp, atrans
  logical :: recovery
  integer :: debug_level, debug_unit
  integer :: ictxt,np,me,err_act
  integer :: nrow, ncol, n_ne
  character(len=20)  :: name, ch_err

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  name = 'mld_aggrmap_bld'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  !
  ictxt=psb_cd_get_context(desc_a)
  call psb_info(ictxt,me,np)
  nrow  = psb_cd_get_local_rows(desc_a)
  ncol  = psb_cd_get_local_cols(desc_a)

  select case (aggr_type)
  case (mld_dec_aggr_)  
    call mld_dec_map_bld(theta,a,desc_a,nlaggr,ilaggr,info)    

  case (mld_sym_dec_aggr_)  
    nr = a%get_nrows()
    call a%csclip(atmp,info,imax=nr,jmax=nr,&
         & rscale=.false.,cscale=.false.)
    call atmp%set_nrows(nr)
    call atmp%set_ncols(nr)
    if (info == psb_success_) call atrans%transp(atmp)
    if (info == psb_success_) call atrans%cscnv(info,type='COO')
    if (info == psb_success_) call psb_rwextd(nr,atmp,info,b=atrans,rowscale=.false.) 
    call atmp%set_nrows(nr)
    call atmp%set_ncols(nr)
    if (info == psb_success_) call atrans%free()
    if (info == psb_success_) call atmp%cscnv(info,type='CSR')
    if (info == psb_success_) call mld_dec_map_bld(theta,atmp,desc_a,nlaggr,ilaggr,info)    
    if (info == psb_success_) call atmp%free()

  case default

    info = -1
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=(/1,aggr_type,0,0,0/))
    goto 9999
  end select

  if (info /= psb_success_) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='dec_map_bld')
    goto 9999
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

contains

  subroutine mld_dec_map_bld(theta,a,desc_a,nlaggr,ilaggr,info)

    use psb_sparse_mod
    use mld_c_inner_mod !, mld_protect_name => mld_daggrmap_bld

    implicit none

    ! Arguments
    type(psb_cspmat_type), intent(in) :: a
    type(psb_desc_type), intent(in)   :: desc_a
    real(psb_spk_), intent(in)        :: theta
    integer, allocatable, intent(out) :: ilaggr(:),nlaggr(:)
    integer, intent(out)              :: info

    ! Local variables
    integer, allocatable  :: ils(:), neigh(:), irow(:), icol(:)
    complex(psb_spk_), allocatable  :: val(:), diag(:)
    integer :: icnt,nlp,k,n,ia,isz,nr, naggr,i,j,m, nz
    real(psb_spk_)  :: cpling, tcl
    logical :: recovery
    integer :: debug_level, debug_unit
    integer :: ictxt,np,me,err_act
    integer :: nrow, ncol, n_ne
    character(len=20)  :: name, ch_err

    if(psb_get_errstatus() /= 0) return 
    info=psb_success_
    name = 'mld_dec_map_bld'
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    !
    ictxt=psb_cd_get_context(desc_a)
    call psb_info(ictxt,me,np)
    nrow  = psb_cd_get_local_rows(desc_a)
    ncol  = psb_cd_get_local_cols(desc_a)

    nr = a%get_nrows()
    allocate(ilaggr(nr),neigh(nr),stat=info)
    if(info /= psb_success_) then
      info=psb_err_alloc_request_
      call psb_errpush(info,name,i_err=(/2*nr,0,0,0,0/),&
           & a_err='integer')
      goto 9999
    end if

    allocate(diag(nr),stat=info)
    if(info /= psb_success_) then
      info=psb_err_alloc_request_
      call psb_errpush(info,name,i_err=(/nr,0,0,0,0/),&
           & a_err='complex(psb_spk_)')
      goto 9999
    end if
    call a%get_diag(diag,info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_sp_getdiag')
      goto 9999
    end if

    do i=1, nr
      ilaggr(i) = -(nr+1)
    end do


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
          icnt      = icnt + 1 
          naggr     = naggr + 1 
          ilaggr(i) = naggr

          call a%csget(i,i,nz,irow,icol,val,info)
          if (info /= psb_success_) then 
            info=psb_err_from_subroutine_
            call psb_errpush(info,name,a_err='csget')
            goto 9999
          end if

          do k=1, nz
            j = icol(k)
            if ((1<=j).and.(j<=nr).and.(i /= j)) then 
              if (abs(val(k)) > theta*sqrt(abs(diag(i)*diag(j)))) then 
                ilaggr(j) = naggr
              else 
                ilaggr(j) = -naggr
              endif
            end if
          enddo

          !
          ! 2. Untouched neighbours of these nodes are marked <0.
          !
          call a%get_neigh(i,neigh,n_ne,info,lev=2)
          if (info /= psb_success_) then 
            info=psb_err_from_subroutine_
            call psb_errpush(info,name,a_err='psb_neigh')
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
    if (debug_level >= psb_debug_outer_) then 
      write(debug_unit,*) me,' ',trim(name),&
           & ' Check 1:',count(ilaggr == -(nr+1))
    end if

    !
    ! Phase two: sweep over leftovers. 
    !
    allocate(ils(nr),stat=info) 
    if(info /= psb_success_) then
      info=psb_err_alloc_request_
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
          info=psb_err_internal_error_
          call psb_errpush(info,name,a_err='loc_Aggregate: n > naggr 1 ?')
          goto 9999        
        else
          ils(n) = ils(n) + 1 
        end if

      end if
    end do
    if (debug_level >= psb_debug_outer_) then 
      write(debug_unit,*) me,' ',trim(name),&
           & 'Phase 1: number of aggregates ',naggr
      write(debug_unit,*) me,' ',trim(name),&
           & 'Phase 1: nodes aggregated     ',sum(ils)
    end if

    recovery=.false.
    do i=1, nr
      if (ilaggr(i) < 0) then 
        !
        ! Now some silly rule to break ties:
        ! Group with adjacent aggregate. 
        !
        isz  = nr+1
        ia   = -1
        cpling = szero
        call a%csget(i,i,nz,irow,icol,val,info)
        if (info /= psb_success_) then 
          info=psb_err_from_subroutine_
          call psb_errpush(info,name,a_err='psb_sp_getrow')
          goto 9999
        end if

        do j=1, nz
          k = icol(j)
          if ((1<=k).and.(k<=nr).and.(k /= i))  then 
            tcl = abs(val(j)) / sqrt(abs(diag(i)*diag(k)))
            if (abs(val(j)) > theta*sqrt(abs(diag(i)*diag(k)))) then 
!!$            if (tcl > theta) then 
              n = ilaggr(k) 
              if (n>0) then 
                if (n>naggr) then 
                  info=psb_err_internal_error_
                  call psb_errpush(info,name,a_err='loc_Aggregate: n > naggr 2 ?')
                  goto 9999        
                end if

                if ((abs(val(j))>cpling) .or. &
                     & ((abs(val(j)) == cpling).and. (ils(n) < isz))) then 
!!$                if ((tcl > cpling) .or. ((tcl  == cpling).and. (ils(n) < isz))) then
                  ia     = n
                  isz    = ils(n)
                  cpling = abs(val(j))
!!$                  cpling = tcl
                endif
              endif
            endif
          end if
        enddo

        if (ia == -1) then 
          ! At this point, the easiest thing is to start a new aggregate
          naggr          = naggr + 1
          ilaggr(i)      = naggr 
          ils(ilaggr(i)) = 1

        else

          ilaggr(i) = ia

          if (ia>naggr) then 
            info=psb_err_internal_error_
            call psb_errpush(info,name,a_err='loc_Aggregate: n > naggr ? ')
            goto 9999        
          end if
          ils(ia)  = ils(ia) + 1
        endif

      end if
    end do
    if (debug_level >= psb_debug_outer_) then 
      if (recovery) then 
        write(debug_unit,*) me,' ',trim(name),&
             & 'Had to recover from strange situation in loc_aggregate.'
        write(debug_unit,*) me,' ',trim(name),&
             & 'Perhaps an unsymmetric pattern?'
      endif
      write(debug_unit,*) me,' ',trim(name),&
           & 'Phase 2: number of aggregates ',naggr,sum(ils) 
      do i=1, naggr 
        write(debug_unit,*) me,' ',trim(name),&
             & 'Size of aggregate ',i,' :',count(ilaggr == i), ils(i)
      enddo
      write(debug_unit,*) me,' ',trim(name),&
           & maxval(ils(1:naggr))
      write(debug_unit,*) me,' ',trim(name),&
           & 'Leftovers ',count(ilaggr<0), '  in ',nlp,' loops'
    end if

    if (count(ilaggr<0) >0) then 
      info=psb_err_internal_error_
      call psb_errpush(info,name,a_err='Fatal error: some leftovers')
      goto 9999
    endif

    deallocate(ils,neigh,stat=info)
    if (info /= psb_success_) then 
      info=psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    end if
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
      call psb_errpush(info,name,i_err=(/np,0,0,0,0/),&
           & a_err='integer')
      goto 9999
    end if

    nlaggr(:) = 0
    nlaggr(me+1) = naggr
    call psb_sum(ictxt,nlaggr(1:np))

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine mld_dec_map_bld

end subroutine mld_caggrmap_bld
