!!$ 
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010,2012,2013
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
! File: mld_cmlprec_bld.f90
!
! Subroutine: mld_cmlprec_bld
! Version:    complex
!
!  This routine builds the preconditioner according to the requirements made by
!  the user trough the subroutines mld_precinit and mld_precset.
!  
!  A multilevel preconditioner is regarded as an array of 'one-level' data structures,
!  each containing the part of the preconditioner associated to a certain level,
!  (for more details see the description of mld_Tonelev_type in mld_prec_type.f90).
!  The levels are numbered in increasing order starting from the finest one, i.e.
!  level 1 is the finest level. No transfer operators are associated to level 1.
! 
!
! Arguments:
!    a       -  type(psb_cspmat_type).
!               The sparse matrix structure containing the local part of the
!               matrix to be preconditioned.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    p       -  type(mld_cprec_type), input/output.
!               The preconditioner data structure containing the local part
!               of the preconditioner to be built.
!    info    -  integer, output.
!               Error code.              
!
!    amold   -  class(psb_c_base_sparse_mat), input, optional
!               Mold for the inner format of matrices contained in the
!               preconditioner
!
!
!    vmold   -  class(psb_c_base_vect_type), input, optional
!               Mold for the inner format of vectors contained in the
!               preconditioner
!
!
!  
subroutine mld_cmlprec_bld(a,desc_a,p,info,amold,vmold)

  use psb_base_mod
  use mld_c_inner_mod, mld_protect_name => mld_cmlprec_bld
  use mld_c_prec_mod

  Implicit None

  ! Arguments
  type(psb_cspmat_type),intent(in), target           :: a
  type(psb_desc_type), intent(in), target              :: desc_a
  type(mld_cprec_type),intent(inout),target          :: p
  integer(psb_ipk_), intent(out)                       :: info
  class(psb_c_base_sparse_mat), intent(in), optional :: amold
  class(psb_c_base_vect_type), intent(in), optional  :: vmold
!!$  character, intent(in), optional         :: upd

  ! Local Variables
  type(mld_cprec_type) :: t_prec
  integer(psb_ipk_)      :: ictxt, me,np
  integer(psb_ipk_)      :: err,i,k, err_act, iszv, newsz, casize
  integer(psb_ipk_)      :: ipv(mld_ifpsz_), val
  integer(psb_ipk_)      :: int_err(5)
  character    :: upd_
  class(mld_c_base_smoother_type), allocatable :: coarse_sm, base_sm, med_sm
  type(mld_sml_parms) :: baseparms, medparms, coarseparms
  type(mld_c_onelev_node), pointer :: head, tail, newnode, current
  integer(psb_ipk_)            :: debug_level, debug_unit
  character(len=20)  :: name, ch_err

  if (psb_get_errstatus().ne.0) return 
  info=psb_success_
  err=0
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  name = 'mld_cmlprec_bld'
  info = psb_success_
  int_err(1) = 0
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Entering '
  !
  ! For the time being we are commenting out the UPDATE argument
  ! we plan to resurrect it later. 
!!$  if (present(upd)) then 
!!$    if (debug_level >= psb_debug_outer_) &
!!$         & write(debug_unit,*) me,' ',trim(name),'UPD ', upd
!!$
!!$    if ((psb_toupper(upd).eq.'F').or.(psb_toupper(upd).eq.'T')) then
!!$      upd_=psb_toupper(upd)
!!$    else
!!$      upd_='F'
!!$    endif
!!$  else
!!$    upd_='F'
!!$  endif
  upd_ = 'F'

  if (.not.allocated(p%precv)) then 
    !! Error: should have called mld_cprecinit
    info=3111
    call psb_errpush(info,name)
    goto 9999
  end if

  !
  ! Check to ensure all procs have the same 
  !   
  newsz = -1
  casize = p%coarse_aggr_size
  iszv  = size(p%precv)
  call psb_bcast(ictxt,iszv)
  call psb_bcast(ictxt,casize)
  if (casize > 0) then 
    if (casize /= p%coarse_aggr_size) then 
      info=psb_err_internal_error_
      call psb_errpush(info,name,a_err='Inconsistent coarse_aggr_size')
      goto 9999
    end if
  else
    if (iszv /= size(p%precv)) then 
      info=psb_err_internal_error_
      call psb_errpush(info,name,a_err='Inconsistent size of precv')
      goto 9999
    end if
  end if

  if (iszv <= 1) then
    ! We should only ever get here for multilevel.
    info=psb_err_from_subroutine_
    ch_err='size bpv'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  endif



  if (casize>0) then 
    ! 
    ! New strategy to build according to coarse size. 
    !
    coarseparms = p%precv(iszv)%parms
    baseparms   = p%precv(1)%parms
    medparms    = p%precv(2)%parms

    allocate(coarse_sm, source=p%precv(iszv)%sm,stat=info) 
    if (info == psb_success_)  &
         & allocate(med_sm, source=p%precv(2)%sm,stat=info) 
    if (info == psb_success_)  &
         & allocate(base_sm, source=p%precv(1)%sm,stat=info) 
    if (info /= psb_success_) then 
      write(0,*) 'Error in saving smoothers',info
      call psb_errpush(psb_err_internal_error_,name,a_err='Base level precbuild.')
      goto 9999
    end if
    !
    ! Replicated matrix should only ever happen at coarse level.
    !
    call mld_check_def(baseparms%coarse_mat,'Coarse matrix',&
         &   mld_distr_mat_,is_distr_ml_coarse_mat)
    !
    !    Now build a doubly linked list
    !
    allocate(newnode,stat=info) 
    if (info /= psb_success_) then 
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='List start'); goto 9999
    end if
    head => newnode
    tail => newnode
    newnode%item%base_a    => a
    newnode%item%base_desc => desc_a
    newnode%item%parms     = baseparms
    newsz = 1
    current => head
    list_build_loop: do 
      allocate(newnode,stat=info) 
      if (info /= psb_success_) then 
        info = psb_err_internal_error_
        call psb_errpush(info,name,a_err='List start'); goto 9999
      end if
      current%next => newnode
      newnode%prev => current
      newsz        = newsz + 1
      newnode%item%parms             = medparms
      newnode%item%parms%aggr_thresh = current%item%parms%aggr_thresh/2
      call mld_coarse_bld(current%item%base_a, current%item%base_desc, &
           & newnode%item,info)
      if (info /= psb_success_) then 
        info = psb_err_internal_error_
        call psb_errpush(info,name,a_err='build next level'); goto 9999
      end if

      if (newsz>2) then 
        if (all(current%item%map%naggr == newnode%item%map%naggr)) then 
          !
          ! We are not gaining anything
          !
          newsz = newsz-1
          current%next => null()
          call newnode%item%free(info)
          if (info == psb_success_) deallocate(newnode,stat=info)
          if (info /= psb_success_) then 
            info = psb_err_internal_error_
            call psb_errpush(info,name,a_err='Deallocate at list end'); goto 9999
          end if
          exit list_build_loop          
        end if
      end if

      current => current%next
      tail    => current
      if (sum(newnode%item%map%naggr) <= casize) then 
        !
        ! Target reached; but we may need to rebuild. 
        !
        exit list_build_loop
      end if
    end do list_build_loop
    !
    ! At this point, we are at  the list tail,
    ! and it needs to be rebuilt in case the parms were
    ! different.
    !
    ! But the threshold has to be fixed before rebuliding
    coarseparms%aggr_thresh = current%item%parms%aggr_thresh
    current%item%parms      = coarseparms
    call mld_coarse_bld(current%prev%item%base_a,&
         & current%prev%item%base_desc, &
         & current%item,info)
    if (info /= psb_success_) then 
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='build next level'); goto 9999
    end if

    !
    ! Ok, now allocate the output vector and fix items. 
    !
    do i=1,iszv
      if (info == psb_success_) call p%precv(i)%free(info)
    end do
    if (info == psb_success_) deallocate(p%precv,stat=info)
    if (info == psb_success_) allocate(p%precv(newsz),stat=info) 
    if (info /= psb_success_) then 
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='Reallocate precv'); goto 9999
    end if
    newnode => head
    do i=1, newsz
      current => newnode
      ! First do a move_alloc.
      ! This handles the AC, DESC_AC and MAP fields
      if (info == psb_success_) &
           & call mld_move_alloc(current%item,p%precv(i),info)
      ! Now set the smoother/solver parts. 
      if (info == psb_success_) then 
        if (i ==1) then 
          ! This is a workaround for a bug in gfortran 4.7.2
          call doallc(i,p%precv,base_sm,info) 
          ! !$          allocate(p%precv(i)%sm,source=base_sm,stat=info) 
        else if (i < newsz) then 
          call doallc(i,p%precv,med_sm,info) 
          ! !$          allocate(p%precv(i)%sm,source=med_sm,stat=info) 
        else
          call doallc(i,p%precv,coarse_sm,info) 
          ! !$          allocate(p%precv(i)%sm,source=coarse_sm,stat=info) 
        end if
      end if
      if (info /= psb_success_) then 
        info = psb_err_internal_error_
        call psb_errpush(info,name,a_err='list cpy'); goto 9999
      end if
      if (i == 1) then 
        p%precv(i)%base_a       => a
        p%precv(i)%base_desc    => desc_a
      else
        p%precv(i)%base_a       => p%precv(i)%ac
        p%precv(i)%base_desc    => p%precv(i)%desc_ac
        p%precv(i)%map%p_desc_X => p%precv(i-1)%base_desc
        p%precv(i)%map%p_desc_Y => p%precv(i)%base_desc
      end if

      newnode => current%next
      deallocate(current) 
    end do
    call base_sm%free(info)
    if (info == psb_success_) call med_sm%free(info)
    if (info == psb_success_) call coarse_sm%free(info)
    if (info == psb_success_) deallocate(coarse_sm,med_sm,base_sm,stat=info)
    if (info /= psb_success_) then 
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='final cleanup'); goto 9999
    end if
    iszv = newsz

  else 
    ! 
    ! Default, oldstyle
    ! 

    !
    ! Build the matrix and the transfer operators corresponding
    ! to the remaining levels
    !
    !
    ! Check on the iprcparm contents: they should be the same
    ! on all processes.
    !
    call psb_bcast(ictxt,p%precv(1)%parms)
    !
    ! Finest level first; remember to fix base_a and base_desc
    ! 
    p%precv(1)%base_a    => a
    p%precv(1)%base_desc => desc_a


    do i=2, iszv
      !
      ! Check on the iprcparm contents: they should be the same
      ! on all processes.
      !
      call psb_bcast(ictxt,p%precv(i)%parms)

      !
      ! Sanity checks on the parameters
      !
      if (i<iszv) then 
        !
        ! A replicated matrix only makes sense at the coarsest level
        !
        call mld_check_def(p%precv(i)%parms%coarse_mat,'Coarse matrix',&
             &   mld_distr_mat_,is_distr_ml_coarse_mat)
      end if

      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & 'Calling mlprcbld at level  ',i
      !
      ! Build the mapping between levels i-1 and i and the matrix
      ! at level i
      ! 
      if (info == psb_success_) call mld_coarse_bld(p%precv(i-1)%base_a,&
           & p%precv(i-1)%base_desc, p%precv(i),info)

      if (info /= psb_success_) then 
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Init upper level preconditioner')
        goto 9999
      endif

      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & 'Return from ',i,' call to mlprcbld ',info      

      if (i>2) then 
        if (all(p%precv(i)%map%naggr == p%precv(i-1)%map%naggr)) then 
          newsz=i-1
        end if
        call psb_bcast(ictxt,newsz)
        if (newsz > 0) exit
      end if
    end do

    if (newsz > 0) then 
      if (me == 0) then 
        write(debug_unit,*) trim(name),&
             &': Warning: aggregates from level ',&
             & newsz
        write(debug_unit,*) trim(name),&
             &':                       to level ',&
             & iszv,' coincide.'
        write(debug_unit,*) trim(name),&
             &': Number of levels actually used :',newsz
        write(debug_unit,*)
      end if
      allocate(t_prec%precv(newsz),stat=info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_from_subroutine_,name,&
             & a_err='prec reallocation')
        goto 9999
      endif
      do i=1,newsz-1
        call mld_move_alloc(p%precv(i),t_prec%precv(i),info)
      end do
      call mld_move_alloc(p%precv(iszv),t_prec%precv(newsz),info)
      do i=newsz+1, iszv
        call p%precv(i)%free(info)
      end do
      call mld_move_alloc(t_prec,p,info) 
      ! Ignore errors from transfer
      info = psb_success_
      !
      ! Restart
      iszv = newsz
      ! Fix the pointers, but the level 1 should
      ! be already OK
      do i=2, iszv - 1 
        p%precv(i)%base_a       => p%precv(i)%ac
        p%precv(i)%base_desc    => p%precv(i)%desc_ac
        p%precv(i)%map%p_desc_X => p%precv(i-1)%base_desc
        p%precv(i)%map%p_desc_Y => p%precv(i)%base_desc
      end do

      i    = iszv 
      if (info == psb_success_) call mld_coarse_bld(p%precv(i-1)%base_a,&
           & p%precv(i-1)%base_desc, p%precv(i),info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_from_subroutine_,name,a_err='coarse rebuild')
        goto 9999
      endif
    end if
  end if

  !
  ! The coarse space hierarchy has been build. 
  !
  ! Now do the preconditioner build.
  !

  do i=1, iszv
    !
    ! build the base preconditioner at level i
    !
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Calling mlprcbld at level  ',i
    call mld_check_def(p%precv(i)%parms%sweeps,&
         & 'Jacobi sweeps',ione,is_legal_jac_sweeps)
    call mld_check_def(p%precv(i)%parms%sweeps_pre,&
         & 'Jacobi sweeps',ione,is_legal_jac_sweeps)
    call mld_check_def(p%precv(i)%parms%sweeps_post,&
         & 'Jacobi sweeps',ione,is_legal_jac_sweeps)

    if (.not.allocated(p%precv(i)%sm)) then 
      !! Error: should have called mld_dprecinit
      info=3111
      call psb_errpush(info,name)
      goto 9999
    end if
    if (.not.allocated(p%precv(i)%sm%sv)) then 
      !! Error: should have called mld_dprecinit
      info=3111
      call psb_errpush(info,name)
      goto 9999
    end if

    call p%precv(i)%sm%build(p%precv(i)%base_a,p%precv(i)%base_desc,&
         & 'F',info,amold=amold,vmold=vmold)

    if ((info == psb_success_).and.(i>1).and.(present(amold))) then 
      call psb_map_cscnv(p%precv(i)%map,info,mold=amold)
      call p%precv(i)%ac%cscnv(info,mold=amold)
    end if
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='One level preconditioner build.')
      goto 9999
    endif

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Return from ',i,' call to mlprcbld ',info      
  end do


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

  subroutine doallc(i,v,src,info)
    type(mld_c_onelev_type), intent(inout)      :: v(:)
    integer(psb_ipk_), intent(in)                 :: i
    class(mld_c_base_smoother_type), intent(in) :: src
    integer(psb_ipk_), intent(out)                :: info
    
    if ((lbound(v,1)<=i).and.(i<=ubound(v,1))) then 
      allocate(v(i)%sm,source=src,stat=info)
    end if
  end subroutine doallc

end subroutine mld_cmlprec_bld
