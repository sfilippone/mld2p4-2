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
! File: mld_c_hierarchy_bld.f90
!
! Subroutine: mld_c_hierarchy_bld
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
subroutine mld_c_hierarchy_bld(a,desc_a,p,info,amold,vmold,imold)

  use psb_base_mod
  use mld_c_inner_mod
  use mld_c_prec_mod, mld_protect_name => mld_c_hierarchy_bld

  Implicit None

  ! Arguments
  type(psb_cspmat_type),intent(in), target           :: a
  type(psb_desc_type), intent(inout), target           :: desc_a
  type(mld_cprec_type),intent(inout),target          :: p
  integer(psb_ipk_), intent(out)                       :: info
  class(psb_c_base_sparse_mat), intent(in), optional :: amold
  class(psb_c_base_vect_type), intent(in), optional  :: vmold
  class(psb_i_base_vect_type), intent(in), optional  :: imold
!!$  character, intent(in), optional         :: upd

  ! Local Variables
  integer(psb_ipk_)  :: ictxt, me,np
  integer(psb_ipk_)  :: err,i,k, err_act, iszv, newsz, casize, nplevs, mxplevs, iaggsize
  real(psb_spk_)     :: mnaggratio, sizeratio
  class(mld_c_base_smoother_type), allocatable :: coarse_sm, base_sm, med_sm
  type(mld_sml_parms)              :: baseparms, medparms, coarseparms
  integer(psb_ipk_), allocatable   :: ilaggr(:), nlaggr(:)
  type(psb_cspmat_type)            :: op_prol
  type(mld_c_onelev_type), allocatable :: tprecv(:)    
  integer(psb_ipk_)  :: int_err(5)
  character          :: upd_
  integer(psb_ipk_)  :: debug_level, debug_unit
  character(len=20)  :: name, ch_err

  if (psb_get_errstatus().ne.0) return 
  info=psb_success_
  err=0
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  name = 'mld_c_hierarchy_bld'
  info = psb_success_
  int_err(1) = 0
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  p%ictxt = ictxt
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Entering '
  !
  ! For the time being we are commenting out the UPDATE argument
  ! we plan to resurrect it later. 
  ! !$  if (present(upd)) then 
  ! !$    if (debug_level >= psb_debug_outer_) &
  ! !$         & write(debug_unit,*) me,' ',trim(name),'UPD ', upd
  ! !$
  ! !$    if ((psb_toupper(upd).eq.'F').or.(psb_toupper(upd).eq.'T')) then
  ! !$      upd_=psb_toupper(upd)
  ! !$    else
  ! !$      upd_='F'
  ! !$    endif
  ! !$  else
  ! !$    upd_='F'
  ! !$  endif
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
  newsz      = -1
  casize     = p%coarse_aggr_size
  nplevs     = p%n_prec_levs
  mxplevs    = p%max_prec_levs
  mnaggratio = p%min_aggr_ratio
  casize     = p%coarse_aggr_size
  iszv       = size(p%precv)
  call psb_bcast(ictxt,iszv)
  call psb_bcast(ictxt,casize)
  call psb_bcast(ictxt,nplevs)
  call psb_bcast(ictxt,mxplevs)
  call psb_bcast(ictxt,mnaggratio)
  if (casize /= p%coarse_aggr_size) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent coarse_aggr_size')
    goto 9999
  end if
  if (nplevs /= p%n_prec_levs) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent n_prec_levs')
    goto 9999
  end if
  if (mxplevs /= p%max_prec_levs) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent max_prec_levs')
    goto 9999
  end if
  if (mnaggratio /= p%min_aggr_ratio) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent min_aggr_ratio')
    goto 9999
  end if
  if (iszv /= size(p%precv)) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent size of precv')
    goto 9999
  end if

  if (iszv <= 1) then
    ! We should only ever get here for multilevel.
    info=psb_err_from_subroutine_
    ch_err='size bpv'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  endif

  !
  ! The strategy:
  ! 1. The maximum number of levels should be already encoded in the
  !    size of the array;
  ! 2. If the user did not specify anything, then a default coarse size
  !    is generated, and the number of levels is set to the maximum;
  ! 3. If the number of levels has been specified, make sure it's capped
  !    at the maximum;
  ! 4. If the size of the array is different from target number of levels,
  !    reallocate;
  ! 5. Build the matrix hierarchy, stopping early if either the target
  !    coarse size is hit, or the gain fall below the mmin_aggr_ratio
  !    threshold.
  !


  if (nplevs <= 0) then
    if (casize <=0) then
      !
      ! Default to the cubic root of the size at base level.
      ! 
      casize = desc_a%get_global_rows()
      casize = int((sone*casize)**(sone/(sone*3)),psb_ipk_)
      casize = max(casize,ione)
      casize = casize*40_psb_ipk_ 
    end if
    nplevs = mxplevs    
  end if

  nplevs = max(itwo,min(nplevs,mxplevs))

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
  ! First set desired number of levels
  !
  if (iszv /= nplevs) then
    allocate(tprecv(nplevs),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,&
           & a_err='prec reallocation')
      goto 9999
    endif
    tprecv(1)%parms = baseparms
    allocate(tprecv(1)%sm,source=base_sm,stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,&
           & a_err='prec reallocation')
      goto 9999
    endif
    do i=2,nplevs-1
      tprecv(i)%parms = medparms
      allocate(tprecv(i)%sm,source=med_sm,stat=info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_from_subroutine_,name,&
             & a_err='prec reallocation')
        goto 9999
      endif
    end do
    tprecv(nplevs)%parms = coarseparms
    allocate(tprecv(nplevs)%sm,source=coarse_sm,stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,&
           & a_err='prec reallocation')
      goto 9999
    endif
    do i=1,iszv
      call p%precv(i)%free(info)
    end do
    call move_alloc(tprecv,p%precv)
    iszv = size(p%precv)
  end if

  !
  ! Finest level first; remember to fix base_a and base_desc
  ! 
  p%precv(1)%base_a    => a
  p%precv(1)%base_desc => desc_a
  newsz = 0
  array_build_loop: do i=2, iszv
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
    if (info == psb_success_) call mld_aggrmap_bld(p%precv(i),&
         & p%precv(i-1)%base_a,p%precv(i-1)%base_desc,&
         & ilaggr,nlaggr,op_prol,info)

    if (info /= psb_success_) then 
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Map build')
      goto 9999
    endif

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Return from ',i,' call to mlprcbld ',info      


    !
    ! Check for early termination of aggregation loop. 
    !      
    iaggsize = sum(nlaggr)
    if (iaggsize <= casize) then
      newsz = i
    end if

    if (i>2) then
      sizeratio = iaggsize
      sizeratio = sum(p%precv(i-1)%map%naggr)/sizeratio
      if (sizeratio < mnaggratio) then
        if (sizeratio > 1) then
          newsz = i
        else
          !
          ! We are not gaining 
          !
          newsz = newsz-1
        end if
      end if

      if (all(nlaggr == p%precv(i-1)%map%naggr)) then 
        newsz=i-1
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
      end if
    end if
    call psb_bcast(ictxt,newsz)
    if (newsz > 0) &
         &  call p%precv(i)%parms%get_coarse(p%precv(iszv)%parms)

    if (info == psb_success_) call mld_lev_mat_asb(p%precv(i),&
         & p%precv(i-1)%base_a,p%precv(i-1)%base_desc,&
         & ilaggr,nlaggr,op_prol,info)

    if (info /= psb_success_) then 
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Map build')
      goto 9999
    endif

    if (newsz > 0) exit array_build_loop
  end do array_build_loop

  if (newsz > 0) then
    !
    ! We exited early from the build loop, need to fix
    ! the size.
    !
    allocate(tprecv(newsz),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,&
           & a_err='prec reallocation')
      goto 9999
    endif
    do i=1,newsz
      call p%precv(i)%move_alloc(tprecv(i),info)
    end do
    do i=newsz+1, iszv
      call p%precv(i)%free(info)
    end do
    call move_alloc(tprecv,p%precv) 
    ! Ignore errors from transfer
    info = psb_success_
    !
    ! Restart
    iszv = newsz
    ! Fix the pointers, but the level 1 should
    ! be already OK
    do i=2, iszv 
      p%precv(i)%base_a       => p%precv(i)%ac
      p%precv(i)%base_desc    => p%precv(i)%desc_ac
      p%precv(i)%map%p_desc_X => p%precv(i-1)%base_desc
      p%precv(i)%map%p_desc_Y => p%precv(i)%base_desc
    end do
  end if


  if (info /= psb_success_) then 
    call psb_errpush(psb_err_internal_error_,name,&
         & a_err='Internal hierarchy build' )
    goto 9999
  endif

  iszv = size(p%precv)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Exiting with',iszv,' levels'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine mld_c_hierarchy_bld
