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

subroutine mld_s_bld_mlhier_array(nplevs,a,desc_a,precv,info)
  use psb_base_mod
  use mld_s_inner_mod, mld_protect_name => mld_s_bld_mlhier_array
  use mld_s_prec_mod
  implicit none 
  integer(psb_ipk_), intent(inout) :: nplevs
  type(psb_sspmat_type),intent(in), target           :: a
  type(psb_desc_type), intent(inout), target           :: desc_a
  type(mld_s_onelev_type),intent(inout), allocatable, target :: precv(:)
  integer(psb_ipk_), intent(out)   :: info
  ! Local
  integer(psb_ipk_)      :: ictxt, me,np
  integer(psb_ipk_)      :: err,i,k, err_act, newsz, iszv
  integer(psb_ipk_)      :: ipv(mld_ifpsz_), val
  class(mld_s_base_smoother_type), allocatable :: coarse_sm, base_sm, med_sm
  type(mld_sml_parms)              :: baseparms, medparms, coarseparms
  type(mld_s_onelev_type), allocatable :: tprecv(:)    
  integer(psb_ipk_)      :: int_err(5)
  integer(psb_ipk_)      :: debug_level, debug_unit
  character(len=20)  :: name, ch_err
  name = 'mld_bld_array_hierarchy'
  if (psb_get_errstatus().ne.0) return 
  info=psb_success_
  err=0
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ictxt  = desc_a%get_ctxt()
  call psb_info(ictxt,me,np)
  iszv        = size(precv)

  coarseparms = precv(iszv)%parms
  baseparms   = precv(1)%parms
  medparms    = precv(2)%parms

  allocate(coarse_sm, source=precv(iszv)%sm,stat=info) 
  if (info == psb_success_)  &
       & allocate(med_sm, source=precv(2)%sm,stat=info) 
  if (info == psb_success_)  &
       & allocate(base_sm, source=precv(1)%sm,stat=info) 
  if (info /= psb_success_) then 
    write(0,*) 'Error in saving smoothers',info
    call psb_errpush(psb_err_internal_error_,name,a_err='Base level precbuild.')
    goto 9999
  end if

  !
  !
  ! Build the matrix and the transfer operators corresponding
  ! to the remaining levels
  !
  !
  ! Check on the iprcparm contents: they should be the same
  ! on all processes.
  !
  call psb_bcast(ictxt,precv(1)%parms)
  iszv = size(precv)
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
      call precv(i)%free(info)
    end do
    call move_alloc(tprecv,precv)
    iszv = size(precv)
  end if
    
  !
  ! Finest level first; remember to fix base_a and base_desc
  ! 
  precv(1)%base_a    => a
  precv(1)%base_desc => desc_a

  array_build_loop: do i=2, iszv
    !
    ! Check on the iprcparm contents: they should be the same
    ! on all processes.
    !
    call psb_bcast(ictxt,precv(i)%parms)

    !
    ! Sanity checks on the parameters
    !
    if (i<iszv) then 
      !
      ! A replicated matrix only makes sense at the coarsest level
      !
      call mld_check_def(precv(i)%parms%coarse_mat,'Coarse matrix',&
           &   mld_distr_mat_,is_distr_ml_coarse_mat)
    end if

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Calling mlprcbld at level  ',i
    !
    ! Build the mapping between levels i-1 and i and the matrix
    ! at level i
    ! 
    if (info == psb_success_) call mld_coarse_bld(precv(i-1)%base_a,&
         & precv(i-1)%base_desc, precv(i),info)

    if (info /= psb_success_) then 
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Init upper level preconditioner')
      goto 9999
    endif

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Return from ',i,' call to mlprcbld ',info      

    if (i>2) then 
      if (all(precv(i)%map%naggr == precv(i-1)%map%naggr)) then 
        newsz=i-1
      end if
      call psb_bcast(ictxt,newsz)
      if (newsz > 0) exit array_build_loop
    end if
  end do array_build_loop

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
    allocate(tprecv(newsz),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,&
           & a_err='prec reallocation')
      goto 9999
    endif
    do i=1,newsz-1
      call precv(i)%move_alloc(tprecv(i),info)
    end do
    call precv(iszv)%move_alloc(tprecv(newsz),info)
    do i=newsz+1, iszv
      call precv(i)%free(info)
    end do
    call move_alloc(tprecv,precv) 
    ! Ignore errors from transfer
    info = psb_success_
    !
    ! Restart
    iszv = newsz
    ! Fix the pointers, but the level 1 should
    ! be already OK
    do i=2, iszv - 1 
      precv(i)%base_a       => precv(i)%ac
      precv(i)%base_desc    => precv(i)%desc_ac
      precv(i)%map%p_desc_X => precv(i-1)%base_desc
      precv(i)%map%p_desc_Y => precv(i)%base_desc
    end do

    i    = iszv 
    if (info == psb_success_) call mld_coarse_bld(precv(i-1)%base_a,&
         & precv(i-1)%base_desc, precv(i),info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='coarse rebuild')
      goto 9999
    endif
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine mld_s_bld_mlhier_array
