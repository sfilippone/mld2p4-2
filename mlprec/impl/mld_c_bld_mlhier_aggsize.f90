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
!
! Build an aggregation hierarchy with a target aggregation size
!
!
subroutine mld_c_bld_mlhier_aggsize(casize,mxplevs,mnaggratio,a,desc_a,precv,info)
  use psb_base_mod
  use mld_c_inner_mod, mld_protect_name => mld_c_bld_mlhier_aggsize
  use mld_c_prec_mod
  implicit none 
  integer(psb_ipk_), intent(in)    :: casize,mxplevs
  real(psb_spk_)                   :: mnaggratio
  type(psb_cspmat_type),intent(in), target           :: a
  type(psb_desc_type), intent(inout), target         :: desc_a
  type(mld_c_onelev_type), allocatable, target, intent(inout)  :: precv(:)
  integer(psb_ipk_), intent(out)   :: info
  ! Local
  integer(psb_ipk_) :: ictxt, me,np
  integer(psb_ipk_) :: err,i,k, err_act, iszv, newsz, iaggsize
  integer(psb_ipk_) :: ipv(mld_ifpsz_), val
  integer(psb_ipk_) :: int_err(5)
  character         :: upd_
  class(mld_c_base_smoother_type), allocatable :: coarse_sm, base_sm, med_sm
  type(mld_sml_parms)              :: baseparms, medparms, coarseparms
  type(mld_c_onelev_node), pointer :: head, tail, newnode, current
  real(psb_spk_)     :: sizeratio
  integer(psb_ipk_)  :: debug_level, debug_unit
  character(len=20)  :: name, ch_err
  name = 'mld_bld_mlhier_aggsize'
  if (psb_get_errstatus().ne.0) return 
  info=psb_success_
  err=0
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ictxt  = desc_a%get_ctxt()
  call psb_info(ictxt,me,np)
  ! 
  ! New strategy to build according to coarse size. 
  !
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

    current => current%next
    tail    => current
    iaggsize = sum(current%item%map%naggr)

    if (iaggsize <= casize) then 
      !
      ! Target reached; but we may need to rebuild. 
      !
      exit list_build_loop
    end if
    if (newsz>2) then
      sizeratio = iaggsize
      sizeratio = sum(current%prev%item%map%naggr)/sizeratio

      if (sizeratio < mnaggratio) then 
        !
        ! We are not gaining 
        !
        newsz = newsz-1
        current => current%prev
        current%next =>null()
        call newnode%item%free(info)
        if (info == psb_success_) deallocate(newnode,stat=info)
        if (info /= psb_success_) then 
          info = psb_err_internal_error_
          call psb_errpush(info,name,a_err='Deallocate at list end'); goto 9999
        end if
        exit list_build_loop          
      end if
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
    if (info == psb_success_) call precv(i)%free(info)
  end do
  if (info == psb_success_) deallocate(precv,stat=info)
  if (info == psb_success_) allocate(precv(newsz),stat=info) 
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
         & call current%item%move_alloc(precv(i),info)
    ! Now set the smoother/solver parts. 
    if (info == psb_success_) then 
      if (i ==1) then 
        ! This is a workaround for a bug in gfortran 4.7.2
        allocate(precv(i)%sm,source=base_sm,stat=info) 
      else if (i < newsz) then 
        allocate(precv(i)%sm,source=med_sm,stat=info) 
      else
        allocate(precv(i)%sm,source=coarse_sm,stat=info) 
      end if
    end if
    if (info /= psb_success_) then 
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='list cpy'); goto 9999
    end if
    if (i == 1) then 
      precv(i)%base_a       => a
      precv(i)%base_desc    => desc_a
    else
      precv(i)%base_a       => precv(i)%ac
      precv(i)%base_desc    => precv(i)%desc_ac
      precv(i)%map%p_desc_X => precv(i-1)%base_desc
      precv(i)%map%p_desc_Y => precv(i)%base_desc
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

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine mld_c_bld_mlhier_aggsize
