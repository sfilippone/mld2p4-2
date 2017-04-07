!  
!   
!                             MLD2P4  version 2.1
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.4)
!    
!    (C) Copyright 2008, 2010, 2012, 2015, 2017 
!  
!        Salvatore Filippone    Cranfield University
!        Ambra Abdullahi Hassan University of Rome Tor Vergata
!        Alfredo Buttari        CNRS-IRIT, Toulouse
!        Pasqua D'Ambra         ICAR-CNR, Naples
!        Daniela di Serafino    University of Campania "L. Vanvitelli", Caserta
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the MLD2P4 group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!  
subroutine mld_d_base_onelev_build(lv,info,amold,vmold,imold)
  use psb_base_mod
  use mld_d_onelev_mod, mld_protect_name => mld_d_base_onelev_build
  implicit none
  class(mld_d_onelev_type), target, intent(inout) :: lv
  integer(psb_ipk_), intent(out) :: info
  class(psb_d_base_sparse_mat), intent(in), optional :: amold
  class(psb_d_base_vect_type), intent(in), optional  :: vmold
  class(psb_i_base_vect_type), intent(in), optional  :: imold
  ! Local
  integer(psb_ipk_)  :: err,i,k, err_act, iszv, newsz, casize
  integer(psb_ipk_)  :: ictxt, me, np
  integer(psb_ipk_)  :: ipv(mld_ifpsz_), val
  integer(psb_ipk_)  :: int_err(5)
  integer(psb_ipk_)  :: debug_level, debug_unit
  character(len=20)  :: name, ch_err

  name = 'mld_onelev_build'
  if (psb_get_errstatus().ne.0) return 
  info=psb_success_
  err=0
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  if (.not.associated(lv%base_desc)) then
    info = psb_err_internal_error_
    call psb_errpush(info,name,&
         & a_err='Unassociated base DESC')
    goto 9999
  end if
  info = psb_success_
  int_err(1) = 0
  ictxt = lv%base_desc%get_ctxt()
  call psb_info(ictxt,me,np)

  if (.not.allocated(lv%sm)) then 
    !! Error: should have called mld_dprecinit
    info=3111
    call psb_errpush(info,name)
    goto 9999
  end if
  if (.not.allocated(lv%sm%sv)) then 
    !! Error: should have called mld_dprecinit
    info=3111
    call psb_errpush(info,name)
    goto 9999
  end if
  lv%ac_nz_loc = lv%ac%get_nzeros()
  lv%ac_nz_tot = lv%ac_nz_loc
  call psb_sum(ictxt,lv%ac_nz_tot)
  !  end do
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Calling mlprcbld at level  ',i
  call mld_check_def(lv%parms%sweeps,&
       & 'Jacobi sweeps',izero,is_int_non_negative)
  call mld_check_def(lv%parms%sweeps_pre,&
       & 'Jacobi sweeps',izero,is_int_non_negative)
  call mld_check_def(lv%parms%sweeps_post,&
       & 'Jacobi sweeps',izero,is_int_non_negative)

  call lv%sm%build(lv%base_a,lv%base_desc,&
       & 'F',info,amold=amold,vmold=vmold,imold=imold)
  if (info == 0) then
    if (allocated(lv%sm2a)) then 
      call lv%sm2a%build(lv%base_a,lv%base_desc,'F',info,&
           & amold=amold,vmold=vmold,imold=imold)
      lv%sm2 => lv%sm2a
    else
      lv%sm2 => lv%sm
    end if
  end if
  if (info /=0) then
    info = psb_err_internal_error_
    call psb_errpush(info,name,&
         & a_err='Smoother bld error')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine mld_d_base_onelev_build
