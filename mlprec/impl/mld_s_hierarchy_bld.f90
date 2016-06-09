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
! File: mld_s_hierarchy_bld.f90
!
! Subroutine: mld_s_hierarchy_bld
! Version:    real
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
!    a       -  type(psb_sspmat_type).
!               The sparse matrix structure containing the local part of the
!               matrix to be preconditioned.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    p       -  type(mld_sprec_type), input/output.
!               The preconditioner data structure containing the local part
!               of the preconditioner to be built.
!    info    -  integer, output.
!               Error code.              
!
!    amold   -  class(psb_s_base_sparse_mat), input, optional
!               Mold for the inner format of matrices contained in the
!               preconditioner
!
!
!    vmold   -  class(psb_s_base_vect_type), input, optional
!               Mold for the inner format of vectors contained in the
!               preconditioner
!
!
!  
subroutine mld_s_hierarchy_bld(a,desc_a,p,info,amold,vmold,imold)

  use psb_base_mod
  use mld_s_inner_mod
  use mld_s_prec_mod, mld_protect_name => mld_s_hierarchy_bld

  Implicit None

  ! Arguments
  type(psb_sspmat_type),intent(in), target           :: a
  type(psb_desc_type), intent(inout), target           :: desc_a
  type(mld_sprec_type),intent(inout),target          :: p
  integer(psb_ipk_), intent(out)                       :: info
  class(psb_s_base_sparse_mat), intent(in), optional :: amold
  class(psb_s_base_vect_type), intent(in), optional  :: vmold
  class(psb_i_base_vect_type), intent(in), optional  :: imold
!!$  character, intent(in), optional         :: upd

  ! Local Variables
  integer(psb_ipk_)  :: ictxt, me,np
  integer(psb_ipk_)  :: err,i,k, err_act, iszv, newsz, casize, nplevs, mxplevs
  real(psb_spk_)     :: mnaggratio
  integer(psb_ipk_)  :: ipv(mld_ifpsz_), val
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

  name = 'mld_s_hierarchy_bld'
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
    !! Error: should have called mld_sprecinit
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

  if (nplevs <= 0) then
    !
    ! This should become the default strategy, we specify a target aggregation size.
    !
    if (casize <=0) then
      !
      ! Default to the cubic root of the size at base level.
      ! 
      casize = desc_a%get_global_rows()
      casize = int((sone*casize)**(sone/(sone*3)),psb_ipk_)
      casize = max(casize,ione)
    end if
    call mld_bld_mlhier_aggsize(casize,mxplevs,mnaggratio,a,desc_a,p%precv,info)
  else 
    ! 
    ! Oldstyle with fixed number of levels. 
    !
    nplevs = max(itwo,min(nplevs,mxplevs))
    call mld_bld_mlhier_array(nplevs,a,desc_a,p%precv,info)
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

end subroutine mld_s_hierarchy_bld
