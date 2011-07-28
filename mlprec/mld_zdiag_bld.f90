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
! File: mld_zdiag_bld.f90
!
! Subroutine: mld_zdiag_bld
! Version:    complex
!
!  This routine builds the diagonal preconditioner corresponding to a given
!  sparse matrix A.    
!
!
! Arguments:
!    a       -  type(psb_zspmat_type), input.
!               The sparse matrix structure containing the local part of the
!               matrix A to be preconditioned.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor associated to the sparse matrix A.
!    p       -  type(mld_zbaseprec_type), input/output.
!               The 'base preconditioner' data structure containing the local 
!               part of the diagonal preconditioner.
!    info    -  integer, output.
!               Error code.
!  
subroutine mld_zdiag_bld(a,desc_a,p,info)

  use psb_base_mod
  use mld_z_inner_mod, mld_protect_name => mld_zdiag_bld

  Implicit None

! Arguments
  type(psb_zspmat_type),intent(in), target :: a
  type(psb_desc_type), intent(in)          :: desc_a
  type(mld_zbaseprec_type),intent(inout)   :: p
  integer, intent(out)                     :: info

! Local variables
  Integer           :: err_act,ictxt, me, np, n_row, n_col,i
  integer           :: debug_level, debug_unit
  character(len=20) :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info  = psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  name  = 'mld_zdiag_bld'
  info  = psb_success_
  ictxt = desc_a%get_context()
  n_row = desc_a%get_local_rows()
  n_col = desc_a%get_local_cols()
  call psb_info(ictxt, me, np)

  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),' Enter'

  call psb_realloc(n_col,p%d,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_realloc')
    goto 9999
  end if

  !
  ! Retrieve the diagonal entries of the matrix A
  !
  call psb_sp_getdiag(a,p%d,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_sp_getdiag'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  !
  ! Copy into p%desc_data the descriptor associated to A
  !
  call psb_cdcpy(desc_a,p%desc_Data,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_cdcpy')
    goto 9999
  end if

  !
  ! The i-th diagonal entry of the preconditioner is set to one if the
  ! corresponding entry a_ii of the sparse matrix A is zero; otherwise 
  ! it is set to one/a_ii
  !
  do i=1,n_row
    if (p%d(i) == zzero) then
      p%d(i) = zone
    else
      p%d(i) = zone/p%d(i)
    endif
  end do

  if (a%pl(1) /= 0) then
    !
    ! Apply the same row permutation as in the sparse matrix A
    !
    call  psb_gelp('n',a%pl,p%d,info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_gelp'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  endif

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),'Done'

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_zdiag_bld

