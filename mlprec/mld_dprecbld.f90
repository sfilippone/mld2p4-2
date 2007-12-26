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
! File: mld_dprecbld.f90
!
! Subroutine: mld_dprecbld
! Version:    real
! Contains:   subroutine init_baseprc_av
!
!  This routine builds the preconditioner according to the requirements made by
!  the user trough the subroutines mld_precinit and mld_precset.
!  
!  A multilevel preconditioner is regarded as an array of 'base preconditioners',
!  each representing the part of the preconditioner associated to a certain level.
!  The levels are numbered in increasing order starting from the finest      one, i.e.
!  level 1 is the finest level. 
! 
!
! Arguments:
!    a       -  type(psb_dspmat_type).
!               The sparse matrix structure containing the local part of the
!               matrix to be preconditioned.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    p       -  type(mld_dprec_type), input/output.
!               The preconditioner data structure containing the local part
!               of the preconditioner to be built.
!    info    -  integer, output.
!               Error code.              
!  
subroutine mld_dprecbld(a,desc_a,p,info,upd)

  use psb_base_mod
  use mld_prec_mod, protect => mld_dprecbld
  
  Implicit None

  ! Arguments
  type(psb_dspmat_type), target           :: a
  type(psb_desc_type), intent(in), target :: desc_a
  type(mld_dprec_type),intent(inout)      :: p
  integer, intent(out)                    :: info
  character, intent(in), optional         :: upd

  ! Local Variables
  Integer      :: err,i,k,ictxt, me,np, err_act, iszv
  integer      :: int_err(5)
  character    :: iupd
  integer            :: debug_level, debug_unit
  character(len=20)  :: name, ch_err

  if (psb_get_errstatus().ne.0) return 
  info=0
  err=0
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  name = 'mld_dprecbld'
  info = 0
  int_err(1) = 0
  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Entering ',desc_a%matrix_data(:)

  if (present(upd)) then 
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),'UPD ', upd

    if ((toupper(upd).eq.'F').or.(toupper(upd).eq.'T')) then
      iupd=toupper(upd)
    else
      iupd='F'
    endif
  else
    iupd='F'
  endif

  if (.not.allocated(p%baseprecv)) then 
    !! Error: should have called mld_dprecinit
    info=3111
    call psb_errpush(info,name)
    goto 9999
  end if

  !
  ! Check to ensure all procs have the same 
  !   
  iszv = size(p%baseprecv)
  call psb_bcast(ictxt,iszv)
  if (iszv /= size(p%baseprecv)) then 
    info=4001
    call psb_errpush(info,name,a_err='Inconsistent size of baseprecv')
    goto 9999
  end if

  if (iszv >= 1) then
    !
    ! Allocate and build the fine level preconditioner
    ! 
    call init_baseprc_av(p%baseprecv(1),info)
    if (info == 0) call mld_baseprc_bld(a,desc_a,p%baseprecv(1),info,iupd)

    if (info /= 0) then 
      call psb_errpush(4001,name,a_err='Base level precbuild.')
      goto 9999
    end if

  else
    info=4010
    ch_err='size bpv'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  endif

  if (iszv > 1) then

    !
    ! Build the base preconditioners corresponding to the remaining
    ! levels
    !
    do i=2, iszv
      
      !
      ! Allocate the av component of the preconditioner data type
      ! at level i
      !
      if (i<iszv) then 
        !
        ! A replicated matrix only makes sense at the coarsest level
        !
        call mld_check_def(p%baseprecv(i)%iprcparm(mld_coarse_mat_),'Coarse matrix',&
             &   mld_distr_mat_,is_distr_ml_coarse_mat)
      end if

      call init_baseprc_av(p%baseprecv(i),info)

      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & 'Calling mlprcbld at level  ',i

      !
      ! Build the base preconditioner corresponding to level i
      !
      if (info == 0) call mld_mlprec_bld(p%baseprecv(i-1)%base_a,&
           & p%baseprecv(i-1)%base_desc, p%baseprecv(i),info)
      if (info /= 0) then 
        call psb_errpush(4001,name,a_err='Init & build upper level preconditioner')
        goto 9999
      endif

      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & 'Return from ',i,' call to mlprcbld ',info      
    end do

  endif

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

  subroutine init_baseprc_av(p,info)
    type(mld_dbaseprc_type), intent(inout) :: p
    integer                                :: info
    if (allocated(p%av)) then
      if (size(p%av) /= mld_max_avsz_) then 
        deallocate(p%av,stat=info)
        if (info /= 0) return 
      endif
    end if
    if (.not.(allocated(p%av))) then 
      allocate(p%av(mld_max_avsz_),stat=info)
      if (info /= 0) return
    end if
    do k=1,size(p%av)
      call psb_nullify_sp(p%av(k))
    end do
    
  end subroutine init_baseprc_av

end subroutine mld_dprecbld

