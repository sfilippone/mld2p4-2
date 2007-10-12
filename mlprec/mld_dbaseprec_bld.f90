!!$ 
!!$ 
!!$                    MD2P4
!!$    Multilevel Domain Decomposition Parallel Preconditioner Package for PSBLAS
!!$                      for 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$                       Daniela di Serafino    Second University of Naples
!!$                       Pasqua D'Ambra         ICAR-CNR                      
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the MD2P4 group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MD2P4 GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
subroutine mld_dbaseprc_bld(a,desc_a,p,info,upd)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_dbaseprc_bld

  Implicit None

  type(psb_dspmat_type), target           :: a
  type(psb_desc_type), intent(in), target :: desc_a
  type(mld_dbaseprc_type),intent(inout)   :: p
  integer, intent(out)                    :: info
  character, intent(in), optional         :: upd

  ! Local scalars
  Integer      :: err, nnzero, n_row, n_col,I,j,k,ictxt,&
       & me,mycol,np,npcol,mglob,lw, mtype, nrg, nzg, err_act
  real(kind(1.d0))         :: temp, real_err(5)
  real(kind(1.d0)),pointer :: gd(:), work(:)
  integer      :: int_err(5)
  character    :: iupd

  logical, parameter :: debug=.false.  
  integer,parameter  :: iroot=0,iout=60,ilout=40
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  err=0
  call psb_erractionsave(err_act)
  name = 'mld_dbaseprc_bld'

  if (debug) write(0,*) 'Entering baseprc_bld'
  info = 0
  int_err(1) = 0
  ictxt   = psb_cd_get_context(desc_a)
  n_row   = psb_cd_get_local_rows(desc_a)
  n_col   = psb_cd_get_local_cols(desc_a)
  mglob   = psb_cd_get_global_rows(desc_a)

  if (debug) write(0,*) 'Preconditioner Blacs_gridinfo'
  call psb_info(ictxt, me, np)

  if (present(upd)) then 
    if (debug) write(0,*) 'UPD ', upd
    if ((UPD.eq.'F').or.(UPD.eq.'T')) then
      IUPD=UPD
    else
      IUPD='F'
    endif
  else
    IUPD='F'
  endif

  !
  ! Should add check to ensure all procs have the same... 
  !

  call mld_check_def(p%iprcparm(mld_prec_type_),'base_prec',&
       &  mld_diag_,is_legal_base_prec)


  call psb_nullify_desc(p%desc_data)

  select case(p%iprcparm(mld_prec_type_)) 
  case (mld_noprec_)
    ! Do nothing. 
    call psb_cdcpy(desc_a,p%desc_data,info)
    if(info /= 0) then
      info=4010
      ch_err='psb_cdcpy'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  case (mld_diag_)

    call mld_diag_bld(a,desc_a,p,iupd,info)
    if(debug) write(0,*)me,': out of mld_diag_bld'
    if(info /= 0) then
      info=4010
      ch_err='mld_diag_bld'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  case(mld_bjac_,mld_as_)

    call mld_check_def(p%iprcparm(mld_n_ovr_),'overlap',&
         &  0,is_legal_n_ovr)
    call mld_check_def(p%iprcparm(mld_sub_restr_),'restriction',&
         &  psb_halo_,is_legal_restrict)
    call mld_check_def(p%iprcparm(mld_sub_prol_),'prolongator',&
         &  psb_none_,is_legal_prolong)
    call mld_check_def(p%iprcparm(mld_sub_ren_),'renumbering',&
         &  mld_renum_none_,is_legal_renum)
    call mld_check_def(p%iprcparm(mld_sub_solve_),'fact',&
         &  mld_ilu_n_,is_legal_ml_fact)
    
    if (p%iprcparm(mld_sub_solve_)==mld_sludist_) then
      p%iprcparm(mld_n_ovr_)      = 0
      p%iprcparm(mld_smooth_sweeps_) = 1
    end if

    if (debug) write(0,*)me, ': Calling mld_bjac_bld'
    if (debug) call psb_barrier(ictxt)

    call mld_bjac_bld(a,desc_a,p,iupd,info)
    if(info /= 0) then
      info=4010
      call psb_errpush(info,name,a_err='mld_bjac_bld')
      goto 9999
    end if

  case default
    info=4010
    ch_err='Unknown mld_prec_type_'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999

  end select

  p%base_a    => a
  p%base_desc => desc_a
  p%iprcparm(mld_prec_status_) = mld_prec_built_
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_dbaseprc_bld

