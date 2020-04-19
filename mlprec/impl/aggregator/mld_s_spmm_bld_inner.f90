!   
!   
!                             MLD2P4  Extensions
!    
!    (C) Copyright 2019
!  
!                        Salvatore Filippone  Cranfield University
!        Pasqua D'Ambra         IAC-CNR, Naples, IT
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
! File: mld_daggrmat_nosmth_bld.F90
!
!
subroutine mld_s_spmm_bld_inner(a_csr,desc_a,nlaggr,parms,ac,&
     & coo_prol,desc_cprol,coo_restr,info)
  use psb_base_mod
  use mld_s_inner_mod
  use mld_s_base_aggregator_mod, mld_protect_name => mld_s_spmm_bld_inner
  implicit none

  ! Arguments
  type(psb_s_csr_sparse_mat), intent(inout) :: a_csr
  type(psb_desc_type), intent(in)            :: desc_a
  integer(psb_lpk_), intent(inout)           :: nlaggr(:)
  type(mld_sml_parms), intent(inout)         :: parms 
  type(psb_ls_coo_sparse_mat), intent(inout)  :: coo_prol, coo_restr
  type(psb_desc_type), intent(inout)         :: desc_cprol
  type(psb_lsspmat_type), intent(out)        :: ac
  integer(psb_ipk_), intent(out)             :: info

  ! Local variables
  integer(psb_ipk_)  :: err_act
  integer(psb_ipk_)  :: ictxt,np,me, icomm, ndx, minfo
  character(len=40)  :: name
  integer(psb_ipk_)  :: ierr(5)
  type(psb_ls_coo_sparse_mat) :: ac_coo, tmpcoo
  type(psb_s_csr_sparse_mat) :: acsr3, csr_prol, ac_csr, csr_restr
  integer(psb_ipk_) :: debug_level, debug_unit, naggr
  integer(psb_lpk_) :: nrow, nglob, ncol, ntaggr, nzl, ip, &
       &  nzt, naggrm1, naggrp1, i, k
  integer(psb_lpk_) ::  nrsave, ncsave, nzsave, nza
  logical, parameter :: do_timings=.true., oldstyle=.false., debug=.false.  
  integer(psb_ipk_), save :: idx_spspmm=-1

  name='mld_spmm_bld_inner'
  if(psb_get_errstatus().ne.0) return 
  info=psb_success_
  call psb_erractionsave(err_act)


  ictxt = desc_a%get_context()
  icomm = desc_a%get_mpic()
  call psb_info(ictxt, me, np)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  nglob = desc_a%get_global_rows()
  nrow  = desc_a%get_local_rows()
  ncol  = desc_a%get_local_cols()

  if ((do_timings).and.(idx_spspmm==-1)) &
       & idx_spspmm = psb_get_timer_idx("SPMM_BLD: par_spspmm")

  naggr   = nlaggr(me+1)
  ntaggr  = sum(nlaggr)
  naggrm1 = sum(nlaggr(1:me))
  naggrp1 = sum(nlaggr(1:me+1)) 
  !write(0,*)me,' ',name,' input sizes',nlaggr(:),':',naggr

  !
  ! COO_PROL should arrive here with local numbering
  !
  if (debug) write(0,*)  me,' ',trim(name),' Size check on entry New: ',&
       & coo_prol%get_fmt(),coo_prol%get_nrows(),coo_prol%get_ncols(),coo_prol%get_nzeros(),&
       & nrow,ntaggr,naggr

  call coo_prol%cp_to_ifmt(csr_prol,info)

  if (debug) write(0,*) me,trim(name),' Product AxPROL ',&
       & a_csr%get_nrows(),a_csr%get_ncols(), csr_prol%get_nrows(), &
       & desc_a%get_local_rows(),desc_a%get_local_cols(),&
       & desc_cprol%get_local_rows(),desc_a%get_local_cols()
  if (debug) flush(0)

  if (do_timings) call psb_tic(idx_spspmm)
  call psb_par_spspmm(a_csr,desc_a,csr_prol,acsr3,desc_cprol,info)
  if (do_timings) call psb_toc(idx_spspmm)  

  if (debug) write(0,*) me,trim(name),' Done AxPROL ',&
       & acsr3%get_nrows(),acsr3%get_ncols(), acsr3%get_nzeros(),&
       & desc_cprol%get_local_rows(),desc_cprol%get_local_cols()

  !
  ! Ok first product done.
  !
  ! Remember that RESTR must be built from PROL after halo extension,
  ! which is done above in psb_par_spspmm
  if (debug) write(0,*)me,' ',name,' No inp_restr, transposing prol ',&
       & csr_prol%get_nrows(),csr_prol%get_ncols(),csr_prol%get_nzeros()
  call csr_prol%cp_to_lcoo(coo_restr,info)
!!$      write(0,*)me,' ',name,' new into transposition ',coo_restr%get_nrows(),&
!!$           & coo_restr%get_ncols(),coo_restr%get_nzeros()
  call coo_restr%transp()
  nzl = coo_restr%get_nzeros()
  call desc_cprol%l2gip(coo_restr%ia(1:nzl),info)
  i=0
  !
  ! Now we have to fix this.  The only rows of the restrictor that are correct 
  ! are those corresponding to "local" aggregates, i.e. indices in ilaggr(:)
  !
  do k=1, nzl
    if ((naggrm1 < coo_restr%ia(k)) .and.(coo_restr%ia(k) <= naggrp1)) then
      i = i+1
      coo_restr%val(i) = coo_restr%val(k)
      coo_restr%ia(i)  = coo_restr%ia(k)
      coo_restr%ja(i)  = coo_restr%ja(k)
    end if
  end do
  call coo_restr%set_nzeros(i)
  call coo_restr%fix(info)
  call coo_restr%cp_to_coo(tmpcoo,info)
!!$      write(0,*)me,' ',name,' after transposition ',coo_restr%get_nrows(),coo_restr%get_ncols(),coo_restr%get_nzeros()

  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='spcnv coo_restr')
    goto 9999
  end if
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'starting sphalo/ rwxtd'
  nzl    = tmpcoo%get_nzeros()    
  call psb_glob_to_loc(tmpcoo%ia(1:nzl),desc_cprol,info,iact='I',owned=.true.)
  call tmpcoo%clean_negidx(info)
  nzl  = tmpcoo%get_nzeros()
  call tmpcoo%set_nrows(desc_cprol%get_local_rows())
  call tmpcoo%set_ncols(desc_a%get_local_cols())
!!$    write(0,*)me,' ',name,' after G2L on rows ',tmpcoo%get_nrows(),tmpcoo%get_ncols(),tmpcoo%get_nzeros()      
  call csr_restr%mv_from_lcoo(tmpcoo,info)

  if (debug) write(0,*) me,trim(name),' Product RESTRxAP ',&
       & csr_restr%get_nrows(),csr_restr%get_ncols(), &
       & desc_cprol%get_local_rows(),desc_a%get_local_cols(),&
       & acsr3%get_nrows(),acsr3%get_ncols()
  if (do_timings) call psb_tic(idx_spspmm)      
  call psb_par_spspmm(csr_restr,desc_a,acsr3,ac_csr,desc_cprol,info)
  if (do_timings) call psb_toc(idx_spspmm)      
  call csr_restr%free()
  call acsr3%free()
  call ac_csr%mv_to_lcoo(ac_coo,info)
  call ac_coo%fix(info)
  nza    = ac_coo%get_nzeros()
  if (debug) write(0,*) me,trim(name),' Fixing ac ',&
       & ac_coo%get_nrows(),ac_coo%get_ncols(), nza
  call desc_cprol%indxmap%l2gip(ac_coo%ia(1:nza),info)
  call desc_cprol%indxmap%l2gip(ac_coo%ja(1:nza),info)
  call ac_coo%set_nrows(ntaggr)
  call ac_coo%set_ncols(ntaggr)
  if (debug) write(0,*)  me,' ',trim(name),' Before mv_from',psb_get_errstatus()
  if (info == 0) call ac%mv_from(ac_coo)
  if (debug) write(0,*)  me,' ',trim(name),' After  mv_from',psb_get_errstatus()
  if (debug) write(0,*)  me,' ',trim(name),' ',ac%get_fmt(),ac%get_nrows(),ac%get_ncols(),ac%get_nzeros(),naggr,ntaggr
  ! write(0,*)  me,' ',trim(name),' Final AC newstyle ',ac%get_fmt(),ac%get_nrows(),ac%get_ncols(),ac%get_nzeros()

  nza = coo_prol%get_nzeros()
  call desc_cprol%indxmap%l2gip(coo_prol%ja(1:nza),info)
  
  if (debug) then
    write(0,*) me,' ',trim(name),' Checkpoint at exit'
    call psb_barrier(ictxt)
    write(0,*) me,' ',trim(name),' Checkpoint through'
  end if

  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='Build ac = coo_restr x am3')
    goto 9999
  end if


  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done spmm_bld_inner '

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
  
end subroutine mld_s_spmm_bld_inner

subroutine mld_ls_spmm_bld_inner(a_csr,desc_a,nlaggr,parms,ac,&
     & coo_prol,desc_cprol,coo_restr,info)
  use psb_base_mod
  use mld_s_inner_mod
  use mld_s_base_aggregator_mod, mld_protect_name => mld_ls_spmm_bld_inner
  implicit none

  ! Arguments
  type(psb_ls_csr_sparse_mat), intent(inout) :: a_csr
  type(psb_desc_type), intent(in)            :: desc_a
  integer(psb_lpk_), intent(inout)           :: nlaggr(:)
  type(mld_sml_parms), intent(inout)         :: parms 
  type(psb_ls_coo_sparse_mat), intent(inout)  :: coo_prol, coo_restr
  type(psb_desc_type), intent(inout)         :: desc_cprol
  type(psb_lsspmat_type), intent(out)        :: ac
  integer(psb_ipk_), intent(out)             :: info

  ! Local variables
  integer(psb_ipk_)  :: err_act
  integer(psb_ipk_)  :: ictxt,np,me, icomm, ndx, minfo
  character(len=40)  :: name
  integer(psb_ipk_)  :: ierr(5)
  type(psb_ls_coo_sparse_mat) :: ac_coo, tmpcoo
  type(psb_ls_csr_sparse_mat) :: acsr3, csr_prol, ac_csr, csr_restr
  integer(psb_ipk_) :: debug_level, debug_unit, naggr
  integer(psb_lpk_) :: nrow, nglob, ncol, ntaggr, nzl, ip, &
       &  nzt, naggrm1, naggrp1, i, k
  integer(psb_lpk_) ::  nrsave, ncsave, nzsave, nza
  logical, parameter :: do_timings=.true., oldstyle=.false., debug=.false.  
  integer(psb_ipk_), save :: idx_spspmm=-1

  name='mld_spmm_bld_inner'
  if(psb_get_errstatus().ne.0) return 
  info=psb_success_
  call psb_erractionsave(err_act)


  ictxt = desc_a%get_context()
  icomm = desc_a%get_mpic()
  call psb_info(ictxt, me, np)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  nglob = desc_a%get_global_rows()
  nrow  = desc_a%get_local_rows()
  ncol  = desc_a%get_local_cols()

  if ((do_timings).and.(idx_spspmm==-1)) &
       & idx_spspmm = psb_get_timer_idx("SPMM_BLD: par_spspmm")

  naggr   = nlaggr(me+1)
  ntaggr  = sum(nlaggr)
  naggrm1 = sum(nlaggr(1:me))
  naggrp1 = sum(nlaggr(1:me+1)) 
  !write(0,*)me,' ',name,' input sizes',nlaggr(:),':',naggr

  !
  ! Here COO_PROL should be with GLOBAL indices on the cols
  ! and LOCAL indices on the rows. 
  !
  if (debug) write(0,*)  me,' ',trim(name),' Size check on entry New: ',&
       & coo_prol%get_fmt(),coo_prol%get_nrows(),coo_prol%get_ncols(),coo_prol%get_nzeros(),&
       & nrow,ntaggr,naggr

  call coo_prol%cp_to_fmt(csr_prol,info)

  if (debug) write(0,*) me,trim(name),' Product AxPROL ',&
       & a_csr%get_nrows(),a_csr%get_ncols(), csr_prol%get_nrows(), &
       & desc_a%get_local_rows(),desc_a%get_local_cols(),&
       & desc_cprol%get_local_rows(),desc_a%get_local_cols()
  if (debug) flush(0)

  if (do_timings) call psb_tic(idx_spspmm)
  call psb_par_spspmm(a_csr,desc_a,csr_prol,acsr3,desc_cprol,info)
  if (do_timings) call psb_toc(idx_spspmm)  

  if (debug) write(0,*) me,trim(name),' Done AxPROL ',&
       & acsr3%get_nrows(),acsr3%get_ncols(), acsr3%get_nzeros(),&
       & desc_cprol%get_local_rows(),desc_cprol%get_local_cols()

  !
  ! Ok first product done.
  !
  ! Remember that RESTR must be built from PROL after halo extension,
  ! which is done above in psb_par_spspmm
  if (debug) write(0,*)me,' ',name,' No inp_restr, transposing prol ',&
       & csr_prol%get_nrows(),csr_prol%get_ncols(),csr_prol%get_nzeros()
  call csr_prol%cp_to_fmt(coo_restr,info)
!!$      write(0,*)me,' ',name,' new into transposition ',coo_restr%get_nrows(),&
!!$           & coo_restr%get_ncols(),coo_restr%get_nzeros()
  call coo_restr%transp()
  nzl = coo_restr%get_nzeros()
  call desc_cprol%l2gip(coo_restr%ia(1:nzl),info)
  i=0
  !
  ! Now we have to fix this.  The only rows of the restrictor that are correct 
  ! are those corresponding to "local" aggregates, i.e. indices in ilaggr(:)
  !
  do k=1, nzl
    if ((naggrm1 < coo_restr%ia(k)) .and.(coo_restr%ia(k) <= naggrp1)) then
      i = i+1
      coo_restr%val(i) = coo_restr%val(k)
      coo_restr%ia(i)  = coo_restr%ia(k)
      coo_restr%ja(i)  = coo_restr%ja(k)
    end if
  end do
  call coo_restr%set_nzeros(i)
  call coo_restr%fix(info)
  call coo_restr%cp_to_coo(tmpcoo,info)
!!$      write(0,*)me,' ',name,' after transposition ',coo_restr%get_nrows(),coo_restr%get_ncols(),coo_restr%get_nzeros()

  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='spcnv coo_restr')
    goto 9999
  end if
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'starting sphalo/ rwxtd'
  nzl    = tmpcoo%get_nzeros()    
  call psb_glob_to_loc(tmpcoo%ia(1:nzl),desc_cprol,info,iact='I',owned=.true.)
  call tmpcoo%clean_negidx(info)
  nzl  = tmpcoo%get_nzeros()
  call tmpcoo%set_nrows(desc_cprol%get_local_rows())
  call tmpcoo%set_ncols(desc_a%get_local_cols())
!!$    write(0,*)me,' ',name,' after G2L on rows ',tmpcoo%get_nrows(),tmpcoo%get_ncols(),tmpcoo%get_nzeros()      
  call csr_restr%mv_from_coo(tmpcoo,info)

  if (debug) write(0,*) me,trim(name),' Product RESTRxAP ',&
       & csr_restr%get_nrows(),csr_restr%get_ncols(), &
       & desc_cprol%get_local_rows(),desc_a%get_local_cols(),&
       & acsr3%get_nrows(),acsr3%get_ncols()
  if (do_timings) call psb_tic(idx_spspmm)      
  call psb_par_spspmm(csr_restr,desc_a,acsr3,ac_csr,desc_cprol,info)
  if (do_timings) call psb_toc(idx_spspmm)      
  call csr_restr%free()
  call ac_csr%mv_to_coo(ac_coo,info)
  nza    = ac_coo%get_nzeros()
  if (debug) write(0,*) me,trim(name),' Fixing ac ',&
       & ac_coo%get_nrows(),ac_coo%get_ncols(), nza
  call ac_coo%fix(info)
  call desc_cprol%indxmap%l2gip(ac_coo%ia(1:nza),info)
  call desc_cprol%indxmap%l2gip(ac_coo%ja(1:nza),info)
  call ac_coo%set_nrows(ntaggr)
  call ac_coo%set_ncols(ntaggr)
  if (debug) write(0,*)  me,' ',trim(name),' Before mv_from',psb_get_errstatus()
  if (info == 0) call ac%mv_from(ac_coo)
  if (debug) write(0,*)  me,' ',trim(name),' After  mv_from',psb_get_errstatus()
  if (debug) write(0,*)  me,' ',trim(name),' ',ac%get_fmt(),ac%get_nrows(),ac%get_ncols(),ac%get_nzeros(),naggr,ntaggr
  ! write(0,*)  me,' ',trim(name),' Final AC newstyle ',ac%get_fmt(),ac%get_nrows(),ac%get_ncols(),ac%get_nzeros()

  nza = coo_prol%get_nzeros()
  call desc_cprol%indxmap%l2gip(coo_prol%ja(1:nza),info)
  
  if (debug) then
    write(0,*) me,' ',trim(name),' Checkpoint at exit'
    call psb_barrier(ictxt)
    write(0,*) me,' ',trim(name),' Checkpoint through'
  end if

  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='Build ac = coo_restr x am3')
    goto 9999
  end if


  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done spmm_bld_inner '

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
  
end subroutine mld_ls_spmm_bld_inner
