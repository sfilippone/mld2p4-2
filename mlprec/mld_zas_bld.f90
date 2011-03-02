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
! File: mld_zas_bld.f90
!
! Subroutine: mld_zas_bld
! Version:    complex
!
!  This routine builds Additive Schwarz (AS) preconditioners. If the AS
!  preconditioner is actually the block-Jacobi one, the routine makes only a
!  copy of the descriptor of the original matrix and then calls mld_fact_bld
!  to perform an LU or ILU factorization of the diagonal blocks of the
!  distributed matrix.
!    
!
! Arguments:
!    a          -  type(psb_zspmat_type), input.
!                  The sparse matrix structure containing the local part of the
!                  matrix to be preconditioned.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the sparse matrix a.
!    p          -  type(mld_zbaseprec_type), input/output.
!                  The 'base preconditioner' data structure containing the local
!                  part of the preconditioner or solver to be built.
!    upd        -  character, input.
!                  If upd='F' then the preconditioner is built from scratch;
!                  if upd=T' then the matrix to be preconditioned has the same
!                  sparsity pattern of a matrix that has been previously
!                  preconditioned, hence some information is reused in building
!                  the new preconditioner.
!    info       -  integer, output.
!                  Error code.
!  
subroutine mld_zas_bld(a,desc_a,p,upd,info)

  use psb_sparse_mod
  use mld_z_inner_mod, mld_protect_name => mld_zas_bld

  Implicit None

  ! Arguments
  type(psb_zspmat_type), intent(in), target :: a
  Type(psb_desc_type), Intent(in)           :: desc_a 
  type(mld_zbaseprec_type), intent(inout)    :: p
  character, intent(in)                     :: upd
  integer, intent(out)                      :: info

  ! Local variables
  integer :: ptype,novr
  integer :: icomm
  Integer ::  np,me,nnzero,ictxt, int_err(5),&
       &  tot_recv, n_row,n_col,nhalo, err_act, data_
  type(psb_zspmat_type) :: blck
  integer           :: debug_level, debug_unit
  character(len=20) :: name, ch_err

  name='mld_as_bld'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  If (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ' start ', upd
  ictxt = psb_cd_get_context(desc_a)
  icomm = psb_cd_get_mpic(desc_a)

  Call psb_info(ictxt, me, np)

  tot_recv=0

  n_row = psb_cd_get_local_rows(desc_a)
  n_col  = psb_cd_get_local_cols(desc_a)
  nnzero = psb_sp_get_nnzeros(a)
  nhalo  = n_col-n_row
  ptype  = p%iprcparm(mld_smoother_type_)
  novr   = p%iprcparm(mld_sub_ovr_)

  select case (ptype)

  case(mld_bjac_) 
    !
    ! Block Jacobi
    !
    data_ = psb_no_comm_
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Calling desccpy'
    if (upd == 'F') then 
      call psb_cdcpy(desc_a,p%desc_data,info)
      If(debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & '  done cdcpy'
      if(info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='psb_cdcpy'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & 'Early return: P>=3 N_OVR=0'
    endif
    call psb_sp_all(0,0,blck,1,info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    blck%fida            = 'COO'
    blck%infoa(psb_nnz_) = 0
    
    call mld_fact_bld(a,p,upd,info,blck=blck)

    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='mld_fact_bld'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if


  case(mld_as_) 
    !
    ! Additive Schwarz 
    !
    if (novr < 0) then
      info=psb_err_invalid_ovr_num_
      int_err(1)=novr
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    endif

    if ((novr == 0).or.(np == 1)) then 
      !
      ! Actually, this is just block Jacobi
      !
      data_ = psb_no_comm_
      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & 'Calling desccpy'
      if (upd == 'F') then 
        call psb_cdcpy(desc_a,p%desc_data,info)
        If(debug_level >= psb_debug_outer_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & '  done cdcpy'
        if(info /= psb_success_) then
          info=psb_err_from_subroutine_
          ch_err='psb_cdcpy'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        end if
        if (debug_level >= psb_debug_outer_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & 'Early return: P>=3 N_OVR=0'
      endif
      call psb_sp_all(0,0,blck,1,info)
      if(info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='psb_sp_all'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
      blck%fida            = 'COO'
      blck%infoa(psb_nnz_) = 0

    else

      If (upd == 'F') Then
        !
        ! Build the auxiliary descriptor desc_p%matrix_data(psb_n_row_).
        ! This is done by psb_cdbldext (interface to psb_cdovr), which is
        ! independent of CSR, and has been placed in the tools directory
        ! of PSBLAS, instead of the mlprec directory of MLD2P4, because it
        ! might be used independently of the AS preconditioner, to build
        ! a descriptor for an extended stencil in a PDE solver. 
        !
        call psb_cdbldext(a,desc_a,novr,p%desc_data,info,extype=psb_ovt_asov_)
        if(debug_level >= psb_debug_outer_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & ' From cdbldext _:',psb_cd_get_local_rows(p%desc_data),&
             & psb_cd_get_local_cols(p%desc_data)
        
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_
          ch_err='psb_cdbldext'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        end if
      Endif

      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & 'Before sphalo ',blck%fida,blck%m,psb_nnz_,blck%infoa(psb_nnz_)

      !
      ! Retrieve the remote sparse matrix rows required for the AS extended
      ! matrix
      data_ = psb_comm_ext_
      Call psb_sphalo(a,p%desc_data,blck,info,data=data_,rowscale=.true.)
      
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='psb_sphalo'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
      
      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & 'After psb_sphalo ',&
           & blck%fida,blck%m,psb_nnz_,blck%infoa(psb_nnz_)

    End if


    call mld_fact_bld(a,p,upd,info,blck=blck)

    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='mld_fact_bld'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  case default

    info=psb_err_internal_error_
    ch_err='Invalid ptype'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999

  End select
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
  Return

End Subroutine mld_zas_bld

