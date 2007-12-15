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
! File: mld_dasmat_bld.f90.
!
! Subroutine: mld_dasmat_bld.
! Version:    real.
!
!  This routine builds the communication descriptor associated to the extended
!  matrices that form the Additive Schwarz (AS) preconditioner and retrieves
!  the remote pieces needed to build the local extended matrix. If the
!  preconditioner is the block-Jacobi one, the routine makes only a copy of
!  the descriptor of the original matrix.
!    
!
! Arguments:
!    ptype      -  integer, input.
!                  The type of preconditioner to be built. Only the values
!                  mld_bjac_ and mld_as_ (see mld_prec_type.f90) are allowed.
!    novr       -  integer, input.
!                  The number of overlap layers in the AS preconditioner.
!    a          -  type(psb_dspmat_type), input.
!                  The sparse matrix structure containing the local part of the
!                  matrix to be preconditioned.
!    blk        -  type(psb_dspmat_type), output.
!                  The sparse matrix structure containing the remote rows that
!                  extend the local matrix according to novr. If novr = 0 then
!                  blk does not contain any row.
!    desc_data  -  type(psb_desc_type), input.
!                  The communication descriptor of the sparse matrix a.
!       desc_p  -  type(psb_desc_type), output.
!                  The communication descriptor associated to the extended 
!                  matrices that form the AS preconditioner.
!    info       -  integer, output.
!                  Error code.
!    outfmt     -  character(len=5), optional.
!                  The storage format of the local extended matrix for the AS
!                  preconditioner. Currently outfmt is set to 'CSR' by the
!                  calling routine mld_bjac_bld.                 
!  
subroutine mld_dasmat_bld(ptype,novr,a,blk,desc_data,upd,desc_p,info,outfmt)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_dasmat_bld
  
  Implicit None

! Arguments
  integer, intent(in)                  :: ptype,novr
  Type(psb_dspmat_type), Intent(in)    :: a
  Type(psb_dspmat_type), Intent(inout) :: blk
  integer, intent(out)                 :: info
  Type(psb_desc_type), Intent(inout)   :: desc_p
  Type(psb_desc_type), Intent(in)      :: desc_data 
  Character, Intent(in)                :: upd
  character(len=5), optional           :: outfmt

! Local variables
  real(kind(1.d0)) :: t1,t2,t3
  integer   icomm
  Integer ::  np,me,nnzero,&
       &  ictxt, n_col,int_err(5),&
       &  tot_recv, n_row,nhalo, nrow_a,err_act
  Logical,Parameter :: debug=.false.
  character(len=20) :: name, ch_err

  name='mld_dasmat_bld'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  If(debug) Write(0,*)'IN DASMATBLD  ', upd
  ictxt = psb_cd_get_context(desc_data)
  icomm = psb_cd_get_mpic(desc_data)

  Call psb_info(ictxt, me, np)

  tot_recv=0

  nrow_a = desc_data%matrix_data(psb_n_row_)
  nnzero = Size(a%aspk)
  n_col  = desc_data%matrix_data(psb_n_col_)
  nhalo  = n_col-nrow_a

  select case (ptype)
  case (mld_bjac_) 
    !
    ! Block-Jacobi preconditioner. Copy the descriptor, just in case 
    ! we want to renumber the rows and columns of the matrix. 
    !
    If(debug) Write(0,*)' asmatbld calling allocate '
    call psb_sp_all(0,0,blk,1,info)
    if(info /= 0) then
      info=4010
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    blk%fida        = 'COO'
    blk%infoa(psb_nnz_) = 0
    If(debug) Write(0,*)' asmatbld done spallocate'
    If (upd == 'F') Then
      call psb_cdcpy(desc_data,desc_p,info)
      If(debug) Write(0,*)' asmatbld done cdcpy'
      if(info /= 0) then
        info=4010
        ch_err='psb_cdcpy'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    endif

  case(mld_as_) 


    !
    ! Additive Schwarz 
    !

    if (novr < 0) then
      info=3
      int_err(1)=novr
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    endif

    if ((novr == 0).or.(np==1)) then 
      !
      ! Actually, this is just block Jacobi
      !
      If(debug) Write(0,*)' asmatbld calling allocate novr=0'
      call psb_sp_all(0,0,blk,1,info)
      if(info /= 0) then
        info=4010
        ch_err='psb_sp_all'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
      blk%fida            = 'COO'
      blk%infoa(psb_nnz_) = 0
      if (debug) write(0,*) 'Calling desccpy'
      if (upd == 'F') then 
        call psb_cdcpy(desc_data,desc_p,info)
        If(debug) Write(0,*)' asmatbld done cdcpy'
        if(info /= 0) then
          info=4010
          ch_err='psb_cdcpy'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        end if
        if (debug) write(0,*) 'Early return from asmatbld: P>=3 N_OVR=0'
      endif
      return
    endif


    If(debug)Write(0,*)'BEGIN dasmatbld',me,upd,novr
    t1 = psb_wtime()

    If (upd == 'F') Then
      !
      ! Build the auxiliary descriptor desc_p%matrix_data(psb_n_row_).
      ! This is done by psb_cdbldext (interface to psb_cdovr), which is
      ! independent of CSR, and has been placed in the tools directory
      ! of PSBLAS, instead of the mlprec directory of MLD2P4, because it
      ! might be used independently of the AS preconditioner, to build
      ! a descriptor for an extended stencil in a PDE solver. 
      !
      call psb_cdbldext(a,desc_data,novr,desc_p,info,extype=psb_ovt_asov_)
      if(info /= 0) then
        info=4010
        ch_err='psb_cdbldext'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    Endif

    if(debug) write(0,*) me,' From cdbldext _:',desc_p%matrix_data(psb_n_row_),&
         & desc_p%matrix_data(psb_n_col_)

    !
    ! Retrieve the remote pieces needed to build the local extended matrix
    !

    n_row = desc_p%matrix_data(psb_n_row_)
    t2 = psb_wtime()

    if (debug) write(0,*) 'Before sphalo ',blk%fida,blk%m,psb_nnz_,blk%infoa(psb_nnz_)
    Call psb_sphalo(a,desc_p,blk,info,&
         & outfmt=outfmt,data=psb_comm_ext_,rowscale=.true.)

    if(info /= 0) then
      info=4010
      ch_err='psb_sphalo'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if (debug) write(0,*) 'After psb_sphalo ',&
         & blk%fida,blk%m,psb_nnz_,blk%infoa(psb_nnz_)
  case default
    if(info /= 0) then
      info=4000
      ch_err='Invalid ptype'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    
  End select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  Return

End Subroutine mld_dasmat_bld

