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
!*****************************************************************************
!*                                                                           *
!* This routine does two things:                                             *
!*    1. Builds the auxiliary descriptor. This is always done even for       *
!*       Block Jacobi.                                                       *
!*    2. Retrieves the remote matrix pieces.                                 *
!*                                                                           *
!*    All of 1. is done under psb_cdovr, which is independent of CSR, and    *
!*    has been placed in the TOOLS directory because it might be used for    *
!*    building a descriptor for an extended stencil in a PDE solver without  *
!*    necessarily applying AS precond.                                       *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*****************************************************************************
Subroutine mld_zasmat_bld(ptype,novr,a,blk,desc_data,upd,desc_p,info,outfmt)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_zasmat_bld

  Implicit None

  !     .. Array Arguments ..
  integer, intent(in)                  :: ptype,novr
  Type(psb_zspmat_type), Intent(in)    :: a
  Type(psb_zspmat_type), Intent(inout) :: blk
  integer, intent(out)                 :: info
  Type(psb_desc_type), Intent(inout)   :: desc_p
  Type(psb_desc_type), Intent(in)      :: desc_data 
  Character, Intent(in)                :: upd
  character(len=5), optional           :: outfmt


  real(kind(1.d0)) :: t1,t2,t3
  integer   icomm

  !     .. Local Scalars ..
  Integer ::  k, np,me,m,nnzero,&
       &  ictxt, n_col,ier,n,int_err(5),&
       &  tot_recv, ircode, n_row,nhalo, nrow_a,err_act
  Logical,Parameter :: debug=.false., debugprt=.false.
  character(len=20) :: name, ch_err

  name='mld_zasmat_bld'
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


  If (ptype == mld_bjac_) Then
    !
    ! Block Jacobi. Copy the descriptor, just in case we want to
    ! do the renumbering. 
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

  Else If (ptype == mld_as_) Then


    !
    ! Additive Schwarz variant. 
    !
    !


    if (novr < 0) then
      info=3
      int_err(1)=novr
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    endif

    if ((novr == 0).or.(np==1)) then 
      !
      ! This is really just Block Jacobi.....
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
      !  Build the  auiliary descriptor',desc_p%matrix_data(psb_n_row_)
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


    n_row = desc_p%matrix_data(psb_n_row_)
    t2 = psb_wtime()

    if (debug) write(0,*) 'Before sphalo ',blk%fida,blk%m,psb_nnz_,blk%infoa(psb_nnz_)

    if (present(outfmt)) then 
      if(debug) write(0,*) me,': Calling outfmt SPHALO with ',size(blk%ia2)
      Call psb_sphalo(a,desc_p,blk,info,&
           & outfmt=outfmt,data=psb_comm_ext_,rowscale=.true.)
    else
      if(debug) write(0,*) me,': Calling SPHALO with ',size(blk%ia2)
      Call psb_sphalo(a,desc_p,blk,info,&
           & data=psb_comm_ext_,rowscale=.true.)
    end if


    if(info /= 0) then
      info=4010
      ch_err='psb_sphalo'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if (debug) write(0,*) 'After psb_sphalo ',&
         & blk%fida,blk%m,psb_nnz_,blk%infoa(psb_nnz_)

    t3 = psb_wtime()
    if (debugprt) then 
      open(40+me) 
      call psb_csprt(40+me,blk,head='% Ovrlap rows')
      close(40+me)
    endif


  End If

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  Return

End Subroutine mld_zasmat_bld

