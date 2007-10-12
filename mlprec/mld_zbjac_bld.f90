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
!* This is where the action takes place.                                     *
!* ASMATBLD does the setup: building the prec descriptor plus retrieving     *
!*                           matrix rows if needed                           *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!* some open code does the renumbering                                       *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*****************************************************************************
subroutine mld_zbjac_bld(a,desc_a,p,upd,info)
  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_zbjac_bld

  implicit none
  !                                                                               
  !     .. Scalar Arguments ..                                                    
  integer, intent(out)                      :: info
  !     .. array Arguments ..                                                     
  type(psb_zspmat_type), intent(in), target :: a
  type(mld_zbaseprc_type), intent(inout)    :: p
  type(psb_desc_type), intent(in)           :: desc_a
  character, intent(in)                     :: upd

  !     .. Local Scalars ..                                                       
  integer  ::    i, j, jj, k, kk, m
  integer  ::    int_err(5)
  character ::        trans, unitd
  type(psb_zspmat_type) :: blck, atmp
  real(kind(1.d0)) :: t1,t2,t3,t4,t5,t6, t7, t8
  logical, parameter :: debugprt=.false., debug=.false., aggr_dump=.false.
  integer   nztota, nztotb, nztmp, nzl, nnr, ir, err_act,&
       & n_row, nrow_a,n_col, nhalo, ind, iind
  integer :: ictxt,np,me
  character(len=20)      :: name, ch_err
  character(len=5), parameter :: coofmt='COO', csrfmt='CSR'

  if(psb_get_errstatus().ne.0) return 
  info=0
  name='mld_zbjac_bld'
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)

  m = a%m
  if (m < 0) then
    info = 10
    int_err(1) = 1
    int_err(2) = m
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif
  trans = 'N'
  unitd = 'U'
  if (p%iprcparm(mld_n_ovr_) < 0) then
    info = 11
    int_err(1) = 1
    int_err(2) = p%iprcparm(mld_n_ovr_)
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  call psb_nullify_sp(blck)
  call psb_nullify_sp(atmp)

  t1= psb_wtime()

  if(debug) write(0,*)me,': calling mld_asmat_bld',&
       & p%iprcparm(mld_prec_type_),p%iprcparm(mld_n_ovr_)
  if (debug) call psb_barrier(ictxt)
  call mld_asmat_bld(p%iprcparm(mld_prec_type_),p%iprcparm(mld_n_ovr_),a,&
       & blck,desc_a,upd,p%desc_data,info,outfmt=csrfmt)

  if (debugprt) then 
    open(60+me)
    call psb_csprt(60+me,a,head='% A')
    close(60+me)
    open(70+me)
    call psb_csprt(70+me,blck,head='% BLCK')
    close(70+me)
  endif
  
  if(info/=0) then
    call psb_errpush(4010,name,a_err='mld_asmat_bld')
    goto 9999
  end if

  t2= psb_wtime()
  if (debug) write(0,*)me,': out of mld_asmat_bld'
  if (debug) call psb_barrier(ictxt)


  select case(p%iprcparm(mld_sub_ren_)) 

  case (1:)

    !
    ! Here we allocate a full copy to hold local A and received BLK
    ! Done inside sp_renum.
    !

    call  mld_sp_renum(a,desc_a,blck,p,atmp,info)

    if (info/=0) then
      call psb_errpush(4010,name,a_err='mld_sp_renum')
      goto 9999
    end if

    !------------------------------------------------------------------
    ! Split AC=M+N  N off-diagonal part
    ! Output in COO format. 
    call psb_sp_clip(atmp,p%av(mld_ap_nd_),info,&
         & jmin=atmp%m+1,rscale=.false.,cscale=.false.)

    call psb_spcnv(p%av(mld_ap_nd_),info,afmt='csr',dupl=psb_dupl_add_)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psb_spcnv csr 1')
      goto 9999
    end if

    k = psb_sp_get_nnzeros(p%av(mld_ap_nd_))
    call psb_sum(ictxt,k)

    if (k == 0) then 
      ! If the off diagonal part is emtpy, there's no point 
      ! in doing multiple  Jacobi sweeps. This is certain 
      ! to happen when running on a single processor.
      p%iprcparm(mld_smooth_sweeps_) = 1
    end if


    t3 = psb_wtime()
    if (debugprt) then 
      call psb_barrier(ictxt)
      open(40+me) 
      call psb_csprt(40+me,atmp,head='% Local matrix')
      close(40+me)
    endif
    if (debug) write(0,*) me,' Factoring rows ',&
         &atmp%m,a%m,blck%m,atmp%ia2(atmp%m+1)-1

    select case(p%iprcparm(mld_sub_solve_))

    case(mld_ilu_n_,mld_milu_n_,mld_ilu_t_) 

      call psb_spcnv(atmp,info,afmt='csr',dupl=psb_dupl_add_)
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='psb_spcnv csr 2')
        goto 9999
      end if

      call mld_ilu_bld(atmp,p%desc_data,p,upd,info)

      if (info/=0) then
        call psb_errpush(4010,name,a_err='mld_ilu_bld')
        goto 9999
      end if


      if (debugprt) then 

        open(80+me)

        call psb_csprt(80+me,p%av(mld_l_pr_),head='% Local L factor')
        write(80+me,*) '% Diagonal: ',p%av(mld_l_pr_)%m
        do i=1,p%av(mld_l_pr_)%m
          write(80+me,*) i,i,p%d(i)
        enddo
        call psb_csprt(80+me,p%av(mld_u_pr_),head='% Local U factor')

        close(80+me)
      endif


    case(mld_slu_)

      call psb_spcnv(atmp,info,afmt='csr',dupl=psb_dupl_add_)
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='psb_spcnv csr 3')
        goto 9999
      end if

      call mld_slu_bld(atmp,p%desc_data,p,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='mld_slu_bld')
        goto 9999
      end if

    case(mld_umf_)

      call psb_spcnv(atmp,info,afmt='csc',dupl=psb_dupl_add_)
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='psb_spcnv csc')
        goto 9999
      end if

      call mld_umf_bld(atmp,p%desc_data,p,info)
      if(debug) write(0,*)me,': Done mld_umf_bld ',info
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='mld_umf_bld')
        goto 9999
      end if

    case(mld_f_none_) 
      info=4010
      call psb_errpush(info,name,a_err='Inconsistent prec  mld_f_none_')
      goto 9999

    case default
      info=4010
      call psb_errpush(info,name,a_err='Unknown mld_sub_solve_')
      goto 9999
    end select



    call psb_sp_free(atmp,info) 

    if(info/=0) then
      call psb_errpush(4010,name,a_err='psb_sp_free')
      goto 9999
    end if



  case(0)  ! No renumbering

    select case(p%iprcparm(mld_sub_solve_))

    case(mld_ilu_n_,mld_milu_n_,mld_ilu_t_) 


      if (p%iprcparm(mld_smooth_sweeps_) > 1) then 
        n_row = psb_cd_get_local_rows(p%desc_data)
        n_col = psb_cd_get_local_cols(p%desc_data)
        nrow_a = a%m 
        ! The following is known to work 
        ! given that the output from CLIP is in COO. 
        call psb_sp_clip(a,p%av(mld_ap_nd_),info,&
             & jmin=nrow_a+1,rscale=.false.,cscale=.false.)
        call psb_sp_clip(blck,atmp,info,&
             & jmin=nrow_a+1,rscale=.false.,cscale=.false.)
        call psb_rwextd(n_row,p%av(mld_ap_nd_),info,b=atmp) 

        call psb_spcnv(p%av(mld_ap_nd_),info,afmt='csr',dupl=psb_dupl_add_)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_spcnv csr 4')
          goto 9999
        end if
        
        k = psb_sp_get_nnzeros(p%av(mld_ap_nd_))
        call psb_sum(ictxt,k)

        if (k == 0) then 
          ! If the off diagonal part is emtpy, there's no point 
          ! in doing multiple  Jacobi sweeps. This is certain 
          ! to happen when running on a single processor.
          p%iprcparm(mld_smooth_sweeps_) = 1
        end if
        call psb_sp_free(atmp,info) 
      end if

      call mld_ilu_bld(a,desc_a,p,upd,info,blck=blck)

      if(info/=0) then
        call psb_errpush(4010,name,a_err='mld_ilu_bld')
        goto 9999
      end if


      if (debugprt) then 

        open(80+me)

        call psb_csprt(80+me,p%av(mld_l_pr_),head='% Local L factor')
        write(80+me,*) '% Diagonal: ',p%av(mld_l_pr_)%m
        do i=1,p%av(mld_l_pr_)%m
          write(80+me,*) i,i,p%d(i)
        enddo
        call psb_csprt(80+me,p%av(mld_u_pr_),head='% Local U factor')

        close(80+me)
      endif


    case(mld_slu_)

      call psb_spcnv(a,atmp,info,afmt='coo')
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='psb_spcnv')
        goto 9999
      end if

      n_row = psb_cd_get_local_rows(p%desc_data)
      n_col = psb_cd_get_local_cols(p%desc_data)
      call psb_rwextd(n_row,atmp,info,b=blck) 

      if (p%iprcparm(mld_smooth_sweeps_) > 1) then 
        !------------------------------------------------------------------
        ! Split AC=M+N  N off-diagonal part
        ! Output in COO format. 
        call psb_sp_clip(atmp,p%av(mld_ap_nd_),info,&
             & jmin=atmp%m+1,rscale=.false.,cscale=.false.)

        call psb_spcnv(p%av(mld_ap_nd_),info,afmt='csr',dupl=psb_dupl_add_)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_spcnv csr 6')
          goto 9999
        end if

        k = psb_sp_get_nnzeros(p%av(mld_ap_nd_))
        call psb_sum(ictxt,k)

        if (k == 0) then 
          ! If the off diagonal part is emtpy, there's no point 
          ! in doing multiple  Jacobi sweeps. This is certain 
          ! to happen when running on a single processor.
          p%iprcparm(mld_smooth_sweeps_) = 1
        end if
      endif

      if (info == 0) call psb_spcnv(atmp,info,afmt='csr',dupl=psb_dupl_add_)
      if (info == 0) call mld_slu_bld(atmp,p%desc_data,p,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='mld_slu_bld')
        goto 9999
      end if

      call psb_sp_free(atmp,info) 
      if(info/=0) then
        call psb_errpush(4010,name,a_err='psb_sp_free')
        goto 9999
      end if


    case(mld_sludist_)

      call psb_spcnv(a,atmp,info,afmt='coo')
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='psb_spcnv')
        goto 9999
      end if
      
      n_row = psb_cd_get_local_rows(p%desc_data)
      n_col = psb_cd_get_local_cols(p%desc_data)
      call psb_rwextd(n_row,atmp,info,b=blck) 

      if (p%iprcparm(mld_smooth_sweeps_) > 1) then 
        !------------------------------------------------------------------
        ! Split AC=M+N  N off-diagonal part
        ! Output in COO format. 
        call psb_sp_clip(atmp,p%av(mld_ap_nd_),info,&
             & jmin=atmp%m+1,rscale=.false.,cscale=.false.)

        call psb_spcnv(p%av(mld_ap_nd_),info,afmt='csr',dupl=psb_dupl_add_)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_spcnv csr 7')
          goto 9999
        end if

        k = psb_sp_get_nnzeros(p%av(mld_ap_nd_))
        call psb_sum(ictxt,k)

        if (k == 0) then 
          ! If the off diagonal part is emtpy, there's no point 
          ! in doing multiple  Jacobi sweeps. This is certain 
          ! to happen when running on a single processor.
          p%iprcparm(mld_smooth_sweeps_) = 1
        end if
      endif
      
      if (info == 0) call psb_spcnv(atmp,info,afmt='csr',dupl=psb_dupl_add_)
      if (info == 0) call mld_sludist_bld(atmp,p%desc_data,p,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='mld_slu_bld')
        goto 9999
      end if

      call psb_sp_free(atmp,info) 
      if(info/=0) then
        call psb_errpush(4010,name,a_err='psb_sp_free')
        goto 9999
      end if

    case(mld_umf_)


      call psb_spcnv(a,atmp,info,afmt='coo')
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='psb_spcnv')
        goto 9999
      end if

      n_row = psb_cd_get_local_rows(p%desc_data)
      n_col = psb_cd_get_local_cols(p%desc_data)
      call psb_rwextd(n_row,atmp,info,b=blck) 

      if (p%iprcparm(mld_smooth_sweeps_) > 1) then 
        !------------------------------------------------------------------
        ! Split AC=M+N  N off-diagonal part
        ! Output in COO format. 
!!$        write(0,*) 'mld_bjac_bld:' size(p%av),mld_ap_nd_
        call psb_sp_clip(atmp,p%av(mld_ap_nd_),info,&
             & jmin=atmp%m+1,rscale=.false.,cscale=.false.)

        call psb_spcnv(p%av(mld_ap_nd_),info,afmt='csr',dupl=psb_dupl_add_)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_spcnv csr 8')
          goto 9999
        end if

        k = psb_sp_get_nnzeros(p%av(mld_ap_nd_))
        call psb_sum(ictxt,k)

        if (k == 0) then 
          ! If the off diagonal part is emtpy, there's no point 
          ! in doing multiple  Jacobi sweeps. This is certain 
          ! to happen when running on a single processor.
          p%iprcparm(mld_smooth_sweeps_) = 1
        end if
      endif

      if (info == 0) call psb_ipcoo2csc(atmp,info,clshr=.true.)
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='psb_ipcoo2csc')
        goto 9999
      end if

      call mld_umf_bld(atmp,p%desc_data,p,info)
      if(debug) write(0,*)me,': Done mld_umf_bld ',info
      if (info /= 0) then
        call psb_errpush(4010,name,a_err='mld_umf_bld')
        goto 9999
      end if

      call psb_sp_free(atmp,info) 
      if(info/=0) then
        call psb_errpush(4010,name,a_err='psb_sp_free')
        goto 9999
      end if


    case(mld_f_none_) 
      info=4010
      call psb_errpush(info,name,a_err='Inconsistent prec  mld_f_none_')
      goto 9999

    case default
      info=4010
      call psb_errpush(info,name,a_err='Unknown mld_sub_solve_')
      goto 9999
    end select

  case default
    info=4010
    call psb_errpush(info,name,a_err='Invalid renum_')
    goto 9999

  end select

  t6 = psb_wtime()

  call psb_sp_free(blck,info)
  if(info/=0) then
    call psb_errpush(4010,name,a_err='psb_sp_free')
    goto 9999
  end if

  if (debug) write(0,*) me,'End of ilu_bld'

  call psb_erractionrestore(err_act)

  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return


end subroutine mld_zbjac_bld


