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
subroutine mld_dbjac_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
  !
  !  Compute   Y <-  beta*Y + alpha*K^-1 X 
  !  where K is a a Block Jacobi  preconditioner stored in prec
  !  Note that desc_data may or may not be the same as prec%desc_data,
  !  but since both are INTENT(IN) this should be legal. 
  ! 

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_dbjac_aply

  implicit none 

  type(psb_desc_type), intent(in)       :: desc_data
  type(mld_dbaseprc_type), intent(in)   :: prec
  real(kind(0.d0)),intent(inout)        :: x(:), y(:)
  real(kind(0.d0)),intent(in)           :: alpha,beta
  character(len=1)                      :: trans
  real(kind(0.d0)),target               :: work(:)
  integer, intent(out)                  :: info

  ! Local variables
  integer :: n_row,n_col
  real(kind(1.d0)), pointer :: ww(:), aux(:), tx(:),ty(:),tb(:)
  character     ::diagl, diagu
  integer :: ictxt,np,me,i, nrg, err_act, int_err(5)
  real(kind(1.d0)) :: t1, t2, t3, t4, t5, t6, t7
  logical,parameter   :: debug=.false., debugprt=.false.
  character(len=20)   :: name, ch_err

  name='mld_dbjac_aply'
  info = 0
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_data)
  call psb_info(ictxt, me, np)

  select case(toupper(trans))
  case('N')
  case('T','C')
  case default
    call psb_errpush(40,name)
    goto 9999
  end select


  n_row = psb_cd_get_local_rows(desc_data)
  n_col = psb_cd_get_local_cols(desc_data)

  if (n_col <= size(work)) then 
    ww => work(1:n_col)
    if ((4*n_col+n_col) <= size(work)) then 
      aux => work(n_col+1:)
    else
      allocate(aux(4*n_col),stat=info)
      if (info /= 0) then 
        info=4025
        call psb_errpush(info,name,i_err=(/4*n_col,0,0,0,0/),&
             & a_err='real(kind(1.d0))')
        goto 9999      
      end if

    endif
  else
    allocate(ww(n_col),aux(4*n_col),stat=info)
    if (info /= 0) then 
      info=4025
      call psb_errpush(info,name,i_err=(/5*n_col,0,0,0,0/),&
           & a_err='real(kind(1.d0))')
      goto 9999      
    end if
  endif

  if (debug) then 
    write(0,*) me,' mld_bjac_APLY: ',prec%iprcparm(mld_sub_solve_),prec%iprcparm(mld_smooth_sweeps_)
  end if
  
  if (prec%iprcparm(mld_smooth_sweeps_) == 1) then 

    select case(prec%iprcparm(mld_sub_solve_))
    case(mld_ilu_n_,mld_milu_n_,mld_ilu_t_) 

      select case(toupper(trans))
      case('N')

        call psb_spsm(done,prec%av(mld_l_pr_),x,dzero,ww,desc_data,info,&
             & trans='N',unit='L',diag=prec%d,choice=psb_none_,work=aux)
        if(info /=0) goto 9999
        call psb_spsm(alpha,prec%av(mld_u_pr_),ww,beta,y,desc_data,info,&
             & trans='N',unit='U',choice=psb_none_, work=aux)
        if(info /=0) goto 9999

      case('T','C')
        call psb_spsm(done,prec%av(mld_u_pr_),x,dzero,ww,desc_data,info,&
             & trans=trans,unit='L',diag=prec%d,choice=psb_none_,work=aux)
        if(info /=0) goto 9999
        call psb_spsm(alpha,prec%av(mld_l_pr_),ww,beta,y,desc_data,info,&
             & trans=trans,unit='U',choice=psb_none_,work=aux)
        if(info /=0) goto 9999

      end select

    case(mld_slu_)

      ww(1:n_row) = x(1:n_row)

      select case(toupper(trans))
      case('N')
        call mld_dslu_solve(0,n_row,1,ww,n_row,prec%iprcparm(mld_slu_ptr_),info)
      case('T','C')
        call mld_dslu_solve(1,n_row,1,ww,n_row,prec%iprcparm(mld_slu_ptr_),info)
      end select

      if(info /=0) goto 9999
      call psb_geaxpby(alpha,ww,beta,y,desc_data,info)

    case(mld_sludist_)

!!$      write(0,*) 'Calling mld_sludist_solve ',n_row
      ww(1:n_row) = x(1:n_row)

      select case(toupper(trans))
      case('N')
        call mld_dsludist_solve(0,n_row,1,ww,n_row,prec%iprcparm(mld_slud_ptr_),info)
      case('T','C')
        call mld_dsludist_solve(1,n_row,1,ww,n_row,prec%iprcparm(mld_slud_ptr_),info)
      end select

      if(info /=0) goto 9999
      call psb_geaxpby(alpha,ww,beta,y,desc_data,info)

    case (mld_umf_) 


      select case(toupper(trans))
      case('N')
        call mld_dumf_solve(0,n_row,ww,x,n_row,prec%iprcparm(mld_umf_numptr_),info)
      case('T','C')
        call mld_dumf_solve(1,n_row,ww,x,n_row,prec%iprcparm(mld_umf_numptr_),info)
      end select

      if(info /=0) goto 9999

      call psb_geaxpby(alpha,ww,beta,y,desc_data,info)

    case default
      write(0,*) 'Unknown factorization type in mld_bjac_aply',prec%iprcparm(mld_sub_solve_)
    end select
    if (debugprt) write(0,*)' Y: ',y(:)

  else if (prec%iprcparm(mld_smooth_sweeps_) > 1) then 

    ! Note: we have to add TRANS to this one !!!!!!!!! 

    if (size(prec%av) < mld_ap_nd_) then 
      info = 4011
      goto 9999
    endif

    allocate(tx(n_col),ty(n_col),stat=info)
    if (info /= 0) then 
      info=4025
      call psb_errpush(info,name,i_err=(/2*n_col,0,0,0,0/),&
           & a_err='real(kind(1.d0))')
      goto 9999      
    end if

    tx = dzero
    ty = dzero
    select case(prec%iprcparm(mld_sub_solve_)) 
    case(mld_ilu_n_,mld_milu_n_,mld_ilu_t_) 
      do i=1, prec%iprcparm(mld_smooth_sweeps_) 
        !   X(k+1) = M^-1*(b-N*X(k))
        ty(1:n_row) = x(1:n_row)
        call psb_spmm(-done,prec%av(mld_ap_nd_),tx,done,ty,&
             &   prec%desc_data,info,work=aux)
        if(info /=0) goto 9999
        call psb_spsm(done,prec%av(mld_l_pr_),ty,dzero,ww,&
             & prec%desc_data,info,&
             & trans='N',unit='L',diag=prec%d,choice=psb_none_,work=aux)
        if(info /=0) goto 9999
        call psb_spsm(done,prec%av(mld_u_pr_),ww,dzero,tx,&
             & prec%desc_data,info,&
             & trans='N',unit='U',choice=psb_none_,work=aux)
        if(info /=0) goto 9999
      end do

    case(mld_sludist_) 
      write(0,*) 'No sense in having SLUDist with Jmld_ac_SWEEPS >1'
      info=4010
      goto 9999
    case(mld_slu_) 
      do i=1, prec%iprcparm(mld_smooth_sweeps_) 
        !   X(k+1) = M^-1*(b-N*X(k))
        ty(1:n_row) = x(1:n_row)
        call psb_spmm(-done,prec%av(mld_ap_nd_),tx,done,ty,&
             &   prec%desc_data,info,work=aux)
        if(info /=0) goto 9999

        call mld_dslu_solve(0,n_row,1,ty,n_row,prec%iprcparm(mld_slu_ptr_),info)
        if(info /=0) goto 9999
        tx(1:n_row) = ty(1:n_row)        
      end do
    case(mld_umf_) 
      do i=1, prec%iprcparm(mld_smooth_sweeps_) 
        !   X(k+1) = M^-1*(b-N*X(k))
        ty(1:n_row) = x(1:n_row)
        call psb_spmm(-done,prec%av(mld_ap_nd_),tx,done,ty,&
             &   prec%desc_data,info,work=aux)
        if(info /=0) goto 9999

        call mld_dumf_solve(0,n_row,ww,ty,n_row,&
             & prec%iprcparm(mld_umf_numptr_),info)
        if(info /=0) goto 9999
        tx(1:n_row) = ww(1:n_row)        
      end do

    end select
    
    call psb_geaxpby(alpha,tx,beta,y,desc_data,info)
    

    deallocate(tx,ty)


  else
    info = 10
    call psb_errpush(info,name,&
         & i_err=(/2,prec%iprcparm(mld_smooth_sweeps_),0,0,0/))
    goto 9999

  endif

  if (n_col <= size(work)) then 
    if ((4*n_col+n_col) <= size(work)) then 
    else
      deallocate(aux)
    endif
  else
    deallocate(ww,aux)
  endif


  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name,i_err=int_err,a_err=ch_err)
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_dbjac_aply

