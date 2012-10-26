subroutine mld_s_ilu_solver_apply_vect(alpha,sv,x,beta,y,desc_data,trans,work,info)
  
  use psb_base_mod
  use mld_s_ilu_solver, mld_protect_name => mld_s_ilu_solver_apply_vect
  implicit none 
  type(psb_desc_type), intent(in)             :: desc_data
  class(mld_s_ilu_solver_type), intent(inout) :: sv
  type(psb_s_vect_type),intent(inout)         :: x
  type(psb_s_vect_type),intent(inout)         :: y
  real(psb_spk_),intent(in)                   :: alpha,beta
  character(len=1),intent(in)                 :: trans
  real(psb_spk_),target, intent(inout)        :: work(:)
  integer, intent(out)                        :: info

  integer    :: n_row,n_col
  type(psb_s_vect_type) :: wv, wv1
  real(psb_spk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer    :: ictxt,np,me,i, err_act
  character          :: trans_
  character(len=20)  :: name='s_ilu_solver_apply'

  call psb_erractionsave(err_act)

  info = psb_success_

  trans_ = psb_toupper(trans)
  select case(trans_)
  case('N')
  case('T')
  case('C')
  case default
    call psb_errpush(psb_err_iarg_invalid_i_,name)
    goto 9999
  end select

  n_row = desc_data%get_local_rows()
  n_col = desc_data%get_local_cols()


  if (x%get_nrows() < n_row) then 
    info = 36
    call psb_errpush(info,name,i_err=(/2,n_row,0,0,0/))
    goto 9999
  end if
  if (y%get_nrows() < n_row) then 
    info = 36
    call psb_errpush(info,name,i_err=(/3,n_row,0,0,0/))
    goto 9999
  end if
  if (.not.allocated(sv%dv%v)) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner: D")
    goto 9999
  end if
  if (sv%dv%get_nrows() < n_row) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner: DV")
    goto 9999
  end if



  if (n_col <= size(work)) then 
    ww => work(1:n_col)
    if ((4*n_col+n_col) <= size(work)) then 
      aux => work(n_col+1:)
    else
      allocate(aux(4*n_col),stat=info)
    endif
  else
    allocate(ww(n_col),aux(4*n_col),stat=info)
  endif

  if (info /= psb_success_) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/5*n_col,0,0,0,0/),&
         & a_err='real(psb_spk_)')
    goto 9999      
  end if

  call psb_geasb(wv,desc_data,info,mold=x%v,scratch=.true.) 
  call psb_geasb(wv1,desc_data,info,mold=x%v,scratch=.true.) 

  select case(trans_)
  case('N')
    call psb_spsm(sone,sv%l,x,szero,wv,desc_data,info,&
         & trans=trans_,scale='L',diag=sv%dv,choice=psb_none_,work=aux)

    if (info == psb_success_) call psb_spsm(alpha,sv%u,wv,beta,y,desc_data,info,&
         & trans=trans_,scale='U',choice=psb_none_, work=aux)

  case('T')
    call psb_spsm(sone,sv%u,x,szero,wv,desc_data,info,&
         & trans=trans_,scale='L',diag=sv%dv,choice=psb_none_,work=aux)
    if (info == psb_success_) call psb_spsm(alpha,sv%l,wv,beta,y,desc_data,info,&
         & trans=trans_,scale='U',choice=psb_none_,work=aux)

  case('C')

    call psb_spsm(sone,sv%u,x,szero,wv,desc_data,info,&
         & trans=trans_,scale='U',choice=psb_none_,work=aux)

    call wv1%mlt(sone,sv%dv,wv,szero,info,conjgx=trans_)

    if (info == psb_success_) call psb_spsm(alpha,sv%l,wv1,beta,y,desc_data,info,&
         & trans=trans_,scale='U',choice=psb_none_,work=aux)

  case default
    call psb_errpush(psb_err_internal_error_,name,a_err='Invalid TRANS in ILU subsolve')
    goto 9999
  end select


  if (info /= psb_success_) then

    call psb_errpush(psb_err_internal_error_,name,a_err='Error in subsolve')
    goto 9999
  endif
  call wv%free(info)
  call wv1%free(info)
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
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_s_ilu_solver_apply_vect
