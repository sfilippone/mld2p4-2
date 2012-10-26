subroutine mld_d_jac_smoother_apply_vect(alpha,sm,x,beta,y,desc_data,trans,sweeps,work,info)
  
  use psb_base_mod
  use mld_d_jac_smoother, mld_protect_name => mld_d_jac_smoother_apply_vect
  implicit none 
  type(psb_desc_type), intent(in)               :: desc_data
  class(mld_d_jac_smoother_type), intent(inout) :: sm
  type(psb_d_vect_type),intent(inout)           :: x
  type(psb_d_vect_type),intent(inout)           :: y
  real(psb_dpk_),intent(in)                     :: alpha,beta
  character(len=1),intent(in)                   :: trans
  integer, intent(in)                           :: sweeps
  real(psb_dpk_),target, intent(inout)          :: work(:)
  integer, intent(out)                          :: info

  integer    :: n_row,n_col
  type(psb_d_vect_type)  :: tx, ty
  real(psb_dpk_), pointer :: ww(:), aux(:)
  integer    :: ictxt,np,me,i, err_act
  character          :: trans_
  character(len=20)  :: name='d_jac_smoother_apply'

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

  if (.not.allocated(sm%sv)) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  end if

  n_row = desc_data%get_local_rows()
  n_col = desc_data%get_local_cols()

  if (n_col <= size(work)) then 
    ww => work(1:n_col)
    if ((4*n_col+n_col) <= size(work)) then 
      aux => work(n_col+1:)
    else
      allocate(aux(4*n_col),stat=info)
      if (info /= psb_success_) then 
        info=psb_err_alloc_request_
        call psb_errpush(info,name,i_err=(/4*n_col,0,0,0,0/),&
             & a_err='real(psb_dpk_)')
        goto 9999      
      end if
    endif
  else
    allocate(ww(n_col),aux(4*n_col),stat=info)
    if (info /= psb_success_) then 
      info=psb_err_alloc_request_
      call psb_errpush(info,name,i_err=(/5*n_col,0,0,0,0/),&
           & a_err='real(psb_dpk_)')
      goto 9999      
    end if
  endif

  if ((sweeps == 1).or.(sm%nnz_nd_tot==0)) then 

    call sm%sv%apply(alpha,x,beta,y,desc_data,trans_,aux,info) 

    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,&
           & name,a_err='Error in sub_aply Jacobi Sweeps = 1')
      goto 9999
    endif

  else if (sweeps  > 1) then 

    !
    !
    ! Apply multiple sweeps of a block-Jacobi solver
    ! to compute an approximate solution of a linear system.
    !
    !
    call tx%bld(x%get_nrows(),mold=x%v)
    call tx%set(dzero)
    call ty%bld(x%get_nrows(),mold=x%v)

    do i=1, sweeps
      !
      ! Compute Y(j+1) = D^(-1)*(X-ND*Y(j)), where D and ND are the
      ! block diagonal part and the remaining part of the local matrix
      ! and Y(j) is the approximate solution at sweep j.
      !
      call psb_geaxpby(done,x,dzero,ty,desc_data,info)
      call psb_spmm(-done,sm%nd,tx,done,ty,desc_data,info,work=aux,trans=trans_)

      if (info /= psb_success_) exit

      call sm%sv%apply(done,ty,dzero,tx,desc_data,trans_,aux,info) 

      if (info /= psb_success_) exit
    end do

    if (info == psb_success_) call psb_geaxpby(alpha,tx,beta,y,desc_data,info)

    if (info /= psb_success_) then 
      info=psb_err_internal_error_
      call psb_errpush(info,name,a_err='subsolve with Jacobi sweeps > 1')
      goto 9999      
    end if

    call tx%free(info) 
    if (info == psb_success_) call ty%free(info) 
    if (info /= psb_success_) then 
      info=psb_err_internal_error_
      call psb_errpush(info,name,a_err='final cleanup with Jacobi sweeps > 1')
      goto 9999      
    end if

  else

    info = psb_err_iarg_neg_
    call psb_errpush(info,name,&
         & i_err=(/2,sweeps,0,0,0/))
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
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_d_jac_smoother_apply_vect
