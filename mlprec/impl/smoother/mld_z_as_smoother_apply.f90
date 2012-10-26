subroutine mld_z_as_smoother_apply(alpha,sm,x,beta,y,desc_data,trans,sweeps,work,info)
  
  use psb_base_mod
  use mld_z_as_smoother, mld_protect_nam => mld_z_as_smoother_apply
  implicit none 
  type(psb_desc_type), intent(in)      :: desc_data
  class(mld_z_as_smoother_type), intent(in) :: sm
  complex(psb_dpk_),intent(inout)         :: x(:)
  complex(psb_dpk_),intent(inout)         :: y(:)
  complex(psb_dpk_),intent(in)            :: alpha,beta
  character(len=1),intent(in)          :: trans
  integer, intent(in)                  :: sweeps
  complex(psb_dpk_),target, intent(inout) :: work(:)
  integer, intent(out)                 :: info

  integer    :: n_row,n_col, nrow_d, i
  complex(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer    :: ictxt,np,me, err_act,isz,int_err(5)
  character          :: trans_
  character(len=20)  :: name='z_as_smoother_apply', ch_err

  call psb_erractionsave(err_act)

  info = psb_success_
  ictxt = desc_data%get_context()
  call psb_info (ictxt,me,np)

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


  n_row = sm%desc_data%get_local_rows()
  n_col = sm%desc_data%get_local_cols()
  nrow_d = desc_data%get_local_rows()
  isz=max(n_row,N_COL)
  if ((6*isz) <= size(work)) then 
    ww => work(1:isz)
    tx => work(isz+1:2*isz)
    ty => work(2*isz+1:3*isz)
    aux => work(3*isz+1:)
  else if ((4*isz) <= size(work)) then 
    aux => work(1:)
    allocate(ww(isz),tx(isz),ty(isz),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_alloc_request_,name,i_err=(/3*isz,0,0,0,0/),&
           & a_err='complex(psb_dpk_)')
      goto 9999      
    end if
  else if ((3*isz) <= size(work)) then 
    ww => work(1:isz)
    tx => work(isz+1:2*isz)
    ty => work(2*isz+1:3*isz)
    allocate(aux(4*isz),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_alloc_request_,name,i_err=(/4*isz,0,0,0,0/),&
           & a_err='complex(psb_dpk_)')
      goto 9999      
    end if
  else 
    allocate(ww(isz),tx(isz),ty(isz),&
         &aux(4*isz),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_alloc_request_,name,i_err=(/4*isz,0,0,0,0/),&
           & a_err='complex(psb_dpk_)')
      goto 9999      
    end if

  endif

  if ((sm%novr == 0).and.(sweeps == 1)) then 
    !
    ! Shortcut: in this case it's just the same
    ! as Block Jacobi.
    !
    call sm%sv%apply(alpha,x,beta,y,desc_data,trans_,aux,info) 

    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Error in sub_aply Jacobi Sweeps = 1')
      goto 9999
    endif

  else 


    tx(1:nrow_d)     = x(1:nrow_d) 
    tx(nrow_d+1:isz) = zzero


    if (sweeps == 1) then 

      select case(trans_)
      case('N')
        !
        ! Get the overlap entries of tx (tx == x)
        ! 
        if (sm%restr == psb_halo_) then 
          call psb_halo(tx,sm%desc_data,info,work=aux,data=psb_comm_ext_)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='psb_halo'
            goto 9999
          end if
        else if (sm%restr /= psb_none_) then 
          call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_restr_')
          goto 9999
        end if


      case('T','C')
        !
        ! With transpose, we have to do it here
        ! 

        select case (sm%prol) 

        case(psb_none_)
          ! 
          ! Do nothing

        case(psb_sum_) 
          !
          ! The transpose of sum is halo
          !
          call psb_halo(tx,sm%desc_data,info,work=aux,data=psb_comm_ext_)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='psb_halo'
            goto 9999
          end if

        case(psb_avg_) 
          !
          ! Tricky one: first we have to scale the overlap entries,
          ! which we can do by assignind mode=0, i.e. no communication
          ! (hence only scaling), then we do the halo
          !
          call psb_ovrl(tx,sm%desc_data,info,&
               & update=psb_avg_,work=aux,mode=0)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='psb_ovrl'
            goto 9999
          end if
          call psb_halo(tx,sm%desc_data,info,work=aux,data=psb_comm_ext_)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='psb_halo'
            goto 9999
          end if

        case default
          call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_prol_')
          goto 9999
        end select


      case default
        info=psb_err_iarg_invalid_i_
        int_err(1)=6
        ch_err(2:2)=trans
        goto 9999
      end select

      call sm%sv%apply(zone,tx,zzero,ty,sm%desc_data,trans_,aux,info) 

      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error in sub_aply Jacobi Sweeps = 1')
        goto 9999
      endif

      select case(trans_)
      case('N')

        select case (sm%prol) 

        case(psb_none_)
          ! 
          ! Would work anyway, but since it is supposed to do nothing ...
          !        call psb_ovrl(ty,sm%desc_data,info,&
          !             & update=sm%prol,work=aux)


        case(psb_sum_,psb_avg_) 
          !
          ! Update the overlap of ty
          !
          call psb_ovrl(ty,sm%desc_data,info,&
               & update=sm%prol,work=aux)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='psb_ovrl'
            goto 9999
          end if

        case default
          call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_prol_')
          goto 9999
        end select

      case('T','C')
        !
        ! With transpose, we have to do it here
        ! 
        if (sm%restr == psb_halo_) then 
          call psb_ovrl(ty,sm%desc_data,info,&
               & update=psb_sum_,work=aux)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='psb_ovrl'
            goto 9999
          end if
        else if (sm%restr /= psb_none_) then 
          call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_restr_')
          goto 9999
        end if

      case default
        info=psb_err_iarg_invalid_i_
        int_err(1)=6
        ch_err(2:2)=trans
        goto 9999
      end select



    else if (sweeps  > 1) then 

      !
      !
      ! Apply prec%iprcparm(mld_smoother_sweeps_) sweeps of a block-Jacobi solver
      ! to compute an approximate solution of a linear system.
      !
      !
      ty = zzero
      do i=1, sweeps
        select case(trans_)
        case('N')
          !
          ! Get the overlap entries of tx (tx == x)
          ! 
          if (sm%restr == psb_halo_) then 
            call psb_halo(tx,sm%desc_data,info,work=aux,data=psb_comm_ext_)
            if(info /= psb_success_) then
              info=psb_err_from_subroutine_
              ch_err='psb_halo'
              goto 9999
            end if
          else if (sm%restr /= psb_none_) then 
            call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_restr_')
            goto 9999
          end if


        case('T','C')
          !
          ! With transpose, we have to do it here
          ! 

          select case (sm%prol) 

          case(psb_none_)
            ! 
            ! Do nothing

          case(psb_sum_) 
            !
            ! The transpose of sum is halo
            !
            call psb_halo(tx,sm%desc_data,info,work=aux,data=psb_comm_ext_)
            if(info /= psb_success_) then
              info=psb_err_from_subroutine_
              ch_err='psb_halo'
              goto 9999
            end if

          case(psb_avg_) 
            !
            ! Tricky one: first we have to scale the overlap entries,
            ! which we can do by assignind mode=0, i.e. no communication
            ! (hence only scaling), then we do the halo
            !
            call psb_ovrl(tx,sm%desc_data,info,&
                 & update=psb_avg_,work=aux,mode=0)
            if(info /= psb_success_) then
              info=psb_err_from_subroutine_
              ch_err='psb_ovrl'
              goto 9999
            end if
            call psb_halo(tx,sm%desc_data,info,work=aux,data=psb_comm_ext_)
            if(info /= psb_success_) then
              info=psb_err_from_subroutine_
              ch_err='psb_halo'
              goto 9999
            end if

          case default
            call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_prol_')
            goto 9999
          end select


        case default
          info=psb_err_iarg_invalid_i_
          int_err(1)=6
          ch_err(2:2)=trans
          goto 9999
        end select
        !
        ! Compute Y(j+1) = D^(-1)*(X-ND*Y(j)), where D and ND are the
        ! block diagonal part and the remaining part of the local matrix
        ! and Y(j) is the approximate solution at sweep j.
        !
        ww(1:n_row) = tx(1:n_row)
        call psb_spmm(-zone,sm%nd,ty,zone,ww,sm%desc_data,info,work=aux,trans=trans_)

        if (info /= psb_success_) exit

        call sm%sv%apply(zone,ww,zzero,ty,sm%desc_data,trans_,aux,info) 

        if (info /= psb_success_) exit


        select case(trans_)
        case('N')

          select case (sm%prol) 

          case(psb_none_)
            ! 
            ! Would work anyway, but since it is supposed to do nothing ...
            !        call psb_ovrl(ty,sm%desc_data,info,&
            !             & update=sm%prol,work=aux)


          case(psb_sum_,psb_avg_) 
            !
            ! Update the overlap of ty
            !
            call psb_ovrl(ty,sm%desc_data,info,&
                 & update=sm%prol,work=aux)
            if(info /= psb_success_) then
              info=psb_err_from_subroutine_
              ch_err='psb_ovrl'
              goto 9999
            end if

          case default
            call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_prol_')
            goto 9999
          end select

        case('T','C')
          !
          ! With transpose, we have to do it here
          ! 
          if (sm%restr == psb_halo_) then 
            call psb_ovrl(ty,sm%desc_data,info,&
                 & update=psb_sum_,work=aux)
            if(info /= psb_success_) then
              info=psb_err_from_subroutine_
              ch_err='psb_ovrl'
              goto 9999
            end if
          else if (sm%restr /= psb_none_) then 
            call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_restr_')
            goto 9999
          end if

        case default
          info=psb_err_iarg_invalid_i_
          int_err(1)=6
          ch_err(2:2)=trans
          goto 9999
        end select
      end do

      if (info /= psb_success_) then 
        info=psb_err_internal_error_
        call psb_errpush(info,name,a_err='subsolve with Jacobi sweeps > 1')
        goto 9999      
      end if


    else

      info = psb_err_iarg_neg_
      call psb_errpush(info,name,&
           & i_err=(/2,sweeps,0,0,0/))
      goto 9999


    end if

    !
    ! Compute y = beta*y + alpha*ty (ty == K^(-1)*tx)
    !
    call psb_geaxpby(alpha,ty,beta,y,desc_data,info) 

  end if


  if ((6*isz) <= size(work)) then 
  else if ((4*isz) <= size(work)) then 
    deallocate(ww,tx,ty)
  else if ((3*isz) <= size(work)) then 
    deallocate(aux)
  else 
    deallocate(ww,aux,tx,ty)
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

end subroutine mld_z_as_smoother_apply
