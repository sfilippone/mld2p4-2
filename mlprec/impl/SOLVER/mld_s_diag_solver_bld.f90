subroutine mld_s_diag_solver_bld(a,desc_a,sv,upd,info,b,amold,vmold)
  
  use psb_base_mod
  use mld_s_diag_solver, mld_protect_name => mld_s_diag_solver_bld

  Implicit None

  ! Arguments
  type(psb_sspmat_type), intent(in), target           :: a
  Type(psb_desc_type), Intent(in)                     :: desc_a 
  class(mld_s_diag_solver_type), intent(inout)        :: sv
  character, intent(in)                               :: upd
  integer, intent(out)                                :: info
  type(psb_sspmat_type), intent(in), target, optional :: b
  class(psb_s_base_sparse_mat), intent(in), optional  :: amold
  class(psb_s_base_vect_type), intent(in), optional   :: vmold
  ! Local variables
  integer :: n_row,n_col, nrow_a, nztota
  real(psb_spk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer :: ictxt,np,me,i, err_act, debug_unit, debug_level
  character(len=20)  :: name='s_diag_solver_bld', ch_err

  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ictxt       = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'


  n_row  = desc_a%get_local_rows()
  nrow_a = a%get_nrows()
  if (allocated(sv%d)) then 
    if (size(sv%d) < n_row) then 
      deallocate(sv%d)
    endif
  endif
  if (.not.allocated(sv%d)) then 
    allocate(sv%d(n_row),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999      
    end if

  endif

  call a%get_diag(sv%d,info)
  if (present(b)) then 
    if (info == psb_success_) call b%get_diag(sv%d(nrow_a+1:), info)
  end if
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='get_diag')
    goto 9999      
  end if

  do i=1,n_row
    if (sv%d(i) == szero) then 
      sv%d(i) = sone
    else
      sv%d(i) = sone/sv%d(i)
    end if
  end do
  allocate(sv%dv,stat=info) 
  if (info == psb_success_) then 
    if (present(vmold)) then 
      allocate(sv%dv%v,mold=vmold,stat=info) 
    else
      allocate(psb_s_base_vect_type :: sv%dv%v,stat=info) 
    end if
  end if
  if (info == psb_success_) then 
    call sv%dv%bld(sv%d)
  else
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate sv%dv')
    goto 9999      
  end if

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' end'

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine mld_s_diag_solver_bld
