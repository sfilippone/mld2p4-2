subroutine mld_c_jac_smoother_bld(a,desc_a,sm,upd,info,amold,vmold)
  

  use psb_base_mod
  use mld_c_diag_solver
  use mld_c_jac_smoother, mld_protect_name => mld_c_jac_smoother_bld
  Implicit None

  ! Arguments
  type(psb_cspmat_type), intent(in), target          :: a
  Type(psb_desc_type), Intent(in)                    :: desc_a 
  class(mld_c_jac_smoother_type), intent(inout)      :: sm
  character, intent(in)                              :: upd
  integer, intent(out)                               :: info
  class(psb_c_base_sparse_mat), intent(in), optional :: amold
  class(psb_c_base_vect_type), intent(in), optional  :: vmold
  ! Local variables
  integer :: n_row,n_col, nrow_a, nztota, nzeros
  complex(psb_spk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer :: ictxt,np,me,i, err_act, debug_unit, debug_level
  character(len=20)  :: name='c_jac_smoother_bld', ch_err

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
  nztota = a%get_nzeros()
  select type (smsv => sm%sv)
  type is (mld_c_diag_solver_type)
    call a%clip_diag(sm%nd,info)
    class default
    call a%csclip(sm%nd,info,&
         & jmin=nrow_a+1,rscale=.false.,cscale=.false.)
  end select
  if (info == psb_success_) then 
    if (present(amold)) then 
      call sm%nd%cscnv(info,&
           & mold=amold,dupl=psb_dupl_add_)
    else
      call sm%nd%cscnv(info,&
           & type='csr',dupl=psb_dupl_add_)
    endif
  end if
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,&
         & a_err='clip & psb_spcnv csr 4')
    goto 9999
  end if

  call sm%sv%build(a,desc_a,upd,info,amold=amold,vmold=vmold)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,&
         & a_err='solver build')
    goto 9999
  end if
  nzeros = sm%nd%get_nzeros()
  call psb_sum(ictxt,nzeros)
  sm%nnz_nd_tot = nzeros
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

end subroutine mld_c_jac_smoother_bld
