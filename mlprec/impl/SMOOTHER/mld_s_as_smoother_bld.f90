subroutine mld_s_as_smoother_bld(a,desc_a,sm,upd,info,amold,vmold)
  
  use psb_base_mod
  use mld_s_as_smoother, mld_protect_nam => mld_s_as_smoother_bld
  Implicit None

  ! Arguments
  type(psb_sspmat_type), intent(in), target          :: a
  Type(psb_desc_type), Intent(in)                    :: desc_a 
  class(mld_s_as_smoother_type), intent(inout)       :: sm
  character, intent(in)                              :: upd
  integer, intent(out)                               :: info
  class(psb_s_base_sparse_mat), intent(in), optional :: amold
  class(psb_s_base_vect_type), intent(in), optional  :: vmold

  ! Local variables
  type(psb_sspmat_type) :: blck, atmp
  integer :: n_row,n_col, nrow_a, nhalo, novr, data_, nzeros
  real(psb_spk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer :: ictxt,np,me,i, err_act, debug_unit, debug_level
  character(len=20)  :: name='s_as_smoother_bld', ch_err

  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ictxt       = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'


  novr   = sm%novr
  if (novr < 0) then
    info=psb_err_invalid_ovr_num_
    call psb_errpush(info,name,i_err=(/novr,0,0,0,0,0/))
    goto 9999
  endif

  if ((novr == 0).or.(np == 1)) then 
    if (psb_toupper(upd) == 'F') then 
      call psb_cdcpy(desc_a,sm%desc_data,info)
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
    call blck%csall(0,0,info,1)
  else

    If (psb_toupper(upd) == 'F') Then
      !
      ! Build the auxiliary descriptor desc_p%matrix_data(psb_n_row_).
      ! This is done by psb_cdbldext (interface to psb_cdovr), which is
      ! independent of CSR, and has been placed in the tools directory
      ! of PSBLAS, instead of the mlprec directory of MLD2P4, because it
      ! might be used independently of the AS preconditioner, to build
      ! a descriptor for an extended stencil in a PDE solver. 
      !
      call psb_cdbldext(a,desc_a,novr,sm%desc_data,info,extype=psb_ovt_asov_)
      if(debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & ' From cdbldext _:',sm%desc_data%get_local_rows(),&
           & sm%desc_data%get_local_cols()

      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='psb_cdbldext'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    Endif

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Before sphalo '

    !
    ! Retrieve the remote sparse matrix rows required for the AS extended
    ! matrix
    data_ = psb_comm_ext_
    Call psb_sphalo(a,sm%desc_data,blck,info,data=data_,rowscale=.true.)

    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_sphalo'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if (debug_level >=psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'After psb_sphalo ',&
         & blck%get_nrows(), blck%get_nzeros()

  End if
  if (info == psb_success_) &
       & call sm%sv%build(a,sm%desc_data,upd,info,&
       &  blck,amold=amold,vmold=vmold)

  nrow_a = a%get_nrows()
  n_row  = sm%desc_data%get_local_rows()
  n_col  = sm%desc_data%get_local_cols()

  if (info == psb_success_) call a%csclip(sm%nd,info,&
       & jmin=nrow_a+1,rscale=.false.,cscale=.false.)
  if (info == psb_success_) call blck%csclip(atmp,info,&
       & jmin=nrow_a+1,rscale=.false.,cscale=.false.)
  if (info == psb_success_) call psb_rwextd(n_row,sm%nd,info,b=atmp) 

  if (info == psb_success_) then 
    if (present(amold)) then 
      call sm%nd%cscnv(info,&
           & mold=amold,dupl=psb_dupl_add_)
    else
      call sm%nd%cscnv(info,&
           & type='csr',dupl=psb_dupl_add_)
    end if
  end if
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='clip & psb_spcnv csr 4')
    goto 9999
  end if
  nzeros = sm%nd%get_nzeros()
!!$    write(0,*) me,' ND nzeors ',nzeros
  call psb_sum(ictxt,nzeros)
  sm%nd_nnz_tot = nzeros

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
end subroutine mld_s_as_smoother_bld
