  subroutine mld_s_id_solver_apply_vect(alpha,sv,x,beta,y,desc_data,trans,work,info)
    
    use psb_base_mod
    use mld_s_id_solver, mld_protect_name => mld_s_id_solver_apply_vect
    implicit none 
    type(psb_desc_type), intent(in)            :: desc_data
    class(mld_s_id_solver_type), intent(inout) :: sv
    type(psb_s_vect_type),intent(inout)        :: x
    type(psb_s_vect_type),intent(inout)        :: y
    real(psb_spk_),intent(in)                  :: alpha,beta
    character(len=1),intent(in)                :: trans
    real(psb_spk_),target, intent(inout)       :: work(:)
    integer, intent(out)                       :: info

    integer    :: n_row,n_col
    integer    :: ictxt,np,me,i, err_act
    character          :: trans_
    character(len=20)  :: name='s_id_solver_apply_vect'

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

    call psb_geaxpby(alpha,x,beta,y,desc_data,info)    

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine mld_s_id_solver_apply_vect
