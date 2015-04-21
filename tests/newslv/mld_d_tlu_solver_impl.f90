!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.3)
!!$  
!!$  (C) Copyright 2008, 2010, 2012, 2015
!!$
!!$                      Salvatore Filippone  University of Rome Tor Vergata
!!$                      Alfredo Buttari      CNRS-IRIT, Toulouse
!!$                      Pasqua D'Ambra       ICAR-CNR, Naples
!!$                      Daniela di Serafino  Second University of Naples
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the MLD2P4 group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$
!
!
!
!


subroutine mld_d_tlu_solver_apply_vect(alpha,sv,x,beta,y,desc_data,trans,work,info)
  use psb_base_mod
  use mld_d_tlu_solver, mld_protect_name => mld_d_tlu_solver_apply_vect
  implicit none 
  type(psb_desc_type), intent(in)             :: desc_data
  class(mld_d_tlu_solver_type), intent(inout) :: sv
  type(psb_d_vect_type),intent(inout)         :: x
  type(psb_d_vect_type),intent(inout)         :: y
  real(psb_dpk_),intent(in)                   :: alpha,beta
  character(len=1),intent(in)                 :: trans
  real(psb_dpk_),target, intent(inout)        :: work(:)
  integer, intent(out)                        :: info

  integer    :: n_row,n_col
  type(psb_d_vect_type) :: wv, wv1
  real(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer    :: ictxt,np,me,i, err_act
  character          :: trans_
  character(len=20)  :: name='d_tlu_solver_apply'

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
         & a_err='real(psb_dpk_)')
    goto 9999      
  end if

  call psb_geasb(wv,desc_data,info,mold=x%v,scratch=.true.) 
  call psb_geasb(wv1,desc_data,info,mold=x%v,scratch=.true.) 

  select case(trans_)
  case('N')
    call psb_spsm(done,sv%l,x,dzero,wv,desc_data,info,&
         & trans=trans_,scale='L',diag=sv%dv,choice=psb_none_,work=aux)

    if (info == psb_success_) call psb_spsm(alpha,sv%u,wv,beta,y,desc_data,info,&
         & trans=trans_,scale='U',choice=psb_none_, work=aux)

  case('T')
    call psb_spsm(done,sv%u,x,dzero,wv,desc_data,info,&
         & trans=trans_,scale='L',diag=sv%dv,choice=psb_none_,work=aux)
    if (info == psb_success_) call psb_spsm(alpha,sv%l,wv,beta,y,desc_data,info,&
         & trans=trans_,scale='U',choice=psb_none_,work=aux)

  case('C')

    call psb_spsm(done,sv%u,x,dzero,wv,desc_data,info,&
         & trans=trans_,scale='U',choice=psb_none_,work=aux)

    call wv1%mlt(done,sv%dv,wv,dzero,info,conjgx=trans_)

    if (info == psb_success_) call psb_spsm(alpha,sv%l,wv1,beta,y,desc_data,info,&
         & trans=trans_,scale='U',choice=psb_none_,work=aux)

  case default
    call psb_errpush(psb_err_internal_error_,name,a_err='Invalid TRANS in TLU subsolve')
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

9999 call psb_error_handler(err_act)

  return

end subroutine mld_d_tlu_solver_apply_vect


subroutine mld_d_tlu_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info)
  use psb_base_mod
  use mld_d_tlu_solver, mld_protect_name => mld_d_tlu_solver_apply
  implicit none 
  type(psb_desc_type), intent(in)      :: desc_data
  class(mld_d_tlu_solver_type), intent(inout) :: sv
  real(psb_dpk_),intent(inout)         :: x(:)
  real(psb_dpk_),intent(inout)         :: y(:)
  real(psb_dpk_),intent(in)            :: alpha,beta
  character(len=1),intent(in)          :: trans
  real(psb_dpk_),target, intent(inout) :: work(:)
  integer, intent(out)                 :: info

  integer    :: n_row,n_col
  real(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer    :: ictxt,np,me,i, err_act
  character          :: trans_
  character(len=20)  :: name='d_tlu_solver_apply'

  call psb_erractionsave(err_act)

  info = psb_success_

  trans_ = psb_toupper(trans)
  select case(trans_)
  case('N')
  case('T','C')
  case default
    call psb_errpush(psb_err_iarg_invalid_i_,name)
    goto 9999
  end select

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

  select case(trans_)
  case('N')
    call psb_spsm(done,sv%l,x,dzero,ww,desc_data,info,&
         & trans=trans_,scale='L',diag=sv%d,choice=psb_none_,work=aux)

    if (info == psb_success_) call psb_spsm(alpha,sv%u,ww,beta,y,desc_data,info,&
         & trans=trans_,scale='U',choice=psb_none_, work=aux)

  case('T')
    call psb_spsm(done,sv%u,x,dzero,ww,desc_data,info,&
         & trans=trans_,scale='L',diag=sv%d,choice=psb_none_,work=aux)
    if (info == psb_success_) call psb_spsm(alpha,sv%l,ww,beta,y,desc_data,info,&
         & trans=trans_,scale='U',choice=psb_none_,work=aux)
  case('C')
    call psb_spsm(done,sv%u,x,dzero,ww,desc_data,info,&
         & trans=trans_,scale='L',diag=sv%d,choice=psb_none_,work=aux)
    if (info == psb_success_) call psb_spsm(alpha,sv%l,ww,beta,y,desc_data,info,&
         & trans=trans_,scale='U',choice=psb_none_,work=aux)
  case default
    call psb_errpush(psb_err_internal_error_,name,a_err='Invalid TRANS in TLU subsolve')
    goto 9999
  end select


  if (info /= psb_success_) then

    call psb_errpush(psb_err_internal_error_,name,a_err='Error in subsolve')
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

9999 call psb_error_handler(err_act)

  return

end subroutine mld_d_tlu_solver_apply

subroutine mld_d_tlu_solver_bld(a,desc_a,sv,upd,info,b,amold,vmold)

  use psb_base_mod
  use mld_d_tlu_solver, mld_protect_name => mld_d_tlu_solver_bld

  Implicit None

  ! Arguments
  type(psb_dspmat_type), intent(in), target           :: a
  Type(psb_desc_type), Intent(in)                     :: desc_a 
  class(mld_d_tlu_solver_type), intent(inout)         :: sv
  character, intent(in)                               :: upd
  integer, intent(out)                                :: info
  type(psb_dspmat_type), intent(in), target, optional :: b
  class(psb_d_base_sparse_mat), intent(in), optional  :: amold
  class(psb_d_base_vect_type), intent(in), optional   :: vmold
  ! Local variables
  integer :: n_row,n_col, nrow_a, nztota
!!$    real(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer :: ictxt,np,me,i, err_act, debug_unit, debug_level
  character(len=20)  :: name='d_tlu_solver_bld', ch_err

  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ictxt       = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'


  n_row  = desc_a%get_local_rows()

  if (psb_toupper(upd) == 'F') then 
    nrow_a = a%get_nrows()
    nztota = a%get_nzeros()
    if (present(b)) then 
      nztota = nztota + b%get_nzeros()
    end if

    call sv%l%csall(n_row,n_row,info,nztota)
    if (info == psb_success_) call sv%u%csall(n_row,n_row,info,nztota)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if (allocated(sv%d)) then 
      if (size(sv%d) < n_row) then 
        deallocate(sv%d)
      endif
    endif
    if (.not.allocated(sv%d))  allocate(sv%d(n_row),stat=info)

    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999      
    endif


    select case(sv%fact_type)

    case (mld_ilu_t_)
      !
      ! ILU(k,t)
      !
      select case(sv%fill_in)

      case(:-1) 
        ! Error: fill-in <= -1
        call psb_errpush(psb_err_input_value_invalid_i_,&
             & name,i_err=(/3,sv%fill_in,0,0,0/))
        goto 9999

      case(0:)
        ! Fill-in >= 0
        call mld_ilut_fact(sv%fill_in,sv%thresh,&
             & a, sv%l,sv%u,sv%d,info,blck=b)
      end select
      if(info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='mld_ilut_fact'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

    case(mld_ilu_n_,mld_milu_n_) 
      !
      ! ILU(k) and MILU(k)
      !
      select case(sv%fill_in)
      case(:-1) 
        ! Error: fill-in <= -1
        call psb_errpush(psb_err_input_value_invalid_i_,&
             & name,i_err=(/3,sv%fill_in,0,0,0/))
        goto 9999
      case(0)
        ! Fill-in 0
        ! Separate implementation of ILU(0) for better performance.
        ! There seems to be a problem with the separate implementation of MILU(0),
        ! contained into mld_ilu0_fact. This must be investigated. For the time being,
        ! resort to the implementation of MILU(k) with k=0.
        if (sv%fact_type == mld_ilu_n_) then 
          call mld_ilu0_fact(sv%fact_type,a,sv%l,sv%u,&
               & sv%d,info,blck=b,upd=upd)
        else
          call mld_iluk_fact(sv%fill_in,sv%fact_type,&
               & a,sv%l,sv%u,sv%d,info,blck=b)
        endif
      case(1:)
        ! Fill-in >= 1
        ! The same routine implements both ILU(k) and MILU(k)
        call mld_iluk_fact(sv%fill_in,sv%fact_type,&
             & a,sv%l,sv%u,sv%d,info,blck=b)
      end select
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='mld_iluk_fact'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

    case default
      ! If we end up here, something was wrong up in the call chain. 
      info = psb_err_input_value_invalid_i_
      call psb_errpush(psb_err_input_value_invalid_i_,name,&
           & i_err=(/3,sv%fact_type,0,0,0/))
      goto 9999

    end select
  else
    ! Here we should add checks for reuse of L and U.
    ! For the time being just throw an error. 
    info = 31
    call psb_errpush(info, name,&
         & i_err=(/3,0,0,0,0/),a_err=upd)
    goto 9999 

    !
    ! What is an update of a factorization??
    ! A first attempt could be to reuse EXACTLY the existing indices
    ! as if it was an ILU(0) (since, effectively, the sparsity pattern
    ! should not grow beyond what is already there).
    !  
    call mld_ilu0_fact(sv%fact_type,a,&
         & sv%l,sv%u,&
         & sv%d,info,blck=b,upd=upd)

  end if

  call sv%l%set_asb()
  call sv%l%trim()
  call sv%u%set_asb()
  call sv%u%trim()
  call sv%dv%bld(sv%d,mold=vmold)

  if (present(amold)) then 
    call sv%l%cscnv(info,mold=amold)
    call sv%u%cscnv(info,mold=amold)
  end if

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end subroutine mld_d_tlu_solver_bld


subroutine mld_d_tlu_solver_dmp(sv,ictxt,level,info,prefix,head,solver)
  use psb_base_mod
  use mld_d_tlu_solver, mld_protect_name => mld_d_tlu_solver_dmp
  implicit none 
  class(mld_d_tlu_solver_type), intent(in) :: sv
  integer, intent(in)              :: ictxt,level
  integer, intent(out)             :: info
  character(len=*), intent(in), optional :: prefix, head
  logical, optional, intent(in)    :: solver
  integer :: i, j, il1, iln, lname, lev
  integer :: icontxt,iam, np
  character(len=80)  :: prefix_
  character(len=120) :: fname ! len should be at least 20 more than
  logical :: solver_
  !  len of prefix_ 

  info = 0


  call psb_info(ictxt,iam,np)

  if (present(solver)) then 
    solver_ = solver
  else
    solver_ = .false. 
  end if

  if (solver_) then 
    if (present(prefix)) then 
      prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
    else
      prefix_ = "dump_slv_d"
    end if
    lname = len_trim(prefix_)
    fname = trim(prefix_)
    write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
    lname = lname + 5

    write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_lower.mtx'
    if (sv%l%is_asb()) &
         & call sv%l%print(fname,head=head)
    write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_diag.mtx'
    if (allocated(sv%d)) &
         & call psb_geprt(fname,sv%d,head=head)
    write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_upper.mtx'
    if (sv%u%is_asb()) &
         & call sv%u%print(fname,head=head)

  end if

end subroutine mld_d_tlu_solver_dmp

