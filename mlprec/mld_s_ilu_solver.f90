!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010, 2010
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
!
!

module mld_s_ilu_solver

  use mld_s_prec_type
  use mld_s_ilu_fact_mod

  type, extends(mld_s_base_solver_type) :: mld_s_ilu_solver_type
    type(psb_sspmat_type)       :: l, u
    real(psb_spk_), allocatable :: d(:)
    integer                     :: fact_type, fill_in
    real(psb_spk_)              :: thresh
  contains
    procedure, pass(sv) :: dump  => s_ilu_solver_dmp
    procedure, pass(sv) :: build => s_ilu_solver_bld
    procedure, pass(sv) :: apply => s_ilu_solver_apply
    procedure, pass(sv) :: free  => s_ilu_solver_free
    procedure, pass(sv) :: seti  => s_ilu_solver_seti
    procedure, pass(sv) :: setc  => s_ilu_solver_setc
    procedure, pass(sv) :: setr  => s_ilu_solver_setr
    procedure, pass(sv) :: descr => s_ilu_solver_descr
    procedure, pass(sv) :: sizeof => s_ilu_solver_sizeof
    procedure, pass(sv) :: default => s_ilu_solver_default
  end type mld_s_ilu_solver_type


  private :: s_ilu_solver_bld, s_ilu_solver_apply, &
       &  s_ilu_solver_free,   s_ilu_solver_seti, &
       &  s_ilu_solver_setc,   s_ilu_solver_setr,&
       &  s_ilu_solver_descr,  s_ilu_solver_sizeof, &
       &  s_ilu_solver_default, s_ilu_solver_dmp

  character(len=15), parameter, private :: &
       &  fact_names(0:mld_slv_delta_+4)=(/&
       &  'none          ','none          ',&
       &  'none          ','none          ',&
       &  'none          ','DIAG ??       ',&
       &  'ILU(n)        ',&
       &  'MILU(n)       ','ILU(t,n)      '/)


contains

  subroutine s_ilu_solver_default(sv)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_s_ilu_solver_type), intent(inout) :: sv

    sv%fact_type = mld_ilu_n_
    sv%fill_in   = 0
    sv%thresh    = szero

    return
  end subroutine s_ilu_solver_default

  subroutine s_ilu_solver_check(sv,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_s_ilu_solver_type), intent(inout) :: sv
    integer, intent(out)                   :: info
    Integer           :: err_act
    character(len=20) :: name='s_ilu_solver_check'

    call psb_erractionsave(err_act)
    info = psb_success_

    call mld_check_def(sv%fact_type,&
         & 'Factorization',mld_ilu_n_,is_legal_ilu_fact)

    select case(sv%fact_type)
    case(mld_ilu_n_,mld_milu_n_)      
      call mld_check_def(sv%fill_in,&
           & 'Level',0,is_legal_ml_lev)
    case(mld_ilu_t_)                 
      call mld_check_def(sv%thresh,&
           & 'Eps',szero,is_legal_s_fact_thrs)
    end select
    
    if (info /= psb_success_) goto 9999
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine s_ilu_solver_check


  subroutine s_ilu_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info)
    use psb_base_mod
    type(psb_desc_type), intent(in)      :: desc_data
    class(mld_s_ilu_solver_type), intent(in) :: sv
    real(psb_spk_),intent(inout)         :: x(:)
    real(psb_spk_),intent(inout)         :: y(:)
    real(psb_spk_),intent(in)            :: alpha,beta
    character(len=1),intent(in)          :: trans
    real(psb_spk_),target, intent(inout) :: work(:)
    integer, intent(out)                 :: info

    integer    :: n_row,n_col
    real(psb_spk_), pointer :: ww(:), aux(:), tx(:),ty(:)
    integer    :: ictxt,np,me,i, err_act
    character          :: trans_
    character(len=20)  :: name='s_ilu_solver_apply'

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
               & a_err='real(psb_spk_)')
          goto 9999      
        end if
      endif
    else
      allocate(ww(n_col),aux(4*n_col),stat=info)
      if (info /= psb_success_) then 
        info=psb_err_alloc_request_
        call psb_errpush(info,name,i_err=(/5*n_col,0,0,0,0/),&
             & a_err='real(psb_spk_)')
        goto 9999      
      end if
    endif

    select case(trans_)
    case('N')
      call psb_spsm(sone,sv%l,x,szero,ww,desc_data,info,&
           & trans=trans_,scale='L',diag=sv%d,choice=psb_none_,work=aux)

      if (info == psb_success_) call psb_spsm(alpha,sv%u,ww,beta,y,desc_data,info,&
           & trans=trans_,scale='U',choice=psb_none_, work=aux)

    case('T','C')
      call psb_spsm(sone,sv%u,x,szero,ww,desc_data,info,&
           & trans=trans_,scale='L',diag=sv%d,choice=psb_none_,work=aux)
      if (info == psb_success_) call psb_spsm(alpha,sv%l,ww,beta,y,desc_data,info,&
           & trans=trans_,scale='U',choice=psb_none_,work=aux)
    case default
      call psb_errpush(psb_err_internal_error_,name,a_err='Invalid TRANS in ILU subsolve')
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

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine s_ilu_solver_apply

  subroutine s_ilu_solver_bld(a,desc_a,sv,upd,info,b)

    use psb_base_mod

    Implicit None

    ! Arguments
    type(psb_sspmat_type), intent(in), target  :: a
    Type(psb_desc_type), Intent(in)             :: desc_a 
    class(mld_s_ilu_solver_type), intent(inout) :: sv
    character, intent(in)                       :: upd
    integer, intent(out)                        :: info
    type(psb_sspmat_type), intent(in), target, optional  :: b
    ! Local variables
    integer :: n_row,n_col, nrow_a, nztota
    real(psb_spk_), pointer :: ww(:), aux(:), tx(:),ty(:)
    integer :: ictxt,np,me,i, err_act, debug_unit, debug_level
    character(len=20)  :: name='s_ilu_solver_bld', ch_err
    
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
      if (.not.allocated(sv%d)) then 
        allocate(sv%d(n_row),stat=info)
        if (info /= psb_success_) then 
          call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
          goto 9999      
        end if

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
      call psb_errpush(info, name, i_err=(/3,0,0,0,0/),a_err=upd)
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
  end subroutine s_ilu_solver_bld


  subroutine s_ilu_solver_seti(sv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_s_ilu_solver_type), intent(inout) :: sv 
    integer, intent(in)                    :: what 
    integer, intent(in)                    :: val
    integer, intent(out)                   :: info
    Integer :: err_act
    character(len=20)  :: name='s_ilu_solver_seti'

    info = psb_success_
    call psb_erractionsave(err_act)

    select case(what) 
    case(mld_sub_solve_) 
      sv%fact_type = val
    case(mld_sub_fillin_)
      sv%fill_in   = val
    case default
!!$      write(0,*) name,': Error: invalid WHAT'
!!$      info = -2
    end select

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine s_ilu_solver_seti

  subroutine s_ilu_solver_setc(sv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_s_ilu_solver_type), intent(inout) :: sv
    integer, intent(in)                    :: what 
    character(len=*), intent(in)           :: val
    integer, intent(out)                   :: info
    Integer :: err_act, ival
    character(len=20)  :: name='s_ilu_solver_setc'

    info = psb_success_
    call psb_erractionsave(err_act)


    call mld_stringval(val,ival,info)
    if (info == psb_success_) call sv%set(what,ival,info)
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info, name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine s_ilu_solver_setc
  
  subroutine s_ilu_solver_setr(sv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_s_ilu_solver_type), intent(inout) :: sv 
    integer, intent(in)                    :: what 
    real(psb_spk_), intent(in)             :: val
    integer, intent(out)                   :: info
    Integer :: err_act
    character(len=20)  :: name='s_ilu_solver_setr'

    call psb_erractionsave(err_act)
    info = psb_success_

    select case(what)
    case(mld_sub_iluthrs_) 
      sv%thresh = val
    case default
!!$      write(0,*) name,': Error: invalid WHAT'
!!$      info = -2
!!$      goto 9999
    end select

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine s_ilu_solver_setr

  subroutine s_ilu_solver_free(sv,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_s_ilu_solver_type), intent(inout) :: sv
    integer, intent(out)                       :: info
    Integer :: err_act
    character(len=20)  :: name='s_ilu_solver_free'

    call psb_erractionsave(err_act)
    info = psb_success_
    
    if (allocated(sv%d)) then 
      deallocate(sv%d,stat=info)
      if (info /= psb_success_) then 
        info = psb_err_alloc_dealloc_
        call psb_errpush(info,name)
        goto 9999 
      end if
    end if
    call sv%l%free()
    call sv%u%free()

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine s_ilu_solver_free

  subroutine s_ilu_solver_descr(sv,info,iout,coarse)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_s_ilu_solver_type), intent(in) :: sv
    integer, intent(out)                     :: info
    integer, intent(in), optional            :: iout
    logical, intent(in), optional            :: coarse

    ! Local variables
    integer      :: err_act
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_s_ilu_solver_descr'
    integer :: iout_

    call psb_erractionsave(err_act)
    info = psb_success_
    if (present(iout)) then 
      iout_ = iout 
    else
      iout_ = 6
    endif
    
    write(iout_,*) '  Incomplete factorization solver: ',&
           &  fact_names(sv%fact_type)
    select case(sv%fact_type)
    case(mld_ilu_n_,mld_milu_n_)      
      write(iout_,*) '  Fill level:',sv%fill_in
    case(mld_ilu_t_)         
      write(iout_,*) '  Fill level:',sv%fill_in
      write(iout_,*) '  Fill threshold :',sv%thresh
    end select

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine s_ilu_solver_descr

  function s_ilu_solver_sizeof(sv) result(val)
    use psb_base_mod
    implicit none 
    ! Arguments
    class(mld_s_ilu_solver_type), intent(in) :: sv
    integer(psb_long_int_k_) :: val
    integer             :: i

    val = 2*psb_sizeof_int + psb_sizeof_dp
    if (allocated(sv%d)) val = val + psb_sizeof_dp * size(sv%d)
    val = val + psb_sizeof(sv%l)
    val = val + psb_sizeof(sv%u)

    return
  end function s_ilu_solver_sizeof

  subroutine s_ilu_solver_dmp(sv,ictxt,level,info,prefix,head,solver)
    use psb_base_mod
    implicit none 
    class(mld_s_ilu_solver_type), intent(in) :: sv
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

    if (present(prefix)) then 
      prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
    else
      prefix_ = "dump_slv_d"
    end if

    call psb_info(ictxt,iam,np)

    if (present(solver)) then 
      solver_ = solver
    else
      solver_ = .false. 
    end if
    lname = len_trim(prefix_)
    fname = trim(prefix_)
    write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
    lname = lname + 5

    if (solver_) then 
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

  end subroutine s_ilu_solver_dmp


end module mld_s_ilu_solver
