!!$
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
! Identity solver. Reference for nullprec. 
!
!

module mld_s_id_solver

  use mld_s_prec_type

  type, extends(mld_s_base_solver_type) :: mld_s_id_solver_type
  contains
    procedure, pass(sv) :: build => s_id_solver_bld
    procedure, pass(sv) :: apply => s_id_solver_apply
    procedure, pass(sv) :: free  => s_id_solver_free
    procedure, pass(sv) :: seti  => s_id_solver_seti
    procedure, pass(sv) :: setc  => s_id_solver_setc
    procedure, pass(sv) :: setr  => s_id_solver_setr
    procedure, pass(sv) :: descr => s_id_solver_descr
    procedure, pass(sv) :: sizeof => s_id_solver_sizeof
  end type mld_s_id_solver_type


  private :: s_id_solver_bld, s_id_solver_apply, &
       &  s_id_solver_free,   s_id_solver_seti, &
       &  s_id_solver_setc,   s_id_solver_setr,&
       &  s_id_solver_descr,  s_id_solver_sizeof


contains

  subroutine s_id_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info)
    use psb_base_mod
    type(psb_desc_type), intent(in)      :: desc_data
    class(mld_s_id_solver_type), intent(in) :: sv
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
    character(len=20)  :: name='s_id_solver_apply'

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

  end subroutine s_id_solver_apply

  subroutine s_id_solver_bld(a,desc_a,sv,upd,info,b)

    use psb_base_mod

    Implicit None

    ! Arguments
    type(psb_sspmat_type), intent(in), target  :: a
    Type(psb_desc_type), Intent(in)             :: desc_a 
    class(mld_s_id_solver_type), intent(inout) :: sv
    character, intent(in)                       :: upd
    integer, intent(out)                        :: info
    type(psb_sspmat_type), intent(in), target, optional  :: b
    ! Local variables
    integer :: n_row,n_col, nrow_a, nztota
    real(psb_spk_), pointer :: ww(:), aux(:), tx(:),ty(:)
    integer :: ictxt,np,me,i, err_act, debug_unit, debug_level
    character(len=20)  :: name='s_id_solver_bld', ch_err
    
    info=psb_success_
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    ictxt       = desc_a%get_context()
    call psb_info(ictxt, me, np)
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' start'


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
  end subroutine s_id_solver_bld


  subroutine s_id_solver_seti(sv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_s_id_solver_type), intent(inout) :: sv 
    integer, intent(in)                    :: what 
    integer, intent(in)                    :: val
    integer, intent(out)                   :: info
    Integer :: err_act
    character(len=20)  :: name='s_id_solver_seti'

    info = psb_success_

    return

  end subroutine s_id_solver_seti

  subroutine s_id_solver_setc(sv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_s_id_solver_type), intent(inout) :: sv
    integer, intent(in)                    :: what 
    character(len=*), intent(in)           :: val
    integer, intent(out)                   :: info
    Integer :: err_act, ival
    character(len=20)  :: name='s_id_solver_setc'

    info = psb_success_

    return
  end subroutine s_id_solver_setc
  
  subroutine s_id_solver_setr(sv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_s_id_solver_type), intent(inout) :: sv 
    integer, intent(in)                    :: what 
    real(psb_spk_), intent(in)             :: val
    integer, intent(out)                   :: info
    Integer :: err_act
    character(len=20)  :: name='s_id_solver_setr'

    info = psb_success_

    return

  end subroutine s_id_solver_setr

  subroutine s_id_solver_free(sv,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_s_id_solver_type), intent(inout) :: sv
    integer, intent(out)                       :: info
    Integer :: err_act
    character(len=20)  :: name='s_id_solver_free'

    info = psb_success_

    return
  end subroutine s_id_solver_free

  subroutine s_id_solver_descr(sv,info,iout,coarse)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_s_id_solver_type), intent(in) :: sv
    integer, intent(out)                      :: info
    integer, intent(in), optional             :: iout
    logical, intent(in), optional             :: coarse

    ! Local variables
    integer      :: err_act
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_s_id_solver_descr'
    integer :: iout_

    info = psb_success_
    if (present(iout)) then 
      iout_ = iout 
    else
      iout_ = 6
    endif
    
    write(iout_,*) '  Identity local solver '

    return

  end subroutine s_id_solver_descr

  function s_id_solver_sizeof(sv) result(val)
    use psb_base_mod
    implicit none 
    ! Arguments
    class(mld_s_id_solver_type), intent(in) :: sv
    integer(psb_long_int_k_) :: val
    integer             :: i

    val = 0

    return
  end function s_id_solver_sizeof

end module mld_s_id_solver
