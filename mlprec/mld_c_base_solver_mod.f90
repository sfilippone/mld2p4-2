!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010
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
! File: mld_c_base_solver_mod.f90
!
! Module: mld_c_base_solver_mod
!
!  This module defines: 
!  - the mld_c_base_solver_type data structure containing the
!    basic solver type acting on a subdomain
!
!  It contains routines for
!  - Building and applying; 
!  - checking if the solver is correctly defined;
!  - printing a	description of the solver;
!  - deallocating the data structure.  
!

module mld_c_base_solver_mod

  use mld_base_prec_type
  use psb_base_mod, only : psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type
  !
  ! 
  ! Type: mld_T_base_solver_type.
  ! 
  !  It holds the local solver; it has no mandatory components. 
  !
  !  type  mld_T_base_solver_type
  !  end type mld_T_base_solver_type
  !
  !    build      -   Compute the actual contents of the smoother; includes
  !                   invocation of the build method on the solver component. 
  !    free       -   Release memory
  !    apply      -   Apply the smoother to a vector (or to an array); includes
  !                   invocation of the apply method on the solver component. 
  !    descr      -   Prints a description of the object.
  !    default    -   Set default values
  !    dump       -   Dump to file object contents
  !    set        -   Sets various parameters; when a request is unknown
  !                   it is passed to the smoother object for further processing.
  !    check      -   Sanity checks.
  !    sizeof     -   Total memory occupation in bytes
  !    get_nzeros -   Number of nonzeros 
  !
  !
  !

  type mld_c_base_solver_type
  contains
    procedure, pass(sv) :: check => c_base_solver_check
    procedure, pass(sv) :: dump  => c_base_solver_dmp
    procedure, pass(sv) :: build => c_base_solver_bld
    procedure, pass(sv) :: apply_v => c_base_solver_apply_vect
    procedure, pass(sv) :: apply_a => c_base_solver_apply
    generic, public     :: apply => apply_a, apply_v
    procedure, pass(sv) :: free  => c_base_solver_free
    procedure, pass(sv) :: seti  => c_base_solver_seti
    procedure, pass(sv) :: setc  => c_base_solver_setc
    procedure, pass(sv) :: setr  => c_base_solver_setr
    generic, public     :: set   => seti, setc, setr
    procedure, pass(sv) :: default => c_base_solver_default
    procedure, pass(sv) :: descr   => c_base_solver_descr
    procedure, pass(sv) :: sizeof  => c_base_solver_sizeof
    procedure, pass(sv) :: get_nzeros => c_base_solver_get_nzeros
  end type mld_c_base_solver_type

  private :: c_base_solver_bld,  c_base_solver_apply, &
       &  c_base_solver_free,    c_base_solver_seti, &
       &  c_base_solver_setc,    c_base_solver_setr, &
       &  c_base_solver_descr,   c_base_solver_sizeof, &
       &  c_base_solver_default, c_base_solver_check,&
       &  c_base_solver_dmp, c_base_solver_apply_vect, &
       &  c_base_solver_get_nzeros



contains
  !
  ! Function returning the size of the data structure
  ! in bytes or in number of nonzeros of the operator(s) involved. 
  !

  function c_base_solver_sizeof(sv) result(val)
    implicit none 
    ! Arguments
    class(mld_c_base_solver_type), intent(in) :: sv
    integer(psb_long_int_k_)                  :: val
    integer             :: i
    val = 0

    return
  end function c_base_solver_sizeof

  function c_base_solver_get_nzeros(sv) result(val)
    implicit none 
    class(mld_c_base_solver_type), intent(in) :: sv
    integer(psb_long_int_k_) :: val
    integer             :: i
    val = 0
  end function c_base_solver_get_nzeros

  
  !
  ! Apply: comes in two versions, on plain arrays or on encapsulated
  ! vectors.
  ! The base version throws an error, since it means that no explicit
  ! choice was made. 
  ! Question: would it make sense to transform the base version into
  ! the ID version, i.e. "base_solver" is the identity operator? 
  ! 

  subroutine c_base_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info)
    use psb_base_mod
    type(psb_desc_type), intent(in)           :: desc_data
    class(mld_c_base_solver_type), intent(in) :: sv
    complex(psb_spk_),intent(inout)              :: x(:)
    complex(psb_spk_),intent(inout)              :: y(:)
    complex(psb_spk_),intent(in)                 :: alpha,beta
    character(len=1),intent(in)               :: trans
    complex(psb_spk_),target, intent(inout)      :: work(:)
    integer, intent(out)                      :: info
    
    Integer :: err_act
    character(len=20)  :: name='d_base_solver_apply'

    call psb_erractionsave(err_act)
    
    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine c_base_solver_apply

  subroutine c_base_solver_apply_vect(alpha,sv,x,beta,y,desc_data,trans,work,info)
    use psb_base_mod
    type(psb_desc_type), intent(in)              :: desc_data
    class(mld_c_base_solver_type), intent(inout) :: sv
    type(psb_c_vect_type),intent(inout)          :: x
    type(psb_c_vect_type),intent(inout)          :: y
    complex(psb_spk_),intent(in)                    :: alpha,beta
    character(len=1),intent(in)                  :: trans
    complex(psb_spk_),target, intent(inout)         :: work(:)
    integer, intent(out)                         :: info
    
    Integer :: err_act
    character(len=20)  :: name='d_base_solver_apply'

    call psb_erractionsave(err_act)
    
    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine c_base_solver_apply_vect


  !
  ! Build
  ! The base version throws an error, since it means that no explicit
  ! choice was made. 
  !
  subroutine c_base_solver_bld(a,desc_a,sv,upd,info,b,amold,vmold)

    use psb_base_mod

    Implicit None

    ! Arguments
    type(psb_cspmat_type), intent(in), target           :: a
    Type(psb_desc_type), Intent(in)                     :: desc_a 
    class(mld_c_base_solver_type), intent(inout)        :: sv
    character, intent(in)                               :: upd
    integer, intent(out)                                :: info
    type(psb_cspmat_type), intent(in), target, optional :: b
    class(psb_c_base_sparse_mat), intent(in), optional  :: amold
    class(psb_c_base_vect_type), intent(in), optional   :: vmold

    Integer :: err_act
    character(len=20)  :: name='d_base_solver_bld'

    call psb_erractionsave(err_act)

    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine c_base_solver_bld

  subroutine c_base_solver_check(sv,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_c_base_solver_type), intent(inout) :: sv
    integer, intent(out)                   :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_solver_check'

    call psb_erractionsave(err_act)
    info = psb_success_


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
  end subroutine c_base_solver_check

  !
  ! Set.
  ! The base version does nothing; the principle is that
  ! SET acts on what is known, and delegates what is unknown.
  ! Since we are at the bottom of the hierarchy, there's no one
  ! to delegate, so we do nothing. 
  !
  subroutine c_base_solver_seti(sv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_c_base_solver_type), intent(inout) :: sv 
    integer, intent(in)                          :: what 
    integer, intent(in)                          :: val
    integer, intent(out)                         :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_solver_seti'
    
    ! Correct action here is doing nothing. 
    info = 0
    
    return
  end subroutine c_base_solver_seti

  subroutine c_base_solver_setc(sv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_c_base_solver_type), intent(inout) :: sv
    integer, intent(in)                          :: what 
    character(len=*), intent(in)                 :: val
    integer, intent(out)                         :: info
    Integer           :: err_act, ival 
    character(len=20) :: name='d_base_solver_setc'

    call psb_erractionsave(err_act)

    info = psb_success_

    call mld_stringval(val,ival,info)
    if (info == psb_success_) call sv%set(what,ival,info)

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
  end subroutine c_base_solver_setc
  
  subroutine c_base_solver_setr(sv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_c_base_solver_type), intent(inout) :: sv 
    integer, intent(in)                          :: what 
    real(psb_spk_), intent(in)                   :: val
    integer, intent(out)                         :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_solver_setr'

    
    ! Correct action here is doing nothing. 
    info = 0
    
    return
  end subroutine c_base_solver_setr

  !
  ! Free
  ! The base version throws an error, since it means that no explicit
  ! choice was made. IS THIS CORRECT? I suspect it would be better
  ! to be silent here, to cover reallocation. 
  !
  subroutine c_base_solver_free(sv,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_c_base_solver_type), intent(inout) :: sv
    integer, intent(out)                         :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_solver_free'

    call psb_erractionsave(err_act)

    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine c_base_solver_free

  subroutine c_base_solver_descr(sv,info,iout,coarse)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_c_base_solver_type), intent(in) :: sv
    integer, intent(out)                      :: info
    integer, intent(in), optional             :: iout
    logical, intent(in), optional             :: coarse

    ! Local variables
    integer      :: err_act
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_c_base_solver_descr'
    integer      :: iout_


    call psb_erractionsave(err_act)

    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine c_base_solver_descr

  !
  ! Dump. For debugging purposes. 
  !
  subroutine c_base_solver_dmp(sv,ictxt,level,info,prefix,head,solver)
    use psb_base_mod
    implicit none 
    class(mld_c_base_solver_type), intent(in) :: sv
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

    ! At base level do nothing for the solver

  end subroutine c_base_solver_dmp

  subroutine c_base_solver_default(sv) 
    implicit none 
    ! Arguments
    class(mld_c_base_solver_type), intent(inout) :: sv
    ! Do nothing for base version

    return
  end subroutine c_base_solver_default



end module mld_c_base_solver_mod
