!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010,2010,2012
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

  use mld_s_base_solver_mod

  type, extends(mld_s_base_solver_type) :: mld_s_id_solver_type
  contains
    procedure, pass(sv) :: build   => s_id_solver_bld
    procedure, pass(sv) :: apply_v => mld_s_id_solver_apply_vect
    procedure, pass(sv) :: apply_a => mld_s_id_solver_apply
    procedure, pass(sv) :: free    => s_id_solver_free
    procedure, pass(sv) :: seti    => s_id_solver_seti
    procedure, pass(sv) :: setc    => s_id_solver_setc
    procedure, pass(sv) :: setr    => s_id_solver_setr
    procedure, pass(sv) :: descr   => s_id_solver_descr
  end type mld_s_id_solver_type


  private :: s_id_solver_bld, &
       &  s_id_solver_free,   s_id_solver_seti, &
       &  s_id_solver_setc,   s_id_solver_setr,&
       &  s_id_solver_descr

  interface 
    subroutine mld_s_id_solver_apply_vect(alpha,sv,x,beta,y,desc_data,trans,work,info)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
           & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, & 
           & mld_s_id_solver_type, psb_ipk_
      type(psb_desc_type), intent(in)              :: desc_data
      class(mld_s_id_solver_type), intent(inout) :: sv
      type(psb_s_vect_type),intent(inout)        :: x
      type(psb_s_vect_type),intent(inout)        :: y
      real(psb_spk_),intent(in)                   :: alpha,beta
      character(len=1),intent(in)                  :: trans
      real(psb_spk_),target, intent(inout)        :: work(:)
      integer(psb_ipk_), intent(out)               :: info
    end subroutine mld_s_id_solver_apply_vect
  end interface
  
  interface 
    subroutine mld_s_id_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
           & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
           & mld_s_id_solver_type, psb_ipk_
      type(psb_desc_type), intent(in)      :: desc_data
      class(mld_s_id_solver_type), intent(in) :: sv
      real(psb_spk_),intent(inout)         :: x(:)
      real(psb_spk_),intent(inout)         :: y(:)
      real(psb_spk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)          :: trans
      real(psb_spk_),target, intent(inout) :: work(:)
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine mld_s_id_solver_apply
  end interface

contains


  subroutine s_id_solver_bld(a,desc_a,sv,upd,info,b,amold,vmold)

    Implicit None

    ! Arguments
    type(psb_sspmat_type), intent(in), target           :: a
    Type(psb_desc_type), Intent(in)                       :: desc_a 
    class(mld_s_id_solver_type), intent(inout)          :: sv
    character, intent(in)                                 :: upd
    integer(psb_ipk_), intent(out)                        :: info
    type(psb_sspmat_type), intent(in), target, optional :: b
    class(psb_s_base_sparse_mat), intent(in), optional  :: amold
    class(psb_s_base_vect_type), intent(in), optional   :: vmold
    ! Local variables
    integer(psb_ipk_) :: n_row,n_col, nrow_a, nztota
    real(psb_spk_), pointer :: ww(:), aux(:), tx(:),ty(:)
    integer(psb_ipk_) :: i, err_act, debug_unit, debug_level
    character(len=20)  :: name='s_id_solver_bld', ch_err
    
    info=psb_success_

    return
  end subroutine s_id_solver_bld


  subroutine s_id_solver_seti(sv,what,val,info)

    Implicit None

    ! Arguments
    class(mld_s_id_solver_type), intent(inout) :: sv 
    integer(psb_ipk_), intent(in)                :: what 
    integer(psb_ipk_), intent(in)                :: val
    integer(psb_ipk_), intent(out)               :: info
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='s_id_solver_seti'

    info = psb_success_

    return

  end subroutine s_id_solver_seti

  subroutine s_id_solver_setc(sv,what,val,info)

    Implicit None

    ! Arguments
    class(mld_s_id_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(in)                :: what 
    character(len=*), intent(in)                 :: val
    integer(psb_ipk_), intent(out)               :: info
    integer(psb_ipk_) :: err_act, ival
    character(len=20) :: name='s_id_solver_setc'

    info = psb_success_

    return
  end subroutine s_id_solver_setc
  
  subroutine s_id_solver_setr(sv,what,val,info)

    Implicit None

    ! Arguments
    class(mld_s_id_solver_type), intent(inout) :: sv 
    integer(psb_ipk_), intent(in)                :: what 
    real(psb_spk_), intent(in)                    :: val
    integer(psb_ipk_), intent(out)               :: info
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='s_id_solver_setr'

    info = psb_success_

    return

  end subroutine s_id_solver_setr

  subroutine s_id_solver_free(sv,info)

    Implicit None

    ! Arguments
    class(mld_s_id_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(out)               :: info
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='s_id_solver_free'

    info = psb_success_

    return
  end subroutine s_id_solver_free

  subroutine s_id_solver_descr(sv,info,iout,coarse)

    Implicit None

    ! Arguments
    class(mld_s_id_solver_type), intent(in) :: sv
    integer(psb_ipk_), intent(out)            :: info
    integer(psb_ipk_), intent(in), optional   :: iout
    logical, intent(in), optional       :: coarse

    ! Local variables
    integer(psb_ipk_)      :: err_act
    character(len=20), parameter :: name='mld_s_id_solver_descr'
    integer(psb_ipk_) :: iout_

    info = psb_success_
    if (present(iout)) then 
      iout_ = iout 
    else
      iout_ = 6
    endif
    
    write(iout_,*) '  Identity local solver '

    return

  end subroutine s_id_solver_descr

end module mld_s_id_solver
