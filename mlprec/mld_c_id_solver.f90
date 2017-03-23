!  
!   
!                             MLD2P4  version 2.0
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.3)
!    
!    (C) Copyright 2008, 2010, 2012, 2015
!  
!                        Salvatore Filippone  University of Rome Tor Vergata
!                        Alfredo Buttari      CNRS-IRIT, Toulouse
!                        Pasqua D'Ambra       ICAR-CNR, Naples
!                        Daniela di Serafino  Second University of Naples
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the MLD2P4 group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!  
!
!
!
! Identity solver. Reference for nullprec. 
!
!
module mld_c_id_solver

  use mld_c_base_solver_mod

  type, extends(mld_c_base_solver_type) :: mld_c_id_solver_type
  contains
    procedure, pass(sv) :: build   => c_id_solver_bld
    procedure, pass(sv) :: clone   => mld_c_id_solver_clone
    procedure, pass(sv) :: apply_v => mld_c_id_solver_apply_vect
    procedure, pass(sv) :: apply_a => mld_c_id_solver_apply
    procedure, pass(sv) :: free    => c_id_solver_free
    procedure, pass(sv) :: descr   => c_id_solver_descr
    procedure, nopass   :: get_fmt   => c_id_solver_get_fmt
  end type mld_c_id_solver_type


  private :: c_id_solver_bld, &
       &  c_id_solver_free, c_id_solver_get_fmt, &
       &  c_id_solver_descr

  interface 
    subroutine mld_c_id_solver_apply_vect(alpha,sv,x,beta,y,desc_data,&
         & trans,work,info,init,initu)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat, &
           & psb_c_vect_type, psb_c_base_vect_type, psb_spk_, & 
           & mld_c_id_solver_type, psb_ipk_
      type(psb_desc_type), intent(in)              :: desc_data
      class(mld_c_id_solver_type), intent(inout) :: sv
      type(psb_c_vect_type),intent(inout)        :: x
      type(psb_c_vect_type),intent(inout)        :: y
      complex(psb_spk_),intent(in)                   :: alpha,beta
      character(len=1),intent(in)                  :: trans
      complex(psb_spk_),target, intent(inout)        :: work(:)
      integer(psb_ipk_), intent(out)               :: info
      character, intent(in), optional                :: init
      type(psb_c_vect_type),intent(inout), optional   :: initu
    end subroutine mld_c_id_solver_apply_vect
  end interface
  
  interface 
    subroutine mld_c_id_solver_apply(alpha,sv,x,beta,y,desc_data,&
         & trans,work,info,init,initu)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat, &
           & psb_c_vect_type, psb_c_base_vect_type, psb_spk_, &
           & mld_c_id_solver_type, psb_ipk_
      type(psb_desc_type), intent(in)      :: desc_data
      class(mld_c_id_solver_type), intent(inout) :: sv
      complex(psb_spk_),intent(inout)         :: x(:)
      complex(psb_spk_),intent(inout)         :: y(:)
      complex(psb_spk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)          :: trans
      complex(psb_spk_),target, intent(inout) :: work(:)
      integer(psb_ipk_), intent(out)                 :: info
      character, intent(in), optional       :: init
      complex(psb_spk_),intent(inout), optional :: initu(:)
    end subroutine mld_c_id_solver_apply
  end interface

  interface
    subroutine mld_c_id_solver_clone(sv,svout,info)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat, &
           & psb_c_vect_type, psb_c_base_vect_type, psb_spk_, &
           & mld_c_base_solver_type, mld_c_id_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(mld_c_id_solver_type), intent(inout)                :: sv
      class(mld_c_base_solver_type), allocatable, intent(inout) :: svout
      integer(psb_ipk_), intent(out)             :: info
    end subroutine mld_c_id_solver_clone
  end interface

contains


  subroutine c_id_solver_bld(a,desc_a,sv,upd,info,b,amold,vmold,imold)

    Implicit None

    ! Arguments
    type(psb_cspmat_type), intent(in), target           :: a
    Type(psb_desc_type), Intent(in)                     :: desc_a 
    class(mld_c_id_solver_type), intent(inout)          :: sv
    character, intent(in)                               :: upd
    integer(psb_ipk_), intent(out)                      :: info
    type(psb_cspmat_type), intent(in), target, optional :: b
    class(psb_c_base_sparse_mat), intent(in), optional  :: amold
    class(psb_c_base_vect_type), intent(in), optional   :: vmold
    class(psb_i_base_vect_type), intent(in), optional   :: imold
    ! Local variables
    integer(psb_ipk_) :: n_row,n_col, nrow_a, nztota
    complex(psb_spk_), pointer :: ww(:), aux(:), tx(:),ty(:)
    integer(psb_ipk_) :: i, err_act, debug_unit, debug_level
    character(len=20)  :: name='c_id_solver_bld', ch_err
    
    info=psb_success_

    return
  end subroutine c_id_solver_bld

  subroutine c_id_solver_free(sv,info)

    Implicit None

    ! Arguments
    class(mld_c_id_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(out)               :: info
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='c_id_solver_free'

    info = psb_success_

    return
  end subroutine c_id_solver_free

  subroutine c_id_solver_descr(sv,info,iout,coarse)

    Implicit None

    ! Arguments
    class(mld_c_id_solver_type), intent(in) :: sv
    integer(psb_ipk_), intent(out)            :: info
    integer(psb_ipk_), intent(in), optional   :: iout
    logical, intent(in), optional       :: coarse

    ! Local variables
    integer(psb_ipk_)      :: err_act
    character(len=20), parameter :: name='mld_c_id_solver_descr'
    integer(psb_ipk_) :: iout_

    info = psb_success_
    if (present(iout)) then 
      iout_ = iout 
    else
      iout_ = 6
    endif
    
    write(iout_,*) '  Identity local solver '

    return

  end subroutine c_id_solver_descr

  function c_id_solver_get_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "Identity solver"
  end function c_id_solver_get_fmt


end module mld_c_id_solver
