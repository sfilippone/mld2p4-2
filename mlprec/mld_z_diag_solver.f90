!  
!   
!                             MLD2P4  version 2.2
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!    
!    (C) Copyright 2008-2018 
!  
!        Salvatore Filippone  
!        Pasqua D'Ambra   
!        Daniela di Serafino   
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
! File: mld_z_diag_solver_mod.f90
!
! Module: mld_z_diag_solver_mod
!
!  This module defines: 
!  - the mld_z_diag_solver_type data structure containing the 
!    simple diagonal solver. This extracts the main diagonal of a matrix
!    and precomputes its inverse. Combined with a Jacobi "smoother" generates
!    what are commonly known as the classic Jacobi iterations
!
module mld_z_diag_solver

  use mld_z_base_solver_mod

  type, extends(mld_z_base_solver_type) :: mld_z_diag_solver_type
    type(psb_z_vect_type), allocatable :: dv
    complex(psb_dpk_), allocatable        :: d(:)
  contains
    procedure, pass(sv) :: dump    => mld_z_diag_solver_dmp
    procedure, pass(sv) :: build   => mld_z_diag_solver_bld
    procedure, pass(sv) :: cnv     => mld_z_diag_solver_cnv
    procedure, pass(sv) :: clone   => mld_z_diag_solver_clone
    procedure, pass(sv) :: apply_v => mld_z_diag_solver_apply_vect
    procedure, pass(sv) :: apply_a => mld_z_diag_solver_apply
    procedure, pass(sv) :: free    => z_diag_solver_free
    procedure, pass(sv) :: descr   => z_diag_solver_descr
    procedure, pass(sv) :: sizeof  => z_diag_solver_sizeof
    procedure, pass(sv) :: get_nzeros  => z_diag_solver_get_nzeros
    procedure, nopass   :: get_fmt   => z_diag_solver_get_fmt
    procedure, nopass   :: get_id    => z_diag_solver_get_id
  end type mld_z_diag_solver_type


  private :: z_diag_solver_free,  z_diag_solver_descr, &
       & z_diag_solver_sizeof, z_diag_solver_get_nzeros, &
       & z_diag_solver_get_fmt, z_diag_solver_get_id


  interface 
    subroutine mld_z_diag_solver_apply_vect(alpha,sv,x,beta,y,desc_data,& 
         & trans,work,wv,info,init,initu)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
       & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
       & mld_z_diag_solver_type, psb_ipk_
      type(psb_desc_type), intent(in)                :: desc_data
      class(mld_z_diag_solver_type), intent(inout) :: sv
      type(psb_z_vect_type), intent(inout)         :: x
      type(psb_z_vect_type), intent(inout)         :: y
      complex(psb_dpk_),intent(in)                     :: alpha,beta
      character(len=1),intent(in)                    :: trans
      complex(psb_dpk_),target, intent(inout)          :: work(:)
      type(psb_z_vect_type),intent(inout)          :: wv(:)
      integer(psb_ipk_), intent(out)                 :: info
      character, intent(in), optional                :: init
      type(psb_z_vect_type),intent(inout), optional   :: initu
    end subroutine mld_z_diag_solver_apply_vect
  end interface
  
  interface 
    subroutine mld_z_diag_solver_apply(alpha,sv,x,beta,y,desc_data,&
         & trans,work,info,init,initu)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
       & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
       & mld_z_diag_solver_type, psb_ipk_
      type(psb_desc_type), intent(in)            :: desc_data
      class(mld_z_diag_solver_type), intent(inout) :: sv
      complex(psb_dpk_), intent(inout)             :: x(:)
      complex(psb_dpk_), intent(inout)             :: y(:)
      complex(psb_dpk_),intent(in)                 :: alpha,beta
      character(len=1),intent(in)                :: trans
      complex(psb_dpk_),target, intent(inout)      :: work(:)
      integer(psb_ipk_), intent(out)             :: info
      character, intent(in), optional       :: init
      complex(psb_dpk_),intent(inout), optional :: initu(:)
    end subroutine mld_z_diag_solver_apply
  end interface
  
  interface 
    subroutine mld_z_diag_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & mld_z_diag_solver_type, psb_ipk_, psb_i_base_vect_type      
      type(psb_zspmat_type), intent(in), target           :: a
      Type(psb_desc_type), Intent(inout)                    :: desc_a 
      class(mld_z_diag_solver_type), intent(inout)        :: sv
      integer(psb_ipk_), intent(out)                        :: info
      type(psb_zspmat_type), intent(in), target, optional :: b
      class(psb_z_base_sparse_mat), intent(in), optional  :: amold
      class(psb_z_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional   :: imold
    end subroutine mld_z_diag_solver_bld
  end interface
  
  interface 
    subroutine mld_z_diag_solver_cnv(sv,info,amold,vmold,imold)
      import :: psb_z_base_sparse_mat, psb_z_base_vect_type, psb_dpk_, &
           & mld_z_diag_solver_type, psb_ipk_, psb_i_base_vect_type      
      class(mld_z_diag_solver_type), intent(inout)        :: sv
      integer(psb_ipk_), intent(out)                        :: info
      class(psb_z_base_sparse_mat), intent(in), optional  :: amold
      class(psb_z_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional   :: imold
    end subroutine mld_z_diag_solver_cnv
  end interface

  interface 
    subroutine mld_z_diag_solver_dmp(sv,ictxt,level,info,prefix,head,solver)
      import :: psb_desc_type, mld_z_diag_solver_type, psb_z_vect_type, psb_dpk_, &
           & psb_zspmat_type, psb_z_base_sparse_mat, psb_z_base_vect_type, &
           & psb_ipk_
      implicit none 
      class(mld_z_diag_solver_type), intent(in) :: sv
      integer(psb_ipk_), intent(in)              :: ictxt
      integer(psb_ipk_), intent(in)              :: level
      integer(psb_ipk_), intent(out)             :: info
      character(len=*), intent(in), optional     :: prefix, head
      logical, optional, intent(in)              :: solver
    end subroutine mld_z_diag_solver_dmp
  end interface
  
  interface
    subroutine mld_z_diag_solver_clone(sv,svout,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & mld_z_base_solver_type, mld_z_diag_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(mld_z_diag_solver_type), intent(inout)              :: sv
      class(mld_z_base_solver_type), allocatable, intent(inout) :: svout
      integer(psb_ipk_), intent(out)               :: info
    end subroutine mld_z_diag_solver_clone
  end interface
  
  
contains

  subroutine z_diag_solver_free(sv,info)

    Implicit None

    ! Arguments
    class(mld_z_diag_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(out)                 :: info
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='z_diag_solver_free'

    call psb_erractionsave(err_act)
    info = psb_success_

    if (allocated(sv%dv)) call sv%dv%free(info)
    
    if (allocated(sv%d)) then 
      deallocate(sv%d,stat=info)
      if (info /= psb_success_) then 
        info = psb_err_alloc_dealloc_
        call psb_errpush(info,name)
        goto 9999 
      end if
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine z_diag_solver_free

  subroutine z_diag_solver_descr(sv,info,iout,coarse)

    Implicit None

    ! Arguments
    class(mld_z_diag_solver_type), intent(in) :: sv
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_), intent(in), optional     :: iout
    logical, intent(in), optional               :: coarse

    ! Local variables
    integer(psb_ipk_)      :: err_act
    character(len=20), parameter :: name='mld_z_diag_solver_descr'
    integer(psb_ipk_) :: iout_

    info = psb_success_
    if (present(iout)) then 
      iout_ = iout 
    else
      iout_ = psb_out_unit
    endif
    
    write(iout_,*) '  Diagonal local solver '

    return

  end subroutine z_diag_solver_descr

  function z_diag_solver_sizeof(sv) result(val)
    implicit none 
    ! Arguments
    class(mld_z_diag_solver_type), intent(in) :: sv
    integer(psb_long_int_k_) :: val
    integer(psb_ipk_)             :: i

    val = 0
    if (allocated(sv%dv)) val = val + sv%dv%sizeof()

    return
  end function z_diag_solver_sizeof

  function z_diag_solver_get_nzeros(sv) result(val)
    implicit none 
    ! Arguments
    class(mld_z_diag_solver_type), intent(in) :: sv
    integer(psb_long_int_k_) :: val
    integer(psb_ipk_)             :: i

    val = 0
    if (allocated(sv%dv)) val = val +  sv%dv%get_nrows()

    return
  end function z_diag_solver_get_nzeros

  function z_diag_solver_get_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "Diag solver"
  end function z_diag_solver_get_fmt

  function z_diag_solver_get_id() result(val)
    implicit none 
    integer(psb_ipk_)  :: val

    val = mld_diag_scale_
  end function z_diag_solver_get_id

end module mld_z_diag_solver
