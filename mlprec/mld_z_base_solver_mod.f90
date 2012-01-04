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
! File: mld_z_base_solver_mod.f90
!
! Module: mld_z_base_solver_mod
!
!  This module defines: 
!  - the mld_z_base_solver_type data structure containing the
!    basic solver type acting on a subdomain
!
!  It contains routines for
!  - Building and applying; 
!  - checking if the solver is correctly defined;
!  - printing a	description of the solver;
!  - deallocating the data structure.  
!

module mld_z_base_solver_mod

  use mld_base_prec_type
  use psb_base_mod, only : psb_desc_type, psb_zspmat_type, psb_long_int_k_, &
       & psb_sizeof, psb_free, psb_cdfree, psb_errpush, psb_act_abort_,&
       & psb_erractionsave, psb_erractionrestore, psb_error, psb_get_errstatus, &
       & psb_success_, psb_err_alloc_dealloc_, &
       & psb_z_vect_type, psb_z_base_vect_type, psb_z_base_sparse_mat, psb_dpk_
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

  type mld_z_base_solver_type
  contains
    procedure, pass(sv) :: check => mld_z_base_solver_check
    procedure, pass(sv) :: dump  => mld_z_base_solver_dmp
    procedure, pass(sv) :: build => mld_z_base_solver_bld
    procedure, pass(sv) :: apply_v => mld_z_base_solver_apply_vect
    procedure, pass(sv) :: apply_a => mld_z_base_solver_apply
    generic, public     :: apply => apply_a, apply_v
    procedure, pass(sv) :: free  => mld_z_base_solver_free
    procedure, pass(sv) :: seti  => mld_z_base_solver_seti
    procedure, pass(sv) :: setc  => mld_z_base_solver_setc
    procedure, pass(sv) :: setr  => mld_z_base_solver_setr
    generic, public     :: set   => seti, setc, setr
    procedure, pass(sv) :: default => z_base_solver_default
    procedure, pass(sv) :: descr   => mld_z_base_solver_descr
    procedure, pass(sv) :: sizeof  => z_base_solver_sizeof
    procedure, pass(sv) :: get_nzeros => z_base_solver_get_nzeros
  end type mld_z_base_solver_type

  private :: z_base_solver_sizeof, z_base_solver_default,&
       &  z_base_solver_get_nzeros


  interface  mld_z_base_solver_apply
    subroutine mld_z_base_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
       & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, mld_z_base_solver_type
      type(psb_desc_type), intent(in)           :: desc_data
      class(mld_z_base_solver_type), intent(in) :: sv
      complex(psb_dpk_),intent(inout)              :: x(:)
      complex(psb_dpk_),intent(inout)              :: y(:)
      complex(psb_dpk_),intent(in)                 :: alpha,beta
      character(len=1),intent(in)               :: trans
      complex(psb_dpk_),target, intent(inout)      :: work(:)
      integer, intent(out)                      :: info
    end subroutine mld_z_base_solver_apply
  end interface mld_z_base_solver_apply
  
      
  interface mld_z_base_solver_apply_vect
    subroutine mld_z_base_solver_apply_vect(alpha,sv,x,beta,y,desc_data,trans,work,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, mld_z_base_solver_type
      
      type(psb_desc_type), intent(in)              :: desc_data
      class(mld_z_base_solver_type), intent(inout) :: sv
      type(psb_z_vect_type),intent(inout)          :: x
      type(psb_z_vect_type),intent(inout)          :: y
      complex(psb_dpk_),intent(in)                    :: alpha,beta
      character(len=1),intent(in)                  :: trans
      complex(psb_dpk_),target, intent(inout)         :: work(:)
      integer, intent(out)                         :: info
    end subroutine mld_z_base_solver_apply_vect
  end interface mld_z_base_solver_apply_vect
  
  interface mld_z_base_solver_bld
    subroutine mld_z_base_solver_bld(a,desc_a,sv,upd,info,b,amold,vmold)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
       & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, mld_z_base_solver_type
      
      Implicit None
      
      ! Arguments
      type(psb_zspmat_type), intent(in), target           :: a
      Type(psb_desc_type), Intent(in)                     :: desc_a 
      class(mld_z_base_solver_type), intent(inout)        :: sv
      character, intent(in)                               :: upd
      integer, intent(out)                                :: info
      type(psb_zspmat_type), intent(in), target, optional :: b
      class(psb_z_base_sparse_mat), intent(in), optional  :: amold
      class(psb_z_base_vect_type), intent(in), optional   :: vmold
    end subroutine mld_z_base_solver_bld
  end interface mld_z_base_solver_bld
  
  interface mld_z_base_solver_check
    subroutine mld_z_base_solver_check(sv,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, mld_z_base_solver_type

      Implicit None
      
      ! Arguments
      class(mld_z_base_solver_type), intent(inout) :: sv
      integer, intent(out)                   :: info
    end subroutine mld_z_base_solver_check
  end interface mld_z_base_solver_check
  
  interface mld_z_base_solver_seti
    subroutine mld_z_base_solver_seti(sv,what,val,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, mld_z_base_solver_type
      
      Implicit None
      
      ! Arguments
      class(mld_z_base_solver_type), intent(inout) :: sv 
      integer, intent(in)                          :: what 
      integer, intent(in)                          :: val
      integer, intent(out)                         :: info
    end subroutine mld_z_base_solver_seti
  end interface mld_z_base_solver_seti
  
  interface  mld_z_base_solver_setc
    subroutine mld_z_base_solver_setc(sv,what,val,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, mld_z_base_solver_type
      Implicit None
      
      ! Arguments
      class(mld_z_base_solver_type), intent(inout) :: sv
      integer, intent(in)                          :: what 
      character(len=*), intent(in)                 :: val
      integer, intent(out)                         :: info
    end subroutine mld_z_base_solver_setc
  end interface mld_z_base_solver_setc
  
  interface  mld_z_base_solver_setr
    subroutine mld_z_base_solver_setr(sv,what,val,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, mld_z_base_solver_type
            
      Implicit None
      
      ! Arguments
      class(mld_z_base_solver_type), intent(inout) :: sv 
      integer, intent(in)                          :: what 
      real(psb_dpk_), intent(in)                   :: val
      integer, intent(out)                         :: info
    end subroutine mld_z_base_solver_setr
  end interface mld_z_base_solver_setr
  
  interface  mld_z_base_solver_free
    subroutine mld_z_base_solver_free(sv,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, mld_z_base_solver_type
      Implicit None
      
      ! Arguments
      class(mld_z_base_solver_type), intent(inout) :: sv
      integer, intent(out)                         :: info
    end subroutine mld_z_base_solver_free
  end interface mld_z_base_solver_free
  
  interface  mld_z_base_solver_descr
    subroutine mld_z_base_solver_descr(sv,info,iout,coarse)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, mld_z_base_solver_type
      
      Implicit None
      
      ! Arguments
      class(mld_z_base_solver_type), intent(in) :: sv
      integer, intent(out)                      :: info
      integer, intent(in), optional             :: iout
      logical, intent(in), optional             :: coarse

    end subroutine mld_z_base_solver_descr
  end interface mld_z_base_solver_descr
  
  interface  mld_z_base_solver_dmp
    subroutine mld_z_base_solver_dmp(sv,ictxt,level,info,prefix,head,solver)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, mld_z_base_solver_type
      
      implicit none 
      class(mld_z_base_solver_type), intent(in) :: sv
      integer, intent(in)              :: ictxt,level
      integer, intent(out)             :: info
      character(len=*), intent(in), optional :: prefix, head
      logical, optional, intent(in)    :: solver
    end subroutine mld_z_base_solver_dmp
  end interface mld_z_base_solver_dmp
  



contains
  !
  ! Function returning the size of the data structure
  ! in bytes or in number of nonzeros of the operator(s) involved. 
  !

  function z_base_solver_sizeof(sv) result(val)
    implicit none 
    ! Arguments
    class(mld_z_base_solver_type), intent(in) :: sv
    integer(psb_long_int_k_)                  :: val
    integer             :: i
    val = 0

    return
  end function z_base_solver_sizeof

  function z_base_solver_get_nzeros(sv) result(val)
    implicit none 
    class(mld_z_base_solver_type), intent(in) :: sv
    integer(psb_long_int_k_) :: val
    integer             :: i
    val = 0
  end function z_base_solver_get_nzeros

  subroutine z_base_solver_default(sv) 
    implicit none 
    ! Arguments
    class(mld_z_base_solver_type), intent(inout) :: sv
    ! Do nothing for base version

    return
  end subroutine z_base_solver_default

end module mld_z_base_solver_mod
