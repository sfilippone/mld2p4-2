!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010,2012
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
! File: mld_c_base_smoother_mod.f90
!
! Module: mld_c_base_smoother_mod
!
!  This module defines: 
!  - the mld_c_base_smoother_type data structure containing the
!    smoother  and related data structures;
!
!  It contains routines for
!  - Building and applying; 
!  - checking if the smoother is correctly defined;
!  - printing a	description of the preconditioner;
!  - deallocating the preconditioner data structure.  
!

module mld_c_base_smoother_mod

  use mld_c_base_solver_mod
  use psb_base_mod, only : psb_desc_type, psb_cspmat_type, psb_long_int_k_,&
       & psb_c_vect_type, psb_c_base_vect_type, psb_c_base_sparse_mat, psb_spk_
  
  !
  !
  ! 
  ! Type: mld_T_base_smoother_type.
  ! 
  !  It holds the smoother a single level. Its only mandatory component is a solver
  !  object which holds a local solver; this decoupling allows to have the same solver
  !  e.g ILU to work with Jacobi with multiple sweeps as well as with any AS variant.
  !
  !  type  mld_T_base_smoother_type
  !    class(mld_T_base_solver_type), allocatable :: sv
  !  end type mld_T_base_smoother_type
  !
  !   Methods:  
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
  !                   it is passed to the solver object for further processing.
  !    check      -   Sanity checks.
  !    sizeof     -   Total memory occupation in bytes
  !    get_nzeros -   Number of nonzeros 
  !
  !
  ! 

  type  mld_c_base_smoother_type
    class(mld_c_base_solver_type), allocatable :: sv
  contains
    procedure, pass(sm) :: check => mld_c_base_smoother_check
    procedure, pass(sm) :: dump  => mld_c_base_smoother_dmp
    procedure, pass(sm) :: build => mld_c_base_smoother_bld
    procedure, pass(sm) :: apply_v => mld_c_base_smoother_apply_vect
    procedure, pass(sm) :: apply_a => mld_c_base_smoother_apply
    generic, public     :: apply => apply_a, apply_v
    procedure, pass(sm) :: free  => mld_c_base_smoother_free
    procedure, pass(sm) :: seti  => mld_c_base_smoother_seti
    procedure, pass(sm) :: setc  => mld_c_base_smoother_setc
    procedure, pass(sm) :: setr  => mld_c_base_smoother_setr
    generic, public     :: set   => seti, setc, setr
    procedure, pass(sm) :: default => c_base_smoother_default
    procedure, pass(sm) :: descr =>   mld_c_base_smoother_descr
    procedure, pass(sm) :: sizeof =>  c_base_smoother_sizeof
    procedure, pass(sm) :: get_nzeros => c_base_smoother_get_nzeros
  end type mld_c_base_smoother_type


  private :: c_base_smoother_sizeof, &
       &  c_base_smoother_default, c_base_smoother_get_nzeros



  interface 
    subroutine mld_c_base_smoother_apply(alpha,sm,x,beta,y,desc_data,trans,sweeps,work,info)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat, &
           & psb_c_vect_type, psb_c_base_vect_type, psb_spk_, &
           & mld_c_base_smoother_type, psb_ipk_
      type(psb_desc_type), intent(in)             :: desc_data
      class(mld_c_base_smoother_type), intent(in) :: sm
      complex(psb_spk_),intent(inout)                :: x(:)
      complex(psb_spk_),intent(inout)                :: y(:)
      complex(psb_spk_),intent(in)                   :: alpha,beta
      character(len=1),intent(in)                  :: trans
      integer(psb_ipk_), intent(in)                :: sweeps
      complex(psb_spk_),target, intent(inout)        :: work(:)
      integer(psb_ipk_), intent(out)               :: info
    end subroutine mld_c_base_smoother_apply
  end interface
  
  interface 
    subroutine mld_c_base_smoother_apply_vect(alpha,sm,x,beta,y,desc_data,&
         &  trans,sweeps,work,info)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat, &
           & psb_c_vect_type, psb_c_base_vect_type, psb_spk_, &
           & mld_c_base_smoother_type, psb_ipk_
      type(psb_desc_type), intent(in)                :: desc_data
      class(mld_c_base_smoother_type), intent(inout) :: sm
      type(psb_c_vect_type),intent(inout)            :: x
      type(psb_c_vect_type),intent(inout)            :: y
      complex(psb_spk_),intent(in)                       :: alpha,beta
      character(len=1),intent(in)                      :: trans
      integer(psb_ipk_), intent(in)                    :: sweeps
      complex(psb_spk_),target, intent(inout)            :: work(:)
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine mld_c_base_smoother_apply_vect
  end interface
  
  interface 
    subroutine mld_c_base_smoother_check(sm,info)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat, &
           & psb_c_vect_type, psb_c_base_vect_type, psb_spk_, &
           & mld_c_base_smoother_type, psb_ipk_
      ! Arguments
      class(mld_c_base_smoother_type), intent(inout) :: sm 
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine mld_c_base_smoother_check
  end interface
  
  interface 
    subroutine mld_c_base_smoother_seti(sm,what,val,info)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat, &
           & psb_c_vect_type, psb_c_base_vect_type, psb_spk_, &
           & mld_c_base_smoother_type, psb_ipk_
      ! Arguments
      class(mld_c_base_smoother_type), intent(inout) :: sm 
      integer(psb_ipk_), intent(in)                    :: what 
      integer(psb_ipk_), intent(in)                    :: val
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine mld_c_base_smoother_seti
  end interface
  
  interface 
    subroutine mld_c_base_smoother_setc(sm,what,val,info)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat, &
           & psb_c_vect_type, psb_c_base_vect_type, psb_spk_, &
           & mld_c_base_smoother_type, psb_ipk_
      class(mld_c_base_smoother_type), intent(inout) :: sm 
      integer(psb_ipk_), intent(in)                    :: what 
      character(len=*), intent(in)                     :: val
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine mld_c_base_smoother_setc
  end interface
  
  interface 
    subroutine mld_c_base_smoother_setr(sm,what,val,info)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat, &
           & psb_c_vect_type, psb_c_base_vect_type, psb_spk_, &
           & mld_c_base_smoother_type, psb_ipk_
      ! Arguments
      class(mld_c_base_smoother_type), intent(inout) :: sm 
      integer(psb_ipk_), intent(in)                    :: what 
      real(psb_spk_), intent(in)                        :: val
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine mld_c_base_smoother_setr
  end interface
  
  interface 
    subroutine mld_c_base_smoother_bld(a,desc_a,sm,upd,info,amold,vmold)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat, &
           & psb_c_vect_type, psb_c_base_vect_type, psb_spk_, &
           & mld_c_base_smoother_type, psb_ipk_
      ! Arguments
      type(psb_cspmat_type), intent(in), target     :: a
      Type(psb_desc_type), Intent(in)                 :: desc_a 
      class(mld_c_base_smoother_type), intent(inout) :: sm 
      character, intent(in)                           :: upd
      integer(psb_ipk_), intent(out)                  :: info
      class(psb_c_base_sparse_mat), intent(in), optional :: amold
      class(psb_c_base_vect_type), intent(in), optional  :: vmold
    end subroutine mld_c_base_smoother_bld
  end interface
  
  interface 
    subroutine mld_c_base_smoother_free(sm,info)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat, &
           & psb_c_vect_type, psb_c_base_vect_type, psb_spk_, &
           & mld_c_base_smoother_type, psb_ipk_
      ! Arguments
      class(mld_c_base_smoother_type), intent(inout) :: sm
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine mld_c_base_smoother_free
  end interface
  
  interface 
    subroutine mld_c_base_smoother_descr(sm,info,iout,coarse)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat, &
           & psb_c_vect_type, psb_c_base_vect_type, psb_spk_, &
           & mld_c_base_smoother_type, psb_ipk_
      ! Arguments
      class(mld_c_base_smoother_type), intent(in) :: sm
      integer(psb_ipk_), intent(out)                :: info
      integer(psb_ipk_), intent(in), optional       :: iout
      logical, intent(in), optional                 :: coarse
    end subroutine mld_c_base_smoother_descr
  end interface
  
  interface 
    subroutine mld_c_base_smoother_dmp(sm,ictxt,level,info,prefix,head,smoother,solver)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat, &
           & psb_c_vect_type, psb_c_base_vect_type, psb_spk_, &
           & mld_c_base_smoother_type, psb_ipk_
      class(mld_c_base_smoother_type), intent(in) :: sm
      integer(psb_ipk_), intent(in)              :: ictxt
      integer(psb_ipk_), intent(in)              :: level
      integer(psb_ipk_), intent(out)             :: info
      character(len=*), intent(in), optional :: prefix, head
      logical, optional, intent(in)    :: smoother, solver
    end subroutine mld_c_base_smoother_dmp
  end interface
  
contains
  !
  ! Function returning the size of the mld_prec_type data structure
  ! in bytes or in number of nonzeros of the operator(s) involved. 
  !

  function c_base_smoother_get_nzeros(sm) result(val)
    implicit none 
    class(mld_c_base_smoother_type), intent(in) :: sm
    integer(psb_long_int_k_) :: val
    integer(psb_ipk_)             :: i
    val = 0
    if (allocated(sm%sv)) &
         &  val =  sm%sv%get_nzeros()
  end function c_base_smoother_get_nzeros

  function c_base_smoother_sizeof(sm) result(val)
    implicit none 
    ! Arguments
    class(mld_c_base_smoother_type), intent(in) :: sm
    integer(psb_long_int_k_)                    :: val
    integer(psb_ipk_)             :: i
    
    val = 0
    if (allocated(sm%sv)) then 
      val = sm%sv%sizeof()
    end if

    return
  end function c_base_smoother_sizeof

  !
  ! Set sensible defaults.
  ! To be called immediately after allocation
  !
  subroutine c_base_smoother_default(sm) 
    implicit none 
    ! Arguments
    class(mld_c_base_smoother_type), intent(inout) :: sm
    ! Do nothing for base version

    if (allocated(sm%sv)) call sm%sv%default()

    return
  end subroutine c_base_smoother_default


end module mld_c_base_smoother_mod
