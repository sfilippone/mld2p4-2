!  
!   
!                             MLD2P4  version 2.1
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!    
!    (C) Copyright 2008, 2010, 2012, 2015, 2017 
!  
!        Salvatore Filippone    Cranfield University, UK
!        Ambra Abdullahi Hassan University of Rome Tor Vergata, IT
!        Alfredo Buttari        CNRS-IRIT, Toulouse, FR
!        Pasqua D'Ambra         IAC-CNR, Naples, IT
!        Daniela di Serafino    University of Campania "L. Vanvitelli", Caserta, IT
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
! File: mld_d_tlu_solver_mod.f90
!
! Module: mld_d_tlu_solver_mod
!
!  This module serves as an example of how to define a new solver and integrate
!  it in MLD2P4 via the P%SET(sv,info) method.
!  In this example we are extending the ILU solver by implementing a new factorization algorithm.
!  In actual reality, we are just giving a new name to ILU(0), but this should be sufficient to show
!  the basics.
!
!  The code is divided in two files:
!  1. The interface file (this one)
!  2. The implementation file (mld_d_tlu_solver_impl.f90)
!  
!  The separation between interface and implementation is an essential part of the
!  object-oriented design. The most appropriate tool would be to have the implementation
!  in a SUBMODULE, something which we plan to do at some point in the future but will
!  need to check for compiler support. 
!  
!
!

module mld_d_tlu_solver

  use mld_d_ilu_solver
  !  use mld_d_ilu_fact_mod

  type, extends(mld_d_ilu_solver_type) :: mld_d_tlu_solver_type
    !
    ! These are already defined in the ILU solver type; since we
    ! are supposedly implementing a new factorization strategy, the
    ! data components are the same, and we only need to change
    ! the methods that build the solver. 
    ! 
    !type(psb_dspmat_type)       :: l, u
    !real(psb_dpk_), allocatable :: d(:)
    !type(psb_d_vect_type)      :: dv
    !integer(psb_ipk_)            :: fact_type, fill_in
    !real(psb_dpk_)                :: thresh
  contains
    !
    ! Again, we can freely reuse these methods because they would be
    ! in common among all possible ILU factorizations
    ! 
    ! 
    !procedure, pass(sv) :: dump    => mld_d_tlu_solver_dmp 
    !procedure, pass(sv) :: ccheck  => d_tlu_solver_check
    !procedure, pass(sv) :: clone   => mld_d_tlu_solver_clone
    !procedure, pass(sv) :: cnv     => mld_d_tlu_solver_cnv
    !procedure, pass(sv) :: apply_v => mld_d_tlu_solver_apply_vect
    !procedure, pass(sv) :: apply_a => mld_d_tlu_solver_apply
    !procedure, pass(sv) :: free    => d_tlu_solver_free
    !procedure, pass(sv) :: seti    => d_tlu_solver_seti
    !procedure, pass(sv) :: setc    => d_tlu_solver_setc
    !procedure, pass(sv) :: setr    => d_tlu_solver_setr
    !procedure, pass(sv) :: cseti   => d_tlu_solver_cseti
    !procedure, pass(sv) :: csetc   => d_tlu_solver_csetc
    !procedure, pass(sv) :: csetr   => d_tlu_solver_csetr
    !procedure, pass(sv) :: descr   => d_tlu_solver_descr    
    !procedure, pass(sv) :: sizeof  => d_tlu_solver_sizeof
    !procedure, pass(sv) :: get_nzeros => d_tlu_solver_get_nzeros

    
    !
    ! These methods are specific for the new solver type
    ! and therefore need to be overridden
    ! 
    procedure, pass(sv) :: default => d_tlu_solver_default
    procedure, pass(sv) :: build   => mld_d_tlu_solver_bld
    procedure, nopass   :: get_fmt => d_tlu_solver_get_fmt
    procedure, nopass   :: get_id  => d_tlu_solver_get_id
  end type mld_d_tlu_solver_type


  private ::  d_tlu_solver_get_fmt,  d_tlu_solver_get_id

  interface 
    subroutine mld_d_tlu_solver_bld(a,desc_a,sv,info,b,amold,vmold, imold)
      import :: psb_desc_type, mld_d_tlu_solver_type, psb_d_vect_type, psb_dpk_, &
           & psb_dspmat_type, psb_d_base_sparse_mat, psb_d_base_vect_type,&
           & psb_ipk_, psb_i_base_vect_type
      implicit none 
      type(psb_dspmat_type), intent(in), target           :: a
      Type(psb_desc_type), Intent(in)                     :: desc_a 
      class(mld_d_tlu_solver_type), intent(inout)         :: sv
      integer(psb_ipk_), intent(out)                      :: info
      type(psb_dspmat_type), intent(in), target, optional :: b
      class(psb_d_base_sparse_mat), intent(in), optional  :: amold
      class(psb_d_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional   :: imold
    end subroutine mld_d_tlu_solver_bld
  end interface

contains

  !
  ! This method is a copy of the ILU method, because for this example we are
  ! only repackaging existing code. 
  !
  subroutine d_tlu_solver_default(sv)
    
    Implicit None

    ! Arguments
    class(mld_d_tlu_solver_type), intent(inout) :: sv

    sv%fact_type = mld_ilu_n_
    sv%fill_in   = 0
    sv%thresh    = dzero

    return
  end subroutine d_tlu_solver_default

  function d_tlu_solver_get_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "TLU solver"
  end function d_tlu_solver_get_fmt

  function d_tlu_solver_get_id() result(val)
    implicit none 
    integer(psb_ipk_)  :: val
    
    val = mld_ilu_n_
  end function d_tlu_solver_get_id
  
end module mld_d_tlu_solver
