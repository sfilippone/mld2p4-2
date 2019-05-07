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
! File: mld_z_jac_smoother_mod.f90
!
! Module: mld_z_jac_smoother_mod
!
!  This module defines: 
!    the mld_z_jac_smoother_type data structure containing the
!    smoother for a Jacobi/block Jacobi smoother.
!  The smoother stores in ND the block off-diagonal matrix.
!  One special case is treated separately, when the solver is DIAG
!  then the ND is the entire off-diagonal part of the matrix (including the
!  main diagonal block), so that it becomes possible to implement
!  pure Jacobi global solver. 
! 
module mld_z_jac_smoother

  use mld_z_base_smoother_mod

  type, extends(mld_z_base_smoother_type) :: mld_z_jac_smoother_type
    ! The local solver component is inherited from the
    ! parent type. 
    !    class(mld_z_base_solver_type), allocatable :: sv
    !    
    type(psb_zspmat_type), pointer  :: pa => null()
    type(psb_zspmat_type) :: nd
    integer(psb_ipk_)       :: nnz_nd_tot
  contains
    procedure, pass(sm) :: dump    => mld_z_jac_smoother_dmp
    procedure, pass(sm) :: build   => mld_z_jac_smoother_bld
    procedure, pass(sm) :: cnv     => mld_z_jac_smoother_cnv
    procedure, pass(sm) :: clone   => mld_z_jac_smoother_clone
    procedure, pass(sm) :: apply_v => mld_z_jac_smoother_apply_vect
    procedure, pass(sm) :: apply_a => mld_z_jac_smoother_apply
    procedure, pass(sm) :: free    => z_jac_smoother_free
    procedure, pass(sm) :: descr   => mld_z_jac_smoother_descr
    procedure, pass(sm) :: sizeof  => z_jac_smoother_sizeof
    procedure, pass(sm) :: get_nzeros => z_jac_smoother_get_nzeros
    procedure, pass(sm) :: get_wrksz => z_jac_smoother_get_wrksize
    procedure, nopass   :: get_fmt    => z_jac_smoother_get_fmt
    procedure, nopass   :: get_id     => z_jac_smoother_get_id
  end type mld_z_jac_smoother_type


  private :: z_jac_smoother_free,   &
       & z_jac_smoother_sizeof,  z_jac_smoother_get_nzeros, &
       & z_jac_smoother_get_fmt, z_jac_smoother_get_id, &
       & z_jac_smoother_get_wrksize


  interface 
    subroutine mld_z_jac_smoother_apply_vect(alpha,sm,x,beta,y,desc_data,trans,& 
         & sweeps,work,wv,info,init,initu)
      import :: psb_desc_type, mld_z_jac_smoother_type, psb_z_vect_type, psb_dpk_, &
           & psb_zspmat_type, psb_z_base_sparse_mat, psb_z_base_vect_type,&
           & psb_ipk_
       
      type(psb_desc_type), intent(in)                 :: desc_data
      class(mld_z_jac_smoother_type), intent(inout) :: sm
      type(psb_z_vect_type),intent(inout)           :: x
      type(psb_z_vect_type),intent(inout)           :: y
      complex(psb_dpk_),intent(in)                      :: alpha,beta
      character(len=1),intent(in)                     :: trans
      integer(psb_ipk_), intent(in)                   :: sweeps
      complex(psb_dpk_),target, intent(inout)           :: work(:)
      type(psb_z_vect_type),intent(inout)           :: wv(:)
      integer(psb_ipk_), intent(out)                  :: info
      character, intent(in), optional                :: init
      type(psb_z_vect_type),intent(inout), optional   :: initu
    end subroutine mld_z_jac_smoother_apply_vect
  end interface
  
  interface 
    subroutine mld_z_jac_smoother_apply(alpha,sm,x,beta,y,desc_data,trans,& 
         & sweeps,work,info,init,initu)
      import :: psb_desc_type, mld_z_jac_smoother_type, psb_z_vect_type, psb_dpk_, &
           & psb_zspmat_type, psb_z_base_sparse_mat, psb_z_base_vect_type, &
           & psb_ipk_
      type(psb_desc_type), intent(in)      :: desc_data
      class(mld_z_jac_smoother_type), intent(inout) :: sm
      complex(psb_dpk_),intent(inout)         :: x(:)
      complex(psb_dpk_),intent(inout)         :: y(:)
      complex(psb_dpk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)           :: trans
      integer(psb_ipk_), intent(in)         :: sweeps
      complex(psb_dpk_),target, intent(inout) :: work(:)
      integer(psb_ipk_), intent(out)        :: info
      character, intent(in), optional       :: init
      complex(psb_dpk_),intent(inout), optional :: initu(:)
    end subroutine mld_z_jac_smoother_apply
  end interface
  
  interface 
    subroutine mld_z_jac_smoother_bld(a,desc_a,sm,info,amold,vmold,imold)
      import :: psb_desc_type, mld_z_jac_smoother_type, psb_z_vect_type, psb_dpk_, &
           & psb_zspmat_type, psb_z_base_sparse_mat, psb_z_base_vect_type,&
           & psb_ipk_, psb_i_base_vect_type
      type(psb_zspmat_type), intent(in), target           :: a
      Type(psb_desc_type), Intent(inout)                  :: desc_a 
      class(mld_z_jac_smoother_type), intent(inout)       :: sm
      integer(psb_ipk_), intent(out)                      :: info
      class(psb_z_base_sparse_mat), intent(in), optional :: amold
      class(psb_z_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine mld_z_jac_smoother_bld
  end interface
  
  interface 
    subroutine mld_z_jac_smoother_cnv(sm,info,amold,vmold,imold)
      import :: mld_z_jac_smoother_type, psb_dpk_, &
           & psb_z_base_sparse_mat, psb_z_base_vect_type,&
           & psb_ipk_, psb_i_base_vect_type
      class(mld_z_jac_smoother_type), intent(inout)       :: sm
      integer(psb_ipk_), intent(out)                      :: info
      class(psb_z_base_sparse_mat), intent(in), optional :: amold
      class(psb_z_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine mld_z_jac_smoother_cnv
  end interface
  
  interface 
    subroutine mld_z_jac_smoother_dmp(sm,ictxt,level,info,prefix,head,smoother,solver)
      import :: psb_zspmat_type, psb_z_vect_type, psb_z_base_vect_type, &
           & psb_dpk_, mld_z_jac_smoother_type, psb_long_int_k_, psb_desc_type, &
           & psb_ipk_
      implicit none 
      class(mld_z_jac_smoother_type), intent(in) :: sm
      integer(psb_ipk_), intent(in)               :: ictxt
      integer(psb_ipk_), intent(in)               :: level
      integer(psb_ipk_), intent(out)              :: info
      character(len=*), intent(in), optional :: prefix, head
      logical, optional, intent(in)    :: smoother, solver
    end subroutine mld_z_jac_smoother_dmp
  end interface
  
  interface 
    subroutine mld_z_jac_smoother_clone(sm,smout,info)
      import :: mld_z_jac_smoother_type, psb_dpk_, &
           & mld_z_base_smoother_type, psb_ipk_
      class(mld_z_jac_smoother_type), intent(inout)               :: sm
      class(mld_z_base_smoother_type), allocatable, intent(inout) :: smout
      integer(psb_ipk_), intent(out)                :: info
    end subroutine mld_z_jac_smoother_clone
  end interface

  interface
    subroutine mld_z_jac_smoother_descr(sm,info,iout,coarse)
      import :: mld_z_jac_smoother_type, psb_ipk_
      class(mld_z_jac_smoother_type), intent(in) :: sm
      integer(psb_ipk_), intent(out)               :: info
      integer(psb_ipk_), intent(in), optional      :: iout
      logical, intent(in), optional              :: coarse
    end subroutine mld_z_jac_smoother_descr
  end interface
    
contains


  subroutine z_jac_smoother_free(sm,info)


    Implicit None

    ! Arguments
    class(mld_z_jac_smoother_type), intent(inout) :: sm
    integer(psb_ipk_), intent(out)                  :: info
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='z_jac_smoother_free'

    call psb_erractionsave(err_act)
    info = psb_success_



    if (allocated(sm%sv)) then 
      call sm%sv%free(info)
      if (info == psb_success_) deallocate(sm%sv,stat=info)
      if (info /= psb_success_) then 
        info = psb_err_alloc_dealloc_
        call psb_errpush(info,name)
        goto 9999 
      end if
    end if
    call sm%nd%free()
    sm%pa => null()
    
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine z_jac_smoother_free

  function z_jac_smoother_sizeof(sm) result(val)

    implicit none 
    ! Arguments
    class(mld_z_jac_smoother_type), intent(in) :: sm
    integer(psb_long_int_k_) :: val
    integer(psb_ipk_)        :: i

    val = psb_sizeof_int 
    if (allocated(sm%sv)) val = val + sm%sv%sizeof()
    val = val + sm%nd%sizeof()

    return
  end function z_jac_smoother_sizeof

  function z_jac_smoother_get_nzeros(sm) result(val)

    implicit none 
    ! Arguments
    class(mld_z_jac_smoother_type), intent(in) :: sm
    integer(psb_long_int_k_) :: val
    integer(psb_ipk_)        :: i

    val = 0
    if (allocated(sm%sv)) val = val + sm%sv%get_nzeros()
    val = val + sm%nd%get_nzeros()

    return
  end function z_jac_smoother_get_nzeros

  function z_jac_smoother_get_wrksize(sm) result(val)
    implicit none 
    class(mld_z_jac_smoother_type), intent(inout) :: sm
    integer(psb_ipk_)  :: val

    val = 2
    if (allocated(sm%sv)) val = val + sm%sv%get_wrksz()
    
  end function z_jac_smoother_get_wrksize
  
  function z_jac_smoother_get_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "Jacobi smoother"
  end function z_jac_smoother_get_fmt

  function z_jac_smoother_get_id() result(val)
    implicit none 
    integer(psb_ipk_)  :: val

    val = mld_jac_
  end function z_jac_smoother_get_id
  
end module mld_z_jac_smoother
