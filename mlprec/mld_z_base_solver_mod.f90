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
  use psb_base_mod, only : psb_zspmat_type, &
       & psb_z_vect_type, psb_z_base_vect_type, psb_z_base_sparse_mat, &
       & psb_dpk_, psb_i_base_vect_type, psb_erractionsave, psb_error_handler
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
  !    stringval  -   convert string to val for internal parms
  !    get_fmt    -   short string descriptor
  !    get_id     -   numeric id descriptro
  !    get_wrksz  -   How many workspace vector does apply_vect need
  !
  !

  type mld_z_base_solver_type
  contains
    procedure, pass(sv) :: check => mld_z_base_solver_check
    procedure, pass(sv) :: dump  => mld_z_base_solver_dmp
    procedure, pass(sv) :: clone => mld_z_base_solver_clone
    procedure, pass(sv) :: build => mld_z_base_solver_bld
    procedure, pass(sv) :: cnv   => mld_z_base_solver_cnv
    procedure, pass(sv) :: apply_v => mld_z_base_solver_apply_vect
    procedure, pass(sv) :: apply_a => mld_z_base_solver_apply
    generic, public     :: apply => apply_a, apply_v
    procedure, pass(sv) :: free  => mld_z_base_solver_free
    procedure, pass(sv) :: seti  => mld_z_base_solver_seti
    procedure, pass(sv) :: setc  => mld_z_base_solver_setc
    procedure, pass(sv) :: setr  => mld_z_base_solver_setr
    procedure, pass(sv) :: cseti => mld_z_base_solver_cseti
    procedure, pass(sv) :: csetc => mld_z_base_solver_csetc
    procedure, pass(sv) :: csetr => mld_z_base_solver_csetr
    generic, public     :: set   => seti, setc, setr, cseti, csetc, csetr
    procedure, pass(sv) :: default => z_base_solver_default
    procedure, pass(sv) :: descr   => mld_z_base_solver_descr
    procedure, pass(sv) :: sizeof  => z_base_solver_sizeof
    procedure, pass(sv) :: get_nzeros => z_base_solver_get_nzeros
    procedure, nopass   :: get_wrksz => z_base_solver_get_wrksize
    procedure, nopass   :: stringval => mld_stringval
    procedure, nopass   :: get_fmt   => z_base_solver_get_fmt
    procedure, nopass   :: get_id    => z_base_solver_get_id
    procedure, nopass   :: is_iterative => z_base_solver_is_iterative
  end type mld_z_base_solver_type

  private :: z_base_solver_sizeof, z_base_solver_default,&
       &  z_base_solver_get_nzeros, z_base_solver_get_fmt, &
       &  z_base_solver_is_iterative, z_base_solver_get_id, &
       &  z_base_solver_get_wrksize


  interface  
    subroutine mld_z_base_solver_apply(alpha,sv,x,beta,y,desc_data,&
         & trans,work,info,init,initu)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
       & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
       & mld_z_base_solver_type, psb_ipk_
      implicit none 
      type(psb_desc_type), intent(in)           :: desc_data
      class(mld_z_base_solver_type), intent(inout) :: sv
      complex(psb_dpk_),intent(inout)              :: x(:)
      complex(psb_dpk_),intent(inout)              :: y(:)
      complex(psb_dpk_),intent(in)                 :: alpha,beta
      character(len=1),intent(in)                :: trans
      complex(psb_dpk_),target, intent(inout)      :: work(:)
      integer(psb_ipk_), intent(out)             :: info
      character, intent(in), optional       :: init
      complex(psb_dpk_),intent(inout), optional :: initu(:)
    end subroutine mld_z_base_solver_apply
  end interface 
  
      
  interface 
    subroutine mld_z_base_solver_apply_vect(alpha,sv,x,beta,y,desc_data,&
         & trans,work,wv,info,init,initu)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & mld_z_base_solver_type, psb_ipk_
      implicit none 
      type(psb_desc_type), intent(in)              :: desc_data
      class(mld_z_base_solver_type), intent(inout) :: sv
      type(psb_z_vect_type),intent(inout)          :: x
      type(psb_z_vect_type),intent(inout)          :: y
      complex(psb_dpk_),intent(in)                     :: alpha,beta
      character(len=1),intent(in)                    :: trans
      complex(psb_dpk_),target, intent(inout)          :: work(:)
      type(psb_z_vect_type),intent(inout)            :: wv(:)
      integer(psb_ipk_), intent(out)                 :: info
      character, intent(in), optional                :: init
      type(psb_z_vect_type),intent(inout), optional   :: initu
    end subroutine mld_z_base_solver_apply_vect
  end interface
  
  interface 
    subroutine mld_z_base_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
       & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
       & mld_z_base_solver_type, psb_ipk_, psb_i_base_vect_type      
      Implicit None
      
      ! Arguments
      type(psb_zspmat_type), intent(in), target             :: a
      Type(psb_desc_type), Intent(in)                       :: desc_a 
      class(mld_z_base_solver_type), intent(inout)          :: sv
      integer(psb_ipk_), intent(out)                        :: info
      type(psb_zspmat_type), intent(in), target, optional   :: b
      class(psb_z_base_sparse_mat), intent(in), optional    :: amold
      class(psb_z_base_vect_type), intent(in), optional     :: vmold
      class(psb_i_base_vect_type), intent(in), optional     :: imold
    end subroutine mld_z_base_solver_bld
  end interface
  
  interface 
    subroutine mld_z_base_solver_cnv(sv,info,amold,vmold,imold)
      import :: psb_z_base_sparse_mat, psb_z_base_vect_type, psb_dpk_, &
       & mld_z_base_solver_type, psb_ipk_, psb_i_base_vect_type      
      Implicit None
      
      ! Arguments
      class(mld_z_base_solver_type), intent(inout)          :: sv
      integer(psb_ipk_), intent(out)                        :: info
      class(psb_z_base_sparse_mat), intent(in), optional    :: amold
      class(psb_z_base_vect_type), intent(in), optional     :: vmold
      class(psb_i_base_vect_type), intent(in), optional     :: imold
    end subroutine mld_z_base_solver_cnv
  end interface
  
  interface 
    subroutine mld_z_base_solver_check(sv,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & mld_z_base_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(mld_z_base_solver_type), intent(inout)   :: sv
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine mld_z_base_solver_check
  end interface
  
  interface 
    subroutine mld_z_base_solver_seti(sv,what,val,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & mld_z_base_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(mld_z_base_solver_type), intent(inout)   :: sv 
      integer(psb_ipk_), intent(in)                  :: what 
      integer(psb_ipk_), intent(in)                  :: val
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine mld_z_base_solver_seti
  end interface
  
  interface 
    subroutine mld_z_base_solver_setc(sv,what,val,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, & 
           & mld_z_base_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(mld_z_base_solver_type), intent(inout)   :: sv
      integer(psb_ipk_), intent(in)                  :: what 
      character(len=*), intent(in)                   :: val
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine mld_z_base_solver_setc
  end interface 
  
  interface 
    subroutine mld_z_base_solver_setr(sv,what,val,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & mld_z_base_solver_type, psb_ipk_            
      Implicit None      
      ! Arguments
      class(mld_z_base_solver_type), intent(inout)   :: sv 
      integer(psb_ipk_), intent(in)                  :: what 
      real(psb_dpk_), intent(in)                      :: val
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine mld_z_base_solver_setr
  end interface 
  
  interface 
    subroutine mld_z_base_solver_cseti(sv,what,val,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & mld_z_base_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(mld_z_base_solver_type), intent(inout) :: sv 
      character(len=*), intent(in)                   :: what 
      integer(psb_ipk_), intent(in)                  :: val
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine mld_z_base_solver_cseti
  end interface
  
  interface 
    subroutine mld_z_base_solver_csetc(sv,what,val,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, & 
           & mld_z_base_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(mld_z_base_solver_type), intent(inout) :: sv
      character(len=*), intent(in)                   :: what 
      character(len=*), intent(in)                   :: val
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine mld_z_base_solver_csetc
  end interface 
  
  interface 
    subroutine mld_z_base_solver_csetr(sv,what,val,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & mld_z_base_solver_type, psb_ipk_            
      Implicit None      
      ! Arguments
      class(mld_z_base_solver_type), intent(inout) :: sv 
      character(len=*), intent(in)                   :: what 
      real(psb_dpk_), intent(in)                      :: val
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine mld_z_base_solver_csetr
  end interface
  
  interface
    subroutine mld_z_base_solver_free(sv,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & mld_z_base_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(mld_z_base_solver_type), intent(inout) :: sv
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine mld_z_base_solver_free
  end interface
  
  interface
    subroutine mld_z_base_solver_descr(sv,info,iout,coarse)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & mld_z_base_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(mld_z_base_solver_type), intent(in) :: sv
      integer(psb_ipk_), intent(out)              :: info
      integer(psb_ipk_), intent(in), optional     :: iout
      logical, intent(in), optional               :: coarse

    end subroutine mld_z_base_solver_descr
  end interface 
  
  interface 
    subroutine mld_z_base_solver_dmp(sv,ictxt,level,info,prefix,head,solver)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & mld_z_base_solver_type, psb_ipk_      
      implicit none 
      class(mld_z_base_solver_type), intent(in) :: sv
      integer(psb_ipk_), intent(in)               :: ictxt
      integer(psb_ipk_), intent(in)               :: level
      integer(psb_ipk_), intent(out)              :: info
      character(len=*), intent(in), optional :: prefix, head
      logical, optional, intent(in)    :: solver
    end subroutine mld_z_base_solver_dmp
  end interface
 
   
  interface
    subroutine mld_z_base_solver_clone(sv,svout,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & mld_z_base_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(mld_z_base_solver_type), intent(inout)              :: sv
      class(mld_z_base_solver_type), allocatable, intent(inout) :: svout
      integer(psb_ipk_), intent(out)               :: info
    end subroutine mld_z_base_solver_clone
  end interface



contains
  !
  ! Function returning the size of the data structure
  ! in bytes or in number of nonzeros of the operator(s) involved. 
  !

  function z_base_solver_sizeof(sv) result(val)
    implicit none 
    ! Arguments
    class(mld_z_base_solver_type), intent(in) :: sv
    integer(psb_epk_)                  :: val
    integer(psb_ipk_)             :: i
    val = 0

    return
  end function z_base_solver_sizeof

  function z_base_solver_get_nzeros(sv) result(val)
    implicit none 
    class(mld_z_base_solver_type), intent(in) :: sv
    integer(psb_epk_) :: val
    integer(psb_ipk_)             :: i
    val = 0
  end function z_base_solver_get_nzeros

  subroutine z_base_solver_default(sv) 
    implicit none 
    ! Arguments
    class(mld_z_base_solver_type), intent(inout) :: sv
    ! Do nothing for base version

    return
  end subroutine z_base_solver_default

  function z_base_solver_get_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "Base solver"
  end function z_base_solver_get_fmt

  !
  ! If this is true, then the solver needs a starting
  ! guess. Currently only handled in JAC smoother. 
  ! 
  function z_base_solver_is_iterative() result(val)
    implicit none 
    logical  :: val

    val = .false.
  end function z_base_solver_is_iterative

  function z_base_solver_get_id() result(val)
    implicit none 
    integer(psb_ipk_)  :: val

    val = mld_f_none_
  end function z_base_solver_get_id

  function z_base_solver_get_wrksize() result(val)
    implicit none 
    integer(psb_ipk_)  :: val

    val = 0
  end function z_base_solver_get_wrksize

end module mld_z_base_solver_mod
