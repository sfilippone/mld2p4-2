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
! File: mld_c_as_smoother_mod.f90
!
! Module: mld_c_as_smoother_mod
!
!  This module defines: 
!    the mld_c_as_smoother_type data structure containing the
!    smoother for an Additive Schwarz smoother.
!
!  To begin with, the build procedure constructs the extended
!  matrix A and its corresponding descriptor (this has multiple
!  halo layers duplicated across different processes); it then 
!  stores in ND the block off-diagonal matrix, and builds the solver
!  on the (extended) block diagonal matrix.
!
!  The code allows for the variations of Additive Schwartz, Restricted
!  Additive Schwartz and Additive Schwartz with Harmonic Extensions.
!  From an implementation point of view, these are handled by
!  combining application/non-application of the prolongator/restrictor
!  operators.
!  
module mld_c_as_smoother

  use mld_c_base_smoother_mod
  
  type, extends(mld_c_base_smoother_type) :: mld_c_as_smoother_type
    ! The local solver component is inherited from the
    ! parent type. 
    !    class(mld_c_base_solver_type), allocatable :: sv
    !    
    type(psb_cspmat_type) :: nd
    type(psb_desc_type)     :: desc_data 
    integer(psb_ipk_)       :: novr, restr, prol, nd_nnz_tot
  contains
    procedure, pass(sm) :: check   => mld_c_as_smoother_check
    procedure, pass(sm) :: dump    => mld_c_as_smoother_dmp
    procedure, pass(sm) :: build   => mld_c_as_smoother_bld
    procedure, pass(sm) :: cnv     => mld_c_as_smoother_cnv
    procedure, pass(sm) :: clone   => mld_c_as_smoother_clone
    procedure, pass(sm) :: apply_v => mld_c_as_smoother_apply_vect
    procedure, pass(sm) :: apply_a => mld_c_as_smoother_apply
    procedure, pass(sm) :: restr_a => mld_c_as_smoother_restr_a
    procedure, pass(sm) :: prol_a  => mld_c_as_smoother_prol_a
    procedure, pass(sm) :: restr_v => mld_c_as_smoother_restr_v
    procedure, pass(sm) :: prol_v  => mld_c_as_smoother_prol_v
    generic, public     :: apply_restr   => restr_v, restr_a
    generic, public     :: apply_prol    => prol_v, prol_a
    procedure, pass(sm) :: free    => mld_c_as_smoother_free
    procedure, pass(sm) :: cseti   => mld_c_as_smoother_cseti
    procedure, pass(sm) :: csetc   => mld_c_as_smoother_csetc
    procedure, pass(sm) :: descr   => c_as_smoother_descr
    procedure, pass(sm) :: sizeof  => c_as_smoother_sizeof
    procedure, pass(sm) :: default => c_as_smoother_default
    procedure, pass(sm) :: get_nzeros => c_as_smoother_get_nzeros
    procedure, pass(sm) :: get_wrksz => c_as_smoother_get_wrksize
    procedure, nopass   :: get_fmt    => c_as_smoother_get_fmt
    procedure, nopass   :: get_id     => c_as_smoother_get_id
  end type mld_c_as_smoother_type
  
  
  private :: c_as_smoother_descr,  c_as_smoother_sizeof, &
       &  c_as_smoother_default, c_as_smoother_get_nzeros, &
       &  c_as_smoother_get_fmt, c_as_smoother_get_id, &
       &  c_as_smoother_get_wrksize

  character(len=6), parameter, private :: &
       &  restrict_names(0:4)=(/'none ','halo ','     ','     ','     '/)
  character(len=12), parameter, private :: &
       &  prolong_names(0:3)=(/'none       ','sum        ','average    ','square root'/)


  interface 
    subroutine mld_c_as_smoother_check(sm,info)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_spk_, mld_c_as_smoother_type, psb_epk_, psb_desc_type, psb_ipk_
      implicit none 
      class(mld_c_as_smoother_type), intent(inout) :: sm 
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine mld_c_as_smoother_check
  end interface
  
  interface 
    subroutine mld_c_as_smoother_restr_v(sm,x,trans,work,info,data)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_spk_, mld_c_as_smoother_type, psb_epk_, &
           & psb_desc_type, psb_ipk_
      implicit none 
      class(mld_c_as_smoother_type), intent(inout) :: sm
      type(psb_c_vect_type),intent(inout)          :: x
      character(len=1),intent(in)                    :: trans
      complex(psb_spk_),target, intent(inout)          :: work(:)
      integer(psb_ipk_), intent(out)                 :: info
      integer(psb_ipk_), optional, intent(in)        :: data
    end subroutine mld_c_as_smoother_restr_v
  end interface
  
  interface 
    subroutine mld_c_as_smoother_restr_a(sm,x,trans,work,info,data)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_spk_, mld_c_as_smoother_type, psb_epk_, &
           & psb_desc_type, psb_ipk_
      implicit none 
      class(mld_c_as_smoother_type), intent(inout) :: sm
      complex(psb_spk_), intent(inout)                 :: x(:)
      character(len=1),intent(in)                    :: trans
      complex(psb_spk_),target, intent(inout)          :: work(:)
      integer(psb_ipk_), intent(out)                 :: info
      integer(psb_ipk_), optional, intent(in)        :: data
    end subroutine mld_c_as_smoother_restr_a
  end interface

  interface 
    subroutine mld_c_as_smoother_prol_v(sm,x,trans,work,info,data)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_spk_, mld_c_as_smoother_type, psb_epk_, &
           & psb_desc_type, psb_ipk_
      implicit none 
      class(mld_c_as_smoother_type), intent(inout) :: sm
      type(psb_c_vect_type),intent(inout)          :: x
      character(len=1),intent(in)                    :: trans
      complex(psb_spk_),target, intent(inout)          :: work(:)
      integer(psb_ipk_), intent(out)                 :: info
      integer(psb_ipk_), optional, intent(in)        :: data
    end subroutine mld_c_as_smoother_prol_v
  end interface
  
  interface 
    subroutine mld_c_as_smoother_prol_a(sm,x,trans,work,info,data)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_spk_, mld_c_as_smoother_type, psb_epk_, &
           & psb_desc_type, psb_ipk_
      implicit none 
      class(mld_c_as_smoother_type), intent(inout) :: sm
      complex(psb_spk_), intent(inout)                 :: x(:)
      character(len=1),intent(in)                    :: trans
      complex(psb_spk_),target, intent(inout)          :: work(:)
      integer(psb_ipk_), intent(out)                 :: info
      integer(psb_ipk_), optional, intent(in)        :: data
    end subroutine mld_c_as_smoother_prol_a
  end interface
  
  
  interface 
    subroutine mld_c_as_smoother_apply_vect(alpha,sm,x,beta,y,desc_data,&
      & trans,sweeps,work,wv,info,init,initu)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_spk_, mld_c_as_smoother_type, psb_epk_, &
           & psb_desc_type, psb_ipk_
      implicit none 
      type(psb_desc_type), intent(in)              :: desc_data
      class(mld_c_as_smoother_type), intent(inout) :: sm
      type(psb_c_vect_type),intent(inout)          :: x
      type(psb_c_vect_type),intent(inout)          :: y
      complex(psb_spk_),intent(in)                     :: alpha,beta
      character(len=1),intent(in)                    :: trans
      integer(psb_ipk_), intent(in)                  :: sweeps
      complex(psb_spk_),target, intent(inout)          :: work(:)
      type(psb_c_vect_type),intent(inout)          :: wv(:)
      integer(psb_ipk_), intent(out)                 :: info
      character, intent(in), optional                :: init
      type(psb_c_vect_type),intent(inout), optional   :: initu
    end subroutine mld_c_as_smoother_apply_vect
  end interface
  
  interface
    subroutine mld_c_as_smoother_apply(alpha,sm,x,beta,y,desc_data,& 
         & trans,sweeps,work,info,init,initu)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_spk_, mld_c_as_smoother_type, psb_epk_,&
           & psb_desc_type, psb_ipk_
      implicit none 
      type(psb_desc_type), intent(in)      :: desc_data
      class(mld_c_as_smoother_type), intent(inout) :: sm
      complex(psb_spk_),intent(inout)         :: x(:)
      complex(psb_spk_),intent(inout)         :: y(:)
      complex(psb_spk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)           :: trans
      integer(psb_ipk_), intent(in)         :: sweeps
      complex(psb_spk_),target, intent(inout) :: work(:)
      integer(psb_ipk_), intent(out)        :: info
      character, intent(in), optional       :: init
      complex(psb_spk_),intent(inout), optional :: initu(:)
    end subroutine mld_c_as_smoother_apply
  end interface
  
  interface
    subroutine mld_c_as_smoother_bld(a,desc_a,sm,info,amold,vmold,imold)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_spk_, mld_c_as_smoother_type, psb_epk_, &
           & psb_desc_type, psb_c_base_sparse_mat, psb_ipk_,&
           & psb_i_base_vect_type
      implicit none 
      type(psb_cspmat_type), intent(in), target        :: a
      Type(psb_desc_type), Intent(inout)                 :: desc_a 
      class(mld_c_as_smoother_type), intent(inout)       :: sm
      integer(psb_ipk_), intent(out)                     :: info
      class(psb_c_base_sparse_mat), intent(in), optional :: amold
      class(psb_c_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine mld_c_as_smoother_bld
  end interface
  
  interface
    subroutine mld_c_as_smoother_cnv(sm,info,amold,vmold,imold)
      import :: psb_c_base_vect_type, &
           & psb_spk_, mld_c_as_smoother_type, psb_epk_, &
           & psb_c_base_sparse_mat, psb_ipk_, psb_i_base_vect_type
      implicit none 
      class(mld_c_as_smoother_type), intent(inout)       :: sm
      integer(psb_ipk_), intent(out)                     :: info
      class(psb_c_base_sparse_mat), intent(in), optional :: amold
      class(psb_c_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine mld_c_as_smoother_cnv
  end interface
  
  interface 
    subroutine mld_c_as_smoother_cseti(sm,what,val,info,idx)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_spk_, mld_c_as_smoother_type, psb_epk_, psb_desc_type, psb_ipk_
      implicit none 
      class(mld_c_as_smoother_type), intent(inout) :: sm 
      character(len=*), intent(in)                   :: what 
      integer(psb_ipk_), intent(in)                  :: val
      integer(psb_ipk_), intent(out)                 :: info
      integer(psb_ipk_), intent(in), optional        :: idx
    end subroutine mld_c_as_smoother_cseti
  end interface
  
  interface 
    subroutine mld_c_as_smoother_csetc(sm,what,val,info,idx)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_spk_, mld_c_as_smoother_type, psb_epk_, psb_desc_type, psb_ipk_
      implicit none 
      class(mld_c_as_smoother_type), intent(inout) :: sm
      character(len=*), intent(in)                   :: what 
      character(len=*), intent(in)                   :: val
      integer(psb_ipk_), intent(out)                 :: info
      integer(psb_ipk_), intent(in), optional        :: idx
    end subroutine mld_c_as_smoother_csetc
  end interface
  
  interface 
    subroutine mld_c_as_smoother_free(sm,info)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_spk_, mld_c_as_smoother_type, psb_epk_, psb_desc_type, psb_ipk_
      implicit none 
      class(mld_c_as_smoother_type), intent(inout) :: sm
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine mld_c_as_smoother_free
  end interface
  
  interface 
    subroutine mld_c_as_smoother_dmp(sm,ictxt,level,info,prefix,head,smoother,solver)
      import :: psb_cspmat_type, psb_c_vect_type, psb_c_base_vect_type, &
           & psb_spk_, mld_c_as_smoother_type, psb_epk_, psb_desc_type, &
           & psb_ipk_
      implicit none 
      class(mld_c_as_smoother_type), intent(in) :: sm
      integer(psb_ipk_), intent(in)               :: ictxt
      integer(psb_ipk_), intent(in)               :: level
      integer(psb_ipk_), intent(out)              :: info
      character(len=*), intent(in), optional :: prefix, head
      logical, optional, intent(in)    :: smoother, solver
    end subroutine mld_c_as_smoother_dmp
  end interface
  
  interface 
    subroutine mld_c_as_smoother_clone(sm,smout,info)
      import :: mld_c_as_smoother_type, psb_spk_, &
           & mld_c_base_smoother_type, psb_ipk_
      class(mld_c_as_smoother_type), intent(inout)                :: sm
      class(mld_c_base_smoother_type), allocatable, intent(inout) :: smout
      integer(psb_ipk_), intent(out)               :: info
    end subroutine mld_c_as_smoother_clone
  end interface
  
contains

  function c_as_smoother_sizeof(sm) result(val)
    implicit none 
    ! Arguments
    class(mld_c_as_smoother_type), intent(in) :: sm
    integer(psb_epk_) :: val
    integer(psb_ipk_)             :: i

    val = psb_sizeof_ip 
    if (allocated(sm%sv)) val = val + sm%sv%sizeof()
    val = val + sm%nd%sizeof()

    return
  end function c_as_smoother_sizeof

  function c_as_smoother_get_nzeros(sm) result(val)
    implicit none 
    class(mld_c_as_smoother_type), intent(in) :: sm
    integer(psb_epk_) :: val
    integer(psb_ipk_)             :: i
    val = 0
    if (allocated(sm%sv)) &
         &  val =  sm%sv%get_nzeros()
    val = val + sm%nd%get_nzeros()

  end function c_as_smoother_get_nzeros

  subroutine c_as_smoother_default(sm)

    use psb_base_mod, only : psb_halo_, psb_none_

    Implicit None

    ! Arguments
    class(mld_c_as_smoother_type), intent(inout) :: sm 

    !
    ! Default: AS with 1 overlap layer
    ! 
    sm%restr = psb_halo_
    sm%prol  = psb_sum_
    sm%novr  = 1

    if (allocated(sm%sv)) then 
      call sm%sv%default()
    end if

    return
  end subroutine c_as_smoother_default


  subroutine c_as_smoother_descr(sm,info,iout,coarse)

    Implicit None

    ! Arguments
    class(mld_c_as_smoother_type), intent(in) :: sm
    integer(psb_ipk_), intent(out)                      :: info
    integer(psb_ipk_), intent(in), optional             :: iout
    logical, intent(in), optional             :: coarse

    ! Local variables
    integer(psb_ipk_)      :: err_act
    character(len=20), parameter :: name='mld_c_as_smoother_descr'
    integer(psb_ipk_) :: iout_
    logical      :: coarse_

    call psb_erractionsave(err_act)
    info = psb_success_
    if (present(coarse)) then 
      coarse_ = coarse
    else
      coarse_ = .false.
    end if
    if (present(iout)) then 
      iout_ = iout 
    else
      iout_ = psb_out_unit
    endif

    if (.not.coarse_) then 
      write(iout_,*) '  Additive Schwarz with  ',&
           &  sm%novr, ' overlap layers.'
      write(iout_,*) '  Restrictor:  ',restrict_names(sm%restr)
      write(iout_,*) '  Prolongator: ',prolong_names(sm%prol)
      write(iout_,*) '  Local solver:'
    endif
    if (allocated(sm%sv)) then 
      call sm%sv%descr(info,iout_,coarse=coarse)
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine c_as_smoother_descr

  function c_as_smoother_get_wrksize(sm) result(val)
    implicit none 
    class(mld_c_as_smoother_type), intent(inout) :: sm
    integer(psb_ipk_)  :: val

    val = 3
    if (allocated(sm%sv)) val = val + sm%sv%get_wrksz()
    
  end function c_as_smoother_get_wrksize
  
  function c_as_smoother_get_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "Additive Schwarz"
  end function c_as_smoother_get_fmt

  function c_as_smoother_get_id() result(val)
    implicit none 
    integer(psb_ipk_)  :: val
    
    val = mld_as_
  end function c_as_smoother_get_id

end module mld_c_as_smoother
