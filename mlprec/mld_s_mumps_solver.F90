!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010,2012,2013
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
!
!
!

module mld_s_mumps_solver
#if defined(HAVE_MUMPS_)
  use smumps_struc_def
#endif
  use mld_s_base_solver_mod
  
#if defined(LONG_INTEGERS)
  
  type, extends(mld_s_base_solver_type) :: mld_s_mumps_solver_type
    
  end type mld_s_mumps_solver_type
#else
  type, extends(mld_s_base_solver_type) :: mld_s_mumps_solver_type
#if defined(HAVE_MUMPS_)
    type(smumps_struc), allocatable  :: id
#else
    integer, allocatable :: id
#endif
    integer(psb_ipk_),dimension(2)   :: ipar
    logical                          :: built=.false.
  contains
    procedure, pass(sv) :: build   => s_mumps_solver_bld
    procedure, pass(sv) :: apply_a => s_mumps_solver_apply
    procedure, pass(sv) :: apply_v => s_mumps_solver_apply_vect
    procedure, pass(sv) :: free    => s_mumps_solver_free
    procedure, pass(sv) :: descr   => s_mumps_solver_descr
    procedure, pass(sv) :: sizeof  => s_mumps_solver_sizeof
    procedure, pass(sv) :: seti    => s_mumps_solver_seti
    procedure, pass(sv) :: setr    => s_mumps_solver_setr
    procedure, pass(sv) :: cseti    =>s_mumps_solver_cseti
    procedure, pass(sv) :: csetr    => s_mumps_solver_csetr
    procedure, pass(sv) :: default  => s_mumps_solver_default
#if defined(HAVE_FINAL) 

    final               :: s_mumps_solver_finalize
#endif
  end type mld_s_mumps_solver_type


  private :: s_mumps_solver_bld, s_mumps_solver_apply, &
       &  s_mumps_solver_free,   s_mumps_solver_descr, &
       &  s_mumps_solver_sizeof, s_mumps_solver_apply_vect,&
       &  s_mumps_solver_seti,   s_mumps_solver_setr,    &
       &  s_mumps_solver_cseti, s_mumps_solver_csetri,   &
       &  s_mumps_solver_default
#if defined(HAVE_FINAL) 
  private :: s_mumps_solver_finalize
#endif

  interface 
    subroutine s_mumps_solver_apply_vect(alpha,sv,x,beta,y,desc_data,trans,work,info)
      import :: psb_desc_type, mld_s_mumps_solver_type, psb_s_vect_type, psb_dpk_, psb_spk_, &
           & psb_sspmat_type, psb_s_base_sparse_mat, psb_s_base_vect_type, psb_ipk_
      implicit none 
      type(psb_desc_type), intent(in)      :: desc_data
      class(mld_s_mumps_solver_type), intent(inout) :: sv
      type(psb_s_vect_type),intent(inout)  :: x
      type(psb_s_vect_type),intent(inout)  :: y
      real(psb_spk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)          :: trans
      real(psb_spk_),target, intent(inout) :: work(:)
      integer, intent(out)                 :: info

      integer(psb_ipk_)    :: err_act
      character(len=20)  :: name='s_mumps_solver_apply_vect'
    end subroutine s_mumps_solver_apply_vect
  end interface

  interface
    subroutine s_mumps_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info)
      import :: psb_desc_type, mld_s_mumps_solver_type, psb_s_vect_type, psb_dpk_, psb_spk_, &
           & psb_sspmat_type, psb_s_base_sparse_mat, psb_s_base_vect_type, psb_ipk_
      implicit none 
      type(psb_desc_type), intent(in)      :: desc_data
      class(mld_s_mumps_solver_type), intent(inout) :: sv
      real(psb_spk_),intent(inout)         :: x(:)
      real(psb_spk_),intent(inout)         :: y(:)
      real(psb_spk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)          :: trans
      real(psb_spk_),target, intent(inout) :: work(:)
      integer(psb_ipk_), intent(out)                 :: info

      integer(psb_ipk_)    :: n_row, n_col, nglob
      real(psb_spk_), pointer     :: ww(:)
      real(psb_spk_), allocatable, target :: gx(:)
      integer(psb_ipk_)  :: ictxt,np,me,i, err_act
      character          :: trans_
      character(len=20)  :: name='s_mumps_solver_apply'
    end subroutine s_mumps_solver_apply
  end interface

  interface
    subroutine s_mumps_solver_bld(a,desc_a,sv,upd,info,b,amold,vmold,imold)

      use mpi    
      import :: psb_desc_type, mld_s_mumps_solver_type, psb_s_vect_type, psb_spk_, &
           & psb_sspmat_type, psb_s_base_sparse_mat, psb_s_base_vect_type,&
           & psb_ipk_, psb_i_base_vect_type

      Implicit None

      ! Arguments
      type(psb_sspmat_type), intent(in), target           :: a
      Type(psb_desc_type), Intent(in)                     :: desc_a 
      class(mld_s_mumps_solver_type), intent(inout)       :: sv
      character, intent(in)                               :: upd
      integer(psb_ipk_), intent(out)                      :: info
      type(psb_sspmat_type), intent(in), target, optional :: b
      class(psb_s_base_sparse_mat), intent(in), optional  :: amold
      class(psb_s_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional   :: imold
    end subroutine s_mumps_solver_bld
  end interface

contains

  subroutine s_mumps_solver_free(sv,info)

    Implicit None

    ! Arguments
    class(mld_s_mumps_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(out)                :: info
    Integer(psb_ipk_) :: err_act
    character(len=20)  :: name='s_mumps_solver_free'

    call psb_erractionsave(err_act)
#if defined(HAVE_MUMPS_)
    if (allocated(sv%id)) then      
      if (sv%built) then 
        sv%id%job = -2
        call smumps(sv%id)
        info = sv%id%infog(1)
        if (info /= psb_success_) goto 9999
      end if
      deallocate(sv%id)
      sv%built=.false.
    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
#endif
  end subroutine s_mumps_solver_free

#if defined(HAVE_FINAL)
  subroutine s_mumps_solver_finalize(sv)

    Implicit None

    ! Arguments
    type(mld_s_mumps_solver_type), intent(inout) :: sv 
    integer :: info
    Integer :: err_act
    character(len=20)  :: name='s_mumps_solver_finalize'

    call sv%free(info) 

    return

  end subroutine s_mumps_solver_finalize
#endif

  subroutine s_mumps_solver_descr(sv,info,iout,coarse)

    Implicit None

    ! Arguments
    class(mld_s_mumps_solver_type), intent(in) :: sv
    integer(psb_ipk_), intent(out)                     :: info
    integer(psb_ipk_), intent(in), optional            :: iout
    logical, intent(in), optional       :: coarse

    ! Local variables
    integer(psb_ipk_)      :: err_act
    integer(psb_ipk_)      :: ictxt, me, np
    character(len=20), parameter :: name='mld_z_mumps_solver_descr'
    integer(psb_ipk_) :: iout_

    call psb_erractionsave(err_act)
    info = psb_success_
    if (present(iout)) then 
      iout_ = iout 
    else
      iout_ = 6
    endif

    write(iout_,*) '  MUMPS  Solver. '

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine s_mumps_solver_descr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$ WARNING: OTHERS PARAMETERS OF MUMPS COULD BE ADDED. FOR THIS, ADD AN     !!$
!!$ INTEGER IN MLD_BASE_PREC_TYPE.F90 AND MODIFY SUBROUTINE SET              !!$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine s_mumps_solver_seti(sv,what,val,info)

    Implicit None

    ! Arguments
    class(mld_s_mumps_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(in)                 :: what
    integer(psb_ipk_), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_)  :: err_act
    character(len=20)  :: name='z_mumps_solver_seti'

    info = psb_success_
    call psb_erractionsave(err_act)
    select case(what)
#if defined(HAVE_MUMPS_)
    case(mld_as_sequential_)   
      sv%ipar(1)=val
    case(mld_mumps_print_err_)
      sv%ipar(2)=val
    !case(mld_print_stat_)
    !  sv%id%icntl(2)=val
    !  sv%ipar(2)=val
    !case(mld_print_glob_)
    !  sv%id%icntl(3)=val
    !  sv%ipar(3)=val
#endif
    case default
      call sv%mld_s_base_solver_type%set(what,val,info)
    end select

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine s_mumps_solver_seti


  subroutine s_mumps_solver_setr(sv,what,val,info)

    Implicit None

    ! Arguments
    class(mld_s_mumps_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(in)                 :: what
    real(psb_spk_), intent(in)                    :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_)  :: err_act
    character(len=20)  :: name='z_mumps_solver_setr'

    info = psb_success_
    call psb_erractionsave(err_act)

    select case(what)
    case default
      call sv%mld_s_base_solver_type%set(what,val,info)
    end select

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine s_mumps_solver_setr

  subroutine s_mumps_solver_cseti(sv,what,val,info)

    Implicit None

    ! Arguments
    class(mld_s_mumps_solver_type), intent(inout) :: sv
    character(len=*), intent(in)                  :: what
    integer(psb_ipk_), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_)  :: err_act, iwhat
    character(len=20)  :: name='s_mumps_solver_cseti'

    info = psb_success_
    call psb_erractionsave(err_act)

    select case(psb_toupper(what))
#if defined(HAVE_MUMPS_)
    case('SET_AS_SEQUENTIAL')
      iwhat=mld_as_sequential_
    case('SET_MUMPS_PRINT_ERR')
      iwhat=mld_mumps_print_err_
#endif
    case default
      iwhat=-1
    end select

    if (iwhat >=0 ) then 
      call sv%set(iwhat,val,info)
    else
      call sv%mld_s_base_solver_type%set(what,val,info)
    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine s_mumps_solver_cseti

  subroutine s_mumps_solver_csetr(sv,what,val,info)

    Implicit None

    ! Arguments
    class(mld_s_mumps_solver_type), intent(inout) :: sv
    character(len=*), intent(in)                  :: what
    real(psb_spk_), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_)  :: err_act, iwhat
    character(len=20)  :: name='s_mumps_solver_csetr'

    info = psb_success_
    call psb_erractionsave(err_act)

    select case(psb_toupper(what))
    case default
      call sv%mld_s_base_solver_type%set(what,val,info)
    end select

    if (iwhat >=0 ) then 
      call sv%set(iwhat,val,info)
    else
      call sv%mld_s_base_solver_type%set(what,val,info)
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine s_mumps_solver_csetr

  !!NOTE: BY DEFAULT BLR is activated with a dropping parameter to 1d-4       !!
  subroutine s_mumps_solver_default(sv)

    Implicit none

    !Argument
    class(mld_s_mumps_solver_type),intent(inout) :: sv
    integer(psb_ipk_) :: info
    integer(psb_ipk_)  :: err_act,ictx,icomm
    character(len=20)  :: name='s_mumps_default'

    info = psb_success_
    call psb_erractionsave(err_act)

#if defined(HAVE_MUMPS_)
    if (.not.allocated(sv%id)) then 
      allocate(sv%id,stat=info)
      if (info /= psb_success_) then
        info=psb_err_alloc_dealloc_
        call psb_errpush(info,name,a_err='mld_smumps_default')
        goto 9999
      end if
      sv%built=.false.
    end if

    ! INSTANTIATION OF sv%id needed to set parmater but mpi communicator needed
    ! sv%id%job = -1
    ! sv%id%par=1
    ! call dmumps(sv%id)    
    sv%ipar(1)=2
    !sv%ipar(10)=6
    !sv%ipar(11)=0
    !sv%ipar(12)=6

#endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine s_mumps_solver_default

  function s_mumps_solver_sizeof(sv) result(val)

    implicit none 
    ! Arguments
    class(mld_s_mumps_solver_type), intent(in) :: sv
    integer(psb_long_int_k_) :: val
    integer             :: i
#if defined(HAVE_MUMPS_)
    val = (sv%id%INFOG(22)+sv%id%INFOG(32))*1d+6
#else
    val = 0 
#endif
    ! val = 2*psb_sizeof_int + psb_sizeof_dp
    ! val = val + sv%symbsize
    ! val = val + sv%numsize
    return
  end function s_mumps_solver_sizeof
#endif
end module mld_s_mumps_solver

