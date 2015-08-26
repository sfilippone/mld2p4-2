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
! File: mld_z_prec_type.f90
!
! Module: mld_z_prec_type
!
!  This module defines: 
!  - the mld_z_prec_type data structure containing the preconditioner and related
!    data structures;
!
!  It contains routines for
!  - Building and applying; 
!  - checking if the preconditioner is correctly defined;
!  - printing a	description of the preconditioner;
!  - deallocating the preconditioner data structure.  
!

module mld_z_prec_type

  use mld_base_prec_type
  use mld_z_base_solver_mod
  use mld_z_base_smoother_mod
  use mld_z_onelev_mod
  use psb_prec_mod, only : psb_zprec_type

  !
  ! Type: mld_Tprec_type.
  !
  !  This is the data type containing all the information about the multilevel
  !  preconditioner (here and in the following 'T' denotes 'd', 's', 'c' and
  !  'z', according to the real/complex, single/double precision version of
  !  MLD2P4). It consists of an array of 'one-level' intermediate data structures
  !  of type mld_Tonelev_type, each containing the information needed to apply
  !  the smoothing and the coarse-space correction at a generic level. RT is the
  !  real data type, i.e. S for both S and C, and D for both D and Z. 
  !
  !  type mld_Tprec_type
  !    type(mld_Tonelev_type), allocatable :: precv(:) 
  !  end type mld_Tprec_type
  ! 
  !  Note that the levels are numbered in increasing order starting from
  !  the finest one and the number of levels is given by size(precv(:)).
  !
  !

  type, extends(psb_zprec_type)        :: mld_zprec_type
    integer(psb_ipk_)                  :: ictxt
    integer(psb_ipk_)                  :: coarse_aggr_size
    real(psb_dpk_)                        :: op_complexity=dzero
    type(mld_z_onelev_type), allocatable :: precv(:) 
  contains
    procedure, pass(prec)               :: psb_z_apply2_vect => mld_z_apply2_vect
    procedure, pass(prec)               :: psb_z_apply1_vect => mld_z_apply1_vect
    procedure, pass(prec)               :: psb_z_apply2v => mld_z_apply2v
    procedure, pass(prec)               :: psb_z_apply1v => mld_z_apply1v
    procedure, pass(prec)               :: dump           => mld_z_dump
    procedure, pass(prec)               :: clone          => mld_z_clone
    procedure, pass(prec)               :: free           => mld_z_prec_free
    procedure, pass(prec)               :: get_complexity => mld_z_get_compl
    procedure, pass(prec)               :: cmp_complexity => mld_z_cmp_compl
    procedure, pass(prec)               :: get_nzeros => mld_z_get_nzeros
    procedure, pass(prec)               :: sizeof => mld_zprec_sizeof
    procedure, pass(prec)               :: setsm  => mld_zprecsetsm
    procedure, pass(prec)               :: setsv  => mld_zprecsetsv
    procedure, pass(prec)               :: seti   => mld_zprecseti
    procedure, pass(prec)               :: setc   => mld_zprecsetc
    procedure, pass(prec)               :: setr   => mld_zprecsetr
    procedure, pass(prec)               :: cseti  => mld_zcprecseti
    procedure, pass(prec)               :: csetc  => mld_zcprecsetc
    procedure, pass(prec)               :: csetr  => mld_zcprecsetr
    generic, public                     :: set => seti, setc, setr, & 
         &       cseti, csetc, csetr, setsm, setsv 
    procedure, pass(prec)               :: get_smoother => mld_z_get_smootherp
    procedure, pass(prec)               :: get_solver   => mld_z_get_solverp
  end type mld_zprec_type

  private :: mld_z_dump, mld_z_get_compl,  mld_z_cmp_compl,&
       &  mld_z_get_nzeros


  !
  ! Interfaces to routines for checking the definition of the preconditioner,
  ! for printing its description and for deallocating its data structure
  !

  interface mld_precfree
    module procedure mld_zprecfree
  end interface


  interface mld_precdescr
    module procedure mld_zfile_prec_descr
  end interface

  interface mld_sizeof
    module procedure mld_zprec_sizeof
  end interface

  interface mld_precaply
    subroutine mld_zprecaply2_vect(prec,x,y,desc_data,info,trans,work)
      import :: psb_zspmat_type, psb_desc_type, &
           & psb_dpk_, psb_z_vect_type, mld_zprec_type, psb_ipk_
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_zprec_type), intent(inout) :: prec
      type(psb_z_vect_type),intent(inout) :: x
      type(psb_z_vect_type),intent(inout) :: y
      integer(psb_ipk_), intent(out)                :: info
      character(len=1), optional          :: trans
      complex(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine mld_zprecaply2_vect
    subroutine mld_zprecaply1_vect(prec,x,desc_data,info,trans,work)
      import :: psb_zspmat_type, psb_desc_type, &
           & psb_dpk_, psb_z_vect_type, mld_zprec_type, psb_ipk_
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_zprec_type), intent(inout) :: prec
      type(psb_z_vect_type),intent(inout) :: x
      integer(psb_ipk_), intent(out)                :: info
      character(len=1), optional          :: trans
      complex(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine mld_zprecaply1_vect
    subroutine mld_zprecaply(prec,x,y,desc_data,info,trans,work)
      import :: psb_zspmat_type, psb_desc_type, psb_dpk_, mld_zprec_type, psb_ipk_
      type(psb_desc_type),intent(in)   :: desc_data
      type(mld_zprec_type), intent(inout) :: prec
      complex(psb_dpk_),intent(inout)     :: x(:)
      complex(psb_dpk_),intent(inout)     :: y(:)
      integer(psb_ipk_), intent(out)             :: info
      character(len=1), optional       :: trans
      complex(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine mld_zprecaply
    subroutine mld_zprecaply1(prec,x,desc_data,info,trans)
      import :: psb_zspmat_type, psb_desc_type, psb_dpk_, mld_zprec_type, psb_ipk_
      type(psb_desc_type),intent(in)   :: desc_data
      type(mld_zprec_type), intent(inout) :: prec
      complex(psb_dpk_),intent(inout)     :: x(:)
      integer(psb_ipk_), intent(out)             :: info
      character(len=1), optional       :: trans
    end subroutine mld_zprecaply1
  end interface

  interface 
    subroutine mld_zprecsetsm(prec,val,info,ilev)
      import :: psb_zspmat_type, psb_desc_type, psb_dpk_, &
           & mld_zprec_type, mld_z_base_smoother_type, psb_ipk_
      class(mld_zprec_type), intent(inout)        :: prec
      class(mld_z_base_smoother_type), intent(in) :: val
      integer(psb_ipk_), intent(out)                :: info
      integer(psb_ipk_), optional, intent(in)       :: ilev
    end subroutine mld_zprecsetsm
    subroutine mld_zprecsetsv(prec,val,info,ilev)
      import :: psb_zspmat_type, psb_desc_type, psb_dpk_, &
           & mld_zprec_type, mld_z_base_solver_type, psb_ipk_
      class(mld_zprec_type), intent(inout)      :: prec
      class(mld_z_base_solver_type), intent(in) :: val
      integer(psb_ipk_), intent(out)              :: info
      integer(psb_ipk_), optional, intent(in)     :: ilev
    end subroutine mld_zprecsetsv
    subroutine mld_zprecseti(prec,what,val,info,ilev)
      import :: psb_zspmat_type, psb_desc_type, psb_dpk_, &
           & mld_zprec_type, psb_ipk_
      class(mld_zprec_type), intent(inout)   :: prec
      integer(psb_ipk_), intent(in)            :: what 
      integer(psb_ipk_), intent(in)            :: val
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), optional, intent(in)  :: ilev
    end subroutine mld_zprecseti
    subroutine mld_zprecsetr(prec,what,val,info,ilev)
      import :: psb_zspmat_type, psb_desc_type, psb_dpk_, &
           & mld_zprec_type, psb_ipk_
      class(mld_zprec_type), intent(inout)   :: prec
      integer(psb_ipk_), intent(in)            :: what 
      real(psb_dpk_), intent(in)                :: val
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), optional, intent(in)  :: ilev
    end subroutine mld_zprecsetr
    subroutine mld_zprecsetc(prec,what,string,info,ilev)
      import :: psb_zspmat_type, psb_desc_type, psb_dpk_, &
           & mld_zprec_type, psb_ipk_
      class(mld_zprec_type), intent(inout)   :: prec
      integer(psb_ipk_), intent(in)            :: what 
      character(len=*), intent(in)             :: string
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), optional, intent(in)  :: ilev
    end subroutine mld_zprecsetc
    subroutine mld_zcprecseti(prec,what,val,info,ilev)
      import :: psb_zspmat_type, psb_desc_type, psb_dpk_, &
           & mld_zprec_type, psb_ipk_
      class(mld_zprec_type), intent(inout)   :: prec
      character(len=*), intent(in)             :: what 
      integer(psb_ipk_), intent(in)            :: val
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), optional, intent(in)  :: ilev
    end subroutine mld_zcprecseti
    subroutine mld_zcprecsetr(prec,what,val,info,ilev)
      import :: psb_zspmat_type, psb_desc_type, psb_dpk_, &
           & mld_zprec_type, psb_ipk_
      class(mld_zprec_type), intent(inout)   :: prec
      character(len=*), intent(in)             :: what 
      real(psb_dpk_), intent(in)                :: val
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), optional, intent(in)  :: ilev
    end subroutine mld_zcprecsetr
    subroutine mld_zcprecsetc(prec,what,string,info,ilev)
      import :: psb_zspmat_type, psb_desc_type, psb_dpk_, &
           & mld_zprec_type, psb_ipk_
      class(mld_zprec_type), intent(inout)   :: prec
      character(len=*), intent(in)             :: what 
      character(len=*), intent(in)             :: string
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), optional, intent(in)  :: ilev
    end subroutine mld_zcprecsetc
  end interface


  interface mld_move_alloc
    module procedure  mld_zprec_move_alloc
  end interface

contains
  !
  ! Function returning a pointer to the smoother
  !
  function mld_z_get_smootherp(prec,ilev) result(val)
    implicit none 
    class(mld_zprec_type), target, intent(in) :: prec
    integer(psb_ipk_), optional                 :: ilev
    class(mld_z_base_smoother_type), pointer  :: val
    integer(psb_ipk_)        :: ilev_
    
    val => null()
    if (present(ilev)) then 
      ilev_ = ilev
    else
      ! What is a good default? 
      ilev_ = 1
    end if
    if (allocated(prec%precv)) then 
      if ((1<=ilev_).and.(ilev_<=size(prec%precv))) then 
        if (allocated(prec%precv(ilev_)%sm)) then 
          val => prec%precv(ilev_)%sm
        end if
      end if
    end if
  end function mld_z_get_smootherp
  !
  ! Function returning a pointer to the solver
  !
  function mld_z_get_solverp(prec,ilev) result(val)
    implicit none 
    class(mld_zprec_type), target, intent(in) :: prec
    integer(psb_ipk_), optional                 :: ilev
    class(mld_z_base_solver_type), pointer  :: val
    integer(psb_ipk_)        :: ilev_
    
    val => null()
    if (present(ilev)) then 
      ilev_ = ilev
    else
      ! What is a good default? 
      ilev_ = 1
    end if
    if (allocated(prec%precv)) then 
      if ((1<=ilev_).and.(ilev_<=size(prec%precv))) then 
        if (allocated(prec%precv(ilev_)%sm)) then 
          if (allocated(prec%precv(ilev_)%sm%sv)) then 
            val => prec%precv(ilev_)%sm%sv
          end if
        end if
      end if
    end if
  end function mld_z_get_solverp
  !
  ! Function returning the size of the mld_prec_type data structure
  ! in bytes or in number of nonzeros of the operator(s) involved. 
  !
  function mld_z_get_nzeros(prec) result(val)
    implicit none 
    class(mld_zprec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer(psb_ipk_)        :: i
    val = 0
    if (allocated(prec%precv)) then 
      do i=1, size(prec%precv)
        val = val + prec%precv(i)%get_nzeros()
      end do
    end if
  end function mld_z_get_nzeros

  function mld_zprec_sizeof(prec) result(val)
    implicit none 
    class(mld_zprec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer(psb_ipk_)        :: i
    val = 0
    val = val + psb_sizeof_int
    if (allocated(prec%precv)) then 
      do i=1, size(prec%precv)
        val = val + prec%precv(i)%sizeof()
      end do
    end if
  end function mld_zprec_sizeof

  !
  ! Operator complexity: ratio of total number
  ! of nonzeros in the aggregated matrices at the
  ! various level to the nonzeroes at the fine level
  ! (original matrix)
  !
  
  function mld_z_get_compl(prec) result(val)
    implicit none 
    class(mld_zprec_type), intent(in) :: prec
    complex(psb_dpk_)  :: val
    
    val = prec%op_complexity

  end function mld_z_get_compl
  
  subroutine mld_z_cmp_compl(prec) 

    implicit none 
    class(mld_zprec_type), intent(inout) :: prec
    
    real(psb_dpk_) :: num,den
    integer(psb_ipk_) :: ictxt 
    integer(psb_ipk_)  :: il 

    num = -done
    den = done
    ictxt = prec%ictxt
    if (allocated(prec%precv)) then 
      il  = 1
      num = prec%precv(il)%base_a%get_nzeros()
      if (num >= dzero) then
        den = num 
        do il=2,size(prec%precv)
          num = num + max(0,prec%precv(il)%base_a%get_nzeros())
        end do
      end if
    end if
    call psb_min(ictxt,num) 
    if (num < dzero) then 
      den = done
    else
      call psb_sum(ictxt,num)
      call psb_sum(ictxt,den)
    end if
    prec%op_complexity = num/den
  end subroutine mld_z_cmp_compl
  
  !
  ! Subroutine: mld_file_prec_descr
  ! Version: complex
  !
  !  This routine prints a description of the preconditioner to the standard 
  !  output or to a file. It must be called after the preconditioner has been
  !  built by mld_precbld.
  !
  ! Arguments:
  !  p       -  type(mld_Tprec_type), input.
  !             The preconditioner data structure to be printed out.
  !  info    -  integer, output.
  !             error code.
  !  iout    -  integer, input, optional.
  !             The id of the file where the preconditioner description
  !             will be printed. If iout is not present, then the standard
  !             output is condidered.
  !  root    -  integer, input, optional.
  !             The id of the process printing the message; -1 acts as a wildcard.
  !             Default is psb_root_
  !
  subroutine mld_zfile_prec_descr(p,info,iout,root)
    implicit none 
    ! Arguments
    type(mld_zprec_type), intent(in)      :: p
    integer(psb_ipk_), intent(out)          :: info
    integer(psb_ipk_), intent(in), optional :: iout
    integer(psb_ipk_), intent(in), optional :: root

    ! Local variables
    integer(psb_ipk_)  :: ilev, nlev
    integer(psb_ipk_) :: ictxt, me, np
    character(len=20), parameter :: name='mld_file_prec_descr'
    integer(psb_ipk_)  :: iout_
    integer(psb_ipk_)  :: root_

    info = psb_success_
    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if
    if (iout_ < 0) iout_ = 6 

    ictxt = p%ictxt

    if (allocated(p%precv)) then

      call psb_info(ictxt,me,np)
      if (present(root)) then 
        root_ = root
      else
        root_ = psb_root_
      end if
      if (root_ == -1) root_ = me

      !
      ! The preconditioner description is printed by processor psb_root_.
      ! This agrees with the fact that all the parameters defining the
      ! preconditioner have the same values on all the procs (this is
      ! ensured by mld_precbld).
      !
      if (me == root_) then
        nlev = size(p%precv)
        do ilev = 1, nlev 
          if (.not.allocated(p%precv(ilev)%sm)) then 
            info = 3111
            write(iout_,*) ' ',name,&
                 & ': error: inconsistent MLPREC part, should call MLD_PRECINIT'
            return
          endif
        end do

        write(iout_,*) 
        write(iout_,'(a)') 'Preconditioner description'
        if (nlev >= 1) then
          !
          ! Print description of base preconditioner
          !
          if (nlev > 1) then
            write(iout_,*) 'Multilevel Schwarz'
            write(iout_,*) 
            write(iout_,*) 'Base preconditioner (smoother) details'
          endif
          call p%precv(1)%sm%descr(info,iout=iout_)
          if (nlev == 1) then 
            if (p%precv(1)%parms%sweeps > 1) then 
              write(iout_,*) '  Number of sweeps : ',&
                   & p%precv(1)%parms%sweeps 
            end if
            write(iout_,*) 
            return 
          end if
        end if

        !
        ! Print multilevel details
        !
        write(iout_,*) 
        write(iout_,*) 'Multilevel details'
        write(iout_,*) ' Number of levels   : ',nlev
        write(iout_,*) ' Operator complexity: ',p%get_complexity()
        do ilev=2,nlev
          call p%precv(ilev)%descr(ilev,nlev,info,iout=iout_)
        end do
        write(iout_,*) 
          
      end if

    else
      write(iout_,*) trim(name), &
           & ': Error: no base preconditioner available, something is wrong!'
      info = -2
      return
    endif
    

  end subroutine mld_zfile_prec_descr


  !
  ! Subroutines: mld_Tprec_free
  ! Version: complex
  !
  !  These routines deallocate the mld_Tprec_type data structures.
  !
  ! Arguments:
  !  p       -  type(mld_Tprec_type), input.
  !             The data structure to be deallocated.
  !  info    -  integer, output.
  !             error code.
  !
  subroutine mld_zprecfree(p,info)
  
    implicit none
    
    ! Arguments
    type(mld_zprec_type), intent(inout) :: p
    integer(psb_ipk_), intent(out)        :: info
    
    ! Local variables
    integer(psb_ipk_)   :: me,err_act,i
    character(len=20)   :: name
    
    if(psb_get_errstatus().ne.0) return 
    info=psb_success_
    name = 'mld_zprecfree'
    call psb_erractionsave(err_act)
    
    me=-1
    
    call p%free(info)

    
    return
    
  end subroutine mld_zprecfree

  subroutine mld_z_prec_free(prec,info)
  
    implicit none
    
    ! Arguments
    class(mld_zprec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out)        :: info
    
    ! Local variables
    integer(psb_ipk_)   :: me,err_act,i
    character(len=20)   :: name
    
    if(psb_get_errstatus().ne.0) return 
    info=psb_success_
    name = 'mld_zprecfree'
    call psb_erractionsave(err_act)
    
    me=-1
    
    if (allocated(prec%precv)) then 
      do i=1,size(prec%precv) 
        call prec%precv(i)%free(info)
      end do
      deallocate(prec%precv,stat=info)
    end if
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine mld_z_prec_free

  

  !
  ! Top level methods. 
  !
  subroutine mld_z_apply2_vect(prec,x,y,desc_data,info,trans,work)
    implicit none 
    type(psb_desc_type),intent(in)        :: desc_data
    class(mld_zprec_type), intent(inout)  :: prec
    type(psb_z_vect_type),intent(inout)   :: x
    type(psb_z_vect_type),intent(inout)   :: y
    integer(psb_ipk_), intent(out)          :: info
    character(len=1), optional              :: trans
    complex(psb_dpk_),intent(inout), optional, target :: work(:)
    Integer(psb_ipk_) :: err_act
    character(len=20) :: name='d_prec_apply'

    call psb_erractionsave(err_act)

    select type(prec) 
    type is (mld_zprec_type)
      call mld_precaply(prec,x,y,desc_data,info,trans,work)
    class default
      info = psb_err_missing_override_method_
      call psb_errpush(info,name)
      goto 9999 
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

  end subroutine mld_z_apply2_vect

  subroutine mld_z_apply1_vect(prec,x,desc_data,info,trans,work)
    implicit none 
    type(psb_desc_type),intent(in)          :: desc_data
    class(mld_zprec_type), intent(inout)  :: prec
    type(psb_z_vect_type),intent(inout)   :: x
    integer(psb_ipk_), intent(out)          :: info
    character(len=1), optional            :: trans
    complex(psb_dpk_),intent(inout), optional, target :: work(:)
    Integer(psb_ipk_) :: err_act
    character(len=20) :: name='d_prec_apply'

    call psb_erractionsave(err_act)

    select type(prec) 
    type is (mld_zprec_type)
      call mld_precaply(prec,x,desc_data,info,trans,work)
    class default
      info = psb_err_missing_override_method_
      call psb_errpush(info,name)
      goto 9999 
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

  end subroutine mld_z_apply1_vect


  subroutine mld_z_apply2v(prec,x,y,desc_data,info,trans,work)
    implicit none 
    type(psb_desc_type),intent(in)    :: desc_data
    class(mld_zprec_type), intent(inout) :: prec
    complex(psb_dpk_),intent(inout)      :: x(:)
    complex(psb_dpk_),intent(inout)      :: y(:)
    integer(psb_ipk_), intent(out)     :: info
    character(len=1), optional        :: trans
    complex(psb_dpk_),intent(inout), optional, target :: work(:)
    Integer(psb_ipk_) :: err_act
    character(len=20) :: name='d_prec_apply'

    call psb_erractionsave(err_act)

    select type(prec) 
    type is (mld_zprec_type)
      call mld_precaply(prec,x,y,desc_data,info,trans,work)
    class default
      info = psb_err_missing_override_method_
      call psb_errpush(info,name)
      goto 9999 
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

  end subroutine mld_z_apply2v

  subroutine mld_z_apply1v(prec,x,desc_data,info,trans)
    implicit none 
    type(psb_desc_type),intent(in)    :: desc_data
    class(mld_zprec_type), intent(inout) :: prec
    complex(psb_dpk_),intent(inout)      :: x(:)
    integer(psb_ipk_), intent(out)     :: info
    character(len=1), optional         :: trans
    Integer(psb_ipk_) :: err_act
    character(len=20) :: name='d_prec_apply'

    call psb_erractionsave(err_act)

    select type(prec) 
    type is (mld_zprec_type)
      call mld_precaply(prec,x,desc_data,info,trans)
    class default
      info = psb_err_missing_override_method_
      call psb_errpush(info,name)
      goto 9999 
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

  end subroutine mld_z_apply1v


  subroutine mld_z_dump(prec,info,istart,iend,prefix,head,ac,rp,smoother,solver)
    
    implicit none 
    class(mld_zprec_type), intent(in)     :: prec
    integer(psb_ipk_), intent(out)          :: info
    integer(psb_ipk_), intent(in), optional :: istart, iend
    character(len=*), intent(in), optional  :: prefix, head
    logical, optional, intent(in)    :: smoother, solver,ac, rp
    integer(psb_ipk_)  :: i, j, il1, iln, lname, lev
    integer(psb_ipk_)  :: icontxt,iam, np
    character(len=80)  :: prefix_
    character(len=120) :: fname ! len should be at least 20 more than
    !  len of prefix_ 

    info = 0

    iln = size(prec%precv)
    if (present(istart)) then 
      il1 = max(1,istart)
    else
      il1 = min(2,iln)
    end if
    if (present(iend)) then 
      iln = min(iln, iend)
    end if

    do lev=il1, iln
      call prec%precv(lev)%dump(lev,info,prefix=prefix,head=head,&
           & ac=ac,smoother=smoother,solver=solver,rp=rp)
    end do

  end subroutine mld_z_dump


  subroutine mld_z_clone(prec,precout,info)

    implicit none 
    class(mld_zprec_type), intent(inout) :: prec
    class(psb_zprec_type), intent(inout) :: precout
    integer(psb_ipk_), intent(out)       :: info
    
    call precout%free(info)
    if (info == 0) call mld_z_inner_clone(prec,precout,info) 

  end subroutine mld_z_clone

  subroutine mld_z_inner_clone(prec,precout,info)

    implicit none 
    class(mld_zprec_type), intent(inout)         :: prec
    class(psb_zprec_type), target, intent(inout) :: precout
    integer(psb_ipk_), intent(out)             :: info
    ! Local vars
    integer(psb_ipk_)  :: i, j, il1, ln, lname, lev
    integer(psb_ipk_)  :: icontxt,iam, np

    info = psb_success_
    select type(pout => precout)
    class is (mld_zprec_type)
      pout%ictxt            = prec%ictxt
      pout%coarse_aggr_size = prec%coarse_aggr_size
      pout%op_complexity    = prec%op_complexity
      if (allocated(prec%precv)) then 
        ln = size(prec%precv) 
        allocate(pout%precv(ln),stat=info)
        if (info /= psb_success_) goto 9999
        if (ln >= 1) then 
          call prec%precv(1)%clone(pout%precv(1),info)
        end if
        do lev=2, ln
          if (info /= psb_success_) exit
          call prec%precv(lev)%clone(pout%precv(lev),info)
          if (info == psb_success_) then 
            pout%precv(lev)%base_a       => pout%precv(lev)%ac
            pout%precv(lev)%base_desc    => pout%precv(lev)%desc_ac
            pout%precv(lev)%map%p_desc_X => pout%precv(lev-1)%base_desc
            pout%precv(lev)%map%p_desc_Y => pout%precv(lev)%base_desc
          end if
        end do
      end if
    class default 
      write(0,*) 'Error: wrong out type'
      info = psb_err_invalid_input_
    end select
9999 continue
  end subroutine mld_z_inner_clone

  subroutine mld_zprec_move_alloc(a, b,info)
    use psb_base_mod
    implicit none
    type(mld_zprec_type), intent(inout) :: a
    type(mld_zprec_type), intent(inout), target :: b
    integer(psb_ipk_), intent(out) :: info 
    integer(psb_ipk_) :: i
    
    if (allocated(b%precv)) then 
      ! This might not be required if FINAL procedures are available.
      call mld_precfree(b,info)
      if (info /= psb_success_) then 
        !       ?????
    !!$        return
      endif
    end if

    call move_alloc(a%precv,b%precv)
    ! Fix the pointers except on level 1.
    do i=2, size(b%precv)
      b%precv(i)%base_a    => b%precv(i)%ac
      b%precv(i)%base_desc => b%precv(i)%desc_ac
      b%precv(i)%map%p_desc_X => b%precv(i-1)%base_desc
      b%precv(i)%map%p_desc_Y => b%precv(i)%base_desc
    end do
  end subroutine mld_zprec_move_alloc
  
end module mld_z_prec_type
