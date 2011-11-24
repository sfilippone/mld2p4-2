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
! File: mld_prec_type.f90
!
! Module: mld_prec_type
!
!  This module defines: 
!  - the mld_prec_type data structure containing the preconditioner and related
!    data structures;
!
!  It contains routines for
!  - converting character constants defining the preconditioner into integer
!    constants; 
!  - checking if the preconditioner is correctly defined;
!  - printing a	description of the preconditioner;
!  - deallocating the preconditioner data structure.  
!

module mld_d_prec_type

  use mld_base_prec_type
  use psb_base_mod, only : psb_d_vect_type, psb_d_base_vect_type
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
  ! Type: mld_Tonelev_type.
  !
  !  It is the data type containing the necessary items for the	current
  !  level (essentially, the smoother, the current-level matrix
  !  and the restriction and prolongation operators).
  !
  !  type mld_Tonelev_type
  !    class(mld_T_base_smoother_type), allocatable :: sm
  !    type(mld_RTml_parms)            :: parms 
  !    type(psb_Tspmat_type)           :: ac
  !    type(psb_Tesc_type)             :: desc_ac
  !    type(psb_Tspmat_type), pointer  :: base_a    => null() 
  !    type(psb_Tesc_type), pointer    :: base_desc => null() 
  !    type(psb_Tlinmap_type)          :: map
  !  end type mld_Tonelev_type
  !
  !  Note that psb_Tpk denotes the kind of the real data type to be chosen
  !  according to single/double precision version of MLD2P4.
  !
  !   sm           -  class(mld_T_base_smoother_type), allocatable
  !                   The current level preconditioner (aka smoother).
  !   parms        -  type(mld_RTml_parms)
  !                   The parameters defining the multilevel strategy.
  !   ac           -  The local part of the current-level matrix, built by
  !                   coarsening the previous-level matrix.
  !   desc_ac      -  type(psb_desc_type).
  !                   The communication descriptor associated to the matrix
  !                   stored in ac.
  !   base_a       -  type(psb_Tspmat_type), pointer.
  !                   Pointer (really a pointer!) to the local part of the current 
  !                   matrix (so we have a unified treatment of residuals).
  !                   We need this to avoid passing explicitly the current matrix
  !                   to the routine which applies the preconditioner.
  !   base_desc    -  type(psb_desc_type), pointer.
  !                   Pointer to the communication descriptor associated to the
  !                   matrix pointed by base_a.
  !   map          -  Stores the maps (restriction and prolongation) between the
  !                   vector spaces associated to the index spaces of the previous
  !                   and current levels.
  !
  !   Methods:  
  !     Most methods follow the encapsulation hierarchy: they take whatever action
  !     is appropriate for the current object, then call the corresponding method for
  !     the contained object.
  !     As an example: the descr() method prints out a description of the
  !     level. It starts by invoking the descr() method of the parms object,
  !     then calls the descr() method of the smoother object. 
  !
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
  !                   it is passed to the smoother object for further processing.
  !    check      -   Sanity checks.
  !    sizeof     -   Total memory occupation in bytes
  !    get_nzeros -   Number of nonzeros 
  !
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

  type mld_d_base_solver_type
  contains
    procedure, pass(sv) :: check => d_base_solver_check
    procedure, pass(sv) :: dump  => d_base_solver_dmp
    procedure, pass(sv) :: build => d_base_solver_bld
    procedure, pass(sv) :: apply_v => d_base_solver_apply_vect
    procedure, pass(sv) :: apply_a => d_base_solver_apply
    generic, public     :: apply => apply_a, apply_v
    procedure, pass(sv) :: free  => d_base_solver_free
    procedure, pass(sv) :: seti  => d_base_solver_seti
    procedure, pass(sv) :: setc  => d_base_solver_setc
    procedure, pass(sv) :: setr  => d_base_solver_setr
    generic, public     :: set   => seti, setc, setr
    procedure, pass(sv) :: default => d_base_solver_default
    procedure, pass(sv) :: descr   => d_base_solver_descr
    procedure, pass(sv) :: sizeof  => d_base_solver_sizeof
    procedure, pass(sv) :: get_nzeros => d_base_solver_get_nzeros
  end type mld_d_base_solver_type

  type  mld_d_base_smoother_type
    class(mld_d_base_solver_type), allocatable :: sv
  contains
    procedure, pass(sm) :: check => d_base_smoother_check
    procedure, pass(sm) :: dump  => d_base_smoother_dmp
    procedure, pass(sm) :: build => d_base_smoother_bld
    procedure, pass(sm) :: apply_v => d_base_smoother_apply_vect
    procedure, pass(sm) :: apply_a => d_base_smoother_apply
    generic, public     :: apply => apply_a, apply_v
    procedure, pass(sm) :: free  => d_base_smoother_free
    procedure, pass(sm) :: seti  => d_base_smoother_seti
    procedure, pass(sm) :: setc  => d_base_smoother_setc
    procedure, pass(sm) :: setr  => d_base_smoother_setr
    generic, public     :: set   => seti, setc, setr
    procedure, pass(sm) :: default => d_base_smoother_default
    procedure, pass(sm) :: descr =>   d_base_smoother_descr
    procedure, pass(sm) :: sizeof =>  d_base_smoother_sizeof
    procedure, pass(sm) :: get_nzeros => d_base_smoother_get_nzeros
  end type mld_d_base_smoother_type

  type mld_donelev_type
    class(mld_d_base_smoother_type), allocatable :: sm
    type(mld_dml_parms)             :: parms 
    type(psb_dspmat_type)           :: ac
    type(psb_desc_type)             :: desc_ac
    type(psb_dspmat_type), pointer  :: base_a    => null() 
    type(psb_desc_type), pointer    :: base_desc => null() 
    type(psb_dlinmap_type)          :: map
  contains
    procedure, pass(lv) :: descr   => d_base_onelev_descr
    procedure, pass(lv) :: default => d_base_onelev_default
    procedure, pass(lv) :: free    => d_base_onelev_free
    procedure, pass(lv) :: nullify => d_base_onelev_nullify
    procedure, pass(lv) :: check => d_base_onelev_check
    procedure, pass(lv) :: dump  => d_base_onelev_dump
    procedure, pass(lv) :: seti  => d_base_onelev_seti
    procedure, pass(lv) :: setr  => d_base_onelev_setr
    procedure, pass(lv) :: setc  => d_base_onelev_setc
    generic, public     :: set   => seti, setr, setc
    procedure, pass(lv) :: sizeof => d_base_onelev_sizeof
    procedure, pass(lv) :: get_nzeros => d_base_onelev_get_nzeros
  end type mld_donelev_type

  type, extends(psb_dprec_type)         :: mld_dprec_type
    integer                             :: ictxt
    real(psb_dpk_)                      :: op_complexity=dzero
    type(mld_donelev_type), allocatable :: precv(:) 
  contains
    procedure, pass(prec)               :: d_apply2_vect => mld_d_apply2_vect
    procedure, pass(prec)               :: d_apply2v => mld_d_apply2v
    procedure, pass(prec)               :: d_apply1v => mld_d_apply1v
    procedure, pass(prec)               :: dump      => mld_d_dump
    procedure, pass(prec)               :: get_complexity => mld_d_get_compl
    procedure, pass(prec)               :: cmp_complexity => mld_d_cmp_compl
    procedure, pass(prec)               :: get_nzeros => mld_d_get_nzeros
  end type mld_dprec_type

  private :: d_base_solver_bld,  d_base_solver_apply, &
       &  d_base_solver_free,    d_base_solver_seti, &
       &  d_base_solver_setc,    d_base_solver_setr, &
       &  d_base_solver_descr,   d_base_solver_sizeof, &
       &  d_base_solver_default, d_base_solver_check,&
       &  d_base_solver_dmp, d_base_solver_apply_vect, &
       &  d_base_smoother_bld,   d_base_smoother_apply, &
       &  d_base_smoother_free,  d_base_smoother_seti, &
       &  d_base_smoother_setc,  d_base_smoother_setr,&
       &  d_base_smoother_descr, d_base_smoother_sizeof, &
       &  d_base_smoother_default, d_base_smoother_check, &
       &  d_base_smoother_dmp, d_base_smoother_apply_vect, &
       &  d_base_onelev_seti, d_base_onelev_setc, &
       &  d_base_onelev_setr, d_base_onelev_check, &
       &  d_base_onelev_default, d_base_onelev_dump, &
       &  d_base_onelev_descr,  d_base_onelev_sizeof, &
       &  d_base_onelev_free,  d_base_onelev_nullify,&
       &  mld_d_dump, &
       &  mld_d_get_compl,  mld_d_cmp_compl,&
       &  mld_d_get_nzeros, d_base_onelev_get_nzeros, &
       &  d_base_smoother_get_nzeros, d_base_solver_get_nzeros


  !
  ! Interfaces to routines for checking the definition of the preconditioner,
  ! for printing its description and for deallocating its data structure
  !

  interface mld_precfree
    module procedure mld_dprec_free
  end interface

  interface mld_nullify_onelevprec
    module procedure  mld_nullify_d_onelevprec
  end interface

  interface mld_precdescr
    module procedure mld_dfile_prec_descr
  end interface

  interface mld_sizeof
    module procedure mld_dprec_sizeof
  end interface

  interface mld_precaply
    subroutine mld_dprecaply_vect(prec,x,y,desc_data,info,trans,work)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, &
           & psb_dpk_, psb_d_vect_type
      import mld_dprec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_dprec_type), intent(inout) :: prec
      type(psb_d_vect_type),intent(inout) :: x
      type(psb_d_vect_type),intent(inout) :: y
      integer, intent(out)                :: info
      character(len=1), optional          :: trans
      real(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine mld_dprecaply_vect
    subroutine mld_dprecaply(prec,x,y,desc_data,info,trans,work)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      import mld_dprec_type
      type(psb_desc_type),intent(in)   :: desc_data
      type(mld_dprec_type), intent(in) :: prec
      real(psb_dpk_),intent(in)        :: x(:)
      real(psb_dpk_),intent(inout)     :: y(:)
      integer, intent(out)             :: info
      character(len=1), optional       :: trans
      real(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine mld_dprecaply
    subroutine mld_dprecaply1(prec,x,desc_data,info,trans)
      use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_dpk_
      import mld_dprec_type
      type(psb_desc_type),intent(in)   :: desc_data
      type(mld_dprec_type), intent(in) :: prec
      real(psb_dpk_),intent(inout)     :: x(:)
      integer, intent(out)             :: info
      character(len=1), optional       :: trans
    end subroutine mld_dprecaply1
  end interface

contains
  !
  ! Function returning the size of the mld_prec_type data structure
  ! in bytes or in number of nonzeros of the operator(s) involved. 
  !

  function d_base_solver_get_nzeros(sv) result(val)
    implicit none 
    class(mld_d_base_solver_type), intent(in) :: sv
    integer(psb_long_int_k_) :: val
    integer             :: i
    val = 0
  end function d_base_solver_get_nzeros

  function d_base_smoother_get_nzeros(sm) result(val)
    implicit none 
    class(mld_d_base_smoother_type), intent(in) :: sm
    integer(psb_long_int_k_) :: val
    integer             :: i
    val = 0
    if (allocated(sm%sv)) &
         &  val =  sm%sv%get_nzeros()
  end function d_base_smoother_get_nzeros

  function d_base_onelev_get_nzeros(lv) result(val)
    implicit none 
    class(mld_donelev_type), intent(in) :: lv
    integer(psb_long_int_k_) :: val
    integer             :: i
    val = 0
    if (allocated(lv%sm)) &
         &  val =  lv%sm%get_nzeros()
  end function d_base_onelev_get_nzeros

  function mld_d_get_nzeros(prec) result(val)
    implicit none 
    class(mld_dprec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    val = 0
    if (allocated(prec%precv)) then 
      do i=1, size(prec%precv)
        val = val + prec%precv(i)%get_nzeros()
      end do
    end if
  end function mld_d_get_nzeros


  function mld_dprec_sizeof(prec) result(val)
    implicit none 
    type(mld_dprec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    val = 0
    val = val + psb_sizeof_int
    if (allocated(prec%precv)) then 
      do i=1, size(prec%precv)
        val = val + prec%precv(i)%sizeof()
      end do
    end if
  end function mld_dprec_sizeof

  function d_base_onelev_sizeof(lv) result(val)
    implicit none 
    class(mld_donelev_type), intent(in) :: lv
    integer(psb_long_int_k_) :: val
    integer             :: i
    
    val = 0
    val = val + lv%desc_ac%sizeof()
    val = val + lv%ac%sizeof()
    val = val + lv%map%sizeof() 
    if (allocated(lv%sm))  val = val + lv%sm%sizeof()
  end function d_base_onelev_sizeof


  !
  ! Operator complexity: ratio of total number
  ! of nonzeros in the aggregated matrices at the
  ! various level to the nonzeroes at the fine level
  ! (original matrix)
  !
  
  function mld_d_get_compl(prec) result(val)
    implicit none 
    class(mld_dprec_type), intent(in) :: prec
    real(psb_dpk_)  :: val
    
    val = prec%op_complexity

  end function mld_d_get_compl
  
  subroutine mld_d_cmp_compl(prec) 
    use psb_base_mod, only : psb_min, psb_sum
    implicit none 
    class(mld_dprec_type), intent(inout) :: prec
    
    real(psb_dpk_) :: num,den
    integer  :: ictxt, il 

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
  end subroutine mld_d_cmp_compl
  
  !
  ! Subroutine: mld_file_prec_descr
  ! Version: real
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
  !
  subroutine mld_dfile_prec_descr(p,info,iout)
    implicit none 
    ! Arguments
    type(mld_dprec_type), intent(in) :: p
    integer, intent(out)             :: info
    integer, intent(in), optional    :: iout

    ! Local variables
    integer      :: ilev, nlev
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_file_prec_descr'
    integer :: iout_

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

      !
      ! The preconditioner description is printed by processor psb_root_.
      ! This agrees with the fact that all the parameters defining the
      ! preconditioner have the same values on all the procs (this is
      ! ensured by mld_precbld).
      !
      if (me == psb_root_) then
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
    

  end subroutine mld_dfile_prec_descr

  subroutine d_base_onelev_descr(lv,il,nl,info,iout)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_donelev_type), intent(in) :: lv
    integer, intent(in)                 :: il,nl
    integer, intent(out)                :: info
    integer, intent(in), optional       :: iout

    ! Local variables
    integer      :: err_act
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_d_base_onelev_descr'
    integer      :: iout_
    logical      :: coarse


    call psb_erractionsave(err_act)


    coarse = (il==nl)

    if (present(iout)) then 
      iout_ = iout
    else 
      iout_ = 6
    end if
    
    write(iout_,*) 
    if (il == 2) then 
      call lv%parms%mldescr(iout_,info)
      write(iout_,*) 
    end if

    if (coarse)  then 
      write(iout_,*) ' Level ',il,' (coarsest)'
    else
      write(iout_,*) ' Level ',il
    end if

    call lv%parms%descr(iout_,info,coarse=coarse)
        
    if (nl > 1) then 
      if (allocated(lv%map%naggr)) then
        write(iout_,*) '  Size of coarse matrix: ', &
             &  sum(lv%map%naggr(:))
        write(iout_,*) '  Sizes of aggregates: ', &
             &  lv%map%naggr(:)
      end if
    end if

    if (coarse.and.allocated(lv%sm)) &
         & call lv%sm%descr(info,iout=iout_,coarse=coarse)

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_onelev_descr


  !
  ! Subroutines: mld_Tbase_precfree, mld_T_onelev_precfree, mld_Tprec_free
  ! Version: real/complex
  !
  !  These routines deallocate the mld_Tbaseprec_type, mld_Tonelev_type and
  !  mld_Tprec_type data structures.
  !
  ! Arguments:
  !  p       -  type(mld_Tbaseprec_type/mld_Tonelev_type/mld_Tprec_type), input.
  !             The data structure to be deallocated.
  !  info    -  integer, output.
  !             error code.
  !
  subroutine d_base_onelev_free(lv,info)
    use psb_base_mod
    implicit none 

    class(mld_donelev_type), intent(inout) :: lv
    integer, intent(out)                :: info
    integer :: i

    info = psb_success_

    ! Actually we might just deallocate the top level array, except 
    ! for the inner UMFPACK or SLU stuff.
    ! We really need FINALs. 
    call lv%sm%free(info)
    
    call lv%ac%free()
    if (psb_is_ok_desc(lv%desc_ac)) &
         & call psb_cdfree(lv%desc_ac,info)

    ! This is a pointer to something else, must not free it here. 
    nullify(lv%base_a) 
    ! This is a pointer to something else, must not free it here. 
    nullify(lv%base_desc) 

    !
    ! free explicitly map???
    ! For now thanks to allocatable semantics
    ! works anyway. 
    !
    
    call lv%nullify()

  end subroutine d_base_onelev_free


  subroutine d_base_onelev_nullify(lv)
    implicit none 

    class(mld_donelev_type), intent(inout) :: lv

    nullify(lv%base_a) 
    nullify(lv%base_desc) 

  end subroutine d_base_onelev_nullify


  subroutine mld_nullify_d_onelevprec(p)
    implicit none 

    type(mld_donelev_type), intent(inout) :: p

    nullify(p%base_a) 
    nullify(p%base_desc) 

  end subroutine mld_nullify_d_onelevprec

  subroutine mld_dprec_free(p,info)
  
    use psb_base_mod
    
    implicit none
    
    ! Arguments
    type(mld_dprec_type), intent(inout) :: p
    integer, intent(out)                :: info
    
    ! Local variables
    integer             :: me,err_act,i
    character(len=20)   :: name
    
    if(psb_get_errstatus().ne.0) return 
    info=psb_success_
    name = 'mld_dprecfree'
    call psb_erractionsave(err_act)
    
    me=-1
    
    if (allocated(p%precv)) then 
      do i=1,size(p%precv) 
        call p%precv(i)%free(info)
      end do
      deallocate(p%precv)
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
    
  end subroutine mld_dprec_free

  
  !
  ! Smoother related routines/methods. 
  !


  !
  ! Apply: comes in two versions, on plain arrays or on encapsulated
  ! vectors.
  ! This basic version just applies the local solver, whatever that
  ! is. 
  !

  subroutine d_base_smoother_apply(alpha,sm,x,beta,y,desc_data,trans,sweeps,work,info)
    use psb_base_mod
    type(psb_desc_type), intent(in)             :: desc_data
    class(mld_d_base_smoother_type), intent(in) :: sm
    real(psb_dpk_),intent(inout)                :: x(:)
    real(psb_dpk_),intent(inout)                :: y(:)
    real(psb_dpk_),intent(in)                   :: alpha,beta
    character(len=1),intent(in)                 :: trans
    integer, intent(in)                         :: sweeps
    real(psb_dpk_),target, intent(inout)        :: work(:)
    integer, intent(out)                        :: info
    
    Integer           :: err_act
    character(len=20) :: name='d_base_smoother_apply'

    call psb_erractionsave(err_act)
    info = psb_success_
    if (allocated(sm%sv)) then 
      call sm%sv%apply(alpha,x,beta,y,desc_data,trans,work,info)
    else
      info = 1121
    endif
    if (info /= psb_success_) then 
      call psb_errpush(info,name)
      goto 9999 
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
    
  end subroutine d_base_smoother_apply

  subroutine d_base_smoother_apply_vect(alpha,sm,x,beta,y,desc_data,&
       &  trans,sweeps,work,info)
    use psb_base_mod
    type(psb_desc_type), intent(in)                :: desc_data
    class(mld_d_base_smoother_type), intent(inout) :: sm
    type(psb_d_vect_type),intent(inout)            :: x
    type(psb_d_vect_type),intent(inout)            :: y
    real(psb_dpk_),intent(in)                      :: alpha,beta
    character(len=1),intent(in)                    :: trans
    integer, intent(in)                            :: sweeps
    real(psb_dpk_),target, intent(inout)           :: work(:)
    integer, intent(out)                           :: info
    
    Integer           :: err_act
    character(len=20) :: name='d_base_smoother_apply'

    call psb_erractionsave(err_act)
    info = psb_success_
    if (allocated(sm%sv)) then 
      call sm%sv%apply(alpha,x,beta,y,desc_data,trans,work,info)
    else
      info = 1121
    endif
    if (info /= psb_success_) then 
      call psb_errpush(info,name)
      goto 9999 
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
    
  end subroutine d_base_smoother_apply_vect

  !
  ! Check:
  ! 1. Check that we do have a solver object
  ! 2. Call its check method
  !

  subroutine d_base_smoother_check(sm,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_base_smoother_type), intent(inout) :: sm 
    integer, intent(out)                   :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_smoother_check'

    call psb_erractionsave(err_act)
    info = psb_success_

    if (allocated(sm%sv)) then 
      call sm%sv%check(info)
    else 
      info=3111
      call psb_errpush(info,name)
      goto 9999
    end if

    if (info /= psb_success_) goto 9999
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_smoother_check

  !
  ! Set methods: the come in multiple versions according
  ! to whether we are setting with integer, real or character
  ! input.
  ! The basic rule is: if the input refers to a parameter
  ! of the smoother, use it, otherwise pass it to the
  ! solver object for further processing.
  ! Since there are no parameters in the base smoother
  ! we just pass everything to the solver object. 
  !

  subroutine d_base_smoother_seti(sm,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_base_smoother_type), intent(inout) :: sm 
    integer, intent(in)                            :: what 
    integer, intent(in)                            :: val
    integer, intent(out)                           :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_smoother_seti'

    call psb_erractionsave(err_act)
    info = psb_success_

    if (allocated(sm%sv)) then 
      call sm%sv%set(what,val,info)
    end if
    if (info /= psb_success_) goto 9999
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_smoother_seti

  subroutine d_base_smoother_setc(sm,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_base_smoother_type), intent(inout) :: sm 
    integer, intent(in)                            :: what 
    character(len=*), intent(in)                   :: val
    integer, intent(out)                           :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_smoother_setc'

    call psb_erractionsave(err_act)

    info = psb_success_

    if (allocated(sm%sv)) then 
      call sm%sv%set(what,val,info)
    end if
    if (info /= psb_success_) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_smoother_setc
  
  subroutine d_base_smoother_setr(sm,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_base_smoother_type), intent(inout) :: sm 
    integer, intent(in)                            :: what 
    real(psb_dpk_), intent(in)                     :: val
    integer, intent(out)                           :: info
    Integer :: err_act
    character(len=20)  :: name='d_base_smoother_setr'

    call psb_erractionsave(err_act)


    info = psb_success_

    if (allocated(sm%sv)) then 
      call sm%sv%set(what,val,info)
    end if
    if (info /= psb_success_) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_smoother_setr



  !
  ! Build method.
  ! At base level we only have to pass data to the inner solver. 
  ! AMOLD/VMOLD allow to have any relevant sparse matrix or vector
  ! to be stored in a given format. This is essential e.g.
  ! when dealing  with GPUs. 
  !
  subroutine d_base_smoother_bld(a,desc_a,sm,upd,info,amold,vmold)

    use psb_base_mod

    Implicit None

    ! Arguments
    type(psb_dspmat_type), intent(in), target      :: a
    Type(psb_desc_type), Intent(in)                :: desc_a 
    class(mld_d_base_smoother_type), intent(inout) :: sm 
    character, intent(in)                          :: upd
    integer, intent(out)                           :: info
    class(psb_d_base_sparse_mat), intent(in), optional :: amold
    class(psb_d_base_vect_type), intent(in), optional  :: vmold
    Integer           :: err_act
    character(len=20) :: name='d_base_smoother_bld'

    call psb_erractionsave(err_act)

    info = psb_success_
    if (allocated(sm%sv)) then 
      call sm%sv%build(a,desc_a,upd,info,amold=amold,vmold=vmold)
    else
      info = 1121
      call psb_errpush(info,name)
    endif
    if (info /= psb_success_) goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_smoother_bld

  !
  ! Free method (aka destructor).
  ! For this one actually we could do without; however
  ! for cases where there are data objects allocated outside
  ! of the Fortran RTE we need to free them explicitly.
  !
  ! Even in that case, we could do without this if FINAL
  ! subroutines were supported, which is not the case
  ! in GNU Fortran up to 4.7. 
  !
  subroutine d_base_smoother_free(sm,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_base_smoother_type), intent(inout) :: sm
    integer, intent(out)                           :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_smoother_free'

    call psb_erractionsave(err_act)
    info = psb_success_
    
    if (allocated(sm%sv)) then 
      call sm%sv%free(info)
    end if
    if (info == psb_success_) deallocate(sm%sv,stat=info) 
    if (info /= psb_success_) then 
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
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
  end subroutine d_base_smoother_free

  !
  ! Print a description
  !

  subroutine d_base_smoother_descr(sm,info,iout,coarse)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_base_smoother_type), intent(in) :: sm
    integer, intent(out)                        :: info
    integer, intent(in), optional               :: iout
    logical, intent(in), optional               :: coarse

    ! Local variables
    integer      :: err_act
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_d_base_smoother_descr'
    integer :: iout_
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
      iout_ = 6
    end if

    if (.not.coarse_) &
         &  write(iout_,*) 'Base smoother with local solver'
    if (allocated(sm%sv)) then 
      call sm%sv%descr(info,iout,coarse)
      if (info /= psb_success_) then 
        info = psb_err_from_subroutine_ 
        call psb_errpush(info,name,a_err='Local solver')
        goto 9999
      end if
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
  end subroutine d_base_smoother_descr

  !
  ! Dump 
  ! to file, for debugging purposes.
  !
  subroutine d_base_smoother_dmp(sm,ictxt,level,info,prefix,head,smoother,solver)
    use psb_base_mod
    implicit none 
    class(mld_d_base_smoother_type), intent(in) :: sm
    integer, intent(in)              :: ictxt,level
    integer, intent(out)             :: info
    character(len=*), intent(in), optional :: prefix, head
    logical, optional, intent(in)    :: smoother, solver
    integer :: i, j, il1, iln, lname, lev
    integer :: icontxt,iam, np
    character(len=80)  :: prefix_
    character(len=120) :: fname ! len should be at least 20 more than
    logical :: smoother_
    !  len of prefix_ 

    info = 0

    if (present(prefix)) then 
      prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
    else
      prefix_ = "dump_smth_d"
    end if

    call psb_info(ictxt,iam,np)

    if (present(smoother)) then 
      smoother_ = smoother
    else
      smoother_ = .false. 
    end if
    lname = len_trim(prefix_)
    fname = trim(prefix_)
    write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
    lname = lname + 5

    ! At base level do nothing for the smoother
    if (allocated(sm%sv)) &
         & call sm%sv%dump(ictxt,level,info,solver=solver)

  end subroutine d_base_smoother_dmp

  function d_base_smoother_sizeof(sm) result(val)
    implicit none 
    ! Arguments
    class(mld_d_base_smoother_type), intent(in) :: sm
    integer(psb_long_int_k_)                    :: val
    integer             :: i
    
    val = 0
    if (allocated(sm%sv)) then 
      val = sm%sv%sizeof()
    end if

    return
  end function d_base_smoother_sizeof


  !
  ! Set sensible defaults.
  ! To be called immediately after allocation
  !
  subroutine d_base_smoother_default(sm) 
    implicit none 
    ! Arguments
    class(mld_d_base_smoother_type), intent(inout) :: sm
    ! Do nothing for base version

    if (allocated(sm%sv)) call sm%sv%default()

    return
  end subroutine d_base_smoother_default

  
  !
  ! Local solver  related routines/methods. 
  !


  !
  ! Apply: comes in two versions, on plain arrays or on encapsulated
  ! vectors.
  ! The base version throws an error, since it means that no explicit
  ! choice was made. 
  ! Question: would it make sense to transform the base version into
  ! the ID version, i.e. "solver" is the identity operator? 
  ! 


  subroutine d_base_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info)
    use psb_base_mod
    type(psb_desc_type), intent(in)           :: desc_data
    class(mld_d_base_solver_type), intent(in) :: sv
    real(psb_dpk_),intent(inout)              :: x(:)
    real(psb_dpk_),intent(inout)              :: y(:)
    real(psb_dpk_),intent(in)                 :: alpha,beta
    character(len=1),intent(in)               :: trans
    real(psb_dpk_),target, intent(inout)      :: work(:)
    integer, intent(out)                      :: info
    
    Integer :: err_act
    character(len=20)  :: name='d_base_solver_apply'

    call psb_erractionsave(err_act)
    
    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine d_base_solver_apply

  subroutine d_base_solver_apply_vect(alpha,sv,x,beta,y,desc_data,trans,work,info)
    use psb_base_mod
    type(psb_desc_type), intent(in)              :: desc_data
    class(mld_d_base_solver_type), intent(inout) :: sv
    type(psb_d_vect_type),intent(inout)          :: x
    type(psb_d_vect_type),intent(inout)          :: y
    real(psb_dpk_),intent(in)                    :: alpha,beta
    character(len=1),intent(in)                  :: trans
    real(psb_dpk_),target, intent(inout)         :: work(:)
    integer, intent(out)                         :: info
    
    Integer :: err_act
    character(len=20)  :: name='d_base_solver_apply'

    call psb_erractionsave(err_act)
    
    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine d_base_solver_apply_vect


  !
  ! Build
  ! The base version throws an error, since it means that no explicit
  ! choice was made. 
  !
  subroutine d_base_solver_bld(a,desc_a,sv,upd,info,b,amold,vmold)

    use psb_base_mod

    Implicit None

    ! Arguments
    type(psb_dspmat_type), intent(in), target           :: a
    Type(psb_desc_type), Intent(in)                     :: desc_a 
    class(mld_d_base_solver_type), intent(inout)        :: sv
    character, intent(in)                               :: upd
    integer, intent(out)                                :: info
    type(psb_dspmat_type), intent(in), target, optional :: b
    class(psb_d_base_sparse_mat), intent(in), optional  :: amold
    class(psb_d_base_vect_type), intent(in), optional   :: vmold

    Integer :: err_act
    character(len=20)  :: name='d_base_solver_bld'

    call psb_erractionsave(err_act)

    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_solver_bld

  subroutine d_base_solver_check(sv,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_base_solver_type), intent(inout) :: sv
    integer, intent(out)                   :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_solver_check'

    call psb_erractionsave(err_act)
    info = psb_success_


    if (info /= psb_success_) goto 9999
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_solver_check

  subroutine d_base_solver_seti(sv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_base_solver_type), intent(inout) :: sv 
    integer, intent(in)                          :: what 
    integer, intent(in)                          :: val
    integer, intent(out)                         :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_solver_seti'
    
    ! Correct action here is doing nothing. 
    info = 0
    
    return
  end subroutine d_base_solver_seti

  subroutine d_base_solver_setc(sv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_base_solver_type), intent(inout) :: sv
    integer, intent(in)                          :: what 
    character(len=*), intent(in)                 :: val
    integer, intent(out)                         :: info
    Integer           :: err_act, ival 
    character(len=20) :: name='d_base_solver_setc'

    call psb_erractionsave(err_act)

    info = psb_success_

    call mld_stringval(val,ival,info)
    if (info == psb_success_) call sv%set(what,ival,info)

    if (info /= psb_success_) goto 9999


    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_solver_setc
  
  subroutine d_base_solver_setr(sv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_base_solver_type), intent(inout) :: sv 
    integer, intent(in)                          :: what 
    real(psb_dpk_), intent(in)                   :: val
    integer, intent(out)                         :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_solver_setr'

    
    ! Correct action here is doing nothing. 
    info = 0
    
    return
  end subroutine d_base_solver_setr

  !
  ! Free
  ! The base version throws an error, since it means that no explicit
  ! choice was made. IS THIS CORRECT? I suspect it would be better
  ! to be silent here, to cover reallocation. 
  !

  subroutine d_base_solver_free(sv,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_base_solver_type), intent(inout) :: sv
    integer, intent(out)                         :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_solver_free'

    call psb_erractionsave(err_act)

    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_solver_free

  subroutine d_base_solver_descr(sv,info,iout,coarse)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_base_solver_type), intent(in) :: sv
    integer, intent(out)                      :: info
    integer, intent(in), optional             :: iout
    logical, intent(in), optional             :: coarse

    ! Local variables
    integer      :: err_act
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_d_base_solver_descr'
    integer      :: iout_


    call psb_erractionsave(err_act)

    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_solver_descr

  subroutine d_base_solver_dmp(sv,ictxt,level,info,prefix,head,solver)
    use psb_base_mod
    implicit none 
    class(mld_d_base_solver_type), intent(in) :: sv
    integer, intent(in)              :: ictxt,level
    integer, intent(out)             :: info
    character(len=*), intent(in), optional :: prefix, head
    logical, optional, intent(in)    :: solver
    integer :: i, j, il1, iln, lname, lev
    integer :: icontxt,iam, np
    character(len=80)  :: prefix_
    character(len=120) :: fname ! len should be at least 20 more than
    logical :: solver_
    !  len of prefix_ 

    info = 0

    if (present(prefix)) then 
      prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
    else
      prefix_ = "dump_slv_d"
    end if

    call psb_info(ictxt,iam,np)

    if (present(solver)) then 
      solver_ = solver
    else
      solver_ = .false. 
    end if
    lname = len_trim(prefix_)
    fname = trim(prefix_)
    write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
    lname = lname + 5

    ! At base level do nothing for the solver

  end subroutine d_base_solver_dmp

  function d_base_solver_sizeof(sv) result(val)
    implicit none 
    ! Arguments
    class(mld_d_base_solver_type), intent(in) :: sv
    integer(psb_long_int_k_)                  :: val
    integer             :: i
    val = 0

    return
  end function d_base_solver_sizeof

  subroutine d_base_solver_default(sv) 
    implicit none 
    ! Arguments
    class(mld_d_base_solver_type), intent(inout) :: sv
    ! Do nothing for base version

    return
  end subroutine d_base_solver_default

  !
  ! Onelevel checks. 
  ! The number of Jacobi sweeps to be applied is not
  ! tied to the Jacobi smoother: logically, you have
  ! a smoother and you can choose to apply it any number
  ! of times you like. 
  !
  subroutine d_base_onelev_check(lv,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_donelev_type), intent(inout) :: lv 
    integer, intent(out)                   :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_onelev_check'

    call psb_erractionsave(err_act)
    info = psb_success_

    call mld_check_def(lv%parms%sweeps,&
         & 'Jacobi sweeps',1,is_legal_jac_sweeps)
    call mld_check_def(lv%parms%sweeps_pre,&
         & 'Jacobi sweeps',1,is_legal_jac_sweeps)
    call mld_check_def(lv%parms%sweeps_post,&
         & 'Jacobi sweeps',1,is_legal_jac_sweeps)

    
    if (allocated(lv%sm)) then 
      call lv%sm%check(info)
    else 
      info=3111
      call psb_errpush(info,name)
      goto 9999
    end if

    if (info /= psb_success_) goto 9999
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_onelev_check

  !
  ! Multilevel defaults: 
  !  multiplicative vs. additive ML framework;
  !  Smoothed decoupled aggregation with zero threshold; 
  !  distributed coarse matrix;
  !  damping omega  computed with the max-norm estimate of the
  !  dominant eigenvalue;
  !  two-sided smoothing (i.e. V-cycle) with 1 smoothing sweep;
  !

  subroutine d_base_onelev_default(lv)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_donelev_type), intent(inout) :: lv 

    lv%parms%sweeps          = 1
    lv%parms%sweeps_pre      = 1
    lv%parms%sweeps_post     = 1
    lv%parms%ml_type         = mld_mult_ml_
    lv%parms%aggr_alg        = mld_dec_aggr_
    lv%parms%aggr_kind       = mld_smooth_prol_
    lv%parms%coarse_mat      = mld_distr_mat_
    lv%parms%smoother_pos    = mld_twoside_smooth_
    lv%parms%aggr_omega_alg  = mld_eig_est_
    lv%parms%aggr_eig        = mld_max_norm_
    lv%parms%aggr_filter     = mld_no_filter_mat_
    lv%parms%aggr_omega_val  = dzero
    lv%parms%aggr_thresh     = dzero
    
    if (allocated(lv%sm)) call lv%sm%default()

    return

  end subroutine d_base_onelev_default

  ! 
  ! Set routines:
  ! Parameters belonging here are:
  !   Number of smoothing sweeps;
  !   Smoother position;
  !   Aggregation related parameters
  !   Record request on coarse level solver,
  !    for checks on solver vs. smoother nomenclature
  !    reconciliation. 
  ! 
  subroutine d_base_onelev_seti(lv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_donelev_type), intent(inout) :: lv 
    integer, intent(in)                          :: what 
    integer, intent(in)                          :: val
    integer, intent(out)                         :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_onelev_seti'

    call psb_erractionsave(err_act)
    info = psb_success_

    select case (what) 

    case (mld_smoother_sweeps_)
      lv%parms%sweeps      = val
      lv%parms%sweeps_pre  = val
      lv%parms%sweeps_post = val

    case (mld_smoother_sweeps_pre_)
      lv%parms%sweeps_pre  = val

    case (mld_smoother_sweeps_post_)
      lv%parms%sweeps_post = val

    case (mld_ml_type_)
      lv%parms%ml_type       = val

    case (mld_aggr_alg_)
      lv%parms%aggr_alg      = val

    case (mld_aggr_kind_)
      lv%parms%aggr_kind     = val

    case (mld_coarse_mat_)
      lv%parms%coarse_mat    = val

    case (mld_smoother_pos_)
      lv%parms%smoother_pos  = val

    case (mld_aggr_omega_alg_)
      lv%parms%aggr_omega_alg= val

    case (mld_aggr_eig_)
      lv%parms%aggr_eig      = val

    case (mld_aggr_filter_)
      lv%parms%aggr_filter   = val

    case (mld_coarse_solve_)
      lv%parms%coarse_solve    = val

    case default
      if (allocated(lv%sm)) then 
        call lv%sm%set(what,val,info)
      end if
      if (info /= psb_success_) goto 9999
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
  end subroutine d_base_onelev_seti

  subroutine d_base_onelev_setc(lv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_donelev_type), intent(inout) :: lv 
    integer, intent(in)                            :: what 
    character(len=*), intent(in)                   :: val
    integer, intent(out)                           :: info
    Integer           :: err_act
    character(len=20) :: name='d_base_onelev_setc'
    integer :: ival 

    call psb_erractionsave(err_act)

    info = psb_success_

    call mld_stringval(val,ival,info)
    if (info == psb_success_) call lv%set(what,ival,info)

    if (info /= psb_success_) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_onelev_setc
  
  subroutine d_base_onelev_setr(lv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_donelev_type), intent(inout) :: lv 
    integer, intent(in)                            :: what 
    real(psb_dpk_), intent(in)                     :: val
    integer, intent(out)                           :: info
    Integer :: err_act
    character(len=20)  :: name='d_base_onelev_setr'

    call psb_erractionsave(err_act)


    info = psb_success_
    
    select case (what) 

    case (mld_aggr_omega_val_)
      lv%parms%aggr_omega_val= val
   
    case (mld_aggr_thresh_)
      lv%parms%aggr_thresh   = val

    case default
      if (allocated(lv%sm)) then 
        call lv%sm%set(what,val,info)
      end if
      if (info /= psb_success_) goto 9999
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
  end subroutine d_base_onelev_setr

  !
  ! Dump on file: can be fine-tuned to include the (aggregated) matrix
  ! as well as smoother and solver. 
  !
  subroutine d_base_onelev_dump(lv,level,info,prefix,head,ac,rp,smoother,solver)
    use psb_base_mod
    implicit none 
    class(mld_donelev_type), intent(in) :: lv
    integer, intent(in)              :: level
    integer, intent(out)             :: info
    character(len=*), intent(in), optional :: prefix, head
    logical, optional, intent(in)    :: ac, rp, smoother, solver
    integer :: i, j, il1, iln, lname, lev
    integer :: icontxt,iam, np
    character(len=80)  :: prefix_
    character(len=120) :: fname ! len should be at least 20 more than
    logical :: ac_, rp_
    !  len of prefix_ 

    info = 0

    if (present(prefix)) then 
      prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
    else
      prefix_ = "dump_lev_d"
    end if

    if (associated(lv%base_desc)) then 
      icontxt = lv%base_desc%get_context()
      call psb_info(icontxt,iam,np)
    else 
      icontxt = -1 
      iam     = -1
    end if
    if (present(ac)) then 
      ac_ = ac
    else
      ac_ = .false. 
    end if
    if (present(rp)) then 
      rp_ = rp
    else
      rp_ = .false. 
    end if
    lname = len_trim(prefix_)
    fname = trim(prefix_)
    write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
    lname = lname + 5

    if (level >= 2) then 
      if (ac_) then 
        write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_ac.mtx'
        write(0,*) 'Filename ',fname
        call lv%ac%print(fname,head=head)
      end if
      if (rp_) then 
        write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_r.mtx'
        write(0,*) 'Filename ',fname
        call lv%map%map_X2Y%print(fname,head=head)
        write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_p.mtx'
        write(0,*) 'Filename ',fname
        call lv%map%map_Y2X%print(fname,head=head)
      end if
    end if
    if (allocated(lv%sm)) &
         & call lv%sm%dump(icontxt,level,info,smoother=smoother,solver=solver)

  end subroutine d_base_onelev_dump


  !
  ! Top level methods. 
  !
  subroutine mld_d_apply2_vect(prec,x,y,desc_data,info,trans,work)
    use psb_base_mod
    type(psb_desc_type),intent(in)        :: desc_data
    class(mld_dprec_type), intent(inout)  :: prec
    type(psb_d_vect_type),intent(inout)   :: x
    type(psb_d_vect_type),intent(inout)   :: y
    integer, intent(out)                  :: info
    character(len=1), optional            :: trans
    real(psb_dpk_),intent(inout), optional, target :: work(:)
    Integer           :: err_act
    character(len=20) :: name='d_prec_apply'

    call psb_erractionsave(err_act)

    select type(prec) 
    type is (mld_dprec_type)
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

  end subroutine mld_d_apply2_vect


  subroutine mld_d_apply2v(prec,x,y,desc_data,info,trans,work)
    use psb_base_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(mld_dprec_type), intent(in) :: prec
    real(psb_dpk_),intent(inout)      :: x(:)
    real(psb_dpk_),intent(inout)      :: y(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans
    real(psb_dpk_),intent(inout), optional, target :: work(:)
    Integer           :: err_act
    character(len=20) :: name='d_prec_apply'

    call psb_erractionsave(err_act)

    select type(prec) 
    type is (mld_dprec_type)
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

  end subroutine mld_d_apply2v

  subroutine mld_d_apply1v(prec,x,desc_data,info,trans)
    use psb_base_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(mld_dprec_type), intent(in) :: prec
    real(psb_dpk_),intent(inout)      :: x(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans
    Integer           :: err_act
    character(len=20) :: name='d_prec_apply'

    call psb_erractionsave(err_act)

    select type(prec) 
    type is (mld_dprec_type)
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

  end subroutine mld_d_apply1v


  subroutine mld_d_dump(prec,info,istart,iend,prefix,head,ac,rp,smoother,solver)
    use psb_base_mod
    implicit none 
    class(mld_dprec_type), intent(in) :: prec
    integer, intent(out)             :: info
    integer, intent(in), optional    :: istart, iend
    character(len=*), intent(in), optional :: prefix, head
    logical, optional, intent(in)    :: smoother, solver,ac, rp
    integer :: i, j, il1, iln, lname, lev
    integer :: icontxt,iam, np
    character(len=80)  :: prefix_
    character(len=120) :: fname ! len should be at least 20 more than
    !  len of prefix_ 

    info = 0

    iln = size(prec%precv)
    if (present(istart)) then 
      il1 = max(1,istart)
    else
      il1 = 2
    end if
    if (present(iend)) then 
      iln = min(iln, iend)
    end if

    do lev=il1, iln
      call prec%precv(lev)%dump(lev,info,prefix=prefix,head=head,&
           & ac=ac,smoother=smoother,solver=solver,rp=rp)
    end do

  end subroutine mld_d_dump
  


end module mld_d_prec_type
