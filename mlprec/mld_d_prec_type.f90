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
! File: mld_d_prec_type.f90
!
! Module: mld_d_prec_type
!
!  This module defines: 
!  - the mld_d_prec_type data structure containing the preconditioner and related
!    data structures;
!
!  It contains routines for
!  - Building and applying; 
!  - checking if the preconditioner is correctly defined;
!  - printing a	description of the preconditioner;
!  - deallocating the preconditioner data structure.  
!

module mld_d_prec_type

  use mld_base_prec_type
  use mld_d_base_solver_mod
  use mld_d_base_smoother_mod
  use mld_d_onelev_mod
  use psb_prec_mod, only : psb_dprec_type

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

  type, extends(psb_dprec_type)         :: mld_dprec_type
    integer                             :: ictxt
    integer                             :: coarse_aggr_size
    real(psb_dpk_)                      :: op_complexity=dzero
    type(mld_d_onelev_type), allocatable :: precv(:) 
  contains
    procedure, pass(prec)               :: psb_d_apply2_vect => mld_d_apply2_vect
    procedure, pass(prec)               :: psb_d_apply1_vect => mld_d_apply1_vect
    procedure, pass(prec)               :: psb_d_apply2v => mld_d_apply2v
    procedure, pass(prec)               :: psb_d_apply1v => mld_d_apply1v
    procedure, pass(prec)               :: dump      => mld_d_dump
    procedure, pass(prec)               :: get_complexity => mld_d_get_compl
    procedure, pass(prec)               :: cmp_complexity => mld_d_cmp_compl
    procedure, pass(prec)               :: get_nzeros => mld_d_get_nzeros
    procedure, pass(prec)               :: sizeof => mld_dprec_sizeof
  end type mld_dprec_type

  private :: mld_d_dump, mld_d_get_compl,  mld_d_cmp_compl,&
       &  mld_d_get_nzeros


  !
  ! Interfaces to routines for checking the definition of the preconditioner,
  ! for printing its description and for deallocating its data structure
  !

  interface mld_precfree
    module procedure mld_dprec_free
  end interface


  interface mld_precdescr
    module procedure mld_dfile_prec_descr
  end interface

  interface mld_sizeof
    module procedure mld_dprec_sizeof
  end interface

  interface mld_precaply
    subroutine mld_dprecaply2_vect(prec,x,y,desc_data,info,trans,work)
      import :: psb_dspmat_type, psb_desc_type, &
           & psb_dpk_, psb_d_vect_type, mld_dprec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_dprec_type), intent(inout) :: prec
      type(psb_d_vect_type),intent(inout) :: x
      type(psb_d_vect_type),intent(inout) :: y
      integer, intent(out)                :: info
      character(len=1), optional          :: trans
      real(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine mld_dprecaply2_vect
    subroutine mld_dprecaply1_vect(prec,x,desc_data,info,trans,work)
      import :: psb_dspmat_type, psb_desc_type, &
           & psb_dpk_, psb_d_vect_type, mld_dprec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_dprec_type), intent(inout) :: prec
      type(psb_d_vect_type),intent(inout) :: x
      integer, intent(out)                :: info
      character(len=1), optional          :: trans
      real(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine mld_dprecaply1_vect
    subroutine mld_dprecaply(prec,x,y,desc_data,info,trans,work)
      import :: psb_dspmat_type, psb_desc_type, psb_dpk_, mld_dprec_type
      type(psb_desc_type),intent(in)   :: desc_data
      type(mld_dprec_type), intent(in) :: prec
      real(psb_dpk_),intent(inout)     :: x(:)
      real(psb_dpk_),intent(inout)     :: y(:)
      integer, intent(out)             :: info
      character(len=1), optional       :: trans
      real(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine mld_dprecaply
    subroutine mld_dprecaply1(prec,x,desc_data,info,trans)
      import :: psb_dspmat_type, psb_desc_type, psb_dpk_, mld_dprec_type
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
    class(mld_dprec_type), intent(in) :: prec
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


  !
  ! Subroutines: mld_Tprec_free
  ! Version: real
  !
  !  These routines deallocate the mld_Tprec_type data structures.
  !
  ! Arguments:
  !  p       -  type(mld_Tprec_type), input.
  !             The data structure to be deallocated.
  !  info    -  integer, output.
  !             error code.
  !
  subroutine mld_dprec_free(p,info)
  
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
  ! Top level methods. 
  !
  subroutine mld_d_apply2_vect(prec,x,y,desc_data,info,trans,work)
    implicit none 
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

  subroutine mld_d_apply1_vect(prec,x,desc_data,info,trans,work)
    implicit none 
    type(psb_desc_type),intent(in)        :: desc_data
    class(mld_dprec_type), intent(inout)  :: prec
    type(psb_d_vect_type),intent(inout)   :: x
    integer, intent(out)                  :: info
    character(len=1), optional            :: trans
    real(psb_dpk_),intent(inout), optional, target :: work(:)
    Integer           :: err_act
    character(len=20) :: name='d_prec_apply'

    call psb_erractionsave(err_act)

    select type(prec) 
    type is (mld_dprec_type)
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

  end subroutine mld_d_apply1_vect


  subroutine mld_d_apply2v(prec,x,y,desc_data,info,trans,work)
    implicit none 
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
    implicit none 
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
