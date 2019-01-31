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
! File: mld_c_prec_type.f90
!
! Module: mld_c_prec_type
!
!  This module defines: 
!  - the mld_c_prec_type data structure containing the preconditioner and related
!    data structures;
!
!  It contains routines for
!  - Building and applying; 
!  - checking if the preconditioner is correctly defined;
!  - printing a	description of the preconditioner;
!  - deallocating the preconditioner data structure.  
!

module mld_c_prec_type

  use mld_base_prec_type
  use mld_c_base_solver_mod
  use mld_c_base_smoother_mod
  use mld_c_base_aggregator_mod
  use mld_c_onelev_mod
  use psb_prec_mod, only : psb_cprec_type

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
  !  the finest one and the number of levels is given by size(precv(:)),
  !  and that is the id of the coarsest level. 
  !  In the multigrid literature authors often number the levels in decreasing
  !  order, with level 0 being the id of the coarsest level.
  !
  !
  integer, parameter, private :: wv_size_=4
  
  type, extends(psb_cprec_type)        :: mld_cprec_type
    ! integer(psb_ipk_)                  :: ictxt ! Now it's in the PSBLAS prec. 
    !
    ! Aggregation defaults:
    !
    ! 1. min_coarse_size = 0        Default target size will be computed  as  40*(N_fine)**(1./3.)
    integer(psb_ipk_)                  :: min_coarse_size = izero
    ! 2. maximum number of levels.   Defaults to  20 
    integer(psb_ipk_)                  :: max_levs    = 20_psb_ipk_
    ! 3. min_cr_ratio   = 1.5     
    real(psb_spk_)                     :: min_cr_ratio   = 1.5_psb_spk_
    real(psb_spk_)                     :: op_complexity  = szero
    real(psb_spk_)                     :: avg_cr         = szero
    !
    ! Number of outer sweeps. Sometimes  2 V-cycles may be better than 1 W-cycle. 
    !
    integer(psb_ipk_)                  :: outer_sweeps = 1
    !
    ! Coarse solver requires some tricky checks, and this needs we record the
    ! choice in the format given by the user, to keep track against what
    ! is put later in the multilevel array
    !
    integer(psb_ipk_)                  :: coarse_solver = -1
    
    !
    ! The multilevel hierarchy
    !
    type(mld_c_onelev_type), allocatable :: precv(:) 
  contains
    procedure, pass(prec)               :: psb_c_apply2_vect => mld_c_apply2_vect
    procedure, pass(prec)               :: psb_c_apply1_vect => mld_c_apply1_vect
    procedure, pass(prec)               :: psb_c_apply2v => mld_c_apply2v
    procedure, pass(prec)               :: psb_c_apply1v => mld_c_apply1v
    procedure, pass(prec)               :: dump           => mld_c_dump
    procedure, pass(prec)               :: cnv            => mld_c_cnv
    procedure, pass(prec)               :: clone          => mld_c_clone
    procedure, pass(prec)               :: free           => mld_c_prec_free
    procedure, pass(prec)               :: allocate_wrk   => mld_c_allocate_wrk
    procedure, pass(prec)               :: free_wrk       => mld_c_free_wrk
    procedure, pass(prec)               :: get_complexity => mld_c_get_compl
    procedure, pass(prec)               :: cmp_complexity => mld_c_cmp_compl
    procedure, pass(prec)               :: get_avg_cr => mld_c_get_avg_cr
    procedure, pass(prec)               :: cmp_avg_cr => mld_c_cmp_avg_cr
    procedure, pass(prec)               :: get_nlevs  => mld_c_get_nlevs
    procedure, pass(prec)               :: get_nzeros => mld_c_get_nzeros
    procedure, pass(prec)               :: sizeof => mld_cprec_sizeof
    procedure, pass(prec)               :: setsm  => mld_cprecsetsm
    procedure, pass(prec)               :: setsv  => mld_cprecsetsv
    procedure, pass(prec)               :: setag  => mld_cprecsetag
    procedure, pass(prec)               :: cseti  => mld_ccprecseti
    procedure, pass(prec)               :: csetc  => mld_ccprecsetc
    procedure, pass(prec)               :: csetr  => mld_ccprecsetr
    generic, public                     :: set => cseti, csetc, csetr, setsm, setsv, setag 
    procedure, pass(prec)               :: get_smoother => mld_c_get_smootherp
    procedure, pass(prec)               :: get_solver   => mld_c_get_solverp
    procedure, pass(prec)               :: move_alloc   => c_prec_move_alloc
    procedure, pass(prec)               :: init         => mld_cprecinit
    procedure, pass(prec)               :: build        => mld_cprecbld
    procedure, pass(prec)               :: hierarchy_build  => mld_c_hierarchy_bld
    procedure, pass(prec)               :: smoothers_build  => mld_c_smoothers_bld
    procedure, pass(prec)               :: descr        =>  mld_cfile_prec_descr
  end type mld_cprec_type

  private :: mld_c_dump, mld_c_get_compl,  mld_c_cmp_compl,&
       & mld_c_get_avg_cr,  mld_c_cmp_avg_cr,&
       & mld_c_get_nzeros, mld_c_get_nlevs, c_prec_move_alloc


  !
  ! Interfaces to routines for checking the definition of the preconditioner,
  ! for printing its description and for deallocating its data structure
  !

  interface mld_precfree
    module procedure mld_cprecfree
  end interface


  interface mld_precdescr
    subroutine mld_cfile_prec_descr(prec,iout,root)
      import :: mld_cprec_type, psb_ipk_
      implicit none 
      ! Arguments
      class(mld_cprec_type), intent(in)      :: prec
      integer(psb_ipk_), intent(in), optional :: iout
      integer(psb_ipk_), intent(in), optional :: root
    end subroutine mld_cfile_prec_descr
  end interface

  interface mld_sizeof
    module procedure mld_cprec_sizeof
  end interface

  interface mld_precapply
    subroutine mld_cprecaply2_vect(prec,x,y,desc_data,info,trans,work)
      import :: psb_cspmat_type, psb_desc_type, &
           & psb_spk_, psb_c_vect_type, mld_cprec_type, psb_ipk_
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_cprec_type), intent(inout) :: prec
      type(psb_c_vect_type),intent(inout) :: x
      type(psb_c_vect_type),intent(inout) :: y
      integer(psb_ipk_), intent(out)                :: info
      character(len=1), optional          :: trans
      complex(psb_spk_),intent(inout), optional, target :: work(:)
    end subroutine mld_cprecaply2_vect
    subroutine mld_cprecaply1_vect(prec,x,desc_data,info,trans,work)
      import :: psb_cspmat_type, psb_desc_type, &
           & psb_spk_, psb_c_vect_type, mld_cprec_type, psb_ipk_
      type(psb_desc_type),intent(in)      :: desc_data
      type(mld_cprec_type), intent(inout) :: prec
      type(psb_c_vect_type),intent(inout) :: x
      integer(psb_ipk_), intent(out)                :: info
      character(len=1), optional          :: trans
      complex(psb_spk_),intent(inout), optional, target :: work(:)
    end subroutine mld_cprecaply1_vect
    subroutine mld_cprecaply(prec,x,y,desc_data,info,trans,work)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, mld_cprec_type, psb_ipk_
      type(psb_desc_type),intent(in)   :: desc_data
      type(mld_cprec_type), intent(inout) :: prec
      complex(psb_spk_),intent(inout)     :: x(:)
      complex(psb_spk_),intent(inout)     :: y(:)
      integer(psb_ipk_), intent(out)             :: info
      character(len=1), optional       :: trans
      complex(psb_spk_),intent(inout), optional, target :: work(:)
    end subroutine mld_cprecaply
    subroutine mld_cprecaply1(prec,x,desc_data,info,trans)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, mld_cprec_type, psb_ipk_
      type(psb_desc_type),intent(in)   :: desc_data
      type(mld_cprec_type), intent(inout) :: prec
      complex(psb_spk_),intent(inout)     :: x(:)
      integer(psb_ipk_), intent(out)             :: info
      character(len=1), optional       :: trans
    end subroutine mld_cprecaply1
  end interface

  interface 
    subroutine mld_cprecsetsm(prec,val,info,ilev,ilmax,pos)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, &
           & mld_cprec_type, mld_c_base_smoother_type, psb_ipk_
      class(mld_cprec_type), target, intent(inout):: prec
      class(mld_c_base_smoother_type), intent(in) :: val
      integer(psb_ipk_), intent(out)              :: info
      integer(psb_ipk_), optional, intent(in)     :: ilev,ilmax
      character(len=*), optional, intent(in)      :: pos
    end subroutine mld_cprecsetsm
    subroutine mld_cprecsetsv(prec,val,info,ilev,ilmax,pos)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, &
           & mld_cprec_type, mld_c_base_solver_type, psb_ipk_
      class(mld_cprec_type), intent(inout)      :: prec
      class(mld_c_base_solver_type), intent(in) :: val
      integer(psb_ipk_), intent(out)              :: info
      integer(psb_ipk_), optional, intent(in)     :: ilev,ilmax
      character(len=*), optional, intent(in)      :: pos
    end subroutine mld_cprecsetsv
    subroutine mld_cprecsetag(prec,val,info,ilev,ilmax,pos)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, &
           & mld_cprec_type, mld_c_base_aggregator_type, psb_ipk_
      class(mld_cprec_type), intent(inout)      :: prec
      class(mld_c_base_aggregator_type), intent(in) :: val
      integer(psb_ipk_), intent(out)              :: info
      integer(psb_ipk_), optional, intent(in)     :: ilev,ilmax
      character(len=*), optional, intent(in)      :: pos
    end subroutine mld_cprecsetag
    subroutine mld_ccprecseti(prec,what,val,info,ilev,ilmax,pos,idx)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, &
           & mld_cprec_type, psb_ipk_
      class(mld_cprec_type), intent(inout)   :: prec
      character(len=*), intent(in)             :: what 
      integer(psb_ipk_), intent(in)            :: val
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), optional, intent(in)  :: ilev,ilmax,idx
      character(len=*), optional, intent(in)      :: pos
    end subroutine mld_ccprecseti
    subroutine mld_ccprecsetr(prec,what,val,info,ilev,ilmax,pos,idx)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, &
           & mld_cprec_type, psb_ipk_
      class(mld_cprec_type), intent(inout)   :: prec
      character(len=*), intent(in)             :: what 
      real(psb_spk_), intent(in)                :: val
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), optional, intent(in)  :: ilev,ilmax,idx
      character(len=*), optional, intent(in)      :: pos
    end subroutine mld_ccprecsetr
    subroutine mld_ccprecsetc(prec,what,string,info,ilev,ilmax,pos,idx)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, &
           & mld_cprec_type, psb_ipk_
      class(mld_cprec_type), intent(inout)   :: prec
      character(len=*), intent(in)             :: what 
      character(len=*), intent(in)             :: string
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), optional, intent(in)  :: ilev,ilmax,idx
      character(len=*), optional, intent(in)      :: pos
    end subroutine mld_ccprecsetc
  end interface

  interface mld_precinit
    subroutine mld_cprecinit(ictxt,prec,ptype,info)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, &
           & mld_cprec_type, psb_ipk_
      integer(psb_ipk_), intent(in)            :: ictxt
      class(mld_cprec_type), intent(inout)    :: prec
      character(len=*), intent(in)             :: ptype
      integer(psb_ipk_), intent(out)           :: info
    end subroutine mld_cprecinit
  end interface mld_precinit

  interface mld_precbld
    subroutine mld_cprecbld(a,desc_a,prec,info,amold,vmold,imold)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, &
           & psb_c_base_sparse_mat, psb_c_base_vect_type, &
           & psb_i_base_vect_type, mld_cprec_type, psb_ipk_
      implicit none
      type(psb_cspmat_type), intent(in), target          :: a
      type(psb_desc_type), intent(inout), target           :: desc_a
      class(mld_cprec_type), intent(inout), target       :: prec
      integer(psb_ipk_), intent(out)                       :: info
      class(psb_c_base_sparse_mat), intent(in), optional :: amold
      class(psb_c_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
      !      character, intent(in),optional             :: upd
    end subroutine mld_cprecbld
  end interface mld_precbld

  interface mld_hierarchy_bld
    subroutine mld_c_hierarchy_bld(a,desc_a,prec,info)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, &
           & mld_cprec_type, psb_ipk_
      implicit none
      type(psb_cspmat_type), intent(in), target          :: a
      type(psb_desc_type), intent(inout), target           :: desc_a
      class(mld_cprec_type), intent(inout), target        :: prec
      integer(psb_ipk_), intent(out)                       :: info
      !      character, intent(in),optional             :: upd
    end subroutine mld_c_hierarchy_bld
  end interface mld_hierarchy_bld

  interface mld_smoothers_bld
    subroutine mld_c_smoothers_bld(a,desc_a,prec,info,amold,vmold,imold)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, &
           & psb_c_base_sparse_mat, psb_c_base_vect_type, &
           & psb_i_base_vect_type, mld_cprec_type, psb_ipk_
      implicit none
      type(psb_cspmat_type), intent(in), target          :: a
      type(psb_desc_type), intent(inout), target           :: desc_a
      class(mld_cprec_type), intent(inout), target        :: prec
      integer(psb_ipk_), intent(out)                       :: info
      class(psb_c_base_sparse_mat), intent(in), optional :: amold
      class(psb_c_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
      !      character, intent(in),optional             :: upd
    end subroutine mld_c_smoothers_bld
  end interface mld_smoothers_bld
  
contains
  !
  ! Function returning a pointer to the smoother
  !
  function mld_c_get_smootherp(prec,ilev) result(val)
    implicit none 
    class(mld_cprec_type), target, intent(in) :: prec
    integer(psb_ipk_), optional                 :: ilev
    class(mld_c_base_smoother_type), pointer  :: val
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
  end function mld_c_get_smootherp
  !
  ! Function returning a pointer to the solver
  !
  function mld_c_get_solverp(prec,ilev) result(val)
    implicit none 
    class(mld_cprec_type), target, intent(in) :: prec
    integer(psb_ipk_), optional                 :: ilev
    class(mld_c_base_solver_type), pointer  :: val
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
  end function mld_c_get_solverp
  !
  ! Function returning the size of the precv(:) array
  !
  function mld_c_get_nlevs(prec) result(val)
    implicit none 
    class(mld_cprec_type), intent(in) :: prec
    integer(psb_ipk_) :: val
    val = 0
    if (allocated(prec%precv)) then 
      val = size(prec%precv)
    end if
  end function mld_c_get_nlevs
  !
  ! Function returning the size of the mld_prec_type data structure
  ! in bytes or in number of nonzeros of the operator(s) involved. 
  !
  function mld_c_get_nzeros(prec) result(val)
    implicit none 
    class(mld_cprec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer(psb_ipk_)        :: i
    val = 0
    if (allocated(prec%precv)) then 
      do i=1, size(prec%precv)
        val = val + prec%precv(i)%get_nzeros()
      end do
    end if
  end function mld_c_get_nzeros

  function mld_cprec_sizeof(prec) result(val)
    implicit none 
    class(mld_cprec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer(psb_ipk_)        :: i
    val = 0
    val = val + psb_sizeof_int
    if (allocated(prec%precv)) then 
      do i=1, size(prec%precv)
        val = val + prec%precv(i)%sizeof()
      end do
    end if
  end function mld_cprec_sizeof

  !
  ! Operator complexity: ratio of total number
  ! of nonzeros in the aggregated matrices at the
  ! various level to the nonzeroes at the fine level
  ! (original matrix)
  !
  
  function mld_c_get_compl(prec) result(val)
    implicit none 
    class(mld_cprec_type), intent(in) :: prec
    complex(psb_spk_)  :: val
    
    val = prec%op_complexity

  end function mld_c_get_compl
  
  subroutine mld_c_cmp_compl(prec) 

    implicit none 
    class(mld_cprec_type), intent(inout) :: prec
    
    real(psb_spk_) :: num, den, nmin
    integer(psb_ipk_) :: ictxt 
    integer(psb_ipk_)  :: il 

    num = -sone
    den = sone
    ictxt = prec%ictxt
    if (allocated(prec%precv)) then 
      il  = 1
      num = prec%precv(il)%base_a%get_nzeros()
      if (num >= szero) then
        den = num 
        do il=2,size(prec%precv)
          num = num + max(0,prec%precv(il)%base_a%get_nzeros())
        end do
      end if
    end if
    nmin = num
    call psb_min(ictxt,nmin) 
    if (nmin < szero) then
      num = szero
      den = sone
    else
      call psb_sum(ictxt,num)
      call psb_sum(ictxt,den)
    end if
    prec%op_complexity = num/den
  end subroutine mld_c_cmp_compl
  
  !
  ! Average coarsening ratio
  !
  
  function mld_c_get_avg_cr(prec) result(val)
    implicit none 
    class(mld_cprec_type), intent(in) :: prec
    complex(psb_spk_)  :: val
    
    val = prec%avg_cr

  end function mld_c_get_avg_cr
  
  subroutine mld_c_cmp_avg_cr(prec) 

    implicit none 
    class(mld_cprec_type), intent(inout) :: prec
    
    real(psb_spk_)  :: avgcr
    integer(psb_ipk_) :: ictxt 
    integer(psb_ipk_) :: il, nl, iam, np


    avgcr = szero
    ictxt = prec%ictxt
    call psb_info(ictxt,iam,np)
    if (allocated(prec%precv)) then
      nl = size(prec%precv)
      do il=2,nl
        avgcr = avgcr + max(szero,prec%precv(il)%szratio)
      end do
      avgcr = avgcr / (nl-1)      
    end if
    call psb_sum(ictxt,avgcr) 
    prec%avg_cr = avgcr/np
  end subroutine mld_c_cmp_avg_cr
  
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
  subroutine mld_cprecfree(p,info)
  
    implicit none
    
    ! Arguments
    type(mld_cprec_type), intent(inout) :: p
    integer(psb_ipk_), intent(out)        :: info
    
    ! Local variables
    integer(psb_ipk_)   :: me,err_act,i
    character(len=20)   :: name
    
    if(psb_get_errstatus().ne.0) return 
    info=psb_success_
    name = 'mld_cprecfree'
    call psb_erractionsave(err_act)
    
    me=-1
    
    call p%free(info)

    
    return
    
  end subroutine mld_cprecfree

  subroutine mld_c_prec_free(prec,info)
  
    implicit none
    
    ! Arguments
    class(mld_cprec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out)        :: info
    
    ! Local variables
    integer(psb_ipk_)   :: me,err_act,i
    character(len=20)   :: name
    
    if(psb_get_errstatus().ne.0) return 
    info=psb_success_
    name = 'mld_cprecfree'
    call psb_erractionsave(err_act)
    
    me=-1
    call prec%free_wrk(info)
    if (allocated(prec%precv)) then 
      do i=1,size(prec%precv) 
        call prec%precv(i)%free(info)
      end do
      deallocate(prec%precv,stat=info)
    end if
    call psb_erractionrestore(err_act)
    return
    
9999 call psb_error_handler(err_act)
    return
    
  end subroutine mld_c_prec_free

  

  !
  ! Top level methods. 
  !
  subroutine mld_c_apply2_vect(prec,x,y,desc_data,info,trans,work)
    implicit none 
    type(psb_desc_type),intent(in)        :: desc_data
    class(mld_cprec_type), intent(inout)  :: prec
    type(psb_c_vect_type),intent(inout)   :: x
    type(psb_c_vect_type),intent(inout)   :: y
    integer(psb_ipk_), intent(out)          :: info
    character(len=1), optional              :: trans
    complex(psb_spk_),intent(inout), optional, target :: work(:)
    Integer(psb_ipk_) :: err_act
    character(len=20) :: name='d_prec_apply'

    call psb_erractionsave(err_act)

    select type(prec) 
    type is (mld_cprec_type)
      call mld_precapply(prec,x,y,desc_data,info,trans,work)
    class default
      info = psb_err_missing_override_method_
      call psb_errpush(info,name)
      goto 9999 
    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine mld_c_apply2_vect

  subroutine mld_c_apply1_vect(prec,x,desc_data,info,trans,work)
    implicit none 
    type(psb_desc_type),intent(in)          :: desc_data
    class(mld_cprec_type), intent(inout)  :: prec
    type(psb_c_vect_type),intent(inout)   :: x
    integer(psb_ipk_), intent(out)          :: info
    character(len=1), optional            :: trans
    complex(psb_spk_),intent(inout), optional, target :: work(:)
    Integer(psb_ipk_) :: err_act
    character(len=20) :: name='d_prec_apply'

    call psb_erractionsave(err_act)

    select type(prec) 
    type is (mld_cprec_type)
      call mld_precapply(prec,x,desc_data,info,trans,work)
    class default
      info = psb_err_missing_override_method_
      call psb_errpush(info,name)
      goto 9999 
    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine mld_c_apply1_vect


  subroutine mld_c_apply2v(prec,x,y,desc_data,info,trans,work)
    implicit none 
    type(psb_desc_type),intent(in)    :: desc_data
    class(mld_cprec_type), intent(inout) :: prec
    complex(psb_spk_),intent(inout)      :: x(:)
    complex(psb_spk_),intent(inout)      :: y(:)
    integer(psb_ipk_), intent(out)     :: info
    character(len=1), optional        :: trans
    complex(psb_spk_),intent(inout), optional, target :: work(:)
    Integer(psb_ipk_) :: err_act
    character(len=20) :: name='d_prec_apply'

    call psb_erractionsave(err_act)

    select type(prec) 
    type is (mld_cprec_type)
      call mld_precapply(prec,x,y,desc_data,info,trans,work)
    class default
      info = psb_err_missing_override_method_
      call psb_errpush(info,name)
      goto 9999 
    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine mld_c_apply2v

  subroutine mld_c_apply1v(prec,x,desc_data,info,trans)
    implicit none 
    type(psb_desc_type),intent(in)    :: desc_data
    class(mld_cprec_type), intent(inout) :: prec
    complex(psb_spk_),intent(inout)      :: x(:)
    integer(psb_ipk_), intent(out)     :: info
    character(len=1), optional         :: trans
    Integer(psb_ipk_) :: err_act
    character(len=20) :: name='d_prec_apply'

    call psb_erractionsave(err_act)

    select type(prec) 
    type is (mld_cprec_type)
      call mld_precapply(prec,x,desc_data,info,trans)
    class default
      info = psb_err_missing_override_method_
      call psb_errpush(info,name)
      goto 9999 
    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
  return

  end subroutine mld_c_apply1v


  subroutine mld_c_dump(prec,info,istart,iend,prefix,head,ac,rp,smoother,solver,tprol,&
       & global_num)
    
    implicit none 
    class(mld_cprec_type), intent(in)     :: prec
    integer(psb_ipk_), intent(out)          :: info
    integer(psb_ipk_), intent(in), optional :: istart, iend
    character(len=*), intent(in), optional  :: prefix, head
    logical, optional, intent(in)    :: smoother, solver,ac, rp, tprol, global_num
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
           & ac=ac,smoother=smoother,solver=solver,rp=rp,tprol=tprol, &
           & global_num=global_num)
    end do

  end subroutine mld_c_dump


  subroutine mld_c_cnv(prec,info,amold,vmold,imold)

    implicit none 
    class(mld_cprec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out)       :: info
    class(psb_c_base_sparse_mat), intent(in), optional :: amold
    class(psb_c_base_vect_type), intent(in), optional  :: vmold
    class(psb_i_base_vect_type), intent(in), optional  :: imold

    integer(psb_ipk_) :: i
    
    info = psb_success_
    if (allocated(prec%precv)) then
      do i=1,size(prec%precv)
        if (info == psb_success_ ) &
             & call prec%precv(i)%cnv(info,amold=amold,vmold=vmold,imold=imold)
      end do
    end if
    
  end subroutine mld_c_cnv

  subroutine mld_c_clone(prec,precout,info)

    implicit none 
    class(mld_cprec_type), intent(inout) :: prec
    class(psb_cprec_type), intent(inout) :: precout
    integer(psb_ipk_), intent(out)       :: info
    
    call precout%free(info)
    if (info == 0) call mld_c_inner_clone(prec,precout,info) 

  end subroutine mld_c_clone

  subroutine mld_c_inner_clone(prec,precout,info)

    implicit none 
    class(mld_cprec_type), intent(inout)         :: prec
    class(psb_cprec_type), target, intent(inout) :: precout
    integer(psb_ipk_), intent(out)             :: info
    ! Local vars
    integer(psb_ipk_)  :: i, j, il1, ln, lname, lev
    integer(psb_ipk_)  :: icontxt,iam, np

    info = psb_success_
    select type(pout => precout)
    class is (mld_cprec_type)
      pout%ictxt            = prec%ictxt
      pout%max_levs         = prec%max_levs
      pout%min_coarse_size  = prec%min_coarse_size
      pout%min_cr_ratio     = prec%min_cr_ratio
      pout%outer_sweeps     = prec%outer_sweeps
      pout%op_complexity    = prec%op_complexity
      pout%avg_cr           = prec%avg_cr
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
      if (allocated(prec%precv(1)%wrk)) &
           & call pout%allocate_wrk(info,vmold=prec%precv(1)%wrk%vx2l%v)

    class default 
      write(0,*) 'Error: wrong out type'
      info = psb_err_invalid_input_
    end select
9999 continue
  end subroutine mld_c_inner_clone

  subroutine c_prec_move_alloc(prec, b,info)
    use psb_base_mod
    implicit none
    class(mld_cprec_type), intent(inout) :: prec
    class(mld_cprec_type), intent(inout), target :: b
    integer(psb_ipk_), intent(out) :: info 
    integer(psb_ipk_) :: i
    
    if (same_type_as(prec,b)) then 
      if (allocated(b%precv)) then 
        ! This might not be required if FINAL procedures are available.
        call b%free(info)
        if (info /= psb_success_) then 
          !?????
!!$        return
        endif
      end if
      
      call move_alloc(prec%precv,b%precv)
      ! Fix the pointers except on level 1.
      do i=2, size(b%precv)
        b%precv(i)%base_a    => b%precv(i)%ac
        b%precv(i)%base_desc => b%precv(i)%desc_ac
        b%precv(i)%map%p_desc_X => b%precv(i-1)%base_desc
        b%precv(i)%map%p_desc_Y => b%precv(i)%base_desc
      end do
            
    else
      write(0,*) 'Warning: PREC%move_alloc onto different type?'
      info = psb_err_internal_error_
    end if
  end subroutine c_prec_move_alloc

  subroutine mld_c_allocate_wrk(prec,info,vmold)
    use psb_base_mod
    implicit none
    
    ! Arguments
    class(mld_cprec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out)        :: info
    class(psb_c_base_vect_type), intent(in), optional  :: vmold

    ! Local variables
    integer(psb_ipk_)   :: me,err_act,i,j,level,nlev, nc2l
    character(len=20)   :: name
    
    if(psb_get_errstatus().ne.0) return 
    info=psb_success_
    name = 'mld_c_allocate_wrk'
    call psb_erractionsave(err_act)
    nlev   = size(prec%precv)  
    level = 1
    do level = 1, nlev
      call prec%precv(level)%allocate_wrk(info,vmold=vmold)
      if (psb_errstatus_fatal()) then 
        nc2l = prec%precv(level)%base_desc%get_local_cols()
        info=psb_err_alloc_request_
        call psb_errpush(info,name,i_err=(/2*nc2l,izero,izero,izero,izero/),&
             & a_err='complex(psb_spk_)')
        goto 9999      
      end if
    end do

    call psb_erractionrestore(err_act)
    return
    
9999 call psb_error_handler(err_act)
    return
    
  end subroutine mld_c_allocate_wrk



  subroutine mld_c_free_wrk(prec,info)
    use psb_base_mod
    implicit none

    ! Arguments
    class(mld_cprec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out)        :: info

    ! Local variables
    integer(psb_ipk_)   :: me,err_act,i,j,level, nlev, nc2l
    character(len=20)   :: name

    if(psb_get_errstatus().ne.0) return 
    info=psb_success_
    name = 'mld_c_free_wrk'
    call psb_erractionsave(err_act)

    nlev   = size(prec%precv)  
    do level = 1, nlev
      call prec%precv(level)%free_wrk(info)
    end do

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine mld_c_free_wrk

end module mld_c_prec_type
