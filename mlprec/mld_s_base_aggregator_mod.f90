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
!  The aggregator object hosts the aggregation method for building
!  the multilevel hierarchy. 
!    
module mld_s_base_aggregator_mod

  use mld_base_prec_type, only : mld_sml_parms, mld_saggr_data 
  use psb_base_mod, only : psb_sspmat_type, psb_s_vect_type, &
       & psb_s_base_vect_type, psb_slinmap_type, psb_spk_, &
       & psb_ipk_, psb_long_int_k_, psb_desc_type, psb_i_base_vect_type, &
       & psb_erractionsave, psb_error_handler, psb_success_, psb_toupper
  !
  !  
  !
  !> \class mld_s_base_aggregator_type
  !!
  !!  It is the data type containing the basic interface definition for
  !!  building a multigrid hierarchy by aggregation. The base object has no attributes,
  !!  it is intended to be essentially an abstract type. 
  !!  
  !!
  !!  type mld_s_base_aggregator_type
  !!  end type
  !!  
  !!  
  !!  Methods:
  !!
  !!  bld_tprol   -   Build a tentative prolongator
  !!  
  !!  mat_bld     -   Build prolongator/restrictor and coarse matrix ac
  !!  
  !!  mat_asb     -   Convert  prolongator/restrictor/coarse matrix
  !!                  and fix their descriptor(s)
  !!                  
  !!  update_next -   Transfer information to the next level; default is
  !!                  to do nothing, i.e. aggregators at different
  !!                  levels are independent.
  !!  
  !!  default     -   Apply defaults
  !!  set_aggr_type - For aggregator that have internal options.
  !!  fmt         - Return a short string description
  !!  descr       - Print a more detailed description
  !!
  !!  cseti, csetr, csetc  - Set internal parameters, if any
  !  
  type mld_s_base_aggregator_type
    ! Do we want to purge explicit zeros when aggregating? 
    logical :: do_clean_zeros
  contains
    procedure, pass(ag) :: bld_tprol   => mld_s_base_aggregator_build_tprol
    procedure, pass(ag) :: mat_bld     => mld_s_base_aggregator_mat_bld
    procedure, pass(ag) :: mat_asb     => mld_s_base_aggregator_mat_asb
    procedure, pass(ag) :: bld_map     => mld_s_base_aggregator_bld_map
    procedure, pass(ag) :: update_next => mld_s_base_aggregator_update_next
    procedure, pass(ag) :: clone       => mld_s_base_aggregator_clone
    procedure, pass(ag) :: free        => mld_s_base_aggregator_free
    procedure, pass(ag) :: default     => mld_s_base_aggregator_default
    procedure, pass(ag) :: descr       => mld_s_base_aggregator_descr
    procedure, pass(ag) :: sizeof      => mld_s_base_aggregator_sizeof   
    procedure, pass(ag) :: set_aggr_type => mld_s_base_aggregator_set_aggr_type
    procedure, nopass   :: fmt         => mld_s_base_aggregator_fmt
    procedure, pass(ag) :: cseti       => mld_s_base_aggregator_cseti
    procedure, pass(ag) :: csetr       => mld_s_base_aggregator_csetr
    procedure, pass(ag) :: csetc       => mld_s_base_aggregator_csetc
    generic, public     :: set         => cseti, csetr, csetc
    procedure, nopass   :: xt_desc     => mld_s_base_aggregator_xt_desc
    procedure, pass(ag) :: backfix     => mld_s_base_aggregator_backfix
  end type mld_s_base_aggregator_type

  abstract interface  
    subroutine mld_s_soc_map_bld(iorder,theta,clean_zeros,a,desc_a,nlaggr,ilaggr,info)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      implicit none 
      integer(psb_ipk_), intent(in)     :: iorder
      logical, intent(in)               :: clean_zeros
      type(psb_sspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)    :: desc_a
      real(psb_spk_), intent(in)         :: theta
      integer(psb_ipk_), allocatable, intent(out)  :: ilaggr(:),nlaggr(:)
      integer(psb_ipk_), intent(out)               :: info
    end subroutine mld_s_soc_map_bld
  end interface

contains

  subroutine mld_s_base_aggregator_cseti(ag,what,val,info,idx)

    Implicit None

    ! Arguments
    class(mld_s_base_aggregator_type), intent(inout) :: ag
    character(len=*), intent(in)                  :: what
    integer(psb_ipk_), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), intent(in), optional       :: idx
    ! Do nothing
    info = 0
  end subroutine mld_s_base_aggregator_cseti

  subroutine mld_s_base_aggregator_csetr(ag,what,val,info,idx)

    Implicit None

    ! Arguments
    class(mld_s_base_aggregator_type), intent(inout) :: ag
    character(len=*), intent(in)                  :: what
    real(psb_spk_), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), intent(in), optional       :: idx
    ! Do nothing
    info = 0
  end subroutine mld_s_base_aggregator_csetr

  subroutine mld_s_base_aggregator_csetc(ag,what,val,info,idx)

    Implicit None

    ! Arguments
    class(mld_s_base_aggregator_type), intent(inout) :: ag
    character(len=*), intent(in)                  :: what
    character(len=*), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), intent(in), optional       :: idx
    ! Set clean zeros, or do nothing. 
    select case (psb_toupper(trim(what)))
    case('AGGR_CLEAN_ZEROS')
      select case (psb_toupper(trim(val)))
      case('TRUE','T')
        ag%do_clean_zeros = .true.
      case('FALSE','F')
        ag%do_clean_zeros = .false.
      end select
    end select
    info = 0
  end subroutine mld_s_base_aggregator_csetc

  
  subroutine  mld_s_base_aggregator_update_next(ag,agnext,info)
    implicit none 
    class(mld_s_base_aggregator_type), target, intent(inout) :: ag, agnext
    integer(psb_ipk_), intent(out)       :: info

    !
    ! Base version does nothing. 
    !
    info = 0 
  end subroutine mld_s_base_aggregator_update_next
  
  subroutine  mld_s_base_aggregator_clone(ag,agnext,info)
    implicit none 
    class(mld_s_base_aggregator_type), intent(inout) :: ag
    class(mld_s_base_aggregator_type), allocatable, intent(inout) :: agnext
    integer(psb_ipk_), intent(out)       :: info

    info = 0 
    if (allocated(agnext)) then
      call agnext%free(info)
      if (info == 0) deallocate(agnext,stat=info)
    end if
    if (info /= 0) return
    allocate(agnext,source=ag,stat=info)
    
  end subroutine mld_s_base_aggregator_clone

  subroutine  mld_s_base_aggregator_free(ag,info)
    implicit none 
    class(mld_s_base_aggregator_type), intent(inout) :: ag
    integer(psb_ipk_), intent(out)       :: info
    
    info = psb_success_
    return
  end subroutine mld_s_base_aggregator_free
  
  subroutine  mld_s_base_aggregator_default(ag)
    implicit none 
    class(mld_s_base_aggregator_type), intent(inout) :: ag
    ! Only one default setting
    ag%do_clean_zeros = .true.
    
    return
  end subroutine mld_s_base_aggregator_default

  function mld_s_base_aggregator_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "Default aggregator "
  end function mld_s_base_aggregator_fmt

  function mld_s_base_aggregator_sizeof(ag) result(val)
    implicit none
    class(mld_s_base_aggregator_type), intent(in)  :: ag
    integer(psb_long_int_k_)  :: val

    val = 1
  end function mld_s_base_aggregator_sizeof


  function mld_s_base_aggregator_xt_desc() result(val)
    implicit none 
    logical  :: val

    val = .false.
  end function mld_s_base_aggregator_xt_desc

  subroutine  mld_s_base_aggregator_backfix(ag,base_a,ac,base_desc,desc_ac,info)
    implicit none
    class(mld_s_base_aggregator_type), intent(inout)  :: ag
    type(psb_sspmat_type), pointer               :: base_a
    type(psb_sspmat_type), intent(inout), target :: ac
    type(psb_desc_type), pointer                   :: base_desc
    type(psb_desc_type), intent(inout), target     :: desc_ac
    integer(psb_ipk_), intent(out)      :: info

    ! Base version is to do nothing.
    info = psb_success_
    
  end subroutine mld_s_base_aggregator_backfix
  
  subroutine  mld_s_base_aggregator_descr(ag,parms,iout,info)
    implicit none 
    class(mld_s_base_aggregator_type), intent(in) :: ag
    type(mld_sml_parms), intent(in)   :: parms
    integer(psb_ipk_), intent(in)  :: iout
    integer(psb_ipk_), intent(out) :: info

    write(iout,*) 'Aggregator object type: ',ag%fmt()
    call parms%mldescr(iout,info)
    
    return
  end subroutine mld_s_base_aggregator_descr
  
  subroutine  mld_s_base_aggregator_set_aggr_type(ag,parms,info)
    implicit none 
    class(mld_s_base_aggregator_type), intent(inout) :: ag
    type(mld_sml_parms), intent(in)   :: parms
    integer(psb_ipk_), intent(out) :: info

    ! Do nothing
    
    return
  end subroutine mld_s_base_aggregator_set_aggr_type

  !
  !> Function  bld_tprol:
  !! \memberof  mld_s_base_aggregator_type
  !! \brief  Build a tentative prolongator.           
  !!         The routine will map the local matrix entries to aggregates.
  !!         The mapping is store in ILAGGR; for each local row index I,      
  !!         ILAGGR(I) contains the index of the aggregate to which index I
  !!         will contribute, in global numbering. 
  !!         Many aggregation produce a binary tentative prolongator, but some
  !!         do not, hence we also need the OP_PROL output.
  !!         AG_DATA is passed here just in case some of the
  !!         aggregators need it internally, most of them will ignore.
  !!         
  !!  \param ag      The input aggregator object
  !!  \param parms   The auxiliary parameters object
  !!  \param ag_data Auxiliary global aggregation info
  !!  \param a       The local matrix part
  !!  \param desc_a  The descriptor
  !!  \param ilaggr  Output aggregation map
  !!  \param nlaggr  Sizes of ilaggr on all processes
  !!  \param op_prol The tentative prolongator operator
  !!  \param info    Return code
  !!           
  !
  subroutine  mld_s_base_aggregator_build_tprol(ag,parms,ag_data,&
       & a,desc_a,ilaggr,nlaggr,op_prol,info)
    use psb_base_mod
    implicit none
    class(mld_s_base_aggregator_type), target, intent(inout) :: ag
    type(mld_sml_parms), intent(inout)  :: parms
    type(mld_saggr_data), intent(in)    :: ag_data
    type(psb_sspmat_type), intent(inout) :: a
    type(psb_desc_type), intent(inout)     :: desc_a
    integer(psb_ipk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
    type(psb_sspmat_type), intent(out)  :: op_prol
    integer(psb_ipk_), intent(out)      :: info
    
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='s_base_aggregator_build_tprol'
    
    call psb_erractionsave(err_act)
    
    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return
    
9999 call psb_error_handler(err_act)
    
    return

  end subroutine mld_s_base_aggregator_build_tprol

  !
  !> Function   mat_bld    
  !! \memberof  mld_s_base_aggregator_type
  !! \brief     Build prolongator/restrictor/coarse matrix.
  !!
  !!
  !!  \param ag        The input aggregator object
  !!  \param parms     The auxiliary parameters object
  !!  \param a         The local matrix part
  !!  \param desc_a    The descriptor
  !!  \param ilaggr    Aggregation map
  !!  \param nlaggr    Sizes of ilaggr on all processes
  !!  \param ac        On output the coarse matrix
  !!  \param op_prol   On input, the  tentative prolongator operator, on output
  !!                   the final prolongator
  !!  \param op_restr  On output, the restrictor operator;
  !!                   in many cases it is the transpose of the prolongator. 
  !!  \param info    Return code
  !!  
  subroutine  mld_s_base_aggregator_mat_bld(ag,parms,a,desc_a,ilaggr,nlaggr,ac,&
       & op_prol,op_restr,info)
    use psb_base_mod
    implicit none
    class(mld_s_base_aggregator_type), target, intent(inout) :: ag
    type(mld_sml_parms), intent(inout)   :: parms 
    type(psb_sspmat_type), intent(in)    :: a
    type(psb_desc_type), intent(in)      :: desc_a
    integer(psb_ipk_), intent(inout)     :: ilaggr(:), nlaggr(:)
    type(psb_sspmat_type), intent(inout)   :: op_prol
    type(psb_sspmat_type), intent(out)   :: ac,op_restr
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='s_base_aggregator_mat_bld'

    call psb_erractionsave(err_act)

    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  end subroutine mld_s_base_aggregator_mat_bld


  !
  !> Function   mat_asb    
  !! \memberof  mld_s_base_aggregator_type
  !! \brief     Build prolongator/restrictor/coarse matrix.
  !!
  !!
  !!  \param ag        The input aggregator object
  !!  \param parms     The auxiliary parameters object
  !!  \param a         The local matrix part
  !!  \param desc_a    The descriptor
  !!  \param ilaggr    Aggregation map
  !!  \param nlaggr    Sizes of ilaggr on all processes
  !!  \param ac        On output the coarse matrix
  !!  \param op_prol   On input, the  tentative prolongator operator, on output
  !!                   the final prolongator
  !!  \param op_restr  On output, the restrictor operator;
  !!                   in many cases it is the transpose of the prolongator. 
  !!  \param info    Return code
  !!  
  subroutine  mld_s_base_aggregator_mat_asb(ag,parms,a,desc_a,ilaggr,nlaggr,&
       & ac,desc_ac, op_prol,op_restr,info)
    use psb_base_mod
    implicit none
    class(mld_s_base_aggregator_type), target, intent(inout) :: ag
    type(mld_sml_parms), intent(inout)  :: parms 
    type(psb_sspmat_type), intent(in)    :: a
    type(psb_desc_type), intent(in)        :: desc_a
    integer(psb_ipk_), intent(inout)       :: ilaggr(:), nlaggr(:)
    type(psb_sspmat_type), intent(inout) :: op_prol, ac,op_restr
    type(psb_desc_type), intent(inout)     :: desc_ac
    integer(psb_ipk_), intent(out)         :: info
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='s_base_aggregator_mat_asb'

    call psb_erractionsave(err_act)

    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  end subroutine mld_s_base_aggregator_mat_asb

  !
  !> Function   bld_map
  !! \memberof  mld_s_base_aggregator_type
  !! \brief     Build linear map between hierarchy levels 
  !!
  !!
  !!  \param ag        The input aggregator object
  !!  \param desc_a    The fine space descriptor
  !!  \param desc_ac   The coarse space descriptor
  !!  \param ilaggr    Aggregation map vector
  !!  \param nlaggr    Sizes of ilaggr on all processes
  !!  \param op_prol   The  prolongator operator
  !!  \param op_restr  The restrictor operator
  !!  \param map       The output map
  !!  \param info    Return code
  !!  
  subroutine  mld_s_base_aggregator_bld_map(ag,desc_a,desc_ac,ilaggr,nlaggr,&
       & op_restr,op_prol,map,info)
    use psb_base_mod
    implicit none
class(mld_s_base_aggregator_type), target, intent(inout) :: ag
    type(psb_desc_type), intent(in), target :: desc_a, desc_ac
    integer(psb_ipk_), intent(inout)     :: ilaggr(:), nlaggr(:)
    type(psb_sspmat_type), intent(inout)   :: op_restr, op_prol
    type(psb_slinmap_type), intent(out)    :: map
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='s_base_aggregator_bld_map'

    call psb_erractionsave(err_act)

    !
    ! Copy the prolongation/restriction matrices into the descriptor map.
    !  op_restr => PR^T   i.e. restriction  operator
    !  op_prol => PR     i.e. prolongation operator
    !
    !  This default implementation reuses desc_a/desc_ac through
    !  pointers in the map structure.
    !      
    map = psb_linmap(psb_map_aggr_,desc_a,&
         & desc_ac,op_restr,op_prol,ilaggr,nlaggr)
    if (info == psb_success_) call op_prol%free()
    if (info == psb_success_) call op_restr%free()
    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='sp_Free')
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  end subroutine mld_s_base_aggregator_bld_map


end module mld_s_base_aggregator_mod
