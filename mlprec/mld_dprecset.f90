!!$ 
!!$ 
!!$                                MLD2P4
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS v.2.0)
!!$  
!!$  (C) Copyright 2007  Alfredo Buttari      University of Rome Tor Vergata
!!$                      Pasqua D'Ambra       ICAR-CNR, Naples
!!$                      Daniela di Serafino  Second University of Naples
!!$                      Salvatore Filippone  University of Rome Tor Vergata    
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
! File: mld_dprecset.f90
!
! Subroutine: mld_dprecseti
!
!  These routines set the parameters defining the preconditioner. More precisely,
!  the parameter identified by 'what' is assigned the value contained in 'val'.
!  mld_dprecsetc works on string parameters, mld_dprecseti works on integer
!  parameters, while mld_dprecsetd works on real ones.
!  For the multilevel preconditioners, the levels are numbered in increasing
!  order starting from the finest one, i.e. level 1 is the finest level. 
!
!
! Arguments:
!    p       -  type(mld_dprec_type), input/output.
!               The preconditioner data structure.
!    what    -  integer, input.
!               The number identifying the parameter to be set.
!               A mnemonic constant has been associated to each of these
!               numbers, as reported in MLD2P4 user's guide.
!    val     -  integer, input.
!               The value of the parameter to be set. The list of allowed
!               values is reported in MLD2P4 user's      guide.
!    info    -  integer, output.
!               Error code.
!    ilev    -  integer, optional, input.
!               For the multilevel preconditioner, the level at which the
!               preconditioner parameter has to be set. 
!               If nlev is not present, the parameter identified by 'what'
!               is set at all the levels but the coarsest one, except when
!               'what' has the values mld_coarse_mat_, mld_coarse_solve_,
!               mld_coarse_sweeps_, mld_coarse_fill_in_, which refer to the
!               coarsest level.
!   
subroutine mld_dprecseti(p,what,val,info,ilev)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_dprecseti

  implicit none

! Arguments
  type(mld_dprec_type), intent(inout)    :: p
  integer, intent(in)                    :: what 
  integer, intent(in)                    :: val
  integer, intent(out)                   :: info
  integer, optional, intent(in)          :: ilev

! Local variables
  integer                                :: ilev_, nlev_

  info = 0

  if (.not.allocated(p%baseprecv)) then 
    info = 3111
    return 
  endif
  nlev_ = size(p%baseprecv)

  if (present(ilev)) then 
    ilev_ = ilev
  else
    ilev_ = 1 
  end if

  if ((ilev_<1).or.(ilev_ > nlev_)) then 
    info = -1
    return
  endif
  if (.not.allocated(p%baseprecv(ilev_)%iprcparm)) then 
    info = 3111
    return 
  endif

  !
  ! Set preconditioner parameters at level ilev.
  !
  if (present(ilev)) then 
    
    if (ilev_ == 1) then
      ! 
      ! Rules for fine level are slightly different.
      ! 
      select case(what) 
      case(mld_prec_type_,mld_sub_solve_,mld_sub_restr_,mld_sub_prol_,mld_sub_ren_,mld_n_ovr_,mld_sub_fill_in_,mld_smooth_sweeps_)
        p%baseprecv(ilev_)%iprcparm(what)  = val
      case default
        write(0,*) 'Error: trying to call PRECSET with an invalid WHAT'
        info = -2
      end select

    else if (ilev_ > 1) then 
      select case(what) 
      case(mld_prec_type_,mld_sub_solve_,mld_sub_restr_,mld_sub_prol_,mld_sub_ren_,mld_n_ovr_,mld_sub_fill_in_,&
           & mld_smooth_sweeps_,mld_ml_type_,mld_aggr_alg_,mld_aggr_kind_,mld_coarse_mat_,&
           & mld_smooth_pos_,mld_aggr_eig_)
        p%baseprecv(ilev_)%iprcparm(what)  = val
      case(mld_coarse_solve_)
        if (ilev_ /= nlev_) then 
          write(0,*) 'Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if
        p%baseprecv(ilev_)%iprcparm(mld_sub_solve_)  = val
      case(mld_coarse_sweeps_)
        if (ilev_ /= nlev_) then 
          write(0,*) 'Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if
        p%baseprecv(ilev_)%iprcparm(mld_smooth_sweeps_)  = val
      case(mld_coarse_fill_in_)
        if (ilev_ /= nlev_) then 
          write(0,*) 'Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if
        p%baseprecv(ilev_)%iprcparm(mld_sub_fill_in_)  = val
      case default
        write(0,*) 'Error: trying to call PRECSET with an invalid WHAT'
        info = -2
      end select

    endif

  else if (.not.present(ilev)) then 
    !
    ! ilev not specified: set preconditioner parameters at levels 1,...,nlev_-1,
    ! except when 'what' has the values mld_coarse_solve_, mld_coarse_sweeps_,
    ! mld_coarse_fill_in_, which refer to the coarsest level
    !

      select case(what) 
      case(mld_prec_type_,mld_sub_solve_,mld_sub_restr_,mld_sub_prol_,mld_sub_ren_,mld_n_ovr_,mld_sub_fill_in_,&
           & mld_smooth_sweeps_,mld_ml_type_,mld_aggr_alg_,mld_aggr_kind_,mld_coarse_mat_,&
           & mld_smooth_pos_,mld_aggr_eig_)
        do ilev_=1,nlev_-1
          if (.not.allocated(p%baseprecv(ilev_)%iprcparm)) then 
            write(0,*) 'Error: trying to call PRECSET on an uninitialized preconditioner'
            info = -1 
            return 
          endif
          p%baseprecv(ilev_)%iprcparm(what)  = val
        end do
      case(mld_coarse_solve_)
        if (.not.allocated(p%baseprecv(nlev_)%iprcparm)) then 
          write(0,*) 'Error: trying to call PRECSET on an uninitialized preconditioner'
          info = -1 
          return 
        endif
        p%baseprecv(nlev_)%iprcparm(mld_sub_solve_)  = val
      case(mld_coarse_sweeps_)
        if (.not.allocated(p%baseprecv(nlev_)%iprcparm)) then 
          write(0,*) 'Error: trying to call PRECSET on an uninitialized preconditioner'
          info = -1 
          return 
        endif
        p%baseprecv(nlev_)%iprcparm(mld_smooth_sweeps_)  = val
      case(mld_coarse_fill_in_)
        if (.not.allocated(p%baseprecv(nlev_)%iprcparm)) then 
          write(0,*) 'Error: trying to call PRECSET on an uninitialized preconditioner'
          info = -1 
          return 
        endif
        p%baseprecv(nlev_)%iprcparm(mld_sub_fill_in_)  = val
      case default
        write(0,*) 'Error: trying to call PRECSET with an invalid WHAT'
        info = -2
      end select

  endif

end subroutine mld_dprecseti

!
! Subroutine: mld_dprecsetc
!
!
! Arguments:
!    p       -  type(mld_dprec_type), input/output.
!               The preconditioner data structure.
!    what    -  integer, input.
!               The number identifying the parameter to be set.
!               A mnemonic constant has been associated to each of these
!               numbers, as reported in MLD2P4 user's guide.
!    val     -  string(len=*) input.
!               The value of the parameter to be set. The list of allowed
!               values is reported in MLD2P4 user's      guide.
!    info    -  integer, output.
!               Error code.
!    ilev    -  integer, optional, input.
!               For the multilevel preconditioner, the level at which the
!               preconditioner parameter has to be set. 
!               If nlev is not present, the parameter identified by 'what'
!               is set at all the levels but the coarsest one, except when
!               'what' has the values mld_coarse_mat_, mld_coarse_solve_,
!               mld_coarse_sweeps_, mld_coarse_fill_in_, which refer to the
!               coarsest level.
!   
subroutine mld_dprecsetc(p,what,string,info,ilev)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_dprecsetc

  implicit none
  type(mld_dprec_type), intent(inout)    :: p
  integer, intent(in)                    :: what 
  character(len=*), intent(in)           :: string
  integer, intent(out)                   :: info
  integer, optional, intent(in)          :: ilev
  integer                                :: ilev_, nlev_,val

  info = 0

  if (.not.allocated(p%baseprecv)) then 
    info = 3111
    return 
  endif
  nlev_ = size(p%baseprecv)

  if (present(ilev)) then 
    ilev_ = ilev
  else
    ilev_ = 1 
  end if

  if ((ilev_<1).or.(ilev_ > nlev_)) then 
    write(0,*) 'PRECSET ERRROR: ilev out of bounds'
    info = -1
    return
  endif
  if (.not.allocated(p%baseprecv(ilev_)%iprcparm)) then 
    info = 3111
    return 
  endif


  if (present(ilev)) then 
    
    if (ilev_ == 1) then 
      ! Rules for fine level are slightly different. 
      select case(what) 
      case(mld_prec_type_,mld_sub_solve_,mld_sub_restr_,mld_sub_prol_)
        call get_stringval(string,val,info)
        p%baseprecv(ilev_)%iprcparm(what)  = val
      case default
        write(0,*) 'Error: trying to call PRECSET with an invalid WHAT'
        info = -2
      end select
    else if (ilev_ > 1) then 
      select case(what) 
      case(mld_prec_type_,mld_sub_solve_,mld_sub_restr_,mld_sub_prol_,&
           & mld_ml_type_,mld_aggr_alg_,mld_aggr_kind_,mld_coarse_mat_,&
           & mld_smooth_pos_,mld_aggr_eig_)
        call get_stringval(string,val,info)
        p%baseprecv(ilev_)%iprcparm(what)  = val
      case(mld_coarse_solve_)
        call get_stringval(string,val,info)
        if (ilev_ /= nlev_) then 
          write(0,*) 'Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if
        p%baseprecv(ilev_)%iprcparm(mld_sub_solve_)  = val
      case default
        write(0,*) 'Error: trying to call PRECSET with an invalid WHAT'
        info = -2
      end select
    endif

  else if (.not.present(ilev)) then 

      select case(what) 
      case(mld_prec_type_,mld_sub_solve_,mld_sub_restr_,mld_sub_prol_,mld_sub_ren_,&
           & mld_smooth_sweeps_,mld_ml_type_,mld_aggr_alg_,mld_aggr_kind_,mld_coarse_mat_,&
           & mld_smooth_pos_)
        call get_stringval(string,val,info)
        do ilev_=1,nlev_-1
          if (.not.allocated(p%baseprecv(ilev_)%iprcparm)) then 
            write(0,*) 'Error: trying to call PRECSET on an uninitialized preconditioner'
            info = -1 
            return 
          endif
          p%baseprecv(ilev_)%iprcparm(what)  = val
        end do
      case(mld_coarse_solve_)
        if (.not.allocated(p%baseprecv(nlev_)%iprcparm)) then 
          write(0,*) 'Error: trying to call PRECSET on an uninitialized preconditioner'
          info = -1 
          return 
        endif
        call get_stringval(string,val,info)
        p%baseprecv(nlev_)%iprcparm(mld_sub_solve_)  = val
      case default
        write(0,*) 'Error: trying to call PRECSET with an invalid WHAT'
        info = -2
      end select

  endif

contains

  subroutine get_stringval(string,val,info)
    character(len=*), intent(in) :: string
    integer, intent(out) :: val, info
    
    info = 0
    select case(toupper(trim(string)))
    case('NONE')
      val = 0
    case('HALO')
      val = psb_halo_ 
    case('SUM')
      val = psb_sum_
    case('AVG')
      val = psb_avg_
    case('ILU')
      val = mld_ilu_n_
    case('MILU')
      val = mld_milu_n_
    case('ILUT')
      val = mld_ilu_t_
    case('SLU')
      val = mld_slu_
    case('UMFP')
      val = mld_umf_
    case('ADD')
      val = mld_add_ml_
    case('MULT')
      val = mld_mult_ml_
    case('DEC')
      val = mld_dec_aggr_
    case('REPL')
      val = mld_repl_mat_
    case('DIST')
      val = mld_distr_mat_
    case('SYMDEC')
      val = mld_sym_dec_aggr_
    case('GLB')
      val = mld_glb_aggr_
    case('SMOOTH')
      val = mld_smooth_prol_
    case('PRE')
      val = mld_pre_smooth_
    case('POST')
      val = mld_post_smooth_
    case('TWOSIDE','BOTH')
      val = mld_twoside_smooth_
    case default
      val  = -1
      info = -1
    end select
    if (info /= 0) then 
      write(0,*) 'Error in get_Stringval: unknown: "',trim(string),'"'
    end if
  end subroutine get_stringval
end subroutine mld_dprecsetc


!
! Subroutine: mld_dprecsetd
!
!
! Arguments:
!    p       -  type(mld_dprec_type), input/output.
!               The preconditioner data structure.
!    what    -  integer, input.
!               The number identifying the parameter to be set.
!               A mnemonic constant has been associated to each of these
!               numbers, as reported in MLD2P4 user's guide.
!    val     -  real(kind(1.d0))
!               The value of the parameter to be set. The list of allowed
!               values is reported in MLD2P4 user's guide.
!    info    -  integer, output.
!               Error code.
!    ilev    -  integer, optional, input.
!               For the multilevel preconditioner, the level at which the
!               preconditioner parameter has to be set. 
!               If nlev is not present, the parameter identified by 'what'
!               is set at all the levels but the coarsest one, except when
!               'what' has the values mld_coarse_mat_, mld_coarse_solve_,
!               mld_coarse_sweeps_, mld_coarse_fill_in_, which refer to the
!               coarsest level.
!   
subroutine mld_dprecsetd(p,what,val,info,ilev)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_dprecsetd

  implicit none
  type(mld_dprec_type), intent(inout)    :: p
  integer, intent(in)                    :: what 
  real(kind(1.d0)), intent(in)           :: val
  integer, intent(out)                   :: info
  integer, optional, intent(in)          :: ilev

  integer                                :: ilev_,nlev_

  info = 0

  if (present(ilev)) then 
    ilev_ = ilev
  else
    ilev_ = 1 
  end if

  if (.not.allocated(p%baseprecv)) then 
    info = 3111
    return 
  endif
  nlev_ = size(p%baseprecv)

  if ((ilev_<1).or.(ilev_ > nlev_)) then 
    write(0,*) 'PRECSET ERRROR: ilev out of bounds'
    info = -1
    return
  endif
  if (.not.allocated(p%baseprecv(ilev_)%dprcparm)) then 
    info = 3111
    return 
  endif

  if (ilev_ == 1) then 
    ! Rules for fine level are slightly different. 
    select case(what) 
    case(mld_fact_thrs_)
      p%baseprecv(ilev_)%dprcparm(what)  = val
    case default
      write(0,*) 'Error: trying to call PRECSET with an invalid WHAT'
      info = -2
    end select
  else if (ilev_ > 1) then 
    select case(what) 
    case(mld_aggr_damp_,mld_fact_thrs_)
      p%baseprecv(ilev_)%dprcparm(what)  = val
    case default
      write(0,*) 'Error: trying to call PRECSET with an invalid WHAT'
      info = -2
    end select
  endif

end subroutine mld_dprecsetd
