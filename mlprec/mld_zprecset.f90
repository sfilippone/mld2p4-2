!!$ 
!!$ 
!!$                    MD2P4
!!$    Multilevel Domain Decomposition Parallel Preconditioner Package for PSBLAS
!!$                      for 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$                       Daniela di Serafino    Second University of Naples
!!$                       Pasqua D'Ambra         ICAR-CNR                      
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the MD2P4 group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MD2P4 GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
subroutine mld_zprecseti(p,what,val,info,ilev)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_zprecseti

  implicit none

  type(mld_zprec_type), intent(inout)    :: p
  integer, intent(in)                    :: what 
  integer, intent(in)                    :: val
  integer, intent(out)                   :: info
  integer, optional, intent(in)          :: ilev
  integer                                :: ilev_, nlev_

  info = 0

  if (.not.allocated(p%baseprecv)) then 
    write(0,*) 'Error: trying to call PRECSET on an uninitialized preconditioner'
    info = -1 
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
    write(0,*) 'Error: trying to call PRECSET on an uninitialized preconditioner'
    info = -1 
    return 
  endif


  if (present(ilev)) then 
    
    if (ilev_ == 1) then 
      ! Rules for fine level are slightly different. 
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

end subroutine mld_zprecseti
subroutine mld_zprecsetd(p,what,val,info,ilev)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_zprecsetd

  implicit none
  type(mld_zprec_type), intent(inout)    :: p
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
    write(0,*) 'Error: trying to call PRECSET on an uninitialized preconditioner'
    info = -1 
    return 
  endif
  nlev_ = size(p%baseprecv)

  if ((ilev_<1).or.(ilev_ > nlev_)) then 
    write(0,*) 'PRECSET ERRROR: ilev out of bounds'
    info = -1
    return
  endif
  if (.not.allocated(p%baseprecv(ilev_)%dprcparm)) then 
    write(0,*) 'Error: trying to call PRECSET on an uninitialized preconditioner'
    info = -1 
    return 
  endif

  if (ilev_ == 1) then 
    ! Rules for fine level are slightly different. 
    select case(what) 
      ! Right now we don't have any at base level. Will  change when
      ! we implement mld_ilu_t_
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

end subroutine mld_zprecsetd
