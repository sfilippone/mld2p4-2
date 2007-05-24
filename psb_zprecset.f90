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
subroutine psb_zprecseti(p,what,val,info,ilev)

  use psb_base_mod
  use psb_prec_mod, psb_protect_name => psb_zprecseti

  implicit none
  type(psb_zprec_type), intent(inout)    :: p
  integer, intent(in)                    :: what 
  integer, intent(in)                    :: val
  integer, intent(out)                   :: info
  integer, optional, intent(in)          :: ilev
  integer                                :: ilev_, nlev_

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
  if (.not.allocated(p%baseprecv(ilev_)%iprcparm)) then 
    write(0,*) 'Error: trying to call PRECSET on an uninitialized preconditioner'
    info = -1 
    return 
  endif



  if (ilev_ == 1) then 
    ! Rules for fine level are slightly different. 
    select case(what) 
    case(p_type_,f_type_,restr_,prol_,iren_,n_ovr_,ilu_fill_in_,jac_sweeps_)
      p%baseprecv(ilev_)%iprcparm(what)  = val
    case default
      write(0,*) 'Error: trying to call PRECSET with an invalid WHAT'
      info = -2
    end select
  else if (ilev_ > 1) then 
    select case(what) 
    case(p_type_,f_type_,restr_,prol_,iren_,n_ovr_,ilu_fill_in_,jac_sweeps_,&
         & ml_type_,aggr_alg_,smth_kind_,coarse_mat_,smth_pos_,om_choice_)
      p%baseprecv(ilev_)%iprcparm(what)  = val
    case default
      write(0,*) 'Error: trying to call PRECSET with an invalid WHAT'
      info = -2
    end select
  endif

end subroutine psb_zprecseti
subroutine psb_zprecsetd(p,what,val,info,ilev)

  use psb_base_mod
  use psb_prec_mod, psb_protect_name => psb_zprecsetd

  implicit none
  type(psb_zprec_type), intent(inout)    :: p
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
      ! we implement F_ILU_E_
    case default
      write(0,*) 'Error: trying to call PRECSET with an invalid WHAT'
      info = -2
    end select
  else if (ilev_ > 1) then 
    select case(what) 
    case(smooth_omega_)
      p%baseprecv(ilev_)%dprcparm(what)  = val
    case default
      write(0,*) 'Error: trying to call PRECSET with an invalid WHAT'
      info = -2
    end select
  endif

end subroutine psb_zprecsetd
