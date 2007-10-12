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
subroutine mld_zaggrmat_asb(a,desc_a,ac,desc_ac,p,info)
  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_zaggrmat_asb

  implicit none

  type(psb_zspmat_type), intent(in), target  :: a
  type(mld_zbaseprc_type), intent(inout),target     :: p
  type(psb_zspmat_type), intent(inout), target :: ac
  type(psb_desc_type), intent(in)            :: desc_a
  type(psb_desc_type), intent(inout)         :: desc_ac
  integer, intent(out)                       :: info

  logical, parameter :: aggr_dump=.false.
  integer ::ictxt,np,me, err_act,icomm
  character(len=20) :: name, ch_err

  name='mld_aggrmat_asb'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  ictxt = psb_cd_get_context(desc_a)
  icomm = psb_cd_get_mpic(desc_a)

  call psb_info(ictxt, me, np)

  select case (p%iprcparm(mld_aggr_kind_))
  case (mld_no_smooth_) 

    call mld_aggrmat_raw_asb(a,desc_a,ac,desc_ac,p,info)

    if(info /= 0) then
      call psb_errpush(4010,name,a_err='raw_aggregate')
      goto 9999
    end if
    if (aggr_dump) call psb_csprt(90+me,ac,head='% Raw aggregate.')

  case(mld_smooth_prol_,mld_biz_prol_) 
    if (aggr_dump) call psb_csprt(70+me,a,head='% Input matrix')
    call mld_aggrmat_smth_asb(a,desc_a,ac,desc_ac,p,info)

    if(info /= 0) then
      call psb_errpush(4010,name,a_err='smooth_aggregate')
      goto 9999
    end if
    if (aggr_dump) call psb_csprt(90+me,ac,head='% Smooth aggregate.')
  case default
    call psb_errpush(4010,name,a_err=name)
    goto 9999

  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_zaggrmat_asb
