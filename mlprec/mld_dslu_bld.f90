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
subroutine mld_dslu_bld(a,desc_a,p,info)
  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_dslu_bld

  implicit none 

  type(psb_dspmat_type), intent(inout)   :: a
  type(psb_desc_type), intent(in)        :: desc_a
  type(mld_dbaseprc_type), intent(inout) :: p
  integer, intent(out)                   :: info

  integer                  :: i,j,nza,nzb,nzt,ictxt,me,np,err_act
  logical, parameter :: debug=.false.
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  name='mld_dslu_bld'
  call psb_erractionsave(err_act)

  ictxt = psb_cd_get_context(desc_a)

  call psb_info(ictxt, me, np)

  if (toupper(a%fida) /= 'CSR') then 
    write(0,*) 'Unimplemented input to mld_slu_BLD'
    goto 9999
  endif


  nzt = psb_sp_get_nnzeros(a)

  if (Debug) then 
    write(0,*) me,'Calling mld_slu_factor ',nzt,a%m,&
         & a%k,p%desc_data%matrix_data(psb_n_row_)
    call psb_barrier(ictxt)
  endif

  call mld_dslu_factor(a%m,nzt,&
       & a%aspk,a%ia2,a%ia1,p%iprcparm(mld_slu_ptr_),info)

  if (info /= 0) then
    ch_err='mld_slu_fact'
    call psb_errpush(4110,name,a_err=ch_err,i_err=(/info,0,0,0,0/))
    goto 9999
  end if

  if (Debug) then 
    write(0,*) me, 'SPLUBLD: Done mld_slu_Factor',info,p%iprcparm(mld_slu_ptr_)
    call psb_barrier(ictxt)
  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_dslu_bld

