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
subroutine mld_zprec_aply(prec,x,y,desc_data,info,trans, work)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_zprec_aply

  implicit none

  type(psb_desc_type),intent(in)      :: desc_data
  type(mld_zprec_type), intent(in)    :: prec
  complex(kind(0.d0)),intent(inout)   :: x(:), y(:)
  integer, intent(out)                :: info
  character(len=1), optional          :: trans
  complex(kind(0.d0)), optional, target  :: work(:)

  ! Local variables
  character     :: trans_ 
  complex(kind(1.d0)), pointer :: work_(:)
  integer :: ictxt,np,me,err_act,iwsz
  logical,parameter                 :: debug=.false., debugprt=.false.
  character(len=20)   :: name
  
  name='mld_zprec_aply'
  info = 0
  call psb_erractionsave(err_act)

  ictxt = psb_cd_get_context(desc_data)
  call psb_info(ictxt, me, np)

  if (present(trans)) then 
    trans_=trans
  else
    trans_='N'
  end if

  if (present(work)) then 
    work_ => work
  else
    iwsz = max(1,4*psb_cd_get_local_cols(desc_data))
    allocate(work_(iwsz),stat=info)
    if (info /= 0) then 
      info=4025
      call psb_errpush(info,name,i_err=(/iwsz,0,0,0,0/),&
           & a_err='complex(kind(1.d0))')
      goto 9999      
    end if

  end if

  if (.not.(allocated(prec%baseprecv))) then 
    write(0,*) 'Inconsistent preconditioner: neither SMTH nor BASE?'      
  end if
  if (size(prec%baseprecv) >1) then 
    if (debug) write(0,*) 'Into mlprc_aply',size(x),size(y)
    call mld_mlprec_aply(zone,prec%baseprecv,x,zzero,y,desc_data,trans_,work_,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='mld_zmlprec_aply')
      goto 9999
    end if

  else  if (size(prec%baseprecv) == 1) then 
    call mld_baseprec_aply(zone,prec%baseprecv(1),x,zzero,y,desc_data,trans_, work_,info)
  else 
    write(0,*) 'Inconsistent preconditioner: size of baseprecv???' 
  endif

  if (present(work)) then 
  else
    deallocate(work_)
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

end subroutine mld_zprec_aply


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
subroutine mld_zprec_aply1(prec,x,desc_data,info,trans)
  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_zprec_aply1

  implicit none

  type(psb_desc_type),intent(in)    :: desc_data
  type(mld_zprec_type), intent(in)  :: prec
  complex(kind(0.d0)),intent(inout) :: x(:)
  integer, intent(out)              :: info
  character(len=1), optional        :: trans
  logical,parameter                 :: debug=.false., debugprt=.false.

  ! Local variables
  character     :: trans_
  integer :: ictxt,np,me,i, isz, err_act, int_err(5)
  complex(kind(1.d0)), pointer :: WW(:), w1(:)
  character(len=20)   :: name, ch_err
  name='mld_zprec_aply1'
  info = 0
  call psb_erractionsave(err_act)
  

  ictxt = psb_cd_get_context(desc_data)
  call psb_info(ictxt, me, np)
  if (present(trans)) then 
    trans_=trans
  else
    trans_='N'
  end if

  allocate(ww(size(x)),w1(size(x)),stat=info)
  if (info /= 0) then 
    info=4025
    call psb_errpush(info,name,i_err=(/2*size(x),0,0,0,0/),&
         & a_err='complex(kind(1.d0))')
    goto 9999      
  end if
  if (debug) write(0,*) 'Prc_aply1 Size(x) ',size(x), size(ww),size(w1)
  call mld_zprec_aply(prec,x,ww,desc_data,info,trans_,work=w1)
  if(info /=0) goto 9999
  x(:) = ww(:)
  deallocate(ww,W1)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name)
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
     call psb_error()
     return
  end if
  return
end subroutine mld_zprec_aply1
