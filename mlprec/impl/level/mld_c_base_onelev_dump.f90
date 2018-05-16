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
subroutine mld_c_base_onelev_dump(lv,level,info,prefix,head,ac,rp,&
     & smoother,solver,global_num)

  use psb_base_mod
  use mld_c_onelev_mod, mld_protect_name => mld_c_base_onelev_dump
  implicit none 
  class(mld_c_onelev_type), intent(in) :: lv
  integer(psb_ipk_), intent(in)          :: level
  integer(psb_ipk_), intent(out)         :: info
  character(len=*), intent(in), optional :: prefix, head
  logical, optional, intent(in)    :: ac, rp, smoother, solver,global_num
  ! Local variables
  integer(psb_ipk_) :: i, j, il1, iln, lname, lev
  integer(psb_ipk_) :: icontxt,iam, np
  character(len=80)  :: prefix_
  character(len=120) :: fname ! len should be at least 20 more than
  logical :: ac_, rp_,global_num_
  integer(psb_ipk_), allocatable :: ivr(:), ivc(:)

  info = 0

  if (present(prefix)) then 
    prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
  else
    prefix_ = "dump_lev_c"
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
  if (present(global_num)) then 
    global_num_ = global_num
  else
    global_num_ = .false. 
  end if
  lname = len_trim(prefix_)
  fname = trim(prefix_)
  write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
  lname = lname + 5

  if (global_num_) then
    if (level >= 2) then 
      if (ac_) then
        ivr = lv%desc_ac%get_global_indices(owned=.false.)
        write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_ac.mtx'
        call lv%ac%print(fname,head=head,iv=ivr)
      end if
      if (rp_) then 
        ivr = lv%map%p_desc_X%get_global_indices(owned=.false.)
        ivc = lv%map%p_desc_Y%get_global_indices(owned=.false.)
        write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_r.mtx'
        call lv%map%map_X2Y%print(fname,head=head,ivr=ivc,ivc=ivr)
        write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_p.mtx'
        call lv%map%map_Y2X%print(fname,head=head,ivr=ivr,ivc=ivc)
      end if
    end if
  else
    if (level >= 2) then 
      if (ac_) then 
        write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_ac.mtx'
        call lv%ac%print(fname,head=head)
      end if
      if (rp_) then 
        write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_r.mtx'
        call lv%map%map_X2Y%print(fname,head=head)
        write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_p.mtx'
        call lv%map%map_Y2X%print(fname,head=head)
      end if
    end if
  end if
  
  if (allocated(lv%sm)) then
    call lv%sm%dump(icontxt,level,info,smoother=smoother, &
         & solver=solver,prefix=trim(prefix_)//"_sm")
  end if
  if (allocated(lv%sm2a)) then
    call lv%sm2a%dump(icontxt,level,info,smoother=smoother, &
         & solver=solver,prefix=trim(prefix_)//"_sm2a")
  end if
  
end subroutine mld_c_base_onelev_dump
