!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.3)
!!$  
!!$  (C) Copyright 2008, 2010, 2012, 2015
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
subroutine mld_s_base_onelev_setsv(lev,val,info,pos)

  use psb_base_mod
  use mld_s_prec_mod, mld_protect_name => mld_s_base_onelev_setsv

  implicit none

  ! Arguments
  class(mld_s_onelev_type), target, intent(inout) :: lev
  class(mld_s_base_solver_type), intent(in)       :: val
  integer(psb_ipk_), intent(out)                  :: info
  character(len=*), optional, intent(in)          :: pos
  
  ! Local variables
  integer(psb_ipk_)                :: ipos_
  character(len=*), parameter      :: name='mld_base_onelev_setsv'

  info = psb_success_

  if (present(pos)) then
    select case(psb_toupper(trim(pos)))
    case('PRE')
      ipos_ = mld_pre_smooth_
    case('POST')
      ipos_ = mld_post_smooth_
    case default
      ipos_ = mld_pre_smooth_
    end select
  else
    ipos_ = mld_pre_smooth_
  end if
  
  select case(ipos_)
  case(mld_pre_smooth_) 
    if (allocated(lev%sm)) then 
      if (allocated(lev%sm%sv)) then
        if (.not.same_type_as(lev%sm%sv,val))  then
          call lev%sm%sv%free(info)
          deallocate(lev%sm%sv,stat=info)
          if (info /= 0) then
            info = 3111
            return
          end if
        end if
      end if
      
      if (.not.allocated(lev%sm%sv)) then 
#ifdef HAVE_MOLD 
        allocate(lev%sm%sv,mold=val,stat=info) 
#else
        allocate(lev%sm%sv,source=val,stat=info) 
#endif
        if (info /= 0) then
          info = 3111
          return
        end if
      end if
      call lev%sm%sv%default()
    else
      info = 3111
      write(psb_err_unit,*) name,&
           &': Error: uninitialized preconditioner component,',&
           &' should call MLD_PRECINIT/MLD_PRECSET' 
      return 
      
    end if
      
  case(mld_post_smooth_) 

    if (allocated(lev%sm2a)) then 
      if (allocated(lev%sm2a%sv)) then
        if (.not.same_type_as(lev%sm2a%sv,val))  then
          call lev%sm2a%sv%free(info)
          deallocate(lev%sm2a%sv,stat=info)
          if (info /= 0) then
            info = 3111
            return
          end if
        end if
      end if
      if (.not.allocated(lev%sm2a%sv)) then 
#ifdef HAVE_MOLD 
        allocate(lev%sm2a%sv,mold=val,stat=info) 
#else
        allocate(lev%sm2a%sv,source=val,stat=info) 
#endif
        if (info /= 0) then
          info = 3111
          return
        end if
      end if
      call lev%sm2a%sv%default()
      
    else
      info = 3111
      write(psb_err_unit,*) name,&
           &': Error: uninitialized preconditioner component,',&
           &' should call MLD_PRECINIT/MLD_PRECSET' 
      return 
      
    end if
    
  end select
  
end subroutine mld_s_base_onelev_setsv

