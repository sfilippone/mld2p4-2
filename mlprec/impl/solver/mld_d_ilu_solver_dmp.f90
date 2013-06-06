!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010,2012,2013
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
subroutine mld_d_ilu_solver_dmp(sv,ictxt,level,info,prefix,head,solver)
  
  use psb_base_mod
  use mld_d_ilu_solver, mld_protect_name => mld_d_ilu_solver_dmp
  implicit none 
  class(mld_d_ilu_solver_type), intent(in) :: sv
  integer(psb_ipk_), intent(in)              :: ictxt,level
  integer(psb_ipk_), intent(out)             :: info
  character(len=*), intent(in), optional     :: prefix, head
  logical, optional, intent(in)              :: solver
  integer(psb_ipk_)  :: i, j, il1, iln, lname, lev
  integer(psb_ipk_)  :: icontxt,iam, np
  character(len=80)  :: prefix_
  character(len=120) :: fname ! len should be at least 20 more than
  logical :: solver_
  !  len of prefix_ 

  info = 0


  call psb_info(ictxt,iam,np)

  if (present(solver)) then 
    solver_ = solver
  else
    solver_ = .false. 
  end if

  if (solver_) then 
    if (present(prefix)) then 
      prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
    else
      prefix_ = "dump_slv_d"
    end if
    lname = len_trim(prefix_)
    fname = trim(prefix_)
    write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
    lname = lname + 5

    write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_lower.mtx'
    if (sv%l%is_asb()) &
         & call sv%l%print(fname,head=head)
    write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_diag.mtx'
    if (allocated(sv%d)) &
         & call psb_geprt(fname,sv%d,head=head)
    write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_upper.mtx'
    if (sv%u%is_asb()) &
         & call sv%u%print(fname,head=head)

  end if

end subroutine mld_d_ilu_solver_dmp
