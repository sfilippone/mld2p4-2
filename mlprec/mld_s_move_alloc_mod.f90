!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010
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
! File: mld_move_alloc_mod.f90
!
! Module: mld_move_alloc_mod
!
!  This module defines move_alloc-like routines, and related interfaces,
!  for the preconditioner data structures. .   
!

module mld_s_move_alloc_mod

  use mld_s_prec_type

  interface mld_move_alloc
    module procedure  mld_sonelev_prec_move_alloc,&
         & mld_sprec_move_alloc
  end interface

contains

  subroutine mld_sonelev_prec_move_alloc(a, b,info)
    use psb_sparse_mod
    implicit none
    type(mld_sonelev_type), intent(inout) :: a, b
    integer, intent(out) :: info 
    
    call mld_precfree(b,info)
    call move_alloc(a%sm,b%sm)
    if (info == psb_success_) call psb_move_alloc(a%ac,b%ac,info) 
    if (info == psb_success_) call psb_move_alloc(a%desc_ac,b%desc_ac,info) 
    if (info == psb_success_) call psb_move_alloc(a%map,b%map,info) 
    b%base_a    => a%base_a
    b%base_desc => a%base_desc
    
  end subroutine mld_sonelev_prec_move_alloc

  subroutine mld_sprec_move_alloc(a, b,info)
    use psb_sparse_mod
    implicit none
    type(mld_sprec_type), intent(inout) :: a
    type(mld_sprec_type), intent(inout), target :: b
    integer, intent(out) :: info 
    integer :: i,isz
    
    if (allocated(b%precv)) then 
      ! This might not be required if FINAL procedures are available.
      call mld_precfree(b,info)
      if (info /= psb_success_) then 
        !       ?????
    !!$        return
      endif
    end if

    call move_alloc(a%precv,b%precv)
    ! Fix the pointers except on level 1.
    do i=2, isz
      b%precv(i)%base_a    => b%precv(i)%ac
      b%precv(i)%base_desc => b%precv(i)%desc_ac
      b%precv(i)%map%p_desc_X => b%precv(i-1)%base_desc
      b%precv(i)%map%p_desc_Y => b%precv(i)%base_desc
    end do
  end subroutine mld_sprec_move_alloc


end module mld_s_move_alloc_mod
