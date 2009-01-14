!!$
!!$ 
!!$                           MLD2P4  version 1.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 2.2)
!!$  
!!$  (C) Copyright 2008
!!$
!!$                      Salvatore Filippone  University of Rome Tor Vergata       
!!$                      Alfredo Buttari      University of Rome Tor Vergata
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
! File: mld_transfer_mod.f90
!
! Package: mld_prec_type
!          Data structure(s) for sparse matrices
!
!  This module defines: 
!  - the mld_prec_type data structure containing the preconditioner;
!  - character constants describing the preconditioner, used by the routine
!    printing out a preconditioner description;
!  - the interfaces to the routines for the management of the preconditioner
!    data structure (see below).
!
!  It contains:
!  - mld_dprec_sizeof, mld_dbaseprc_sizeof, mld_out_prec_descr, mld_file_prec_descr,
!    mld_icheck_def, mld_dcheck_def, mld_dbase_precfree, mld_nullify_dbaseprec,
!    is_legal_..._..., and their complex versions (if applicable).
!  These routines check if the preconditioner is correctly defined,      print a
!  description of the preconditioner, and deallocate its data structure.  
!

module mld_transfer_mod

  use mld_prec_type

  interface mld_transfer
    module procedure  mld_sbaseprec_transfer, mld_sonelev_prec_transfer,&
         & mld_sprec_transfer,&
         & mld_dbaseprec_transfer, mld_donelev_prec_transfer,&
         & mld_dprec_transfer,&
         & mld_cbaseprec_transfer, mld_conelev_prec_transfer,&
         & mld_cprec_transfer,&
         & mld_zbaseprec_transfer, mld_zonelev_prec_transfer,&
         & mld_zprec_transfer
  end interface

contains


  subroutine mld_sbaseprec_transfer(a, b,info)
    use psb_base_mod
    implicit none
    type(mld_sbaseprc_type), intent(inout) :: a, b
    integer, intent(out) :: info 
    integer :: i, isz
    
    call mld_precfree(b,info)
    if (info == 0) call psb_transfer(a%iprcparm,b%iprcparm,info) 
    if (info == 0) call psb_transfer(a%rprcparm,b%rprcparm,info) 
    if (info == 0) call psb_transfer(a%desc_data,b%desc_data,info) 
    if (info == 0) call psb_transfer(a%perm,b%perm,info) 
    if (info == 0) call psb_transfer(a%invperm,b%invperm,info) 
#ifdef HAVE_MOVE_ALLOC
    call move_alloc(a%av,b%av)
#else
    if (allocated(a%av)) then 
      isz = size(a%av)
      allocate(b%av(isz),stat=info) 
      do i=1,isz
        if (info == 0) call psb_transfer(a%av(i), b%av(i), info)
      end do
      
      if (info == 0) deallocate(a%av,stat=info)
    end if
#endif
    if (info /= 0) then
      write(0,*) 'Error in baseprec_:transfer',info
    end if

  end subroutine mld_sbaseprec_transfer

  subroutine mld_sonelev_prec_transfer(a, b,info)
    use psb_base_mod
    implicit none
    type(mld_s_interlev_prec_type), intent(inout) :: a, b
    integer, intent(out) :: info 
    
    call mld_precfree(b,info)
    if (info == 0) call mld_transfer(a%prec,b%prec,info) 
    if (info == 0) call psb_transfer(a%iprcparm,b%iprcparm,info) 
    if (info == 0) call psb_transfer(a%rprcparm,b%rprcparm,info) 
    if (info == 0) call psb_transfer(a%ac,b%ac,info) 
    if (info == 0) call psb_transfer(a%desc_ac,b%desc_ac,info) 
    if (info == 0) call psb_transfer(a%mlia,b%mlia,info) 
    if (info == 0) call psb_transfer(a%nlaggr,b%nlaggr,info) 
    if (info == 0) call psb_transfer(a%map_desc,b%map_desc,info) 
    b%base_a    => a%base_a
    b%base_desc => a%base_desc
    
  end subroutine mld_sonelev_prec_transfer

  subroutine mld_sprec_transfer(a, b,info)
    use psb_base_mod
    implicit none
    type(mld_sprec_type), intent(inout) :: a
    type(mld_sprec_type), intent(inout), target :: b
    integer, intent(out) :: info 
    integer :: i,isz
    
    if (allocated(b%precv)) then 
      ! This might not be required if FINAL procedures are available.
      call mld_precfree(b,info)
      if (info /= 0) then 
        !       ?????
    !!$        return
      endif
    end if

#ifdef HAVE_MOVE_ALLOC
    call move_alloc(a%precv,b%precv)
#else 
    if (.not.allocated(a%precv)) return
    isz = size(a%precv)
    allocate(b%precv(isz),stat=info) 
    if (info /= 0) then 
      write(0,*) 'Memory allocation failure in prec_transfer'
      return
    end if
    do i=1,isz
      call mld_transfer(a%precv(i),b%precv(i),info)
    end do
    deallocate(a%precv,stat=info)
#endif
    ! Fix the pointers except on level 1.
    do i=2, isz
      b%precv(i)%base_a    => b%precv(i)%ac
      b%precv(i)%base_desc => b%precv(i)%desc_ac
      b%precv(i)%map_desc%p_desc_X => b%precv(i-1)%base_desc
      b%precv(i)%map_desc%p_desc_Y => b%precv(i)%base_desc
    end do
  end subroutine mld_sprec_transfer


  subroutine mld_dbaseprec_transfer(a, b,info)
    use psb_base_mod
    implicit none
    type(mld_dbaseprc_type), intent(inout) :: a, b
    integer, intent(out) :: info 
    integer :: i, isz
    
    call mld_precfree(b,info)
    if (info == 0) call psb_transfer(a%iprcparm,b%iprcparm,info) 
    if (info == 0) call psb_transfer(a%rprcparm,b%rprcparm,info) 
    if (info == 0) call psb_transfer(a%desc_data,b%desc_data,info) 
    if (info == 0) call psb_transfer(a%perm,b%perm,info) 
    if (info == 0) call psb_transfer(a%invperm,b%invperm,info) 
#ifdef HAVE_MOVE_ALLOC
    call move_alloc(a%av,b%av)
#else
    if (allocated(a%av)) then 
      isz = size(a%av)
      allocate(b%av(isz),stat=info) 
      do i=1,isz
        if (info == 0) call psb_transfer(a%av(i), b%av(i), info)
      end do
      
      if (info == 0) deallocate(a%av,stat=info)
    end if
#endif
    if (info /= 0) then
      write(0,*) 'Error in baseprec_:transfer',info
    end if

  end subroutine mld_dbaseprec_transfer

  subroutine mld_donelev_prec_transfer(a, b,info)
    use psb_base_mod
    implicit none
    type(mld_d_interlev_prec_type), intent(inout) :: a, b
    integer, intent(out) :: info 
    
    call mld_precfree(b,info)
    if (info == 0) call mld_transfer(a%prec,b%prec,info) 
    if (info == 0) call psb_transfer(a%iprcparm,b%iprcparm,info) 
    if (info == 0) call psb_transfer(a%rprcparm,b%rprcparm,info) 
    if (info == 0) call psb_transfer(a%ac,b%ac,info) 
    if (info == 0) call psb_transfer(a%desc_ac,b%desc_ac,info) 
    if (info == 0) call psb_transfer(a%mlia,b%mlia,info) 
    if (info == 0) call psb_transfer(a%nlaggr,b%nlaggr,info) 
    if (info == 0) call psb_transfer(a%map_desc,b%map_desc,info) 
    b%base_a    => a%base_a
    b%base_desc => a%base_desc
    
  end subroutine mld_donelev_prec_transfer

  subroutine mld_dprec_transfer(a, b,info)
    use psb_base_mod
    implicit none
    type(mld_dprec_type), intent(inout) :: a
    type(mld_dprec_type), intent(inout), target :: b
    integer, intent(out) :: info 
    integer :: i,isz
    
    if (allocated(b%precv)) then 
      ! This might not be required if FINAL procedures are available.
      call mld_precfree(b,info)
      if (info /= 0) then 
        !       ?????
    !!$        return
      endif
    end if

#ifdef HAVE_MOVE_ALLOC
    call move_alloc(a%precv,b%precv)
#else 
    if (.not.allocated(a%precv)) return
    isz = size(a%precv)
    allocate(b%precv(isz),stat=info) 
    if (info /= 0) then 
      write(0,*) 'Memory allocation failure in prec_transfer'
      return
    end if
    do i=1,isz
      call mld_transfer(a%precv(i),b%precv(i),info)
    end do
    deallocate(a%precv,stat=info)
#endif
    ! Fix the pointers except on level 1.
    do i=2, isz
      b%precv(i)%base_a    => b%precv(i)%ac
      b%precv(i)%base_desc => b%precv(i)%desc_ac
      b%precv(i)%map_desc%p_desc_X => b%precv(i-1)%base_desc
      b%precv(i)%map_desc%p_desc_Y => b%precv(i)%base_desc
    end do
  end subroutine mld_dprec_transfer


  subroutine mld_cbaseprec_transfer(a, b,info)
    use psb_base_mod
    implicit none
    type(mld_cbaseprc_type), intent(inout) :: a, b
    integer, intent(out) :: info 
    integer :: i, isz
    
    call mld_precfree(b,info)
    if (info == 0) call psb_transfer(a%iprcparm,b%iprcparm,info) 
    if (info == 0) call psb_transfer(a%rprcparm,b%rprcparm,info) 
    if (info == 0) call psb_transfer(a%desc_data,b%desc_data,info) 
    if (info == 0) call psb_transfer(a%perm,b%perm,info) 
    if (info == 0) call psb_transfer(a%invperm,b%invperm,info) 
#ifdef HAVE_MOVE_ALLOC
    call move_alloc(a%av,b%av)
#else
    if (allocated(a%av)) then 
      isz = size(a%av)
      allocate(b%av(isz),stat=info) 
      do i=1,isz
        if (info == 0) call psb_transfer(a%av(i), b%av(i), info)
      end do
      
      if (info == 0) deallocate(a%av,stat=info)
    end if
#endif
    if (info /= 0) then
      write(0,*) 'Error in baseprec_:transfer',info
    end if

  end subroutine mld_cbaseprec_transfer

  subroutine mld_conelev_prec_transfer(a, b,info)
    use psb_base_mod
    implicit none
    type(mld_c_interlev_prec_type), intent(inout) :: a, b
    integer, intent(out) :: info 
    
    call mld_precfree(b,info)
    if (info == 0) call mld_transfer(a%prec,b%prec,info) 
    if (info == 0) call psb_transfer(a%iprcparm,b%iprcparm,info) 
    if (info == 0) call psb_transfer(a%rprcparm,b%rprcparm,info) 
    if (info == 0) call psb_transfer(a%ac,b%ac,info) 
    if (info == 0) call psb_transfer(a%desc_ac,b%desc_ac,info) 
    if (info == 0) call psb_transfer(a%mlia,b%mlia,info) 
    if (info == 0) call psb_transfer(a%nlaggr,b%nlaggr,info) 
    if (info == 0) call psb_transfer(a%map_desc,b%map_desc,info) 
    b%base_a    => a%base_a
    b%base_desc => a%base_desc
    
  end subroutine mld_conelev_prec_transfer

  subroutine mld_cprec_transfer(a, b,info)
    use psb_base_mod
    implicit none
    type(mld_cprec_type), intent(inout) :: a
    type(mld_cprec_type), intent(inout), target :: b
    integer, intent(out) :: info 
    integer :: i,isz
    
    if (allocated(b%precv)) then 
      ! This might not be required if FINAL procedures are available.
      call mld_precfree(b,info)
      if (info /= 0) then 
        !       ?????
    !!$        return
      endif
    end if

#ifdef HAVE_MOVE_ALLOC
    call move_alloc(a%precv,b%precv)
#else 
    if (.not.allocated(a%precv)) return
    isz = size(a%precv)
    allocate(b%precv(isz),stat=info) 
    if (info /= 0) then 
      write(0,*) 'Memory allocation failure in prec_transfer'
      return
    end if
    do i=1,isz
      call mld_transfer(a%precv(i),b%precv(i),info)
    end do
    deallocate(a%precv,stat=info)
#endif
    ! Fix the pointers except on level 1.
    do i=2, isz
      b%precv(i)%base_a    => b%precv(i)%ac
      b%precv(i)%base_desc => b%precv(i)%desc_ac
      b%precv(i)%map_desc%p_desc_X => b%precv(i-1)%base_desc
      b%precv(i)%map_desc%p_desc_Y => b%precv(i)%base_desc
    end do
  end subroutine mld_cprec_transfer


  subroutine mld_zbaseprec_transfer(a, b,info)
    use psb_base_mod
    implicit none
    type(mld_zbaseprc_type), intent(inout) :: a, b
    integer, intent(out) :: info 
    integer :: i, isz
    
    call mld_precfree(b,info)
    if (info == 0) call psb_transfer(a%iprcparm,b%iprcparm,info) 
    if (info == 0) call psb_transfer(a%rprcparm,b%rprcparm,info) 
    if (info == 0) call psb_transfer(a%desc_data,b%desc_data,info) 
    if (info == 0) call psb_transfer(a%perm,b%perm,info) 
    if (info == 0) call psb_transfer(a%invperm,b%invperm,info) 
#ifdef HAVE_MOVE_ALLOC
    call move_alloc(a%av,b%av)
#else
    if (allocated(a%av)) then 
      isz = size(a%av)
      allocate(b%av(isz),stat=info) 
      do i=1,isz
        if (info == 0) call psb_transfer(a%av(i), b%av(i), info)
      end do
      
      if (info == 0) deallocate(a%av,stat=info)
    end if
#endif
    if (info /= 0) then
      write(0,*) 'Error in baseprec_:transfer',info
    end if

  end subroutine mld_zbaseprec_transfer

  subroutine mld_zonelev_prec_transfer(a, b,info)
    use psb_base_mod
    implicit none
    type(mld_z_interlev_prec_type), intent(inout) :: a, b
    integer, intent(out) :: info 
    
    call mld_precfree(b,info)
    if (info == 0) call mld_transfer(a%prec,b%prec,info) 
    if (info == 0) call psb_transfer(a%iprcparm,b%iprcparm,info) 
    if (info == 0) call psb_transfer(a%rprcparm,b%rprcparm,info) 
    if (info == 0) call psb_transfer(a%ac,b%ac,info) 
    if (info == 0) call psb_transfer(a%desc_ac,b%desc_ac,info) 
    if (info == 0) call psb_transfer(a%mlia,b%mlia,info) 
    if (info == 0) call psb_transfer(a%nlaggr,b%nlaggr,info) 
    if (info == 0) call psb_transfer(a%map_desc,b%map_desc,info) 
    b%base_a    => a%base_a
    b%base_desc => a%base_desc
    
  end subroutine mld_zonelev_prec_transfer

  subroutine mld_zprec_transfer(a, b,info)
    use psb_base_mod
    implicit none
    type(mld_zprec_type), intent(inout) :: a
    type(mld_zprec_type), intent(inout), target :: b
    integer, intent(out) :: info 
    integer :: i,isz
    
    if (allocated(b%precv)) then 
      ! This might not be required if FINAL procedures are available.
      call mld_precfree(b,info)
      if (info /= 0) then 
        !       ?????
    !!$        return
      endif
    end if

#ifdef HAVE_MOVE_ALLOC
    call move_alloc(a%precv,b%precv)
#else 
    if (.not.allocated(a%precv)) return
    isz = size(a%precv)
    allocate(b%precv(isz),stat=info) 
    if (info /= 0) then 
      write(0,*) 'Memory allocation failure in prec_transfer'
      return
    end if
    do i=1,isz
      call mld_transfer(a%precv(i),b%precv(i),info)
    end do
    deallocate(a%precv,stat=info)
#endif
    ! Fix the pointers except on level 1.
    do i=2, isz
      b%precv(i)%base_a    => b%precv(i)%ac
      b%precv(i)%base_desc => b%precv(i)%desc_ac
      b%precv(i)%map_desc%p_desc_X => b%precv(i-1)%base_desc
      b%precv(i)%map_desc%p_desc_Y => b%precv(i)%base_desc
    end do
  end subroutine mld_zprec_transfer


end module mld_transfer_mod
