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
! File: mld_cslu_bld.f90
!
! Subroutine: mld_cslu_bld
! Version:    complex
!
!  This routine computes the LU factorization of the local part of the matrix
!  stored into a, by using SuperLU. 
!  
!  The local matrix to be factorized is
!  - either a submatrix of the distributed matrix corresponding to any level
!    of a multilevel preconditioner, and its factorization is used to build
!    the 'base preconditioner' corresponding to that level,
!  - or a copy of the whole matrix corresponding to the coarsest level of
!    a multilevel preconditioner, and its factorization is used to build
!    the 'base preconditioner' corresponding to the coarsest level.
!
!  The data structure allocated by SuperLU to store the L and U factors is
!  pointed by p%iprcparm(mld_slu_ptr_).
!
!
! Arguments:
!    a       -  type(psb_zspmat_type), input/output.
!               The sparse matrix structure containing the local submatrix to
!               be factorized.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor associated to a.
!    p       -  type(mld_cbaseprec_type), input/output.
!               The 'base preconditioner' data structure containing the pointer, 
!               p%iprcparm(mld_slu_ptr_), to the data structure used by SuperLU
!               to store the L and U factors.
!    info    -  integer, output.                                                             
!               Error code.
! 
subroutine mld_cslu_bld(a,desc_a,p,info)

  use psb_base_mod
  use mld_c_inner_mod, mld_protect_name => mld_cslu_bld

  implicit none 

  ! Arguments
  type(psb_cspmat_type), intent(inout)   :: a
  type(psb_desc_type), intent(in)        :: desc_a
  type(mld_cbaseprec_type), intent(inout) :: p
  integer, intent(out)                   :: info

  ! Local variables
  integer           :: ictxt,me,np,err_act
  character(len=20) :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=psb_success_
  name='mld_cslu_bld'
  call psb_erractionsave(err_act)

  ictxt = desc_a%get_context()

  call psb_info(ictxt, me, np)

  !
  ! Compute the LU factorization
  !
  select type(aa=>a%a)
  type is (psb_c_csr_sparse_mat) 
    call mld_cslu_fact(aa%get_nrows(),aa%get_nzeros(),&
       & aa%val,aa%ja,aa%irp,p%iprcparm(mld_slu_ptr_),info)

    if (info /= psb_success_) then
      ch_err='mld_slu_fact'
      call psb_errpush(4110,name,a_err=ch_err,i_err=(/info,0,0,0,0/))
      goto 9999
    end if

  class default 
    info=psb_err_unsupported_format_
    ch_err=a%get_fmt()
    call psb_errpush(info,name,a_err=ch_err)
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

end subroutine mld_cslu_bld

