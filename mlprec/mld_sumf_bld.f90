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
! File: mld_sumf_bld.f90
!
! Subroutine: mld_sumf_bld
! Version:    real
!
!  This routine computes the LU factorization of the local part of the matrix
!  stored into a, by using UMFPACK. 
!  
!  The local matrix to be factorized is
!  - either a submatrix of the distributed matrix corresponding to any level
!    of a multilevel preconditioner, and its factorization is used to build
!    the 'base preconditioner' corresponding to that level,
!  - or a copy of the whole matrix corresponding to the coarsest level of
!    a multilevel preconditioner, and its factorization is used to build
!    the 'base preconditioner' corresponding to the coarsest level.
!
!  The data structures allocated by UMFPACK to compute the symbolic and the
!  numeric factorization are pointed by p%iprcparm(mld_umf_symptr_) and
!  p%iprcparm(mld_umf_numptr_).
!
!
! Arguments:
!    a       -  type(psb_sspmat_type), input/output.
!               The sparse matrix structure containing the local submatrix
!               to be factorized. Note that a is intent(inout), and not only
!               intent(in), since the row and column indices of the      matrix
!               stored in a are shifted by -1, and then again by +1, by the
!               routine mld_sumf_fact, which is an interface to the UMFPACK
!               C code performing the factorization.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor associated to a.
!    p       -  type(mld_sbaseprec_type), input/output.
!               The 'base preconditioner' data structure containing the pointers,
!               p%iprcparm(mld_umf_symptr_) and p%iprcparm(mld_umf_numptr_),
!               to the data structures used by UMFPACK for computing the LU
!               factorization.
!    info    -  integer, output.                                                
!               Error code.
!
subroutine mld_sumf_bld(a,desc_a,p,info)

  use psb_sparse_mod
  use mld_s_inner_mod, mld_protect_name => mld_sumf_bld

  implicit none 

! Arguments
  type(psb_sspmat_type), intent(inout)      :: a
  type(psb_desc_type), intent(in)        :: desc_a
  type(mld_sbaseprec_type), intent(inout) :: p
  integer, intent(out)                   :: info

  ! Local variables
  integer            :: ictxt,me,np,err_act
  character(len=20)  :: name, ch_err

  info=psb_success_
  name='mld_sumf_bld'
  call psb_erractionsave(err_act)
  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)

  !
  ! Compute the LU factorization
  !
  select type(aa=>a%a)
!!$  type is (psb_s_csc_sparse_mat)
!!$    call mld_sumf_fact(aa%m,aa%get_nzeros(),&
!!$         & aa%val,aa%ia,aa%icp,&
!!$         & p%iprcparm(mld_umf_symptr_),p%iprcparm(mld_umf_numptr_),info)
!!$    
!!$    if (info /= psb_success_) then
!!$      info=4110
!!$      call psb_errpush(info,name,a_err='mld_umf_fact',i_err=(/info,0,0,0,0/))
!!$      goto 9999
!!$    end if
  class default 
    info=psb_err_unsupported_format_
    ch_err = aa%get_fmt()
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

end subroutine mld_sumf_bld



