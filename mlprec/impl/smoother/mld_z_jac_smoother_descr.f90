!  
!   
!                             MLD2P4  version 2.1
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!    
!    (C) Copyright 2008, 2010, 2012, 2015, 2017 
!  
!        Salvatore Filippone    Cranfield University, UK
!        Ambra Abdullahi Hassan University of Rome Tor Vergata, IT
!        Alfredo Buttari        CNRS-IRIT, Toulouse, FR
!        Pasqua D'Ambra         IAC-CNR, Naples, IT
!        Daniela di Serafino    University of Campania "L. Vanvitelli", Caserta, IT
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
subroutine mld_z_jac_smoother_descr(sm,info,iout,coarse)

  use psb_base_mod
  use mld_z_diag_solver
  use mld_z_jac_smoother, mld_protect_name => mld_z_jac_smoother_descr
  use mld_z_diag_solver
  use mld_z_gs_solver

  Implicit None

  ! Arguments
  class(mld_z_jac_smoother_type), intent(in) :: sm
  integer(psb_ipk_), intent(out)               :: info
  integer(psb_ipk_), intent(in), optional      :: iout
  logical, intent(in), optional              :: coarse

  ! Local variables
  integer(psb_ipk_)      :: err_act
  character(len=20), parameter :: name='mld_z_jac_smoother_descr'
  integer(psb_ipk_)      :: iout_
  logical      :: coarse_

  call psb_erractionsave(err_act)
  info = psb_success_
  if (present(coarse)) then 
    coarse_ = coarse
  else
    coarse_ = .false.
  end if
  if (present(iout)) then 
    iout_ = iout 
  else
    iout_ = psb_out_unit
  endif

  if (.not.coarse_) then
    if (allocated(sm%sv)) then
      select type(smv=>sm%sv)
      class is (mld_z_diag_solver_type)
        write(iout_,*) '  Point Jacobi smoother '
      class is (mld_z_bwgs_solver_type)
        write(iout_,*) '  Hybrid Backward Gauss-Seidel smoother '
!!$        write(iout_,*) '       Local Gauss-Seidel solver details:'
!!$        call smv%descr(info,iout_,coarse=coarse)
      class is (mld_z_gs_solver_type)
        write(iout_,*) '  Hybrid Forward Gauss-Seidel smoother '
!!$        write(iout_,*) '       Local Gauss-Seidel solver details:'
!!$        call smv%descr(info,iout_,coarse=coarse)
      class default
        write(iout_,*) '  Block Jacobi smoother '
        write(iout_,*) '       Local solver details:'
        call smv%descr(info,iout_,coarse=coarse)
      end select
      
    else
      write(iout_,*) '  Block Jacobi smoother '
    end if
  else
    if (allocated(sm%sv)) then 
      call sm%sv%descr(info,iout_,coarse=coarse)
    end if
  end if
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine mld_z_jac_smoother_descr
