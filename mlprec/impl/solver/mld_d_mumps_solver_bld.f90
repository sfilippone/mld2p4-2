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


  subroutine d_mumps_solver_bld(a,desc_a,sv,upd,info,b,amold,vmold,imold)

    use psb_base_mod
    use mpi
    use mld_d_mumps_solver
    Implicit None

    ! Arguments
    type(psb_dspmat_type)                               :: c 
    type(psb_dspmat_type), intent(in), target           :: a
    Type(psb_desc_type), Intent(in)                     :: desc_a 
    class(mld_d_mumps_solver_type), intent(inout)       :: sv
    character, intent(in)                               :: upd
    integer, intent(out)                                :: info
    type(psb_dspmat_type), intent(in), target, optional :: b
    class(psb_d_base_sparse_mat), intent(in), optional  :: amold
    class(psb_d_base_vect_type), intent(in), optional   :: vmold
    class(psb_i_base_vect_type), intent(in), optional   :: imold
    ! Local variables
    type(psb_dspmat_type)      :: atmp
    type(psb_d_coo_sparse_mat), target :: acoo
    integer                    :: n_row,n_col, nrow_a, nztota, nglob, nzt, npr, npc
    integer                    :: ifrst, ibcheck
    integer                    :: ictxt, ictxt1, icomm, np, me, i, err_act, debug_unit, debug_level
    character(len=20)          :: name='d_mumps_solver_bld', ch_err

    info=psb_success_

    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    ictxt       = desc_a%get_context()
    if (sv%ipar(9) < 0 ) then
        call psb_info(ictxt, me, np)
    	call psb_init(ictxt1,np=1,basectxt=ictxt,ids=(/me/))
    	call psb_get_mpicomm(ictxt1, icomm)
    	write(*,*)'mumps_bld: +++++>',icomm,ictxt1,mpi_comm_world
    	call psb_info(ictxt1, me, np)
    	npr  = np
    else 
        call psb_get_mpicomm(ictxt,icomm)
    	write(*,*)'mumps_bld: +++++>',icomm,ictxt,mpi_comm_world
    	call psb_info(ictxt, me, np)
    	npr  = np
    end if
    npc  = 1
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' start'
!    if (allocated(sv%id)) then 
!      call sv%free(info)

 !     deallocate(sv%id)
 !   end if
     if(.not.allocated(sv%id)) then
      allocate(sv%id,stat=info)
      if (info /= psb_success_) then
        info=psb_err_alloc_dealloc_
        call psb_errpush(info,name,a_err='mld_dmumps_default')
        goto 9999
      end if        
     end if


    if (psb_toupper(upd) == 'F') then 

      sv%id%comm    =  icomm
      sv%id%job = -1
      sv%id%par=1
      call dmumps(sv%id)   
      !WARNING: CALLING DMUMPS WITH JOB=-1 DESTROY THE SETTING OF DEFAULT:TO FIX
      sv%id%icntl(14)=sv%ipar(1)
      sv%id%keep(487)=sv%ipar(2)
      sv%id%keep(486)=sv%ipar(3)
      sv%id%keep(488)=sv%ipar(4)
      sv%id%keep(489)=sv%ipar(5)
      sv%id%keep(490)=sv%ipar(6)
      sv%id%keep(491)=sv%ipar(7)
      sv%id%keep(492)=sv%ipar(8)
      sv%id%dkeep(8)=sv%rpar(1)
      sv%id%cntl(2)=sv%rpar(2)
      sv%id%icntl(1)=-1
      !sv%id%icntl(2)=1
      !sv%id%icntl(3)=1
      !sv%id%icntl(4)=2
      nglob  = desc_a%get_global_rows()
      if (sv%ipar(9) < 0) then
	nglob=desc_a%get_local_rows()
      	call a%csclip(c,info,jmax=a%get_nrows())
      	call c%cp_to(acoo)
        nglob = c%get_nrows()
      else
      	call a%cp_to(acoo)
      end if
      nztota = acoo%get_nzeros()
      
      ! switch to global numbering
      if (sv%ipar(9) >= 0 ) then
      	call psb_loc_to_glob(acoo%ja(1:nztota), desc_a, info, iact='I')
      	call psb_loc_to_glob(acoo%ia(1:nztota), desc_a, info, iact='I')
      end if
      sv%id%irn_loc=> acoo%ia
      sv%id%jcn_loc=> acoo%ja
      sv%id%a_loc=> acoo%val
      sv%id%icntl(18)=3
      if(acoo%is_upper() .or. acoo%is_lower()) then
         sv%id%sym = 2
      else
         sv%id%sym = 0
      end if
      sv%id%n       =  nglob
      ! there should be a better way for this
      sv%id%nz_loc  =  acoo%get_nzeros()
      sv%id%nz      =  acoo%get_nzeros()
      sv%id%job = 4
      call psb_barrier(ictxt)
      write(*,*)'calling mumps N,nz,nz_loc',sv%id%n,sv%id%nz,sv%id%nz_loc
      call dmumps(sv%id)
      call psb_barrier(ictxt)
      info = sv%id%infog(1)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='mld_dmumps_fact '
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
      nullify(sv%id%irn)
      nullify(sv%id%jcn)
      nullify(sv%id%a)

      call acoo%free()
      sv%built=.true.
    else
      ! ? 
        info=psb_err_internal_error_
        call psb_errpush(info,name)
        goto 9999
      
    end if

    if (mld_as_sequential_ < 0) then
    	call psb_exit(ictxt1,close=.false.)	
    end if
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' end'

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_mumps_solver_bld
