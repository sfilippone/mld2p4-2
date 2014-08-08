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
!
!
!
!
!
!

module mld_d_mumps_solver
  use dmumps_struc_def
  use mld_d_base_solver_mod

#if defined(LONG_INTEGERS)

  type, extends(mld_d_base_solver_type) :: mld_d_mumps_solver_type

  end type mld_d_mumps_solver_type
#else
  type, extends(mld_d_base_solver_type) :: mld_d_mumps_solver_type
     type(dmumps_struc), allocatable  :: id
  contains
    procedure, pass(sv) :: build   => d_mumps_solver_bld
    procedure, pass(sv) :: apply_a => d_mumps_solver_apply
    procedure, pass(sv) :: apply_v => d_mumps_solver_apply_vect
    procedure, pass(sv) :: free    => d_mumps_solver_free
    procedure, pass(sv) :: descr   => d_mumps_solver_descr
    procedure, pass(sv) :: sizeof  => d_mumps_solver_sizeof
    procedure, pass(sv) :: seti    => d_mumps_solver_seti
    procedure, pass(sv) :: setr    => d_mumps_solver_setr
    procedure, pass(sv) :: cseti    => d_mumps_solver_cseti
    procedure, pass(sv) :: csetr    => d_mumps_solver_csetr
    procedure, pass(sv) :: default  => d_mumps_solver_default
#if defined(HAVE_FINAL) 

    final               :: d_mumps_solver_finalize
#endif
  end type mld_d_mumps_solver_type


  private :: d_mumps_solver_bld, d_mumps_solver_apply, &
       &  d_mumps_solver_free,   d_mumps_solver_descr, &
       &  d_mumps_solver_sizeof, d_mumps_solver_apply_vect,&
       &  d_mumps_solver_seti,   d_mumps_solver_setr,    &
       &  d_mumps_solver_cseti, d_mumps_solver_csetri,   &
       &  d_mumps_solver_default
#if defined(HAVE_FINAL) 
  private :: d_mumps_solver_finalize
#endif


contains

  subroutine d_mumps_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info)
    use psb_base_mod
    implicit none 
    type(psb_desc_type), intent(in)      :: desc_data
    class(mld_d_mumps_solver_type), intent(inout) :: sv
    real(psb_dpk_),intent(inout)         :: x(:)
    real(psb_dpk_),intent(inout)         :: y(:)
    real(psb_dpk_),intent(in)            :: alpha,beta
    character(len=1),intent(in)          :: trans
    real(psb_dpk_),target, intent(inout) :: work(:)
    integer, intent(out)                 :: info

    integer    :: n_row, n_col, nglob
    real(psb_dpk_), pointer     :: ww(:)
    real(psb_dpk_), allocatable, target :: gx(:)
    integer    :: ictxt,np,me,i, err_act
    character          :: trans_
    character(len=20)  :: name='d_mumps_solver_apply'

    call psb_erractionsave(err_act)

    info = psb_success_
    trans_ = psb_toupper(trans)
    select case(trans_)
    case('N')
    case('T')
    case default
      call psb_errpush(psb_err_iarg_invalid_i_,name)
      goto 9999
    end select

    nglob = desc_data%get_global_rows()
    n_row = desc_data%get_local_rows()
    n_col = desc_data%get_local_cols()

    if (n_col <= size(work)) then 
      ww => work(1:n_col)
    else
      allocate(ww(n_col),stat=info)
      if (info /= psb_success_) then 
        info=psb_err_alloc_request_
        call psb_errpush(info,name,i_err=(/n_col,0,0,0,0/),&
             & a_err='complex(psb_spk_)')
        goto 9999      
      end if
    end if
    allocate(gx(nglob),stat=info)
    if (info /= psb_success_) then 
       info=psb_err_alloc_request_
       call psb_errpush(info,name,i_err=(/nglob,0,0,0,0/),&
             & a_err='complex(psb_spk_)')
       goto 9999      
    end if
    call psb_gather(gx, x, desc_data, info, root=0)
    select case(trans_)
    case('N')
      sv%id%icntl(9) = 1
    case('T')
      sv%id%icntl(9) = 2
    case default
      call psb_errpush(psb_err_internal_error_,&
           & name,a_err='Invalid TRANS in subsolve')
      goto 9999
    end select

    sv%id%rhs  => gx
    sv%id%nrhs =  1
    sv%id%icntl(1)=-1
    sv%id%icntl(2)=-1
    sv%id%icntl(3)=-1
    sv%id%icntl(4)=-1
    sv%id%job = 3
    call dmumps(sv%id)
    call psb_scatter(gx, ww, desc_data, info, root=0)
    
    if (info == psb_success_) then
          call psb_geaxpby(alpha,ww,beta,y,desc_data,info)
      end if

    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,&
           & name,a_err='Error in subsolve')
      goto 9999
    endif

    if (nglob > size(work)) then 
      deallocate(ww)
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine d_mumps_solver_apply

  subroutine d_mumps_solver_apply_vect(alpha,sv,x,beta,y,desc_data,trans,work,info)
    use psb_base_mod
    implicit none 
    type(psb_desc_type), intent(in)      :: desc_data
    class(mld_d_mumps_solver_type), intent(inout) :: sv
    type(psb_d_vect_type),intent(inout)  :: x
    type(psb_d_vect_type),intent(inout)  :: y
    real(psb_dpk_),intent(in)            :: alpha,beta
    character(len=1),intent(in)          :: trans
    real(psb_dpk_),target, intent(inout) :: work(:)
    integer, intent(out)                 :: info

    integer    :: err_act
    character(len=20)  :: name='d_mumps_solver_apply_vect'

    call psb_erractionsave(err_act)

    info = psb_success_

    call x%v%sync()
    call y%v%sync()
    call sv%apply(alpha,x%v%v,beta,y%v%v,desc_data,trans,work,info)
    call y%v%set_host()
    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine d_mumps_solver_apply_vect

  subroutine d_mumps_solver_bld(a,desc_a,sv,upd,info,b,amold,vmold,imold)

    use psb_base_mod
    use mpi
    Implicit None

    ! Arguments
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
    integer                    :: ictxt, icomm, np, me, i, err_act, debug_unit, debug_level
    character(len=20)          :: name='d_mumps_solver_bld', ch_err

    info=psb_success_

    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    ictxt       = desc_a%get_context()
    call psb_get_mpicomm(ictxt, icomm)
    write(*,*)'mumps_bld: +++++>',icomm,ictxt,mpi_comm_world
    call psb_info(ictxt, me, np)
    npr  = np
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
      call sv%default

       nglob  = desc_a%get_global_rows()
      
      call a%cp_to(acoo)
      nztota = acoo%get_nzeros()
      
      ! switch to global numbering
      call psb_loc_to_glob(acoo%ja(1:nztota), desc_a, info, iact='I')
      call psb_loc_to_glob(acoo%ia(1:nztota), desc_a, info, iact='I')

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
      write(*,*)'calling mumps N,nz,nz_loc',sv%id%n,sv%id%nz,sv%id%nz_loc
      call dmumps(sv%id)
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
    else
      ! ? 
        info=psb_err_internal_error_
        call psb_errpush(info,name)
        goto 9999
      
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

  subroutine d_mumps_solver_free(sv,info)

    Implicit None

    ! Arguments
    class(mld_d_mumps_solver_type), intent(inout) :: sv
    integer, intent(out)                       :: info
    Integer :: err_act
    character(len=20)  :: name='d_mumps_solver_free'

    call psb_erractionsave(err_act)
    if (allocated(sv%id)) then 
      sv%id%job = -2
      call dmumps(sv%id)
      info = sv%id%infog(1)
      if (info /= psb_success_) goto 9999
      deallocate(sv%id)
    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_mumps_solver_free

#if defined(HAVE_FINAL)
  subroutine d_mumps_solver_finalize(sv)

    Implicit None

    ! Arguments
    type(mld_d_mumps_solver_type), intent(inout) :: sv 
    integer :: info
    Integer :: err_act
    character(len=20)  :: name='d_mumps_solver_finalize'

    call sv%free(info) 

    return
  
  end subroutine d_mumps_solver_finalize
#endif

  subroutine d_mumps_solver_descr(sv,info,iout,coarse)

    Implicit None

    ! Arguments
    class(mld_d_mumps_solver_type), intent(in) :: sv
    integer, intent(out)                     :: info
    integer, intent(in), optional            :: iout
    logical, intent(in), optional       :: coarse

    ! Local variables
    integer      :: err_act
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_d_mumps_solver_descr'
    integer :: iout_

    call psb_erractionsave(err_act)
    info = psb_success_
    if (present(iout)) then 
      iout_ = iout 
    else
      iout_ = 6
    endif
    
    write(iout_,*) '  MUMPS  Solver. '

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_mumps_solver_descr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$ LIST OF PARAMETERS FOR SUBROUTINE SET:ESEENTIALLY PARAMETER FOR THE USE  $!!
!!$ OF LOW RANK :                                                            $!!
!!$ WORKINGSPACE matching ICNTL(14);controls the percentage increase in the  $!!
!!$ estimated working space.Possible value:>0                                $!!
!!$ CLUSTERING matching KEEP(486):describes the clustering strategy and
!!$ activates BLR.Possible values:0 BLR is not activated; 1 BLR is activated
!!$  with inherit clustering                                                 !!$
!!$ HALO_DEPTH matching keep(487):the halo depth used in the clustering      $!!
!!$ operation. Possible values:                                              !!$
!!$ 0 : the clustering is deactivated, i.e. variables are kept in the original
!!$ order given by the global ordering of the matrix                         !!$
!!$ >0:user value                                                            !!$
!!$ TARGET_CLUSTER_SIZE matching dkeep(488):describes the clustering strategy!!$
!!$ and activates BLR. Posisble value >0                                     !!$
!!$ ALGORITHM matching KEEP(489):the factorization algorithm used within     !!$
!!$ fronts performed with BLR methods. Possible values:                      !!$
!!$ 0 : FSUUC (not yet available since it will be profitable only when the   !!$
!solution phase based on BLR blocks will be implemented).                    !!$
!!$ 1 : FCSUU (not implemented because not compatible with pivoting).        !!$
!!$ 2 : FSCUU without CB compression.                                        !!$
!!$ 3 : FSCUU with CB compression.                                           !!$
!!$ NASS_MIN matching with KEEP(490):the minimum number of assembled variables!$
!!s (NASS MIN) for a node to be selected for BLR. Possible values:>0         !!$
!!$ NFRONT_MIN matching KEEP(491):the minimum front size (NFRONT MIN) for a  !!$
!!$ node to be selected for BLR.Posisble values:>0                           !!$
!!$ SELECT_FRONT matching KEEP(492):describes how fronts are selected for BLR!!$
!!$ Possible values :                                                        !!$
!!$ < 0 : only front number |KEEP(492)| is selected                          !!$
!!$ 0 : none of the fronts are processed with BLR                            !!$
!!$ > 0 : all the fronts matching the criteria defined by KEEP(490) and      !!$
!!$ KEEP(491) are selected.                                                  !!$
!!$ QR_EPSILON matching DKEEP(8):the dropping parameter used for the QR      !!$
!!$ compression (ε) expressed with a double precision, real value            !!$
!!$ Possible values :                                                        !!$
!!$ > 0 : the dropping parameter is DKEEP(8)                                 !!$
!!$ < 0 : the dropping parameter is |DKEEP(8)| × ||A||                       !!$
!!$ STOPPING_CRITERION matching with CNTL(2):the stopping criterion for      !!$
!!$ iterative refinement                                                     !!$
!!$ WARNING: OTHERS PARAMETERS OF MUMPS COULD BE ADDED. FOR THIS, ADD AN     !!$
!!$ INTEGER IN MLD_BASE_PREC_TYPE.F90 AND MODIFY SUBROUTINE SET              !!$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  subroutine d_mumps_solver_seti(sv,what,val,info)

    Implicit None

    ! Arguments
    class(mld_d_mumps_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(in)                 :: what
    integer(psb_ipk_), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_)  :: err_act
    character(len=20)  :: name='d_mumps_solver_seti'

    info = psb_success_
    call psb_erractionsave(err_act)
    select case(what)
    case(mld_workspace_)
     sv%id%icntl(14)=val
    case(mld_halo_depth_)
     sv%id%keep(487)=val
    case(mld_clustering_)
     sv%id%keep(486)=val
    case(mld_cluster_size_)
     sv%id%keep(488)=val
    case(mld_algorithm_)
     sv%id%keep(489)=val
    case(mld_nass_min_)
     sv%id%keep(490)=val
    case(mld_nfront_min_)
     sv%id%keep(491)=val
    case(mld_select_front_)
     sv%id%keep(492)=val
    case default
      call sv%mld_d_base_solver_type%set(what,val,info)
    end select

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_mumps_solver_seti


  subroutine d_mumps_solver_setr(sv,what,val,info)

    Implicit None

    ! Arguments
    class(mld_d_mumps_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(in)                 :: what
    real(psb_dpk_), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_)  :: err_act
    character(len=20)  :: name='d_mumps_solver_setr'

    info = psb_success_
    call psb_erractionsave(err_act)

    select case(what)
    case(mld_qr_eps_)
     sv%id%dkeep(8)=val
    case(mld_stop_criterion_)
     sv%id%cntl(2)=val
    case default
      call sv%mld_d_base_solver_type%set(what,val,info)
    end select

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_mumps_solver_setr

  subroutine d_mumps_solver_cseti(sv,what,val,info)

    Implicit None

    ! Arguments
    class(mld_d_mumps_solver_type), intent(inout) :: sv
    character(len=*), intent(in)                  :: what
    integer(psb_ipk_), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_)  :: err_act, iwhat
    character(len=20)  :: name='d_mumps_solver_cseti'

    info = psb_success_
    call psb_erractionsave(err_act)

      select case(psb_toupper(what))
      case('WORKINGSPACE')
        iwhat=mld_workspace_
      case('HALO_DEPTH')
        iwhat=mld_halo_depth_
      case('TARGET_CLUSTER_SIZE')
        iwhat=mld_cluster_size_
      case('CLUSTERING')
        iwhat=mld_clustering_
      case('ALGORITHM')
        iwhat=mld_algorithm_
      case('NASS_MIN')
        iwhat=mld_nass_min_
      case('NFRONT_MIN')
        iwhat=mld_nfront_min_
      case('SELECT_FRONT')
        iwhat=mld_select_front_
      case default
        iwhat=-1
      end select

      if (iwhat >=0 ) then 
        call sv%set(iwhat,val,info)
      else
        call sv%mld_d_base_solver_type%set(what,val,info)
      end if
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_mumps_solver_cseti

subroutine d_mumps_solver_csetr(sv,what,val,info)

    Implicit None

    ! Arguments
    class(mld_d_mumps_solver_type), intent(inout) :: sv
    character(len=*), intent(in)                  :: what
    real(psb_dpk_), intent(in)                    :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_)  :: err_act, iwhat
    character(len=20)  :: name='d_mumps_solver_csetr'

    info = psb_success_
    call psb_erractionsave(err_act)

    select case(psb_toupper(what))
    case('QR_EPSILON')
      iwhat=mld_qr_eps_
    case('STOPPING_CRITERION')
      iwhat=mld_stop_criterion_
    case default
     call sv%mld_d_base_solver_type%set(what,val,info)
    end select
     
    if (iwhat >=0 ) then 
        call sv%set(iwhat,val,info)
    else
        call sv%mld_d_base_solver_type%set(what,val,info)
     end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_mumps_solver_csetr

!!NOTE: BY DEFAULT BLR is activated with a dropping parameter to 1d-4       !!
  subroutine d_mumps_solver_default(sv)

    Implicit none

    !Argument
    class(mld_d_mumps_solver_type),intent(inout) :: sv
    integer(psb_ipk_) :: info
    integer(psb_ipk_)  :: err_act,ictx,icomm
    character(len=20)  :: name='d_mumps_default'

    info = psb_success_
    call psb_erractionsave(err_act)
    
    if (.not.allocated(sv%id)) then 
      allocate(sv%id,stat=info)
      if (info /= psb_success_) then
        info=psb_err_alloc_dealloc_
        call psb_errpush(info,name,a_err='mld_dmumps_default')
        goto 9999
      end if
    end if

   !INSTANCIATION OF sv%id needed to set parmater but mpi communicator needed
   ! sv%id%job = -1
   ! sv%id%par=1
   ! call dmumps(sv%id)    
    !activation of Block Low rank factorization 
    sv%id%keep(486)=1
    !dropping paramater 
    sv%id%dkeep(8)=1d-4

    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine d_mumps_solver_default 

  function d_mumps_solver_sizeof(sv) result(val)

    implicit none 
    ! Arguments
    class(mld_d_mumps_solver_type), intent(in) :: sv
    integer(psb_long_int_k_) :: val
    integer             :: i
    val = 0
    ! val = 2*psb_sizeof_int + psb_sizeof_dp
    ! val = val + sv%symbsize
    ! val = val + sv%numsize
    return
  end function d_mumps_solver_sizeof
#endif
end module mld_d_mumps_solver
