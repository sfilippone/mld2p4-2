!!$
!!$ 
!!$                                MLD2P4
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS v.2.0)
!!$  
!!$  (C) Copyright 2007  Alfredo Buttari      University of Rome Tor Vergata
!!$                      Pasqua D'Ambra       ICAR-CNR, Naples
!!$                      Daniela di Serafino  Second University of Naples
!!$                      Salvatore Filippone  University of Rome Tor Vergata       
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
! File: mld_prec_type.f90.
!
! Module: mld_prec_type.
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

module mld_prec_type

  !
  ! This reduces the size of .mod file. Without the ONLY clause compilation 
  ! blows up on some systems.
  !
  use psb_base_mod, only : psb_dspmat_type, psb_zspmat_type, psb_desc_type,&
       & psb_sizeof

  !
  !  mld_dprec_type and mld_zprec_type are the real and complex      preconditioner
  !  data structures. In the following description 'd' and 'z' are omitted.
  !
  ! The multilevel preconditioner data structure, mld_prec_type, consists
  ! of an array of 'base preconditioner' data structures, mld_dbaseprc_type,
  ! each containing the local part of the preconditioner associated      to a
  ! certain level. For each level ilev, the base preconditioner K(ilev) is 
  ! built from a matrix A(ilev), which is obtained by 'tranferring' the 
  ! original matrix A (i.e. the matrix to be preconditioned) to level ilev,
  ! through smoothed aggregation.
  !
  ! The levels are numbered in increasing order      starting from the finest
  ! one, i.e. level 1 is the finest level and A(1) is the matrix A.
  !
  ! Details on mld_prec_type:
  !
  !         baseprecv  -  type(mld_dbaseprc_type), dimension(:), allocatable.
  !                  baseprecv(ilev) is the base preconditioner      at level ilev.
  !
  !    Note that number of levels = size(baseprecv(:)).
  !                 
  ! Details on mld_baseprc_type:
  !
  !    av         -  type(psb_dspmat_type), dimension(:), allocatable(:).
  !                  The sparse matrices needed to apply the preconditioner at
  !                  the current level ilev. 
  !      av(mld_l_pr_)     -  The L factor of the ILU factorization of the local
  !                           diagonal block of A(ilev).
  !      av(mld_u_pr_)     -  The U factor of the ILU factorization of the local
  !                           diagonal block of A(ilev), except its diagonal entries
  !                           (stored in d).
  !      av(mld_ap_nd_)    -  The entries of the local part of A(ilev) outside
  !                           the diagonal block, for block-Jacobi sweeps.
  !      av(mld_ac_)       -  The local part of the matrix A(ilev).
  !      av(mld_sm_pr_)    -  The smoother prolongator.   
  !                           It maps vectors (ilev) ---> (ilev-1).
  !      av(mld_sm_pr_t_)  -  The smoother prolongator transpose.   
  !                           It maps vectors (ilev-1) ---> (ilev).
  !      Shouldn't we keep just one of the last two items and handle the transpose
  !      in the Sparse BLAS? Maybe.
  !
  !   d            -  real(kind(1.d0)), dimension(:), allocatable.
  !                   The diagonal entries of the U factor in the ILU factorization
  !                   of A(ilev).
  !   desc_data    -  type(psb_desc_type).
  !                   The communication descriptor associated to the base preconditioner,
  !                   i.e. to the sparse matrices needed to apply the base preconditioner
  !                   at the current level.
  !   desc_ac      -  type(psb_desc_type).
  !                   The communication descriptor associated to the sparse matrix
  !                   A(ilev), stored in av(mld_ac_).
  !   iprcparm     -  integer, dimension(:), allocatable.
  !                   The integer parameters defining the base preconditioner K(ilev)
  !                   (the iprcparm entries and values are specified below).
  !   dprcparm     -  real(kind(1.d0)), dimension(:), allocatable.
  !                   The real parameters defining the base preconditioner K(ilev)
  !                   (the dprcparm entries and values are specified below).
  !   perm         -  integer, dimension(:), allocatable.
  !                   The row and column permutations applied to the local part of
  !                   A(ilev) (defined only if iprcparm(mld_sub_ren_)>0). 
  !   invperm      -  integer, dimension(:), allocatable.
  !                   The inverse of the permutation stored in perm.
  !   mlia         -  integer, dimension(:), allocatable.
  !                   The aggregation map (ilev-1) --> (ilev).
  !                   In case of non-smoothed aggregation, it is used instead of
  !                   mld_sm_pr_.
  !   nlaggr       -  integer, dimension(:), allocatable.
  !                   The number of aggregates (rows of A(ilev)) on the
  !                   various processes. 
  !   base_a       -  type(psb_zspmat_type), pointer.
  !                   Pointer (really a pointer!) to the local part of the base matrix 
  !                   of the current level, i.e. A(ilev); so we have a unified treatment
  !                   of residuals. We need this to avoid passing explicitly the matrix
  !                   A(ilev) to the routine which applies the preconditioner.
  !   base_desc    -  type(psb_desc_type), pointer.
  !                   Pointer to the communication descriptor associated to the sparse
  !                   matrix pointed by base_a. 
  !   dorig        -  real(kind(1.d0)), dimension(:), allocatable.
  !                   Diagonal entries of the matrix pointed by base_a.
  !
  !   Note that when the LU factorization of the matrix A(ilev) is computed instead of
  !   the ILU one, by using UMFPACK or SuperLU_dist, the corresponding L and U factors
  !   are stored in data structures provided by UMFPACK or SuperLU_dist and pointed by
  !   iprcparm(mld_umf_ptr) or iprcparm(mld_slu_ptr), respectively.
  !

  !
  ! Preconditioner data types
  !

  type mld_dbaseprc_type

    type(psb_dspmat_type), allocatable :: av(:) 
    real(kind(1.d0)), allocatable      :: d(:)  
    type(psb_desc_type)                :: desc_data , desc_ac
    integer, allocatable               :: iprcparm(:) 
    real(kind(1.d0)), allocatable      :: dprcparm(:) 
    integer, allocatable               :: perm(:),  invperm(:) 
    integer, allocatable               :: mlia(:), nlaggr(:) 
    type(psb_dspmat_type), pointer     :: base_a    => null() 
    type(psb_desc_type), pointer       :: base_desc => null() 
    real(kind(1.d0)), allocatable      :: dorig(:) 

  end type mld_dbaseprc_type

  type mld_dprec_type
    type(mld_dbaseprc_type), allocatable  :: baseprecv(:) 
  end type mld_dprec_type

  type mld_zbaseprc_type

    type(psb_zspmat_type), allocatable :: av(:) 
    complex(kind(1.d0)), allocatable   :: d(:)  
    type(psb_desc_type)                :: desc_data , desc_ac
    integer, allocatable               :: iprcparm(:) 
    real(kind(1.d0)), allocatable      :: dprcparm(:) 
    integer, allocatable               :: perm(:),  invperm(:) 
    integer, allocatable               :: mlia(:), nlaggr(:) 
    type(psb_zspmat_type), pointer     :: base_a    => null() 
    type(psb_desc_type), pointer       :: base_desc => null() 
    complex(kind(1.d0)), allocatable   :: dorig(:)

  end type mld_zbaseprc_type

  type mld_zprec_type
    type(mld_zbaseprc_type), allocatable  :: baseprecv(:) 
  end type mld_zprec_type

  !
  ! Entries in iprcparm
  !
  integer, parameter :: mld_prec_type_=1         
  integer, parameter :: mld_sub_solve_=2
  integer, parameter :: mld_sub_restr_=3
  integer, parameter :: mld_sub_prol_=4
  integer, parameter :: mld_sub_ren_=5
  integer, parameter :: mld_n_ovr_=6
  integer, parameter :: mld_sub_fill_in_=8
  integer, parameter :: mld_smooth_sweeps_=9
  integer, parameter :: mld_ml_type_=10
  integer, parameter :: mld_smooth_pos_=11
  integer, parameter :: mld_aggr_alg_=12
  integer, parameter :: mld_aggr_kind_=13
  integer, parameter :: mld_aggr_eig_=14
  integer, parameter :: mld_coarse_mat_=16
  !! 2 ints for 64 bit versions
  integer, parameter :: mld_slu_ptr_=17
  integer, parameter :: mld_umf_symptr_=17
  integer, parameter :: mld_umf_numptr_=19
  integer, parameter :: mld_slud_ptr_=21
  integer, parameter :: mld_prec_status_=24
  integer, parameter :: mld_coarse_solve_  =25 
  integer, parameter :: mld_coarse_sweeps_ =26
  integer, parameter :: mld_coarse_fill_in_=27
  integer, parameter :: mld_ifpsz_=32

  !
  ! Legal values for entry: mld_prec_type_
  ! 
  integer, parameter :: mld_min_prec_=0, mld_noprec_=0, mld_diag_=1, mld_bjac_=2,&
       & mld_as_=3, mld_max_prec_=3
  !
  ! Legal values for entry: mld_ml_type_
  !
  integer, parameter :: mld_no_ml_=0, mld_add_ml_=1, mld_mult_ml_=2
  integer, parameter :: mld_new_ml_prec_=3, mld_max_ml_=mld_new_ml_prec_
  !
  ! Legal values for entry: mld_smooth_pos_
  !
  integer, parameter :: mld_pre_smooth_=1, mld_post_smooth_=2,&
       &  mld_twoside_smooth_=3, mld_max_smooth_=mld_twoside_smooth_
  !
  ! Legal values for entry: mld_sub_solve_
  !
  integer, parameter :: mld_f_none_=0,mld_ilu_n_=1,mld_milu_n_=2, mld_ilu_t_=3
  integer, parameter :: mld_slu_=4, mld_umf_=5, mld_sludist_=6
  !  
  ! Legal values for entry: mld_aggr_alg_
  !
  integer, parameter :: mld_dec_aggr_=0, mld_sym_dec_aggr_=1
  integer, parameter :: mld_glb_aggr_=2, mld_new_dec_aggr_=3
  integer, parameter :: mld_new_glb_aggr_=4, mld_max_aggr_=mld_new_glb_aggr_
  !
  ! Legal values for entry: mld_aggr_kind_
  !
  integer, parameter :: mld_no_smooth_=0, mld_smooth_prol_=1, mld_biz_prol_=2
  !
  ! Legal values for entry: mld_aggr_eig_
  !
  integer, parameter :: mld_max_norm_=0, mld_user_choice_=999
  !
  ! Legal values for entry: mld_coarse_mat_
  !
  integer, parameter :: mld_distr_mat_=0, mld_repl_mat_=1
  !
  ! Legal values for entry: mld_prec_status_
  !
  integer, parameter :: mld_prec_built_=98765
  !
  ! Legal values for entry: mld_sub_ren_
  !
  integer, parameter :: mld_renum_none_=0, mld_renum_glb_=1, mld_renum_gps_=2

  !
  ! Entries in dprcparm: ILU(k,t) threshold, smoothed aggregation omega
  !
  integer, parameter :: mld_fact_thrs_=1
  integer, parameter :: mld_aggr_damp_=2
  integer, parameter :: mld_aggr_thresh_=3
  integer, parameter :: mld_dfpsz_=4

  !
  ! Fields for sparse matrices ensembles stored in av()
  ! 
  integer, parameter :: mld_l_pr_=1, mld_u_pr_=2, mld_bp_ilu_avsz_=2
  integer, parameter :: mld_ap_nd_=3, mld_ac_=4, mld_sm_pr_t_=5, mld_sm_pr_=6
  integer, parameter :: mld_smth_avsz_=6, mld_max_avsz_=mld_smth_avsz_ 

  !
  ! Character constants used by      mld_file_prec_descr
  !
  character(len=15), parameter, private :: &
       &  smooth_names(1:3)=(/'Pre-smoothing ','Post-smoothing',&
       & 'Smooth both   '/)
  character(len=15), parameter, private :: &
       &  smooth_kinds(0:2)=(/'No  smoother  ','Omega smoother',&
       &           'Bizr. smoother'/)
  character(len=15), parameter, private :: &
       &  matrix_names(0:1)=(/'Distributed   ','Replicated    '/)
  character(len=18), parameter, private :: &
       &  aggr_names(0:4)=(/'Local aggregation ','Sym. local aggr.  ',&
       &     'Global aggregation', 'New local aggr.   ','New global aggr.  '/)
  character(len=6), parameter, private :: &
       &  restrict_names(0:4)=(/'None ','Halo ','     ','     ','     '/)
  character(len=12), parameter, private :: &
       &  prolong_names(0:3)=(/'None       ','Sum        ','Average    ','Square root'/)
  character(len=15), parameter, private :: &
       &  ml_names(0:3)=(/'None          ','Additive      ','Multiplicative',&
       & 'New ML        '/)
  character(len=15), parameter, private :: &
       &  fact_names(0:6)=(/'None          ','ILU(n)        ',&
       &  'MILU(n)       ','ILU(T)        ',&
       &  'Sparse SuperLU','UMFPACK Sp. LU',&
       &  'SuperLU_Dist  '/)

  !
  ! Interfaces to routines for checking the definition of the preconditioner,
  ! for printing its description and for deallocating its data structure
  !

  interface mld_base_precfree
    module procedure mld_dbase_precfree, mld_zbase_precfree
  end interface

  interface mld_nullify_baseprec
    module procedure mld_nullify_dbaseprec, mld_nullify_zbaseprec
  end interface

  interface mld_check_def
    module procedure mld_icheck_def, mld_dcheck_def
  end interface

  interface mld_prec_descr
    module procedure mld_out_prec_descr, mld_file_prec_descr, &
         &  mld_zout_prec_descr, mld_zfile_prec_descr
  end interface

  interface mld_prec_short_descr
    module procedure mld_prec_short_descr, mld_zprec_short_descr
  end interface

  interface mld_sizeof
    module procedure mld_dprec_sizeof, mld_zprec_sizeof, &
         & mld_dbaseprc_sizeof, mld_zbaseprc_sizeof
  end interface

contains

  !
  ! Function returning the size of the mld_prec_type data structure
  !

  function mld_dprec_sizeof(prec)
    use psb_base_mod
    type(mld_dprec_type), intent(in) :: prec
    integer             :: mld_dprec_sizeof
    integer             :: val,i
    val = 8
    if (allocated(prec%baseprecv)) then 
      do i=1, size(prec%baseprecv)
        val = val + mld_sizeof(prec%baseprecv(i))
      end do
    end if
    mld_dprec_sizeof = val
  end function mld_dprec_sizeof

  function mld_zprec_sizeof(prec)
    use psb_base_mod
    type(mld_zprec_type), intent(in) :: prec
    integer             :: mld_zprec_sizeof
    integer             :: val,i
    val = 8
    if (allocated(prec%baseprecv)) then 
      do i=1, size(prec%baseprecv)
        val = val + mld_sizeof(prec%baseprecv(i))
      end do
    end if
    mld_zprec_sizeof = val
  end function mld_zprec_sizeof

  !
  ! Function returning the size of the mld_baseprc_type data structure
  !

  function mld_dbaseprc_sizeof(prec)
    use psb_base_mod
    type(mld_dbaseprc_type), intent(in) :: prec
    integer             :: mld_dbaseprc_sizeof
    integer             :: val,i
    
    val = 0
    if (allocated(prec%iprcparm)) then 
      val = val + 4 * size(prec%iprcparm)
      if (prec%iprcparm(mld_prec_status_) == mld_prec_built_) then 
        select case(prec%iprcparm(mld_sub_solve_)) 
        case(mld_ilu_n_,mld_ilu_t_)
          ! do nothing
        case(mld_slu_)
          write(0,*) 'Should implement check for size of SuperLU data structs'
        case(mld_umf_)
          write(0,*) 'Should implement check for size of UMFPACK data structs'
        case(mld_sludist_)
          write(0,*) 'Should implement check for size of SuperLUDist data structs'
        case default
        end select
        
      end if
    end if
    if (allocated(prec%dprcparm)) val = val + 8 * size(prec%dprcparm)
    if (allocated(prec%d))        val = val + 8 * size(prec%d)
    if (allocated(prec%perm))     val = val + 4 * size(prec%perm)
    if (allocated(prec%invperm))  val = val + 4 * size(prec%invperm)
                                  val = val + psb_sizeof(prec%desc_data)
    if (allocated(prec%av))  then 
      do i=1,size(prec%av)
        val = val + psb_sizeof(prec%av(i))
      end do
    end if

    mld_dbaseprc_sizeof = val 
    
  end function mld_dbaseprc_sizeof

  function mld_zbaseprc_sizeof(prec)
    use psb_base_mod
    type(mld_zbaseprc_type), intent(in) :: prec
    integer             :: mld_zbaseprc_sizeof
    integer             :: val,i
    
    val = 0
    if (allocated(prec%iprcparm)) then 
      val = val + 4 * size(prec%iprcparm)
      if (prec%iprcparm(mld_prec_status_) == mld_prec_built_) then 
        select case(prec%iprcparm(mld_sub_solve_)) 
        case(mld_ilu_n_,mld_ilu_t_)
          ! do nothing
        case(mld_slu_)
          write(0,*) 'Should implement check for size of SuperLU data structs'
        case(mld_umf_)
          write(0,*) 'Should implement check for size of UMFPACK data structs'
        case(mld_sludist_)
          write(0,*) 'Should implement check for size of SuperLUDist data structs'
        case default
        end select
        
      end if
    end if
    if (allocated(prec%dprcparm)) val = val + 8 * size(prec%dprcparm)
    if (allocated(prec%d))        val = val + 16 * size(prec%d)
    if (allocated(prec%perm))     val = val + 4 * size(prec%perm)
    if (allocated(prec%invperm))  val = val + 4 * size(prec%invperm)
                                  val = val + psb_sizeof(prec%desc_data)
    if (allocated(prec%av))  then 
      do i=1,size(prec%av)
        val = val + psb_sizeof(prec%av(i))
      end do
    end if
    
    mld_zbaseprc_sizeof = val 
    
  end function mld_zbaseprc_sizeof
    
  !
  ! Routines printing out a description of the preconditioner
  !

  !
  ! Subroutine: mld_out_prec_descr.
  ! Version: real.
  !
  !  This routine prints to the standard output a description of the
  !  preconditioner.
  !
  ! Arguments:
  !       p       -  type(mld_dprec_type), input.
  !             The preconditioner data structure to be printed out.
  !
  subroutine mld_out_prec_descr(p)
    use psb_base_mod
    type(mld_dprec_type), intent(in) :: p
    call mld_file_prec_descr(6,p)
  end subroutine mld_out_prec_descr

  subroutine mld_zout_prec_descr(p)
    use psb_base_mod
    type(mld_zprec_type), intent(in) :: p
    call mld_zfile_prec_descr(6,p)
  end subroutine mld_zout_prec_descr

  !
  ! Subroutine: mld_file_prec_descr.
  ! Version: real.
  !
  !  This routine prints to a file a description of the preconditioner.
  !
  ! Arguments:
  !  iout    -  integer, input.
  !             The id of the file where the preconditioner description
  !             will be printed.
  !       p       -  type(mld_dprec_type), input.
  !             The preconditioner data structure to be printed out.
  !
  subroutine mld_file_prec_descr(iout,p)
    use psb_base_mod
    integer, intent(in)              :: iout
    type(mld_dprec_type), intent(in) :: p
    integer  :: ilev

    write(iout,*) 'Preconditioner description'
    if (allocated(p%baseprecv)) then 
      if (size(p%baseprecv)>=1) then 
        ilev = 1
        write(iout,*) 'Base preconditioner'
        select case(p%baseprecv(ilev)%iprcparm(mld_prec_type_))
        case(mld_noprec_)
          write(iout,*) 'No preconditioning'
        case(mld_diag_)
          write(iout,*) 'Diagonal scaling'
        case(mld_bjac_)
          write(iout,*) 'Block Jacobi with: ',&
               &  fact_names(p%baseprecv(ilev)%iprcparm(mld_sub_solve_))
          select case(p%baseprecv(ilev)%iprcparm(mld_sub_solve_))
          case(mld_ilu_n_,mld_milu_n_)      
            write(iout,*) 'Fill level:',p%baseprecv(ilev)%iprcparm(mld_sub_fill_in_)
          case(mld_ilu_t_)         
            write(iout,*) 'Fill threshold :',p%baseprecv(ilev)%dprcparm(mld_fact_thrs_)
          case(mld_slu_,mld_umf_,mld_sludist_) 
          case default
            write(iout,*) 'Should never get here!'
          end select
        case(mld_as_)
          write(iout,*) 'Additive Schwarz with: ',&
               &  fact_names(p%baseprecv(ilev)%iprcparm(mld_sub_solve_))
          select case(p%baseprecv(ilev)%iprcparm(mld_sub_solve_))
          case(mld_ilu_n_,mld_milu_n_)      
            write(iout,*) 'Fill level:',p%baseprecv(ilev)%iprcparm(mld_sub_fill_in_)
          case(mld_ilu_t_)         
            write(iout,*) 'Fill threshold :',p%baseprecv(ilev)%dprcparm(mld_fact_thrs_)
          case(mld_slu_,mld_umf_,mld_sludist_) 
          case default
            write(iout,*) 'Should never get here!'
          end select
          write(iout,*) 'Overlap:',&
               &  p%baseprecv(ilev)%iprcparm(mld_n_ovr_)
          write(iout,*) 'Restriction: ',&
               &  restrict_names(p%baseprecv(ilev)%iprcparm(mld_sub_restr_))
          write(iout,*) 'Prolongation: ',&
               &  prolong_names(p%baseprecv(ilev)%iprcparm(mld_sub_prol_))
        end select
      end if
      if (size(p%baseprecv)>=2) then 
        do ilev = 2, size(p%baseprecv) 
          if (.not.allocated(p%baseprecv(ilev)%iprcparm)) then 
            write(iout,*) 'Inconsistent MLPREC part!'
            return
          endif
          
          write(iout,*) 'Multilevel: Level No', ilev
          write(iout,*) 'Multilevel type: ',&
               &   ml_names(p%baseprecv(ilev)%iprcparm(mld_ml_type_))
          if (p%baseprecv(ilev)%iprcparm(mld_ml_type_)>mld_no_ml_) then 
            write(iout,*) 'Multilevel aggregation: ', &
                 &   aggr_names(p%baseprecv(ilev)%iprcparm(mld_aggr_alg_))
            write(iout,*) 'Aggregation smoothing: ', &
                 &  smooth_kinds(p%baseprecv(ilev)%iprcparm(mld_aggr_kind_))
            if (p%baseprecv(ilev)%iprcparm(mld_aggr_kind_) /= mld_no_smooth_) then 
              write(iout,*) 'Damping omega: ', &
                   & p%baseprecv(ilev)%dprcparm(mld_aggr_damp_)
              write(iout,*) 'Multilevel smoother position: ',&
                   & smooth_names(p%baseprecv(ilev)%iprcparm(mld_smooth_pos_))
            end if
            write(iout,*) 'Coarse matrix: ',&
                 & matrix_names(p%baseprecv(ilev)%iprcparm(mld_coarse_mat_))
            if (allocated(p%baseprecv(ilev)%nlaggr)) then 
              write(iout,*) 'Sizes of aggregates: ', &
                   &  sum( p%baseprecv(ilev)%nlaggr(:)),' : ',p%baseprecv(ilev)%nlaggr(:)
            end if
            write(iout,*) 'Factorization type: ',&
                 & fact_names(p%baseprecv(ilev)%iprcparm(mld_sub_solve_))
            select case(p%baseprecv(ilev)%iprcparm(mld_sub_solve_))
            case(mld_ilu_n_,mld_milu_n_)      
              write(iout,*) 'Fill level:',p%baseprecv(ilev)%iprcparm(mld_sub_fill_in_)
            case(mld_ilu_t_)         
              write(iout,*) 'Fill threshold :',p%baseprecv(ilev)%dprcparm(mld_fact_thrs_)
            case(mld_slu_,mld_umf_,mld_sludist_) 
            case default
              write(iout,*) 'Should never get here!'
            end select
            write(iout,*) 'Number of Jacobi sweeps: ', &
                 &   (p%baseprecv(ilev)%iprcparm(mld_smooth_sweeps_))
          end if
        end do
      end if

    else
      write(iout,*) 'No Base preconditioner available, something is wrong!'
      return
    endif

  end subroutine mld_file_prec_descr

  function  mld_prec_short_descr(p)
    use psb_base_mod
    type(mld_dprec_type), intent(in) :: p
    character(len=20) :: mld_prec_short_descr
    mld_prec_short_descr = ' '
!!$    write(iout,*) 'Preconditioner description'
!!$    if (associated(p%baseprecv)) then 
!!$      if (size(p%baseprecv)>=1) then 
!!$        write(iout,*) 'Base preconditioner'
!!$        select case(p%baseprecv(1)%iprcparm(mld_prec_type_))
!!$        case(mld_noprec_)
!!$          write(iout,*) 'No preconditioning'
!!$        case(mld_diag_)
!!$          write(iout,*) 'Diagonal scaling'
!!$        case(mld_bjac_)
!!$          write(iout,*) 'Block Jacobi with: ',&
!!$               &  fact_names(p%baseprecv(1)%iprcparm(mld_sub_solve_))
!!$        case(mld_as_,rmld_as_,ash_,rash_)
!!$          write(iout,*) 'Additive Schwarz with: ',&
!!$               &  fact_names(p%baseprecv(1)%iprcparm(mld_sub_solve_))
!!$          write(iout,*) 'Overlap:',&
!!$               &  p%baseprecv(1)%iprcparm(mld_n_ovr_)
!!$          write(iout,*) 'Restriction: ',&
!!$               &  restrict_names(p%baseprecv(1)%iprcparm(mld_sub_restr_))
!!$          write(iout,*) 'Prolongation: ',&
!!$               &  prolong_names(p%baseprecv(1)%iprcparm(mld_sub_prol_))
!!$        end select
!!$      end if
!!$      if (size(p%baseprecv)>=2) then 
!!$        if (.not.associated(p%baseprecv(2)%iprcparm)) then 
!!$          write(iout,*) 'Inconsistent MLPREC part!'
!!$          return
!!$        endif
!!$        write(iout,*) 'Multilevel: ',ml_names(p%baseprecv(2)%iprcparm(mld_ml_type_))
!!$        if (p%baseprecv(2)%iprcparm(mld_ml_type_)>mld_no_ml_) then 
!!$          write(iout,*) 'Multilevel aggregation: ', &
!!$               &   aggr_names(p%baseprecv(2)%iprcparm(mld_aggr_alg_))
!!$          write(iout,*) 'Multilevel smoothing: ', &
!!$               &  smooth_kinds(p%baseprecv(2)%iprcparm(mld_aggr_kind_))
!!$          write(iout,*) 'damping omega: ', p%baseprecv(2)%dprcparm(mld_aggr_damp_)
!!$          write(iout,*) 'Multilevel smoother position: ',&
!!$               & smooth_names(p%baseprecv(2)%iprcparm(mld_smooth_pos_))
!!$          write(iout,*) 'Coarse matrix: ',&
!!$               & matrix_names(p%baseprecv(2)%iprcparm(mld_coarse_mat_))
!!$          write(iout,*) 'Factorization type: ',&
!!$               & fact_names(p%baseprecv(2)%iprcparm(mld_sub_solve_))
!!$          select case(p%baseprecv(2)%iprcparm(mld_sub_solve_))
!!$          case(mld_ilu_n_)      
!!$            write(iout,*) 'Fill level:',p%baseprecv(2)%iprcparm(mld_sub_fill_in_)
!!$          case(mld_ilu_t_)         
!!$            write(iout,*) 'Fill threshold :',p%baseprecv(2)%dprcparm(mld_fact_thrs_)
!!$          case(mld_slu_,mld_umf_,mld_sludist_)         
!!$          case default
!!$            write(iout,*) 'Should never get here!'
!!$          end select
!!$          write(iout,*) 'Number of Jacobi sweeps: ', &
!!$               &   (p%baseprecv(2)%iprcparm(mld_smooth_sweeps_))
!!$
!!$        end if
!!$      end if
!!$
!!$    else
!!$      write(iout,*) 'No Base preconditioner available, something is wrong!'
!!$      return
!!$    endif

  end function mld_prec_short_descr


  subroutine mld_zfile_prec_descr(iout,p)
    use psb_base_mod
    integer, intent(in)              :: iout
    type(mld_zprec_type), intent(in) :: p
    integer  :: ilev

    write(iout,*) 'Preconditioner description'
    if (allocated(p%baseprecv)) then 
      if (size(p%baseprecv)>=1) then 
        write(iout,*) 'Base preconditioner'
        ilev=1
        select case(p%baseprecv(ilev)%iprcparm(mld_prec_type_))
        case(mld_noprec_)
          write(iout,*) 'No preconditioning'
        case(mld_diag_)
          write(iout,*) 'Diagonal scaling'
        case(mld_bjac_)
          write(iout,*) 'Block Jacobi with: ',&
               &  fact_names(p%baseprecv(ilev)%iprcparm(mld_sub_solve_))
            select case(p%baseprecv(ilev)%iprcparm(mld_sub_solve_))
            case(mld_ilu_n_,mld_milu_n_)      
              write(iout,*) 'Fill level:',p%baseprecv(ilev)%iprcparm(mld_sub_fill_in_)
            case(mld_ilu_t_)         
              write(iout,*) 'Fill threshold :',p%baseprecv(ilev)%dprcparm(mld_fact_thrs_)
            case(mld_slu_,mld_umf_,mld_sludist_) 
            case default
              write(iout,*) 'Should never get here!'
            end select
        case(mld_as_)
          write(iout,*) 'Additive Schwarz with: ',&
               &  fact_names(p%baseprecv(ilev)%iprcparm(mld_sub_solve_))
            select case(p%baseprecv(ilev)%iprcparm(mld_sub_solve_))
            case(mld_ilu_n_,mld_milu_n_)      
              write(iout,*) 'Fill level:',p%baseprecv(ilev)%iprcparm(mld_sub_fill_in_)
            case(mld_ilu_t_)         
              write(iout,*) 'Fill threshold :',p%baseprecv(ilev)%dprcparm(mld_fact_thrs_)
            case(mld_slu_,mld_umf_,mld_sludist_) 
            case default
              write(iout,*) 'Should never get here!'
            end select
          write(iout,*) 'Overlap:',&
               &  p%baseprecv(ilev)%iprcparm(mld_n_ovr_)
          write(iout,*) 'Restriction: ',&
               &  restrict_names(p%baseprecv(ilev)%iprcparm(mld_sub_restr_))
          write(iout,*) 'Prolongation: ',&
               &  prolong_names(p%baseprecv(ilev)%iprcparm(mld_sub_prol_))
        end select
      end if
      if (size(p%baseprecv)>=2) then 
        do ilev = 2, size(p%baseprecv) 
          if (.not.allocated(p%baseprecv(ilev)%iprcparm)) then 
            write(iout,*) 'Inconsistent MLPREC part!'
            return
          endif
          
          write(iout,*) 'Multilevel: Level No', ilev
          write(iout,*) 'Multilevel type: ',&
               &   ml_names(p%baseprecv(ilev)%iprcparm(mld_ml_type_))
          if (p%baseprecv(ilev)%iprcparm(mld_ml_type_)>mld_no_ml_) then 
            write(iout,*) 'Multilevel aggregation: ', &
                 &   aggr_names(p%baseprecv(ilev)%iprcparm(mld_aggr_alg_))
            write(iout,*) 'Smoother:               ', &
                 &  smooth_kinds(p%baseprecv(ilev)%iprcparm(mld_aggr_kind_))
            if (p%baseprecv(ilev)%iprcparm(mld_aggr_kind_) /= mld_no_smooth_) then 
              write(iout,*) 'Smoothing omega: ', &
                   & p%baseprecv(ilev)%dprcparm(mld_aggr_damp_)
              write(iout,*) 'Smoothing position: ',&
                   & smooth_names(p%baseprecv(ilev)%iprcparm(mld_smooth_pos_))
            end if
            write(iout,*) 'Coarse matrix: ',&
                 & matrix_names(p%baseprecv(ilev)%iprcparm(mld_coarse_mat_))
            if (allocated(p%baseprecv(ilev)%nlaggr)) then 
              write(iout,*) 'Aggregation sizes: ', &
                   &  sum( p%baseprecv(ilev)%nlaggr(:)),' : ',p%baseprecv(ilev)%nlaggr(:)
            end if
            write(iout,*) 'Factorization type: ',&
                 & fact_names(p%baseprecv(ilev)%iprcparm(mld_sub_solve_))
            select case(p%baseprecv(ilev)%iprcparm(mld_sub_solve_))
            case(mld_ilu_n_,mld_milu_n_)      
              write(iout,*) 'Fill level:',p%baseprecv(ilev)%iprcparm(mld_sub_fill_in_)
            case(mld_ilu_t_)         
              write(iout,*) 'Fill threshold :',p%baseprecv(ilev)%dprcparm(mld_fact_thrs_)
            case(mld_slu_,mld_umf_,mld_sludist_) 
            case default
              write(iout,*) 'Should never get here!'
            end select
            write(iout,*) 'Number of Jacobi sweeps: ', &
                 &   (p%baseprecv(ilev)%iprcparm(mld_smooth_sweeps_))
          end if
        end do
      end if

    else
      write(iout,*) 'No Base preconditioner available, something is wrong!'
      return
    endif

  end subroutine mld_zfile_prec_descr

  function  mld_zprec_short_descr(p)
    use psb_base_mod
    type(mld_zprec_type), intent(in) :: p
    character(len=20) :: mld_zprec_short_descr
    mld_zprec_short_descr = ' '
!!$    write(iout,*) 'Preconditioner description'
!!$    if (associated(p%baseprecv)) then 
!!$      if (size(p%baseprecv)>=1) then 
!!$        write(iout,*) 'Base preconditioner'
!!$        select case(p%baseprecv(1)%iprcparm(mld_prec_type_))
!!$        case(mld_noprec_)
!!$          write(iout,*) 'No preconditioning'
!!$        case(mld_diag_)
!!$          write(iout,*) 'Diagonal scaling'
!!$        case(mld_bjac_)
!!$          write(iout,*) 'Block Jacobi with: ',&
!!$               &  fact_names(p%baseprecv(1)%iprcparm(mld_sub_solve_))
!!$        case(mld_as_,rmld_as_,ash_,rash_)
!!$          write(iout,*) 'Additive Schwarz with: ',&
!!$               &  fact_names(p%baseprecv(1)%iprcparm(mld_sub_solve_))
!!$          write(iout,*) 'Overlap:',&
!!$               &  p%baseprecv(1)%iprcparm(mld_n_ovr_)
!!$          write(iout,*) 'Restriction: ',&
!!$               &  restrict_names(p%baseprecv(1)%iprcparm(mld_sub_restr_))
!!$          write(iout,*) 'Prolongation: ',&
!!$               &  prolong_names(p%baseprecv(1)%iprcparm(mld_sub_prol_))
!!$        end select
!!$      end if
!!$      if (size(p%baseprecv)>=2) then 
!!$        if (.not.associated(p%baseprecv(2)%iprcparm)) then 
!!$          write(iout,*) 'Inconsistent MLPREC part!'
!!$          return
!!$        endif
!!$        write(iout,*) 'Multilevel: ',ml_names(p%baseprecv(2)%iprcparm(mld_ml_type_))
!!$        if (p%baseprecv(2)%iprcparm(mld_ml_type_)>mld_no_ml_) then 
!!$          write(iout,*) 'Multilevel aggregation: ', &
!!$               &   aggr_names(p%baseprecv(2)%iprcparm(mld_aggr_alg_))
!!$          write(iout,*) 'Smoother:               ', &
!!$               &  smooth_kinds(p%baseprecv(2)%iprcparm(mld_aggr_kind_))
!!$          write(iout,*) 'Smoothing omega: ', p%baseprecv(2)%dprcparm(mld_aggr_damp_)
!!$          write(iout,*) 'Smoothing position: ',&
!!$               & smooth_names(p%baseprecv(2)%iprcparm(mld_smooth_pos_))
!!$          write(iout,*) 'Coarse matrix: ',&
!!$               & matrix_names(p%baseprecv(2)%iprcparm(mld_coarse_mat_))
!!$          write(iout,*) 'Factorization type: ',&
!!$               & fact_names(p%baseprecv(2)%iprcparm(mld_sub_solve_))
!!$          select case(p%baseprecv(2)%iprcparm(mld_sub_solve_))
!!$          case(mld_ilu_n_)      
!!$            write(iout,*) 'Fill level:',p%baseprecv(2)%iprcparm(mld_sub_fill_in_)
!!$          case(mld_ilu_t_)         
!!$            write(iout,*) 'Fill threshold :',p%baseprecv(2)%dprcparm(mld_fact_thrs_)
!!$          case(mld_slu_,mld_umf_,mld_sludist_)         
!!$          case default
!!$            write(iout,*) 'Should never get here!'
!!$          end select
!!$          write(iout,*) 'Number of Jacobi sweeps: ', &
!!$               &   (p%baseprecv(2)%iprcparm(mld_smooth_sweeps_))
!!$
!!$        end if
!!$      end if
!!$
!!$    else
!!$      write(iout,*) 'No Base preconditioner available, something is wrong!'
!!$      return
!!$    endif

  end function mld_zprec_short_descr


  !
  ! Functions/subroutines checking if the preconditioner is correctly defined
  !

  function is_legal_base_prec(ip)
    use psb_base_mod
    integer, intent(in) :: ip
    logical             :: is_legal_base_prec

    is_legal_base_prec = ((ip>=mld_noprec_).and.(ip<=mld_max_prec_))
    return
  end function is_legal_base_prec
  function is_legal_n_ovr(ip)
    use psb_base_mod
    integer, intent(in) :: ip
    logical             :: is_legal_n_ovr

    is_legal_n_ovr = (ip >=0) 
    return
  end function is_legal_n_ovr
  function is_legal_renum(ip)
    use psb_base_mod
    integer, intent(in) :: ip
    logical             :: is_legal_renum
    ! For the time being we are disabling renumbering options. 
    is_legal_renum = (ip ==0) 
    return
  end function is_legal_renum
  function is_legal_jac_sweeps(ip)
    use psb_base_mod
    integer, intent(in) :: ip
    logical             :: is_legal_jac_sweeps

    is_legal_jac_sweeps = (ip >= 1) 
    return
  end function is_legal_jac_sweeps
  function is_legal_prolong(ip)
    use psb_base_mod
    integer, intent(in) :: ip
    logical             :: is_legal_prolong
    is_legal_prolong = ((ip>=psb_none_).and.(ip<=psb_square_root_))
    return
  end function is_legal_prolong
  function is_legal_restrict(ip)
    use psb_base_mod
    integer, intent(in) :: ip
    logical             :: is_legal_restrict
    is_legal_restrict = ((ip==psb_nohalo_).or.(ip==psb_halo_))
    return
  end function is_legal_restrict
  function is_legal_ml_type(ip)
    use psb_base_mod
    integer, intent(in) :: ip
    logical             :: is_legal_ml_type

    is_legal_ml_type = ((ip>=mld_no_ml_).and.(ip<=mld_max_ml_))
    return
  end function is_legal_ml_type
  function is_legal_ml_aggr_kind(ip)
    use psb_base_mod
    integer, intent(in) :: ip
    logical             :: is_legal_ml_aggr_kind

    is_legal_ml_aggr_kind = ((ip>=mld_dec_aggr_).and.(ip<=mld_max_aggr_))
    return
  end function is_legal_ml_aggr_kind
  function is_legal_ml_smooth_pos(ip)
    use psb_base_mod
    integer, intent(in) :: ip
    logical             :: is_legal_ml_smooth_pos

    is_legal_ml_smooth_pos = ((ip>=mld_pre_smooth_).and.(ip<=mld_max_smooth_))
    return
  end function is_legal_ml_smooth_pos
  function is_legal_ml_smth_kind(ip)
    use psb_base_mod
    integer, intent(in) :: ip
    logical             :: is_legal_ml_smth_kind

    is_legal_ml_smth_kind = ((ip>=mld_no_smooth_).and.(ip<=mld_biz_prol_))
    return
  end function is_legal_ml_smth_kind
  function is_legal_ml_coarse_mat(ip)
    use psb_base_mod
    integer, intent(in) :: ip
    logical             :: is_legal_ml_coarse_mat

    is_legal_ml_coarse_mat = ((ip>=mld_distr_mat_).and.(ip<=mld_repl_mat_))
    return
  end function is_legal_ml_coarse_mat
  function is_distr_ml_coarse_mat(ip)
    use psb_base_mod
    integer, intent(in) :: ip
    logical             :: is_distr_ml_coarse_mat

    is_distr_ml_coarse_mat = (ip==mld_distr_mat_)
    return
  end function is_distr_ml_coarse_mat
  function is_legal_ml_fact(ip)
    use psb_base_mod
    integer, intent(in) :: ip
    logical             :: is_legal_ml_fact

    is_legal_ml_fact = ((ip>=mld_ilu_n_).and.(ip<=mld_sludist_))
    return
  end function is_legal_ml_fact
  function is_legal_ml_lev(ip)
    use psb_base_mod
    integer, intent(in) :: ip
    logical             :: is_legal_ml_lev

    is_legal_ml_lev = (ip>=0)
    return
  end function is_legal_ml_lev
  function is_legal_omega(ip)
    use psb_base_mod
    real(kind(1.d0)), intent(in) :: ip
    logical             :: is_legal_omega

    is_legal_omega = ((ip>=0.0d0).and.(ip<=2.0d0))
    return
  end function is_legal_omega
  function is_legal_fact_thrs(ip)
    use psb_base_mod
    real(kind(1.d0)), intent(in) :: ip
    logical             :: is_legal_fact_thrs

    is_legal_fact_thrs = (ip>=0.0d0)
    return
  end function is_legal_fact_thrs


  subroutine mld_icheck_def(ip,name,id,is_legal)
    use psb_base_mod
    integer, intent(inout) :: ip
    integer, intent(in)    :: id
    character(len=*), intent(in) :: name
    interface 
      function is_legal(i)
        integer, intent(in) :: i
        logical             :: is_legal
      end function is_legal
    end interface

    if (.not.is_legal(ip)) then     
      write(0,*) 'Illegal value for ',name,' :',ip, '. defaulting to ',id
      ip = id
    end if
  end subroutine mld_icheck_def

  subroutine mld_dcheck_def(ip,name,id,is_legal)
    use psb_base_mod
    real(kind(1.d0)), intent(inout) :: ip
    real(kind(1.d0)), intent(in)    :: id
    character(len=*), intent(in) :: name
    interface 
      function is_legal(i)
        real(kind(1.d0)), intent(in) :: i
        logical             :: is_legal
      end function is_legal
    end interface

    if (.not.is_legal(ip)) then     
      write(0,*) 'Illegal value for ',name,' :',ip, '. defaulting to ',id
      ip = id
    end if
  end subroutine mld_dcheck_def

  subroutine mld_dbase_precfree(p,info)
    use psb_base_mod

    type(mld_dbaseprc_type), intent(inout) :: p
    integer, intent(out)                :: info
    integer :: i

    info = 0

    ! Actually we might just deallocate the top level array, except 
    ! for the inner UMFPACK or SLU stuff

    if (allocated(p%d)) then 
      deallocate(p%d,stat=info)
    end if

    if (allocated(p%av))  then 
      do i=1,size(p%av) 
        call psb_sp_free(p%av(i),info)
        if (info /= 0) then 
          ! Actually, we don't care here about this.
          ! Just let it go.
          ! return
        end if
      enddo
      deallocate(p%av,stat=info)
    end if

    if (allocated(p%desc_data%matrix_data)) &
         & call psb_cdfree(p%desc_data,info)
    if (allocated(p%desc_ac%matrix_data)) &
         & call psb_cdfree(p%desc_ac,info)
    
    if (allocated(p%dprcparm)) then 
      deallocate(p%dprcparm,stat=info)
    end if
    ! This is a pointer to something else, must not free it here. 
    nullify(p%base_a) 
    ! This is a pointer to something else, must not free it here. 
    nullify(p%base_desc) 

    if (allocated(p%dorig)) then 
      deallocate(p%dorig,stat=info)
    endif

    if (allocated(p%mlia)) then 
      deallocate(p%mlia,stat=info)
    endif

    if (allocated(p%nlaggr)) then 
      deallocate(p%nlaggr,stat=info)
    endif

    if (allocated(p%perm)) then 
      deallocate(p%perm,stat=info)
    endif

    if (allocated(p%invperm)) then 
      deallocate(p%invperm,stat=info)
    endif

    if (allocated(p%iprcparm)) then 
      if (p%iprcparm(mld_sub_solve_)==mld_slu_) then 
        call mld_dslu_free(p%iprcparm(mld_slu_ptr_),info)
      end if
      if (p%iprcparm(mld_sub_solve_)==mld_sludist_) then 
        call mld_dsludist_free(p%iprcparm(mld_slud_ptr_),info)
      end if
      if (p%iprcparm(mld_sub_solve_)==mld_umf_) then 
        call mld_dumf_free(p%iprcparm(mld_umf_symptr_),&
             & p%iprcparm(mld_umf_numptr_),info)
      end if
      deallocate(p%iprcparm,stat=info)
    end if
    call mld_nullify_baseprec(p)
  end subroutine mld_dbase_precfree

  subroutine mld_nullify_dbaseprec(p)
    use psb_base_mod

    type(mld_dbaseprc_type), intent(inout) :: p

    nullify(p%base_a) 
    nullify(p%base_desc) 
!!$    nullify(p%av,p%d,p%iprcparm,p%dprcparm,p%perm,p%invperm,p%mlia,&
!!$         & p%nlaggr,p%base_a,p%base_desc,p%dorig,p%desc_data, p%desc_ac)

  end subroutine mld_nullify_dbaseprec

  subroutine mld_zbase_precfree(p,info)
    use psb_base_mod
    type(mld_zbaseprc_type), intent(inout) :: p
    integer, intent(out)                :: info
    integer :: i

    info = 0

    if (allocated(p%d)) then 
      deallocate(p%d,stat=info)
    end if

    if (allocated(p%av))  then 
      do i=1,size(p%av) 
        call psb_sp_free(p%av(i),info)
        if (info /= 0) then 
          ! Actually, we don't care here about this.
          ! Just let it go.
          ! return
        end if
      enddo
      deallocate(p%av,stat=info)

    end if
    ! call psb_cdfree(p%desc_data,info)
    ! call psb_cdfree(p%desc_ac,info)

    if (allocated(p%dprcparm)) then 
      deallocate(p%dprcparm,stat=info)
    end if
    ! This is a pointer to something else, must not free it here. 
    nullify(p%base_a) 
    ! This is a pointer to something else, must not free it here. 
    nullify(p%base_desc) 

    if (allocated(p%dorig)) then 
      deallocate(p%dorig,stat=info)
    endif

    if (allocated(p%mlia)) then 
      deallocate(p%mlia,stat=info)
    endif

    if (allocated(p%nlaggr)) then 
      deallocate(p%nlaggr,stat=info)
    endif

    if (allocated(p%perm)) then 
      deallocate(p%perm,stat=info)
    endif

    if (allocated(p%invperm)) then 
      deallocate(p%invperm,stat=info)
    endif

    if (allocated(p%iprcparm)) then 
      if (p%iprcparm(mld_sub_solve_)==mld_slu_) then 
        call mld_zslu_free(p%iprcparm(mld_slu_ptr_),info)
      end if
      if (p%iprcparm(mld_sub_solve_)==mld_umf_) then 
        call mld_zumf_free(p%iprcparm(mld_umf_symptr_),&
             & p%iprcparm(mld_umf_numptr_),info)
      end if
      deallocate(p%iprcparm,stat=info)
    end if
    call mld_nullify_baseprec(p)
  end subroutine mld_zbase_precfree

  subroutine mld_nullify_zbaseprec(p)
    use psb_base_mod

    type(mld_zbaseprc_type), intent(inout) :: p


    nullify(p%base_a) 
    nullify(p%base_desc) 

  end subroutine mld_nullify_zbaseprec


  function pr_to_str(iprec)
    use psb_base_mod

    integer, intent(in)  :: iprec
    character(len=10)     :: pr_to_str

    select case(iprec)
    case(mld_noprec_)
      pr_to_str='NOPREC'
    case(mld_diag_)         
      pr_to_str='DIAG'
    case(mld_bjac_)         
      pr_to_str='BJAC'
    case(mld_as_)      
      pr_to_str='AS'
    end select

  end function pr_to_str

end module mld_prec_type
