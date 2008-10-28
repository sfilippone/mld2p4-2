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
! File: mld_cmlprec_aply.f90
!
! Subroutine: mld_cmlprec_aply
! Version:    complex
!
!  This routine computes
!  
!                        Y = beta*Y + alpha*op(M^(-1))*X,
!  where 
!  - M is a multilevel domain decomposition (Schwarz) preconditioner associated
!    to a certain matrix A and stored in the array precv,
!  - op(M^(-1)) is M^(-1) or its transpose, according to the value of trans,
!  - X and Y are vectors,
!  - alpha and beta are scalars.
!
!  For each level we have as many submatrices as processes (except for the coarsest
!  level where we might have a replicated index space) and each process takes care
!  of one submatrix.
!
!  The multilevel preconditioner M is regarded as an array of 'one-level preconditioners',
!  each representing the part of the preconditioner associated to a certain level.
!  For each level ilev, the  preconditioner K(ilev) is stored in precv(ilev)
!  and is associated to a matrix A(ilev), obtained by 'tranferring' the original
!  matrix A (i.e. the matrix to be preconditioned) to the level ilev, through smoothed
!  aggregation.
!
!  The levels are numbered in increasing order starting from the finest one, i.e.
!  level 1 is the finest level and A(1) is the matrix A.
!
!  For a general description of (parallel) multilevel preconditioners see
!    -  B.F. Smith, P.E. Bjorstad & W.D. Gropp,
!       Domain decomposition: parallel multilevel methods for elliptic partial
!       differential equations,
!       Cambridge University Press, 1996.
!    -  K. Stuben,
!       Algebraic Multigrid (AMG): An Introduction with Applications,
!       GMD Report N. 70, 1999.
!
!
! Arguments:
!  alpha       -   complex(psb_spk_), input.
!                  The scalar alpha.
!   precv      -   type(mld_c_onelev_prec_type), dimension(:), input.
!                  The array of one-level preconditioner data structures containing the
!                  local parts of the preconditioners to be applied at each level.
!      Note that nlev = size(precv) = number of levels.
!      precv(ilev)%prec  -  type(psb_cbaseprc_type)
!                           The "base" preconditioner for the current level
!      precv(ilev)%ac        -  type(psb_cspmat_type) 
!                                The local part of the matrix A(ilev).
!      precv(ilev)%desc_ac   -  type(psb_desc_type).
!                               The communication descriptor associated to the sparse
!                               matrix A(ilev)
!      precv(ilev)%map_desc  -  type(psb_inter_desc_type)
!                               Stores the linear operators mapping between levels
!                               (ilev-1) and (ilev). These are the restriction and
!                               prolongation operators described in the sequel. 
!      precv(ilev)%iprcparm  -  integer, dimension(:), allocatable.
!                               The integer parameters defining the multilevel
!                               strategy
!      precv(ilev)%rprcparm  -  real(psb_spk_), dimension(:), allocatable.
!                               The real parameters defining the multilevel strategy
!      precv(ilev)%mlia      -  integer, dimension(:), allocatable.
!                               The aggregation map (ilev-1) --> (ilev).
!                               In case of non-smoothed aggregation, it is used
!                               instead of mld_sm_pr_.
!      precv(ilev)%nlaggr    -  integer, dimension(:), allocatable.
!                               The number of aggregates (rows of A(ilev)) on the
!                               various processes. 
!      precv(ilev)%base_a    -  type(psb_cspmat_type), pointer.
!                               Pointer (really a pointer!) to the base matrix of
!                               the current level, i.e. the local part of A(ilev);
!                               so we have a unified treatment of residuals. We
!                               need this to avoid passing explicitly the matrix
!                               A(ilev) to the routine which applies the
!                               preconditioner.
!      precv(ilev)%base_desc -  type(psb_desc_type), pointer.
!                               Pointer to the communication descriptor associated
!                               to the sparse matrix pointed by base_a.  
!                  
!   x          -  complex(psb_spk_), dimension(:), input.
!                 The local part of the vector X.
!   beta       -  complex(psb_spk_), input.
!                 The scalar beta.
!   y          -  complex(psb_spk_), dimension(:), input/output.
!                 The local part of the vector Y.
!   desc_data  -  type(psb_desc_type), input.
!                 The communication descriptor associated to the matrix to be
!                 preconditioned.
!   trans      -  character, optional.
!                 If trans='N','n' then op(M^(-1)) = M^(-1);
!                 if trans='T','t' then op(M^(-1)) = M^(-T) (transpose of M^(-1)).
!   work       -  complex(psb_spk_), dimension (:), optional, target.
!                 Workspace. Its size must be at least 4*psb_cd_get_local_cols(desc_data).
!   info       -  integer, output.
!                 Error code.
!
!   Note that when the LU factorization of the matrix A(ilev) is computed instead of
!   the ILU one, by using UMFPACK or SuperLU, the corresponding L and U factors
!   are stored in data structures provided by UMFPACK or SuperLU and pointed by
!   precv(ilev)%prec%iprcparm(mld_umf_ptr) or precv(ilev)%prec%iprcparm(mld_slu_ptr),
!   respectively.
!  
subroutine mld_cmlprec_aply(alpha,precv,x,beta,y,desc_data,trans,work,info)

  use psb_base_mod
  use mld_inner_mod, mld_protect_name => mld_cmlprec_aply

  implicit none

  ! Arguments
  type(psb_desc_type),intent(in)      :: desc_data
  type(mld_c_onelev_prec_type), intent(in) :: precv(:)
  complex(psb_spk_),intent(in)      :: alpha,beta
  complex(psb_spk_),intent(in)      :: x(:)
  complex(psb_spk_),intent(inout)   :: y(:)
  character, intent(in)               :: trans
  complex(psb_spk_),target          :: work(:)
  integer, intent(out)                :: info

  ! Local variables
  integer           :: ictxt, np, me, err_act
  integer           :: debug_level, debug_unit
  character(len=20) :: name
  character         :: trans_

  name = 'mld_cmlprec_aply'
  info = 0
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = psb_cd_get_context(desc_data)
  call psb_info(ictxt, me, np)

  if (debug_level >= psb_debug_inner_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ' Entry  ', size(precv)

  trans_ = psb_toupper(trans)

  select case(precv(2)%iprcparm(mld_ml_type_)) 

  case(mld_no_ml_)
    !
    ! No preconditioning, should not really get here
    ! 
    call psb_errpush(4001,name,a_err='mld_no_ml_ in mlprc_aply?')
    goto 9999      

  case(mld_add_ml_)
    !
    ! Additive multilevel
    !

    call add_ml_aply(alpha,precv,x,beta,y,desc_data,trans_,work,info)

  case(mld_mult_ml_)
    ! 
    !  Multiplicative multilevel (multiplicative among the levels, additive inside
    !  each level)
    !
    !  Pre/post-smoothing versions.
    !  Note that the transpose switches pre <-> post.
    !

    select case(precv(2)%iprcparm(mld_smoother_pos_))

    case(mld_post_smooth_)

      select case (trans_) 
      case('N')
        call mlt_post_ml_aply(alpha,precv,x,beta,y,desc_data,trans_,work,info)
      case('T','C')
        call mlt_pre_ml_aply(alpha,precv,x,beta,y,desc_data,trans_,work,info)
      case default
        info = 4001
        call psb_errpush(info,name,a_err='invalid trans')
        goto 9999      
      end select

    case(mld_pre_smooth_)

      select case (trans_) 
      case('N')
        call mlt_pre_ml_aply(alpha,precv,x,beta,y,desc_data,trans_,work,info)
      case('T','C')
        call mlt_post_ml_aply(alpha,precv,x,beta,y,desc_data,trans_,work,info)
      case default
        info = 4001
        call psb_errpush(info,name,a_err='invalid trans')
        goto 9999      
      end select

    case(mld_twoside_smooth_)

      call mlt_twoside_ml_aply(alpha,precv,x,beta,y,desc_data,trans_,work,info)

    case default
      info = 4013
      call psb_errpush(info,name,a_err='invalid smooth_pos',&
           &  i_Err=(/precv(2)%iprcparm(mld_smoother_pos_),0,0,0,0/))
      goto 9999      

    end select

  case default
    info = 4013
    call psb_errpush(info,name,a_err='invalid mltype',&
         &  i_Err=(/precv(2)%iprcparm(mld_ml_type_),0,0,0,0/))
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

contains
  !
  ! Subroutine: add_ml_aply
  ! Version:    complex
  ! Note: internal subroutine of mld_dmlprec_aply.
  ! 
  ! This routine computes
  !  
  !                        Y = beta*Y + alpha*op(M^(-1))*X,
  !  where 
  !  - M is an additive multilevel domain decomposition (Schwarz) preconditioner
  !    associated to a certain matrix A and stored in the array precv,
  !  - op(M^(-1)) is M^(-1) or its (conjugate) transpose, according to
  !    the value of trans,
  !  - X and Y are vectors,
  !  - alpha and beta are scalars.
  !
  !  The preconditioner M is additive both through the levels and inside each
  !  level.
  !
  !  For each level we have as many submatrices as processes (except for the coarsest
  !  level where we might have a replicated index space) and each process takes care
  !  of one submatrix. 
  !
  !  The multilevel preconditioner M is regarded as an array of 'one-level preconditioners',
  !  each representing the part of the preconditioner associated to a certain level.
  !  For each level ilev, the base preconditioner K(ilev) is stored in precv(ilev)
  !  and is associated to a matrix A(ilev), obtained by 'tranferring' the original
  !  matrix A (i.e. the matrix to be preconditioned) to the level ilev, through smoothed
  !  aggregation.
  !
  !  The levels are numbered in increasing order starting from the finest one, i.e.
  !  level 1 is the finest level and A(1) is the matrix A. 
  !
  !  For details on the additive multilevel Schwarz preconditioner see the
  !  Algorithm 3.1.1 in the book:
  !    B.F. Smith, P.E. Bjorstad & W.D. Gropp,
  !    Domain decomposition: parallel multilevel methods for elliptic partial
  !    differential equations, Cambridge University Press, 1996.
  !
  !  For a description of the arguments see mld_cmlprec_aply.
  !
  !  A sketch of the algorithm implemented in this routine is provided below
  !  (AV(ilev; sm_pr_) denotes the smoothed prolongator from level ilev to
  !  level ilev-1, while AV(ilev; sm_pr_t_) denotes its transpose, i.e. the
  !  corresponding restriction operator from level ilev-1 to level ilev).
  !
  !   1. ! Apply the base preconditioner at level 1.
  !      ! The sum over the subdomains is carried out in the
  !      ! application of K(1).
  !        X(1) = Xest
  !        Y(1) = (K(1)^(-1))*X(1)
  !
  !    2.  DO ilev=2,nlev
  !
  !         ! Transfer X(ilev-1) to the next coarser level.
  !           X(ilev) = AV(ilev; sm_pr_t_)*X(ilev-1)
  !
  !         ! Apply the base preconditioner at the current level.
  !         ! The sum over the subdomains is carried out in the
  !         ! application of K(ilev).
  !           Y(ilev) = (K(ilev)^(-1))*X(ilev)
  !
  !        ENDDO
  !
  !    3.  DO ilev=nlev-1,1,-1
  !
  !         ! Transfer Y(ilev+1) to the next finer level.
  !           Y(ilev) = AV(ilev+1; sm_pr_)*Y(ilev+1)
  !
  !        ENDDO
  !     
  !    4.  Yext = beta*Yext + alpha*Y(1)
  !
  subroutine add_ml_aply(alpha,precv,x,beta,y,desc_data,trans,work,info)

    implicit none 

    ! Arguments
    type(psb_desc_type),intent(in)      :: desc_data
    type(mld_c_onelev_prec_type), intent(in) :: precv(:)
    complex(psb_spk_),intent(in)      :: alpha,beta
    complex(psb_spk_),intent(in)      :: x(:)
    complex(psb_spk_),intent(inout)   :: y(:)
    character, intent(in)               :: trans
    complex(psb_spk_),target          :: work(:)
    integer, intent(out)                :: info

    ! Local variables
    integer            :: ictxt,np,me,i, nr2l,nc2l,err_act
    integer            :: debug_level, debug_unit
    integer            :: nlev, ilev
    character(len=20)  :: name

    type psb_mlprec_wrk_type
      complex(psb_spk_), allocatable  :: tx(:),ty(:),x2l(:),y2l(:)
    end type psb_mlprec_wrk_type
    type(psb_mlprec_wrk_type), allocatable :: mlprec_wrk(:)

    name = 'add_ml_aply'
    info = 0
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    ictxt = psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)

    if (debug_level >= psb_debug_inner_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' Entry  ', size(precv)

    nlev = size(precv)
    allocate(mlprec_wrk(nlev),stat=info) 
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    !
    ! STEP 1
    !
    ! Apply the base preconditioner at the finest level
    !
    allocate(mlprec_wrk(1)%x2l(size(x)),mlprec_wrk(1)%y2l(size(y)), stat=info)
    if (info /= 0) then 
      info=4025
      call psb_errpush(info,name,i_err=(/size(x)+size(y),0,0,0,0/),&
           & a_err='complex(psb_spk_)')
      goto 9999      
    end if

    mlprec_wrk(1)%x2l(:) = x(:) 
    mlprec_wrk(1)%y2l(:) = czero 

    call mld_baseprec_aply(alpha,precv(1)%prec,x,beta,y,&
         & precv(1)%base_desc,trans,work,info)
    if (info /=0) then 
      call psb_errpush(4010,name,a_err='baseprec_aply')
      goto 9999
    end if
    !
    ! STEP 2
    !
    ! For each level except the finest one ...
    !
    do ilev = 2, nlev
      nc2l  = psb_cd_get_local_cols(precv(ilev)%base_desc)
      nr2l  = psb_cd_get_local_rows(precv(ilev)%base_desc)
      allocate(mlprec_wrk(ilev)%x2l(nc2l),mlprec_wrk(ilev)%y2l(nc2l),&
           & stat=info)
      if (info /= 0) then 
        info=4025
        call psb_errpush(info,name,i_err=(/2*nc2l,0,0,0,0/),&
             & a_err='complex(psb_spk_)')
        goto 9999      
      end if
      
      ! Apply prolongator transpose, i.e. restriction
      call psb_forward_map(cone,mlprec_wrk(ilev-1)%x2l,&
           & czero,mlprec_wrk(ilev)%x2l,&
           & precv(ilev)%map_desc,info,work=work)
      
      if (info /=0) then
        call psb_errpush(4001,name,a_err='Error during restriction')
        goto 9999
      end if

      !
      ! Apply the base preconditioner
      !
      call mld_baseprec_aply(cone,precv(ilev)%prec,&
           & mlprec_wrk(ilev)%x2l,czero,mlprec_wrk(ilev)%y2l,&
           & precv(ilev)%base_desc,trans,work,info)

    enddo

    !
    ! STEP 3
    !
    ! For each level except the finest one ...
    !
    do ilev =nlev,2,-1

      nc2l  = psb_cd_get_local_cols(precv(ilev)%base_desc)
      nr2l  = psb_cd_get_local_rows(precv(ilev)%base_desc)

      !
      ! Apply prolongator
      !  
      call psb_backward_map(cone,mlprec_wrk(ilev)%y2l,&
           & cone,mlprec_wrk(ilev-1)%y2l,&
           & precv(ilev)%map_desc,info,work=work)

      if (info /=0) then
        call psb_errpush(4001,name,a_err='Error during prolongation')
        goto 9999
      end if
    end do

    !
    ! STEP 4
    !
    ! Compute the output vector Y
    !
    call psb_geaxpby(alpha,mlprec_wrk(1)%y2l,cone,y,precv(1)%base_desc,info)
    if (info /= 0) then
      call psb_errpush(4001,name,a_err='Error on final update')
      goto 9999
    end if

    deallocate(mlprec_wrk,stat=info)
    if (info /= 0) then 
      call psb_errpush(4000,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine add_ml_aply
  !
  ! Subroutine: mlt_pre_ml_aply
  ! Version:    complex
  ! Note: internal subroutine of mld_cmlprec_aply.
  ! 
  !  This routine computes
  !  
  !                        Y = beta*Y + alpha*op(M^(-1))*X,
  !  where 
  !  - M is a hybrid multilevel domain decomposition (Schwarz) preconditioner
  !    associated to a certain matrix A and stored in the array precv,
  !  - op(M^(-1)) is M^(-1) or its (conjugate) transpose, according to
  !    the value of trans,
  !  - X and Y are vectors,
  !  - alpha and beta are scalars.
  !
  !  The preconditioner M is hybrid in the sense that it is multiplicative through the
  !  levels and additive inside a level; pre-smoothing only is applied at each level.
  !
  !  For each level we have as many submatrices as processes (except for the coarsest
  !  level where we might have a replicated index space) and each process takes care
  !  of one submatrix. 
  !
  !  The multilevel preconditioner M is regarded as an array of 'one-level preconditioners',
  !  each representing the part of the preconditioner associated to a certain level.
  !  For each level ilev, the base preconditioner K(ilev) is stored in precv(ilev)
  !  and is associated to a matrix A(ilev), obtained by 'tranferring' the original
  !  matrix A (i.e. the matrix to be preconditioned) to the level ilev, through smoothed
  !  aggregation.
  !
  !  The levels are numbered in increasing order starting from the finest one, i.e.
  !  level 1 is the finest level and A(1) is the matrix A. 
  !
  !  For details on the pre-smoothed hybrid multiplicative multilevel Schwarz
  !  preconditioner, see the Algorithm 3.2.1 in the book:
  !    B.F. Smith, P.E. Bjorstad & W.D. Gropp,
  !    Domain decomposition: parallel multilevel methods for elliptic partial
  !    differential equations, Cambridge University Press, 1996.
  !
  !  For a description of the arguments see mld_cmlprec_aply.
  ! 
  !  A sketch of the algorithm implemented in this routine is provided below
  !  (AV(ilev; sm_pr_) denotes the smoothed prolongator from level ilev to
  !  level ilev-1, while AV(ilev; sm_pr_t_) denotes its transpose, i.e. the
  !  corresponding restriction operator from level ilev-1 to level ilev).
  !
  !    1.   X(1) = Xext
  !
  !    2. ! Apply the base preconditioner at the finest level.
  !         Y(1) = (K(1)^(-1))*X(1)
  !
  !    3. ! Compute the residual at the finest level.
  !         TX(1) = X(1) - A(1)*Y(1)
  !
  !    4.   DO ilev=2, nlev
  !
  !          ! Transfer the residual to the current (coarser) level.
  !            X(ilev) = AV(ilev; sm_pr_t_)*TX(ilev-1)
  !
  !          ! Apply the base preconditioner at the current level.
  !          ! The sum over the subdomains is carried out in the
  !          ! application of K(ilev).     
  !            Y(ilev) = (K(ilev)^(-1))*X(ilev)
  !
  !          ! Compute the residual at the current level (except at
  !          ! the coarsest level).
  !            IF (ilev < nlev)
  !               TX(ilev) = (X(ilev)-A(ilev)*Y(ilev))
  !
  !         ENDDO
  !
  !    5.   DO ilev=nlev-1,1,-1
  !
  !          ! Transfer Y(ilev+1) to the next finer level
  !            Y(ilev) = Y(ilev) + AV(ilev+1; sm_pr_)*Y(ilev+1)
  !
  !         ENDDO
  !     
  !    6.  Yext = beta*Yext + alpha*Y(1)
  ! 
  !
  subroutine mlt_pre_ml_aply(alpha,precv,x,beta,y,desc_data,trans,work,info)

    implicit none 

    ! Arguments
    type(psb_desc_type),intent(in)      :: desc_data
    type(mld_c_onelev_prec_type), intent(in) :: precv(:)
    complex(psb_spk_),intent(in)      :: alpha,beta
    complex(psb_spk_),intent(in)      :: x(:)
    complex(psb_spk_),intent(inout)   :: y(:)
    character, intent(in)               :: trans
    complex(psb_spk_),target          :: work(:)
    integer, intent(out)                :: info

    ! Local variables
    integer            :: ictxt,np,me,i, nr2l,nc2l,err_act
    integer            :: debug_level, debug_unit
    integer            :: nlev, ilev
    character(len=20)  :: name

    type psb_mlprec_wrk_type
      complex(psb_spk_), allocatable  :: tx(:),ty(:),x2l(:),y2l(:)
    end type psb_mlprec_wrk_type
    type(psb_mlprec_wrk_type), allocatable :: mlprec_wrk(:)

    name = 'mlt_pre_ml_aply'
    info = 0
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    ictxt = psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)

    if (debug_level >= psb_debug_inner_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' Entry  ', size(precv)

    nlev = size(precv)
    allocate(mlprec_wrk(nlev),stat=info) 
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    !
    ! STEP 1
    !
    ! Copy the input vector X
    !
    nc2l  = psb_cd_get_local_cols(precv(1)%base_desc)

    allocate(mlprec_wrk(1)%x2l(nc2l),mlprec_wrk(1)%y2l(nc2l), &
         & mlprec_wrk(1)%tx(nc2l), stat=info)
    if (info /= 0) then 
      info=4025
      call psb_errpush(info,name,i_err=(/4*nc2l,0,0,0,0/),&
           & a_err='complex(psb_spk_)')
      goto 9999
    end if

    mlprec_wrk(1)%x2l(:) = x
    !
    ! STEP 2
    !
    ! Apply the base preconditioner at the finest level
    !
    call mld_baseprec_aply(cone,precv(1)%prec,mlprec_wrk(1)%x2l,&
         &  czero,mlprec_wrk(1)%y2l,precv(1)%base_desc,&
         &  trans,work,info)
    if (info /=0) then
      call psb_errpush(4010,name,a_err=' baseprec_aply')
      goto 9999
    end if

    !
    ! STEP 3
    !
    ! Compute the residual at the finest level
    !
    mlprec_wrk(1)%tx = mlprec_wrk(1)%x2l

    call psb_spmm(-cone,precv(1)%base_a,mlprec_wrk(1)%y2l,&
         & cone,mlprec_wrk(1)%tx,precv(1)%base_desc,info,&
         & work=work,trans=trans)
    if (info /=0) then
      call psb_errpush(4001,name,a_err=' fine level residual')
      goto 9999
    end if

    !
    ! STEP 4
    !
    ! For each level but the finest one ...
    !
    do ilev = 2, nlev

      nc2l  = psb_cd_get_local_cols(precv(ilev)%base_desc)
      nr2l  = psb_cd_get_local_rows(precv(ilev)%base_desc)

      allocate(mlprec_wrk(ilev)%tx(nc2l),mlprec_wrk(ilev)%y2l(nc2l),&
           &   mlprec_wrk(ilev)%x2l(nc2l), stat=info)
      if (info /= 0) then 
        info=4025
        call psb_errpush(info,name,i_err=(/4*nc2l,0,0,0,0/),&
             & a_err='complex(psb_spk_)')
        goto 9999      
      end if

      ! Apply prolongator transpose, i.e. restriction      
      call psb_forward_map(cone,mlprec_wrk(ilev-1)%tx,&
           & czero,mlprec_wrk(ilev)%x2l,&
           & precv(ilev)%map_desc,info,work=work)
      
      if (info /=0) then
        call psb_errpush(4001,name,a_err='Error during restriction')
        goto 9999
      end if

      !
      ! Apply the base preconditioner
      !
      call mld_baseprec_aply(cone,precv(ilev)%prec,mlprec_wrk(ilev)%x2l,&
           & czero,mlprec_wrk(ilev)%y2l,precv(ilev)%base_desc,trans,work,info)

      !
      ! Compute the residual (at all levels but the coarsest one)
      !
      if (ilev < nlev) then
        mlprec_wrk(ilev)%tx = mlprec_wrk(ilev)%x2l
        if (info == 0) call psb_spmm(-cone,precv(ilev)%base_a,&
             & mlprec_wrk(ilev)%y2l,cone,mlprec_wrk(ilev)%tx,&
             & precv(ilev)%base_desc,info,work=work,trans=trans)
      endif
      if (info /=0) then
        call psb_errpush(4001,name,a_err='Error on up sweep residual')
        goto 9999
      end if
    enddo

    !
    ! STEP 5
    !
    ! For each level but the coarsest one ...
    !
    do ilev = nlev-1, 1, -1
      !
      ! Apply prolongator
      !  
      call psb_backward_map(cone,mlprec_wrk(ilev+1)%y2l,&
           & cone,mlprec_wrk(ilev)%y2l,&
           & precv(ilev+1)%map_desc,info,work=work)

      if (info /=0) then
        call psb_errpush(4001,name,a_err='Error during prolongation')
        goto 9999
      end if
    enddo

    !
    ! STEP 6
    !
    ! Compute the output vector Y
    !
    call psb_geaxpby(alpha,mlprec_wrk(1)%y2l,beta,y,&
         &  precv(1)%base_desc,info)
    if (info /=0) then
      call psb_errpush(4001,name,a_err='Error on final update')
      goto 9999
    end if

    deallocate(mlprec_wrk,stat=info)
    if (info /= 0) then 
      call psb_errpush(4000,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine mlt_pre_ml_aply
  !
  ! Subroutine: mlt_post_ml_aply
  ! Version:    complex
  ! Note:       internal subroutine of mld_cmlprec_aply.
  ! 
  !  This routine computes
  !  
  !                        Y = beta*Y + alpha*op(M^(-1))*X,
  !  where 
  !  - M is a hybrid multilevel domain decomposition (Schwarz) preconditioner
  !    associated to a certain matrix A and stored in the array precv,
  !  - op(M^(-1)) is M^(-1) or its (conjugate) transpose, according to
  !    the value of trans,
  !  - X and Y are vectors,
  !  - alpha and beta are scalars.
  !
  !  The preconditioner M is hybrid in the sense that it is multiplicative through the
  !  levels and additive inside a level; post-smoothing only is applied at each level.
  !
  !  For each level we have as many submatrices as processes (except for the coarsest
  !  level where we might have a replicated index space) and each process takes care
  !  of one submatrix. 
  !
  !  The multilevel preconditioner M is regarded as an array of 'one-level preconditioners',
  !  each representing the part of the preconditioner associated to a certain level.
  !  For each level ilev, the base preconditioner K(ilev) is stored in precv(ilev)
  !  and is associated to a matrix A(ilev), obtained by 'tranferring' the original
  !  matrix A (i.e. the matrix to be preconditioned) to the level ilev, through smoothed
  !  aggregation.
  !
  !  The levels are numbered in increasing order starting from the finest one, i.e.
  !  level 1 is the finest level and A(1) is the matrix A. 
  !
  !  For details on hybrid multiplicative multilevel Schwarz preconditioners, see
  !    B.F. Smith, P.E. Bjorstad & W.D. Gropp,
  !    Domain decomposition: parallel multilevel methods for elliptic partial
  !    differential equations, Cambridge University Press, 1996.
  !
  !  For a description of the arguments see mld_cmlprec_aply.
  !
  !  A sketch of the algorithm implemented in this routine is provided below.
  !  (AV(ilev; sm_pr_) denotes the smoothed prolongator from level ilev to
  !  level ilev-1, while AV(ilev; sm_pr_t_) denotes its transpose, i.e. the
  !  corresponding restriction operator from level ilev-1 to level ilev).
  !
  !    1.  X(1) = Xext
  !
  !    2.  DO ilev=2, nlev
  !
  !         ! Transfer X(ilev-1) to the next coarser level.
  !           X(ilev) = AV(ilev; sm_pr_t_)*X(ilev-1) 
  !
  !        ENDDO
  !   
  !    3.! Apply the preconditioner at the coarsest level.
  !        Y(nlev) = (K(nlev)^(-1))*X(nlev)
  !
  !    4.  DO ilev=nlev-1,1,-1
  !
  !         ! Transfer Y(ilev+1) to the next finer level.
  !           Y(ilev) = AV(ilev+1; sm_pr_)*Y(ilev+1)
  !
  !         ! Compute the residual at the current level and apply to it the
  !         ! base preconditioner. The sum over the subdomains is carried out
  !         ! in the application of K(ilev).
  !           Y(ilev) = Y(ilev) + (K(ilev)^(-1))*(X(ilev)-A(ilev)*Y(ilev))
  !
  !        ENDDO
  !
  !    5.  Yext = beta*Yext + alpha*Y(1)
  ! 
  !
  subroutine mlt_post_ml_aply(alpha,precv,x,beta,y,desc_data,trans,work,info)

    implicit none 

    ! Arguments
    type(psb_desc_type),intent(in)      :: desc_data
    type(mld_c_onelev_prec_type), intent(in) :: precv(:)
    complex(psb_spk_),intent(in)      :: alpha,beta
    complex(psb_spk_),intent(in)      :: x(:)
    complex(psb_spk_),intent(inout)   :: y(:)
    character, intent(in)               :: trans
    complex(psb_spk_),target          :: work(:)
    integer, intent(out)                :: info

    ! Local variables
    integer            :: ictxt,np,me,i, nr2l,nc2l,err_act
    integer            :: debug_level, debug_unit
    integer            :: nlev, ilev
    character(len=20)  :: name

    type psb_mlprec_wrk_type
      complex(psb_spk_), allocatable  :: tx(:),ty(:),x2l(:),y2l(:)
    end type psb_mlprec_wrk_type
    type(psb_mlprec_wrk_type), allocatable :: mlprec_wrk(:)

    name = 'mlt_post_ml_aply'
    info = 0
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    ictxt = psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)

    if (debug_level >= psb_debug_inner_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' Entry  ', size(precv)

    nlev = size(precv)
    allocate(mlprec_wrk(nlev),stat=info) 
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    !
    ! STEP 1
    !
    ! Copy the input vector X
    !
    if (debug_level >= psb_debug_inner_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' desc_data status',allocated(desc_data%matrix_data)    

    nc2l  = psb_cd_get_local_cols(precv(1)%base_desc)

    allocate(mlprec_wrk(1)%x2l(nc2l),mlprec_wrk(1)%y2l(nc2l), &
         & mlprec_wrk(1)%tx(nc2l), stat=info)
    if (info /= 0) then 
      info=4025
      call psb_errpush(info,name,i_err=(/4*nc2l,0,0,0,0/),&
           & a_err='complex(psb_spk_)')
      goto 9999
    end if

    call psb_geaxpby(cone,x,czero,mlprec_wrk(1)%tx,&
         & precv(1)%base_desc,info)
    call psb_geaxpby(cone,x,czero,mlprec_wrk(1)%x2l,&
         & precv(1)%base_desc,info)

    !
    ! STEP 2
    !
    ! For each level but the finest one ...
    !
    do ilev=2, nlev

      nc2l  = psb_cd_get_local_cols(precv(ilev)%base_desc)
      nr2l  = psb_cd_get_local_rows(precv(ilev)%base_desc)

      if (debug_level >= psb_debug_inner_) &
           & write(debug_unit,*) me,' ',trim(name), &
           & ' starting up sweep ',&
           & ilev,allocated(precv(ilev)%iprcparm),nc2l, nr2l

      allocate(mlprec_wrk(ilev)%tx(nc2l),mlprec_wrk(ilev)%y2l(nc2l),&
           &   mlprec_wrk(ilev)%x2l(nc2l), stat=info)

      if (info /= 0) then 
        info=4025
        call psb_errpush(info,name,i_err=(/4*nc2l,0,0,0,0/),&
             & a_err='complex(psb_spk_)')
        goto 9999      
      end if

      ! Apply prolongator transpose, i.e. restriction
      call psb_forward_map(cone,mlprec_wrk(ilev-1)%x2l,&
           & czero,mlprec_wrk(ilev)%x2l,&
           & precv(ilev)%map_desc,info,work=work)
      
      if (info /=0) then
        call psb_errpush(4001,name,a_err='Error during restriction')
        goto 9999
      end if

      !
      ! update x2l
      !
      call psb_geaxpby(cone,mlprec_wrk(ilev)%x2l,czero,mlprec_wrk(ilev)%tx,&
           & precv(ilev)%base_desc,info)
      if (info /= 0) then
        call psb_errpush(4001,name,a_err='Error in update')
        goto 9999
      end if

      if (debug_level >= psb_debug_inner_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & ' done up sweep ', ilev

    enddo

    !
    ! STEP 3
    !
    ! Apply the base preconditioner at the coarsest level
    !
    call mld_baseprec_aply(cone,precv(nlev)%prec,mlprec_wrk(nlev)%x2l, &
         & czero, mlprec_wrk(nlev)%y2l,precv(nlev)%base_desc,trans,work,info)
    if (info /=0) then
      call psb_errpush(4010,name,a_err='baseprec_aply')
      goto 9999
    end if

    if (debug_level >= psb_debug_inner_) write(debug_unit,*) &
         & me,' ',trim(name), ' done baseprec_aply ', nlev

    !
    ! STEP 4
    !
    ! For each level but the coarsest one ...
    !
    do ilev=nlev-1, 1, -1

      if (debug_level >= psb_debug_inner_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & ' starting down sweep',ilev

      !
      ! Apply prolongator
      !  
      call psb_backward_map(cone,mlprec_wrk(ilev+1)%y2l,&
           & czero,mlprec_wrk(ilev)%y2l,&
           & precv(ilev+1)%map_desc,info,work=work)

      if (info /=0) then
        call psb_errpush(4001,name,a_err='Error during prolongation')
        goto 9999
      end if

      !
      ! Compute the residual
      !
      call psb_spmm(-cone,precv(ilev)%base_a,mlprec_wrk(ilev)%y2l,&
           & cone,mlprec_wrk(ilev)%tx,precv(ilev)%base_desc,info,&
           & work=work,trans=trans)

      !
      ! Apply the base preconditioner
      !
      if (info == 0) call mld_baseprec_aply(cone,precv(ilev)%prec,&
           & mlprec_wrk(ilev)%tx,cone,mlprec_wrk(ilev)%y2l,precv(ilev)%base_desc,&
           & trans,work,info)
      if (info /=0) then
        call psb_errpush(4001,name,a_err=' spmm/baseprec_aply')
        goto 9999
      end if

      if (debug_level >= psb_debug_inner_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & ' done down sweep',ilev
    enddo

    !
    ! STEP 5
    !
    ! Compute the output vector Y
    !
    call psb_geaxpby(alpha,mlprec_wrk(1)%y2l,beta,y,precv(1)%base_desc,info)

    if (info /=0) then
      call psb_errpush(4001,name,a_err=' Final update')
      goto 9999
    end if



    deallocate(mlprec_wrk,stat=info)
    if (info /= 0) then 
      call psb_errpush(4000,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine mlt_post_ml_aply
  !
  ! Subroutine: mlt_twoside_ml_aply
  ! Version:    complex
  ! Note:       internal subroutine of mld_cmlprec_aply.
  ! 
  !  This routine computes
  !  
  !                        Y = beta*Y + alpha*op(M^(-1))*X,
  !  where 
  !  - M is a symmetrized hybrid multilevel domain decomposition (Schwarz)
  !    preconditioner associated to a certain matrix A and stored in the array
  !    precv,
  !  - op(M^(-1)) is M^(-1) or its (conjugate) transpose, according to
  !    the value of trans,
  !  - X and Y are vectors,               
  !  - alpha and beta are scalars.
  !
  !  The preconditioner M is hybrid in the sense that it is multiplicative through
  !  the levels and additive inside a level; it is symmetrized since pre-smoothing
  !  and post-smoothing are applied at each level.
  !
  !  For each level we have as many submatrices as processes (except for the coarsest
  !  level where we might have a replicated index space) and each process takes care
  !  of one submatrix.   
  !
  !  The multilevel preconditioner M is regarded as an array of 'one-level preconditioners',
  !  each representing the part of the preconditioner associated to a certain level.
  !  For each level ilev, the base preconditioner K(ilev) is stored in precv(ilev)
  !  and is associated to a matrix A(ilev), obtained by 'tranferring' the original
  !  matrix A (i.e. the matrix to be preconditioned) to the level ilev, through smoothed
  !  aggregation.
  !
  !  The levels are numbered in increasing order starting from the finest one, i.e.
  !  level 1 is the finest level and A(1) is the matrix A. 
  !
  !  For details on the symmetrized hybrid multiplicative multilevel Schwarz
  !  preconditioner, see the Algorithm 3.2.2 of the book:
  !    B.F. Smith, P.E. Bjorstad & W.D. Gropp,
  !    Domain decomposition: parallel multilevel methods for elliptic partial
  !    differential equations, Cambridge University Press, 1996.
  !
  !  For a description of the arguments see mld_cmlprec_aply.
  !
  !  A sketch of the algorithm implemented in this routine is provided below.
  !  (AV(ilev; sm_pr_) denotes the smoothed prolongator from level ilev to
  !  level ilev-1, while AV(ilev; sm_pr_t_) denotes its transpose, i.e. the
  !  corresponding restriction operator from level ilev-1 to level ilev).
  !
  !    1.   X(1)  = Xext
  !
  !    2. ! Apply the base peconditioner at the finest level
  !         Y(1)  = (K(1)^(-1))*X(1)
  !
  !    3. ! Compute the residual at the finest level
  !         TX(1) = X(1) - A(1)*Y(1)
  !
  !    4.   DO ilev=2, nlev
  !
  !          ! Transfer the residual to the current (coarser) level
  !            X(ilev) = AV(ilev; sm_pr_t)*TX(ilev-1)
  !    
  !          ! Apply the base preconditioner at the current level.
  !          ! The sum over the subdomains is carried out in the
  !          ! application of K(ilev)
  !            Y(ilev) = (K(ilev)^(-1))*X(ilev)
  !
  !          ! Compute the residual at the current level
  !            TX(ilev) = (X(ilev)-A(ilev)*Y(ilev))
  !
  !         ENDDO
  !
  !    5.   DO ilev=NLEV-1,1,-1
  !
  !          ! Transfer Y(ilev+1) to the next finer level
  !            Y(ilev) = Y(ilev) + AV(ilev+1; sm_pr_)*Y(ilev+1)
  !
  !          ! Compute the residual at the current level and apply to it the
  !          ! base preconditioner. The sum over the subdomains is carried out
  !          ! in the application of K(ilev)     
  !            Y(ilev) = Y(ilev) + (K(ilev)**(-1))*(X(ilev)-A(ilev)*Y(ilev))
  !
  !         ENDDO
  !
  !    6.  Yext = beta*Yext + alpha*Y(1)
  !
  subroutine mlt_twoside_ml_aply(alpha,precv,x,beta,y,desc_data,trans,work,info)

    implicit none 

    ! Arguments
    type(psb_desc_type),intent(in)      :: desc_data
    type(mld_c_onelev_prec_type), intent(in) :: precv(:)
    complex(psb_spk_),intent(in)      :: alpha,beta
    complex(psb_spk_),intent(in)      :: x(:)
    complex(psb_spk_),intent(inout)   :: y(:)
    character, intent(in)               :: trans
    complex(psb_spk_),target          :: work(:)
    integer, intent(out)                :: info

    ! Local variables
    integer            :: ictxt,np,me,i, nr2l,nc2l,err_act
    integer            :: debug_level, debug_unit
    integer            :: nlev, ilev
    character(len=20)  :: name

    type psb_mlprec_wrk_type
      complex(psb_spk_), allocatable  :: tx(:),ty(:),x2l(:),y2l(:)
    end type psb_mlprec_wrk_type
    type(psb_mlprec_wrk_type), allocatable :: mlprec_wrk(:)

    name = 'mlt_twoside_ml_aply'
    info = 0
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    ictxt = psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)

    if (debug_level >= psb_debug_inner_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' Entry  ', size(precv)

    nlev = size(precv)
    allocate(mlprec_wrk(nlev),stat=info) 
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    ! STEP 1
    !
    ! Copy the input vector X
    !
    nc2l  = psb_cd_get_local_cols(precv(1)%base_desc)

    allocate(mlprec_wrk(1)%x2l(nc2l),mlprec_wrk(1)%y2l(nc2l), &
         & mlprec_wrk(1)%ty(nc2l), mlprec_wrk(1)%tx(nc2l), stat=info)

    if (info /= 0) then 
      info=4025
      call psb_errpush(info,name,i_err=(/4*nc2l,0,0,0,0/),&
           & a_err='complex(psb_spk_)')
      goto 9999
    end if

    call psb_geaxpby(cone,x,czero,mlprec_wrk(1)%x2l,&
         & precv(1)%base_desc,info)
    call psb_geaxpby(cone,x,czero,mlprec_wrk(1)%tx,&
         & precv(1)%base_desc,info)

    !
    ! STEP 2
    !
    ! Apply the base preconditioner at the finest level
    !
    call mld_baseprec_aply(cone,precv(1)%prec,mlprec_wrk(1)%x2l,&
         &  czero,mlprec_wrk(1)%y2l,precv(1)%base_desc,&
         &  trans,work,info)
    !
    ! STEP 3
    !
    ! Compute the residual at the finest level
    !
    mlprec_wrk(1)%ty = mlprec_wrk(1)%x2l
    if (info == 0) call psb_spmm(-cone,precv(1)%base_a,mlprec_wrk(1)%y2l,&
         & cone,mlprec_wrk(1)%ty,precv(1)%base_desc,info,&
         & work=work,trans=trans)
    if (info /=0) then
      call psb_errpush(4010,name,a_err='Fine level baseprec/residual')
      goto 9999
    end if

    !
    ! STEP 4
    !
    ! For each level but the finest one ...
    !
    do ilev = 2, nlev

      nc2l  = psb_cd_get_local_cols(precv(ilev)%base_desc)
      nr2l  = psb_cd_get_local_rows(precv(ilev)%base_desc)

      allocate(mlprec_wrk(ilev)%tx(nc2l),mlprec_wrk(ilev)%ty(nc2l),&
           &  mlprec_wrk(ilev)%y2l(nc2l),mlprec_wrk(ilev)%x2l(nc2l), stat=info)

      if (info /= 0) then 
        info=4025
        call psb_errpush(info,name,i_err=(/4*nc2l,0,0,0,0/),&
             & a_err='complex(psb_spk_)')
        goto 9999      
      end if

      ! Apply prolongator transpose, i.e. restriction
      call psb_forward_map(cone,mlprec_wrk(ilev-1)%ty,&
           & czero,mlprec_wrk(ilev)%x2l,&
           & precv(ilev)%map_desc,info,work=work)
      
      if (info /=0) then
        call psb_errpush(4001,name,a_err='Error during restriction')
        goto 9999
      end if

      call psb_geaxpby(cone,mlprec_wrk(ilev)%x2l,czero,mlprec_wrk(ilev)%tx,&
           & precv(ilev)%base_desc,info)
      !
      ! Apply the base preconditioner
      !
      if (info == 0) call mld_baseprec_aply(cone,precv(ilev)%prec,&
           & mlprec_wrk(ilev)%x2l,czero,mlprec_wrk(ilev)%y2l,&
           &precv(ilev)%base_desc,trans,work,info)
      !
      ! Compute the residual (at all levels but the coarsest one)
      !
      if(ilev < nlev) then
        mlprec_wrk(ilev)%ty = mlprec_wrk(ilev)%x2l
        if (info == 0) call psb_spmm(-cone,precv(ilev)%base_a,&
             & mlprec_wrk(ilev)%y2l,cone,mlprec_wrk(ilev)%ty,&
             & precv(ilev)%base_desc,info,work=work,trans=trans)
      endif
      if (info /=0) then
        call psb_errpush(4001,name,a_err='baseprec_aply/residual')
        goto 9999
      end if

    enddo

    !
    ! STEP 5
    !
    ! For each level but the coarsest one ...
    !
    do ilev=nlev-1, 1, -1

      !
      ! Apply prolongator
      !  
      call psb_backward_map(cone,mlprec_wrk(ilev+1)%y2l,&
           & cone,mlprec_wrk(ilev)%y2l,&
           & precv(ilev+1)%map_desc,info,work=work)

      if (info /=0 ) then
        call psb_errpush(4001,name,a_err='Error during restriction')
        goto 9999
      end if

      !
      ! Compute the residual
      !
      call psb_spmm(-cone,precv(ilev)%base_a,mlprec_wrk(ilev)%y2l,&
           & cone,mlprec_wrk(ilev)%tx,precv(ilev)%base_desc,info,&
           & work=work,trans=trans)
      !
      ! Apply the base preconditioner
      !
      if (info == 0) call mld_baseprec_aply(cone,precv(ilev)%prec,mlprec_wrk(ilev)%tx,&
           & cone,mlprec_wrk(ilev)%y2l,precv(ilev)%base_desc, trans, work,info)
      if (info /= 0) then
        call psb_errpush(4001,name,a_err='Error: residual/baseprec_aply')
        goto 9999
      end if
    enddo

    !
    ! STEP 6
    !
    ! Compute the output vector Y
    !
    call psb_geaxpby(alpha,mlprec_wrk(1)%y2l,beta,y,&
         &   precv(1)%base_desc,info)

    if (info /= 0) then
      call psb_errpush(4001,name,a_err='Error final update')
      goto 9999
    end if



    deallocate(mlprec_wrk,stat=info)
    if (info /= 0) then 
      call psb_errpush(4000,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine mlt_twoside_ml_aply

end subroutine mld_cmlprec_aply

