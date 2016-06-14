!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.3)
!!$  
!!$  (C) Copyright 2008, 2010, 2012, 2015
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
! File: mld_dmlprec_aply.f90
!
! Subroutine: mld_dmlprec_aply
! Version:    real
!
!  This routine computes
!  
!                        Y = beta*Y + alpha*op(M^(-1))*X,
!  where 
!  - M is a multilevel domain decomposition (Schwarz) preconditioner associated
!    to a certain matrix A and stored in p,
!  - op(M^(-1)) is M^(-1) or its transpose, according to the value of trans,
!  - X and Y are vectors,
!  - alpha and beta are scalars.
!
!  For each level we have as many submatrices as processes (except for the coarsest
!  level where we might have a replicated index space) and each process takes care
!  of one submatrix.
!
!  A multilevel preconditioner is regarded as an array of 'one-level' data structures,
!  each containing the part of the preconditioner associated to a certain level
!  (for more details see the description of mld_Tonelev_type in mld_prec_type.f90).
!  For each level ilev, the 'base preconditioner' K(ilev) is stored in
!   p%precv(ilev)%prec
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
!   alpha      -   real(psb_dpk_), input.
!                  The scalar alpha.
!   p          -   type(mld_dprec_type), input.
!                  The multilevel preconditioner data structure containing the
!                  local part of the preconditioner to be applied.
!      Note that nlev = size(p%precv) = number of levels.
!      p%precv(ilev)%prec      -  type(psb_dbaseprec_type)
!                                 The 'base preconditioner' for the current level
!      p%precv(ilev)%ac        -  type(psb_dspmat_type) 
!                                 The local part of the matrix A(ilev).
!      p%precv(ilev)%desc_ac   -  type(psb_desc_type).
!                                 The communication descriptor associated to the sparse
!                                 matrix A(ilev)
!      p%precv(ilev)%map       -  type(psb_inter_desc_type)
!                                 Stores the linear operators mapping level (ilev-1)
!                                 to (ilev) and vice versa. These are the restriction
!                                 and prolongation operators described in the sequel. 
!      p%precv(ilev)%iprcparm  -  integer, dimension(:), allocatable.
!                                 The integer parameters defining the multilevel
!                                 strategy 
!      p%precv(ilev)%rprcparm  -  real(psb_dpk_), dimension(:), allocatable.
!                                 The real parameters defining the multilevel strategy
!      p%precv(ilev)%mlia      -  integer, dimension(:), allocatable.
!                                 The aggregation map (ilev-1) --> (ilev).
!      p%precv(ilev)%nlaggr    -  integer, dimension(:), allocatable.
!                                 The number of aggregates (rows of A(ilev)) on the
!                                 various processes. 
!      p%precv(ilev)%base_a    -  type(psb_dspmat_type), pointer.
!                                 Pointer (really a pointer!) to the base matrix of
!                                 the current level, i.e. the local part of A(ilev);
!                                 so we have a unified treatment of residuals. We
!                                 need this to avoid passing explicitly the matrix
!                                 A(ilev) to the routine which applies the
!                                 preconditioner.
!      p%precv(ilev)%base_desc -  type(psb_desc_type), pointer.
!                                 Pointer to the communication descriptor associated
!                                 to the sparse matrix pointed by base_a.  
!                  
!   x          -  real(psb_dpk_), dimension(:), input.
!                 The local part of the vector X.
!   beta       -  real(psb_dpk_), input.
!                 The scalar beta.
!   y          -  real(psb_dpk_), dimension(:), input/output.
!                 The local part of the vector Y.
!   desc_data  -  type(psb_desc_type), input.
!                 The communication descriptor associated to the matrix to be
!                 preconditioned.
!   trans      -  character, optional.
!                 If trans='N','n' then op(M^(-1)) = M^(-1);
!                 if trans='T','t' then op(M^(-1)) = M^(-T) (transpose of M^(-1)).
!   work       -  real(psb_dpk_), dimension (:), optional, target.
!                 Workspace. Its size must be at least 4*desc_data%get_local_cols().
!   info       -  integer, output.
!                 Error code.
!
!   Note that when the LU factorization of the matrix A(ilev) is computed instead of
!   the ILU one, by using UMFPACK or SuperLU, the corresponding L and U factors
!   are stored in data structures provided by UMFPACK or SuperLU and pointed by
!   p%precv(ilev)%prec%iprcparm(mld_umf_ptr) or p%precv(ilev)%prec%iprcparm(mld_slu_ptr),
!   respectively.
!
!   This routine is formulated in a recursive way, so it is very compact.
!   In the original code the recursive formulation was explicitly unrolled.
!   The description of the various alternatives is given below in the explicit
!   formulation, hopefully it will be clear enough when related to the
!   recursive formulation. 
!   
!   This routine computes
!                        Y = beta*Y + alpha*op(M^(-1))*X,
!  where 
!  - M is a multilevel domain decomposition (Schwarz) preconditioner
!    associated to a certain matrix A and stored in p,
!  - op(M^(-1)) is M^(-1) or its transpose, according to the value of trans,
!  - X and Y are vectors,
!  - alpha and beta are scalars.
!
!  For each level we have as many submatrices as processes (except for the coarsest
!  level where we might have a replicated index space) and each process takes care
!  of one submatrix. 
!
!  The multilevel preconditioner is regarded as an array of 'one-level' data structures,
!  each containing the part of the preconditioner associated to a certain level
!  (for more details see the description of mld_Tonelev_type in mld_prec_type.f90).
!  For each level ilev, the 'base preconditioner' K(ilev) is stored in
!  p%precv(ilev)%prec
!  and is associated to a matrix A(ilev), obtained by 'tranferring' the original
!  matrix A (i.e. the matrix to be preconditioned) to the level ilev, through smoothed
!  aggregation.
!  The levels are numbered in increasing order starting from the finest one, i.e.
!  level 1 is the finest level and A(1) is the matrix A. 
!
!
! Additive multilevel
!    This is additive both within the levels and among levels.
!
!  For details on the additive multilevel Schwarz preconditioner, see
!  Algorithm 3.1.1 in the book:
!    B.F. Smith, P.E. Bjorstad & W.D. Gropp,
!    Domain decomposition: parallel multilevel methods for elliptic partial
!    differential equations, Cambridge University Press, 1996.
!
!  (P(ilev) denotes the smoothed prolongator from level ilev to level
!  ilev-1, while PT(ilev) denotes its transpose, i.e. the corresponding
!  restriction operator from level ilev-1 to level ilev).
!
!   0. Transfer the outer vector Xest to x(1) (inner X at level 1)
!   
!   1. If ilev > 1   Transfer x(ilev-1) to the current level:
!           x(ilev) = PT(ilev)*x(ilev-1)
!
!   2. Apply the base preconditioner at the current level:
!         ! The sum over the subdomains is carried out in the
!         ! application of K(ilev)
!          y(ilev) = (K(ilev)^(-1))*x(ilev)
!
!   3. If ilev < nlevel
!         a. Call recursively itself
!         b. Transfer y(ilev+1) to the current level:
!           y(ilev) = y(ilev) + P(ilev+1)*y(ilev+1)
!           
!    4. if ilev == 1  Transfer the inner y to the external:
!         Yext = beta*Yext + alpha*y(1)
!
!
!
!  Hybrid multiplicative---pre-smoothing
!  
!  The preconditioner M is hybrid in the sense that it is multiplicative through the
!  levels and additive inside a level. 
!
!  For details on the pre-smoothed hybrid multiplicative multilevel Schwarz
!  preconditioner, see Algorithm 3.2.1 in the book:
!    B.F. Smith, P.E. Bjorstad & W.D. Gropp,
!    Domain decomposition: parallel multilevel methods for elliptic partial
!    differential equations, Cambridge University Press, 1996.
!
!
!   0. Transfer the outer vector Xest to x(1) (inner X at level 1)
!   
!   1. If ilev >1   Transfer x(ilev-1) to the current level:
!           x(ilev) = PT(ilev)*x(ilev-1)
!
!   2. Apply the base preconditioner at the current level:
!         ! The sum over the subdomains is carried out in the
!         ! application of K(ilev).
!          y(ilev) = (K(ilev)^(-1))*x(ilev)
!
!   3. If ilev < nlevel
!         a. Compute the residual:
!            r(ilev) = x(ilev) - A(ilev)*y(ilev)
!         b. Call recursively itself passing
!            r(ilev) for transfer to the next level
!            (r(ilev) matches x(ilev-1) in step 1)
!
!         c. Transfer y(ilev+1) to the current level:
!            y(ilev) = y(ilev) + P(ilev+1)*y(ilev+1)
!           
!    4. if ilev == 1  Transfer the inner y to the external:
!         Yext = beta*Yext + alpha*y(1)
!
!
!
!  Hybrid multiplicative, post-smoothing variant
!
!   0. Transfer the outer vector Xest to x(1) (inner X at level 1)
!   
!   1. If ilev > 1   Transfer x(ilev-1) to the current level:
!           x(ilev) = PT(ilev)*x(ilev-1)
!
!   2.  If ilev < nlev 
!         a. Call recursively itself passing
!            x(ilev) for transfer to the next level
!         b. Transfer y(ilev+1) to the current level:
!            y(ilev) = P(ilev+1)*y(ilev+1)
!         c. Compute the residual:
!            x(ilev) = x(ilev) - A(ilev)*y(ilev)
!         d. Apply the base preconditioner to the residual at the current level:
!            ! The sum over the subdomains is carried out in the
!            ! application of K(ilev)
!            y(ilev) = y(ilev) + (K(ilev)^(-1))*x(ilev)
!       Else
!         Apply the base preconditioner to the residual at the current level:
!         ! The sum over the subdomains is carried out in the
!         ! application of K(ilev)
!         y(ilev) = (K(ilev)^(-1))*x(ilev)
!      
!    4. if ilev == 1 Transfer the inner Y to the external:
!         Yext = beta*Yext + alpha*Y(1)
!
!
!
!  Hybrid multiplicative, pre- and post-smoothing (two-side) variant
!
!  For details on the symmetrized hybrid multiplicative multilevel Schwarz
!  preconditioner, see Algorithm 3.2.2 in the book:
!    B.F. Smith, P.E. Bjorstad & W.D. Gropp,
!    Domain decomposition: parallel multilevel methods for elliptic partial
!    differential equations, Cambridge University Press, 1996.
!
!
!   0. Transfer the outer vector Xest to x(1) (inner X at level 1)
!   
!   1. If ilev > 1   Transfer x(ilev-1) to the current level:
!           x(ilev) = PT(ilev)*x(ilev-1)
!
!   2. Apply the base preconditioner at the current level:
!         ! The sum over the subdomains is carried out in the
!         ! application of K(ilev)
!          y(ilev) = (K(ilev)^(-1))*x(ilev)
!
!   3. If ilev < nlevel
!         a. Compute the residual:
!            r(ilev) = x(ilev) - A(ilev)*y(ilev)
!         b. Call recursively itself passing
!            r(ilev) for transfer to the next level
!            (r(ilev) matches x(ilev-1) in step 1) 
!         c. Transfer y(ilev+1) to the current level:
!            y(ilev) = y(ilev) + P(ilev+1)*y(ilev+1)
!         d. Compute the residual:
!            r(ilev) = x(ilev) - A(ilev)*y(ilev)
!         e. Apply the base preconditioner at the current level to the residual:
!            ! The sum over the subdomains is carried out in the
!            ! application of K(ilev)
!            y(ilev) = y(ilev) + (K(ilev)^(-1))*r(ilev)
!           
!    5. if ilev == 1 Transfer the inner Y to the external:
!         Yext = beta*Yext + alpha*Y(1)
!
!
subroutine mld_dmlprec_aply(alpha,p,x,beta,y,desc_data,trans,work,info)

  use psb_base_mod
  use mld_d_inner_mod, mld_protect_name => mld_dmlprec_aply

  implicit none

  ! Arguments
  type(psb_desc_type),intent(in)    :: desc_data
  type(mld_dprec_type), intent(inout)  :: p
  real(psb_dpk_),intent(in)         :: alpha,beta
  real(psb_dpk_),intent(inout)      :: x(:)
  real(psb_dpk_),intent(inout)      :: y(:)
  character, intent(in)             :: trans
  real(psb_dpk_),target             :: work(:)
  integer(psb_ipk_), intent(out)              :: info

  ! Local variables
  integer(psb_ipk_)  :: ictxt, np, me
  integer(psb_ipk_)  :: err_act
  integer(psb_ipk_)  :: debug_level, debug_unit, nlev,nc2l,nr2l,level
  character(len=20)  :: name
  character          :: trans_
  type mld_mlprec_wrk_type
    real(psb_dpk_), allocatable  :: tx(:), ty(:), x2l(:), y2l(:)
  end type mld_mlprec_wrk_type
  type(mld_mlprec_wrk_type), allocatable, target  :: mlprec_wrk(:)

  name='mld_dmlprec_aply'
  info = psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = desc_data%get_context()
  call psb_info(ictxt, me, np)

  if (debug_level >= psb_debug_inner_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ' Entry  ', size(p%precv)

  trans_ = psb_toupper(trans)

  nlev = size(p%precv)
  allocate(mlprec_wrk(nlev),stat=info) 
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999      
  end if
  level = 1

  nc2l  = p%precv(level)%base_desc%get_local_cols()
  nr2l  = p%precv(level)%base_desc%get_local_rows()
  allocate(mlprec_wrk(level)%x2l(nc2l),mlprec_wrk(level)%y2l(nc2l),&
       & stat=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,&
         & i_err=(/ione*(size(x)+size(y)),izero,izero,izero,izero/),&
         & a_err='real(psb_dpk_)')
    goto 9999      
  end if

  mlprec_wrk(level)%x2l(:) = x(:) 
  mlprec_wrk(level)%y2l(:) = dzero 

  call inner_ml_aply(level,p,mlprec_wrk,trans_,work,info)    

  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,&
         & a_err='Inner prec aply')
    goto 9999
  end if

  call psb_geaxpby(alpha,mlprec_wrk(level)%y2l,beta,y,&
       &   p%precv(level)%base_desc,info)

  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,&
         & a_err='Error final update')
    goto 9999
  end if


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

contains

  recursive subroutine inner_ml_aply(level,p,mlprec_wrk,trans,work,info)    

    implicit none 

    ! Arguments
    integer(psb_ipk_)                           :: level 
    type(mld_dprec_type), target, intent(inout) :: p
    type(mld_mlprec_wrk_type), intent(inout)    :: mlprec_wrk(:)
    character, intent(in)                       :: trans
    real(psb_dpk_),target            :: work(:)
    integer(psb_ipk_), intent(out)    :: info

    ! Local variables
    integer(psb_ipk_)  :: ictxt,np,me
    integer(psb_ipk_)  :: i, nr2l,nc2l,err_act
    integer(psb_ipk_)  :: debug_level, debug_unit
    integer(psb_ipk_)  :: nlev, ilev, sweeps
    character(len=20)  :: name

    name = 'inner_ml_aply'
    info = psb_success_
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    nlev = size(p%precv)
    if ((level < 1) .or. (level > nlev)) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='wrong call level to inner_ml')
      goto 9999      
    end if
    ictxt = p%precv(level)%base_desc%get_context()
    call psb_info(ictxt, me, np)

    if (level > 1) then 
      nc2l  = p%precv(level)%base_desc%get_local_cols()
      nr2l  = p%precv(level)%base_desc%get_local_rows()
      allocate(mlprec_wrk(level)%x2l(nc2l),&
           & mlprec_wrk(level)%y2l(nc2l),&
           & stat=info)
      if (info /= psb_success_) then 
        info=psb_err_alloc_request_
        call psb_errpush(info,name,i_err=(/2*nc2l,izero,izero,izero,izero/),&
             & a_err='real(psb_dpk_)')
        goto 9999      
      end if
    end if

    select case(p%precv(level)%parms%ml_type) 

    case(mld_no_ml_)
      !
      ! No preconditioning, should not really get here
      ! 
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='mld_no_ml_ in mlprc_aply?')
      goto 9999      

    case(mld_add_ml_)
      !
      ! Additive multilevel
      !

      if (level > 1) then 
        ! Apply the restriction
        call psb_map_X2Y(done,mlprec_wrk(level-1)%x2l,&
             & dzero,mlprec_wrk(level)%x2l,&
             & p%precv(level)%map,info,work=work)

        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='Error during restriction')
          goto 9999
        end if

      end if

      sweeps = p%precv(level)%parms%sweeps 
      call p%precv(level)%sm%apply(done,&
           & mlprec_wrk(level)%x2l,dzero,mlprec_wrk(level)%y2l,&
           & p%precv(level)%base_desc, trans,&
           & sweeps,work,info)

      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during smoother_apply')
        goto 9999
      end if

      if (level < nlev) then
        call inner_ml_aply(level+1,p,mlprec_wrk,trans,work,info)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='Error in recursive call')
          goto 9999
        end if
        !
        ! Apply the prolongator
        !  
        call psb_map_Y2X(done,mlprec_wrk(level+1)%y2l,&
             & done,mlprec_wrk(level)%y2l,&
             & p%precv(level+1)%map,info,work=work)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='Error during prolongation')
          goto 9999
        end if

      end if

    case(mld_mult_ml_)
      ! 
      !  Multiplicative multilevel (multiplicative among the levels, additive inside
      !  each level)
      !
      !  Pre/post-smoothing versions.
      !  Note that the transpose switches pre <-> post.
      !
      select case(p%precv(level)%parms%smoother_pos)

      case(mld_post_smooth_)

        select case (trans_) 
        case('N')
          if (level > 1) then 
            ! Apply the restriction
            call psb_map_X2Y(done,mlprec_wrk(level-1)%x2l,&
                 & dzero,mlprec_wrk(level)%x2l,&
                 & p%precv(level)%map,info,work=work)

            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error during restriction')
              goto 9999
            end if
          end if

          ! This is one step of post-smoothing 
          if (level < nlev) then 
            call inner_ml_aply(level+1,p,mlprec_wrk,trans,work,info) 
            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error in recursive call')
              goto 9999
            end if
            !
            ! Apply the prolongator
            !  
            call psb_map_Y2X(done,mlprec_wrk(level+1)%y2l,&
                 & dzero,mlprec_wrk(level)%y2l,&
                 & p%precv(level+1)%map,info,work=work)

            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error during prolongation')
              goto 9999
            end if
            !
            ! Compute the residual
            !
            call psb_spmm(-done,p%precv(level)%base_a,mlprec_wrk(level)%y2l,&
                 & done,mlprec_wrk(level)%x2l,p%precv(level)%base_desc,info,&
                 & work=work,trans=trans)
            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error during residue')
              goto 9999
            end if

            sweeps = p%precv(level)%parms%sweeps_post 
            call p%precv(level)%sm2%apply(done,&
                 & mlprec_wrk(level)%x2l,done,mlprec_wrk(level)%y2l,&
                 & p%precv(level)%base_desc, trans,&
                 & sweeps,work,info)
            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error during smoother_apply')
              goto 9999
            end if
           
          else
            sweeps = p%precv(level)%parms%sweeps 
            call p%precv(level)%sm%apply(done,&
                 & mlprec_wrk(level)%x2l,dzero,mlprec_wrk(level)%y2l,&
                 & p%precv(level)%base_desc, trans,&
                 & sweeps,work,info)
            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error during smoother_apply')
              goto 9999
            end if
            
          end if

        case('T','C')

          ! Post-smoothing transpose is pre-smoothing


          if (level > 1) then 
            ! Apply the restriction
            call psb_map_X2Y(done,mlprec_wrk(level-1)%x2l,&
                 & dzero,mlprec_wrk(level)%x2l,&
                 & p%precv(level)%map,info,work=work)

            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error during restriction')
              goto 9999
            end if


          end if

          !
          ! Apply the base preconditioner
          !
          if (level < nlev) then 
            sweeps = p%precv(level)%parms%sweeps_post
            call p%precv(level)%sm2%apply(done,&
                 & mlprec_wrk(level)%x2l,dzero,mlprec_wrk(level)%y2l,&
                 & p%precv(level)%base_desc, trans,&
                 & sweeps,work,info)
          else
            sweeps = p%precv(level)%parms%sweeps
            call p%precv(level)%sm%apply(done,&
                 & mlprec_wrk(level)%x2l,dzero,mlprec_wrk(level)%y2l,&
                 & p%precv(level)%base_desc, trans,&
                 & sweeps,work,info)
          end if

          if (info /= psb_success_) then
            call psb_errpush(psb_err_internal_error_,name,&
                 & a_err='Error during smoother_apply')
            goto 9999
          end if
          
          !
          ! Compute the residual (at all levels but the coarsest one)
          !
          if (level < nlev) then
            call psb_spmm(-done,p%precv(level)%base_a,&
                 & mlprec_wrk(level)%y2l,done,mlprec_wrk(level)%x2l,&
                 & p%precv(level)%base_desc,info,work=work,trans=trans)

            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error during residue')
              goto 9999
            end if

            call inner_ml_aply(level+1,p,mlprec_wrk,trans,work,info) 
            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error in recursive call')
              goto 9999
            end if


            call psb_map_Y2X(done,mlprec_wrk(level+1)%y2l,&
                 & done,mlprec_wrk(level)%y2l,&
                 & p%precv(level+1)%map,info,work=work)

            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error during prolongation')
              goto 9999
            end if


          end if

        case default
          info = psb_err_internal_error_
          call psb_errpush(info,name,a_err='invalid trans')
          goto 9999      
        end select

      case(mld_pre_smooth_)

        select case (trans_) 
        case('N')
          ! One step of pre-smoothing 


          if (level > 1) then 
            ! Apply the restriction
            call psb_map_X2Y(done,mlprec_wrk(level-1)%x2l,&
                 & dzero,mlprec_wrk(level)%x2l,&
                 & p%precv(level)%map,info,work=work)

            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error during restriction')
              goto 9999
            end if

          end if

          !
          ! Apply the base preconditioner
          !
          if (level < nlev) then 
            sweeps = p%precv(level)%parms%sweeps_pre
          else
            sweeps = p%precv(level)%parms%sweeps
          end if
          call p%precv(level)%sm%apply(done,&
               & mlprec_wrk(level)%x2l,dzero,mlprec_wrk(level)%y2l,&
               & p%precv(level)%base_desc, trans,&
               & sweeps,work,info)

          if (info /= psb_success_) then
            call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during smoother_apply')
            goto 9999
          end if

          !
          ! Compute the residual (at all levels but the coarsest one)
          !
          if (level < nlev) then
            call psb_spmm(-done,p%precv(level)%base_a,&
                 & mlprec_wrk(level)%y2l,done,mlprec_wrk(level)%x2l,&
                 & p%precv(level)%base_desc,info,work=work,trans=trans)
            if (info /= psb_success_) goto 9999
            call inner_ml_aply(level+1,p,mlprec_wrk,trans,work,info) 
            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error in recursive call')
              goto 9999
            end if


            call psb_map_Y2X(done,mlprec_wrk(level+1)%y2l,&
                 & done,mlprec_wrk(level)%y2l,&
                 & p%precv(level+1)%map,info,work=work)

            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error during prolongation')
              goto 9999
            end if


          end if


        case('T','C')

          ! pre-smooth transpose is post-smoothing 


          if (level > 1) then 
            ! Apply the restriction
            call psb_map_X2Y(done,mlprec_wrk(level-1)%x2l,&
                 & dzero,mlprec_wrk(level)%x2l,&
                 & p%precv(level)%map,info,work=work)

            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error during restriction')
              goto 9999
            end if

          end if

          if (level < nlev) then 
            call inner_ml_aply(level+1,p,mlprec_wrk,trans,work,info) 
            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error in recursive call')
              goto 9999
            end if

            !
            ! Apply the prolongator
            !  
            call psb_map_Y2X(done,mlprec_wrk(level+1)%y2l,&
                 & dzero,mlprec_wrk(level)%y2l,&
                 & p%precv(level+1)%map,info,work=work)
            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error during prolongation')
              goto 9999
            end if

            !
            ! Compute the residual
            !
            call psb_spmm(-done,p%precv(level)%base_a,mlprec_wrk(level)%y2l,&
                 & done,mlprec_wrk(level)%x2l,p%precv(level)%base_desc,info,&
                 & work=work,trans=trans)
            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error during residue')
              goto 9999
            end if


            sweeps = p%precv(level)%parms%sweeps_pre 
            call p%precv(level)%sm%apply(done,&
                 & mlprec_wrk(level)%x2l,done,mlprec_wrk(level)%y2l,&
                 & p%precv(level)%base_desc, trans,&
                 & sweeps,work,info)
            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error during smoother_apply')
              goto 9999
            end if
          else
            sweeps = p%precv(level)%parms%sweeps 
            call p%precv(level)%sm%apply(done,&
                 & mlprec_wrk(level)%x2l,dzero,mlprec_wrk(level)%y2l,&
                 & p%precv(level)%base_desc, trans,&
                 & sweeps,work,info)
            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error during smoother_apply')
              goto 9999
            end if

          end if

        case default
          info = psb_err_internal_error_
          call psb_errpush(info,name,a_err='invalid trans')
          goto 9999      
        end select

      case(mld_twoside_smooth_)

        ! CHECK
        if (.not.(associated(p%precv(level)%sm2,p%precv(level)%sm2a))) then
          write(0,*) 'inner_ml_aply: unassociated sm2 at level ',level
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='Error during restriction')
          goto 9999
        end if
        nc2l  = p%precv(level)%base_desc%get_local_cols()
        nr2l  = p%precv(level)%base_desc%get_local_rows()
        allocate(mlprec_wrk(level)%ty(nc2l), mlprec_wrk(level)%tx(nc2l), stat=info)
        if (info /= psb_success_) then 
          info=psb_err_alloc_request_
          call psb_errpush(info,name,i_err=(/2*nc2l,izero,izero,izero,izero/),&
               & a_err='real(psb_dpk_)')
          goto 9999      
        end if

        if (level > 1) then 
          ! Apply the restriction
          call psb_map_X2Y(done,mlprec_wrk(level-1)%ty,&
               & dzero,mlprec_wrk(level)%x2l,&
               & p%precv(level)%map,info,work=work)

          if (info /= psb_success_) then
            call psb_errpush(psb_err_internal_error_,name,&
                 & a_err='Error during restriction')
            goto 9999
          end if
        end if
        call psb_geaxpby(done,mlprec_wrk(level)%x2l,dzero,mlprec_wrk(level)%tx,&
             & p%precv(level)%base_desc,info)
        !
        ! Apply the base preconditioner
        !
        if (level < nlev) then 
          if (trans == 'N') then 
            sweeps = p%precv(level)%parms%sweeps_pre
          else
            sweeps = p%precv(level)%parms%sweeps_post
          end if
        else
          sweeps = p%precv(level)%parms%sweeps
        end if

        if (trans == 'N') then 
          if (info == psb_success_) call p%precv(level)%sm%apply(done,&
               & mlprec_wrk(level)%x2l,dzero,mlprec_wrk(level)%y2l,&
               & p%precv(level)%base_desc, trans,&
               & sweeps,work,info)
        else
          if (info == psb_success_) call p%precv(level)%sm2%apply(done,&
               & mlprec_wrk(level)%x2l,dzero,mlprec_wrk(level)%y2l,&
               & p%precv(level)%base_desc, trans,&
               & sweeps,work,info)          
        end if
          
        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='Error during smoother_apply')
          goto 9999
        end if
        !
        ! Compute the residual (at all levels but the coarsest one)
        ! and call recursively
        !
        if(level < nlev) then
          mlprec_wrk(level)%ty = mlprec_wrk(level)%x2l
          if (info == psb_success_) call psb_spmm(-done,p%precv(level)%base_a,&
               & mlprec_wrk(level)%y2l,done,mlprec_wrk(level)%ty,&
               & p%precv(level)%base_desc,info,work=work,trans=trans)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_internal_error_,name,&
                 & a_err='Error during residue')
            goto 9999
          end if

          call inner_ml_aply(level+1,p,mlprec_wrk,trans,work,info)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_internal_error_,name,&
                 & a_err='Error in recursive call')
            goto 9999
          end if


          !
          ! Apply the prolongator
          !  
          call psb_map_Y2X(done,mlprec_wrk(level+1)%y2l,&
               & done,mlprec_wrk(level)%y2l,&
               & p%precv(level+1)%map,info,work=work)
          
          if (info /= psb_success_) then
            call psb_errpush(psb_err_internal_error_,name,&
                 & a_err='Error during prolongation')
            goto 9999
          end if

          !
          ! Compute the residual
          !
          call psb_spmm(-done,p%precv(level)%base_a,mlprec_wrk(level)%y2l,&
               & done,mlprec_wrk(level)%tx,p%precv(level)%base_desc,info,&
               & work=work,trans=trans)
            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error during residue')
              goto 9999
            end if
          !
          ! Apply the base preconditioner
          !
          if (trans == 'N') then 
            sweeps = p%precv(level)%parms%sweeps_post
          else
            sweeps = p%precv(level)%parms%sweeps_pre
          end if
          if (trans == 'N') then 
            if (info == psb_success_) call p%precv(level)%sm2%apply(done,&
                 & mlprec_wrk(level)%tx,done,mlprec_wrk(level)%y2l,&
                 & p%precv(level)%base_desc, trans,&
                 & sweeps,work,info)
          else
            if (info == psb_success_) call p%precv(level)%sm%apply(done,&
                 & mlprec_wrk(level)%tx,done,mlprec_wrk(level)%y2l,&
                 & p%precv(level)%base_desc, trans,&
                 & sweeps,work,info)
          end if
          
          if (info /= psb_success_) then
            call psb_errpush(psb_err_internal_error_,name,&
                 & a_err='Error during smoother_apply')
            goto 9999
          end if
          
        endif

      case default
        info = psb_err_from_subroutine_ai_
        call psb_errpush(info,name,a_err='invalid smooth_pos',&
             &  i_Err=(/p%precv(level)%parms%smoother_pos,izero,izero,izero,izero/))
        goto 9999      

      end select

    case default
      info = psb_err_from_subroutine_ai_
      call psb_errpush(info,name,a_err='invalid mltype',&
           &  i_Err=(/p%precv(level)%parms%ml_type,izero,izero,izero,izero/))
      goto 9999      

    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine inner_ml_aply

end subroutine mld_dmlprec_aply

subroutine mld_dmlprec_aply_vect(alpha,p,x,beta,y,desc_data,trans,work,info)

  use psb_base_mod
  use mld_d_inner_mod, mld_protect_name => mld_dmlprec_aply_vect

  implicit none

  ! Arguments
  type(psb_desc_type),intent(in)        :: desc_data
  type(mld_dprec_type), intent(inout) :: p
  real(psb_dpk_),intent(in)            :: alpha,beta
  type(psb_d_vect_type),intent(inout) :: x
  type(psb_d_vect_type),intent(inout) :: y
  character, intent(in)                 :: trans
  real(psb_dpk_),target                :: work(:)
  integer(psb_ipk_), intent(out)        :: info

  ! Local variables
  integer(psb_ipk_)  :: ictxt, np, me
  integer(psb_ipk_)  :: debug_level, debug_unit, nlev,nc2l,nr2l,level, err_act
  character(len=20)  :: name
  character          :: trans_
  real(psb_dpk_)     :: par
  type mld_mlprec_wrk_type
    real(psb_dpk_), allocatable  :: tx(:), ty(:), x2l(:), y2l(:)
    type(psb_d_vect_type)  :: vtx, vty, vx2l, vy2l
  end type mld_mlprec_wrk_type
  type(mld_mlprec_wrk_type), allocatable, target  :: mlprec_wrk(:)

  name='mld_dmlprec_aply'
  info = psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = desc_data%get_context()
  call psb_info(ictxt, me, np)

  if (debug_level >= psb_debug_inner_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ' Entry  ', size(p%precv)

  trans_ = psb_toupper(trans)
  nlev = size(p%precv)
  allocate(mlprec_wrk(nlev),stat=info) 
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999      
  end if
  level = 1
  do level = 1, nlev
    call psb_geasb(mlprec_wrk(level)%vx2l,&
         & p%precv(level)%base_desc,info,&
         & scratch=.true.,mold=x%v)
    call psb_geasb(mlprec_wrk(level)%vy2l,&
         & p%precv(level)%base_desc,info,&
         & scratch=.true.,mold=x%v)
    call psb_geasb(mlprec_wrk(level)%vtx,&
         & p%precv(level)%base_desc,info,&
         & scratch=.true.,mold=x%v)
    call psb_geasb(mlprec_wrk(level)%vty,&
         & p%precv(level)%base_desc,info,&
         & scratch=.true.,mold=x%v)
    if (psb_errstatus_fatal()) then 
      nc2l = p%precv(level)%base_desc%get_local_cols()
      info=psb_err_alloc_request_
      call psb_errpush(info,name,i_err=(/2*nc2l,izero,izero,izero,izero/),&
           & a_err='real(psb_dpk_)')
      goto 9999      
    end if
  end do
  level = 1

  call psb_geaxpby(done,x,dzero,mlprec_wrk(level)%vx2l,p%precv(level)%base_desc,info)
  call mlprec_wrk(level)%vy2l%zero() 


  call inner_ml_aply(level,p,mlprec_wrk,trans_,work,info)    

  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,&
         & a_err='Inner prec aply')
    goto 9999
  end if

  call psb_geaxpby(alpha,mlprec_wrk(level)%vy2l,beta,y,&
       &   p%precv(level)%base_desc,info)
  do level = 1, nlev

    call mlprec_wrk(level)%vx2l%free(info)
    call mlprec_wrk(level)%vy2l%free(info)
    call mlprec_wrk(level)%vtx%free(info)
    call mlprec_wrk(level)%vty%free(info)
    if (psb_errstatus_fatal()) then 
      info=psb_err_alloc_request_
      nc2l = p%precv(level)%base_desc%get_local_cols()
      call psb_errpush(info,name,i_err=(/2*nc2l,izero,izero,izero,izero/),&
           & a_err='real(psb_dpk_)')
      goto 9999      
    end if
  end do

  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,&
         & a_err='Error final update')
    goto 9999
  end if


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

contains

  recursive subroutine inner_ml_aply(level,p,mlprec_wrk,trans,work,info, U)    

    implicit none 

    ! Arguments
    integer(psb_ipk_)                           :: level 
    type(mld_dprec_type), target, intent(inout) :: p
    type(mld_mlprec_wrk_type), intent(inout), target    :: mlprec_wrk(:)
    character, intent(in)                       :: trans
    real(psb_dpk_),target                      :: work(:)
    integer(psb_ipk_), intent(out)              :: info
    type(psb_d_vect_type),intent(inout), optional :: u

    type(psb_d_vect_type) :: res
    type(psb_d_vect_type), pointer :: current
    integer(psb_ipk_) :: sweeps_post, sweeps_pre
    ! Local variables
    integer(psb_ipk_)  :: ictxt,np,me
    integer(psb_ipk_)  :: i, nr2l,nc2l,err_act
    integer(psb_ipk_)  :: debug_level, debug_unit
    integer(psb_ipk_)  :: nlev, ilev, sweeps
    logical            :: pre, post
    character(len=20)  :: name



    name = 'inner_ml_aply'
    info = psb_success_
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    nlev = size(p%precv)
    if ((level < 1) .or. (level > nlev)) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='wrong call level to inner_ml')
      goto 9999      
    end if
    ictxt = p%precv(level)%base_desc%get_context()
    call psb_info(ictxt, me, np)

    nc2l  = p%precv(level)%base_desc%get_local_cols()
    nr2l  = p%precv(level)%base_desc%get_local_rows()
    if(debug_level > 1) then
      write(debug_unit,*) me,' inner_ml_aply at level ',level
    end if

    select case(p%precv(level)%parms%ml_type) 

    case(mld_no_ml_)
      !
      ! No preconditioning, should not really get here
      ! 
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='mld_no_ml_ in mlprc_aply?')
      goto 9999      

    case(mld_add_ml_)

      call mld_d_inner_add(p, mlprec_wrk, level, trans, work)


    case(mld_mult_ml_)
      ! 
      !  Multiplicative multilevel (multiplicative among the levels, additive inside
      !  each level)
      !
      !  Pre/post-smoothing versions.
      !  Note that the transpose switches pre <-> post.
      !
      select case(p%precv(level)%parms%smoother_pos)

      case(mld_post_smooth_)
        p%precv(level)%parms%sweeps_pre = 0
        call mld_d_inner_mult(p, mlprec_wrk, level, trans, work)        
        

      case(mld_pre_smooth_)
        p%precv(level)%parms%sweeps_post = 0
        call mld_d_inner_mult(p, mlprec_wrk, level, trans, work)        

      case(mld_twoside_smooth_)
        call mld_d_inner_mult(p, mlprec_wrk, level, trans, work)        

      case default
        info = psb_err_from_subroutine_ai_
        call psb_errpush(info,name,a_err='invalid smooth_pos',&
             &  i_Err=(/p%precv(level)%parms%smoother_pos,izero,izero,izero,izero/))
        goto 9999      

      end select


    case(mld_mult_dev_ml_) 

      call mld_d_inner_mult(p, mlprec_wrk, level, trans, work)


    case(mld_vcycle_ml_, mld_wcycle_ml_)

      call mld_d_inner_vw_cycle(p, mlprec_wrk, level, trans, work, u=u)

    case(mld_kcycle_ml_, mld_kcyclesym_ml_)

      call mld_d_inner_k_cycle(p, mlprec_wrk, level, trans, work, u=u)
      
    case default
      info = psb_err_from_subroutine_ai_
      call psb_errpush(info,name,a_err='invalid mltype',&
           &  i_Err=(/p%precv(level)%parms%ml_type,izero,izero,izero,izero/))
      goto 9999      

    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine inner_ml_aply


  recursive subroutine mld_d_inner_add(p, mlprec_wrk, level, trans, work)
    use psb_base_mod
    use mld_prec_mod

    implicit none

    !Input/Oputput variables
    type(mld_dprec_type), intent(inout)  :: p

    type(mld_mlprec_wrk_type), target,  intent(inout) :: mlprec_wrk(:)
    integer(psb_ipk_), intent(in) :: level
    character, intent(in)             :: trans
    real(psb_dpk_),target            :: work(:)
    type(psb_d_vect_type) :: res
    type(psb_d_vect_type), pointer :: current
    integer(psb_ipk_) :: sweeps_post, sweeps_pre
    ! Local variables
    integer(psb_ipk_)  :: ictxt,np,me
    integer(psb_ipk_)  :: i, nr2l,nc2l,err_act
    integer(psb_ipk_)  :: debug_level, debug_unit
    integer(psb_ipk_)  :: nlev, ilev, sweeps
    logical            :: pre, post
    character(len=20)  :: name



    name = 'inner_inner_add'
    info = psb_success_
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    nlev = size(p%precv)
    if ((level < 1) .or. (level > nlev)) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='wrong call level to inner_add')
      goto 9999      
    end if
    ictxt = p%precv(level)%base_desc%get_context()
    call psb_info(ictxt, me, np)

    nc2l  = p%precv(level)%base_desc%get_local_cols()
    nr2l  = p%precv(level)%base_desc%get_local_rows()
    if(debug_level > 1) then
      write(debug_unit,*) me,' inner_add at level ',level
    end if

    if ((level<1).or.(level>nlev)) then
      info = psb_err_internal_error_ 
      call psb_errpush(info,name,&
           & a_err='Invalid LEVEL>NLEV')
      goto 9999
    end if
    
    
    if (level > 1) then 
      ! Apply the restriction
      call psb_map_X2Y(done,mlprec_wrk(level-1)%vx2l,&
           & dzero,mlprec_wrk(level)%vx2l,&
           & p%precv(level)%map,info,work=work)

      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during restriction')
        goto 9999
      end if

    end if
    
    sweeps = p%precv(level)%parms%sweeps 
    call p%precv(level)%sm%apply(done,&
         & mlprec_wrk(level)%vx2l,dzero,mlprec_wrk(level)%vy2l,&
         & p%precv(level)%base_desc, trans,&
         & sweeps,work,info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Error during ADD smoother_apply')
      goto 9999
    end if

    if (level < nlev) then
      call inner_ml_aply(level+1,p,mlprec_wrk,trans,work,info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error in recursive call')
        goto 9999
      end if

      !
      ! Apply the prolongator
      !  
      call psb_map_Y2X(done,mlprec_wrk(level+1)%vy2l,&
           & done,mlprec_wrk(level)%vy2l,&
           & p%precv(level+1)%map,info,work=work)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during prolongation')
        goto 9999
      end if


    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine mld_d_inner_add

  recursive subroutine mld_d_inner_mult(p, mlprec_wrk, level, trans, work)
    use psb_base_mod
    use mld_prec_mod

    implicit none

    !Input/Oputput variables
    type(mld_dprec_type), intent(inout)  :: p

    type(mld_mlprec_wrk_type), target,  intent(inout) :: mlprec_wrk(:)
    integer(psb_ipk_), intent(in) :: level
    character, intent(in)             :: trans
    real(psb_dpk_),target            :: work(:)
    type(psb_d_vect_type) :: res
    type(psb_d_vect_type), pointer :: current
    integer(psb_ipk_) :: sweeps_post, sweeps_pre
    ! Local variables
    integer(psb_ipk_)  :: ictxt,np,me
    integer(psb_ipk_)  :: i, nr2l,nc2l,err_act
    integer(psb_ipk_)  :: debug_level, debug_unit
    integer(psb_ipk_)  :: nlev, ilev, sweeps
    logical            :: pre, post
    character(len=20)  :: name



    name = 'inner_inner_mult'
    info = psb_success_
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    nlev = size(p%precv)
    if ((level < 1) .or. (level > nlev)) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='wrong call level to inner_mult')
      goto 9999      
    end if
    ictxt = p%precv(level)%base_desc%get_context()
    call psb_info(ictxt, me, np)

    nc2l  = p%precv(level)%base_desc%get_local_cols()
    nr2l  = p%precv(level)%base_desc%get_local_rows()
    if(debug_level > 1) then
      write(debug_unit,*) me,' inner_mult at level ',level
    end if


    if ((level < nlev).or.(nlev == 1)) then
      sweeps_post = p%precv(level)%parms%sweeps_post
      sweeps_pre  = p%precv(level)%parms%sweeps_pre
    else
      sweeps_post = p%precv(level-1)%parms%sweeps_post
      sweeps_pre  = p%precv(level-1)%parms%sweeps_pre
    endif

    pre  = ((sweeps_pre>0).and.(trans=='N')).or.((sweeps_post>0).and.(trans/='N'))
    post = ((sweeps_post>0).and.(trans=='N')).or.((sweeps_pre>0).and.(trans/='N'))


    if (level > 1) then 
      ! Apply the restriction
      if (pre) then
        current => mlprec_wrk(level-1)%vty
      else
        current => mlprec_wrk(level-1)%vx2l
      endif
      call psb_map_X2Y(done,current,&
           & dzero,mlprec_wrk(level)%vx2l,&
           & p%precv(level)%map,info,work=work)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during restriction')
        goto 9999
      end if
    end if

    
    if (level < nlev) then 

      !
      ! Apply the base preconditioner
      !

      if (pre) then
        if (trans == 'N') then 
          sweeps = p%precv(level)%parms%sweeps_pre
          if (info == psb_success_) call p%precv(level)%sm%apply(done,&
               & mlprec_wrk(level)%vx2l,dzero,mlprec_wrk(level)%vy2l,&
               & p%precv(level)%base_desc, trans,&
               & sweeps,work,info)
        else
          sweeps = p%precv(level)%parms%sweeps_post
          if (info == psb_success_) call p%precv(level)%sm2%apply(done,&
               & mlprec_wrk(level)%vx2l,dzero,mlprec_wrk(level)%vy2l,&
               & p%precv(level)%base_desc, trans,&
               & sweeps,work,info)
        end if
        
        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='Error during 2-PRE smoother_apply')
          goto 9999
        end if
      endif

      !
      ! Compute the residual and call recursively
      !
      if (pre) then
        call psb_geaxpby(done,mlprec_wrk(level)%vx2l,&
             & dzero,mlprec_wrk(level)%vty,&
             & p%precv(level)%base_desc,info)
        
        if (info == psb_success_) call psb_spmm(-done,p%precv(level)%base_a,&
             & mlprec_wrk(level)%vy2l,done,mlprec_wrk(level)%vty,&
             & p%precv(level)%base_desc,info,work=work,trans=trans)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='Error during residue')
          goto 9999
        end if
      endif
      
      call inner_ml_aply(level+1,p,mlprec_wrk,trans,work,info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error in recursive call')
        goto 9999
      end if
      
      
      !
      ! Apply the prolongator
      !  
      
      if (pre) then
        par = done
      else
        par = dzero
      endif
      
      call psb_map_Y2X(done,mlprec_wrk(level+1)%vy2l,&
           & par,mlprec_wrk(level)%vy2l,&
           & p%precv(level+1)%map,info,work=work)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during prolongation')
        goto 9999
      end if
      
      !
      ! Compute the residual
      !
      if (post) then
        call psb_geaxpby(done,mlprec_wrk(level)%vx2l,&
             & dzero,mlprec_wrk(level)%vtx,&
             & p%precv(level)%base_desc,info)
        call psb_spmm(-done,p%precv(level)%base_a,mlprec_wrk(level)%vy2l,&
             & done,mlprec_wrk(level)%vtx,p%precv(level)%base_desc,info,&
             & work=work,trans=trans)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='Error during residue')
          goto 9999
        end if
        !
        ! Apply the base preconditioner
        !
        if (trans == 'N') then
          sweeps = p%precv(level)%parms%sweeps_post
          if (info == psb_success_) call p%precv(level)%sm2%apply(done,&
               & mlprec_wrk(level)%vtx,done,mlprec_wrk(level)%vy2l,&
               & p%precv(level)%base_desc, trans,&
               & sweeps,work,info)
        else 
          sweeps = p%precv(level)%parms%sweeps_pre
          if (info == psb_success_) call p%precv(level)%sm%apply(done,&
               & mlprec_wrk(level)%vtx,done,mlprec_wrk(level)%vy2l,&
               & p%precv(level)%base_desc, trans,&
               & sweeps,work,info)
        end if
        
        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='Error during 2-POST smoother_apply')
          goto 9999
        end if
      endif            
      
    else if (level == nlev) then
      
      sweeps = p%precv(level)%parms%sweeps
      if (info == psb_success_) call p%precv(level)%sm%apply(done,&
           & mlprec_wrk(level)%vx2l,dzero,mlprec_wrk(level)%vy2l,&
           & p%precv(level)%base_desc, trans,&
           & sweeps,work,info)

    else

      info = psb_err_internal_error_ 
      call psb_errpush(info,name,&
           & a_err='Invalid LEVEL>NLEV')
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine mld_d_inner_mult
    
  recursive subroutine mld_d_inner_vw_cycle(p, mlprec_wrk, level, trans, work,u)
    use psb_base_mod
    use mld_prec_mod

    implicit none

    !Input/Oputput variables
    type(mld_dprec_type), intent(inout)  :: p
    type(mld_mlprec_wrk_type), target,  intent(inout) :: mlprec_wrk(:)
    integer(psb_ipk_), intent(in) :: level
    character, intent(in)             :: trans
    real(psb_dpk_),target            :: work(:)
    type(psb_d_vect_type),intent(inout), optional :: u



    type(psb_d_vect_type) :: res
    type(psb_d_vect_type), pointer :: current
    integer(psb_ipk_) :: sweeps_post, sweeps_pre
    ! Local variables
    integer(psb_ipk_)  :: ictxt,np,me
    integer(psb_ipk_)  :: i, nr2l,nc2l,err_act
    integer(psb_ipk_)  :: debug_level, debug_unit
    integer(psb_ipk_)  :: nlev, ilev, sweeps
    logical            :: pre, post
    character(len=20)  :: name



    name = 'inner_inner_add'
    info = psb_success_
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    nlev = size(p%precv)
    if ((level < 1) .or. (level > nlev)) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='wrong call level to inner_add')
      goto 9999      
    end if
    ictxt = p%precv(level)%base_desc%get_context()
    call psb_info(ictxt, me, np)

    nc2l  = p%precv(level)%base_desc%get_local_cols()
    nr2l  = p%precv(level)%base_desc%get_local_rows()
    if(debug_level > 1) then
      write(debug_unit,*) me,' inner_add at level ',level
    end if

    if ((level<1).or.(level>nlev)) then
      info = psb_err_internal_error_ 
      call psb_errpush(info,name,&
           & a_err='Invalid LEVEL>NLEV')
      goto 9999
    end if
    call psb_geasb(res,p%precv(level)%base_desc,info,&
         & scratch=.true., mold=mlprec_wrk(level)%vx2l%v)

    !V/W cycle
    if (level > 1) then 
      ! Apply the restriction
      call psb_map_X2Y(done,mlprec_wrk(level-1)%vty,&
           & dzero,mlprec_wrk(level)%vx2l,&
           & p%precv(level)%map,info,work=work)

      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during restriction')
        goto 9999
      end if
    end if

    call psb_geaxpby(done,mlprec_wrk(level)%vx2l,&
         & dzero,mlprec_wrk(level)%vtx,&
         & p%precv(level)%base_desc,info)
    !
    ! Apply the base preconditioner
    !
    if (level < nlev) then 

      if (present(u)) then
        !    call mlprec_wrk(level)%vy2l%set(u%get_vect())
        call psb_geaxpby(done,u,&
             & dzero,mlprec_wrk(level)%vy2l,&
             & p%precv(level)%base_desc,info)

      else
        call mlprec_wrk(level)%vy2l%zero()
      endif
      call psb_geaxpby(done,mlprec_wrk(level)%vx2l,dzero,res,&
           & p%precv(level)%base_desc,info)

      call psb_spmm(-done,p%precv(level)%base_a,mlprec_wrk(level)%vy2l,& 
           done, res, p%precv(level)%base_desc, info, work=work, trans=trans) 

      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during residue')
        goto 9999
      end if

      if (trans == 'N') then 
        sweeps = p%precv(level)%parms%sweeps_pre
        if (info == psb_success_) call p%precv(level)%sm%apply(done,&
             & mlprec_wrk(level)%vx2l,dzero,mlprec_wrk(level)%vy2l,&
             & p%precv(level)%base_desc, trans,&
             & sweeps,work,info)
      else
        sweeps = p%precv(level)%parms%sweeps_post
        if (info == psb_success_) call p%precv(level)%sm2%apply(done,&
             & mlprec_wrk(level)%vx2l,dzero,mlprec_wrk(level)%vy2l,&
             & p%precv(level)%base_desc, trans,&
             & sweeps,work,info)
      end if
    else
      sweeps = p%precv(level)%parms%sweeps
      if (info == psb_success_) call p%precv(level)%sm%apply(done,&
           & mlprec_wrk(level)%vx2l,dzero,mlprec_wrk(level)%vy2l,&
           & p%precv(level)%base_desc, trans,&
           & sweeps,work,info)
    end if

    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Error during 2-PRE smoother_apply')
      goto 9999
    end if


    !
    ! Compute the residual (at all levels but the coarsest one)
    ! and call recursively
    !
    if(level < nlev) then

      call psb_geaxpby(done,mlprec_wrk(level)%vx2l,&
           & dzero,mlprec_wrk(level)%vty,&
           & p%precv(level)%base_desc,info)

      if (info == psb_success_) call psb_spmm(-done,p%precv(level)%base_a,&
           & mlprec_wrk(level)%vy2l,done,mlprec_wrk(level)%vty,&
           & p%precv(level)%base_desc,info,work=work,trans=trans)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during residue')
        goto 9999
      end if

      call inner_ml_aply(level+1,p,mlprec_wrk,trans,work,info)

      if (p%precv(level)%parms%ml_type == mld_wcycle_ml_) then
        call inner_ml_aply(level+1,p,mlprec_wrk,trans,work,info, u=mlprec_wrk(level+1)%vy2l)
      endif

      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error in recursive call')
        goto 9999
      end if


      !
      ! Apply the prolongator
      !  
      call psb_map_Y2X(done,mlprec_wrk(level+1)%vy2l,&
           & done,mlprec_wrk(level)%vy2l,&
           & p%precv(level+1)%map,info,work=work)

      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during prolongation')
        goto 9999
      end if

      !
      ! Compute the residual
      !
      call psb_spmm(-done,p%precv(level)%base_a,mlprec_wrk(level)%vy2l,&
           & done,mlprec_wrk(level)%vtx,p%precv(level)%base_desc,info,&
           & work=work,trans=trans)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during residue')
        goto 9999
      end if
      !
      ! Apply the base preconditioner
      !
      if (trans == 'N') then 
        sweeps = p%precv(level)%parms%sweeps_post
        if (info == psb_success_) call p%precv(level)%sm2%apply(done,&
             & mlprec_wrk(level)%vtx,done,mlprec_wrk(level)%vy2l,&
             & p%precv(level)%base_desc, trans,&
             & sweeps,work,info)
      else
        sweeps = p%precv(level)%parms%sweeps_pre
        if (info == psb_success_) call p%precv(level)%sm%apply(done,&
             & mlprec_wrk(level)%vtx,done,mlprec_wrk(level)%vy2l,&
             & p%precv(level)%base_desc, trans,&
             & sweeps,work,info)
      end if

      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during 2-POST smoother_apply')
        goto 9999
      end if

    endif


    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine mld_d_inner_vw_cycle

  recursive subroutine mld_d_inner_k_cycle(p, mlprec_wrk, level, trans, work,u)
    use psb_base_mod
    use mld_prec_mod

    implicit none

    !Input/Oputput variables
    type(mld_dprec_type), intent(inout)  :: p
    type(mld_mlprec_wrk_type), target,  intent(inout) :: mlprec_wrk(:)
    integer(psb_ipk_), intent(in) :: level
    character, intent(in)             :: trans
    real(psb_dpk_),target            :: work(:)
    type(psb_d_vect_type),intent(inout), optional :: u



    type(psb_d_vect_type) :: res
    type(psb_d_vect_type), pointer :: current
    integer(psb_ipk_) :: sweeps_post, sweeps_pre
    ! Local variables
    integer(psb_ipk_)  :: ictxt,np,me
    integer(psb_ipk_)  :: i, nr2l,nc2l,err_act
    integer(psb_ipk_)  :: debug_level, debug_unit
    integer(psb_ipk_)  :: nlev, ilev, sweeps
    logical            :: pre, post
    character(len=20)  :: name



    name = 'inner_inner_add'
    info = psb_success_
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    nlev = size(p%precv)
    if ((level < 1) .or. (level > nlev)) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='wrong call level to inner_add')
      goto 9999      
    end if
    ictxt = p%precv(level)%base_desc%get_context()
    call psb_info(ictxt, me, np)

    nc2l  = p%precv(level)%base_desc%get_local_cols()
    nr2l  = p%precv(level)%base_desc%get_local_rows()
    if(debug_level > 1) then
      write(debug_unit,*) me,' inner_add at level ',level
    end if

    if ((level<1).or.(level>nlev)) then
      info = psb_err_internal_error_ 
      call psb_errpush(info,name,&
           & a_err='Invalid LEVEL>NLEV')
      goto 9999
    end if

    !K cycle

    call psb_geaxpby(done,mlprec_wrk(level)%vx2l,&
         & dzero,mlprec_wrk(level)%vtx,&
         & p%precv(level)%base_desc,info)
    !
    ! Apply the base preconditioner
    !
    if (level < nlev) then 

      if (present(u)) then
        call psb_geaxpby(done,u,&
             & dzero,mlprec_wrk(level)%vy2l,&
             & p%precv(level)%base_desc,info)
      else
        call mlprec_wrk(level)%vy2l%zero()
      endif
      res = mlprec_wrk(level)%vx2l

      call psb_spmm(-done,p%precv(level)%base_a,mlprec_wrk(level)%vy2l,& 
           done, res, p%precv(level)%base_desc, info, work=work, trans=trans) 

      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during residue')
        goto 9999
      end if

      if (trans == 'N') then 
        sweeps = p%precv(level)%parms%sweeps_pre
        if (info == psb_success_) call p%precv(level)%sm%apply(done,&
             & mlprec_wrk(level)%vx2l,dzero,mlprec_wrk(level)%vy2l,&
             & p%precv(level)%base_desc, trans,&
             & sweeps,work,info)
      else
        sweeps = p%precv(level)%parms%sweeps_post
        if (info == psb_success_) call p%precv(level)%sm2%apply(done,&
             & mlprec_wrk(level)%vx2l,dzero,mlprec_wrk(level)%vy2l,&
             & p%precv(level)%base_desc, trans,&
             & sweeps,work,info)
      end if
    else
      sweeps = p%precv(level)%parms%sweeps
      if (info == psb_success_) call p%precv(level)%sm%apply(done,&
           & mlprec_wrk(level)%vx2l,dzero,mlprec_wrk(level)%vy2l,&
           & p%precv(level)%base_desc, trans,&
           & sweeps,work,info)
    end if

    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Error during 2-PRE smoother_apply')
      goto 9999
    end if


    !
    ! Compute the residual (at all levels but the coarsest one)
    ! and call recursively
    !
    if(level < nlev) then

      call psb_geaxpby(done,mlprec_wrk(level)%vx2l,&
           & dzero,mlprec_wrk(level)%vty,&
           & p%precv(level)%base_desc,info)

      if (info == psb_success_) call psb_spmm(-done,p%precv(level)%base_a,&
           & mlprec_wrk(level)%vy2l,done,mlprec_wrk(level)%vty,&
           & p%precv(level)%base_desc,info,work=work,trans=trans)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during residue')
        goto 9999
      end if

      ! Apply the restriction
      call psb_map_X2Y(done,mlprec_wrk(level)%vty,&
           & dzero,mlprec_wrk(level + 1)%vx2l,&
           & p%precv(level + 1)%map,info,work=work)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during restriction')
        goto 9999
      end if

      !Set the preconditioner


      if ((level < nlev - 2)) then
        if (p%precv(level)%parms%ml_type == mld_kcyclesym_ml_) then
          call mld_dinneritkcycle(p, mlprec_wrk, level + 1, trans, work, 'FCG')
        elseif (p%precv(level)%parms%ml_type == mld_kcycle_ml_) then
          call mld_dinneritkcycle(p, mlprec_wrk, level + 1, trans, work, 'CGR') 
        endif
      else
        call inner_ml_aply(level + 1 ,p,mlprec_wrk,trans,work,info)
      endif


      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error in recursive call')
        goto 9999
      end if


      !
      ! Apply the prolongator
      !  
      call psb_map_Y2X(done,mlprec_wrk(level+1)%vy2l,&
           & done,mlprec_wrk(level)%vy2l,&
           & p%precv(level+1)%map,info,work=work)

      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during prolongation')
        goto 9999
      end if

      !
      ! Compute the residual
      !
      call psb_spmm(-done,p%precv(level)%base_a,mlprec_wrk(level)%vy2l,&
           & done,mlprec_wrk(level)%vtx,p%precv(level)%base_desc,info,&
           & work=work,trans=trans)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during residue')
        goto 9999
      end if
      !
      ! Apply the base preconditioner
      !
      if (trans == 'N') then 
        sweeps = p%precv(level)%parms%sweeps_post
        if (info == psb_success_) call p%precv(level)%sm2%apply(done,&
             & mlprec_wrk(level)%vtx,done,mlprec_wrk(level)%vy2l,&
             & p%precv(level)%base_desc, trans,&
             & sweeps,work,info)
      else
        sweeps = p%precv(level)%parms%sweeps_pre
        if (info == psb_success_) call p%precv(level)%sm%apply(done,&
             & mlprec_wrk(level)%vtx,done,mlprec_wrk(level)%vy2l,&
             & p%precv(level)%base_desc, trans,&
             & sweeps,work,info)
      end if

      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during 2-POST smoother_apply')
        goto 9999
      end if

    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine mld_d_inner_k_cycle


  recursive subroutine mld_dinneritkcycle(p, mlprec_wrk, level, trans, work, innersolv)
    use psb_base_mod
    use mld_prec_mod
    use mld_d_inner_mod, mld_protect_name => mld_dmlprec_aply

    implicit none

    !Input/Oputput variables
    type(mld_dprec_type), intent(inout)  :: p

    type(mld_mlprec_wrk_type), intent(inout) :: mlprec_wrk(:)
    integer(psb_ipk_), intent(in) :: level
    character, intent(in)             :: trans, innersolv
    real(psb_dpk_),target            :: work(:)

    !Other variables
    type(psb_d_vect_type)  :: v, w, rhs, v1, x
    type(psb_d_vect_type), dimension(0:1) ::  d
    real(psb_dpk_) :: delta_old, rhs_norm, alpha, tau, tau1, tau2, tau3, tau4, beta

    real(psb_dpk_) :: l2_norm, delta, rtol=0.25
    real(psb_dpk_), allocatable :: temp_v(:)
    integer(psb_ipk_) :: info, nlev, i, iter, max_iter=2, idx

    !Assemble rhs, w, v, v1, x

    call psb_geasb(rhs,&
         & p%precv(level)%base_desc,info,&
         & scratch=.true.,mold=mlprec_wrk(level)%vx2l%v)
    call psb_geasb(w,&
         & p%precv(level)%base_desc,info,&
         & scratch=.true.,mold=mlprec_wrk(level)%vx2l%v)
    call psb_geasb(v,&
         & p%precv(level)%base_desc,info,&
         & scratch=.true.,mold=mlprec_wrk(level)%vx2l%v)
    call psb_geasb(v1,&
         & p%precv(level)%base_desc,info,&
         & scratch=.true.,mold=mlprec_wrk(level)%vx2l%v)
    call psb_geasb(x,&
         & p%precv(level)%base_desc,info,&
         & scratch=.true.,mold=mlprec_wrk(level)%vx2l%v)

    call x%set(dzero)

    ! rhs=vx2l and w=rhs
    call psb_geaxpby(done,mlprec_wrk(level)%vx2l,dzero,rhs,&
         &   p%precv(level)%base_desc,info)
    call psb_geaxpby(done,mlprec_wrk(level)%vx2l,dzero,w,&
         &   p%precv(level)%base_desc,info)

    if (psb_errstatus_fatal()) then 
      nc2l = p%precv(level)%base_desc%get_local_cols()
      info=psb_err_alloc_request_
      call psb_errpush(info,name,i_err=(/2*nc2l,izero,izero,izero,izero/),&
           & a_err='TYPE@(psb_dpk_)')
      goto 9999      
    end if

    delta = psb_gedot(w, w, p%precv(level)%base_desc, info)

    !Apply the preconditioner

    call mlprec_wrk(level)%vy2l%set(dzero) 

    idx=0
    call inner_ml_aply(level,p,mlprec_wrk,trans,work,info)

    !Assemble d(0) and d(1)
    call psb_geasb(d(0),&
         & p%precv(level)%base_desc,info,&
         & scratch=.true.,mold=mlprec_wrk(level)%vy2l%v)
    call psb_geasb(d(1),&
         & p%precv(level)%base_desc,info,&
         & scratch=.true.,mold=mlprec_wrk(level)%vy2l%v)

    call psb_geaxpby(done,mlprec_wrk(level)%vy2l,dzero,d(idx),p%precv(level)%base_desc,info)


    call psb_spmm(done,p%precv(level)%base_a,d(idx),dzero,v,p%precv(level)%base_desc,info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Error during residue')
      goto 9999
    end if

    !FCG
    if (innersolv == 'FCG') then
      delta_old = psb_gedot(d(idx), w, p%precv(level)%base_desc, info)
      tau = psb_gedot(d(idx), v, p%precv(level)%base_desc, info) 
      !CGR
    else
      delta_old = psb_gedot(v, w, p%precv(level)%base_desc, info)
      tau = psb_gedot(v, v, p%precv(level)%base_desc, info) 
    endif

    alpha = delta_old/tau
    !Update residual w
    call psb_geaxpby(-alpha, v, done, w, p%precv(level)%base_desc, info) 

    l2_norm = psb_gedot(w, w, p%precv(level)%base_desc, info)
    iter = 0 


    if (l2_norm <= rtol*delta) then
      !Update solution x
      call psb_geaxpby(alpha, d(idx), done, x, p%precv(level)%base_desc, info)   
    else
      iter = iter + 1
      idx=mod(iter,2)

      !Apply preconditioner
      call psb_geaxpby(done,w,dzero,mlprec_wrk(level)%vx2l,p%precv(level)%base_desc,info)    
      call inner_ml_aply(level,p,mlprec_wrk,trans,work,info)
      call psb_geaxpby(done,mlprec_wrk(level)%vy2l,dzero,d(idx),p%precv(level)%base_desc,info)

      !Sparse matrix vector product

      call psb_spmm(done,p%precv(level)%base_a,d(idx),dzero,v1,p%precv(level)%base_desc,info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during residue')
        goto 9999
      end if

      !tau1, tau2, tau3, tau4
      !FCG
      if (innersolv == 'FCG') then
        tau1= psb_gedot(d(idx), v, p%precv(level)%base_desc, info)
        tau2= psb_gedot(d(idx), v1, p%precv(level)%base_desc, info)
        tau3= psb_gedot(d(idx), w, p%precv(level)%base_desc, info)
        tau4= tau2 - (tau1*tau1)/tau
        !CGR
      else
        tau1= psb_gedot(v1, v, p%precv(level)%base_desc, info)
        tau2= psb_gedot(v1, v1, p%precv(level)%base_desc, info)
        tau3= psb_gedot(v1, w, p%precv(level)%base_desc, info)
        tau4= tau2 - (tau1*tau1)/tau
      endif

      !Update solution
      alpha=alpha-(tau1*tau3)/(tau*tau4)
      call psb_geaxpby(alpha,d(idx - 1),done,x,p%precv(level)%base_desc,info)
      alpha=tau3/tau4
      call psb_geaxpby(alpha,d(idx),done,x,p%precv(level)%base_desc,info)
    endif

    !Free vectors
    call psb_geaxpby(done,x,dzero,mlprec_wrk(level)%vy2l,p%precv(level)%base_desc,info)
    call psb_gefree(v, p%precv(level)%base_desc, info)
    call psb_gefree(v1, p%precv(level)%base_desc, info)
    call psb_gefree(w, p%precv(level)%base_desc, info)
    call psb_gefree(x, p%precv(level)%base_desc, info)
    call psb_gefree(d(0), p%precv(level)%base_desc, info)
    call psb_gefree(d(1), p%precv(level)%base_desc, info)

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine mld_dinneritkcycle

end subroutine mld_dmlprec_aply_vect

