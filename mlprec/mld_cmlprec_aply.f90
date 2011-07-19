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
!  alpha       -   complex(psb_spk_), input.
!                  The scalar alpha.
!   p          -   type(mld_cprec_type), input.
!                  The multilevel preconditioner data structure containing the
!                  local part of the preconditioner to be applied.
!      Note that nlev = size(p%precv) = number of levels.
!      p%precv(ilev)%prec      -  type(psb_cbaseprec_type)
!                                 The 'base preconditioner' for the current level
!      p%precv(ilev)%ac        -  type(psb_cspmat_type) 
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
!      p%precv(ilev)%rprcparm  -  real(psb_spk_), dimension(:), allocatable.
!                                 The real parameters defining the multilevel strategy
!      p%precv(ilev)%mlia      -  integer, dimension(:), allocatable.
!                                 The aggregation map (ilev-1) --> (ilev).
!      p%precv(ilev)%nlaggr    -  integer, dimension(:), allocatable.
!                                 The number of aggregates (rows of A(ilev)) on the
!                                 various processes. 
!      p%precv(ilev)%base_a    -  type(psb_cspmat_type), pointer.
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
!
! Additive multilevel
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
subroutine mld_cmlprec_aply(alpha,p,x,beta,y,desc_data,trans,work,info)

  use psb_base_mod
  use mld_c_inner_mod, mld_protect_name => mld_cmlprec_aply

  implicit none

  ! Arguments
  type(psb_desc_type),intent(in)    :: desc_data
  type(mld_cprec_type), intent(in)  :: p
  complex(psb_spk_),intent(in)      :: alpha,beta
  complex(psb_spk_),intent(inout)   :: x(:)
  complex(psb_spk_),intent(inout)   :: y(:)
  character, intent(in)             :: trans
  complex(psb_spk_),target          :: work(:)
  integer, intent(out)              :: info

  ! Local variables
  integer           :: ictxt, np, me, err_act
  integer           :: debug_level, debug_unit, nlev,nc2l,nr2l,level
  character(len=20) :: name
  character         :: trans_
  type psb_mlprec_wrk_type
    complex(psb_spk_), allocatable  :: tx(:), ty(:), x2l(:), y2l(:)
  end type psb_mlprec_wrk_type
  type(psb_mlprec_wrk_type), allocatable  :: mlprec_wrk(:)

  name = 'mld_cmlprec_aply'
  info = psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = psb_cd_get_context(desc_data)
  call psb_info(ictxt, me, np)

  if (debug_level >= psb_debug_inner_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ' Entry  ', size(p%precv)

  trans_ = psb_toupper(trans)

  nlev = size(p%precv)
  allocate(mlprec_wrk(nlev),stat=info) 
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,&
         & a_err='Allocate')
    goto 9999      
  end if
  level = 1

  nc2l  = psb_cd_get_local_cols(p%precv(level)%base_desc)
  nr2l  = psb_cd_get_local_rows(p%precv(level)%base_desc)
  allocate(mlprec_wrk(level)%x2l(nc2l),mlprec_wrk(level)%y2l(nc2l),&
       & stat=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/size(x)+size(y),0,0,0,0/),&
         & a_err='real(psb_spk_)')
    goto 9999      
  end if

  mlprec_wrk(level)%x2l(:) = x(:) 
  mlprec_wrk(level)%y2l(:) = czero 

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

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

contains

  recursive subroutine inner_ml_aply(level,p,mlprec_wrk,trans,work,info)    

    implicit none 

    ! Arguments
    integer                          :: level 
    type(mld_cprec_type), intent(in) :: p
    type(psb_mlprec_wrk_type), intent(inout) :: mlprec_wrk(:)
    character, intent(in)            :: trans
    complex(psb_spk_),target         :: work(:)
    integer, intent(out)             :: info

    ! Local variables
    integer            :: ictxt,np,me,i, nr2l,nc2l,err_act
    integer            :: debug_level, debug_unit
    integer            :: nlev, ilev, sweeps

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
    ictxt = psb_cd_get_context(p%precv(level)%base_desc)
    call psb_info(ictxt, me, np)


    if (level > 1) then 
      nc2l  = psb_cd_get_local_cols(p%precv(level)%base_desc)
      nr2l  = psb_cd_get_local_rows(p%precv(level)%base_desc)
      allocate(mlprec_wrk(level)%x2l(nc2l),mlprec_wrk(level)%y2l(nc2l),&
           & stat=info)
      if (info /= psb_success_) then 
        info=psb_err_alloc_request_
        call psb_errpush(info,name,i_err=(/2*nc2l,0,0,0,0/),&
             & a_err='complex(psb_spk_)')
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
        call psb_map_X2Y(cone,mlprec_wrk(level-1)%x2l,&
             & czero,mlprec_wrk(level)%x2l,&
             & p%precv(level)%map,info,work=work)

        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='Error during restriction')
          goto 9999
        end if

      end if

      sweeps = p%precv(level)%parms%sweeps 
      call p%precv(level)%sm%apply(cone,&
           & mlprec_wrk(level)%x2l,czero,mlprec_wrk(level)%y2l,&
           & p%precv(level)%base_desc, trans,&
           & sweeps,work,info)

      if (info /= psb_success_) goto 9999
      if (level < nlev) then
        call inner_ml_aply(level+1,p,mlprec_wrk,trans,work,info)
        if (info /= psb_success_) goto 9999      
        !
        ! Apply the prolongator
        !  
        call psb_map_Y2X(cone,mlprec_wrk(level+1)%y2l,&
             & cone,mlprec_wrk(level)%y2l,&
             & p%precv(level+1)%map,info,work=work)
        if (info /= psb_success_) goto 9999

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
            call psb_map_X2Y(cone,mlprec_wrk(level-1)%x2l,&
                 & czero,mlprec_wrk(level)%x2l,&
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
            if (info /= psb_success_) goto 9999      
            !
            ! Apply the prolongator
            !  
            call psb_map_Y2X(cone,mlprec_wrk(level+1)%y2l,&
                 & czero,mlprec_wrk(level)%y2l,&
                 & p%precv(level+1)%map,info,work=work)
            if (info /= psb_success_) goto 9999
            !
            ! Compute the residual
            !
            call psb_spmm(-cone,p%precv(level)%base_a,mlprec_wrk(level)%y2l,&
                 & cone,mlprec_wrk(level)%x2l,p%precv(level)%base_desc,info,&
                 & work=work,trans=trans)
            if (info /= psb_success_) goto 9999

            sweeps = p%precv(level)%parms%sweeps_post 
            call p%precv(level)%sm%apply(cone,&
                 & mlprec_wrk(level)%x2l,cone,mlprec_wrk(level)%y2l,&
                 & p%precv(level)%base_desc, trans,&
                 & sweeps,work,info)
          else
            sweeps = p%precv(level)%parms%sweeps 
            call p%precv(level)%sm%apply(cone,&
                 & mlprec_wrk(level)%x2l,czero,mlprec_wrk(level)%y2l,&
                 & p%precv(level)%base_desc, trans,&
                 & sweeps,work,info)

          end if

        case('T','C')

          ! Post-smoothing transpose is pre-smoothing


          if (level > 1) then 
            ! Apply the restriction
            call psb_map_X2Y(cone,mlprec_wrk(level-1)%x2l,&
                 & czero,mlprec_wrk(level)%x2l,&
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
          else
            sweeps = p%precv(level)%parms%sweeps
          end if
          call p%precv(level)%sm%apply(cone,&
               & mlprec_wrk(level)%x2l,czero,mlprec_wrk(level)%y2l,&
               & p%precv(level)%base_desc, trans,&
               & sweeps,work,info)

          if (info /= psb_success_) goto 9999

          !
          ! Compute the residual (at all levels but the coarsest one)
          !
          if (level < nlev) then
            call psb_spmm(-cone,p%precv(level)%base_a,&
                 & mlprec_wrk(level)%y2l,cone,mlprec_wrk(level)%x2l,&
                 & p%precv(level)%base_desc,info,work=work,trans=trans)
            if (info /= psb_success_) goto 9999
            call inner_ml_aply(level+1,p,mlprec_wrk,trans,work,info) 
            if (info /= psb_success_) goto 9999

            call psb_map_Y2X(cone,mlprec_wrk(level+1)%y2l,&
                 & cone,mlprec_wrk(level)%y2l,&
                 & p%precv(level+1)%map,info,work=work)
            if (info /= psb_success_) goto 9999

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
            call psb_map_X2Y(cone,mlprec_wrk(level-1)%x2l,&
                 & czero,mlprec_wrk(level)%x2l,&
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
          call p%precv(level)%sm%apply(cone,&
               & mlprec_wrk(level)%x2l,czero,mlprec_wrk(level)%y2l,&
               & p%precv(level)%base_desc, trans,&
               & sweeps,work,info)

          if (info /= psb_success_) goto 9999

          !
          ! Compute the residual (at all levels but the coarsest one)
          !
          if (level < nlev) then
            call psb_spmm(-cone,p%precv(level)%base_a,&
                 & mlprec_wrk(level)%y2l,cone,mlprec_wrk(level)%x2l,&
                 & p%precv(level)%base_desc,info,work=work,trans=trans)
            if (info /= psb_success_) goto 9999
            call inner_ml_aply(level+1,p,mlprec_wrk,trans,work,info) 
            if (info /= psb_success_) goto 9999

            call psb_map_Y2X(cone,mlprec_wrk(level+1)%y2l,&
                 & cone,mlprec_wrk(level)%y2l,&
                 & p%precv(level+1)%map,info,work=work)
            if (info /= psb_success_) goto 9999

          end if


        case('T','C')

          ! pre-smooth transpose is post-smoothing 


          if (level > 1) then 
            ! Apply the restriction
            call psb_map_X2Y(cone,mlprec_wrk(level-1)%x2l,&
                 & czero,mlprec_wrk(level)%x2l,&
                 & p%precv(level)%map,info,work=work)

            if (info /= psb_success_) then
              call psb_errpush(psb_err_internal_error_,name,&
                   & a_err='Error during restriction')
              goto 9999
            end if

          end if

          if (level < nlev) then 
            call inner_ml_aply(level+1,p,mlprec_wrk,trans,work,info) 
            if (info /= psb_success_) goto 9999      
            !
            ! Apply the prolongator
            !  
            call psb_map_Y2X(cone,mlprec_wrk(level+1)%y2l,&
                 & czero,mlprec_wrk(level)%y2l,&
                 & p%precv(level+1)%map,info,work=work)
            if (info /= psb_success_) goto 9999
            !
            ! Compute the residual
            !
            call psb_spmm(-cone,p%precv(level)%base_a,mlprec_wrk(level)%y2l,&
                 & cone,mlprec_wrk(level)%x2l,p%precv(level)%base_desc,info,&
                 & work=work,trans=trans)
            if (info /= psb_success_) goto 9999

            sweeps = p%precv(level)%parms%sweeps_pre 
            call p%precv(level)%sm%apply(cone,&
                 & mlprec_wrk(level)%x2l,cone,mlprec_wrk(level)%y2l,&
                 & p%precv(level)%base_desc, trans,&
                 & sweeps,work,info)
          else
            sweeps = p%precv(level)%parms%sweeps 
            call p%precv(level)%sm%apply(cone,&
                 & mlprec_wrk(level)%x2l,czero,mlprec_wrk(level)%y2l,&
                 & p%precv(level)%base_desc, trans,&
                 & sweeps,work,info)
          end if

        case default
          info = psb_err_internal_error_
          call psb_errpush(info,name,a_err='invalid trans')
          goto 9999      
        end select

      case(mld_twoside_smooth_)

        nc2l  = psb_cd_get_local_cols(p%precv(level)%base_desc)
        nr2l  = psb_cd_get_local_rows(p%precv(level)%base_desc)
        allocate(mlprec_wrk(level)%ty(nc2l), mlprec_wrk(level)%tx(nc2l), stat=info)
        if (info /= psb_success_) then 
          info=psb_err_alloc_request_
          call psb_errpush(info,name,i_err=(/2*nc2l,0,0,0,0/),&
               & a_err='real(psb_spk_)')
          goto 9999      
        end if

        if (level > 1) then 
          ! Apply the restriction
          call psb_map_X2Y(cone,mlprec_wrk(level-1)%ty,&
               & czero,mlprec_wrk(level)%x2l,&
               & p%precv(level)%map,info,work=work)

          if (info /= psb_success_) then
            call psb_errpush(psb_err_internal_error_,name,&
                 & a_err='Error during restriction')
            goto 9999
          end if
        end if
        call psb_geaxpby(cone,mlprec_wrk(level)%x2l,czero,&
             & mlprec_wrk(level)%tx,&
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
        if (info == psb_success_) call p%precv(level)%sm%apply(cone,&
             & mlprec_wrk(level)%x2l,czero,mlprec_wrk(level)%y2l,&
             & p%precv(level)%base_desc, trans,&
             & sweeps,work,info)
        !
        ! Compute the residual (at all levels but the coarsest one)
        ! and call recursively
        !
        if(level < nlev) then
          mlprec_wrk(level)%ty = mlprec_wrk(level)%x2l
          if (info == psb_success_) call psb_spmm(-cone,p%precv(level)%base_a,&
               & mlprec_wrk(level)%y2l,cone,mlprec_wrk(level)%ty,&
               & p%precv(level)%base_desc,info,work=work,trans=trans)

          call inner_ml_aply(level+1,p,mlprec_wrk,trans,work,info)


          !
          ! Apply the prolongator
          !  
          call psb_map_Y2X(cone,mlprec_wrk(level+1)%y2l,&
               & cone,mlprec_wrk(level)%y2l,&
               & p%precv(level+1)%map,info,work=work)

          if (info /= psb_success_ ) then
            call psb_errpush(psb_err_internal_error_,name,&
                 & a_err='Error during restriction')
            goto 9999
          end if

          !
          ! Compute the residual
          !
          call psb_spmm(-cone,p%precv(level)%base_a,mlprec_wrk(level)%y2l,&
               & cone,mlprec_wrk(level)%tx,p%precv(level)%base_desc,info,&
               & work=work,trans=trans)
          !
          ! Apply the base preconditioner
          !
          if (trans == 'N') then 
            sweeps = p%precv(level)%parms%sweeps_post
          else
            sweeps = p%precv(level)%parms%sweeps_pre
          end if
          if (info == psb_success_) call p%precv(level)%sm%apply(cone,&
               & mlprec_wrk(level)%tx,cone,mlprec_wrk(level)%y2l,&
               & p%precv(level)%base_desc, trans,&
               & sweeps,work,info)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_internal_error_,name,&
                 & a_err='Error: residual/baseprec_aply')
            goto 9999
          end if

        endif

      case default
        info = psb_err_from_subroutine_ai_
        call psb_errpush(info,name,a_err='invalid smooth_pos',&
             &  i_Err=(/p%precv(level)%parms%smoother_pos,0,0,0,0/))
        goto 9999      

      end select

    case default
      info = psb_err_from_subroutine_ai_
      call psb_errpush(info,name,a_err='invalid mltype',&
           &  i_Err=(/p%precv(level)%parms%ml_type,0,0,0,0/))
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

  end subroutine inner_ml_aply

end subroutine mld_cmlprec_aply

