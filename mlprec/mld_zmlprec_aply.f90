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
! File: mld_zmlprec_aply.f90
!
! Subroutine: mld_zmlprec_aply
! Version:    complex
!
!  This routine computes
!  
!                        Y = beta*Y + alpha*op(M^(-1))*X,
!  where 
!  - M is a multilevel domain decomposition (Schwarz) preconditioner associated
!    to a certain matrix A and stored in the array baseprecv,
!  - op(M^(-1)) is M^(-1) or its transpose, according to the value of trans,
!  - X and Y are vectors,
!  - alpha and beta are scalars.
!
!  For each level we have as many subdomains as processes (except for the coarsest
!  level where we might have a replicated index space) and each process      takes care
!  of one subdomain.
!
!  The multilevel preconditioner M is regarded as an array of 'base preconditioners',
!  each representing the part of the preconditioner associated to a certain level.
!  For each level ilev, the base preconditioner K(ilev) is stored in baseprecv(ilev)
!  and is associated to a matrix A(ilev), obtained by 'tranferring' the original
!  matrix A (i.e. the matrix to be preconditioned) to level ilev, through smoothed
!  aggregation.
!
!  The levels are numbered in increasing order starting from the finest one, i.e.
!  level 1 is the finest level and A(1) is the matrix A.
!
!  For a general description of (parallel) multilevel preconditioners see
!    1. B.F. Smith, P.E. Bjorstad & W.D. Gropp,
!       Domain decomposition: parallel multilevel methods for elliptic partial
!       differential equations,
!       Cambridge University Press, 1996.
!    2. K. Stuben,
!       Algebraic Multigrid (AMG): An Introduction with Applications,
!       GMD Report N. 70, 1999.
!
!
! Arguments:
!  alpha       -   complex(kind(0.d0)), input.
!                  The scalar alpha.
!   baseprecv  -   type(mld_zbaseprc_type), dimension(:), input.
!                  The array of base preconditioner data structures containing the
!                  local parts of the preconditioners to be applied at each level.
!      Note that nlev = size(baseprecv) = number of levels.
!      baseprecv(ilev)%av  -  type(psb_zspmat_type), dimension(:), allocatable(:).
!                             The sparse matrices needed to apply the preconditioner 
!                             at level ilev. 
!        baseprecv(ilev)%av(mld_l_pr_)    -  The L factor of the ILU factorization of the 
!                                            local diagonal block of A(ilev).
!        baseprecv(ilev)%av(mld_u_pr_)    -  The U factor of the ILU factorization of the
!                                            local diagonal block of A(ilev), except its
!                                            diagonal entries (stored in baseprecv(ilev)%d).
!        baseprecv(ilev)%av(mld_ap_nd_)   -  The entries of the local part of A(ilev)
!                                            outside the diagonal block, for block-Jacobi
!                                            sweeps.
!        baseprecv(ilev)%av(mld_ac_)      -  The local part of the matrix A(ilev).
!        baseprecv(ilev)%av(mld_sm_pr_)   -  The smoother prolongator.   
!                                            It maps vectors (ilev) ---> (ilev-1).
!        baseprecv(ilev)%av(mld_sm_pr_t_) -  The smoother prolongator transpose.   
!                                            It maps vectors (ilev-1) ---> (ilev).
!      baseprecv(ilev)%d         -  complex(kind(1.d0)), dimension(:), allocatable.
!                                              The diagonal entries of the U factor in the ILU
!                                   factorization of A(ilev).
!      baseprecv(ilev)%desc_data -  type(psb_desc_type).
!                                   The communication descriptor associated to the base
!                                   preconditioner,      i.e. to the sparse matrices needed
!                                   to apply the base preconditioner at the current level.
!      baseprecv(ilev)%desc_ac   -  type(psb_desc_type).
!                                              The communication descriptor associated to the sparse
!                                   matrix A(ilev), stored in baseprecv(ilev)%av(mld_ac_).
!      baseprecv(ilev)%iprcparm  -  integer, dimension(:), allocatable.
!                                              The integer parameters defining the base preconditioner
!                                   K(ilev).
!      baseprecv(ilev)%dprcparm  -  complex(kind(1.d0)), dimension(:), allocatable.
!                                   The real parameters defining the base preconditioner
!                                   K(ilev).
!      baseprecv(ilev)%perm      -  integer, dimension(:), allocatable.
!                                   The row and column permutations applied to the local
!                                   part of      A(ilev) (defined only if baseprecv(ilev)%
!                                   iprcparm(mld_sub_ren_)>0). 
!      baseprecv(ilev)%invperm   -  integer, dimension(:), allocatable.
!                                   The inverse of the permutation stored in
!                                   baseprecv(ilev)%perm.
!      baseprecv(ilev)%mlia      -  integer, dimension(:), allocatable.
!                                   The aggregation map (ilev-1) --> (ilev).
!                                   In case of non-smoothed aggregation, it is used
!                                   instead of mld_sm_pr_.
!      baseprecv(ilev)%nlaggr    -  integer, dimension(:), allocatable.
!                                   The number of aggregates (rows of A(ilev)) on the
!                                   various processes. 
!      baseprecv(ilev)%base_a    -  type(psb_zspmat_type), pointer.
!                                   Pointer (really a pointer!) to the base matrix of
!                                   the current level, i.e. the local part of A(ilev);
!                                   so we have a unified treatment of residuals. We
!                                   need this to avoid passing explicitly the matrix
!                                   A(ilev) to the routine which applies the
!                                   preconditioner.
!      baseprecv(ilev)%base_desc -  type(psb_desc_type), pointer.
!                                              Pointer to the communication descriptor associated
!                                   to the sparse matrix pointed by base_a.  
!      baseprecv(ilev)%dorig     -  complex(kind(1.d0)), dimension(:), allocatable.
!                                              Diagonal entries of the matrix pointed by base_a.
!                  
!   x          -  complex(kind(0.d0)), dimension(:), input.
!                 The local part of the vector X.
!   beta       -  complex(kind(0.d0)), input.
!                 The scalar beta.
!   y          -  complex(kind(0.d0)), dimension(:), input/output.
!                 The local part of the vector Y.
!   desc_data  -  type(psb_desc_type), input.
!                 The communication descriptor associated to the matrix to be
!                 preconditioned.
!   trans      -  character, optional.
!                 If trans='N','n' then op(M^(-1)) = M^(-1);
!                 if trans='T','t' then op(M^(-1)) = M^(-T) (transpose of M^(-1)).
!   work       -  complex(kind(0.d0)), dimension (:), optional, target.
!                 Workspace. Its size must be at least 4*psb_cd_get_local_cols(desc_data).
!   info       -  integer, output.
!                 Error code.
!
!   Note that when the LU factorization of the matrix A(ilev) is computed instead of
!   the ILU one, by using UMFPACK or SuperLU_dist, the corresponding L and U factors
!   are stored in data structures provided by UMFPACK or SuperLU_dist and pointed by
!   baseprecv(ilev)%iprcparm(mld_umf_ptr) or baseprecv(ilev)%iprcparm(mld_slu_ptr),
!   respectively.
!  
subroutine mld_zmlprec_aply(alpha,baseprecv,x,beta,y,desc_data,trans,work,info)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_zmlprec_aply

  implicit none

  ! Arguments
  type(psb_desc_type),intent(in)      :: desc_data
  type(mld_zbaseprc_type), intent(in) :: baseprecv(:)
  complex(kind(1.d0)),intent(in)      :: alpha,beta
  complex(kind(1.d0)),intent(in)      :: x(:)
  complex(kind(1.d0)),intent(inout)   :: y(:)
  character                           :: trans
  complex(kind(1.d0)),target          :: work(:)
  integer, intent(out)                :: info

  ! Local variables
  integer :: n_row,n_col
  integer :: ictxt,np,me,i,  nr2l,nc2l,err_act
  logical, parameter :: debug=.false., debugprt=.false.
  integer            :: ismth, nlev, ilev, icm
  character(len=20)  :: name

  type psb_mlprec_wrk_type
    complex(kind(1.d0)), allocatable  :: tx(:),ty(:),x2l(:),y2l(:)
  end type psb_mlprec_wrk_type
  type(psb_mlprec_wrk_type), allocatable :: mlprec_wrk(:)

  name='mld_zmlprec_aply'
  info = 0
  call psb_erractionsave(err_act)

  ictxt = psb_cd_get_context(desc_data)
  call psb_info(ictxt, me, np)

  if (debug) write(0,*) me,'Entry to mlprec_aply ',&
       & size(baseprecv)

  nlev = size(baseprecv)
  allocate(mlprec_wrk(nlev),stat=info) 
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Allocate')
    goto 9999      
  end if


  select case(baseprecv(2)%iprcparm(mld_ml_type_)) 

  case(mld_no_ml_)
    !
    ! No preconditioning, should not really get here
    ! 
    call psb_errpush(4010,name,a_err='mld_no_ml_ in mlprc_aply?')
    goto 9999      


  case(mld_add_ml_)

    !
    !       Additive multilevel
    !
    !   1. ! Apply the base preconditioner at level 1.
    !      ! The sum over the subdomains is carried      out in the
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


    !
    ! STEP 1
    !
    ! Apply the base preconditioner      at the finest level
    !
    call mld_baseprec_aply(alpha,baseprecv(1),x,beta,y,&
         & baseprecv(1)%base_desc,trans,work,info)
    if(info /=0) goto 9999

    allocate(mlprec_wrk(1)%x2l(size(x)),mlprec_wrk(1)%y2l(size(y)), stat=info)
    if (info /= 0) then 
      info=4025
      call psb_errpush(info,name,i_err=(/size(x)+size(y),0,0,0,0/),&
           & a_err='real(kind(1.d0))')
      goto 9999      
    end if

    mlprec_wrk(1)%x2l(:) = x(:) 


    !
    ! STEP 2
    !
    !
    !      For each level except the finest one ...
    !
    do ilev = 2, nlev
      n_row = psb_cd_get_local_rows(baseprecv(ilev-1)%base_desc)
      n_col = psb_cd_get_local_cols(baseprecv(ilev-1)%desc_data)
      nc2l  = psb_cd_get_local_cols(baseprecv(ilev)%desc_data)
      nr2l  = psb_cd_get_local_rows(baseprecv(ilev)%desc_data)
      allocate(mlprec_wrk(ilev)%x2l(nc2l),mlprec_wrk(ilev)%y2l(nc2l),&
           & mlprec_wrk(ilev)%tx(max(n_row,n_col)),&
           & mlprec_wrk(ilev)%ty(max(n_row,n_col)), stat=info)
      if (info /= 0) then 
        info=4025
        call psb_errpush(info,name,i_err=(/2*(nc2l+max(n_row,n_col)),0,0,0,0/),&
             & a_err='real(kind(1.d0))')
        goto 9999      
      end if

      mlprec_wrk(ilev)%x2l(:) = zzero
      mlprec_wrk(ilev)%y2l(:) = zzero
      mlprec_wrk(ilev)%tx(1:n_row) = mlprec_wrk(ilev-1)%x2l(1:n_row) 
      mlprec_wrk(ilev)%tx(n_row+1:max(n_row,n_col)) = zzero
      mlprec_wrk(ilev)%ty(:) = zzero


      ismth = baseprecv(ilev)%iprcparm(mld_aggr_kind_)
      icm   = baseprecv(ilev)%iprcparm(mld_coarse_mat_)
      if (ismth  /= mld_no_smooth_) then 
        !
        ! Apply the smoothed prolongator transpose
        !
        call psb_halo(mlprec_wrk(ilev-1)%x2l,baseprecv(ilev-1)%base_desc,&
             &  info,work=work) 
        if(info /=0) goto 9999

        call psb_csmm(zone,baseprecv(ilev)%av(mld_sm_pr_t_),mlprec_wrk(ilev-1)%x2l,&
             & zzero,mlprec_wrk(ilev)%x2l,info)
        if(info /=0) goto 9999

      else
        !
        ! Apply the raw aggregation map transpose (take a shortcut)
        !
        do i=1,n_row
          mlprec_wrk(ilev)%x2l(baseprecv(ilev)%mlia(i)) = &
               &  mlprec_wrk(ilev)%x2l(baseprecv(ilev)%mlia(i)) + &
               &  mlprec_wrk(ilev-1)%x2l(i)
        end do

      end if

      if (icm == mld_repl_mat_) then 
        call psb_sum(ictxt,mlprec_wrk(ilev)%x2l(1:nr2l))
      else if (icm /= mld_distr_mat_) then 
        write(0,*) 'Unknown value for baseprecv(2)%iprcparm(mld_coarse_mat_) ',icm 
      endif

      !
      ! Apply the base preconditioner
      !
      call mld_baseprec_aply(zone,baseprecv(ilev),&
           & mlprec_wrk(ilev)%x2l,zzero,mlprec_wrk(ilev)%y2l,&
           & baseprecv(ilev)%desc_data, 'N',work,info)

    enddo

    !
    ! STEP 3
    !
    !
    !      For each level except the finest one ...
    !
    do ilev =nlev,2,-1

      n_row = psb_cd_get_local_rows(baseprecv(ilev-1)%base_desc)
      n_col = psb_cd_get_local_cols(baseprecv(ilev-1)%desc_data)
      nc2l  = psb_cd_get_local_cols(baseprecv(ilev)%desc_data)
      nr2l  = psb_cd_get_local_rows(baseprecv(ilev)%desc_data)
      ismth = baseprecv(ilev)%iprcparm(mld_aggr_kind_)
      icm   = baseprecv(ilev)%iprcparm(mld_coarse_mat_)

      if (ismth  /= mld_no_smooth_) then 
        !
        ! Apply the smoothed prolongator
        !
        call psb_csmm(zone,baseprecv(ilev)%av(mld_sm_pr_),mlprec_wrk(ilev)%y2l,&
             & zone,mlprec_wrk(ilev-1)%y2l,info)
        if(info /=0) goto 9999

      else
        !
        ! Apply the raw aggregation map (take a shortcut)
        !
        do i=1, n_row
          mlprec_wrk(ilev-1)%y2l(i) = mlprec_wrk(ilev-1)%y2l(i) + &
               &   mlprec_wrk(ilev)%y2l(baseprecv(ilev)%mlia(i))
        enddo

      end if
    end do

    !
    ! STEP 4
    !
    ! Compute the output vector Y
    !
    call psb_geaxpby(alpha,mlprec_wrk(1)%y2l,zone,y,baseprecv(1)%base_desc,info)
    if(info /=0) goto 9999


  case(mld_mult_ml_)

    ! 
    !  Multiplicative multilevel (multiplicative among the levels, additive inside
    !  each level)
    !
    !  Pre/post-smoothing versions 
    !

    select case(baseprecv(2)%iprcparm(mld_smooth_pos_))

    case(mld_post_smooth_)

      !
      !    Post-smoothing
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
      !                  ! Transfer Y(ilev+1) to the next finer level.
      !           Y(ilev) = AV(ilev+1; sm_pr_)*Y(ilev+1)
      !
      !         ! Compute the residual at the current level and apply to it      the
      !         ! base preconditioner. The sum over the subdomains is carried out
      !         ! in the application of K(ilev).
      !           Y(ilev) = Y(ilev) + (K(ilev)^(-1))*(X(ilev)-A(ilev)*Y(ilev))
      !
      !        ENDDO
      !
      !    5.  Yext    = beta*Yext + alpha*Y(1)
      ! 

      !
      ! STEP 1
      !
      ! Copy the input vector X
      !
      if (debug) write(0,*) me, 'mlprec_aply desc_data',&
           & allocated(desc_data%matrix_data)    

      n_col = psb_cd_get_local_cols(desc_data)
      nc2l  = psb_cd_get_local_cols(baseprecv(1)%desc_data)

      allocate(mlprec_wrk(1)%x2l(nc2l),mlprec_wrk(1)%y2l(nc2l), &
           & mlprec_wrk(1)%tx(nc2l), stat=info)
      mlprec_wrk(1)%x2l(:) = zzero
      mlprec_wrk(1)%y2l(:) = zzero
      mlprec_wrk(1)%tx(:)  = zzero

      call psb_geaxpby(zone,x,zzero,mlprec_wrk(1)%tx,&
           & baseprecv(1)%base_desc,info)
      call psb_geaxpby(zone,x,zzero,mlprec_wrk(1)%x2l,&
           & baseprecv(1)%base_desc,info)

      !
      ! STEP 2
      !
      !      For each level but the finest one ...
      !
      do ilev=2, nlev

        n_row = psb_cd_get_local_rows(baseprecv(ilev-1)%base_desc)
        n_col = psb_cd_get_local_cols(baseprecv(ilev-1)%desc_data)
        nc2l  = psb_cd_get_local_cols(baseprecv(ilev)%desc_data)
        nr2l  = psb_cd_get_local_rows(baseprecv(ilev)%desc_data)
        ismth = baseprecv(ilev)%iprcparm(mld_aggr_kind_)
        icm   = baseprecv(ilev)%iprcparm(mld_coarse_mat_)

        if (debug) write(0,*) me, 'mlprec_aply starting up sweep ',&
             & ilev,allocated(baseprecv(ilev)%iprcparm),n_row,n_col,&
             & nc2l, nr2l,ismth

        allocate(mlprec_wrk(ilev)%tx(nc2l),mlprec_wrk(ilev)%y2l(nc2l),&
             &   mlprec_wrk(ilev)%x2l(nc2l), stat=info)

        if (info /= 0) then 
          info=4025
          call psb_errpush(info,name,i_err=(/4*nc2l,0,0,0,0/),&
               & a_err='real(kind(1.d0))')
          goto 9999      
        end if

        mlprec_wrk(ilev)%x2l(:) = zzero
        mlprec_wrk(ilev)%y2l(:) = zzero
        mlprec_wrk(ilev)%tx(:) = zzero
        if (ismth  /= mld_no_smooth_) then 
          !
          ! Apply the smoothed prolongator transpose
          !
          if (debug) write(0,*) me, 'mlprec_aply halo in up sweep ', ilev
          call psb_halo(mlprec_wrk(ilev-1)%x2l,&
               &  baseprecv(ilev-1)%base_desc,info,work=work) 
          if(info /=0) goto 9999
          if (debug) write(0,*) me, 'mlprec_aply csmm in up sweep ', ilev
          call psb_csmm(zone,baseprecv(ilev)%av(mld_sm_pr_t_),mlprec_wrk(ilev-1)%x2l, &
               & zzero,mlprec_wrk(ilev)%x2l,info)
          if(info /=0) goto 9999

        else
          !
          ! Apply the raw aggregation map transpose (take a shortcut)
          !
          do i=1,n_row
            mlprec_wrk(ilev)%x2l(baseprecv(ilev)%mlia(i)) = &
                 & mlprec_wrk(ilev)%x2l(baseprecv(ilev)%mlia(i)) + &
                 & mlprec_wrk(ilev-1)%x2l(i)
          end do

        end if

        if (debug) write(0,*) me, 'mlprec_aply possible sum in up sweep ', &
             & ilev,icm,associated(baseprecv(ilev)%base_desc),mld_repl_mat_
        if (debug) write(0,*) me, 'mlprec_aply geaxpby in up sweep X', &
             & ilev,associated(baseprecv(ilev)%base_desc),&
             & baseprecv(ilev)%base_desc%matrix_data(psb_n_row_),&
             & baseprecv(ilev)%base_desc%matrix_data(psb_n_col_),&
             & size(mlprec_wrk(ilev)%tx),size(mlprec_wrk(ilev)%x2l)

        if (icm == mld_repl_mat_) Then 
          if (debug) write(0,*) 'Entering psb_sum ',nr2l
          call psb_sum(ictxt,mlprec_wrk(ilev)%x2l(1:nr2l))
        else if (icm  /= mld_distr_mat_) Then 
          write(0,*) 'Unknown value for baseprecv(2)%iprcparm(mld_coarse_mat_) ', icm 
        endif

        !
        ! update x2l
        !
        call psb_geaxpby(zone,mlprec_wrk(ilev)%x2l,zzero,mlprec_wrk(ilev)%tx,&
             & baseprecv(ilev)%base_desc,info)
        if(info /=0) goto 9999
        if (debug) write(0,*) me, 'mlprec_aply done up sweep ',&
             & ilev

      enddo

      !
      ! STEP 3
      !
      ! Apply the base preconditioner at the coarsest level
      !
      call mld_baseprec_aply(zone,baseprecv(nlev),mlprec_wrk(nlev)%x2l, &
           & zzero, mlprec_wrk(nlev)%y2l,baseprecv(nlev)%desc_data,'N',work,info)

      if(info /=0) goto 9999
      if (debug) write(0,*) me, 'mlprec_aply done prc_apl ',&
           & nlev

      !
      ! STEP 4
      !
      !      For each level but the coarsest one      ...
      !
      do ilev=nlev-1, 1, -1

        if (debug) write(0,*) me, 'mlprec_aply starting down sweep',ilev
        ismth = baseprecv(ilev+1)%iprcparm(mld_aggr_kind_)
        n_row = psb_cd_get_local_rows(baseprecv(ilev)%base_desc)

        if (ismth  /= mld_no_smooth_) then
          !
          ! Apply the smoothed prolongator
          !  
          if (ismth == mld_smooth_prol_) &
               & call psb_halo(mlprec_wrk(ilev+1)%y2l,baseprecv(ilev+1)%desc_data,&
               &  info,work=work) 
          call psb_csmm(zone,baseprecv(ilev+1)%av(mld_sm_pr_),mlprec_wrk(ilev+1)%y2l,&
               &  zzero,mlprec_wrk(ilev)%y2l,info)
          if(info /=0) goto 9999

        else
          !
          ! Apply the raw aggregation map (take a shortcut)
          !
          mlprec_wrk(ilev)%y2l(:) = zzero
          do i=1, n_row
            mlprec_wrk(ilev)%y2l(i) = mlprec_wrk(ilev)%y2l(i) + &
                 & mlprec_wrk(ilev+1)%y2l(baseprecv(ilev+1)%mlia(i))
          enddo

        end if

        !
        ! Compute the residual
        !
        call psb_spmm(-zone,baseprecv(ilev)%base_a,mlprec_wrk(ilev)%y2l,&
             &   zone,mlprec_wrk(ilev)%tx,baseprecv(ilev)%base_desc,info,work=work)

        if(info /=0) goto 9999

        !
        ! Apply the base preconditioner
        !
        call mld_baseprec_aply(zone,baseprecv(ilev),mlprec_wrk(ilev)%tx,&
             & zone,mlprec_wrk(ilev)%y2l,baseprecv(ilev)%base_desc, trans, work,info)

        if(info /=0) goto 9999
        if (debug) write(0,*) me, 'mlprec_aply done down sweep',ilev
      enddo

      !
      ! STEP 5
      !
      ! Compute the output vector Y
      !
      call psb_geaxpby(alpha,mlprec_wrk(1)%y2l,beta,y,baseprecv(1)%base_desc,info)

      if(info /=0) goto 9999


    case(mld_pre_smooth_)

      !
      !    Pre-smoothing
      !
      !    1.   X(1)  = Xext
      !
      !    2. ! Apply the base preconditioner at the finest level.
      !         Y(1)  = (K(1)^(-1))*X(1)
      !
      !    3. ! Compute the residual at the finest level.
      !         TX(1) = X(1) - A(1)*Y(1)
      !
      !    4.   DO ilev=2, nlev
      !
      !          ! Transfer      the residual to the current (coarser) level.
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
      ! STEP 1
      !
      ! Copy the input vector X
      !
      n_col = psb_cd_get_local_cols(desc_data)
      nc2l  = psb_cd_get_local_cols(baseprecv(1)%desc_data)

      allocate(mlprec_wrk(1)%x2l(nc2l),mlprec_wrk(1)%y2l(nc2l), &
           & mlprec_wrk(1)%tx(nc2l), stat=info)
      if (info /= 0) then 
        info=4025
        call psb_errpush(info,name,i_err=(/4*nc2l,0,0,0,0/),&
             & a_err='real(kind(1.d0))')
        goto 9999
      end if

      mlprec_wrk(1)%y2l(:) = zzero
      mlprec_wrk(1)%x2l(:) = x

      !
      ! STEP 2
      !
      ! Apply the base preconditioner at the finest level
      !
      call mld_baseprec_aply(zone,baseprecv(1),mlprec_wrk(1)%x2l,&
           &  zzero,mlprec_wrk(1)%y2l,&
           &  baseprecv(1)%base_desc,&
           &  trans,work,info)

      if(info /=0) goto 9999


      !
      ! STEP 3
      !
      ! Compute the residual at the finest level
      !
      mlprec_wrk(1)%tx = mlprec_wrk(1)%x2l

      call psb_spmm(-zone,baseprecv(1)%base_a,mlprec_wrk(1)%y2l,&
           & zone,mlprec_wrk(1)%tx,baseprecv(1)%base_desc,info,work=work)
      if(info /=0) goto 9999

      !
      ! STEP 4
      !
      ! For each level but the finest one ...
      !
      do ilev = 2, nlev

        n_row = psb_cd_get_local_rows(baseprecv(ilev-1)%base_desc)
        n_col = psb_cd_get_local_cols(baseprecv(ilev-1)%desc_data)
        nc2l  = psb_cd_get_local_cols(baseprecv(ilev)%desc_data)
        nr2l  = psb_cd_get_local_rows(baseprecv(ilev)%desc_data)
        ismth = baseprecv(ilev)%iprcparm(mld_aggr_kind_)
        icm   = baseprecv(ilev)%iprcparm(mld_coarse_mat_)

        allocate(mlprec_wrk(ilev)%tx(nc2l),mlprec_wrk(ilev)%y2l(nc2l),&
             &   mlprec_wrk(ilev)%x2l(nc2l), stat=info)
        if (info /= 0) then 
          info=4025
          call psb_errpush(info,name,i_err=(/4*nc2l,0,0,0,0/),&
               & a_err='real(kind(1.d0))')
          goto 9999      
        end if

        mlprec_wrk(ilev)%x2l(:) = zzero
        mlprec_wrk(ilev)%y2l(:) = zzero
        mlprec_wrk(ilev)%tx(:) = zzero


        if (ismth  /= mld_no_smooth_) then 
          !
          ! Apply the smoothed prolongator transpose
          !
          call psb_halo(mlprec_wrk(ilev-1)%tx,baseprecv(ilev-1)%base_desc,&
               & info,work=work) 
          if(info /=0) goto 9999

          call psb_csmm(zone,baseprecv(ilev)%av(mld_sm_pr_t_),mlprec_wrk(ilev-1)%tx,zzero,&
               & mlprec_wrk(ilev)%x2l,info)
          if(info /=0) goto 9999

        else
          !
          ! Apply the raw aggregation map transpose (take a shortcut)
          !
          mlprec_wrk(ilev)%x2l = zzero
          do i=1,n_row
            mlprec_wrk(ilev)%x2l(baseprecv(ilev)%mlia(i)) = &
                 & mlprec_wrk(ilev)%x2l(baseprecv(ilev)%mlia(i)) + &
                 &  mlprec_wrk(ilev-1)%tx(i)
          end do
        end if

        if (icm ==mld_repl_mat_) then 
          call psb_sum(ictxt,mlprec_wrk(ilev)%x2l(1:nr2l))
        else if (icm  /= mld_distr_mat_) then 
          write(0,*) 'Unknown value for baseprecv(2)%iprcparm(mld_coarse_mat_) ', icm 
        endif

        !
        ! Apply the base preconditioner
        !
        call mld_baseprec_aply(zone,baseprecv(ilev),mlprec_wrk(ilev)%x2l,&
             & zzero,mlprec_wrk(ilev)%y2l,baseprecv(ilev)%desc_data, 'N',work,info)

        if(info /=0) goto 9999

        !
        ! Compute the residual (at all levels but the coarsest one)
        !
        if (ilev < nlev) then
          mlprec_wrk(ilev)%tx = mlprec_wrk(ilev)%x2l
          call psb_spmm(-zone,baseprecv(ilev)%base_a,mlprec_wrk(ilev)%y2l,&
               & zone,mlprec_wrk(ilev)%tx,baseprecv(ilev)%base_desc,info,work=work)
          if(info /=0) goto 9999
        endif

      enddo

      !
      ! STEP 5
      !
      !      For each level but the coarsest one ...
      !
      do ilev = nlev-1, 1, -1

        ismth = baseprecv(ilev+1)%iprcparm(mld_aggr_kind_)
        n_row = psb_cd_get_local_rows(baseprecv(ilev)%base_desc)

        if (ismth  /= mld_no_smooth_) then 
          !
          ! Apply the smoothed prolongator
          !
          if (ismth == mld_smooth_prol_) &
               & call psb_halo(mlprec_wrk(ilev+1)%y2l,&
               & baseprecv(ilev+1)%desc_data,info,work=work) 
          call psb_csmm(zone,baseprecv(ilev+1)%av(mld_sm_pr_),mlprec_wrk(ilev+1)%y2l,&
               & zone,mlprec_wrk(ilev)%y2l,info)

          if(info /=0) goto 9999

        else
          !
          ! Apply the raw aggregation map (take a shortcut)
          !
          do i=1, n_row
            mlprec_wrk(ilev)%y2l(i) = mlprec_wrk(ilev)%y2l(i) + &
                 & mlprec_wrk(ilev+1)%y2l(baseprecv(ilev+1)%mlia(i))
          enddo

        end if

      enddo

      !
      ! STEP 6
      !
      !      Compute the output vector Y
      !
      call psb_geaxpby(alpha,mlprec_wrk(1)%y2l,beta,y,&
           &  baseprecv(1)%base_desc,info)

      if(info /=0) goto 9999


    case(mld_twoside_smooth_)

      !
      !    Pre- and post-smoothing (symmetrized)
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
      !          ! Transfer      the residual to the current (coarser) level
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

      !
      ! STEP 1
      !
      ! Copy the input vector X
      !
      n_col = psb_cd_get_local_cols(desc_data)
      nc2l  = psb_cd_get_local_cols(baseprecv(1)%desc_data)

      allocate(mlprec_wrk(1)%x2l(nc2l),mlprec_wrk(1)%y2l(nc2l), &
           & mlprec_wrk(1)%ty(nc2l), mlprec_wrk(1)%tx(nc2l), stat=info)

      if (info /= 0) then 
        info=4025
        call psb_errpush(info,name,i_err=(/4*nc2l,0,0,0,0/),&
             & a_err='real(kind(1.d0))')
        goto 9999
      end if
      mlprec_wrk(1)%x2l(:) = zzero
      mlprec_wrk(1)%y2l(:) = zzero
      mlprec_wrk(1)%tx(:)  = zzero
      mlprec_wrk(1)%ty(:)  = zzero
      

      call psb_geaxpby(zone,x,zzero,mlprec_wrk(1)%x2l,&
           & baseprecv(1)%base_desc,info)
      call psb_geaxpby(zone,x,zzero,mlprec_wrk(1)%tx,&
           & baseprecv(1)%base_desc,info)

      !
      ! STEP 2
      !
      ! Apply the base preconditioner at the finest level
      !
      call mld_baseprec_aply(zone,baseprecv(1),mlprec_wrk(1)%x2l,&
           &  zzero,mlprec_wrk(1)%y2l,&
           &  baseprecv(1)%base_desc,&
           &  trans,work,info)

      if(info /=0) goto 9999

      !
      ! STEP 3
      !
      ! Compute the residual at the finest level
      !
      mlprec_wrk(1)%ty = mlprec_wrk(1)%x2l
      call psb_spmm(-zone,baseprecv(1)%base_a,mlprec_wrk(1)%y2l,&
           & zone,mlprec_wrk(1)%ty,baseprecv(1)%base_desc,info,work=work)
      if(info /=0) goto 9999

      !
      ! STEP 4
      !
      ! For each level but the finest one ...
      !
      do ilev = 2, nlev

        n_row = psb_cd_get_local_rows(baseprecv(ilev-1)%base_desc)
        n_col = psb_cd_get_local_cols(baseprecv(ilev-1)%desc_data)
        nc2l  = psb_cd_get_local_cols(baseprecv(ilev)%desc_data)
        nr2l  = psb_cd_get_local_rows(baseprecv(ilev)%desc_data)
        ismth = baseprecv(ilev)%iprcparm(mld_aggr_kind_)
        icm   = baseprecv(ilev)%iprcparm(mld_coarse_mat_)
        allocate(mlprec_wrk(ilev)%ty(nc2l),mlprec_wrk(ilev)%y2l(nc2l),&
             &   mlprec_wrk(ilev)%x2l(nc2l), stat=info)

        if (info /= 0) then 
          info=4025
          call psb_errpush(info,name,i_err=(/4*nc2l,0,0,0,0/),&
               & a_err='real(kind(1.d0))')
          goto 9999      
        end if

        mlprec_wrk(ilev)%x2l(:) = zzero
        mlprec_wrk(ilev)%y2l(:) = zzero
        mlprec_wrk(ilev)%tx(:)  = zzero
        mlprec_wrk(ilev)%ty(:)  = zzero
      

        if (ismth  /= mld_no_smooth_) then 
          !
          ! Apply the smoothed prolongator transpose
          !
          call psb_halo(mlprec_wrk(ilev-1)%ty,baseprecv(ilev-1)%base_desc,&
               & info,work=work) 
          if(info /=0) goto 9999
          call psb_csmm(zone,baseprecv(ilev)%av(mld_sm_pr_t_),mlprec_wrk(ilev-1)%ty,zzero,&
               & mlprec_wrk(ilev)%x2l,info)
          if(info /=0) goto 9999

        else
          !
          ! Apply the raw aggregation map transpose (take a shortcut)
          !
          mlprec_wrk(ilev)%x2l = zzero
          do i=1,n_row
            mlprec_wrk(ilev)%x2l(baseprecv(ilev)%mlia(i)) = &
                 & mlprec_wrk(ilev)%x2l(baseprecv(ilev)%mlia(i)) + &
                 &  mlprec_wrk(ilev-1)%ty(i)
          end do
        end if

        if (icm == mld_repl_mat_) then 
          call psb_sum(ictxt,mlprec_wrk(ilev)%x2l(1:nr2l))
        else if (icm /= mld_distr_mat_) then 
          write(0,*) 'Unknown value for baseprecv(2)%iprcparm(mld_coarse_mat_) ', icm 
        endif

        call psb_geaxpby(zone,mlprec_wrk(ilev)%x2l,zzero,mlprec_wrk(ilev)%tx,&
             & baseprecv(ilev)%base_desc,info)
        if(info /=0) goto 9999

        !
        ! Apply the base preconditioner
        !
        call mld_baseprec_aply(zone,baseprecv(ilev),mlprec_wrk(ilev)%x2l,&
             & zzero,mlprec_wrk(ilev)%y2l,baseprecv(ilev)%desc_data, 'N',work,info)

        if(info /=0) goto 9999

        !
        ! Compute the residual (at all levels but the coarsest one)
        !
        if(ilev < nlev) then
          mlprec_wrk(ilev)%ty = mlprec_wrk(ilev)%x2l
          call psb_spmm(-zone,baseprecv(ilev)%base_a,mlprec_wrk(ilev)%y2l,&
               & zone,mlprec_wrk(ilev)%ty,baseprecv(ilev)%base_desc,info,work=work)
          if(info /=0) goto 9999
        endif

      enddo

      !
      ! STEP 5
      !
      ! For each level but the coarsest one ...
      !
      do ilev=nlev-1, 1, -1

        ismth = baseprecv(ilev+1)%iprcparm(mld_aggr_kind_)
        n_row = psb_cd_get_local_rows(baseprecv(ilev)%base_desc)

        if (ismth  /= mld_no_smooth_) then
          !
          ! Apply the smoothed prolongator
          ! 
          if (ismth == mld_smooth_prol_) &
               & call psb_halo(mlprec_wrk(ilev+1)%y2l,baseprecv(ilev+1)%desc_data,&
               &  info,work=work) 
          call psb_csmm(zone,baseprecv(ilev+1)%av(mld_sm_pr_),mlprec_wrk(ilev+1)%y2l,&
               &  zone,mlprec_wrk(ilev)%y2l,info)
          if(info /=0) goto 9999

        else
          !
          ! Apply the raw aggregation map (take a shortcut)
          !
          do i=1, n_row
            mlprec_wrk(ilev)%y2l(i) = mlprec_wrk(ilev)%y2l(i) + &
                 & mlprec_wrk(ilev+1)%y2l(baseprecv(ilev+1)%mlia(i))
          enddo

        end if

        !
        ! Compute the residual
        !
        call psb_spmm(-zone,baseprecv(ilev)%base_a,mlprec_wrk(ilev)%y2l,&
             &   zone,mlprec_wrk(ilev)%tx,baseprecv(ilev)%base_desc,info,work=work)

        if(info /=0) goto 9999

        !
        ! Apply the base preconditioner
        !
        call mld_baseprec_aply(zone,baseprecv(ilev),mlprec_wrk(ilev)%tx,&
             & zone,mlprec_wrk(ilev)%y2l,baseprecv(ilev)%base_desc, trans, work,info)

        if(info /=0) goto 9999

      enddo

      !
      ! STEP 6
      !
      ! Compute the output vector Y
      !
      call psb_geaxpby(alpha,mlprec_wrk(1)%y2l,beta,y,&
           &   baseprecv(1)%base_desc,info)

      if(info /=0) goto 9999

    case default

      call psb_errpush(4013,name,a_err='wrong smooth_pos',&
           &  i_Err=(/baseprecv(2)%iprcparm(mld_smooth_pos_),0,0,0,0/))
      goto 9999      

    end select

  case default
    call psb_errpush(4013,name,a_err='wrong mltype',&
         &  i_Err=(/baseprecv(2)%iprcparm(mld_ml_type_),0,0,0,0/))
    goto 9999      

  end select

  deallocate(mlprec_wrk)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name)
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_zmlprec_aply

