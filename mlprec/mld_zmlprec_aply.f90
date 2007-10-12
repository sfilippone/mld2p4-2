!!$ 
!!$ 
!!$                    MD2P4
!!$    Multilevel Domain Decomposition Parallel Preconditioner Package for PSBLAS
!!$                      for 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$                       Daniela di Serafino    Second University of Naples
!!$                       Pasqua D'Ambra         ICAR-CNR                      
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the MD2P4 group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MD2P4 GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
subroutine mld_zmlprec_aply(alpha,baseprecv,x,beta,y,desc_data,trans,work,info)
  !
  !  Compute   Y <-  beta*Y + alpha*K^-1 X 
  !  where K is a multilevel  preconditioner stored in baseprecv
  ! 
  !  cfr.: Smith, Biorstad & Gropp
  !        Domain Decomposition
  !        Cambridge Univ. Press 
  ! 
  !  To each level I there corresponds a matrix A(I) and a preconditioner K(I)
  !
  !  A notational difference: in the DD reference above the preconditioner for 
  !  a given level K(I) is written out as a sum over the subdomains
  !
  !  SUM_k(R_k^T A_k R_k) 
  !
  !  whereas in this code the sum is implicit in the parallelization, 
  !  i.e. each process takes care of one subdomain, and for each level we have 
  !  as many subdomains as there are processes (except for the coarsest level where 
  !  we might have a replicated index space). Thus the sum apparently disappears 
  !  from our code, but only apparently, because it is implicit in the call 
  !  to mld_baseprec_aply. 
  !
  !  A bit of description of the baseprecv(:) data structure:
  !   1. Number of levels = NLEV = size(baseprecv(:))
  !   2. baseprecv(ilev)%av(:)    sparse matrices needed for the current level. 
  !      Includes:
  !   2.1.:  baseprecv(ilev)%av(mld_l_pr_)    L factor of ILU preconditioners
  !   2.2.:  baseprecv(ilev)%av(mld_u_pr_)    U factor of ILU preconditioners
  !   2.3.:  baseprecv(ilev)%av(mld_ap_nd_)   Off-diagonal part of A for Jacobi sweeps
  !   2.4.:  baseprecv(ilev)%av(mld_ac_)      Aggregated matrix of level ILEV 
  !   2.5.:  baseprecv(ilev)%av(mld_sm_pr_t_) Smoother prolongator transpose; maps vectors  
  !                                          (ilev-1) --->  (ilev) 
  !   2.6.:  baseprecv(ilev)%av(mld_sm_pr_)   Smoother prolongator; maps vectors  
  !                                          (ilev)   --->  (ilev-1) 
  !   Shouldn't we keep just one of them and handle transpose in the sparse BLAS? maybe 
  !
  !   3.    baseprecv(ilev)%desc_data     comm descriptor for level ILEV
  !   4.    baseprecv(ilev)%base_a        Pointer (really a pointer!) to the base matrix 
  !         baseprecv(ilev)%base_desc     of the current level, i.e.: if ILEV=1 then  A
  !                                       else the aggregated matrix av(mld_ac_); so we have 
  !                                       a unified treatment of residuals. Need this to 
  !                                       avoid passing explicitly matrix A to the 
  !                                       outer prec. routine
  !   5.    baseprecv(ilev)%mlia          The aggregation map from (ilev-1)-->(ilev)
  !                                       if no smoother, it is used instead of sm_pr_
  !   6.    baseprecv(ilev)%nlaggr        Number of aggregates on the various procs. 
  !   

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_zmlprec_aply

  implicit none

  type(psb_desc_type),intent(in)      :: desc_data
  type(mld_zbaseprc_type), intent(in) :: baseprecv(:)
  complex(kind(1.d0)),intent(in)      :: alpha,beta
  complex(kind(1.d0)),intent(inout)   :: x(:), y(:)
  character                           :: trans
  complex(kind(1.d0)),target          :: work(:)
  integer, intent(out)                :: info


  ! Local variables
  integer :: n_row,n_col
  character     ::diagl, diagu
  integer :: ictxt,np,me,i, isz, nr2l,nc2l,err_act, iptype, int_err(5)
  real(kind(1.d0)) :: omega
  real(kind(1.d0)) :: t1, t2, t3, t4, t5, t6, t7
  logical, parameter          :: debug=.false., debugprt=.false.
  integer      :: ismth, nlev, ilev, icm
  character(len=20)   :: name, ch_err

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
    ! Should not really get here.
    call psb_errpush(4010,name,a_err='mld_no_ml_ in mlprc_aply?')
    goto 9999      


  case(mld_add_ml_)

    
    !
    !    Additive is very simple. 
    !    1.  X(1) = Xext
    !    2.  DO ILEV=2,NLEV 
    !           X(ILEV) = AV(PR_SM_T_)*X(ILEV-1)    
    !           Y(ILEV) = (K(ILEV)**(-1))*X(ILEV)
    !    3.  DO  ILEV=NLEV-1,1,-1
    !           Y(ILEV) = AV(PR_SM_)*Y(ILEV+1)     
    !    4.  Yext    = beta*Yext + alpha*Y(1)
    !
    !    Note: level numbering reversed wrt ref. DD, i.e. 
    !         1..NLEV <=>  (j) <-> 0
    

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
        ! Smoothed aggregation
        !
        call psb_halo(mlprec_wrk(ilev-1)%x2l,baseprecv(ilev-1)%base_desc,&
             &  info,work=work) 
        if(info /=0) goto 9999

        call psb_csmm(zone,baseprecv(ilev)%av(mld_sm_pr_t_),mlprec_wrk(ilev-1)%x2l,&
             & zzero,mlprec_wrk(ilev)%x2l,info)
        if(info /=0) goto 9999

      else
        !
        ! Raw  aggregation, may take shortcut
        !
        do i=1,n_row
          mlprec_wrk(ilev)%x2l(baseprecv(ilev)%mlia(i)) = &
               &  mlprec_wrk(ilev)%x2l(baseprecv(ilev)%mlia(i)) + &
               &  mlprec_wrk(ilev-1)%x2l(i)
        end do

      end if

      if (icm ==mld_repl_mat_) Then 
        call psb_sum(ictxt,mlprec_wrk(ilev)%x2l(1:nr2l))
      else if (icm/= mld_distr_mat_) Then 
        write(0,*) 'Unknown value for baseprecv(2)%iprcparm(mld_coarse_mat_) ',icm 
      endif

      call mld_baseprec_aply(zone,baseprecv(ilev),&
           & mlprec_wrk(ilev)%x2l,zzero,mlprec_wrk(ilev)%y2l,&
           & baseprecv(ilev)%desc_data, 'N',work,info)

    enddo

    do ilev =nlev,2,-1

      n_row = psb_cd_get_local_rows(baseprecv(ilev-1)%base_desc)
      n_col = psb_cd_get_local_cols(baseprecv(ilev-1)%desc_data)
      nc2l  = psb_cd_get_local_cols(baseprecv(ilev)%desc_data)
      nr2l  = psb_cd_get_local_rows(baseprecv(ilev)%desc_data)
      ismth = baseprecv(ilev)%iprcparm(mld_aggr_kind_)
      icm   = baseprecv(ilev)%iprcparm(mld_coarse_mat_)

      if (ismth  /= mld_no_smooth_) then 

        call psb_csmm(zone,baseprecv(ilev)%av(mld_sm_pr_),mlprec_wrk(ilev)%y2l,&
             & zone,mlprec_wrk(ilev-1)%y2l,info)
        if(info /=0) goto 9999

      else

        do i=1, n_row
          mlprec_wrk(ilev-1)%y2l(i) = mlprec_wrk(ilev-1)%y2l(i) + &
               &   mlprec_wrk(ilev)%y2l(baseprecv(ilev)%mlia(i))
        enddo

      end if
    end do

    call psb_geaxpby(alpha,mlprec_wrk(1)%y2l,zone,y,baseprecv(1)%base_desc,info)
    if(info /=0) goto 9999


  case(mld_mult_ml_)

    ! 
    !  Multiplicative multilevel
    !  Pre/post smoothing versions. 
    !

    select case(baseprecv(2)%iprcparm(mld_smooth_pos_))

    case(mld_post_smooth_)

      !
      !    Post smoothing. 
      !    1.   X(1) = Xext
      !    2.   DO ILEV=2, NLEV :: X(ILEV) = AV(PR_SM_T_,ILEV)*X(ILEV-1)    
      !    3.   Y(NLEV) = (K(NLEV)**(-1))*X(NLEV)
      !    4.   DO  ILEV=NLEV-1,1,-1
      !          Y(ILEV) = AV(PR_SM_,ILEV+1)*Y(ILEV+1)     
      !          Y(ILEV) = Y(ILEV) + (K(ILEV)**(-1))*(X(ILEV)-A(ILEV)*Y(ILEV))
      !
      !    5.  Yext    = beta*Yext + alpha*Y(1)
      !
      !    Note: level numbering reversed wrt ref. DD, i.e. 
      !         1..NLEV <=>  (j) <-> 0
      !
      !    Also: post smoothing in the ref. DD is only presented for NLEV=2. 
      ! 
      ! 

      if (debug) write(0,*) me, 'mlpr_aply desc_data',&
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

      do ilev=2, nlev

      
        n_row = psb_cd_get_local_rows(baseprecv(ilev-1)%base_desc)
        n_col = psb_cd_get_local_cols(baseprecv(ilev-1)%desc_data)
        nc2l  = psb_cd_get_local_cols(baseprecv(ilev)%desc_data)
        nr2l  = psb_cd_get_local_rows(baseprecv(ilev)%desc_data)
        ismth = baseprecv(ilev)%iprcparm(mld_aggr_kind_)
        icm   = baseprecv(ilev)%iprcparm(mld_coarse_mat_)
          
        if (debug) write(0,*) me, 'mlpr_aply starting up sweep ',&
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
          ! Smoothed aggregation
          !
          if (debug) write(0,*) me, 'mlpr_aply halo in up sweep ', ilev
          
          call psb_halo(mlprec_wrk(ilev-1)%x2l,&
               &  baseprecv(ilev-1)%base_desc,info,work=work) 
          if(info /=0) goto 9999

          call psb_csmm(zone,baseprecv(ilev)%av(mld_sm_pr_t_),mlprec_wrk(ilev-1)%x2l, &
               & zzero,mlprec_wrk(ilev)%x2l,info)
          if(info /=0) goto 9999

        else
          !
          ! Raw  aggregation, may take shortcut
          !
          do i=1,n_row
            mlprec_wrk(ilev)%x2l(baseprecv(ilev)%mlia(i)) = &
                 & mlprec_wrk(ilev)%x2l(baseprecv(ilev)%mlia(i)) + &
                 & mlprec_wrk(ilev-1)%x2l(i)
          end do
        end if

        if (debug) write(0,*) me, 'mlpr_aply possible sum in up sweep ', &
             & ilev,icm,associated(baseprecv(ilev)%base_desc),mld_repl_mat_
        if (debug) write(0,*) me, 'mlpr_aply geaxpby in up sweep X', &
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
        call psb_geaxpby(zone,mlprec_wrk(ilev)%x2l,zzero,mlprec_wrk(ilev)%tx,&
             & baseprecv(ilev)%base_desc,info)
        if(info /=0) goto 9999

      enddo


      call mld_baseprec_aply(zone,baseprecv(nlev),mlprec_wrk(nlev)%x2l, &
           & zzero, mlprec_wrk(nlev)%y2l,baseprecv(nlev)%desc_data,'N',work,info)

      if(info /=0) goto 9999


      do ilev=nlev-1, 1, -1
        ismth = baseprecv(ilev+1)%iprcparm(mld_aggr_kind_)
        n_row = psb_cd_get_local_rows(baseprecv(ilev)%base_desc)

        if (ismth  /= mld_no_smooth_) then  
          if (ismth == mld_smooth_prol_) &
               & call psb_halo(mlprec_wrk(ilev+1)%y2l,baseprecv(ilev+1)%desc_data,&
               &  info,work=work) 
          call psb_csmm(zone,baseprecv(ilev+1)%av(mld_sm_pr_),mlprec_wrk(ilev+1)%y2l,&
               &  zzero,mlprec_wrk(ilev)%y2l,info)
          if(info /=0) goto 9999

        else
          mlprec_wrk(ilev)%y2l(:) = zzero
          do i=1, n_row
            mlprec_wrk(ilev)%y2l(i) = mlprec_wrk(ilev)%y2l(i) + &
                 & mlprec_wrk(ilev+1)%y2l(baseprecv(ilev+1)%mlia(i))
          enddo

        end if

        call psb_spmm(-zone,baseprecv(ilev)%base_a,mlprec_wrk(ilev)%y2l,&
             &   zone,mlprec_wrk(ilev)%tx,baseprecv(ilev)%base_desc,info,work=work)

        if(info /=0) goto 9999

        call mld_baseprec_aply(zone,baseprecv(ilev),mlprec_wrk(ilev)%tx,&
             & zone,mlprec_wrk(ilev)%y2l,baseprecv(ilev)%base_desc, trans, work,info)

        if(info /=0) goto 9999

      enddo

      call psb_geaxpby(alpha,mlprec_wrk(1)%y2l,beta,y,baseprecv(1)%base_desc,info)

      if(info /=0) goto 9999


    case(mld_pre_smooth_)

      !
      !    Pre smoothing. 
      !    1.   X(1)  = Xext
      !    2.   Y(1)  = (K(1)**(-1))*X(1)
      !    3.   TX(1) = X(1) - A(1)*Y(1)
      !    4.   DO ILEV=2, NLEV 
      !          X(ILEV) = AV(PR_SM_T_,ILEV)*TX(ILEV-1)     
      !          Y(ILEV) = (K(ILEV)**(-1))*X(ILEV)
      !          TX(ILEV) = (X(ILEV)-A(ILEV)*Y(ILEV))
      !    5.   DO  ILEV=NLEV-1,1,-1
      !          Y(ILEV) = Y(ILEV) + AV(PR_SM_,ILEV+1)*Y(ILEV+1)     
      !    6.  Yext    = beta*Yext + alpha*Y(1)
      !
      !    Note: level numbering reversed wrt ref. DD, i.e. 
      !         1..NLEV <=>  (j) <-> 0
      ! 
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

      call mld_baseprec_aply(zone,baseprecv(1),mlprec_wrk(1)%x2l,&
           &  zzero,mlprec_wrk(1)%y2l,&
           &  baseprecv(1)%base_desc,&
           &  trans,work,info)

      if(info /=0) goto 9999

      mlprec_wrk(1)%tx = mlprec_wrk(1)%x2l

      call psb_spmm(-zone,baseprecv(1)%base_a,mlprec_wrk(1)%y2l,&
           & zone,mlprec_wrk(1)%tx,baseprecv(1)%base_desc,info,work=work)
      if(info /=0) goto 9999

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
          !Smoothed Aggregation
          !
          call psb_halo(mlprec_wrk(ilev-1)%tx,baseprecv(ilev-1)%base_desc,&
               & info,work=work) 
          if(info /=0) goto 9999

          call psb_csmm(zone,baseprecv(ilev)%av(mld_sm_pr_t_),mlprec_wrk(ilev-1)%tx,zzero,&
               & mlprec_wrk(ilev)%x2l,info)
          if(info /=0) goto 9999

        else
          !
          ! Raw  aggregation, may take shortcuts
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


        call mld_baseprec_aply(zone,baseprecv(ilev),mlprec_wrk(ilev)%x2l,&
             & zzero,mlprec_wrk(ilev)%y2l,baseprecv(ilev)%desc_data, 'N',work,info)

        if(info /=0) goto 9999

        if(ilev < nlev) then
          mlprec_wrk(ilev)%tx = mlprec_wrk(ilev)%x2l
          call psb_spmm(-zone,baseprecv(ilev)%base_a,mlprec_wrk(ilev)%y2l,&
               & zone,mlprec_wrk(ilev)%tx,baseprecv(ilev)%base_desc,info,work=work)
          if(info /=0) goto 9999
        endif

      enddo

      do ilev = nlev-1, 1, -1

        ismth = baseprecv(ilev+1)%iprcparm(mld_aggr_kind_)
        n_row = psb_cd_get_local_rows(baseprecv(ilev)%base_desc)

        if (ismth  /= mld_no_smooth_) then 

          if (ismth == mld_smooth_prol_) &
               & call psb_halo(mlprec_wrk(ilev+1)%y2l,&
               & baseprecv(ilev+1)%desc_data,info,work=work) 
          call psb_csmm(zone,baseprecv(ilev+1)%av(mld_sm_pr_),mlprec_wrk(ilev+1)%y2l,&
               & zone,mlprec_wrk(ilev)%y2l,info)

          if(info /=0) goto 9999

        else

          do i=1, n_row
            mlprec_wrk(ilev)%y2l(i) = mlprec_wrk(ilev)%y2l(i) + &
                 & mlprec_wrk(ilev+1)%y2l(baseprecv(ilev+1)%mlia(i))
          enddo

        end if

      enddo

      call psb_geaxpby(alpha,mlprec_wrk(1)%y2l,beta,y,&
           &  baseprecv(1)%base_desc,info)

      if(info /=0) goto 9999



    case(mld_twoside_smooth_)

      !
      !    Symmetrized  smoothing. 
      !    1.   X(1)  = Xext
      !    2.   Y(1)  = (K(1)**(-1))*X(1)
      !    3.   TX(1) = X(1) - A(1)*Y(1)
      !    4.   DO ILEV=2, NLEV 
      !          X(ILEV) = AV(PR_SM_T_,ILEV)*TX(ILEV-1)     
      !          Y(ILEV) = (K(ILEV)**(-1))*X(ILEV)
      !          TX(ILEV) = (X(ILEV)-A(ILEV)*Y(ILEV))
      !    5.   DO  ILEV=NLEV-1,1,-1
      !          Y(ILEV) = Y(ILEV) + AV(PR_SM_,ILEV+1)*Y(ILEV+1)     
      !          Y(ILEV) = Y(ILEV) + (K(ILEV)**(-1))*(X(ILEV)-A(ILEV)*Y(ILEV))
      !    6.  Yext    = beta*Yext + alpha*Y(1)
      !
      !    Note: level numbering reversed wrt ref. DD, i.e. 
      !         1..NLEV <=>  (j) <-> 0
      !
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

      call mld_baseprec_aply(zone,baseprecv(1),mlprec_wrk(1)%x2l,&
           &  zzero,mlprec_wrk(1)%y2l,&
           &  baseprecv(1)%base_desc,&
           &  trans,work,info)

      if(info /=0) goto 9999

      mlprec_wrk(1)%ty = mlprec_wrk(1)%x2l

      call psb_spmm(-zone,baseprecv(1)%base_a,mlprec_wrk(1)%y2l,&
           & zone,mlprec_wrk(1)%ty,baseprecv(1)%base_desc,info,work=work)
      if(info /=0) goto 9999

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
          !Smoothed Aggregation
          !
          call psb_halo(mlprec_wrk(ilev-1)%ty,baseprecv(ilev-1)%base_desc,&
               & info,work=work) 
          if(info /=0) goto 9999

          call psb_csmm(zone,baseprecv(ilev)%av(mld_sm_pr_t_),mlprec_wrk(ilev-1)%ty,zzero,&
               & mlprec_wrk(ilev)%x2l,info)
          if(info /=0) goto 9999

        else
          !
          ! Raw  aggregation, may take shortcuts
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

        call mld_baseprec_aply(zone,baseprecv(ilev),mlprec_wrk(ilev)%x2l,&
             & zzero,mlprec_wrk(ilev)%y2l,baseprecv(ilev)%desc_data, 'N',work,info)

        if(info /=0) goto 9999

        if(ilev < nlev) then
          mlprec_wrk(ilev)%ty = mlprec_wrk(ilev)%x2l
          call psb_spmm(-zone,baseprecv(ilev)%base_a,mlprec_wrk(ilev)%y2l,&
               & zone,mlprec_wrk(ilev)%ty,baseprecv(ilev)%base_desc,info,work=work)
          if(info /=0) goto 9999
        endif

      enddo


      do ilev=nlev-1, 1, -1

        ismth = baseprecv(ilev+1)%iprcparm(mld_aggr_kind_)
        n_row = psb_cd_get_local_rows(baseprecv(ilev)%base_desc)

        if (ismth  /= mld_no_smooth_) then 
          if (ismth == mld_smooth_prol_) &
               & call psb_halo(mlprec_wrk(ilev+1)%y2l,baseprecv(ilev+1)%desc_data,&
               &  info,work=work) 
          call psb_csmm(zone,baseprecv(ilev+1)%av(mld_sm_pr_),mlprec_wrk(ilev+1)%y2l,&
               &  zone,mlprec_wrk(ilev)%y2l,info)
          if(info /=0) goto 9999

        else
          do i=1, n_row
            mlprec_wrk(ilev)%y2l(i) = mlprec_wrk(ilev)%y2l(i) + &
                 & mlprec_wrk(ilev+1)%y2l(baseprecv(ilev+1)%mlia(i))
          enddo

        end if

        call psb_spmm(-zone,baseprecv(ilev)%base_a,mlprec_wrk(ilev)%y2l,&
             &   zone,mlprec_wrk(ilev)%tx,baseprecv(ilev)%base_desc,info,work=work)

        if(info /=0) goto 9999

        call mld_baseprec_aply(zone,baseprecv(ilev),mlprec_wrk(ilev)%tx,&
             & zone,mlprec_wrk(ilev)%y2l,baseprecv(ilev)%base_desc, trans, work,info)

        if(info /=0) goto 9999

      enddo

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

!!$contains
!!$  subroutine mlprec_wrk_free(wrk)
!!$    type(psb_mlprec_wrk_type) :: wrk(:)
!!$    ! This will not be needed when we have allocatables, as 
!!$    ! it is sufficient to deallocate the container, and 
!!$    ! the compiler is supposed to recursively deallocate the 
!!$    ! various components. 
!!$    integer i
!!$
!!$    do i=1, size(wrk)
!!$      if (associated(wrk(i)%tx))  deallocate(wrk(i)%tx)
!!$      if (associated(wrk(i)%ty))  deallocate(wrk(i)%ty)
!!$      if (associated(wrk(i)%x2l)) deallocate(wrk(i)%x2l)
!!$      if (associated(wrk(i)%y2l)) deallocate(wrk(i)%y2l)
!!$      if (associated(wrk(i)%b2l)) deallocate(wrk(i)%b2l)
!!$      if (associated(wrk(i)%tty)) deallocate(wrk(i)%tty)
!!$    end do
!!$  end subroutine mlprec_wrk_free

end subroutine mld_zmlprec_aply

