program df_bench
  use psb_base_mod
  use psb_util_mod
  use mld_prec_mod
  use psb_krylov_mod
  implicit none

  ! input parameters
  character(len=20) :: kmethd
  character(len=80) :: outf1, outf2, outf3
  character(len=20), allocatable :: mtrx(:),rhs(:)
  type precdata
     character(len=10)  :: lv1, lv2    ! First and second level prec type
     integer            :: nlev        !
     integer            :: novr        ! number of overlapping levels
     integer            :: restr       ! restriction  over application of as
     integer            :: prol        ! prolongation over application of as
     integer            :: ftype1      ! Factorization type: ILU, SuperLU, UMFPACK. 
     integer            :: fill1       ! Fill-in for factorization 1
     real(psb_dpk_)     :: thr1        ! Threshold for fact. 1 ILU(T)
     integer            :: mltype      ! additive or multiplicative 2nd level prec
     integer            :: aggr        ! local or global aggregation
     integer            :: smthkind    ! smoothing type
     integer            :: cmat        ! coarse mat
     integer            :: smthpos     ! pre, post, both smoothing
     integer            :: glbsmth     ! global smoothing
     integer            :: ftype2      ! Factorization type: ILU, SuperLU, UMFPACK. 
     integer            :: fill2       ! Fill-in for factorization 1
     real(psb_dpk_)     :: thr2        ! Threshold for fact. 1 ILU(T)
     integer            :: jswp        ! Jacobi sweeps
     real(psb_dpk_)     :: omega       ! smoother omega
     character(len=40)  :: descr       ! verbose description of the prec
  end type precdata
  type(precdata), allocatable  :: precs(:)

  ! sparse matrices
  type(psb_dspmat_type) :: a, aux_a

  ! preconditioner data
  type(mld_dprec_type)  :: pre
  integer               :: igsmth, matop, novr


  ! dense matrices
  real(psb_dpk_), allocatable, target        ::  aux_b(:,:), d(:)
  real(psb_dpk_), allocatable , save :: b_col(:), x_col(:), r_col(:), &
       & x_col_glob(:), r_col_glob(:)
  real(psb_dpk_), pointer  :: b_col_glob(:)
  ! communications data structure
  type(psb_desc_type):: desc_a

  ! blacs variables
  integer               :: ictxt, iam, np
  logical               :: out1, out2

  ! solver paramters
  integer            :: iter, itmax, ierr, itrace, ircode, ipart,nlev,&
       & methd, istopc, iprec, ml, irnum, irst, ntry, nmat, ilev,ipsize,asize,cdsize
  real(psb_dpk_)   :: err, eps

  character(len=5)   :: afmt
  character(len=20)  :: name
  integer   :: iparm(20)

  ! other variables
  integer            :: i,info,j,m_problem, nm, nt
  integer            :: internal, m,ii,nnzero, nprecs, pp
  real(psb_dpk_) :: t1, t2, tprec, r_amax, b_amax,&
       &scale,resmx,resmxp, mttot, mtslv, mtprec
  integer :: nrhs, nrow, n_row, dim, nv, ne
  integer, allocatable  :: ipv(:), neigh(:), ivg(:)


  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif


  name='df_bench'
  if(psb_get_errstatus() /= 0) goto 9999
  info=0
  call psb_set_errverbosity(2)

  !
  !  get parameters
  !
  call get_parms(ictxt,irst,irnum,ntry,nmat,mtrx,rhs,kmethd,nprecs,precs,&
       & ipart,afmt,istopc,itmax,itrace,eps,outf1,outf2)

  if(iam == psb_root_) then
    if(outf1 /= 'NONE') then
      open(8,file=outf1,action='write')
      out1=.true.
    else
      out1=.false.
    end if

    if(outf2 /= 'NONE') then
      open(10,file=outf2,action='write')
      out2=.true.
    else
      out2=.false.
    end if
  end if

  do nm=1, nmat

    if(iam == psb_root_) write(*,'(25("=")," ",a20," ",25("="))')mtrx(nm)
    call psb_barrier(ictxt)
    t1 = psb_wtime()  
    ! read the input matrix to be processed and (possibly) the rhs 
    nrhs = 1

    if (iam == psb_root_) then
      call read_mat(mtrx(nm), aux_a, ictxt)

      m_problem = aux_a%m
      call psb_bcast(ictxt,m_problem,psb_root_)

      if(rhs(nm) /= 'none') then
        !  reading an rhs
        call read_rhs(rhs(nm),aux_b,ictxt)
      end if

      if (allocated(aux_b).and.size(aux_b,1)==m_problem) then
        ! if any rhs were present, broadcast the first one
        b_col_glob =>aux_b(:,1)
      else
        allocate(aux_b(m_problem,1), stat=ircode)
        if (ircode /= 0) then
          call psb_errpush(4000,name)
          goto 9999
        endif

        b_col_glob => aux_b(:,1)
        do i=1, m_problem
          b_col_glob(i) = 1.d0
        enddo
      endif
      call psb_bcast(ictxt,b_col_glob(1:m_problem),psb_root_)
    else
      call psb_bcast(ictxt,m_problem,psb_root_)
      allocate(aux_b(m_problem,1), stat=ircode)
      if (ircode /= 0) then
        call psb_errpush(4000,name)
        goto 9999
      endif
      b_col_glob =>aux_b(:,1)
      call psb_bcast(ictxt,b_col_glob(1:m_problem),psb_root_) 
    end if


    ! switch over different partition types
    if (ipart == 0) then 
      call psb_barrier(ictxt)
      !        if (iam == psb_root_) write(*,'("Partition type: block")')
      allocate(ivg(m_problem),ipv(np))
      if (.true.) then 
        do i=1,m_problem
          call part_block(i,m_problem,np,ipv,nv)
          ivg(i) = ipv(1)
        enddo
        call psb_matdist(aux_a, a, ivg, ictxt, &
             & desc_a,b_col_glob,b_col,info,fmt=afmt)
      else
        call psb_matdist(aux_a, a, part_block, ictxt, &
             & desc_a,b_col_glob,b_col,info,fmt=afmt)
      end if
    else if (ipart == 2) then 
      if (iam == psb_root_) then 
        call build_mtpart(aux_a%m,aux_a%fida,aux_a%ia1,aux_a%ia2,np)
      endif
      call psb_barrier(ictxt)
      call distr_mtpart(0,ictxt)
      call getv_mtpart(ivg)
      call psb_matdist(aux_a, a, ivg, ictxt, &
           & desc_a,b_col_glob,b_col,info,fmt=afmt)
      call free_part(info)
    else 
      call psb_matdist(aux_a, a, part_block, ictxt, &
           & desc_a,b_col_glob,b_col,info,fmt=afmt)
    end if
!!$    open(20+iam)
!!$    call psb_cdprt(20+iam,desc_a,short=.false.)
!!$    close(20+iam)
!!$    write(0,*) iam,'After CDPRT '
!!$    call flush(0)
!!$    call flush(6)
!!$    call psb_barrier(ictxt)

    call psb_geall(x_col,desc_a,info)
    x_col(:) =0.0
    call psb_geasb(x_col,desc_a,info)
    call psb_geall(r_col,desc_a,info)
    r_col(:) =0.0
    call psb_geasb(r_col,desc_a,info)
    t2 = psb_wtime() - t1

!!$    call psb_csprt(10+iam,a,head='% (A)')    
    call psb_sp_free(aux_a,info)

    call psb_amx(ictxt, t2)

!!$    call psb_csprt(20+iam,a,head='% (A)')    
    !
    !  prepare the preconditioning matrix. note the availability
    !  of optional parameters
    !


    do pp=1, nprecs

      mttot=1.d300

      do nt=1,ntry

        if (precs(pp)%lv2(1:2) == 'ml') then
          if (precs(pp)%nlev < 2) then 
            write(0,*) 'Inconsistent number of levels ',precs(pp)%nlev,&
                 & ' forcing 2'
            precs(pp)%nlev = 2
          end if
          nlev = precs(pp)%nlev
          call mld_precinit(pre,precs(pp)%lv2,info,nlev=nlev)
          ! Defaults are OK for all intermediate levels. Only fix last level. 
          if (precs(pp)%omega>=0.0) then 
            call mld_precset(pre,mld_aggr_damp_,precs(pp)%omega,info,ilev=nlev)
          end if
          call mld_precset(pre,mld_ml_type_,       precs(pp)%mltype,   info,ilev=nlev)
          call mld_precset(pre,mld_aggr_alg_,      precs(pp)%aggr,     info,ilev=nlev)
          call mld_precset(pre,mld_coarse_mat_,    precs(pp)%cmat,     info,ilev=nlev)
          call mld_precset(pre,mld_smooth_pos_,    precs(pp)%smthpos,  info,ilev=nlev)
          call mld_precset(pre,mld_sub_solve_,     precs(pp)%ftype2,   info,ilev=nlev)
          call mld_precset(pre,mld_sub_fill_in_,   precs(pp)%fill2,    info,ilev=nlev)
          call mld_precset(pre,mld_fact_thrs_,     precs(pp)%thr2,     info,ilev=nlev)
          call mld_precset(pre,mld_smooth_sweeps_, precs(pp)%jswp,     info,ilev=nlev)
          call mld_precset(pre,mld_aggr_kind_,     precs(pp)%smthkind, info,ilev=nlev)
        else
          call mld_precinit(pre,precs(pp)%lv1,info)
        end if
        call mld_precset(pre,mld_n_ovr_,       precs(pp)%novr,   info,ilev=1)
        call mld_precset(pre,mld_sub_restr_,   precs(pp)%restr,  info,ilev=1)
        call mld_precset(pre,mld_sub_prol_,    precs(pp)%prol,   info,ilev=1)
        call mld_precset(pre,mld_sub_solve_,   precs(pp)%ftype1, info,ilev=1)
        call mld_precset(pre,mld_sub_fill_in_, precs(pp)%fill1,  info,ilev=1)
        call mld_precset(pre,mld_fact_thrs_,   precs(pp)%thr1,   info,ilev=1)


        !  setting initial guess to zero
        call psb_geaxpby(dzero,b_col,dzero,x_col,desc_a,info)
    
        ! building the preconditioner
!!$        write(0,*) 'Check in main program on hasv in:',allocated(desc_a%hashv)
!!$        call flush(0)
        call psb_barrier(ictxt)
        t1 = psb_wtime()
        call mld_precbld(a,desc_a,pre,info)
        tprec = psb_wtime()-t1
        if (info /= 0) then
          write(0,*) 'INFO from precbld ',info
          call psb_error()
          goto 9999
        end if
        if (psb_get_errstatus() /= 0) then 
          write(0,*) 'INFO from precbld ',info
          call psb_error()
          goto 9999
        end if

!!$        write(0,*) 'Check in main program on hasv out:',allocated(desc_a%hashv)
!!$        call flush(0)
        call psb_amx(ictxt,tprec)
!!$        call psb_csprt(40+iam,a,head='% (A)')    
!!$        if (iam == psb_root_) then 
!!$          write(*,'("Matrix : ",a)') mtrx(nm)
!!$          write(*,'("RHS    : ",a)') rhs(nm)
!!$          write(*,'("Method : ",a)') kmethd
!!$          write(*,'("Preconditioner : ",a)') precs(pp)%descr
!!$          call mld_prec_descr(6,pre)
!!$          call flush(6)
!!$        end if
        iparm = 0
        call psb_barrier(ictxt)
        t1 = psb_wtime()
        call  psb_krylov(kmethd,a,pre,b_col,x_col,eps,desc_a,info,& 
             & itmax=itmax,iter=iter,err=err,itrace=itrace,&
             & irst=irst,istop=istopc)     
        call psb_barrier(ictxt)
        t2 = psb_wtime() - t1
        call psb_amx(ictxt,t2)
        if (info/=0) then 
          write(0,*) 'INFO from solver ',info
          call psb_errpush(4010,name,a_err='psb_SOLVER')
          goto 9999
        end if
!!$        write(0,*) iam,'Done Solver'
!!$        call flush(0)
!!$        call flush(6)
!!$        call psb_barrier(ictxt)
        if(iam == psb_root_.and.out2) &
             & write(10,'(a20,2(1x,i3),1x,i5,3(1x,g9.4),1x,a8,1x,a)') &
             & mtrx(nm),np,precs(pp)%novr,iter,tprec,t2,t2+tprec,&
             & trim(kmethd),trim(precs(pp)%descr)
        if(iam == psb_root_) &
             & write(0,'(a20,2(1x,i3),1x,i5,3(1x,g9.4),1x,a8,1x,a)') &
             & mtrx(nm),np,precs(pp)%novr,iter,tprec,t2,t2+tprec,&
             & trim(kmethd),trim(precs(pp)%descr)
        if (nt.lt.ntry) call mld_precfree(pre,info)
        if((t2+tprec).lt.mttot) then
          mtslv=t2
          mtprec=tprec
          mttot=t2+tprec
        end if
      end do
!!$      call psb_csprt(50+iam,a,head='% (A)')    
!!$      write(0,*) 'Check hashv after precfree:',allocated(desc_a%hashv)
!!$      call flush(0)

      call psb_geaxpby(done,b_col,dzero,r_col,desc_a,info)
      call psb_spmm(-done,a,x_col,done,r_col,desc_a,info)
      call psb_genrm2s(resmx,r_col,desc_a,info)
      call psb_geamaxs(resmxp,r_col,desc_a,info)
      
      ipsize = mld_sizeof(pre)
      asize  = psb_sizeof(a)
      cdsize = psb_sizeof(desc_a)
      call psb_sum(ictxt,ipsize)
      call psb_sum(ictxt,asize)
      call psb_sum(ictxt,cdsize)

      if (iam == psb_root_) then 
        write(*,'("Matrix : ",a)') mtrx(nm)
        write(*,'("RHS    : ",a)') rhs(nm)
        write(*,'("Method : ",a)') kmethd
        write(*,'("Preconditioner : ",a)') precs(pp)%descr
        call mld_prec_descr(pre)
        write(*,'("Computed solution on ",i4," processors")')np
        write(*,'(" ")')
        write(*,'("Iterations to convergence: ",i6)')  iter
        write(*,'("Error indicator on exit  : ",g9.4)') err
        write(*,'("Time to buil prec.       : ",es10.4)')mtprec
        write(*,'("Time to solve matrix     : ",es10.4)')mtslv
        write(*,'("Time per iteration       : ",es10.4)')mtslv/(iter)
        write(*,'("Total time               : ",es10.4)')mttot
        write(*,'("Residual norm 2          : ",es10.4)')resmx
        write(*,'("Residual norm inf        : ",es10.4)')resmxp
        write(*,'("Total memory occupation for A:      ",i10)')asize
        write(*,'("Total memory occupation for DESC_A: ",i10)')cdsize
        write(*,'("Total memory occupation for PRE:    ",i10)')ipsize

        write(*,'(72("="))')
        write(*,'(" ")')
        write(*,'(" ")')
        write(*,'(" ")')

        if(out1) write(8,'(a20,2(1x,i3),1x,i5,5(1x,g9.4),1x,a8,1x,a)') mtrx(nm),&
             & np,precs(pp)%novr,&
             & iter,mtprec,mtslv,mttot,resmx,resmxp,&
             & trim(kmethd),trim(precs(pp)%descr)
      end if

      call mld_precfree(pre,info)

      if (.false.) then 
        allocate(x_col_glob(m_problem),r_col_glob(m_problem),stat=ierr)
        if (ierr /= 0) then 
          write(0,*) 'allocation error: no data collection'
        else
          call psb_gather(x_col_glob,x_col,desc_a,info,root=psb_root_)
          call psb_gather(r_col_glob,r_col,desc_a,info,root=psb_root_)
          if (iam == psb_root_) then
            write(0,'(" ")')
            write(0,'("Saving x on file")')
            write(outf3,'(a,a,a)')trim(mtrx(nm)),'.psb_sol.',&
                 & psb_tolower(trim(precs(pp)%descr))
            open(20,file=outf3)
            write(20,*) 'matrix: ',mtrx(nm)
            write(20,*) 'computed solution on ',np,' processors.'
            write(20,*) 'iterations to convergence: ',iter
            write(20,*) 'error indicator (infinity norm) on exit:', &
                 & ' ||r||/(||a||||x||+||b||) = ',err
            write(20,*) 'max residual = ',resmx, resmxp
            do i=1,m_problem
              write(20,998) i,x_col_glob(i),b_col_glob(i)
            enddo
            close(20)
          end if
        end if
      end if
998   format(i8,4(2x,g20.14))
993   format(i6,4(1x,e12.6))
      
      
!!$      if (pp == 1) call psb_csprt(40+iam,a,head='% (A)')    
    end do

    call psb_gefree(b_col, desc_a,info)
    call psb_gefree(x_col, desc_a,info)
    call psb_spfree(a, desc_a,info)
    call psb_cdfree(desc_a,info)
    deallocate(r_col,stat=info)
    deallocate(aux_b,stat=info)
    if (ipart==0) then 
      deallocate(ivg,ipv,stat=info)
    endif
  end do

9999 continue
  if(info /= 0) then
    call psb_error(ictxt)
  end if
  call psb_exit(ictxt)
  stop


contains
  !
  ! get iteration parameters from standard input
  !
  subroutine  get_parms(icontxt,irst,irnum,ntry,nmat,mtrx,rhs,kmethd,nprecs,precs,ipart,&
       & afmt,istopc,itmax,itrace,eps,outf1,outf2)

    use psb_base_mod
    implicit none

    integer           :: icontxt
    character(len=20) :: kmethd
    character(len=80) :: outf1, outf2
    character(len=20),allocatable :: mtrx(:), rhs(:)
    type(precdata),allocatable  :: precs(:)
    integer             :: iret, istopc,itmax,itrace,ipart,nmat,nprecs,irst,irnum,ntry
    character(len=1024) :: charbuf
    real(psb_dpk_)    :: eps, omega,thr1,thr2
    character           :: afmt*5, lv1*10, lv2*10, pdescr*40
    integer             :: iam, nm, np, i, idx
    integer, parameter  :: npparms=14
    integer             :: inparms(40), ip, pparms(npparms)

    call psb_info(icontxt,iam,np)

    if (iam==psb_root_) then
      ! read input parameters
      read(*,*) outf1
      read(*,*) outf2
      read(*,*) kmethd
      read(*,*) eps
      read(*,*) afmt

      call psb_bcast(icontxt,kmethd)
      call psb_bcast(icontxt,eps,0)

      call psb_bcast(icontxt,afmt)

      read(*,*) ipart
      read(*,*) itmax
      read(*,*) itrace
      read(*,*) istopc
      read(*,*) irst
      read(*,*) irnum
      read(*,*) ntry
      read(*,*) nprecs
      ! broadcast parameters to all processors    

      inparms(1) = ipart
      inparms(2) = itmax
      inparms(3) = itrace
      inparms(4) = istopc
      inparms(5) = irst
      inparms(6) = irnum
      inparms(7) = ntry    
      call psb_bcast(icontxt,inparms(1:7))

      call psb_bcast(icontxt,nprecs)

      allocate(precs(nprecs))

      do np=1,nprecs
        read(*,'(a)')charbuf
        charbuf = adjustl(charbuf)
        idx=index(charbuf," ")
        read(charbuf(1:idx-1),'(a)')lv1
        charbuf=adjustl(charbuf(idx:))
        idx=index(charbuf," ")
        read(charbuf(1:idx-1),'(a)')lv2
        charbuf=adjustl(charbuf(idx:))
        do i=1, npparms
          idx=index(charbuf," ")
          read(charbuf(1:idx),*) pparms(i)
          charbuf=adjustl(charbuf(idx:))
        end do
        idx=index(charbuf," ")
        read(charbuf(1:idx),*) omega
        charbuf=adjustl(charbuf(idx:))
        idx=index(charbuf," ")
        read(charbuf(1:idx),*) thr1
        charbuf=adjustl(charbuf(idx:))
        idx=index(charbuf," ")
        read(charbuf(1:idx),*) thr2
        charbuf=adjustl(charbuf(idx:))
        read(charbuf,'(a)') pdescr

        call psb_bcast(icontxt,pdescr)
        precs(np)%descr=pdescr

        call psb_bcast(icontxt,lv1)
        call psb_bcast(icontxt,lv2)
        call psb_bcast(icontxt,pparms(1:npparms))
        call psb_bcast(icontxt,omega)
        call psb_bcast(icontxt,thr1)
        call psb_bcast(icontxt,thr2)

        precs(np)%lv1      = lv1
        precs(np)%lv2      = lv2
        precs(np)%novr     = pparms(1)
        precs(np)%restr    = pparms(2)
        precs(np)%prol     = pparms(3)
        precs(np)%ftype1   = pparms(4)
        precs(np)%fill1    = pparms(5)
        precs(np)%mltype   = pparms(6)
        precs(np)%aggr     = pparms(7)
        precs(np)%smthkind = pparms(8)
        precs(np)%cmat     = pparms(9) 
        precs(np)%smthpos  = pparms(10)
        precs(np)%ftype2   = pparms(11)
        precs(np)%fill2    = pparms(12)
        precs(np)%jswp     = pparms(13)
        precs(np)%nlev     = pparms(14)
        precs(np)%omega    = omega
        precs(np)%thr1     = thr1
        precs(np)%thr2     = thr2
      end do

      read(*,*) nmat
      call psb_bcast(icontxt,nmat)
      allocate(mtrx(nmat),rhs(nmat))

      do nm=1, nmat
        read(*,'(a)') charbuf
        charbuf=adjustl(charbuf)
        idx=index(charbuf," ")
        mtrx(nm)=charbuf(1:idx-1)
        rhs(nm)=adjustl(charbuf(idx+1:))
        call psb_bcast(icontxt,mtrx(nm))
        call psb_bcast(icontxt,rhs(nm))
      end do

    else
      ! receive parameters
      call psb_bcast(icontxt,kmethd)
      call psb_bcast(icontxt,eps)     

      call psb_bcast(icontxt,afmt)

      call psb_bcast(icontxt,inparms(1:7))

      ipart   =  inparms(1) 
      itmax   =  inparms(2) 
      itrace  =  inparms(3) 
      istopc  =  inparms(4) 
      irst    =  inparms(5) 
      irnum   =  inparms(6) 
      ntry    =  inparms(7) 

      call psb_bcast(icontxt,nprecs)
      allocate(precs(nprecs))

      do np=1,nprecs
        call psb_bcast(icontxt,pdescr)
        precs(np)%descr=pdescr

        call psb_bcast(icontxt,lv1)

        call psb_bcast(icontxt,lv2)

        call psb_bcast(icontxt,pparms(1:npparms))
        call psb_bcast(icontxt,omega)     
        call psb_bcast(icontxt,thr1)
        call psb_bcast(icontxt,thr2)

        precs(np)%lv1      = lv1
        precs(np)%lv2      = lv2
        precs(np)%novr     = pparms(1)
        precs(np)%restr    = pparms(2)
        precs(np)%prol     = pparms(3)
        precs(np)%ftype1   = pparms(4)
        precs(np)%fill1    = pparms(5)
        precs(np)%mltype   = pparms(6)
        precs(np)%aggr     = pparms(7)
        precs(np)%smthkind = pparms(8)
        precs(np)%cmat     = pparms(9) 
        precs(np)%smthpos  = pparms(10)
        precs(np)%ftype2   = pparms(11)
        precs(np)%fill2    = pparms(12)
        precs(np)%jswp     = pparms(13)
        precs(np)%nlev     = pparms(14)
        precs(np)%omega    = omega
        precs(np)%thr1     = thr1
        precs(np)%thr2     = thr2
      end do


      call psb_bcast(icontxt,nmat)
      allocate(mtrx(nmat),rhs(nmat))

      do nm=1,nmat

        call psb_bcast(icontxt,mtrx(nm))
        call psb_bcast(icontxt,rhs(nm))

      end do

    end if

  end subroutine get_parms
  subroutine pr_usage(iout)
    integer iout
    write(iout, *) ' number of parameters is incorrect!'
    write(iout, *) ' use: hb_sample mtrx_file methd prec [ptype &
         &itmax istopc itrace]' 
    write(iout, *) ' where:'
    write(iout, *) '     mtrx_file      is stored in hb format'
    write(iout, *) '     methd          may be: cgstab '
    write(iout, *) '     itmax          max iterations [500]        '
    write(iout, *) '     istopc         stopping criterion [1]      '
    write(iout, *) '     itrace         0  (no tracing, default) or '
    write(iout, *) '                    >= 0 do tracing every itrace'
    write(iout, *) '                    iterations ' 
    write(iout, *) '     prec           may be: ilu diagsc none'
    write(iout, *) '     ptype          partition strategy default 0'
    write(iout, *) '                    0: block partition '
  end subroutine pr_usage

end program df_bench





