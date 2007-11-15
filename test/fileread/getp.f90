module getp

  public get_parms
  public pr_usage

contains
  !
  ! get iteration parameters from standard input
  !
  subroutine  get_parms(icontxt,irst,irnum,ntry,nmat,mtrx,rhs,kmethd,nprecs,precs,ipart,&
       & afmt,istopc,itmax,itrace,eps,outf1,outf2)

    use psb_base_mod
    use precd
    implicit none

    integer           :: icontxt
    character(len=20) :: kmethd
    character(len=80) :: outf1, outf2
    character(len=20),allocatable :: mtrx(:), rhs(:)
    type(precdata),allocatable  :: precs(:)
    integer             :: iret, istopc,itmax,itrace,ipart,nmat,nprecs,irst,irnum,ntry
    character(len=1024) :: charbuf
    real(kind(1.d0))    :: eps, omega,thr1,thr2
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
end module getp
