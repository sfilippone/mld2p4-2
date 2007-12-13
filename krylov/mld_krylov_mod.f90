!!$
!!$ 
!!$                                MLD2P4
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS v.2.0)
!!$  
!!$  (C) Copyright 2006  Alfredo Buttari      University of Rome Tor Vergata
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
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
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
! This is a partial duplication of material already present in the Krylov 
! part of the PSBLAS subtree. This is not as bad as it sounds, since 
! 1. There is no actual computation here, these are just wrappers around 
!    the actual methods, so very little code is necessary;
! 2. Changing the wrappers enables issuing error messages with a subroutine 
!    name coherent with what the user sees;
! 3. At the very least we would need a layer renaming the PSBLAS krylov interfaces,
!    so we are doing approx. the same amount of (re)coding. 
!
Module mld_krylov_mod


  interface mld_krylov
    module procedure mld_dkrylov, mld_zkrylov
  end interface

  interface psb_cg
     subroutine psb_dcg(a,prec,b,x,eps,&
          & desc_a,info,itmax,iter,err,itrace,istop)
       use psb_base_mod
       use psb_prec_mod
       type(psb_dspmat_type), intent(in)  :: a
       type(psb_desc_type), intent(in)    :: desc_a
       real(kind(1.d0)), intent(in)       :: b(:)
       real(kind(1.d0)), intent(inout)    :: x(:)
       real(kind(1.d0)), intent(in)       :: eps
       type(psb_dprec_type), intent(in)   :: prec
       integer, intent(out)               :: info
       integer, optional, intent(in)      :: itmax, itrace,istop
       integer, optional, intent(out)     :: iter
       real(kind(1.d0)), optional, intent(out) :: err
     end subroutine psb_dcg
  end interface

  interface psb_bicg
     subroutine psb_dbicg(a,prec,b,x,eps,&
          & desc_a,info,itmax,iter,err,itrace,istop)
       use psb_base_mod
       use psb_prec_mod
       type(psb_dspmat_type), intent(in)  :: a
       type(psb_desc_type), intent(in)    :: desc_a
       real(kind(1.d0)), intent(in)       :: b(:)
       real(kind(1.d0)), intent(inout)    :: x(:)
       real(kind(1.d0)), intent(in)       :: eps
       type(psb_dprec_type), intent(in)   :: prec
       integer, intent(out)               :: info
       integer, optional, intent(in)      :: itmax, itrace,istop
       integer, optional, intent(out)     :: iter
       real(kind(1.d0)), optional, intent(out) :: err
     end subroutine psb_dbicg
  end interface

  interface psb_bicgstab
     subroutine psb_dcgstab(a,prec,b,x,eps,&
          & desc_a,info,itmax,iter,err,itrace,istop)
       use psb_base_mod
       use psb_prec_mod
       type(psb_dspmat_type), intent(in)  :: a
       type(psb_desc_type), intent(in)    :: desc_a
       real(kind(1.d0)), intent(in)       :: b(:)
       real(kind(1.d0)), intent(inout)    :: x(:)
       real(kind(1.d0)), intent(in)       :: eps
       type(psb_dprec_type), intent(in)   :: prec
       integer, intent(out)               :: info
       integer, optional, intent(in)      :: itmax, itrace,istop
       integer, optional, intent(out)     :: iter
       real(kind(1.d0)), optional, intent(out) :: err
     end subroutine psb_dcgstab
     subroutine psb_zcgstab(a,prec,b,x,eps,&
          & desc_a,info,itmax,iter,err,itrace,istop)
       use psb_base_mod
       use psb_prec_mod
       type(psb_zspmat_type), intent(in)  :: a
       type(psb_desc_type), intent(in)    :: desc_a
       complex(kind(1.d0)), intent(in)       :: b(:)
       complex(kind(1.d0)), intent(inout)    :: x(:)
       real(kind(1.d0)), intent(in)       :: eps
       type(psb_zprec_type), intent(in)   :: prec
       integer, intent(out)               :: info
       integer, optional, intent(in)      :: itmax, itrace,istop
       integer, optional, intent(out)     :: iter
       real(kind(1.d0)), optional, intent(out) :: err
     end subroutine psb_zcgstab
  end interface

  interface psb_bicgstabl
    Subroutine psb_dcgstabl(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err, itrace,irst,istop)
      use psb_base_mod
      use psb_prec_mod
      Type(psb_dspmat_type), Intent(in)  :: a
      Type(psb_desc_type), Intent(in)    :: desc_a
      type(psb_dprec_type), intent(in)   :: prec
      Real(Kind(1.d0)), Intent(in)       :: b(:)
      Real(Kind(1.d0)), Intent(inout)    :: x(:)
      Real(Kind(1.d0)), Intent(in)       :: eps
      integer, intent(out)               :: info
      Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
      Integer, Optional, Intent(out)     :: iter
      Real(Kind(1.d0)), Optional, Intent(out) :: err
    end subroutine psb_dcgstabl
  end interface

  interface psb_rgmres
    Subroutine psb_dgmresr(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,irst,istop)
      use psb_base_mod
      use psb_prec_mod
      Type(psb_dspmat_type), Intent(in)  :: a
      Type(psb_desc_type), Intent(in)    :: desc_a
      type(psb_dprec_type), intent(in)   :: prec 
      Real(Kind(1.d0)), Intent(in)       :: b(:)
      Real(Kind(1.d0)), Intent(inout)    :: x(:)
      Real(Kind(1.d0)), Intent(in)       :: eps
      integer, intent(out)               :: info
      Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
      Integer, Optional, Intent(out)     :: iter
      Real(Kind(1.d0)), Optional, Intent(out) :: err
    end subroutine psb_dgmresr
    Subroutine psb_zgmresr(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,irst,istop)
      use psb_base_mod
      use psb_prec_mod
      Type(psb_zspmat_type), Intent(in)  :: a
      Type(psb_desc_type), Intent(in)    :: desc_a
      type(psb_zprec_type), intent(in)   :: prec 
      complex(Kind(1.d0)), Intent(in)    :: b(:)
      complex(Kind(1.d0)), Intent(inout) :: x(:)
      Real(Kind(1.d0)), Intent(in)       :: eps
      integer, intent(out)               :: info
      Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
      Integer, Optional, Intent(out)     :: iter
      Real(Kind(1.d0)), Optional, Intent(out) :: err
    end subroutine psb_zgmresr
  end interface

  interface psb_cgs
    subroutine psb_dcgs(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
      use psb_base_mod
      use psb_prec_mod
      type(psb_dspmat_type), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a 
      type(psb_dprec_type), intent(in)   :: prec 
      real(kind(1.d0)), intent(in)       :: b(:)
      real(kind(1.d0)), intent(inout)    :: x(:)
      real(kind(1.d0)), intent(in)       :: eps
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(kind(1.d0)), optional, intent(out) :: err
    end subroutine psb_dcgs
     subroutine psb_zcgs(a,prec,b,x,eps,&
          & desc_a,info,itmax,iter,err,itrace,istop)
       use psb_base_mod
       use psb_prec_mod
       type(psb_zspmat_type), intent(in)  :: a
       type(psb_desc_type), intent(in)    :: desc_a
       complex(kind(1.d0)), intent(in)       :: b(:)
       complex(kind(1.d0)), intent(inout)    :: x(:)
       real(kind(1.d0)), intent(in)       :: eps
       type(psb_zprec_type), intent(in)   :: prec
       integer, intent(out)               :: info
       integer, optional, intent(in)      :: itmax, itrace,istop
       integer, optional, intent(out)     :: iter
       real(kind(1.d0)), optional, intent(out) :: err
     end subroutine psb_zcgs
  end interface
  
contains


  !
  ! File: mld_krylov_mod.f90
  !
  ! Subroutine: mld_dkrylov
  ! 
  !    Front-end for the Krylov subspace iterations, real version
  !    
  ! Arguments:
  !
  !    methd  -  character                    The specific method; can take the values:
  !                                           CG
  !                                           CGS
  !                                           BICG
  !                                           BICGSTAB
  !                                           BICGSTABL
  !                                           RGMRES
  !                                           
  !    a      -  type(<psb_dspmat_type>)      Input: sparse matrix containing A.
  !    prec   -  type(<mld_dprec_type>)       Input: preconditioner
  !    b      -  real,dimension(:)            Input: vector containing the
  !                                           right hand side B
  !    x      -  real,dimension(:)            Input/Output: vector containing the
  !                                           initial guess and final solution X.
  !    eps    -  real                         Input: Stopping tolerance; the iteration is
  !                                           stopped when the error estimate
  !                                           |err| <= eps
  !    desc_a -  type(<psb_desc_type>).       Input: The communication descriptor.
  !    info   -  integer.                     Output: Return code
  !
  !    itmax  -  integer(optional)            Input: maximum number of iterations to be
  !                                           performed.
  !    iter   -  integer(optional)            Output: how many iterations have been
  !                                           performed.
  !    err    -  real   (optional)            Output: error estimate on exit
  !    itrace -  integer(optional)            Input: print an informational message
  !                                           with the error estimate every itrace
  !                                           iterations
  !    irst   -  integer(optional)            Input: restart parameter for RGMRES and 
  !                                           BICGSTAB(L) methods
  !    istop  -  integer(optional)            Input: stopping criterion, or how
  !                                           to estimate the error. 
  !                                           1: err =  |r|/|b|
  !                                           2: err =  |r|/(|a||x|+|b|)
  !                                           where r is the (preconditioned, recursive
  !                                           estimate of) residual 
  ! 
  Subroutine mld_dkrylov(method,a,prec,b,x,eps,desc_a,info,&
       &itmax,iter,err,itrace,irst,istop)

    use psb_base_mod
    use psb_prec_mod

    character(len=*)                   :: method
    Type(psb_dspmat_type), Intent(in)  :: a
    Type(psb_desc_type), Intent(in)    :: desc_a
    type(psb_dprec_type), intent(in)   :: prec 
    Real(Kind(1.d0)), Intent(in)       :: b(:)
    Real(Kind(1.d0)), Intent(inout)    :: x(:)
    Real(Kind(1.d0)), Intent(in)       :: eps
    integer, intent(out)               :: info
    Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
    Integer, Optional, Intent(out)     :: iter
    Real(Kind(1.d0)), Optional, Intent(out) :: err

    integer                            :: ictxt,me,np,err_act
    character(len=20)             :: name

    info = 0
    name = 'mld_krylov'
    call psb_erractionsave(err_act)

    ictxt=psb_cd_get_context(desc_a)
    
    call psb_info(ictxt, me, np)

    select case(toupper(method))
    case('CG') 
      call  psb_cg(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
    case('CGS') 
      call  psb_cgs(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
    case('BICG') 
      call  psb_bicg(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
    case('BICGSTAB') 
      call  psb_bicgstab(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
    case('RGMRES')
      call  psb_rgmres(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,irst,istop)
    case('BICGSTABL')
      call  psb_bicgstabl(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,irst,istop)
    case default
      if (me==0) write(0,*) 'Warning: Unknown method  ',method,&
           & ' in ',name,'  defaulting to BiCGSTAB'
      call  psb_bicgstab(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
    end select

    if(info/=0) then
      call psb_errpush(info,name)
      goto 9999
    end if
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    
  end subroutine mld_dkrylov



  !
  ! File: mld_krylov_mod.f90
  !
  ! Subroutine: mld_zkrylov
  ! 
  !    Front-end for the Krylov subspace iterations, complexversion
  !    
  ! Arguments:
  !
  !    methd  -  character                    The specific method; can take the values:
  !                                           CGS
  !                                           BICGSTAB
  !                                           RGMRES
  !                                           
  !    a      -  type(<psb_zspmat_type>)      Input: sparse matrix containing A.
  !    prec   -  type(<mld_zprec_type>)       Input: preconditioner
  !    b      -  complex,dimension(:)         Input: vector containing the
  !                                           right hand side B
  !    x      -  complex,dimension(:)         Input/Output: vector containing the
  !                                           initial guess and final solution X.
  !    eps    -  real                         Input: Stopping tolerance; the iteration is
  !                                           stopped when the error estimate
  !                                           |err| <= eps
  !    desc_a -  type(<psb_desc_type>).       Input: The communication descriptor.
  !    info   -  integer.                     Output: Return code
  !
  !    itmax  -  integer(optional)            Input: maximum number of iterations to be
  !                                           performed.
  !    iter   -  integer(optional)            Output: how many iterations have been
  !                                           performed.
  !    err    -  real   (optional)            Output: error estimate on exit
  !    itrace -  integer(optional)            Input: print an informational message
  !                                           with the error estimate every itrace
  !                                           iterations
  !    irst   -  integer(optional)            Input: restart parameter for RGMRES and 
  !                                           BICGSTAB(L) methods
  !    istop  -  integer(optional)            Input: stopping criterion, or how
  !                                           to estimate the error. 
  !                                           1: err =  |r|/|b|
  !                                           2: err =  |r|/(|a||x|+|b|)
  !                                           where r is the (preconditioned, recursive
  !                                           estimate of) residual 
  ! 
  Subroutine mld_zkrylov(method,a,prec,b,x,eps,desc_a,info,&
       &itmax,iter,err,itrace,irst,istop)
    use psb_base_mod
    use psb_prec_mod
    character(len=*)                   :: method
    Type(psb_zspmat_type), Intent(in)  :: a
    Type(psb_desc_type), Intent(in)    :: desc_a
    type(psb_zprec_type), intent(in)   :: prec 
    complex(Kind(1.d0)), Intent(in)    :: b(:)
    complex(Kind(1.d0)), Intent(inout) :: x(:)
    Real(Kind(1.d0)), Intent(in)       :: eps
    integer, intent(out)               :: info
    Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
    Integer, Optional, Intent(out)     :: iter
    Real(Kind(1.d0)), Optional, Intent(out) :: err

    integer                            :: ictxt,me,np,err_act
    character(len=20)             :: name

    info = 0
    name = 'mld_krylov'
    call psb_erractionsave(err_act)


    ictxt=psb_cd_get_context(desc_a)
    
    call psb_info(ictxt, me, np)
    

    select case(toupper(method))
!!$    case('CG') 
!!$      call  psb_cg(a,prec,b,x,eps,desc_a,info,&
!!$         &itmax,iter,err,itrace,istop)
    case('CGS') 
      call  psb_cgs(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
!!$    case('BICG') 
!!$      call  psb_bicg(a,prec,b,x,eps,desc_a,info,&
!!$         &itmax,iter,err,itrace,istop)
    case('BICGSTAB') 
      call  psb_bicgstab(a,prec,b,x,eps,desc_a,info,&
           & itmax,iter,err,itrace,istop)
    case('RGMRES')
      call  psb_rgmres(a,prec,b,x,eps,desc_a,info,&
           & itmax,iter,err,itrace,irst,istop)
!!$    case('BICGSTABL')
!!$      call  psb_bicgstabl(a,prec,b,x,eps,desc_a,info,&
!!$         &itmax,iter,err,itrace,irst,istop)
    case default
      if (me==0) write(0,*) 'Warning: Unknown method ',method,&
           & ' in ',name,' defaulting to BiCGSTAB'
      call  psb_bicgstab(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
    end select
    
    if(info/=0) then
      call psb_errpush(info,name)
      goto 9999
    end if
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if

  end subroutine mld_zkrylov



end module mld_krylov_mod


  
