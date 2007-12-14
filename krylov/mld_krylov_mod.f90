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
! File: mld_krylov_mod.f90.
!
! Module:   mld_krylov_mod.
! Contains: mld_dkrylov, mld_zkrylov.
!
!  This module defines the interfaces to the real and complex versions of the
!  routines implementing Krylov solvers:
!  - mld_krylov    : driver providing a general interface to the Krylov
!                    solvers available in PSBLAS;
!  - psb_cg        : Conjugate Gradient (CG) solver (real version only);
!  - psb_bicg      : BiConjugate Gradient (BiCG) solver (real version only);
!  - psb_bicgstab  : BiConjugate Gradient Stabilized (BiCGStab) solver;
!  - psb_bicgstabl : BiCGStab(l) solver (real version only);
!  - psb_rgmres         : Generalized Minimal Residual (GMRES) solver with
!                    Restarting;
!  - psb_cgs         : Conjugate Gradient Squared (CGS) solver. 
!
!  This is a partial duplication of material already present in the Krylov 
!  part of the PSBLAS subtree. This is not bad since 
!  1. there is no actual computation here, these are just wrappers around 
!     the actual methods, so very little code is necessary;
!  2. these wrappers enable issuing error messages with a subroutine name 
!     coherent with what the user sees;
!  3. at the very least we would need a layer renaming the PSBLAS krylov
!     interfaces, so we are doing approximately the same amount of (re)coding. 
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
  ! Subroutine: mld_dkrylov/mld_zkrylov.
  ! Version:    real/complex.
  ! Note: internal subroutine of module mld_krylov_mod.
  !
  !       This routine solves a linear system
  !
  !                                Ax = b
  !
  !  by using one of following preconditioned Krylov solvers: CG, BiCG, BiCGStab,
  !  BiCGStab(l), Restarted GMRES, CGS. The preconditioning is applied as right
  !  preconditioning.
  ! 
  !  The user can choose among two stopping criteria:
  !
  !  - one based on the normwise backward error in the infinity norm, i.e.
  !
  !        ||r_i||_inf < eps * (||A||_inf * ||x_i||_inf + ||b||_inf);
  !
  !  - the other based on the relative residual reduction in the 2-norm, i.e.
  !
  !                          ||r_i||_2 < eps * ||b||_2.
  !
  !  In the above formulae, x_i is the approximate solution and r_i = b - A*x_i
  !  is the corresponding residual at the i-th iteration.
  !
  !  This routine is a driver that provides a general interface to the Krylov
  !  solvers implemented in PSBLAS. These solvers can use all the preconditioners
  !  available in MLD2P4.
  !
  !
  ! Arguments:
  !    method  -  character(len=*), input.
  !               The Krylov method to be applied. The following values are
  !               are allowed in the real case: 'cg', 'bicg', 'bicgstab',
  !               'bicgstabl', 'rgmres', 'cgs' (and all the upper/lower case
  !               variants). The same values except 'cg', 'bicg' and 'bicgstabl'
  !               are allowed in the complex case.
  !    a       -  type(<psb_dspmat_type>)/type(<psb_zspmat_type>), input.
  !               The sparse matrix structure containing the local part of the
  !               matrix A.
  !    prec    -  type(<mld_dprec_type>)/type(<mld_zprec_type>), input.
  !               The preconditioner data structure containing the local part
  !               of the preconditioner to be applied.
  !    b       -  real(kind(1.d0))/complex(kind(1.d0)), dimension(:), input.
  !               The local part of the right-hand side.
  !    x       -  real(kind(1.d0))/complex(kind(1.d0)), dimension(:), input/output.
  !               The local part of the approximate solution. In input it
  !               contains an initial guess of the solution, while, in output,
  !               the approximation computed by the selected Krylov      solver.
  !    eps     -  real(kind(1.d0)), input.
  !               The tolerance used in the stopping criterion.
  !    desc_a  -  type(<psb_desc_type>), input.
  !               The communication descriptor associated to a.
  !    info    -  integer, output.
  !               Error code.
  !    itmax   -  integer, optional, input.
  !               The maximum number of iterations to be performed. If itmax
  !               is not present, a maximum number of 1000 iterations is
  !               assumed.
  !    iter    -  integer, optional, output.
  !               The number of iterations actually performed.
  !    err     -  real(kind(1.d0)), optional, output.
  !               An estimate of the error in the computed solution, i.e. the
  !               relative backward error ||r|| / (||A||*||x||+||b||), where
  !               || || is the infinity norm, or the relative residual
  !               ||r||/||b||, where || || is the 2-norm, according to the
  !               selected stopping criterion (see istop below).
  !    itrace  -  integer, optional, input.
  !                          If itrace>0 then information about convergence is printed
  !               every itrace iterations; otherwise it is ignored.      If itrace
  !               is not present, itrace=0 is assumed.
  !    irst    -  integer, optional, input.
  !                          The restart parameter used in the GMRES and BiCGSTAB(l)
  !               solvers. It is ignored if other solvers are chosen. If irst
  !               is not present, irst=10 is assumed for Restarded GMRES and
  !               irst=l=1 for BiCGStab(l) (note that BiCGStab(1)=BiCGStab).
  !    istop   -  integer, optional, input.
  !                          The normwise backward error stopping criterion is used if
  !               istop=1; the relative residual one if istop=2. If istop
  !               is not present, istop=1 is assumed.
  !

  Subroutine mld_dkrylov(method,a,prec,b,x,eps,desc_a,info,&
       &itmax,iter,err,itrace,irst,istop)

    use psb_base_mod
    use psb_prec_mod

  ! Arguments
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

  ! Local variables
    integer                            :: ictxt,me,np,err_act
    character(len=20)                  :: name

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

  Subroutine mld_zkrylov(method,a,prec,b,x,eps,desc_a,info,&
       &itmax,iter,err,itrace,irst,istop)

    use psb_base_mod
    use psb_prec_mod

  ! Arguments
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

  ! Local variables
    integer                            :: ictxt,me,np,err_act
    character(len=20)                  :: name

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


  
