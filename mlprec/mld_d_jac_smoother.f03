!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009, 2010
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
!
!
!
!
!
!
module mld_d_jac_smoother

  use mld_d_prec_type

  type, extends(mld_d_base_smoother_type) :: mld_d_jac_smoother_type
    ! The local solver component is inherited from the
    ! parent type. 
    !    class(mld_d_base_solver_type), allocatable :: sv
    !    
    type(psb_d_sparse_mat) :: nd
    integer                :: sweeps 
  contains
    procedure, pass(sm) :: build => d_jac_smoother_bld
    procedure, pass(sm) :: apply => d_jac_smoother_apply
    procedure, pass(sm) :: free  => d_jac_smoother_free
    procedure, pass(sm) :: seti  => d_jac_smoother_seti
    procedure, pass(sm) :: setc  => d_jac_smoother_setc
    procedure, pass(sm) :: setr  => d_jac_smoother_setr
    procedure, pass(sm) :: descr => d_jac_smoother_descr
    procedure, pass(sm) :: sizeof => d_jac_smoother_sizeof
  end type mld_d_jac_smoother_type


  private :: d_jac_smoother_bld, d_jac_smoother_apply, &
       &  d_jac_smoother_free,   d_jac_smoother_seti, &
       &  d_jac_smoother_setc,   d_jac_smoother_setr,&
       &  d_jac_smoother_descr,  d_jac_smoother_sizeof



contains

  subroutine d_jac_smoother_apply(alpha,sm,x,beta,y,desc_data,trans,work,info)
    use psb_sparse_mod
    type(psb_desc_type), intent(in)      :: desc_data
    class(mld_d_jac_smoother_type), intent(in) :: sm
    real(psb_dpk_),intent(in)            :: x(:)
    real(psb_dpk_),intent(inout)         :: y(:)
    real(psb_dpk_),intent(in)            :: alpha,beta
    character(len=1),intent(in)          :: trans
    real(psb_dpk_),target, intent(inout) :: work(:)
    integer, intent(out)                 :: info

    integer    :: n_row,n_col
    real(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
    integer    :: ictxt,np,me,i, err_act
    character          :: trans_
    character(len=20)  :: name='d_jac_smoother_apply'

    call psb_erractionsave(err_act)

    info = 0

    trans_ = psb_toupper(trans)
    select case(trans_)
    case('N')
    case('T','C')
    case default
      call psb_errpush(40,name)
      goto 9999
    end select

    if (.not.allocated(sm%sv)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    end if

    n_row = psb_cd_get_local_rows(desc_data)
    n_col = psb_cd_get_local_cols(desc_data)

    if (n_col <= size(work)) then 
      ww => work(1:n_col)
      if ((4*n_col+n_col) <= size(work)) then 
        aux => work(n_col+1:)
      else
        allocate(aux(4*n_col),stat=info)
        if (info /= 0) then 
          info=4025
          call psb_errpush(info,name,i_err=(/4*n_col,0,0,0,0/),&
               & a_err='real(psb_dpk_)')
          goto 9999      
        end if
      endif
    else
      allocate(ww(n_col),aux(4*n_col),stat=info)
      if (info /= 0) then 
        info=4025
        call psb_errpush(info,name,i_err=(/5*n_col,0,0,0,0/),&
             & a_err='real(psb_dpk_)')
        goto 9999      
      end if
    endif

    if (sm%sweeps == 1) then 

      call sm%sv%apply(alpha,x,beta,y,desc_data,trans_,aux,info) 

      if (info /= 0) then
        call psb_errpush(4001,name,a_err='Error in sub_aply Jacobi Sweeps = 1')
        goto 9999
      endif

    else if (sm%sweeps  > 1) then 

      !
      !
      ! Apply prec%iprcparm(mld_smoother_sweeps_) sweeps of a block-Jacobi solver
      ! to compute an approximate solution of a linear system.
      !
      !
      allocate(tx(n_col),ty(n_col),stat=info)
      if (info /= 0) then 
        info=4025
        call psb_errpush(info,name,i_err=(/2*n_col,0,0,0,0/),&
             & a_err='real(psb_dpk_)')
        goto 9999      
      end if

      tx = dzero
      ty = dzero
      do i=1, sm%sweeps
        !
        ! Compute Y(j+1) = D^(-1)*(X-ND*Y(j)), where D and ND are the
        ! block diagonal part and the remaining part of the local matrix
        ! and Y(j) is the approximate solution at sweep j.
        !
        ty(1:n_row) = x(1:n_row)
        call psb_spmm(-done,sm%nd,tx,done,ty,desc_data,info,work=aux,trans=trans_)

        if (info /=0) exit

        call sm%sv%apply(done,ty,dzero,tx,desc_data,trans_,aux,info) 

        if (info /=0) exit
      end do

      if (info == 0) call psb_geaxpby(alpha,tx,beta,y,desc_data,info)

      if (info /= 0) then 
        info=4001
        call psb_errpush(info,name,a_err='subsolve with Jacobi sweeps > 1')
        goto 9999      
      end if

      deallocate(tx,ty,stat=info)
      if (info /= 0) then 
        info=4001
        call psb_errpush(info,name,a_err='final cleanup with Jacobi sweeps > 1')
        goto 9999      
      end if

    else

      info = 10
      call psb_errpush(info,name,&
           & i_err=(/2,sm%sweeps,0,0,0/))
      goto 9999

    endif


    if (n_col <= size(work)) then 
      if ((4*n_col+n_col) <= size(work)) then 
      else
        deallocate(aux)
      endif
    else
      deallocate(ww,aux)
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine d_jac_smoother_apply

  subroutine d_jac_smoother_bld(a,desc_a,sm,upd,info)

    use psb_sparse_mod
    use mld_d_diag_solver
    Implicit None

    ! Arguments
    type(psb_d_sparse_mat), intent(in), target     :: a
    Type(psb_desc_type), Intent(in)                :: desc_a 
    class(mld_d_jac_smoother_type), intent(inout) :: sm
    character, intent(in)                          :: upd
    integer, intent(out)                           :: info
    ! Local variables
    integer :: n_row,n_col, nrow_a, nztota
    real(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
    integer :: ictxt,np,me,i, err_act, debug_unit, debug_level
    character(len=20)  :: name='d_jac_smoother_bld', ch_err
    
    info=0
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    ictxt       = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' start'


    n_row  = psb_cd_get_local_rows(desc_a)

    nrow_a = a%get_nrows()
    nztota = a%get_nzeros()
    select type (smsv => sm%sv)
    type is (mld_d_diag_solver_type)
      call a%clip_diag(sm%nd,info)
    class default
      call a%csclip(sm%nd,info,&
           & jmin=nrow_a+1,rscale=.false.,cscale=.false.)
    end select 
    if (info == 0) call sm%nd%cscnv(info,&
         & type='csr',dupl=psb_dupl_add_)
    if (info == 0) &
         & call sm%sv%build(a,desc_a,upd,info)
    if (info /= 0) then
      call psb_errpush(4010,name,a_err='clip & psb_spcnv csr 4')
      goto 9999
    end if

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' end'

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine d_jac_smoother_bld


  subroutine d_jac_smoother_seti(sm,what,val,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_d_jac_smoother_type), intent(inout) :: sm 
    integer, intent(in)                    :: what 
    integer, intent(in)                    :: val
    integer, intent(out)                   :: info
    Integer :: err_act
    character(len=20)  :: name='d_jac_smoother_seti'

    info = 0 
    call psb_erractionsave(err_act)

    select case(what) 
    case(mld_smoother_sweeps_) 
      sm%sweeps = val
    case default
      if (allocated(sm%sv)) then 
        call sm%sv%set(what,val,info)
      else
!!$        write(0,*) trim(name),' Missing component, not setting!'
!!$        info = 1121
      end if
    end select

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_jac_smoother_seti

  subroutine d_jac_smoother_setc(sm,what,val,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_d_jac_smoother_type), intent(inout) :: sm
    integer, intent(in)                    :: what 
    character(len=*), intent(in)           :: val
    integer, intent(out)                   :: info
    Integer :: err_act, ival
    character(len=20)  :: name='d_jac_smoother_setc'

    info = 0 
    call psb_erractionsave(err_act)


    call mld_stringval(val,ival,info)
    if (info == 0) call sm%set(what,ival,info)
    if (info /= 0) then
      info = 4010
      call psb_errpush(info, name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_jac_smoother_setc
  
  subroutine d_jac_smoother_setr(sm,what,val,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_d_jac_smoother_type), intent(inout) :: sm 
    integer, intent(in)                    :: what 
    real(psb_dpk_), intent(in)             :: val
    integer, intent(out)                   :: info
    Integer :: err_act
    character(len=20)  :: name='d_jac_smoother_setr'

    call psb_erractionsave(err_act)
    info = 0


    if (allocated(sm%sv)) then 
      call sm%sv%set(what,val,info)
    else
!!$      write(0,*) trim(name),' Missing component, not setting!'
!!$      info = 1121
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_jac_smoother_setr

  subroutine d_jac_smoother_free(sm,info)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_d_jac_smoother_type), intent(inout) :: sm
    integer, intent(out)                       :: info
    Integer :: err_act
    character(len=20)  :: name='d_jac_smoother_free'

    call psb_erractionsave(err_act)
    info = 0

    
    
    if (allocated(sm%sv)) then 
      call sm%sv%free(info)
      if (info == 0) deallocate(sm%sv,stat=info)
      if (info /= 0) then 
        info = 4000
        call psb_errpush(info,name)
        goto 9999 
      end if
    end if
    call sm%nd%free()

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_jac_smoother_free

  subroutine d_jac_smoother_descr(sm,info,iout)

    use psb_sparse_mod

    Implicit None

    ! Arguments
    class(mld_d_jac_smoother_type), intent(in) :: sm
    integer, intent(out)                       :: info
    integer, intent(in), optional              :: iout

    ! Local variables
    integer      :: err_act
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_d_jac_smoother_descr'
    integer :: iout_

    call psb_erractionsave(err_act)
    info = 0
    if (present(iout)) then 
      iout_ = iout 
    else
      iout_ = 6
    endif
    
    write(iout_,*) '  Block Jacobi smoother with  ',&
           &  sm%sweeps,' sweeps.'
    write(iout_,*) '  Local solver:'
    if (allocated(sm%sv)) then 
      call sm%sv%descr(info,iout_)
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_jac_smoother_descr

  function d_jac_smoother_sizeof(sm) result(val)
    use psb_sparse_mod
    implicit none 
    ! Arguments
    class(mld_d_jac_smoother_type), intent(in) :: sm
    integer(psb_long_int_k_) :: val
    integer             :: i

    val = psb_sizeof_int 
    if (allocated(sm%sv)) val = val + sm%sv%sizeof()
    val = val + sm%nd%sizeof()

    return
  end function d_jac_smoother_sizeof

end module mld_d_jac_smoother
