!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010, 2010
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
!
!
!
!
!
!
module mld_c_as_smoother

  use mld_c_prec_type

  type, extends(mld_c_base_smoother_type) :: mld_c_as_smoother_type
    ! The local solver component is inherited from the
    ! parent type. 
    !    class(mld_c_base_solver_type), allocatable :: sv
    !    
    type(psb_cspmat_type) :: nd
    type(psb_desc_type)    :: desc_data 
    integer                :: novr, restr, prol, nd_nnz_tot
  contains
    procedure, pass(sm) :: check => c_as_smoother_check
    procedure, pass(sm) :: dump  => c_as_smoother_dmp
    procedure, pass(sm) :: build => c_as_smoother_bld
    procedure, pass(sm) :: apply => c_as_smoother_apply
    procedure, pass(sm) :: free  => c_as_smoother_free
    procedure, pass(sm) :: seti  => c_as_smoother_seti
    procedure, pass(sm) :: setc  => c_as_smoother_setc
    procedure, pass(sm) :: setr  => c_as_smoother_setr
    procedure, pass(sm) :: descr => c_as_smoother_descr
    procedure, pass(sm) :: sizeof => c_as_smoother_sizeof
    procedure, pass(sm) :: default => c_as_smoother_default
  end type mld_c_as_smoother_type
  
  
  private :: c_as_smoother_bld, c_as_smoother_apply, &
       &  c_as_smoother_free,   c_as_smoother_seti, &
       &  c_as_smoother_setc,   c_as_smoother_setr,&
       &  c_as_smoother_descr,  c_as_smoother_sizeof, &
       &  c_as_smoother_check,  c_as_smoother_default,&
       &  c_as_smoother_dmp
  
  character(len=6), parameter, private :: &
       &  restrict_names(0:4)=(/'none ','halo ','     ','     ','     '/)
  character(len=12), parameter, private :: &
       &  prolong_names(0:3)=(/'none       ','sum        ','average    ','square root'/)


contains

  subroutine c_as_smoother_default(sm)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_c_as_smoother_type), intent(inout) :: sm 

    
    sm%restr = psb_halo_
    sm%prol  = psb_none_
    sm%novr  = 1

    
    if (allocated(sm%sv)) then 
      call sm%sv%default()
    end if

    return
  end subroutine c_as_smoother_default

  subroutine c_as_smoother_check(sm,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_c_as_smoother_type), intent(inout) :: sm 
    integer, intent(out)                   :: info
    Integer           :: err_act
    character(len=20) :: name='c_as_smoother_check'

    call psb_erractionsave(err_act)
    info = psb_success_

    call mld_check_def(sm%restr,&
         & 'Restrictor',psb_halo_,is_legal_restrict)
    call mld_check_def(sm%prol,&
         & 'Prolongator',psb_none_,is_legal_prolong)
    call mld_check_def(sm%novr,&
         & 'Overlap layers ',0,is_legal_n_ovr)

    
    if (allocated(sm%sv)) then 
      call sm%sv%check(info)
    else 
      info=3111
      call psb_errpush(info,name)
      goto 9999
    end if

    if (info /= psb_success_) goto 9999
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine c_as_smoother_check

  subroutine c_as_smoother_apply(alpha,sm,x,beta,y,desc_data,trans,sweeps,work,info)
    use psb_base_mod
    type(psb_desc_type), intent(in)      :: desc_data
    class(mld_c_as_smoother_type), intent(in) :: sm
    complex(psb_spk_),intent(inout)      :: x(:)
    complex(psb_spk_),intent(inout)      :: y(:)
    complex(psb_spk_),intent(in)         :: alpha,beta
    character(len=1),intent(in)          :: trans
    integer, intent(in)                  :: sweeps
    complex(psb_spk_),target, intent(inout) :: work(:)
    integer, intent(out)                 :: info

    integer    :: n_row,n_col, nrow_d, i
    complex(psb_spk_), pointer :: ww(:), aux(:), tx(:),ty(:)
    integer    :: ictxt,np,me, err_act,isz,int_err(5)
    character          :: trans_
    character(len=20)  :: name='c_as_smoother_apply', ch_err

    call psb_erractionsave(err_act)

    info = psb_success_
    ictxt = desc_data%get_context()
    call psb_info (ictxt,me,np)

    trans_ = psb_toupper(trans)
    select case(trans_)
    case('N')
    case('T','C')
    case default
      call psb_errpush(psb_err_iarg_invalid_i_,name)
      goto 9999
    end select

    if (.not.allocated(sm%sv)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    end if


    n_row = sm%desc_data%get_local_rows()
    n_col = sm%desc_data%get_local_cols()
    nrow_d = desc_data%get_local_rows()
    isz=max(n_row,N_COL)
    if ((6*isz) <= size(work)) then 
      ww => work(1:isz)
      tx => work(isz+1:2*isz)
      ty => work(2*isz+1:3*isz)
      aux => work(3*isz+1:)
    else if ((4*isz) <= size(work)) then 
      aux => work(1:)
      allocate(ww(isz),tx(isz),ty(isz),stat=info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_alloc_request_,name,i_err=(/3*isz,0,0,0,0/),&
             & a_err='complex(psb_spk_)')
        goto 9999      
      end if
    else if ((3*isz) <= size(work)) then 
      ww => work(1:isz)
      tx => work(isz+1:2*isz)
      ty => work(2*isz+1:3*isz)
      allocate(aux(4*isz),stat=info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_alloc_request_,name,i_err=(/4*isz,0,0,0,0/),&
             & a_err='complex(psb_spk_)')
        goto 9999      
      end if
    else 
      allocate(ww(isz),tx(isz),ty(isz),&
           &aux(4*isz),stat=info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_alloc_request_,name,i_err=(/4*isz,0,0,0,0/),&
             & a_err='complex(psb_spk_)')
        goto 9999      
      end if

    endif

    if ((sm%novr == 0).and.(sweeps == 1)) then 
      !
      ! Shortcut: in this case it's just the same
      ! as Block Jacobi.
      !
      call sm%sv%apply(alpha,x,beta,y,desc_data,trans_,aux,info) 

      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error in sub_aply Jacobi Sweeps = 1')
        goto 9999
      endif

    else 


      tx(1:nrow_d)     = x(1:nrow_d) 
      tx(nrow_d+1:isz) = czero


      if (sweeps == 1) then 

        select case(trans_)
        case('N')
          !
          ! Get the overlap entries of tx (tx == x)
          ! 
          if (sm%restr == psb_halo_) then 
            call psb_halo(tx,sm%desc_data,info,work=aux,data=psb_comm_ext_)
            if(info /= psb_success_) then
              info=psb_err_from_subroutine_
              ch_err='psb_halo'
              goto 9999
            end if
          else if (sm%restr /= psb_none_) then 
            call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_restr_')
            goto 9999
          end if


        case('T','C')
          !
          ! With transpose, we have to do it here
          ! 

          select case (sm%prol) 

          case(psb_none_)
            ! 
            ! Do nothing

          case(psb_sum_) 
            !
            ! The transpose of sum is halo
            !
            call psb_halo(tx,sm%desc_data,info,work=aux,data=psb_comm_ext_)
            if(info /= psb_success_) then
              info=psb_err_from_subroutine_
              ch_err='psb_halo'
              goto 9999
            end if

          case(psb_avg_) 
            !
            ! Tricky one: first we have to scale the overlap entries,
            ! which we can do by assignind mode=0, i.e. no communication
            ! (hence only scaling), then we do the halo
            !
            call psb_ovrl(tx,sm%desc_data,info,&
                 & update=psb_avg_,work=aux,mode=0)
            if(info /= psb_success_) then
              info=psb_err_from_subroutine_
              ch_err='psb_ovrl'
              goto 9999
            end if
            call psb_halo(tx,sm%desc_data,info,work=aux,data=psb_comm_ext_)
            if(info /= psb_success_) then
              info=psb_err_from_subroutine_
              ch_err='psb_halo'
              goto 9999
            end if

          case default
            call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_prol_')
            goto 9999
          end select


        case default
          info=psb_err_iarg_invalid_i_
          int_err(1)=6
          ch_err(2:2)=trans
          goto 9999
        end select


        call sm%sv%apply(cone,tx,czero,ty,sm%desc_data,trans_,aux,info) 

        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='Error in sub_aply Jacobi Sweeps = 1')
          goto 9999
        endif

        select case(trans_)
        case('N')

          select case (sm%prol) 

          case(psb_none_)
            ! 
            ! Would work anyway, but since it is supposed to do nothing ...
            !        call psb_ovrl(ty,sm%desc_data,info,&
            !             & update=sm%prol,work=aux)


          case(psb_sum_,psb_avg_) 
            !
            ! Update the overlap of ty
            !
            call psb_ovrl(ty,sm%desc_data,info,&
                 & update=sm%prol,work=aux)
            if(info /= psb_success_) then
              info=psb_err_from_subroutine_
              ch_err='psb_ovrl'
              goto 9999
            end if

          case default
            call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_prol_')
            goto 9999
          end select

        case('T','C')
          !
          ! With transpose, we have to do it here
          ! 
          if (sm%restr == psb_halo_) then 
            call psb_ovrl(ty,sm%desc_data,info,&
                 & update=psb_sum_,work=aux)
            if(info /= psb_success_) then
              info=psb_err_from_subroutine_
              ch_err='psb_ovrl'
              goto 9999
            end if
          else if (sm%restr /= psb_none_) then 
            call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_restr_')
            goto 9999
          end if

        case default
          info=psb_err_iarg_invalid_i_
          int_err(1)=6
          ch_err(2:2)=trans
          goto 9999
        end select



      else if (sweeps  > 1) then 

        !
        !
        ! Apply prec%iprcparm(mld_smoother_sweeps_) sweeps of a block-Jacobi solver
        ! to compute an approximate solution of a linear system.
        !
        !
        ty = czero
        do i=1, sweeps
          select case(trans_)
          case('N')
            !
            ! Get the overlap entries of tx (tx == x)
            ! 
            if (sm%restr == psb_halo_) then 
              call psb_halo(tx,sm%desc_data,info,work=aux,data=psb_comm_ext_)
              if(info /= psb_success_) then
                info=psb_err_from_subroutine_
                ch_err='psb_halo'
                goto 9999
              end if
            else if (sm%restr /= psb_none_) then 
              call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_restr_')
              goto 9999
            end if


          case('T','C')
            !
            ! With transpose, we have to do it here
            ! 

            select case (sm%prol) 

            case(psb_none_)
              ! 
              ! Do nothing

            case(psb_sum_) 
              !
              ! The transpose of sum is halo
              !
              call psb_halo(tx,sm%desc_data,info,work=aux,data=psb_comm_ext_)
              if(info /= psb_success_) then
                info=psb_err_from_subroutine_
                ch_err='psb_halo'
                goto 9999
              end if

            case(psb_avg_) 
              !
              ! Tricky one: first we have to scale the overlap entries,
              ! which we can do by assignind mode=0, i.e. no communication
              ! (hence only scaling), then we do the halo
              !
              call psb_ovrl(tx,sm%desc_data,info,&
                   & update=psb_avg_,work=aux,mode=0)
              if(info /= psb_success_) then
                info=psb_err_from_subroutine_
                ch_err='psb_ovrl'
                goto 9999
              end if
              call psb_halo(tx,sm%desc_data,info,work=aux,data=psb_comm_ext_)
              if(info /= psb_success_) then
                info=psb_err_from_subroutine_
                ch_err='psb_halo'
                goto 9999
              end if

            case default
              call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_prol_')
              goto 9999
            end select


          case default
            info=psb_err_iarg_invalid_i_
            int_err(1)=6
            ch_err(2:2)=trans
            goto 9999
          end select
          !
          ! Compute Y(j+1) = D^(-1)*(X-ND*Y(j)), where D and ND are the
          ! block diagonal part and the remaining part of the local matrix
          ! and Y(j) is the approximate solution at sweep j.
          !
          ww(1:n_row) = tx(1:n_row)
          call psb_spmm(-cone,sm%nd,ty,cone,ww,sm%desc_data,info,work=aux,trans=trans_)

          if (info /= psb_success_) exit

          call sm%sv%apply(cone,ww,czero,ty,sm%desc_data,trans_,aux,info) 

          if (info /= psb_success_) exit


          select case(trans_)
          case('N')

            select case (sm%prol) 

            case(psb_none_)
              ! 
              ! Would work anyway, but since it is supposed to do nothing ...
              !        call psb_ovrl(ty,sm%desc_data,info,&
              !             & update=sm%prol,work=aux)


            case(psb_sum_,psb_avg_) 
              !
              ! Update the overlap of ty
              !
              call psb_ovrl(ty,sm%desc_data,info,&
                   & update=sm%prol,work=aux)
              if(info /= psb_success_) then
                info=psb_err_from_subroutine_
                ch_err='psb_ovrl'
                goto 9999
              end if

            case default
              call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_prol_')
              goto 9999
            end select

          case('T','C')
            !
            ! With transpose, we have to do it here
            ! 
            if (sm%restr == psb_halo_) then 
              call psb_ovrl(ty,sm%desc_data,info,&
                   & update=psb_sum_,work=aux)
              if(info /= psb_success_) then
                info=psb_err_from_subroutine_
                ch_err='psb_ovrl'
                goto 9999
              end if
            else if (sm%restr /= psb_none_) then 
              call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_restr_')
              goto 9999
            end if

          case default
            info=psb_err_iarg_invalid_i_
            int_err(1)=6
            ch_err(2:2)=trans
            goto 9999
          end select
        end do

        if (info /= psb_success_) then 
          info=psb_err_internal_error_
          call psb_errpush(info,name,a_err='subsolve with Jacobi sweeps > 1')
          goto 9999      
        end if


      else

        info = psb_err_iarg_neg_
        call psb_errpush(info,name,&
             & i_err=(/2,sweeps,0,0,0/))
        goto 9999


      end if

      !
      ! Compute y = beta*y + alpha*ty (ty == K^(-1)*tx)
      !
      call psb_geaxpby(alpha,ty,beta,y,desc_data,info) 

    end if


    if ((6*isz) <= size(work)) then 
    else if ((4*isz) <= size(work)) then 
      deallocate(ww,tx,ty)
    else if ((3*isz) <= size(work)) then 
      deallocate(aux)
    else 
      deallocate(ww,aux,tx,ty)
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

  end subroutine c_as_smoother_apply

  subroutine c_as_smoother_bld(a,desc_a,sm,upd,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    type(psb_cspmat_type), intent(in), target   :: a
    Type(psb_desc_type), Intent(in)              :: desc_a 
    class(mld_c_as_smoother_type), intent(inout) :: sm
    character, intent(in)                        :: upd
    integer, intent(out)                         :: info
    ! Local variables
    type(psb_cspmat_type) :: blck, atmp
    integer :: n_row,n_col, nrow_a, nhalo, novr, data_, nzeros
    complex(psb_spk_), pointer :: ww(:), aux(:), tx(:),ty(:)
    integer :: ictxt,np,me,i, err_act, debug_unit, debug_level
    character(len=20)  :: name='c_as_smoother_bld', ch_err

    info=psb_success_
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    ictxt       = desc_a%get_context()
    call psb_info(ictxt, me, np)
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' start'


    novr   = sm%novr
    if (novr < 0) then
      info=psb_err_invalid_ovr_num_
      call psb_errpush(info,name,i_err=(/novr,0,0,0,0,0/))
      goto 9999
    endif

    if ((novr == 0).or.(np == 1)) then 
      if (psb_toupper(upd) == 'F') then 
        call psb_cdcpy(desc_a,sm%desc_data,info)
        If(debug_level >= psb_debug_outer_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & '  done cdcpy'
        if(info /= psb_success_) then
          info=psb_err_from_subroutine_
          ch_err='psb_cdcpy'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        end if
        if (debug_level >= psb_debug_outer_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & 'Early return: P>=3 N_OVR=0'
      endif
      call blck%csall(0,0,info,1)
    else

      If (psb_toupper(upd) == 'F') Then
        !
        ! Build the auxiliary descriptor desc_p%matrix_data(psb_n_row_).
        ! This is done by psb_cdbldext (interface to psb_cdovr), which is
        ! independent of CSR, and has been placed in the tools directory
        ! of PSBLAS, instead of the mlprec directory of MLD2P4, because it
        ! might be used independently of the AS preconditioner, to build
        ! a descriptor for an extended stencil in a PDE solver. 
        !
        call psb_cdbldext(a,desc_a,novr,sm%desc_data,info,extype=psb_ovt_asov_)
        if(debug_level >= psb_debug_outer_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & ' From cdbldext _:',sm%desc_data%get_local_rows(),&
             & sm%desc_data%get_local_cols()
        
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_
          ch_err='psb_cdbldext'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        end if
      Endif

      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & 'Before sphalo '

      !
      ! Retrieve the remote sparse matrix rows required for the AS extended
      ! matrix
      data_ = psb_comm_ext_
      Call psb_sphalo(a,sm%desc_data,blck,info,data=data_,rowscale=.true.)
      
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='psb_sphalo'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
      
      if (debug_level >=psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & 'After psb_sphalo ',&
           & blck%get_nrows(), blck%get_nzeros()

    End if
    if (info == psb_success_) &
         & call sm%sv%build(a,sm%desc_data,upd,info,blck)

    nrow_a = a%get_nrows()
    n_row  = sm%desc_data%get_local_rows()
    n_col  = sm%desc_data%get_local_cols()
    
    if (info == psb_success_) call a%csclip(sm%nd,info,&
         & jmin=nrow_a+1,rscale=.false.,cscale=.false.)
    if (info == psb_success_) call blck%csclip(atmp,info,&
         & jmin=nrow_a+1,rscale=.false.,cscale=.false.)
    if (info == psb_success_) call psb_rwextd(n_row,sm%nd,info,b=atmp) 
    if (info == psb_success_) call sm%nd%cscnv(info,&
         & type='csr',dupl=psb_dupl_add_)

    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='clip & psb_spcnv csr 4')
      goto 9999
    end if
    nzeros = sm%nd%get_nzeros()
!!$    write(0,*) me,' ND nzeors ',nzeros
    call psb_sum(ictxt,nzeros)
    sm%nd_nnz_tot = nzeros

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
  end subroutine c_as_smoother_bld


  subroutine c_as_smoother_seti(sm,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_c_as_smoother_type), intent(inout) :: sm 
    integer, intent(in)                    :: what 
    integer, intent(in)                    :: val
    integer, intent(out)                   :: info
    Integer :: err_act
    character(len=20)  :: name='c_as_smoother_seti'

    info = psb_success_
    call psb_erractionsave(err_act)

    select case(what) 
!!$    case(mld_smoother_sweeps_) 
!!$      sm%sweeps = val
    case(mld_sub_ovr_) 
      sm%novr   = val
    case(mld_sub_restr_) 
      sm%restr  = val
    case(mld_sub_prol_) 
      sm%prol   = val
    case default
      if (allocated(sm%sv)) then 
        call sm%sv%set(what,val,info)
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
  end subroutine c_as_smoother_seti

  subroutine c_as_smoother_setc(sm,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_c_as_smoother_type), intent(inout) :: sm
    integer, intent(in)                    :: what 
    character(len=*), intent(in)           :: val
    integer, intent(out)                   :: info
    Integer :: err_act, ival
    character(len=20)  :: name='c_as_smoother_setc'

    info = psb_success_
    call psb_erractionsave(err_act)


    call mld_stringval(val,ival,info)
    if (info == psb_success_) call sm%set(what,ival,info)
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
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
  end subroutine c_as_smoother_setc

  subroutine c_as_smoother_setr(sm,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_c_as_smoother_type), intent(inout) :: sm 
    integer, intent(in)                    :: what 
    real(psb_spk_), intent(in)             :: val
    integer, intent(out)                   :: info
    Integer :: err_act
    character(len=20)  :: name='c_as_smoother_setr'

    call psb_erractionsave(err_act)
    info = psb_success_


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
  end subroutine c_as_smoother_setr

  subroutine c_as_smoother_free(sm,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_c_as_smoother_type), intent(inout) :: sm
    integer, intent(out)                       :: info
    Integer :: err_act
    character(len=20)  :: name='c_as_smoother_free'

    call psb_erractionsave(err_act)
    info = psb_success_



    if (allocated(sm%sv)) then 
      call sm%sv%free(info)
      if (info == psb_success_) deallocate(sm%sv,stat=info)
      if (info /= psb_success_) then 
        info = psb_err_alloc_dealloc_
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
  end subroutine c_as_smoother_free

  subroutine c_as_smoother_descr(sm,info,iout,coarse)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_c_as_smoother_type), intent(in) :: sm
    integer, intent(out)                      :: info
    integer, intent(in), optional             :: iout
    logical, intent(in), optional             :: coarse

    ! Local variables
    integer      :: err_act
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_c_as_smoother_descr'
    integer :: iout_
    logical      :: coarse_

    call psb_erractionsave(err_act)
    info = psb_success_
    if (present(coarse)) then 
      coarse_ = coarse
    else
      coarse_ = .false.
    end if
    if (present(iout)) then 
      iout_ = iout 
    else
      iout_ = 6
    endif

    if (.not.coarse_) then 
      write(iout_,*) '  Additive Schwarz with  ',&
           &  sm%novr, ' overlap layers.'
      write(iout_,*) '  Restrictor:  ',restrict_names(sm%restr)
      write(iout_,*) '  Prolongator: ',prolong_names(sm%prol)
      write(iout_,*) '  Local solver:'
    endif
    if (allocated(sm%sv)) then 
      call sm%sv%descr(info,iout_,coarse=coarse)
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
  end subroutine c_as_smoother_descr

  function c_as_smoother_sizeof(sm) result(val)
    use psb_base_mod
    implicit none 
    ! Arguments
    class(mld_c_as_smoother_type), intent(in) :: sm
    integer(psb_long_int_k_) :: val
    integer             :: i

    val = psb_sizeof_int 
    if (allocated(sm%sv)) val = val + sm%sv%sizeof()
    val = val + sm%nd%sizeof()

    return
  end function c_as_smoother_sizeof

  subroutine c_as_smoother_dmp(sm,ictxt,level,info,prefix,head,smoother,solver)
    use psb_base_mod
    implicit none 
    class(mld_c_as_smoother_type), intent(in) :: sm
    integer, intent(in)              :: ictxt,level
    integer, intent(out)             :: info
    character(len=*), intent(in), optional :: prefix, head
    logical, optional, intent(in)    :: smoother, solver
    integer :: i, j, il1, iln, lname, lev
    integer :: icontxt,iam, np
    character(len=80)  :: prefix_
    character(len=120) :: fname ! len should be at least 20 more than
    logical :: smoother_
    !  len of prefix_ 

    info = 0

    if (present(prefix)) then 
      prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
    else
      prefix_ = "dump_smth_c"
    end if

    call psb_info(ictxt,iam,np)

    if (present(smoother)) then 
      smoother_ = smoother
    else
      smoother_ = .false. 
    end if
    lname = len_trim(prefix_)
    fname = trim(prefix_)
    write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
    lname = lname + 5

    if (smoother_) then 
      write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_nd.mtx'
      if (sm%nd%is_asb()) &
           & call sm%nd%print(fname,head=head)
    end if
    ! At base level do nothing for the smoother
    if (allocated(sm%sv)) &
         & call sm%sv%dump(ictxt,level,info,solver=solver)

  end subroutine c_as_smoother_dmp

end module mld_c_as_smoother
