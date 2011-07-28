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
! File: mld_zas_aply.f90
!
! Subroutine: mld_zas_aply
! Version:    real
!
!  This routine applies the Additive Schwarz preconditioner by computing
!
!                          Y = beta*Y + alpha*op(K^(-1))*X,
!  where
!  - K is the base preconditioner, stored in prec,
!  - op(K^(-1)) is K^(-1) or its transpose, according to the value of trans,
!  - X and Y are vectors,
!  - alpha and beta are scalars.
!
!
! Arguments:
!   alpha      -  real(psb_dpk_), input.
!                 The scalar alpha.
!   prec       -  type(mld_dbaseprec_type), input.
!                 The base preconditioner data structure containing the local part
!                 of the preconditioner K.
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
!                 If trans='N','n' then op(K^(-1)) = K^(-1);
!                 if trans='T','t' then op(K^(-1)) = K^(-T) (transpose of K^(-1)).
!   work       -  real(psb_dpk_), dimension (:), optional, target.
!                 Workspace. Its size must be at least 4*desc_data%get_local_cols().
!   info       -  integer, output.
!                 Error code.
!  
subroutine mld_zas_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)

  use psb_base_mod
  use mld_z_inner_mod, mld_protect_name => mld_zas_aply

  implicit none 

  ! Arguments
  type(psb_desc_type),intent(in)      :: desc_data
  type(mld_zbaseprec_type), intent(in) :: prec
  complex(psb_dpk_),intent(in)         :: x(:)
  complex(psb_dpk_),intent(inout)      :: y(:)
  complex(psb_dpk_),intent(in)         :: alpha,beta
  character(len=1)                    :: trans
  complex(psb_dpk_),target             :: work(:)
  integer, intent(out)                :: info

  ! Local variables
  integer :: n_row,n_col, int_err(5), nrow_d
  complex(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer           :: ictxt,np,me,isz, err_act
  character(len=20) :: name, ch_err
  character         :: trans_

  name='mld_zas_aply'
  info = psb_success_
  call psb_erractionsave(err_act)

  ictxt = desc_data%get_context()

  call psb_info(ictxt, me, np)
  
  trans_ = psb_toupper(trans)

  select case(prec%iprcparm(mld_smoother_type_))

  case(mld_bjac_)
    
    call mld_sub_aply(alpha,prec,x,beta,y,desc_data,trans_,work,info)
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='mld_sub_aply'
      goto 9999
    end if
    
  case(mld_as_)
    !
    ! Additive Schwarz preconditioner
    !

    if ((prec%iprcparm(mld_sub_ovr_) == 0).or.(np==1)) then
      ! 
      ! Shortcut: this fixes performance for RAS(0) == BJA
      !
      call mld_sub_aply(alpha,prec,x,beta,y,desc_data,trans_,work,info)
      if(info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='mld_sub_aply'
        goto 9999
      end if

    else
      !
      ! Overlap > 0
      !

      n_row  = prec%desc_data%get_local_rows()
      n_col  = prec%desc_data%get_local_cols()
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
               & a_err='complex(psb_dpk_)')
          goto 9999      
        end if
      else if ((3*isz) <= size(work)) then 
        ww => work(1:isz)
        tx => work(isz+1:2*isz)
        ty => work(2*isz+1:3*isz)
        allocate(aux(4*isz),stat=info)
        if (info /= psb_success_) then 
          call psb_errpush(psb_err_alloc_request_,name,i_err=(/4*isz,0,0,0,0/),&
               & a_err='complex(psb_dpk_)')
          goto 9999      
        end if
      else 
        allocate(ww(isz),tx(isz),ty(isz),&
             &aux(4*isz),stat=info)
        if (info /= psb_success_) then 
          call psb_errpush(psb_err_alloc_request_,name,i_err=(/4*isz,0,0,0,0/),&
               & a_err='complex(psb_dpk_)')
          goto 9999      
        end if

      endif

      tx(1:nrow_d)     = x(1:nrow_d) 
      tx(nrow_d+1:isz) = dzero

      select case(trans_)
      case('N')
        !
        ! Get the overlap entries of tx (tx == x)
        ! 
        if (prec%iprcparm(mld_sub_restr_) == psb_halo_) then 
          call psb_halo(tx,prec%desc_data,info,work=aux,data=psb_comm_ext_)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='psb_halo'
            goto 9999
          end if
        else if (prec%iprcparm(mld_sub_restr_) /= psb_none_) then 
          call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_restr_')
          goto 9999
        end if

        !
        ! If required, reorder tx according to the row/column permutation of the
        ! local extended matrix, stored into the permutation vector prec%perm
        !
        if (prec%iprcparm(mld_sub_ren_)>0) then 
          call psb_gelp('n',prec%perm,tx,info)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='psb_gelp'
            goto 9999
          end if
        endif

        !
        ! Apply to tx the block-Jacobi preconditioner/solver (multiple sweeps of the
        ! block-Jacobi solver can be applied at the coarsest level of a multilevel
        ! preconditioner). The resulting vector is ty.
        !
        call mld_sub_aply(zone,prec,tx,zzero,ty,prec%desc_data,trans_,aux,info)
        if(info /= psb_success_) then
          info=psb_err_from_subroutine_
          ch_err='mld_sub_aply'
          goto 9999
        end if

        !
        ! Apply to ty the inverse permutation of prec%perm
        !
        if (prec%iprcparm(mld_sub_ren_)>0) then 
          call psb_gelp('n',prec%invperm,ty,info)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='psb_gelp'
            goto 9999
          end if
        endif

        select case (prec%iprcparm(mld_sub_prol_)) 

        case(psb_none_)
          ! 
          ! Would work anyway, but since it is supposed to do nothing ...
          !        call psb_ovrl(ty,prec%desc_data,info,&
          !             & update=prec%iprcparm(mld_sub_prol_),work=aux)


        case(psb_sum_,psb_avg_) 
          !
          ! Update the overlap of ty
          !
          call psb_ovrl(ty,prec%desc_data,info,&
               & update=prec%iprcparm(mld_sub_prol_),work=aux)
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

        select case (prec%iprcparm(mld_sub_prol_)) 

        case(psb_none_)
          ! 
          ! Do nothing

        case(psb_sum_) 
          !
          ! The transpose of sum is halo
          !
          call psb_halo(tx,prec%desc_data,info,work=aux,data=psb_comm_ext_)
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
          call psb_ovrl(tx,prec%desc_data,info,&
               & update=psb_avg_,work=aux,mode=0)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='psb_ovrl'
            goto 9999
          end if
          call psb_halo(tx,prec%desc_data,info,work=aux,data=psb_comm_ext_)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='psb_halo'
            goto 9999
          end if

        case default
          call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_prol_')
          goto 9999
        end select

        !
        ! If required, reorder tx according to the row/column permutation of the
        ! local extended matrix, stored into the permutation vector prec%perm
        !
        if (prec%iprcparm(mld_sub_ren_)>0) then 
          call psb_gelp('n',prec%perm,tx,info)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='psb_gelp'
            goto 9999
          end if
        endif

        !
        ! Apply to tx the block-Jacobi preconditioner/solver (multiple sweeps of the
        ! block-Jacobi solver can be applied at the coarsest level of a multilevel
        ! preconditioner). The resulting vector is ty.
        !
        call mld_sub_aply(zone,prec,tx,zzero,ty,prec%desc_data,trans_,aux,info)
        if(info /= psb_success_) then
          info=psb_err_from_subroutine_
          ch_err='mld_sub_aply'
          goto 9999
        end if

        !
        ! Apply to ty the inverse permutation of prec%perm
        !
        if (prec%iprcparm(mld_sub_ren_)>0) then 
          call psb_gelp('n',prec%invperm,ty,info)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='psb_gelp'
            goto 9999
          end if
        endif

        !
        ! With transpose, we have to do it here
        ! 
        if (prec%iprcparm(mld_sub_restr_) == psb_halo_) then 
          call psb_ovrl(ty,prec%desc_data,info,&
               & update=psb_sum_,work=aux)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='psb_ovrl'
            goto 9999
          end if
        else if (prec%iprcparm(mld_sub_restr_) /= psb_none_) then 
          call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_sub_restr_')
          goto 9999
        end if

      case default
        info=psb_err_iarg_invalid_i_
        int_err(1)=6
        ch_err(2:2)=trans
        goto 9999
      end select

      !
      ! Compute y = beta*y + alpha*ty (ty == K^(-1)*tx)
      !
      call psb_geaxpby(alpha,ty,beta,y,desc_data,info) 


      if ((6*isz) <= size(work)) then 
      else if ((4*isz) <= size(work)) then 
        deallocate(ww,tx,ty)
      else if ((3*isz) <= size(work)) then 
        deallocate(aux)
      else 
        deallocate(ww,aux,tx,ty)
      endif
    end if

  case default

    call psb_errpush(psb_err_internal_error_,name,a_err='Invalid mld_smoother_type_')
    goto 9999

  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name,i_err=int_err,a_err=ch_err)
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_zas_aply

