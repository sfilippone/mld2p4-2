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
! File: mld_das_aply.f90
!
! Subroutine: mld_das_aply
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
!   alpha      -  real(kind(0.d0)), input.
!                 The scalar alpha.
!   prec       -  type(mld_dbaseprc_type), input.
!                 The base preconditioner data structure containing the local part
!                 of the preconditioner K.
!   x          -  real(kind(0.d0)), dimension(:), input.
!                 The local part of the vector X.
!   beta       -  real(kind(0.d0)), input.
!                 The scalar beta.
!   y          -  real(kind(0.d0)), dimension(:), input/output.
!                 The local part of the vector Y.
!   desc_data  -  type(psb_desc_type), input.
!                 The communication descriptor associated to the matrix to be
!                 preconditioned.
!   trans      -  character, optional.
!                 If trans='N','n' then op(K^(-1)) = K^(-1);
!                 if trans='T','t' then op(K^(-1)) = K^(-T) (transpose of K^(-1)).
!   work       -  real(kind(0.d0)), dimension (:), optional, target.
!                 Workspace. Its size must be at least 4*psb_cd_get_local_cols(desc_data).
!   info       -  integer, output.
!                 Error code.
!  
subroutine mld_das_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_das_aply

  implicit none 

  ! Arguments
  type(psb_desc_type),intent(in)      :: desc_data
  type(mld_dbaseprc_type), intent(in) :: prec
  real(kind(0.d0)),intent(in)         :: x(:)
  real(kind(0.d0)),intent(inout)      :: y(:)
  real(kind(0.d0)),intent(in)         :: alpha,beta
  character(len=1)                    :: trans
  real(kind(0.d0)),target             :: work(:)
  integer, intent(out)                :: info

  ! Local variables
  integer :: n_row,n_col, int_err(5), nrow_d
  real(kind(1.d0)), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer           :: ictxt,np,me,isz, err_act
  character(len=20) :: name, ch_err

  name='mld_das_aply'
  info = 0
  call psb_erractionsave(err_act)

  ictxt = psb_cd_get_context(desc_data)

  call psb_info(ictxt, me, np)

  select case(trans)
  case('N','n')
  case('T','t','C','c')
  case default
    info=40
    int_err(1)=6
    ch_err(2:2)=trans
    goto 9999
  end select

  select case(prec%iprcparm(mld_prec_type_))

  case(mld_as_)
    !
    ! Additive Schwarz preconditioner
    !

    if (prec%iprcparm(mld_n_ovr_)==0) then
      ! 
      ! shortcut: this fixes performance for RAS(0) == BJA
      !
      call mld_bjac_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      if(info /= 0) then
        info=4010
        ch_err='psb_bjacaply'
        goto 9999
      end if

    else

      !
      ! Note: currently trans is unused
      !

      n_row  = psb_cd_get_local_rows(prec%desc_data)
      n_col  = psb_cd_get_local_cols(prec%desc_data)
      nrow_d = psb_cd_get_local_rows(desc_data)
      isz=max(n_row,N_COL)
      if ((6*isz) <= size(work)) then 
        ww => work(1:isz)
        tx => work(isz+1:2*isz)
        ty => work(2*isz+1:3*isz)
        aux => work(3*isz+1:)
      else if ((4*isz) <= size(work)) then 
        aux => work(1:)
        allocate(ww(isz),tx(isz),ty(isz),stat=info)
        if (info /= 0) then 
          call psb_errpush(4025,name,i_err=(/3*isz,0,0,0,0/),&
               & a_err='real(kind(1.d0))')
          goto 9999      
        end if
      else if ((3*isz) <= size(work)) then 
        ww => work(1:isz)
        tx => work(isz+1:2*isz)
        ty => work(2*isz+1:3*isz)
        allocate(aux(4*isz),stat=info)
        if (info /= 0) then 
          call psb_errpush(4025,name,i_err=(/4*isz,0,0,0,0/),&
               & a_err='real(kind(1.d0))')
          goto 9999      
        end if
      else 
        allocate(ww(isz),tx(isz),ty(isz),&
             &aux(4*isz),stat=info)
        if (info /= 0) then 
          call psb_errpush(4025,name,i_err=(/4*isz,0,0,0,0/),&
               & a_err='real(kind(1.d0))')
          goto 9999      
        end if

      endif

      tx(1:nrow_d)     = x(1:nrow_d) 
      tx(nrow_d+1:isz) = dzero

      !
      ! Get the overlap entries of tx (tx==x)
      ! 
      if (prec%iprcparm(mld_sub_restr_)==psb_halo_) then 
        call psb_halo(tx,prec%desc_data,info,work=aux,data=psb_comm_ext_)
        if(info /=0) then
          info=4010
          ch_err='psb_halo'
          goto 9999
        end if
      else if (prec%iprcparm(mld_sub_restr_) /= psb_none_) then 
        call psb_errpush(4001,name,a_err='Invalid mld_sub_restr_')
        goto 9999
      end if

      !
      ! If required, reorder tx according to the row/column permutation of the
      ! local extended matrix, stored into the permutation vector prec%perm
      !
      if (prec%iprcparm(mld_sub_ren_)>0) then 
        call psb_gelp('n',prec%perm,tx,info)
        if(info /=0) then
          info=4010
          ch_err='psb_gelp'
          goto 9999
        end if
      endif

      !
      ! Apply to tx the block-Jacobi preconditioner/solver (multiple sweeps of the
      ! block-Jacobi solver can be applied at the coarsest level of a multilevel
      ! preconditioner). The resulting vector is ty.
      !
      call mld_bjac_aply(done,prec,tx,dzero,ty,prec%desc_data,trans,aux,info)
      if(info /= 0) then
        info=4010
        ch_err='mld_bjac_aply'
        goto 9999
      end if

      !
      ! Apply to ty the inverse permutation of prec%perm
      !
      if (prec%iprcparm(mld_sub_ren_)>0) then 
        call psb_gelp('n',prec%invperm,ty,info)
        if(info /= 0) then
          info=4010
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
        if(info /=0) then
          info=4010
          ch_err='psb_ovrl'
          goto 9999
        end if

      case default
        call psb_errpush(4001,name,a_err='Invalid mld_sub_prol_')
        goto 9999
      end select

      !
      ! Compute y = beta*y + alpha*ty (ty==K^(-1)*tx)
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

    call psb_errpush(4001,name,a_err='Invalid mld_prec_type_')
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

end subroutine mld_das_aply

