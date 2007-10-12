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
subroutine mld_dsp_renum(a,desc_a,blck,p,atmp,info)
  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_dsp_renum

  implicit none

  !     .. array Arguments ..                                                     
  type(psb_dspmat_type), intent(in)      :: a,blck
  type(psb_dspmat_type), intent(inout)   :: atmp
  type(mld_dbaseprc_type), intent(inout) :: p
  type(psb_desc_type), intent(in)        :: desc_a
  integer, intent(out)   :: info


  character(len=20)      :: name, ch_err
  integer   nztota, nztotb, nztmp, nzl, nnr, ir, mglob, mtype, n_row, &
       & nrow_a,n_col, nhalo,lovr,  ind, iind, pi,nr,ns,i,j,jj,k,kk
  integer ::ictxt,np,me, err_act
  integer, allocatable  :: itmp(:), itmp2(:)
  real(kind(1.d0)), allocatable  :: rtmp(:)
  real(kind(1.d0)) :: t1,t2,t3,t4,t5,t6, t7, t8

  if (psb_get_errstatus().ne.0) return 
  info=0
  name='mld_dsp_renum'
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)

  !
  ! CHANGE: Start with a COO atmp. Then change if/when necessary. 
  ! Exit with a COO atmp. 
  !
  ! Renumbering type: 
  !     1. Global column indices
  !     (2. GPS band reduction disabled for the time being)
  nztota=psb_sp_get_nnzeros(a)
  nztotb=psb_sp_get_nnzeros(blck)
  call psb_spcnv(a,atmp,info,afmt='coo',dupl=psb_dupl_add_)
  call psb_rwextd(a%m+blck%m,atmp,info,blck)

  if (p%iprcparm(mld_sub_ren_)==mld_renum_glb_) then 

    ! This is the renumbering coherent with global indices..
    mglob = psb_cd_get_global_rows(desc_a)

    !
    !  Remember: we have switched IA1=COLS and IA2=ROWS
    !  Now identify the set of distinct local column indices
    !

    nnr = p%desc_data%matrix_data(psb_n_row_)
    allocate(p%perm(nnr),p%invperm(nnr),itmp2(nnr),stat=info)
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    do k=1,nnr
      itmp2(k) = p%desc_data%loc_to_glob(k)
    enddo
    !
    !  We want:  NEW(I) = OLD(PERM(I))
    ! 
    call psb_msort(itmp2(1:nnr),ix=p%perm)

    do k=1, nnr 
      p%invperm(p%perm(k)) = k
    enddo
    t3 = psb_wtime()

  else if (p%iprcparm(mld_sub_ren_)==mld_renum_gps_) then 
    
    call psb_spcnv(atmp,info,afmt='csr',dupl=psb_dupl_add_)
    nztmp = psb_sp_get_nnzeros(atmp)
    ! This is a renumbering with Gibbs-Poole-Stockmeyer 
    ! band reduction. Switched off for now. To be fixed,
    ! gps_reduction should get p%perm. 

    !          write(0,*) me,' Renumbering: realloc perms',atmp%m
    call psb_realloc(atmp%m,p%perm,info)
    if(info/=0) then
      info=4010
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_realloc(atmp%m,p%invperm,info)
    if(info/=0) then
      info=4010
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    allocate(itmp(max(8,atmp%m+2,nztmp+2)),stat=info)
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    itmp(1:8) = 0
    !          write(0,*) me,' Renumbering: Calling Metis'

    !          write(0,*) size(p%av(mld_u_pr_)%pl),size(p%av(mld_l_pr_)%pr)
    call  gps_reduction(atmp%m,atmp%ia2,atmp%ia1,p%perm,p%invperm,info)
    if(info/=0) then
      info=4010
      ch_err='gps_reduction'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    !      write(0,*) me,' Renumbering: Done GPS'
    !    call psb_barrier(ictxt)
    do i=1, atmp%m 
      if (p%perm(i) /= i) then 
        write(0,*) me,' permutation is not identity '
        exit
      endif
    enddo


    do k=1, nnr 
      p%invperm(p%perm(k)) = k
    enddo
    t3 = psb_wtime()
    
    call psb_spcnv(atmp,info,afmt='coo',dupl=psb_dupl_add_)

  end if

  ! Rebuild  ATMP with new numbering. 
  
  nztmp=psb_sp_get_nnzeros(atmp)
  do i=1,nztmp
    atmp%ia1(i) = p%perm(a%ia1(i))            
    atmp%ia2(i) = p%invperm(a%ia2(i))
  end do
  call psb_spcnv(atmp,info,afmt='coo',dupl=psb_dupl_add_)
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='psb_fixcoo')
    goto 9999      
  end if
  
  t4 = psb_wtime()
  
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

contains


  subroutine gps_reduction(m,ia,ja,perm,iperm,info)
    integer i,j,dgConn,Npnt,m
    integer n,idpth,ideg,ibw2,ipf2
    integer,dimension(:) :: perm,iperm,ia,ja
    integer, intent(out) :: info

    integer,dimension(:,:),allocatable::NDstk
    integer,dimension(:),allocatable::iOld,renum,ndeg,lvl,lvls1,lvls2,ccstor

    character(len=20)      :: name, ch_err

    if(psb_get_errstatus().ne.0) return 
    info=0
    name='gps_reduction'
    call psb_erractionsave(err_act)


    !--- Calcolo il massimo grado di connettivita'.
    npnt = m
    write(6,*) ' GPS su ',npnt
    dgConn=0
    do i=1,m
      dgconn = max(dgconn,(ia(i+1)-ia(i)))
    enddo
    !--- Il max valore di connettivita' e "dgConn"

    !--- Valori della common
    n=Npnt       !--- Numero di righe
    iDeg=dgConn  !--- Massima connettivita'
    !    iDpth=       !--- Numero di livelli non serve settarlo

    allocate(NDstk(Npnt,dgConn),stat=info)
    if (info/=0) then 
      info=4000
      call psb_errpush(info,name)
      goto 9999
    else
      write(0,*) 'gps_reduction first alloc OK'
    endif
    allocate(iOld(Npnt),renum(Npnt+1),ndeg(Npnt),lvl(Npnt),lvls1(Npnt),&
         &lvls2(Npnt),ccstor(Npnt),stat=info)
    if (info/=0) then 
      info=4000
      call psb_errpush(info,name)
      goto 9999
    else
      write(0,*) 'gps_reduction 2nd alloc OK'
    endif

    !--- Prepariamo il grafo della matrice
    Ndstk(:,:)=0
    do i=1,Npnt
      k=0
      do j = ia(i),ia(i+1) - 1 
        if ((1<=ja(j)).and.( ja( j ) /= i ).and.(ja(j)<=npnt)) then
          k = k+1
          Ndstk(i,k)=ja(j)
        endif
      enddo
      ndeg(i)=k
    enddo

    !--- Numerazione.
    do i=1,Npnt
      iOld(i)=i
    enddo
    write(0,*) 'gps_red : Preparation done'
    !--- 
    !--- Chiamiamo funzione reduce.
    call psb_gps_reduce(Ndstk,Npnt,iOld,renum,ndeg,lvl,lvls1, lvls2,ccstor,&
         & ibw2,ipf2,n,idpth,ideg)
    write(0,*) 'gps_red : Done reduce'
    !--- Permutazione
    perm(1:Npnt)=renum(1:Npnt)
    !--- Inversa permutazione
    do i=1,Npnt
      iperm(perm(i))=i
    enddo
    !--- Puliamo tutto.
    deallocate(NDstk,iOld,renum,ndeg,lvl,lvls1,lvls2,ccstor)

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine gps_reduction

end subroutine mld_dsp_renum
