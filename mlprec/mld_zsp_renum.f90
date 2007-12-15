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
! File: mld_zsp_renum.f90.
!
! Subroutine: mld_zsp_renum.
! Version:    real.
! Contains:   gps_reduction.
!
!  This routine reorders the rows and the columns of the local part of a sparse
!  distributed matrix, according to one of the following criteria:
!  1. the numbering of the global column indices,
!  2. the Gibbs-Poole-Stockmeyer (GPS) band reduction algorithm.
!  NOTE: the GPS algorithm is disabled for the time being (see mld_prec_type.f90).
!
!  The matrix to be reordered is stored into a and blck, as specified in the
!  description of the arguments below.
!
!  If required by the user (p%iprcparm(mld_sub_ren_) /= 0), the routine is
!  used by mld_bjac_bld in building the block-Jacobi and Additive Schwarz 
!  'base preconditioners' corresponding to any level of a multilevel
!  preconditioner.
!  
!
! Arguments:
!    a       -  type(psb_zspmat_type), input.
!               The sparse matrix structure containing the 'original' local
!               part of the matrix to be reordered, i.e. the rows of the matrix
!               held by the calling process according to the initial data
!               distribution.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor associated to a.
!    blck    -  type(psb_zspmat_type), input.
!               The sparse matrix structure containing the remote rows of the
!               matrix to be reordered, that have been retrieved by mld_asmat_bld
!               to build an Additive Schwarz base preconditioner with overlap
!               greater than 0.If the overlap is 0, then blck does not contain
!               any row.
!    p       -  type(mld_zbaseprc_type), input/output.
!               The base preconditioner data structure containing the local
!               part of the base preconditioner to be built. In input it 
!               contains information on the type of reordering to be applied
!               and on the matrix to be reordered. In output it contains
!               information on the reordering applied.
!    atmp    -  type(psb_zspmat_type), output.
!               The sparse matrix structure containing the whole local reordered 
!               matrix.
!    info    -  integer, output.                                                             
!               Error code.
! 
subroutine mld_zsp_renum(a,desc_a,blck,p,atmp,info)

  use psb_base_mod
  use mld_prec_mod, mld_protect_name => mld_zsp_renum

  implicit none

  ! Arguments                                                    
  type(psb_zspmat_type), intent(in)      :: a,blck
  type(psb_zspmat_type), intent(inout)   :: atmp
  type(mld_zbaseprc_type), intent(inout) :: p
  type(psb_desc_type), intent(in)        :: desc_a
  integer, intent(out)   :: info

  ! Local variables
  character(len=20)      :: name, ch_err
  integer   nztota, nztotb, nztmp, nnr, mglob, i,k
  integer ::ictxt,np,me, err_act
  integer, allocatable  :: itmp(:), itmp2(:)
  real(kind(1.d0)) :: t3,t4

  if (psb_get_errstatus().ne.0) return 
  info=0
  name='mld_zsp_renum'
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)

  !
  ! NOTE: the matrix to be reordered is converted into the COO format.
  !       If necessary it is converted from the COO to the CSR format.
  !       The output matrix is in COO format.
  !

  !
  ! Convert a into the COO format and extend it up to a%m+blck%m rows
  ! by adding null rows. The converted extended matrix is stored in atmp.
  !
  nztota=psb_sp_get_nnzeros(a)
  nztotb=psb_sp_get_nnzeros(blck)
  call psb_spcnv(a,atmp,info,afmt='coo',dupl=psb_dupl_add_)
  call psb_rwextd(a%m+blck%m,atmp,info,blck)

  if (p%iprcparm(mld_sub_ren_)==mld_renum_glb_) then

    !
    ! Compute the row and column reordering coherent with the
    ! global indices
    !

    mglob = psb_cd_get_global_rows(desc_a)
    !
    !  Remember: we have switched IA1=COLS and IA2=ROWS.
    !  Now identify the set of distinct local column indices.
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
    ! Compute reordering. We want new(i) = old(perm(i)).
    ! 
    call psb_msort(itmp2(1:nnr),ix=p%perm)
    !
    ! Compute the inverse of the permutation stored in perm
    !
    do k=1, nnr 
      p%invperm(p%perm(k)) = k
    enddo
    t3 = psb_wtime()

  else if (p%iprcparm(mld_sub_ren_)==mld_renum_gps_) then

    !
    ! This is a renumbering with Gibbs-Poole-Stockmeyer 
    ! band reduction. Switched off for now. To be fixed,
    ! gps_reduction should get p%perm.
    !

    !
    ! Convert atmp into the CSR format
    !
    call psb_spcnv(atmp,info,afmt='csr',dupl=psb_dupl_add_)
    nztmp = psb_sp_get_nnzeros(atmp)

    !
    ! Realloc the permutation arrays
    !
    ! write(0,*) me,' Renumbering: realloc perms',atmp%m
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

    ! write(0,*) me,' Renumbering: Calling Metis'
    ! write(0,*) size(p%av(mld_u_pr_)%pl),size(p%av(mld_l_pr_)%pr)

    !
    ! Renumber rows and columns according to the GPS algorithm
    !
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

    !
    ! Compute the inverse permutation
    !
    do k=1, atmp%m
      p%invperm(p%perm(k)) = k
    enddo
    t3 = psb_wtime()

    call psb_spcnv(atmp,info,afmt='coo',dupl=psb_dupl_add_)

  end if

  !
  ! Rebuild atmp with the new numbering      (COO format)
  !
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

  !
  ! Subroutine: gps_reduction.
  ! Note: internal subroutine of mld_zsp_renum.
  !
  !  Compute a renumbering of the row and column indices of a sparse matrix
  !  according to the  Gibbs-Poole-Stockmeyer band reduction algorithm. The
  !  matrix is stored in CSR format.
  !
  !  This routine has been obtained by adapting ACM TOMS Algorithm 582.  
  !  
  !
  ! Arguments:
  !    m       -  integer, ...
  !                        The number of rows of the matrix to which the renumbering
  !               is applied.
  !    ia      -  integer, dimension(:), ...
  !               The indices identifying the first nonzero entry of each row
  !               of the matrix, according to the CSR storage format.
  !    ja      -  integer, dimension(:), ...
  !               The column indices of the nonzero entries of the matrix,
  !               according to the CSR storage format.
  !    perm    -  integer, dimension(:), ...
  !               The row/column index permutation corresponding to the
  !               renumbering.
  !    iperm   -  integer, dimension(:),...
  !               The inverse of the row/column permutation stored in perm.
  !    info    -  integer, output.
  !               Error code
  !
  subroutine gps_reduction(m,ia,ja,perm,iperm,info)

    integer i,j,dgConn,Npnt,m
    integer n,idpth,ideg,ibw2,ipf2
    integer,dimension(:) :: perm,iperm,ia,ja
    integer, intent(out) :: info

    integer,dimension(:,:),allocatable::NDstk
    integer,dimension(:),allocatable::iOld,renum,ndeg,lvl,lvls1,lvls2,ccstor

    character(len=20)      :: name

    if(psb_get_errstatus().ne.0) return 
    info=0
    name='gps_reduction'
    call psb_erractionsave(err_act)

    ! Compute the maximum connectivity degree
    npnt = m
    write(6,*) ' GPS su ',npnt
    dgConn=0
    do i=1,m
      dgconn = max(dgconn,(ia(i+1)-ia(i)))
    enddo
    ! The maximum connectivity value is dgConn

    n=Npnt       ! Max number of rows
    iDeg=dgConn  ! Max connectivity
    ! iDpth=     ! Number of level (initialization not needed)

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

    ! Prepare the matrix graph
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

    ! Numbering
    do i=1,Npnt
      iOld(i)=i
    enddo
    ! write(0,*) 'gps_red : Preparation done'

    ! Call gps_reduce
    call psb_gps_reduce(Ndstk,Npnt,iOld,renum,ndeg,lvl,lvls1, lvls2,ccstor,&
         & ibw2,ipf2,n,idpth,ideg)
    write(0,*) 'gps_red : Done reduce'

    ! Build permutation vector
    perm(1:Npnt)=renum(1:Npnt)

    !Build inverse permutation vector
    do i=1,Npnt
      iperm(perm(i))=i
    enddo

    ! Deallocate memory
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

end subroutine mld_zsp_renum
