module mld_d_ilu_fact_mod

  use mld_base_prec_type

  interface mld_ilu0_fact
    subroutine mld_dilu0_fact(ialg,a,l,u,d,info,blck,upd)
      import psb_dspmat_type, psb_dpk_
      integer, intent(in)                 :: ialg
      integer, intent(out)                :: info
      type(psb_dspmat_type),intent(in)    :: a
      type(psb_dspmat_type),intent(inout) :: l,u
      type(psb_dspmat_type),intent(in), optional, target :: blck
      character, intent(in), optional      :: upd
      real(psb_dpk_), intent(inout)     ::  d(:)
    end subroutine mld_dilu0_fact
  end interface

  interface mld_iluk_fact
    subroutine mld_diluk_fact(fill_in,ialg,a,l,u,d,info,blck)
      import psb_dspmat_type, psb_dpk_
      integer, intent(in)                  :: fill_in,ialg
      integer, intent(out)                 :: info
      type(psb_dspmat_type),intent(in)    :: a
      type(psb_dspmat_type),intent(inout) :: l,u
      type(psb_dspmat_type),intent(in), optional, target :: blck
      real(psb_dpk_), intent(inout)        ::  d(:)
    end subroutine mld_diluk_fact
  end interface

  interface mld_ilut_fact
    subroutine mld_dilut_fact(fill_in,thres,a,l,u,d,info,blck)
      import  psb_dspmat_type, psb_dpk_
      integer, intent(in)                  :: fill_in
      real(psb_dpk_), intent(in)           :: thres
      integer, intent(out)                 :: info
      type(psb_dspmat_type),intent(in)    :: a
      type(psb_dspmat_type),intent(inout) :: l,u
      type(psb_dspmat_type),intent(in), optional, target :: blck
      real(psb_dpk_), intent(inout)        ::  d(:)
    end subroutine mld_dilut_fact
  end interface

end module mld_d_ilu_fact_mod
