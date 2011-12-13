module mld_s_ilu_fact_mod

  use psb_base_mod, only : psb_sspmat_type, psb_spk_
  use mld_base_prec_type 

  interface mld_ilu0_fact
    subroutine mld_silu0_fact(ialg,a,l,u,d,info,blck,upd)
      import psb_sspmat_type, psb_spk_
      integer, intent(in)                 :: ialg
      integer, intent(out)                :: info
      type(psb_sspmat_type),intent(in)    :: a
      type(psb_sspmat_type),intent(inout) :: l,u
      type(psb_sspmat_type),intent(in), optional, target :: blck
      character, intent(in), optional      :: upd
      real(psb_spk_), intent(inout)     ::  d(:)
    end subroutine mld_silu0_fact
  end interface

  interface mld_iluk_fact
    subroutine mld_siluk_fact(fill_in,ialg,a,l,u,d,info,blck)
      import psb_sspmat_type, psb_spk_
      integer, intent(in)                  :: fill_in,ialg
      integer, intent(out)                 :: info
      type(psb_sspmat_type),intent(in)    :: a
      type(psb_sspmat_type),intent(inout) :: l,u
      type(psb_sspmat_type),intent(in), optional, target :: blck
      real(psb_spk_), intent(inout)        ::  d(:)
    end subroutine mld_siluk_fact
  end interface

  interface mld_ilut_fact
    subroutine mld_silut_fact(fill_in,thres,a,l,u,d,info,blck,iscale)
      import  psb_sspmat_type, psb_spk_
      integer, intent(in)                 :: fill_in
      real(psb_spk_), intent(in)          :: thres
      integer, intent(out)                :: info
      type(psb_sspmat_type),intent(in)    :: a
      type(psb_sspmat_type),intent(inout) :: l,u
      real(psb_spk_), intent(inout)       :: d(:)
      type(psb_sspmat_type),intent(in), optional, target :: blck
      integer, intent(in), optional       :: iscale
    end subroutine mld_silut_fact
  end interface

end module mld_s_ilu_fact_mod
