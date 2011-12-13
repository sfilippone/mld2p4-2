module mld_z_ilu_fact_mod

  use psb_base_mod, only : psb_zspmat_type, psb_dpk_
  use mld_base_prec_type 

  interface mld_ilu0_fact
    subroutine mld_zilu0_fact(ialg,a,l,u,d,info,blck,upd)
      import psb_zspmat_type, psb_dpk_
      integer, intent(in)                 :: ialg
      integer, intent(out)                :: info
      type(psb_zspmat_type),intent(in)    :: a
      type(psb_zspmat_type),intent(inout) :: l,u
      type(psb_zspmat_type),intent(in), optional, target :: blck
      character, intent(in), optional      :: upd
      complex(psb_dpk_), intent(inout)     ::  d(:)
    end subroutine mld_zilu0_fact
  end interface

  interface mld_iluk_fact
    subroutine mld_ziluk_fact(fill_in,ialg,a,l,u,d,info,blck)
      import psb_zspmat_type, psb_dpk_
      integer, intent(in)                  :: fill_in,ialg
      integer, intent(out)                 :: info
      type(psb_zspmat_type),intent(in)    :: a
      type(psb_zspmat_type),intent(inout) :: l,u
      type(psb_zspmat_type),intent(in), optional, target :: blck
      complex(psb_dpk_), intent(inout)        ::  d(:)
    end subroutine mld_ziluk_fact
  end interface

  interface mld_ilut_fact
    subroutine mld_zilut_fact(fill_in,thres,a,l,u,d,info,blck,iscale)
      import  psb_zspmat_type, psb_dpk_
      integer, intent(in)                 :: fill_in
      real(psb_dpk_), intent(in)          :: thres
      integer, intent(out)                :: info
      type(psb_zspmat_type),intent(in)    :: a
      type(psb_zspmat_type),intent(inout) :: l,u
      complex(psb_dpk_), intent(inout)       :: d(:)
      type(psb_zspmat_type),intent(in), optional, target :: blck
      integer, intent(in), optional       :: iscale
    end subroutine mld_zilut_fact
  end interface

end module mld_z_ilu_fact_mod
