module mld_d_pde_mod
contains
  !
  ! functions parametrizing the differential equation 
  !  
  function b1(x,y,z)
    use psb_base_mod, only : psb_dpk_, done
    real(psb_dpk_) :: b1
    real(psb_dpk_), intent(in) :: x,y,z
    b1=dzero
  end function b1
  function b2(x,y,z)
    use psb_base_mod, only : psb_dpk_, done
    real(psb_dpk_) ::  b2
    real(psb_dpk_), intent(in) :: x,y,z
    b2=dzero
  end function b2
  function b3(x,y,z)
    use psb_base_mod, only : psb_dpk_, done
    real(psb_dpk_) ::  b3
    real(psb_dpk_), intent(in) :: x,y,z      
    b3=dzero
  end function b3
  function c(x,y,z)
    use psb_base_mod, only : psb_dpk_, done
    real(psb_dpk_) ::  c
    real(psb_dpk_), intent(in) :: x,y,z      
    c=dzero
  end function c
  function a1(x,y,z)
    use psb_base_mod, only : psb_dpk_, done
    real(psb_dpk_) ::  a1   
    real(psb_dpk_), intent(in) :: x,y,z
    a1=done
  end function a1
  function a2(x,y,z)
    use psb_base_mod, only : psb_dpk_, done
    real(psb_dpk_) ::  a2
    real(psb_dpk_), intent(in) :: x,y,z
    a2=done
  end function a2
  function a3(x,y,z)
    use psb_base_mod, only : psb_dpk_, done
    real(psb_dpk_) ::  a3
    real(psb_dpk_), intent(in) :: x,y,z
    a3=done
  end function a3
  function g(x,y,z)
    use psb_base_mod, only : psb_dpk_, done, dzero
    real(psb_dpk_) ::  g
    real(psb_dpk_), intent(in) :: x,y,z
    g = dzero
    if (x == done) then
      g = done
    else if (x == dzero) then 
      g = exp(y**2-z**2)
    end if
  end function g
end module mld_d_pde_mod

