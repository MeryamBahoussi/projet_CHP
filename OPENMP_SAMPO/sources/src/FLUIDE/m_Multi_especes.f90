module MULTIESPECES
  use m_struct
  use parametres_fluide

!!!
!!! Attention : ici, on suppose que eps_i = Cv_i T
!!! Si la relation entre eps_i et T est non lineaire, on ne peut plus utiliser la
!!! formule  p = rho*esp*(gamma-1)
!!!
!!! De plus, ici on note cv les cv_figes
!!!

contains
  
!!!-------------------------------------
!!!  Gaz multi especes : relations de fermeture thermo
!!!-------------------------------------
!
!!! p = p(rho, eps)
!!! T = T(rho, p)
  function MESP_p_from_rho_eps(rho, eps, gamma) result (p)
    real*8 :: rho, eps, gamma
    real*8 :: p
    
    p = rho*eps*(gamma-1.d0)
  end function MESP_p_from_rho_eps
  
  function MESP_T_from_p_rho(p, rho, cv, gamma) result (T)
    real*8 :: p, rho, cv, gamma
    real*8 :: T
   
    T = p/(rho*Cv*(gamma-1.d0))
  end function MESP_T_from_p_rho

end module MULTIESPECES
