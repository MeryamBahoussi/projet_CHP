module EquationOfState
  use m_struct
  use parametres_fluide
  use m_MPI

contains
!!!-------------------------------------
!!!       Loi d'etat : Relations 
!!!-------------------------------------
!
!!! p = p(rho,epsilon)
!!! T = T(rho,p)
!
!!! rho = rho(p,T)    ! pour avoir la densite à la paroi dans le bilan d energie pour la fusion
!!! eps = eps(rho,T)  ! pour calculer les derivees de T par rapport a U (premier visqueux)
!!! eps = eps(p,rho)  ! pour l'encodage
!!! eps = eps(z, y, T, rho)
!!!
!!! eps = epsilon = energie interne

  function EOS_p_from_rho_eps(EOS, rho, eps) result (p)
    type(STR_EOS) :: EOS
    real*8        :: rho, eps
    
    real*8 :: p
    
    p = rho*(eps-EOS%q)*(EOS%gamma-1.d0)-EOS%gamma*EOS%pi
  end function EOS_p_from_rho_eps
  
  function EOS_T_from_p_rho(EOS, p, rho) result (T)
    type(STR_EOS) :: EOS
    real*8        :: p, rho
    
    real*8 :: T
    
    T = EOS%T0 + 1.d0/EOS%Cv * ( (p+EOS%gamma*EOS%pi)/(rho*(EOS%gamma-1.d0)) &
         + (EOS%q-EOS%eps0)-deps_dv(rho,EOS) )
 
!!$    if (EOS%nom == stiffened_gas) then
!!$       T = (p+EOS%pi)/(EOS%cv*rho*(EOS%gamma-1.d0))
!!$    end if
  end function EOS_T_from_p_rho

  function EOS_rho_from_p_T(EOS,p,T) result (rho)
    type(STR_EOS) :: EOS
    real*8        :: p, T
    real*8        :: rho
    
    select case(EOS%nom)
    case( GAZ_PARFAIT )
       rho = p/(EOS%Cv * T * (EOS%gamma-1.d0))
    case( STIFFENED_GAS )
       rho = (p+EOS%pi)/(EOS%Cv * T * (EOS%gamma-1.d0))
    case( MIE_GRUNEISEN )
       ! pas bon, le f(rho) de la formule 4.52 de Manu a disparu
       ! pour trouver la bonne densite, il faudrait faire un Newton
       ! gPi = EOS%p0 + EOS%gamma*EOS%pi
       ! rho = EOS%rho0 * (p + EOS%gamma*EOS%pi) / &
       !      (EOS%rho0 * (EOS%gamma-1.d0) * EOS%Cv * (T - EOS%T0) + gPi)

       call newton_rho(EOS, T, p, rho)
    case default
       print*, "EOS non geree : EOS_rho_from_p_T"
       call arret_code
    end select

  end function EOS_rho_from_p_T

  function EOS_eps_from_rho_T(EOS, rho, T) result (eps)
    type(STR_EOS) :: EOS
    real*8        :: rho, T
    
    real*8 :: eps
    
    eps = EOS%eps0 + deps_dv(rho,EOS) + EOS%cv*(T-EOS%T0)   
  end function EOS_eps_from_rho_T
  
  function EOS_eps_from_p_rho(EOS, p, rho) result (eps)
    type(STR_EOS) :: EOS
    real*8        :: p, rho

    real*8 :: eps

    eps = (p + EOS%gamma*EOS%pi)/(rho*(EOS%gamma-1.d0)) + EOS%q
  end function EOS_eps_from_p_rho


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Derivees partielles  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function dT_dp(EOS, rho)
    type(STR_EOS) :: EOS
    real*8        :: rho 
    
    real*8 :: dT_dp
    
    dT_dp = 1.d0/(EOS%Cv*rho*(EOS%gamma-1.d0))
  end function dT_dp

  function deps_dv(rho, EOS)
    real*8,intent(in) :: rho
    real*8 :: deps_dv
    real*8 :: gamma, pi, rho0, T0, p0, Cv
    type(STR_EOS), intent(in) :: EOS

    gamma = EOS%gamma
    pi = EOS%pi
    rho0 = EOS%rho0
    T0 = EOS%T0
    p0 = EOS%p0
    Cv = EOS%Cv 

    select case(EOS%nom)
    case (gaz_parfait)
       deps_dv = 0.d0
    case (stiffened_gas)
       deps_dv = pi/rho
    case (mie_gruneisen)
       deps_dv = pi/rho*(1.d0-(rho/rho0)**gamma) + &
            (Cv*T0-(gamma*pi+p0)/(rho0*(gamma-1.d0)))*(1.d0-(rho/rho0)**(gamma-1.d0))
    case default
       print*, "EOS non geree : deps_dv"
       call arret_code
    end select

  end function deps_dv
  
  subroutine newton_rho(EOS, T, p, rho)
    implicit none
    
    type(STR_EOS), intent(in)  :: EOS
    real*8       , intent(in)  :: T, p
    real*8       , intent(out) :: rho

    real*8 :: rho0, rho_k, rho_kp1
    real*8 :: g, dg
    real*8 :: residu, tol

    rho0=EOS%rho0
    rho_k=rho0
    residu = 1.d70
    tol=1.d-8

    do while (residu>tol)
       g = rho_k*EOS%Cv*T - (p+EOS%pi)/(EOS%gamma-1.d0) - &
            (rho_k/rho0)**EOS%gamma*(rho0*(EOS%Cv*EOS%T0+EOS%q-EOS%eps0)+EOS%pi) 
       
       dg = EOS%Cv*T - EOS%gamma*(rho_k/rho0)**(EOS%gamma-1.d0)*(EOS%Cv*EOS%T0+EOS%q-EOS%eps0+EOS%pi/rho0)
            
       rho_kp1 = rho_k - g/dg
       
       residu = abs(rho_kp1-rho_k)/abs(rho0) 
       rho_k = rho_kp1
    end do
    
    rho = rho_kp1

  end subroutine newton_rho

end module EquationOfState
