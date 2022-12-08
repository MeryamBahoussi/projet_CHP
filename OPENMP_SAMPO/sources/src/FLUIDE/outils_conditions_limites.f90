module outils_conditions_limites
  use m_MPI
  
  implicit none

  integer, private :: choix_p=4
  integer, private :: choix_U=1
  

contains

  subroutine pression_entree_subsonique_2d(p_e, U_i, p_i, U_inf, p_n, n, a, fac)
    implicit none
     
    real*8, intent(inout) :: p_e         ! p de la maille fictive
    real*8, intent(in)    :: U_i(2), p_i ! u, v, p du point interieur
    real*8, intent(in)    :: U_inf(2)    ! u, v de la condition amont
    real*8, intent(in)    :: p_n         ! p de l'instant n
    real*8, intent(in)    :: n(2)        ! normale
    real*8, intent(in)    :: a, fac      ! vitesse du son Lag, fac = c*dt/aire
      

!!! equation caracteristique
!!! dtau p - a dtau u - a ( dm p - a dm u) = 0
!!! a = rho*c
!!! lambda = a_0^n dt/dm

    select case(choix_p)
    case(1) ! Euler implicite (fonctionne bump)
!!! p_0^n+1 = (p_0^n + lambda * (p_1^n+1-a^n_0*(u^n+1_1-u_inf)) )/(1+lambda)
       p_e = ( p_n + fac*(p_i - a*dot_product(U_i(:)-U_inf(:),n(:))) )/(1.d0+fac)
       
    case(3) ! invariance le long de la caracteristique (fonctionne bump)
       p_e = p_i - a*dot_product((U_i-U_inf),n(:))
       
    case(4) ! approche naive (fonctionne bump)
       p_e = p_i
       
    case default
       call arret_code
    end select
  end subroutine pression_entree_subsonique_2d

  subroutine temperature_entree_subsonique_2d(T_e, T_i, U_i, n, infi, a)
    use m_struct
    use decodage, only : dT_dp

    real*8         , intent(inout) :: T_e       ! T de la maille fictive
    real*8         , intent(in)    :: T_i       ! T du point interieur
    real*8         , intent(in)    :: U_i(2)    ! vitesse du point interieur
    real*8         , intent(in)    :: n(2)      ! normale
    type (STR_cell), pointer       :: infi      ! condition amont
    real*8         , intent(in)    :: a         ! vitesse du son Lag 

    select case(choix_p)
       !    case(1)
       ! que faire
    case(3)
       ! je pense que ca ne marche pas si la condition amont n'est pas une phase pure
       T_e = T_i - a*( (U_i(1)-infi%u)*n(1) + (U_i(2)-infi%v)*n(2) ) &
            * dT_dp(infi%fluide%EOS, infi%rho)
    case(4)
       T_e = T_i

    case default
       print*, "temperature dans la maille fictive non definie pour ce choix de pression dans la condition d'entree subsonique"
       call arret_code
    end select

  end subroutine temperature_entree_subsonique_2d

  !!! matrice telle que (u_e, v_e, p_e)^n+1 = M (u_i, v_i, p_i)^n+1
  subroutine mat_pression_entree_subsonique_2d(a, n, fac, mat)
    
    real*8                , intent(in)  :: a, fac
    real*8, dimension(2)  , intent(in)  :: n
    real*8, dimension(3,3), intent(out) :: mat
    mat(:,:) = 0.d0

    select case(choix_p)
    case(1) ! Euler implicite
       mat(3,1) = -a*n(1)*fac/(1.d0+fac)
       mat(3,2) = -a*n(2)*fac/(1.d0+fac)
       mat(3,3) = fac/(1.d0+fac)

    case(2) ! Euler explicite
       ! tout est au second membre

    case(3) ! invariance le long de la caracteristique
       mat(3,1) = -a*n(1)
       mat(3,2) = -a*n(2)
       mat(3,3) = 1.d0

    case(4) ! approche naive  
       mat(3,3) = 1.d0

    end select

  end subroutine mat_pression_entree_subsonique_2d

  subroutine temperature_entree_subsonique(T_e, T_i, u_i, infi, a)
    use m_struct
    use decodage, only : dT_dp

    real*8         , intent(inout) :: T_e       ! T de la maille fictive
    real*8         , intent(in)    :: T_i, u_i  ! T, u du point interieur
    type (STR_cell), pointer       :: infi      ! condition amont
    real*8         , intent(in)    :: a         ! vitesse du son Lag 
  
    select case(choix_p)
!    case(1)
       ! que faire
    case(3)
       ! je pense que ca ne marche pas si la condition amont n'est pas une phase pure
       T_e = T_i - a*(u_i-infi%u) * dT_dp(infi%fluide%EOS, infi%rho)
    case(4)
       T_e = T_i
       
    case default
       print*, "temperature dans la maille fictive non definie pour ce choix de pression dans la condition d'entree subsonique"
       call arret_code
    end select

  end subroutine temperature_entree_subsonique

!!! matrice telle que (u_e, v_e, p_e)^n+1 = M (u_i, v_i, p_i)^n+1
  subroutine mat_pression_entree_subsonique(a, fac, mat)
    
    real*8                , intent(in)  :: a, fac
    real*8, dimension(3,3), intent(out) :: mat
    mat(:,:) = 0.d0

    select case(choix_p)
    case(1) ! Euler implicite
       mat(3,1) = -a*fac/(1.d0+fac)
       mat(3,3) = fac/(1.d0+fac)

    case(2) ! Euler explicite
       ! tout est au second membre

    case(3) ! invariance le long de la caracteristique
       mat(3,1) = -a
       mat(3,3) = 1.d0

    case(4) ! approche naive  
       mat(3,3) = 1.d0

    end select

  end subroutine mat_pression_entree_subsonique

!!! calcul de la vitesse en sortie a partir de l'interieur
  subroutine vitesse_sortie_subsonique(U_e, U_i, U_n, p_i, p_inf, a, fac)
    
    real*8, intent(out)   :: U_e(2)      ! vitesse de la maille fictive
    real*8, intent(in)    :: U_i(2), p_i ! u, v, p du point interieur
    real*8, intent(in)    :: U_n(2)      ! vitesse de l'instant n
    real*8, intent(in)    :: p_inf       ! pression de la condition aval
    real*8, intent(in)    :: a, fac      ! vitesse du son Lag, fac = c*dt/aire
      

!!! equation caracteristique
!!! dtau p + a dtau u + a ( dm p + a dm u) = 0
!!! a = rho*c
!!! lambda = a_0^n dt/dm

    select case (choix_U)
    case(1)
       U_e(1) = ( U_n(1) - fac/a*(p_inf-p_i - a*U_i(1)) )/(1.d0+fac)
       U_e(2) = U_i(2)
              
    case(2)
       U_e(:) = U_i(:)
       
    case(3)
       U_e(1) = (p_inf - p_i)/a + U_i(1)
       U_e(2) = U_i(2)
    end select
    
  end subroutine vitesse_sortie_subsonique

!!! matrice telle que (u_e, v_e, p_e)^n+1 = M (u_i, v_i, p_i)^n+1
  subroutine mat_vitesse_sortie_subsonique(a, fac, mat)
    
    real*8                , intent(in)  :: a, fac
    real*8, dimension(3,3), intent(out) :: mat
    mat(:,:) = 0.d0
    
    select case(choix_U)
    case(1) 
       mat(1,1) = fac/(1.d0+fac)
       mat(2,2) = 1.d0
       mat(1,3) = fac/(a*(1.d0+fac))
    
    case(2) ! approche naive  
       mat(1,1) = 1.d0
       mat(2,2) = 1.d0        
    
    case(3) ! invariance le long de la caracteristique
       mat(1,1) = 1.d0
       mat(2,2) = 1.d0
       mat(1,3) = -1.d0/a   
    end select

  end subroutine mat_vitesse_sortie_subsonique

  function cond_limite_injection_temps(y, t)
    real*8 :: y, t
    real*8 :: cond_limite_injection_temps
    
    cond_limite_injection_temps = 100.d0*t*(y-0.1d0)*(0.2d0-y)
    
  end function cond_limite_injection_temps

end module outils_conditions_limites
