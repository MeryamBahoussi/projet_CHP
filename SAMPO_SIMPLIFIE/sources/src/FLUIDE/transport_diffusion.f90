module m_transport_diffusion
  use parametres_fluide
  use chimie

  implicit none

contains
  
  ! calcul des coefficients de transport (viscosite dynamique, conductivite thermique)
  subroutine calcul_coeff_transport_b(ldeb, lfin, mdeb, mfin, fld, b)
    integer         , intent(in)  :: ldeb, lfin, mdeb, mfin
    type(STR_fluide), intent(in)  :: fld
    type (STR_BLOC) , pointer     :: b

    integer :: l, m, i, j, nb_especes
    real*8  :: temp, Cp
    real*8  :: phi, somme
    real*8  :: x_esp(1:fld%nb_esp)         ! fractions molaires
    real*8, parameter :: T_ref  = 273.15d0

    nb_especes = fld%nb_esp

    ! choix de la loi de viscosite
    select case (fld%loi_viscosite)
    case (CONSTANTE)
       b%mu(ldeb:lfin, mdeb:mfin) = fld%mu_cst 
    case (SUTHERLAND)
       do m = mdeb, mfin
          do l = ldeb, lfin
             temp = b%T(l,m) / T_ref
             b%mu(l,m)   = 1.711d-5 * temp**1.5d0 * (1.d0 + 110.4d0 / T_ref)  / (temp + 110.4d0 / T_ref)
          end do
       end do
    case (WILKE)
!!! cas d'un mélange multiespeces : formule de Wilke modifiée par Bird (Thèse Velghe (Annexe A))
       do m = mdeb, mfin
          do l = ldeb, lfin

             b%mu(l,m) = 0.d0

             call xesp_from_cesp(fld%especes, b%c_esp(:,l,m), x_esp)
             do i = 1, nb_especes
                phi = 0.d0
                somme = 0.d0

                do j = 1, nb_especes
                   phi = ( 1.d0 + sqrt(fld%especes(i)%Mu / fld%especes(j)%Mu) * &
                        (fld%especes(j)%masse / fld%especes(i)%masse)**.25d0 )**2.d0 
                   phi = phi / sqrt( 8.d0 *(1.d0 + (fld%especes(i)%masse/fld%especes(j)%masse) ) )
                   somme = somme + x_esp(j) * phi
                end do
                
                b%mu(l,m) = b%mu(l,m) + x_esp(i) * fld%especes(i)%Mu / somme 
             end do
          end do
       end do
    case default
       print*, "loi de viscosite inconnue", fld%loi_viscosite
       call arret_code
    end select
    
    ! choix de la loi de conductivite thermique
    select case (fld%loi_conductivite)
    case (CONSTANTE)
       b%K(ldeb:lfin, mdeb:mfin) = fld%K_cst
    case (PRANDTL)
       do m = mdeb, mfin
          do l = ldeb, lfin
             Cp = b%gamma_eq(l,m) * b%Cv(l,m) 
             b%K(l,m) = b%mu(l,m) * Cp / fld%Pr
          end do
       end do
    case default
       print*, "loi de conductivite inconnue", fld%loi_conductivite
       call arret_code
    end select
    
    ! choix de la loi de diffusion
    if (nb_especes > 1) then
       select case (fld%loi_diffusion)
       case (LEWIS) 
          ! loi de Fick (Lewis constant)
          do m = mdeb, mfin
             do l = ldeb, lfin
                Cp = b%gamma_eq(l,m) * b%Cv(l,m)
                b%D(l,m) = fld%Le * b%K(l,m) / (b%rho(l,m) * Cp)
             end do
          end do
       case default
          print*, "loi de diffusion inconnue", fld%loi_diffusion
          call arret_code
       end select
    end if

  end subroutine calcul_coeff_transport_b

  subroutine calcul_coeff_transport_loc(fld, T, c_esp, rho, Cp, mu, K, D)
    type(STR_fluide), intent(in)  :: fld
    real*8          , intent(in)  :: T
    real*8          , intent(in)  :: c_esp(1:fld%nb_esp)
    real*8          , intent(in)  :: rho ! rho1 ou rho2 pas le rho du melange
    real*8          , intent(in)  :: Cp
    real*8          , intent(out) :: mu, K, D

    integer :: i, j
    real*8  :: temp
    real*8  :: phi, somme
    real*8  :: x_esp(1:fld%nb_esp)         ! fractions molaires
    real*8, parameter :: T_ref  = 273.15d0

    ! choix de la loi de viscosite
    select case (fld%loi_viscosite)
    case (CONSTANTE)
       mu = fld%mu_cst 
    case (SUTHERLAND)
       temp = T / T_ref
       mu = 1.711d-5 * temp**1.5d0 * (1.d0 + 110.4d0 / T_ref)  / (temp + 110.4d0 / T_ref)
    case (WILKE)
!!! cas d'un mélange multiespeces : formule de Wilke modifiée par Bird (Thèse Velghe (Annexe A))
       mu = 0.d0
       call xesp_from_cesp(fld%especes, c_esp(:), x_esp)
       do i = 1, fld%nb_esp
          phi = 0.d0
          somme = 0.d0
          
          do j = 1, fld%nb_esp
             phi = ( 1.d0 + sqrt(fld%especes(i)%Mu / fld%especes(j)%Mu) * &
                  (fld%especes(j)%masse / fld%especes(i)%masse)**.25d0 )**2.d0 
             phi = phi / sqrt( 8.d0 *(1.d0 + (fld%especes(i)%masse/fld%especes(j)%masse) ) )
             somme = somme + x_esp(j) * phi
          end do
          
          mu = mu + x_esp(i) * fld%especes(i)%Mu / somme 
       end do
    case default
       print*, "loi de viscosite inconnue", fld%loi_viscosite
       call arret_code
    end select
    
    ! choix de la loi de conductivite thermique
    select case (fld%loi_conductivite)
    case (CONSTANTE)
       K = fld%K_cst
    case (PRANDTL)
       K = mu * Cp / fld%Pr
    case default
       print*, "loi de conductivite inconnue", fld%loi_conductivite
       call arret_code
    end select
    
    ! choix de la loi de diffusion
    if (fld%nb_esp > 1) then
       select case (fld%loi_diffusion)
       case (LEWIS)
          ! loi de Fick (Lewis constant)
          D = fld%Le * K / (rho * Cp)
       case default
          print*, "loi de diffusion inconnue", fld%loi_diffusion
          call arret_code
       end select
    end if
  end subroutine calcul_coeff_transport_loc
  
  subroutine calcul_coeff_transport_diphasique(ldeb, lfin, mdeb, mfin, fld1, fld2, b)
    integer         , intent(in)  :: ldeb, lfin, mdeb, mfin
    type(STR_fluide), intent(in)  :: fld1, fld2
    type (STR_BLOC) , pointer     :: b

    integer :: l,m
    real*8  :: Cp
    real*8  :: mu1, K1, D1
    real*8  :: mu2, K2, D2

    do m = mdeb, mfin
       do l = ldeb, lfin
          
          Cp = b%gamma1(l,m) * b%Cv1(l,m)
          call calcul_coeff_transport_loc(fld1, b%T1(l,m), b%c_esp(:,l,m), b%rho1(l,m), Cp, mu1, K1, D1)
          Cp = 1.d80 ! diphasique + chimie non gere
          call calcul_coeff_transport_loc(fld2, b%T2(l,m), b%c_esp(:,l,m), b%rho2(l,m), Cp, mu2, K2, D2)
          
          b%mu(l,m) = b%z(l,m)*mu1 + (1.d0-b%z(l,m))*mu2
          b%K(l,m)  = b%z(l,m)*K1  + (1.d0-b%z(l,m))*K2
          b%D(l,m)  = 1.d80 ! diphasique + chimie non gere
       end do
    end do
  end subroutine calcul_coeff_transport_diphasique
  
end module m_transport_diffusion
