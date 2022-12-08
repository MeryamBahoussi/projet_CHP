module thermo_fluide
  
  use decodage_chimie

  implicit none


contains
  
!!! routine qui ne calcule que la densite et l'enthalpie pour un melange de Gaz Parfaits
!!! utile pour la routine qui calcule le flux impose lorsqu'il y a injection de carbone
  subroutine calc_rho_h_from_p_T_Ci(p,T,Ci,rho,h, Cp)
    implicit none 
    real*8, intent(in)  :: p, T, Ci(1:fluide1%nb_esp)
    real*8, intent(out) :: rho, h, Cp

    real*8 :: Cv, gamma

    rho = 1.d0 ! pour ne pas diviser par 0 dans le premier calcul de D
    call calcul_Cv_Cp_fige(fluide1, fluide1%especes, ci, Cv, Cp)
    gamma = Cp/Cv

    rho =  p / ( (gamma-1.d0) * Cv * T )
    h = Cp * T ! equivalent a ein + p / rho
  end subroutine calc_rho_h_from_p_T_Ci

  
  function eps_from_y_T_rho(z, y, T, rho) result (eps)
    real*8 :: z, y, T, rho
    real*8 :: eps
  
    real*8 :: rho1, f1
    real*8 :: rho2, f2
    real*8 :: Cv
    
    rho1=0.d0; rho2=0.d0
    f1=0.d0; f2=0.d0
    if (abs(z) > tol_z .and. abs(y) > tol_z) then
       rho1 = rho*y/z
       f1 = deps_dv(rho1, fluide1%EOS)
    end if
    if (abs(1.d0-z) > tol_z .and. abs(1.d0-y) > tol_z) then
       rho2 = rho*(1.d0-y)/(1.d0-z)
       f2 = deps_dv(rho2, fluide2%EOS)
    end if

    Cv = y*fluide1%EOS%Cv + (1.d0-y)*fluide2%EOS%Cv
    
    eps = Cv*T + y* (fluide1%EOS%eps0+f1-fluide1%EOS%Cv*fluide1%EOS%T0) + &
         (1.d0-y) * (fluide2%EOS%eps0+f2-fluide2%EOS%Cv*fluide2%EOS%T0)
    
  end function eps_from_y_T_rho
  
  subroutine calc_Ec_from_variable_reconstruite(b, direction, ordre2)
    use parametres_globaux, only : DIR_X, DIR_Y
   
    type (STR_BLOC)            , pointer       :: b
    integer                    , intent(in)    :: direction
    type(variable_reconstruite), intent(inout) :: ordre2
   
    integer :: hl, hm
    integer :: l, m
    real*8  :: rhoyci(1:b%ne)
    real*8  :: rhoy2, ux, uy, p, z 
    real*8  :: E,c
    
    select case(direction)
    case (DIR_X)
       hl = 1; hm = 0
    case (DIR_Y)
       hl = 0; hm = 1
    end select

    ! boucle sur les cellules (+1 maille fictive)
    do m = b%md-hm, b%mf+hm
       do l = b%ld-hl, b%lf+hl
          
          ! valeur de gauche pour une cellule donnee
          rhoyci = ordre2%rhoyCiL(:,l,m)
          rhoy2  = ordre2%rhoy2L(l,m)
          ux = b%ordre2%uxL(l,m)
          uy = b%ordre2%uyL(l,m)
          p = ordre2%pL(l,m)
          z = ordre2%zL(l,m)
          
          call E_c_from_rhoyci_rhoy2_uvpz(rhoyci,rhoy2,ux,uy,p,z,E,c)
          
          ordre2%EL(l,m) = E
          ordre2%cL(l,m) = c
          
          ! valeur de droite pour une cellule donnee
          rhoyci = ordre2%rhoyCiR(:,l,m)
          rhoy2  = ordre2%rhoy2R(l,m)
          ux = b%ordre2%uxR(l,m)
          uy = b%ordre2%uyR(l,m)
          p = ordre2%pR(l,m)
          z = ordre2%zR(l,m)
          
          call E_c_from_rhoyci_rhoy2_uvpz(rhoyci,rhoy2,ux,uy,p,z,E,c)
          
          ordre2%ER(l,m) = E
          ordre2%cR(l,m) = c
       end do
    end do
   
  end subroutine calc_Ec_from_variable_reconstruite

  subroutine E_c_from_rhoyci_rhoy2_uvpz(rhoyci,rhoy2,ux,uy,p,z,E,c)
    real*8, dimension(:), intent(in)  :: rhoyci
    real*8              , intent(in)  :: rhoy2, ux, uy, p, z
    real*8              , intent(out) :: E, c
    
    integer :: ierreur, i, ne
    real*8  :: rhoy, rho, y, rho1, eps
    real*8  :: ci(ubound(rhoyci,1))
    real*8  :: gamma, pi, cv, cp
    type(STR_fluide) :: melange
    
    ne = fluide1%ne
    
    rhoy = sum(rhoyci(:))
    rho = rhoy + rhoy2
    y = rhoy/rho
    if (rhoy>0.d0) then
       do i = 1, ne
          ci(i) = rhoyci(i)/rhoy
       end do
    else if (rhoy == 0.d0) then
       ci(:) = 1.d0
    else
       if (rhoy > -1.d-8) then
          ci(:) = 1.d0
       else
          print*, "E_c_from_rhoyci_rhoy2_uvpz : rhoy<0"
          call arret_code
       end if
    end if
    
    if (z>1.d0 .or. z< 0.d0) then
       print*, "z reconstruit n'est pas entre 0 et 1"
    end if
    
    rho1 = 0.d0
    if (abs(z) > tol_z) rho1 = rhoy/z

    if (fluide1%type == CHIMIE_FIGEE) then
       ! permet de calculer cvi
       call calcul_Cv_Cp_fige(fluide1, fluide1%especes, ci(:), cv, Cp)
       !call calcul_coeff_transport(fluide1, fluide1%especes, ci(:,l,m), b%rho1(l,m), fluide1%mu, fluide1%K, fluide1%D)
    end if
    call grandeurs_melange(fluide1, fluide2, z, y, melange)
    eps = EOS_eps_from_p_rho(melange%EOS, p, rho)
    E = eps + 0.5d0 * (ux**2.d0 + uy**2.d0)
    
    select case (MODELE_DIPHASIQUE)
    case(MONOPHASIQUE, ALLAIRE)
       gamma = melange%EOS%gamma
       pi = melange%EOS%pi
       call vitesse_du_son_melange_Allaire(gamma, pi, p, rho, z, c, ierreur)
    end select

  end subroutine E_c_from_rhoyci_rhoy2_uvpz
  
  subroutine pression_partielle(b, i_espece, p_i)
    type (STR_BLOC)           , pointer     :: b
    integer                   , intent(in)  :: i_espece
    real*8, dimension(-1:,-1:), intent(out) :: p_i
    
    integer :: l, m
    
    !  on pourrait le calculer que sur les deux bandes a la paroi
    do m = b%md, b%mf
       do l = b%ld, b%lf
          p_i(l,m) = b%rho1(l,m)*b%c_esp(i_espece,l,m)*fluide1%especes(i_espece)%R*b%T1(l,m)
       end do
    end do
  end subroutine pression_partielle

end module thermo_fluide
