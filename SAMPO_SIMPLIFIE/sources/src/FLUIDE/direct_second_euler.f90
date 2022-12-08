!!!!!!!!!!!!!!!!
!!!   schéma direct (sans splitting)
!!!   routines de manuel Latige
!!!
!!! Ordre 2 pas top
!!! on reconstruit rhoyCi, rhoy2 et z (stockes dans U) et u, v, et p
!!!
!!!!!!!!!!!!!!!!

module second_euler_direct
  use m_struct
  use m_MPI
  use parametres_fluide
  use parametres_globaux, only : DIR_X, DIR_Y, FLUIDE
  use outils_fluide
  use thermo_fluide
  use m_transport_diffusion
  implicit none

contains
  
  subroutine second_membre_euler_direct(b, sb, sm)

    type (STR_BLOC)                  , pointer :: b
    type (STR_SUPER_BLOC)            , pointer :: sb
    real*8, dimension(1:,b%ld:,b%md:), intent(inout) :: sm
    
    integer :: i, l, m
    integer :: ld, lf, md, mf
    logical :: move
    real*8  :: rayon   
    real*8  :: fff(iz-1)

    ld=b%ld; lf=b%lf
    md=b%md; mf=b%mf 
   
    if (b%maillage%move .and. .not. sb%sol_stat%actif) then 
       move = .true.
       print*, " euler direct + priseEnCompteInterfaceMobileDansFluide non gere"
       call arret_code
    else
       move = .false.
    end if
    
    ! les pentes sont déjà calculées pour le pas de temps
    sm(:,:,:)=0.d0 
    b%F(:,:,:)=0.d0
    call flux_euler_direct(b, DIR_X, sb%solveur_Riemann, move, b%F)
    do m = md, mf
       do l = ld, lf
          sm(1:iz-1,l,m) = sm(1:iz-1,l,m) + (b%F(1:iz-1,l,m)-b%F(1:iz-1,l-1,m))
       end do
    end do
    
    b%F(:,:,:)=0.d0
    call flux_euler_direct(b, DIR_Y, sb%solveur_Riemann, move, b%F)
    do m = md, mf
       do l = ld, lf
          sm(1:iz-1,l,m) = sm(1:iz-1,l,m) + (b%F(1:iz-1,l,m)-b%F(1:iz-1,l,m-1))
       end do
    end do

    !¤¤¤¤¤¤¤¤¤¤¤¤¤¤ Schéma d'advection pour le z
    sm(iz,:,:)=0.d0
    call advection_z_ordre1(b, DIR_X, sb%solveur_Riemann, .false., sm)
    call advection_z_ordre1(b, DIR_Y, sb%solveur_Riemann, .false., sm)
    
  end subroutine second_membre_euler_direct

  subroutine flux_euler_direct( b, direction, solveur_riemann, move, flux )
    use decodage
    
    type (STR_BLOC)                       , pointer       :: b
    integer                               , intent(in)    :: direction 
    integer                               , intent(in)    :: solveur_Riemann
    logical                               , intent(in)    :: move
    real*8, dimension(1:,b%ld-1:, b%md-1:), intent(inout) :: flux
    

    integer :: l,m

    b%ordre2%rhoyCiL = b%U(irho1z1:irho2z2-1,b%ld-1:b%lf+1,b%md-1:b%mf+1)
    b%ordre2%rhoyCiR = b%U(irho1z1:irho2z2-1,b%ld-1:b%lf+1,b%md-1:b%mf+1)
    b%ordre2%rhoy2L = b%U(irho2z2,b%ld-1:b%lf+1,b%md-1:b%mf+1)
    b%ordre2%rhoy2R = b%U(irho2z2,b%ld-1:b%lf+1,b%md-1:b%mf+1)

    !¤¤¤¤¤¤¤¤¤¤¤¤¤¤ Flux centres 
    call flux_centres_direct(b,direction,flux)
    if (direction==DIR_X) then
       do m = b%md, b%mf
          do l = b%ld, b%lf
             b%sm(1:iz-1,l,m) = b%sm(1:iz-1,l,m) + (flux(1:iz-1,l,m)-flux(1:iz-1,l-1,m))
          end do
       end do
    else
       do m = b%md, b%mf
          do l = b%ld, b%lf
             b%sm(1:iz-1,l,m) = b%sm(1:iz-1,l,m) + (flux(1:iz-1,l,m)-flux(1:iz-1,l,m-1))
          end do
       end do
    end if
    flux = 0.d0

    !¤¤¤¤¤¤¤¤¤¤¤¤¤¤ Appel du solveur de Riemann pour les flux decentres
    select case( solveur_riemann )
    case( HLL )
       call solveur_HLL(b,direction,move,flux)
    case( GALLICE )
       call solveur_gallice(b,direction,move,flux)
    case( ROE ) 
       call solveur_ROE(b,direction,move,flux)
    case default
       print*, "solveur_riemann inconnu"
       call arret_code
    end select
    
  end subroutine flux_euler_direct
  
  subroutine flux_centres_direct(b,direction,flux)
 
    type (STR_BLOC), pointer     :: b
    integer        , intent(in)  :: direction 
    real*8, dimension(1:,b%ld-1:,b%md-1:), intent(out) :: flux
    ! on est oblige de specifier le premier indice des tableaux sinon ca commence a 1...
    
    integer :: i, l, m, hl, hm
    integer :: lg, ld, mg, md
    real*8  :: p_g, Un_g, rhog, uxg, uyg, Eg, rhoy2g
    real*8  :: p_d, Un_d, rhod, uxd, uyd, Ed, rhoy2d
    real*8, dimension(1:b%ne) :: rhoyCig, rhoyCid
    real*8, dimension(:,:)  , pointer :: dl
    real*8, dimension(:,:,:), pointer :: n
    real*8, dimension(:,:)  , pointer :: pL, pR, uxL, uxR, uyL, uyR, EL, ER
    real*8, dimension(:,:)  , pointer :: rhoy2L, rhoy2R
    real*8, dimension(:,:,:), pointer :: rhoyCiL, rhoyCiR

    ! Selection de la direction
    select case(direction)
    case (DIR_X) ! direction l
       n => b%n_dir_l
       dl => b%dl_l
       hl=1; hm=0
    case (DIR_Y) ! direction m
       n => b%n_dir_m
       dl => b%dl_m
       hl=0; hm=1
    end select
    
    rhoyCiL => b%ordre2%rhoyCiL; rhoyCiR => b%ordre2%rhoyCiR
    rhoy2L  => b%ordre2%rhoy2L ; rhoy2R  => b%ordre2%rhoy2R
       
    uxL=> b%u_x; uxR => b%u_x
    uyL=> b%u_y; uyR => b%u_y
    pL => b%p  ; pR  => b%p
    EL => b%E  ; ER  => b%E
        
    ! boucle sur les faces
    do m = b%md-hm, b%mf 
       do l = b%ld-hl, b%lf
          lg = l; ld = l+hl
          mg = m; md = m+hm
          
          rhoyCig(:) = rhoyCiR(:,lg,mg)
          rhoyCid(:) = rhoyCiL(:,ld,md)
          rhoy2g = rhoy2R(lg,mg)
          rhoy2d = rhoy2L(ld,md)
          uxg = uxR(lg,mg)
          uxd = uxL(ld,md)
          uyg = uyR(lg,mg)
          uyd = uyL(ld,md)
          p_g = pR(lg,mg)
          p_d = pL(ld,md)
          Eg = ER(lg,mg)
          Ed = EL(ld,md)
          
          rhog = sum(rhoyCig(:)) + rhoy2g
          rhod = sum(rhoyCid(:)) + rhoy2d
          
          Un_g = uxg*n(1,l,m) + uyg*n(2,l,m)
          Un_d = uxd*n(1,l,m) + uyd*n(2,l,m)
          
!!! partie centree (F_i+F_i+1)
          do i = 1, b%ne
             flux(irho1z1+i-1,l,m) = rhoyCig(i)*Un_g + rhoyCid(i)*Un_d
          end do
          flux(irho2z2,l,m) = rhoy2g*Un_g                  + rhoy2d*Un_d
          flux(irhou  ,l,m) = rhog*uxg*Un_g + p_g*n(1,l,m) + rhod*uxd*Un_d +p_d*n(1,l,m)
          flux(irhov  ,l,m) = rhog*uyg*Un_g + p_g*n(2,l,m) + rhod*uyd*Un_d +p_d*n(2,l,m)
          flux(irhoE  ,l,m) = (rhog*Eg+p_g)*Un_g           + (rhod*Ed+p_d)*Un_d
    
          flux(:,l,m)= 0.5d0 * flux(:,l,m) * dl(l,m)
       end do
    end do
    
  end subroutine flux_centres_direct
  
  subroutine solveur_gallice(b,direction,move,flux)

    type (STR_BLOC)                      , pointer       :: b
    integer                              , intent(in)    :: direction 
    logical                              , intent(in)    :: move
    real*8, dimension(1:,b%ld-1:,b%md-1:), intent(inout) :: flux

    integer :: i, nbe        ! (Multi-especes)
    integer :: l, m, hl, hm
    integer :: lg, ld, mg, md
    real*8 :: alpha, p_beta, u_alpha, rho
    real*8 :: p_g, Un_g, Un_s
    real*8 :: p_d, Un_d
    real*8 :: phim, Rm(b%neq)
    real*8 :: phip, Rp(b%neq)
    real*8 :: rho_g, rho_gs, rho_ds, rho_d
    real*8 :: y2_g, y2_gs, y2_ds, y2_d
    real*8, dimension(1:b%neq) :: V_g, V_gs, V_ds, V_d, delta_U
    real*8 :: epsilon, T(2), coef_epsilon
    real*8 :: Lambda1, Lambda2, Lambda3
    real*8, dimension(:,:), pointer :: Cm, Cp, dl, v_mesh
    real*8, dimension(:,:,:), pointer :: n
    real*8, dimension(:,:)  , pointer :: pL, pR, uxL, uxR, uyL, uyR, EL, ER
    real*8, dimension(:,:)  , pointer :: rhoy2L, rhoy2R
    real*8, dimension(:,:,:), pointer :: rhoyCiL, rhoyCiR
    
    select case(direction)
    case (DIR_X) ! direction l
       Cm => b%Cm_l
       Cp => b%Cp_l
       n => b%n_dir_l
       dl => b%dl_l
       hl=1; hm=0
       if (move .eqv. .true.) then
          v_mesh => b%dX_dt
       else
          v_mesh => b%zeros
       end if
    case (DIR_Y) ! direction m
       Cm => b%Cm_m
       Cp => b%Cp_m
       n => b%n_dir_m
       dl => b%dl_m
       hl=0; hm=1
       if (move .eqv. .true.) then
          v_mesh => b%dY_dt
       else
          v_mesh => b%zeros
       end if
    end select
    
    rhoyCiL => b%ordre2%rhoyCiL; rhoyCiR => b%ordre2%rhoyCiR
    rhoy2L  => b%ordre2%rhoy2L ; rhoy2R  => b%ordre2%rhoy2R

    pL => b%p; pR => b%p
    uxL=> b%u_x; uxR => b%u_x
    uyL=> b%u_y; uyR => b%u_y
    EL => b%E ; ER  => b%E
    
    nbe = b%ne
    ! boucle sur les faces
    do m = b%md-hm, b%mf 
       do l = b%ld-hl, b%lf
          ! (l,m) : indice de l'interface
          ! (lg,mg) : indice de la cellule à gauche de l'interface
          ! (ld,md) : indice de la cellule à droite de l'interface
          lg = l; ld = l+hl
          mg = m; md = m+hm
          
          ! initialisation
          V_g = 0.d0 ; V_gs = 0.d0 ; V_ds = 0.d0 ; V_d = 0.d0
                    
          ! U_to_V
          rho = sum(rhoyCiR(:,lg,mg)) + rhoy2R(lg,mg)
          V_g(ivol) = 1.d0 / rho
          do i = 1, nbe
             V_g(iy+i-1) = rhoyCiR(irho1z1+i-1,lg,mg)*V_g(ivol)
          end do
          V_g(iu) = uxR(lg,mg)
          V_g(iv) = uyR(lg,mg)
          V_g(iE) = ER(lg,mg)

          !Droit
          rho = sum(rhoyCiL(:,ld,md)) + rhoy2L(ld,md)
          V_d(ivol) = 1.d0 / rho
          do i = 1, nbe
             V_d(iy+i-1) = rhoyCiL(irho1z1+i-1,ld,md)*V_d(ivol)
          end do
          V_d(iu) = uxL(ld,md)
          V_d(iv) = uyL(ld,md)
          V_d(iE) = EL(ld,md)
          
          ! Constante des moyennes lagrangienne
          alpha = Cm(l,m)/(Cm(l,m)+Cp(l,m))

          p_g = pR(lg,mg); p_d = pL(ld,md)
          Un_g = V_g(iu)*n(1,l,m) + V_g(iv)*n(2,l,m)
          Un_d = V_d(iu)*n(1,l,m) + V_d(iv)*n(2,l,m)

          p_beta = (1.d0-alpha)*p_g + alpha*p_d
          u_alpha = alpha*Un_g + (1.d0-alpha)*Un_d

          !Calcul de R- et R+
          ! R+ et R-
          Rp = 0.d0
          Rp(iy)   = 0.d0
          Rp(ivol) = -1.d0
          Rp(iu)   = Cp(l,m)*n(1,l,m)
          Rp(iv)   = Cp(l,m)*n(2,l,m)
          Rp(iE)   = p_beta + u_alpha*Cp(l,m)

          Rm = 0.d0
          Rm(iy)   = 0.d0
          Rm(ivol) = -1.d0
          Rm(iu)   = -Cm(l,m)*n(1,l,m)
          Rm(iv)   = -Cm(l,m)*n(2,l,m)
          Rm(iE)   = p_beta - u_alpha*Cm(l,m)

          ! phi+ et phi-
          phip = ((p_d-p_g) + Cm(l,m)*(Un_d-Un_g))/(Cm(l,m)*Cp(l,m)+Cp(l,m)**2.d0)
          phim = ((p_d-p_g) - Cp(l,m)*(Un_d-Un_g))/(Cm(l,m)*Cp(l,m)+Cm(l,m)**2.d0)
          
          ! Etats intermediaires de Lagrange
          !Gauche
          V_gs = V_g + Rm * phim
          !Droit
          V_ds = V_d - Rp * phip

          !Calcul des flux décentrés
          Un_s = V_gs(iu)*n(1,l,m) + V_gs(iv)*n(2,l,m)
          
          !Valeurs propres (u-c, u, u+c) + vitesse du maillage
          Lambda1 = Un_g - Cm(l,m) * V_g(ivol) + v_mesh(l,m)
          Lambda2 = Un_s + v_mesh(l,m)
          Lambda3 = Un_d + Cp(l,m) * V_d(ivol) + v_mesh(l,m)
          
!!! TERME DE CORRECTION DES VALEURS ABSOLUES DES VALEURS PROPRES ###################
          T = (/ -n(2,l,m), n(1,l,m) /) ! Vecteur tangente
          ! Au passage lagrange->Euler, les valeurs propres euler ont 
          ! plusieurs formules possibles avec les états intermédiares.
          epsilon =  0.d0
          epsilon = b%delta_epsilon * Max( abs(Lambda1), abs(Lambda2), abs(Lambda3) )
          ! correction du epsilon pour l'annuler dans la couche limite visqueuse ( il est correct dans la partie euler)
          coef_epsilon = 1.d0
          epsilon = epsilon * coef_epsilon
!!! ###############################################################################
          
          !Flux décentré ( - \sum_k |\lambda_k| \delta U_k)
          delta_U = 0.d0
          
!!! passage Lagrange Euler
          rho_g  = 1.d0/V_g(ivol) ; y2_g  = 1.d0-sum(V_g(iy:iy+nbe-1))
          rho_gs = 1.d0/V_gs(ivol); y2_gs = 1.d0-sum(V_gs(iy:iy+nbe-1))
          rho_ds = 1.d0/V_ds(ivol); y2_ds = 1.d0-sum(V_ds(iy:iy+nbe-1))
          rho_d  = 1.d0/V_d(ivol) ; y2_d  = 1.d0-sum(V_d(iy:iy+nbe-1))
          
          do i = 1, nbe
             delta_U(irho1z1+i-1) = rho_gs*V_gs(iy+i-1) - rho_g*V_g(iy+i-1) 
          end do
          delta_U(irho2z2) = rho_gs*y2_gs    - rho_g*y2_g 
          delta_U(irhou)   = rho_gs*V_gs(iu) - rho_g*V_g(iu)
          delta_U(irhov)   = rho_gs*V_gs(iv) - rho_g*V_g(iv)
          delta_U(irhoE)   = rho_gs*V_gs(iE) - rho_g*V_g(iE)
          
          flux(:,l,m) = flux(:,l,m) - val_abs(epsilon,Lambda1) * delta_U
          
          do i = 1, nbe
             delta_U(irho1z1+i-1) = rho_ds*V_ds(iy+i-1) - rho_gs*V_gs(iy+i-1) 
          end do
          delta_U(irho2z2) = rho_ds*y2_ds    - rho_gs*y2_gs
          delta_U(irhou)   = rho_ds*V_ds(iu) - rho_gs*V_gs(iu)
          delta_U(irhov)   = rho_ds*V_ds(iv) - rho_gs*V_gs(iv)
          delta_U(irhoE)   = rho_ds*V_ds(iE) - rho_gs*V_gs(iE)
          
          flux(:,l,m) = flux(:,l,m) - val_abs(epsilon,Lambda2) * delta_U
          
          do i = 1, nbe
             delta_U(irho1z1+i-1) = rho_d*V_d(iy+i-1) - rho_ds*V_ds(iy+i-1) 
          end do
          delta_U(irho2z2) = rho_d*y2_d    - rho_ds*y2_ds
          delta_U(irhou)   = rho_d*V_d(iu) - rho_ds*V_ds(iu)
          delta_U(irhov)   = rho_d*V_d(iv) - rho_ds*V_ds(iv)
          delta_U(irhoE)   = rho_d*V_d(iE) - rho_ds*V_ds(iE)
          
          flux(:,l,m) = flux(:,l,m) - val_abs(epsilon,Lambda3) * delta_U

!!! assemblage final
          flux(:,l,m)= 0.5d0 * flux(:,l,m) * dl(l,m)
       end do
    end do
  end subroutine solveur_gallice
  
  subroutine solveur_ROE(b,direction,move,flux)
    use variables_globales, only : Htot_infini_amont
    
    type (STR_BLOC)                      , pointer       :: b
    integer                              , intent(in)    :: direction 
    logical                              , intent(in)    :: move
    real*8, dimension(1:,b%ld-1:,b%md-1:), intent(inout) :: flux
    
    integer :: i, ne          ! (Multi-especes)
    integer :: l, m, hl, hm
    integer :: lg, ld, mg, md
    real*8 :: rhog, rhod, HtotR, HtotL, xiR, xiL, xi_roe
    real*8 :: rho_moy, y2_roe, Htot_roe, moyenne_c2, c_star
    real*8 :: gamma_1_g, gamma_1_d, h_g, h_d, h_roe, chig, chid
    real*8 :: coef_R, coef_L
    real*8 :: V_g(2), V_d(2), V_roe(2)
    real*8 :: epsilon, coef_epsilon
    real*8 :: Lambda1, Lambda2, Lambda3
    real*8 :: alpha1, alpha5
    real*8 :: R1(b%neq), R5(b%neq), delta_U(b%neq)
    real*8 :: Ta(2), tab_y1_roe(1:b%ne)
    real*8, dimension(:,:), pointer :: dl, v_mesh
    real*8, dimension(:,:,:), pointer :: n
    real*8, dimension(:,:)  , pointer :: pL, pR, uxL, uxR, uyL, uyR, EL, ER, cL, cR
    real*8, dimension(:,:)  , pointer :: rhoy2L, rhoy2R, zL, ZR
    real*8, dimension(:,:,:), pointer :: rhoyCiL, rhoyCiR
    
    select case(direction)
    case (DIR_X) ! direction l
       n => b%n_dir_l
       dl => b%dl_l
       hl=1; hm=0
       if (move .eqv. .true.) then
          v_mesh => b%dX_dt
       else
          v_mesh => b%zeros
       end if
    case (DIR_Y) ! direction m
       n => b%n_dir_m
       dl => b%dl_m
       hl=0; hm=1
       if (move .eqv. .true.) then
          v_mesh => b%dY_dt
       else
          v_mesh => b%zeros
       end if
    end select
    
    rhoyCiL => b%ordre2%rhoyCiL; rhoyCiR => b%ordre2%rhoyCiR
    rhoy2L  => b%ordre2%rhoy2L ; rhoy2R  => b%ordre2%rhoy2R

    pL => b%p; pR => b%p
    uxL=> b%u_x; uxR => b%u_x
    uyL=> b%u_y; uyR => b%u_y
    EL => b%E ; ER  => b%E
    cL => b%c ; cR  => b%c
    zL => b%z ; zR  => b%z

    ne = b%ne
    
    ! boucle sur les faces
    do m = b%md-hm, b%mf 
       do l = b%ld-hl, b%lf
          ! (l,m) : indice de l'interface
          ! (lg,mg) : indice de la cellule à gauche de l'interface
          ! (ld,md) : indice de la cellule à droite de l'interface
          lg = l; ld = l+hl
          mg = m; md = m+hm
          
          ! initialisation
          rhog = sum(rhoyCiR(:,lg,mg)) + rhoy2R(lg,mg)
          rhod = sum(rhoyCiL(:,ld,md)) + rhoy2L(ld,md)
          
          ! Initialisation de H = E + p/rho
          HtotR = ER(lg,mg) + pR(lg,mg)/rhog
          HtotL = EL(ld,md) + pL(ld,md)/rhod                   
          
          ! Determination des vecteurs vitesse
          V_g(1) = uxR(lg,mg)
          V_g(2) = uyR(lg,mg)
          V_d(1) = uxL(ld,md)
          V_d(2) = uyL(ld,md)
          
          ! Determination des moyennes de ROE
          coef_R = sqrt(rhog) / ( sqrt(rhog) + sqrt(rhod) )
          coef_L = sqrt(rhod) / ( sqrt(rhog) + sqrt(rhod) )
          
          do i = 1, ne
             tab_y1_roe(i) = coef_R*rhoyCiR(irho1z1+i-1,lg,mg)/rhog + coef_L*rhoyCiL(irho1z1+i-1,ld,md)/rhod
          end do
          
          y2_roe = coef_R*rhoy2R(lg,mg)/rhog + coef_L*rhoy2L(ld,md)/rhod
          V_roe(1) = coef_R * V_g(1) + coef_L * V_d(1)
          V_roe(2) = coef_R * V_g(2) + coef_L * V_d(2)
          Htot_roe = coef_R * HtotR  + coef_L * HtotL
          
          ! Determination de rho_moy
          rho_moy = coef_R * rhod + coef_L * rhog 
          
          ! Determination de la variable xi = 1/(gamma-1)
          xiR = zR(lg,mg)/(b%gamma1(lg,mg)-1.d0) + (1.d0-zR(lg,mg))/(fluide2%EOS%gamma-1.d0)
          xiL = zL(ld,md)/(b%gamma1(ld,md)-1.d0) + (1.d0-zL(ld,md))/(fluide2%EOS%gamma-1.d0)
             
          xi_roe = coef_R * xiR + coef_L * xiL
             
          ! Determination de la vitesse du son : cf formule (63) p 598 de Allaire
          moyenne_c2 = coef_R * xiR * cR(lg,mg)**2.d0 + coef_L * xiL * cL(ld,md)**2.d0
          c_star = sqrt( moyenne_c2  / xi_roe )

          ! Determination des vecteurs et des coefficients du solveur de ROE
          alpha1 = ( (pL(ld,md)-pR(lg,mg)) - rho_moy * c_star * ( dot_product(V_d,n(:,l,m))-dot_product(V_g,n(:,l,m)) ) ) / ( 2.0d0 * c_star**2.0d0 )
          
          alpha5 = ( (pL(ld,md)-pR(lg,mg)) + rho_moy * c_star * ( dot_product(V_d,n(:,l,m))-dot_product(V_g,n(:,l,m)) ) ) / ( 2.0d0 * c_star**2.0d0 )
          
          do i = 1, ne
             R1(irho1z1+i-1) = tab_y1_roe(i)
          end do
          R1(irho2z2) = y2_roe
          R1(irhou)   = V_roe(1) - c_star * n(1,l,m)
          R1(irhov)   = V_roe(2) - c_star * n(2,l,m)
          R1(irhoE)   = Htot_roe - dot_product(V_roe,n(:,l,m)) * c_star
          
          do i = 1, ne
             R5(irho1z1+i-1) = tab_y1_roe(i)
          end do
          R5(irho2z2) = y2_roe
          R5(irhou)   = V_roe(1) + c_star * n(1,l,m)
          R5(irhov)   = V_roe(2) + c_star * n(2,l,m)
          R5(irhoE)   = Htot_roe + dot_product(V_roe,n(:,l,m))* c_star
          
          do i = 1, ne
             delta_U(irho1z1+i-1) = rhoyCiL(irho1z1+i-1,ld,md)-rhoyCiR(irho1z1+i-1,lg,mg)
          end do
          delta_U(irho2z2) = rhoy2L(ld,md)   - rhoy2R(lg,mg)
          delta_U(irhou)   = rhod*uxL(ld,md) - rhog*uxR(lg,mg)
          delta_U(irhov)   = rhod*uyL(ld,md) - rhog*uyR(lg,mg)
          delta_U(irhoE)   = rhod* EL(ld,md) - rhog* ER(lg,mg)
          
          ! Calcul des flux
          
          ! Valeur propre
          Lambda1 = dot_product(V_roe,n(:,l,m)) - c_star + v_mesh(l,m) 
          Lambda2 = dot_product(V_roe,n(:,l,m)) + v_mesh(l,m)
          Lambda3 = dot_product(V_roe,n(:,l,m)) + c_star + v_mesh(l,m)

!!! TERME DE CORRECTION DES VALEURS ABSOLUES DES VALEURS PROPRES ##################	
          Ta = (/ -n(2,l,m), n(1,l,m) /)
          epsilon = b%delta_epsilon * ( abs(dot_product(V_roe,n(:,l,m))+v_mesh(l,m)) + abs(dot_product(V_roe,Ta)) + c_star )
          
!!! correction du epsilon pour l'annuler dans la couche limite visqueuse 
!!! ( il est correct dans la partie euler)
          coef_epsilon = 5.d0*(Htot_roe/Htot_infini_amont - b%delta_epsilon_couche_lim) + 1.d0
          coef_epsilon = max(0.0d0, min(1.d0, coef_epsilon))
          epsilon = epsilon * coef_epsilon
!!! ###############################################################################		
          flux(:,l,m) = flux(:,l,m) - val_abs(epsilon,Lambda1) * alpha1 * R1(:)
          
          flux(:,l,m) = flux(:,l,m) - val_abs(epsilon,Lambda2) * ( delta_U(:) - alpha1 * R1(:) - alpha5 * R5(:) )
          
          flux(:,l,m) = flux(:,l,m) - val_abs(epsilon,Lambda3) * alpha5 * R5(:)
          
          flux(:,l,m)= 0.5d0 * flux(:,l,m) * dl(l,m)
       end do
    end do
  end subroutine solveur_ROE

  subroutine solveur_HLL(b,direction,move,flux)
    type (STR_BLOC)                      , pointer       :: b
    integer                              , intent(in)    :: direction 
    logical                              , intent(in)    :: move
    real*8, dimension(1:,b%ld-1:,b%md-1:), intent(inout) :: flux

    integer :: l, m, hl, hm, i, j
    integer :: lg, ld, mg, md
    real*8 :: kappa, Un_g, Un_d, rhog, rhod, delta_U(b%neq)
    real*8, dimension(:,:), pointer :: dl, v_mesh
    real*8, dimension(:,:,:), pointer :: n
    real*8, dimension(:,:)  , pointer :: uxL, uxR, uyL, uyR, cL, cR, EL, ER
    real*8, dimension(:,:)  , pointer :: rhoy2L, rhoy2R
    real*8, dimension(:,:,:), pointer :: rhoyCiL, rhoyCiR
    
    select case(direction)
    case (DIR_X) ! direction l
       n => b%n_dir_l
       dl => b%dl_l
       hl=1; hm=0
       if (move .eqv. .true.) then
          v_mesh => b%dX_dt
       else
          v_mesh => b%zeros
       end if
    case (DIR_Y) ! direction m
       n => b%n_dir_m
       dl => b%dl_m
       hl=0; hm=1
       if (move .eqv. .true.) then
          v_mesh => b%dY_dt
       else
          v_mesh => b%zeros
       end if
    end select
    
    rhoyCiL => b%ordre2%rhoyCiL; rhoyCiR => b%ordre2%rhoyCiR
    rhoy2L  => b%ordre2%rhoy2L ; rhoy2R  => b%ordre2%rhoy2R

    uxL=> b%u_x; uxR => b%u_x
    uyL=> b%u_y; uyR => b%u_y
    EL => b%E  ; ER  => b%E
    cL => b%c  ; cR  => b%c

    ! boucle sur les faces
    do m = b%md-hm, b%mf 
       do l = b%ld-hl, b%lf
          ! (l,m) : indice de l'interface
          ! (lg,mg) : indice de la cellule à gauche de l'interface
          ! (ld,md) : indice de la cellule à droite de l'interface
          lg = l; ld = l+hl
          mg = m; md = m+hm
          
          rhog = sum(rhoyCiR(:,lg,mg)) + rhoy2R(lg,mg)
          rhod = sum(rhoyCiL(:,ld,md)) + rhoy2L(ld,md)
          Un_g = uxR(lg,mg)*n(1,l,m) + uyR(lg,mg)*n(2,l,m)
          Un_d = uxL(ld,md)*n(1,l,m) + uyL(ld,md)*n(2,l,m)
          
!!! Cette valeur de kappa (~ max des Q) est en accord avec le livre de Toro p 324
          kappa = 0.d0
          kappa = max( abs(Un_g)+cR(lg,mg), abs(Un_d)+cL(ld,md) ) + v_mesh(l,m)
          
          do i = 1, b%ne
             j = irho1z1+i-1
             delta_U(j) = rhoyCiL(j,ld,md)-rhoyCiR(j,lg,mg)
          end do
          delta_U(irho2z2) = rhoy2L(ld,md)   - rhoy2R(lg,mg)
          delta_U(irhou)   = rhod*uxL(ld,md) - rhog*uxR(lg,mg)
          delta_U(irhov)   = rhod*uyL(ld,md) - rhog*uyR(lg,mg)
          delta_U(irhoE)   = rhod* EL(ld,md) - rhog* ER(lg,mg)
          
          flux(:,l,m) = flux(:,l,m) - kappa * delta_U(:)
          
!!! assemblage final
          flux(:,l,m)= 0.5d0 * flux(:,l,m) * dl(l,m)
       end do
    end do

  end subroutine solveur_HLL
  
  subroutine advection_z_ordre1(b,direction,solveur_Riemann,move,sm)
    type (STR_BLOC)                  , pointer       :: b
    integer                          , intent(in)    :: direction 
    integer                          , intent(in)    :: solveur_Riemann 
    logical                          , intent(in)    :: move
    real*8, dimension(1:,b%ld:,b%md:), intent(inout) :: sm

    integer :: l, m, hl, hm
    integer :: lg, ld, mg, md
    real*8 :: Un_g, Un_d, delta_p, Us, rhog, rhod
    real*8, dimension(1:b%neq) :: Ug, Ud
    real*8, dimension(:,:)  , pointer :: Cm, Cp, dl, v_mesh
    real*8, dimension(:,:,:), pointer :: n, fluxp, fluxm
    
    select case(direction)
    case (DIR_X) ! direction l
       Cm => b%Cm_l
       Cp => b%Cp_l
       n => b%n_dir_l
       dl => b%dl_l
       hl=1; hm=0
       if (move .eqv. .true.) then
          v_mesh => b%dX_dt
       else
          v_mesh => b%zeros
       end if
    case (DIR_Y) ! direction m
       Cm => b%Cm_m
       Cp => b%Cp_m
       n => b%n_dir_m
       dl => b%dl_m
       hl=0; hm=1
       if (move .eqv. .true.) then
          v_mesh => b%dY_dt
       else
          v_mesh => b%zeros
       end if
    end select
   
    fluxp => b%Fg_acoustique
    fluxp(1,:,:) = 0.d0 ! on stocke dans 1 mais c est juste une commodite
    fluxm => b%Fd_acoustique
    fluxm(1,:,:) = 0.d0

    ! boucle sur les faces
    do m = b%md-hm, b%mf 
       do l = b%ld-hl, b%lf
          
          lg = l; ld = l+hl
          mg = m; md = m+hm

          Ug(:) = b%U(:,lg,mg); Ud(:) = b%U(:,ld,md)
          call rho_from_U(Ug, rhog); call rho_from_U(Ud, rhod)
            
          Un_g = Ug(irhou)/rhog*n(1,l,m) + Ug(irhov)/rhog*n(2,l,m)
          Un_d = Ud(irhou)/rhod*n(1,l,m) + Ud(irhov)/rhod*n(2,l,m)
             
          if ( solveur_riemann == GALLICE ) then
             delta_p = b%p(ld,md) - b%p(lg,mg)
             Us = ( Cm(l,m)*Un_g + Cp(l,m)*Un_d - delta_p ) / ( Cm(l,m) + Cp(l,m) ) ! U*
          else
             Us = ( sqrt(rhog)*Un_g + sqrt(rhod)*Un_d ) / ( sqrt(rhog) + sqrt(rhod) ) ! moyenne de roe
          end if
          
          Us = Us + v_mesh(l,m)
          
          fluxp(1,l,m) = .5d0*(Us+abs(Us)) * (b%U(iz,ld,md) - b%U(iz,lg,mg))
          fluxm(1,l,m) = .5d0*(Us-abs(Us)) * (b%U(iz,ld,md) - b%U(iz,lg,mg))
       end do
    end do
    
    do m = b%md, b%mf
       do l = b%ld, b%lf
          sm(iz,l,m) = sm(iz,l,m) + (fluxm(1,l,m)*dl(l,m)+fluxp(1,l-1,m)*dl(l-1,m))
       end do
    end do
   
  end subroutine advection_z_ordre1
  
end module second_euler_direct
