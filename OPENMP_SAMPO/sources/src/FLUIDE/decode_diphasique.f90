module decodage_diphasique
  use EquationOfState
  use parametres_fluide

  implicit none

  integer, parameter :: choix_T = 1
!!! choix_T=1 : temperature de melange
!!! choix_T=2 : temperature liee a l'entropie (ne pas utiliser)

contains

  subroutine decode_diphasique(ldeb, lfin, mdeb, mfin, b, ierreur)
    integer        , intent(in)  :: ldeb, lfin, mdeb, mfin
    type (STR_BLOC), pointer     :: b
    integer        , intent(out) :: ierreur

    integer :: l, m
    real*8  :: z1,z2, eps
    type(STR_fluide) :: melange
    real*8 :: rho, rho1, rho2, c1, c2
    integer :: ierreur_c
    if (b%ne > 1) then
       print*, "diphasique + multiespeces non gere"
       call arret_code
    end if

    b%c_esp(1,ldeb:lfin,mdeb:mfin) = 1.d0
    b%c_ele(1,ldeb:lfin,mdeb:mfin) = 1.d0

    ! fractions volumique et massique
    ! on assure que ces deux quantitees sont bien entre 0 et 1
    ! le schema AD pour le transport donne parfois des z<0 et y<0 !??!
    do m = mdeb, mfin
       do l = ldeb, lfin
          if(b%U(iz,l,m)<=tol_z.or.b%U(irho1z1,l,m)<0.d0)then
            b%U(iz,l,m)      = 0.d0
            b%U(irho1z1,l,m) = 0.d0
            b%U(irho2z2,l,m) = b%rho(l,m)

            b%z(l,m)    = 0.d0
            b%y(l,m)    = 0.d0
            b%rho1(l,m) = 0.d0
            b%rho2(l,m) = b%rho(l,m)
          else if(b%U(iz,l,m)>=(1.d0-tol_z).or.b%U(irho2z2,l,m)<0.d0)then
            b%U(iz,l,m)      = 1.d0
            b%U(irho1z1,l,m) = b%rho(l,m)
            b%U(irho2z2,l,m) = 0.d0

            b%z(l,m)    = 1.d0
            b%y(l,m)    = 1.d0
            b%rho1(l,m) = b%rho(l,m)
            b%rho2(l,m) = 0.d0
          else 
            b%z(l,m)    = b%U(iz,l,m)
            b%y(l,m)    = b%U(irho1z1,l,m)/b%rho(l,m)
            b%rho1(l,m) = b%U(irho1z1,l,m)/b%z(l,m)
            b%rho2(l,m) = b%U(irho2z2,l,m)/(1.d0-b%z(l,m))
          end if

          if(b%rho1(l,m)<0.d0.or.b%rho2(l,m)<0.d0.or.b%y(l,m)<0.d0)then
            print*, 'pb décodage : y ou rho1 ou rho2 <0'
            print*, b%y(l,m), b%rho1(l,m), b%rho2(l,m)
            ierreur = 1
          end if
       end do
    end do

    ! on determine le type de fermeture du mélange en fonction de chaque fluide
    z1 = minval(b%z(ldeb:lfin,mdeb:mfin)) ! z_min
    z2 = maxval(b%z(ldeb:lfin,mdeb:mfin)) ! z_max
    call type_melange(fluide1, fluide2, z1, z2, melange)

    select case(melange%type)
    case(EOS)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! calcul de Cv, Cp et des coefficients de transport K, mu, D
!!! puis calcul des parametres de la loi de melange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       b%cv1(ldeb:lfin,mdeb:mfin) = fluide1%EOS%cv
       b%gamma1(ldeb:lfin,mdeb:mfin) = fluide1%EOS%gamma
       do m = mdeb, mfin
          do l = ldeb, lfin
             call grandeurs_melange(fluide1, fluide2, b%z(l,m), b%y(l,m), melange)
             call stocke_grandeurs_melange(melange, b, l, m)
          end do
       end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! decodage de la pression du melange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       select case (b%eq_energie)
       case (0) ! cas normal ou on a une equation d'energie
          do m = mdeb, mfin
             do l = ldeb, lfin
                melange%EOS%gamma = b%gamma_eq(l,m)
                melange%EOS%pi = b%pi_eq(l,m)
                melange%EOS%q = b%q_eq(l,m)
                eps = b%E(l,m) - 0.5d0*(b%u_x(l,m)**2+b%u_y(l,m)**2)
                b%p(l,m) = EOS_p_from_rho_eps(melange%EOS, b%rho(l,m), eps)
             end do
          end do
       case (1) ! isentropique
          do m = mdeb, mfin
             do l = ldeb, lfin
                b%p(l,m) = b%entropie_ref*b%rho(l,m)**melange%EOS%gamma - melange%EOS%pi
             end do
          end do
       case (2) ! isotherme
          print*, "isotherme non gere"
          call arret_code
       case default
          print*, "decode_diphasique : non gere"
          call arret_code
       end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! vitesse du son du melange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       select case (MODELE_DIPHASIQUE)
       case(ALLAIRE)
          if(b%correction_son==0)then
             do m = mdeb, mfin
                do l = ldeb, lfin
                   call vitesse_du_son_melange_Allaire(b%gamma_eq(l,m), b%pi_eq(l,m), b%p(l,m), b%rho(l,m), b%z(l,m), b%c(l,m),ierreur)
                   if(b%p(l,m)+b%pi_eq(l,m)<0.d0)then
                     print*, 'pb son  ', numproc, l, m
                     print*, 'c & z   ', b%c(l,m), b%z(l,m)
                     print*, 'p & pi  ', b%p(l,m), b%pi_eq(l,m)
                     print*, 'E & kin ', b%E(l,m), 0.5d0*(b%u_x(l,m)**2+b%u_y(l,m)**2)
                   end if
                enddo
             enddo
          else if(b%correction_son==1)then
             do m = mdeb, mfin
                do l = ldeb, lfin
                   ierreur_c = 0
                   call vitesse_du_son_melange_Allaire(b%gamma_eq(l,m), b%pi_eq(l,m), b%p(l,m), b%rho(l,m), b%z(l,m), b%c(l,m),ierreur_c)
                   !! solution 1 sur les pb de vitesse du son négative
                   if(ierreur_c==1)then

                      b%debug_c(l,m) = b%debug_c(l,m) + 1
                      rho1 = b%rho(l,m)*      b%y(l,m) /b%z(l,m)
                      rho2 = b%rho(l,m)*(1.d0-b%y(l,m))/(1.d0-b%z(l,m))

                      call vitesse_du_son_melange_Allaire(fluide2%EOS%gamma, fluide2%EOS%pi, b%p(l,m), rho2, 0.d0, c2, ierreur_c)
                      call vitesse_du_son_melange_Allaire(fluide1%EOS%gamma, fluide1%EOS%pi, b%p(l,m), rho1, 1.d0, c1, ierreur_c)

                      b%c(l,m) = max(c1,c2)
                      if(b%c(l,m)<=0.d0)then
                        print*, 'correction du son négatif', b%c(l,m)
                        ierreur   = 1
                      else if(b%c(l,m)>0.d0)then
                        ! rien
                      else
                        ! si on a des nan
                        print*, 'correction nan'
                        print*, b%p(l,m), b%c(l,m)
                        ierreur   = 1
                      end if
                   end if
                enddo
             enddo
          else if(b%correction_son==2)then
             do m = mdeb, mfin
                do l = ldeb, lfin
                   ierreur_c = 0
                   call vitesse_du_son_melange_Allaire(b%gamma_eq(l,m), b%pi_eq(l,m), b%p(l,m), b%rho(l,m), b%z(l,m), b%c(l,m),ierreur_c)
                   if(b%p(l,m)<0.d0)then
                      b%debug_c(l,m) = b%debug_c(l,m) + 1
                      rho1 = b%rho(l,m)*      b%y(l,m) /b%z(l,m)
                      rho2 = b%rho(l,m)*(1.d0-b%y(l,m))/(1.d0-b%z(l,m))

                      call vitesse_du_son_melange_Allaire(fluide2%EOS%gamma, fluide2%EOS%pi, b%p(l,m), rho2, 0.d0, c2, ierreur_c)
                      call vitesse_du_son_melange_Allaire(fluide1%EOS%gamma, fluide1%EOS%pi, b%p(l,m), rho1, 1.d0, c1, ierreur_c)

                      b%c(l,m) = max(c1,c2)
                      ierreur_c = 0
                      if(b%c(l,m)<=0.d0)then
                        print*, 'correction du son négatif', b%c(l,m)
                        ierreur_c = 0
                        ierreur   = 1
                      else
                        ierreur_c = 0
                      end if
                   end if
                end do
             end do
          end if
       end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Temperature du melange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (b%visqueux) then
       b%T1(ldeb:lfin,mdeb:mfin) = 0.d0
       b%T2(ldeb:lfin,mdeb:mfin) = 0.d0
       b%T (ldeb:lfin,mdeb:mfin) = 0.d0

       do m = mdeb, mfin
          do l = ldeb, lfin
             if (b%rho1(l,m) /=0.d0) then
                b%T1(l,m) = EOS_T_from_p_rho(fluide1%EOS, b%p(l,m), b%rho1(l,m))
             end if
          end do
       end do

       do m = mdeb, mfin
          do l = ldeb, lfin
             if (b%rho2(l,m) /=0.d0) then
                b%T2(l,m) = EOS_T_from_p_rho(fluide2%EOS, b%p(l,m), b%rho2(l,m))
             end if
          end do
       end do

       do m = mdeb, mfin
          do l = ldeb, lfin
             call T_melange(b%Cv1(l,m), b%T1(l,m), fluide2%EOS%Cv, b%T2(l,m), &
                  b%y(l,m), b%Cv(l,m), b%T(l,m), ierreur)
          end do
       end do
    end if
    case(CHIMIE_FIGEE)
       print*, "decode_diphasique : chimie dans le fluide 1 non gere"
       call arret_code
    case default
       print*, "decode_diphasique : type de fluide inconnu"
       call arret_code
    end select

  end subroutine decode_diphasique

  subroutine decode_diphasique_local(U,cell)
    real*8, dimension(:), intent(in)    :: U
    type (STR_cell)     , intent(inout) :: cell

    integer :: ierreur
    real*8  :: z1, z2, T1, T2, eps

    if (fluide1%ne > 1) then
       print*, "diphasique (local) + multiespeces non gere"
       call arret_code
    end if

    ! fraction volumique
    cell%z = max( min(U(iz), 1.d0-epsilon_z), 0.d0+epsilon_z )

    ! fraction massique des especes du melange gazeux
    cell%ci(:) = 1.d0
    if (U(irho1z1) > tol_z) then
       cell%y = U(irho1z1) / cell%rho      ! y = rho*y/rho
    else
       print*, "rho1z1 negatif", U(irho1z1)
       call arret_code
    end if
    cell%y = max( min(cell%y, 1.d0-epsilon_z), 0.d0+epsilon_z )

    z1 = cell%z; z2 = 1.d0 - z1
    call type_melange(fluide1, fluide2, z1, z2, cell%fluide)

!!! densite de chaque fluide
    cell%rho1 = 0.d0; cell%rho2=0.d0
    if (abs(z1) > tol_z) cell%rho1 = U(irho1z1)/z1
    if (abs(z2) > tol_z) cell%rho2 = U(irho2z2)/z2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! calcul des Cv, Cp et des coefficients de transport K, mu, D
!!! puis calcul des parametres de la loi de melange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! on determine gamma, pi et q, mu, k, cv du melange
    select case(fluide1%type)
    case(EOS)
       call grandeurs_melange(fluide1, fluide2, cell%z, cell%y, cell%fluide)
    case(CHIMIE_FIGEE)
       print*, "decode_diphasique_local : chimie dans le fluide 1 non gere"
       call arret_code
    case default
       print*, "decode_diphasique_local : type de fluide inconnu"
       call arret_code
    end select


    select case(cell%fluide%type)
    case(EOS)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! decodage de la pression du melange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       eps = U(irhoE)/cell%rho - 0.5d0*(cell%u**2+cell%v**2)
       cell%p = EOS_p_from_rho_eps(cell%fluide%EOS, cell%rho, eps)

!!! vitesse du son
       select case(MODELE_DIPHASIQUE)
       case(ALLAIRE)
          call vitesse_du_son_melange_Allaire(cell%fluide%EOS%gamma, cell%fluide%EOS%pi, cell%p, cell%rho, cell%z, cell%c, ierreur)
       end select

!!! temperature du melange
       T1 = 0.d0
       if (cell%rho1 /=0.d0) then
          T1 = EOS_T_from_p_rho(fluide1%EOS, cell%p, cell%rho1)
       end if

       T2 = 0.d0
       if (cell%rho2 /=0.d0) then
          T2 = EOS_T_from_p_rho(fluide2%EOS, cell%p, cell%rho2)
       end if

       cell%T = 0.d0
       call T_melange(fluide1%EOS%Cv, T1, fluide2%EOS%Cv, T2, &
            cell%y, cell%fluide%EOS%Cv, cell%T, ierreur)
    case default
       print*, "decode_diphasique_local : non gere"
       call arret_code
    end select

  end subroutine decode_diphasique_local


!!!! temperature de melange
!!! T = yi Cvi Ti,   mais on a pas epsilon = Cv T + pi/rho
!!! on a epsilon = Cv T + f(rho)
  subroutine T_melange(Cv1, T1, Cv2, T2, y, Cv, T, ierreur)
    real*8          , intent(in)  :: Cv1, T1
    real*8          , intent(in)  :: Cv2, T2
    real*8          , intent(in)  :: Cv, y
    real*8          , intent(out) :: T
    integer         , intent(out) :: ierreur

    real*8 :: y1, y2

    y1 = y; y2 = 1.d0 - y

    T = (y1*Cv1*T1 + y2*Cv2*T2) / Cv

    if ( T < 0.d0 ) then
       print*, "temperature negative", T
       print*, "T1", real(T1), "y1", real(y1)
       print*, "T2", real(T2), "y2", real(y2)
       ierreur = 1
    end if

  end subroutine T_melange

!!! Vitesse du son du melange
!!! Modele d'Allaire et al. (ou Massoni et al.)
!!! c^2 = gamma(z)*(p+pi(z))/rho
!!!
  subroutine vitesse_du_son_melange_Allaire(gamma, pi, p, rho, z, c, ierreur)

    real*8          , intent(in)  :: gamma, pi, p, rho, z
    real*8          , intent(out) :: c
    integer         , intent(out) :: ierreur

    real*8 :: c2

    c2 = gamma*(p + pi)/rho

    if (c2 <= 0.d0) then
       ! print*, "vitesse du son negative : ", c2
       ! print*, "p  ", p  , "z  ", z
       ! print*, "rho", rho, "pi ", pi
       c = c2
       ierreur = 1
    else
       c = sqrt(c2) ! vitesse du son
    end if

  end subroutine vitesse_du_son_melange_Allaire

  ! On determine le type du mélange et si le code le supporte
  subroutine type_melange(f1, f2, z_min, z_max, melange)
    type(STR_fluide), intent(in)  :: f1, f2
    real*8          , intent(in)  :: z_min, z_max
    type(STR_fluide), intent(out) :: melange

    ! on regarde si les deux fluides ont le même type de fermeture (EOS, CHIMIE_FIGEE, ...)
    if (f1%type == f2%type) then
       melange%type = f1%type
    else
       if (z_min==z_max) then
          if (z_min == 1.d0) then
             melange%type = f1%type
          else if (z_min == 0.d0) then
             melange%type = f2%type
          else
             melange%type = -1 ! pas gere par le code
          end if
       else
          melange%type = -1 ! pas gere par le code
       end if
    end if

  end subroutine type_melange

  subroutine grandeurs_melange(f1, f2, z, y, melange)
    type(STR_fluide), intent(in)  :: f1, f2
    real*8          , intent(in)  :: z, y
    type(STR_fluide), intent(out) :: melange

    real*8 :: xi

    xi = z/(f1%EOS%gamma-1.d0) + (1.d0-z)/(f2%EOS%gamma-1.d0)
    melange%EOS%gamma = 1.d0 + 1.d0 / xi

    melange%EOS%pi = (melange%EOS%gamma-1.d0)/melange%EOS%gamma * &
         (  z    *f1%EOS%pi*f1%EOS%gamma/(f1%EOS%gamma-1.d0) + &
         (1.d0-z)*f2%EOS%pi*f2%EOS%gamma/(f2%EOS%gamma-1.d0) )

!!! q et cv
    melange%EOS%q  = y*f1%EOS%q  + (1.d0-y)*f2%EOS%q  ! ok coherent avec eps = y_i eps_i
    melange%EOS%cv = y*f1%EOS%cv + (1.d0-y)*f2%EOS%cv

  end subroutine grandeurs_melange

  subroutine stocke_grandeurs_melange(melange, b, l, m)
    type(STR_fluide), intent(in) :: melange
    type (STR_BLOC) , pointer    :: b
    integer         , intent(in) :: l,m

    b%gamma_eq(l,m) = melange%EOS%gamma
    b%pi_eq(l,m) = melange%EOS%pi
    b%q_eq(l,m)  = melange%EOS%q
    b%Cv(l,m) = melange%EOS%Cv

  end subroutine stocke_grandeurs_melange



!!!! temperature de melange (2 options)
!!! Choix 1 :
!!! T = yi Cvi Ti,   mais on a pas epsilon = Cv T + pi/rho
!!! on a epsilon = Cv T + f(rho)
!!!
!!! Choix 2 :
!!! Cv T = (p+pi)/(rho*(gamma-1)), dans ce cas on a epsilon = Cv T + pi/rho
!!! ne marche pas avec la loi d'état de type Mie Gruneisen
  subroutine T_melange_from_p_rho_z_y_old(melange, p, rho, z, y, T, T1, T2, ierreur)

    type(STR_fluide), intent(in)  :: melange
    real*8          , intent(in)  :: p, rho, z, y
    real*8          , intent(out) :: T, T1, T2
    integer         , intent(out) :: ierreur

    real*8 :: Cv, gamma, pi
    real*8 :: z1, y1, rho1
    real*8 :: z2, y2, rho2

    Cv = melange%EOS%Cv
    gamma = melange%EOS%gamma
    pi = melange%EOS%pi

    z1 = z; z2 = 1.d0 - z
    y1 = y; y2 = 1.d0 - y

    ! calcul de la temperature des phases pures
    T1   = 0.d0;   T2 = 0.d0
    rho1 = 0.d0; rho2 = 0.d0
    if (choix_T==1) then
       if (abs(z1) > tol_z .and. abs(y1) > tol_z) then
          rho1 = rho*y1/z1
          T1 = EOS_T_from_p_rho(fluide1%EOS, p, rho1)
       end if
       if (abs(z2) > tol_z .and. abs(y2) > tol_z) then
          rho2 = rho*y2/z2
          T2 = EOS_T_from_p_rho(fluide2%EOS, p, rho2)
       end if
    end if

    if (choix_T==1) then
       T = (y1*fluide1%EOS%Cv*T1 + y2*fluide2%EOS%Cv*T2) / Cv
    else if (choix_T==2) then
       T = (y1*fluide1%EOS%T0 + y2*fluide2%EOS%T0) + (p+pi)/(rho*Cv*(gamma-1.d0))
    end if

    if ( T < 0.d0 ) then
       print*, "temperature negative", T
       print*, "pression", p
       print*, "densite", rho
       print*, "T1", real(T1), "z1", real(z1), "y1", real(y1), "rho1", real(rho1)
       print*, "T2", real(T2), "z2", real(z2), "y2", real(y2), "rho2", real(rho2)
       ierreur = 1
    end if

  end subroutine T_melange_from_p_rho_z_y_old

end module decodage_diphasique
