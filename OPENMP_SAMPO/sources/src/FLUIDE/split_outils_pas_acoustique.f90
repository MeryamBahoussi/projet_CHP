module outils_pas_acoustique
  use m_struct
  use m_MPI
  use parametres_globaux, only : DIR_X, DIR_Y
  use parametres_fluide
  implicit none

contains

   ! Definition des pentes du solveur de Riemann
   ! si cdp_Riemann==.true., on calcule les pentes en tenant compte des conditions
   ! de positivite
   ! sinon on prend la vitesse du son Lagrangienne au centre de la face
   subroutine pentes_solveur_Riemann(sb)

     type (STR_SUPER_BLOC), pointer :: sb

     integer :: i, ib
     type (STR_BLOC), pointer :: b
     
     do i = 1, sb%nb_membre
        ib = sb%membre(i)
        b => acces_bloc(ib)%pointeur
        if (bloc_est_local (ib)) then

           select case (sb%solveur_Riemann)
           case(GALLICE)
              if (cdp_Riemann) then
                 call pentes_solveur_1d(b, DIR_X, sb%egalite_pentes_solveur, b%Cm_l, b%Cp_l)
                 call pentes_solveur_1d(b, DIR_Y, sb%egalite_pentes_solveur, b%Cm_m, b%Cp_m)
              else
                 call vitesse_du_son_Lagrangienne(b, b%a_l, b%a_m)
                 b%Cm_l(b%ld-1:b%lf, b%md:b%mf) = b%a_l(b%ld-1:b%lf, b%md:b%mf)
                 b%Cp_l(b%ld-1:b%lf, b%md:b%mf) = b%a_l(b%ld-1:b%lf, b%md:b%mf)
                 b%Cm_m(b%ld:b%lf, b%md-1:b%mf) = b%a_m(b%ld:b%lf, b%md-1:b%mf)
                 b%Cp_m(b%ld:b%lf, b%md-1:b%mf) = b%a_m(b%ld:b%lf, b%md-1:b%mf)
              end if

           case(HLL)
              !!--- pas vraiment teste
              !call cdp_pentes_HLL(b, DIR_X, b%Cm_l, b%Cp_l)
              !call cdp_pentes_HLL(b, DIR_Y, b%Cm_m, b%Cp_m)
              !b%a_l(b%ld-1:b%lf, b%md:b%mf) = b%Cm_l(b%ld-1:b%lf, b%md:b%mf)
              !b%a_m(b%ld:b%lf, b%md-1:b%mf) = b%Cm_m(b%ld:b%lf, b%md-1:b%mf)
              call vitesse_du_son_Lagrangienne(b, b%a_l, b%a_m)
              b%Cm_l(b%ld-1:b%lf, b%md:b%mf) = b%a_l(b%ld-1:b%lf, b%md:b%mf)
              b%Cp_l(b%ld-1:b%lf, b%md:b%mf) = b%a_l(b%ld-1:b%lf, b%md:b%mf)
              b%Cm_m(b%ld:b%lf, b%md-1:b%mf) = b%a_m(b%ld:b%lf, b%md-1:b%mf)
              b%Cp_m(b%ld:b%lf, b%md-1:b%mf) = b%a_m(b%ld:b%lf, b%md-1:b%mf)
           case(ROE)
              call vitesse_du_son_Lagrangienne(b, b%a_l, b%a_m)
              b%Cm_l(b%ld-1:b%lf, b%md:b%mf) = b%a_l(b%ld-1:b%lf, b%md:b%mf)
              b%Cp_l(b%ld-1:b%lf, b%md:b%mf) = b%a_l(b%ld-1:b%lf, b%md:b%mf)
              b%Cm_m(b%ld:b%lf, b%md-1:b%mf) = b%a_m(b%ld:b%lf, b%md-1:b%mf)
              b%Cp_m(b%ld:b%lf, b%md-1:b%mf) = b%a_m(b%ld:b%lf, b%md-1:b%mf)
           case default
              print*, "solveur de Riemann non defini"
              call arret_code
           end select
        end if
     end do

   end subroutine pentes_solveur_Riemann

   ! calcul de la vitesse du son Lagrangienne au centre des faces
   subroutine vitesse_du_son_Lagrangienne(b, a_l, a_m)
     type (STR_BLOC), pointer :: b
     real*8, dimension(:,:), pointer :: a_l, a_m  ! vitesse du son aux interfaces

     integer :: l, m
     integer :: ld, lf, md, mf
     real*8, dimension(:,:), pointer :: rho, c

     ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
     rho => b%rho; c => b%c

     do m = md-1, mf
        do l = ld-1, lf
           a_l(l,m) = b%K_Riemann*max(rho(l,m)*c(l,m), rho(l+1,m)*c(l+1,m))
           a_m(l,m) = b%K_Riemann*max(rho(l,m)*c(l,m), rho(l,m+1)*c(l,m+1))
        end do
     end do

   end subroutine vitesse_du_son_Lagrangienne

   subroutine cdp_pentes_HLL(b, direction, Cm, Cp)
     type (STR_BLOC), pointer :: b
     integer, intent(in) :: direction
     real*8, dimension(:,:), pointer :: Cm, Cp

     integer :: l, m
     integer :: lg, ld, mg, md, hl, hm
     ! cdp : contraintes de positivite
     real*8 :: cdp_vol, cdp_eps, cdp_c2g, cdp_c2d
     real*8 :: delta, delta_Un, eps_a, p_a, vol_a
     real*8 :: rho_g, p_g, c_g, un_g, eps_g, pi_g
     real*8 :: rho_d, p_d, c_d, un_d, eps_d, pi_d
     real*8, dimension(:,:,:), pointer :: n

     select case(direction)
     case (DIR_X)
        n => b%n_dir_l
        hl=1; hm=0
     case (DIR_Y)
        n => b%n_dir_m
        hl=0; hm=1
     end select

     do m = b%md-hm, b%mf
        do l = b%ld-hl, b%lf

           lg = l; ld = l+hl
           mg = m; md = m+hm

           ! on l'initialise a 0, au cas ou on les calcule pas
           cdp_vol=0.d0
           cdp_eps=0.d0
           cdp_c2g=0.d0
           cdp_c2d=0.d0

           ! grandeurs des deux cotes de l'interface
           rho_g = b%rho(lg,mg) ; rho_d = b%rho(ld,md)
           pi_g = b%pi_eq(lg,mg); pi_d = b%pi_eq(lg,mg)
           p_g = b%p(lg,mg)     ; p_d = b%p(ld,md)
           c_g = b%c(lg,mg)     ; c_d = b%c(ld,md)
           Un_g = b%u_x(lg,mg)*n(1,l,m) + b%u_y(lg,mg)*n(2,l,m)
           Un_d = b%u_x(ld,md)*n(1,l,m) + b%u_y(ld,md)*n(2,l,m)
           eps_g = b%E(lg,mg) - 0.5d0 * (b%u_x(lg,mg)**2 + b%u_y(lg,mg)**2)
           eps_d = b%E(ld,md) - 0.5d0 * (b%u_x(ld,md)**2 + b%u_y(ld,md)**2)

           delta_Un = Un_d-Un_g
           eps_a = .5d0*(eps_g+eps_d)
           p_a = .5d0*(p_g+p_d)
           vol_a = .5d0*(1.d0/rho_g+1.d0/rho_d)  ! volume specifique

!!! condition sur le volume specifique
           cdp_vol = -.5d0*delta_Un/vol_a

!!! condition sur l'energie interne
           delta = .25d0*(delta_Un*p_a/eps_a)**2.d0+.5d0*(p_d-p_g)**2.d0/eps_a
           cdp_eps = .25d0*delta_Un*p_a/eps_a + .5d0*sqrt(delta)

!!! condition sur le carre de la vitesse du son (a gauche)
           delta = .25d0*(delta_Un*(p_a+pi_g)/eps_a)**2.d0 + &
                4.d0/eps_a*(vol_a*pi_g + (p_d-p_g)**2.d0/8.d0)
           cdp_c2g = .5d0*(delta_Un*(p_a+pi_g)/(2.d0*eps_a) + sqrt(delta))

!!! condition sur le carre de la vitesse du son (a droite)
           delta = .25d0*(delta_Un*(p_a+pi_d)/eps_a)**2.d0 + &
                4.d0/eps_a*(vol_a*pi_d + (p_d-p_g)**2.d0/8.d0)
           cdp_c2d = .5d0*(delta_Un*(p_a+pi_d)/(2.d0*eps_a) + sqrt(delta))

           Cm(l,m) = max(c_g*rho_g, c_d*rho_d, cdp_vol, cdp_eps, cdp_c2g, cdp_c2d)

           Cm(l,m) = b%K_Riemann*Cm(l,m)
           Cp(l,m) = Cm(l,m)

        end do
     end do
   end subroutine cdp_pentes_HLL


!!! Calcul des pentes du solveur de Riemann suivant une direction
!!! en tenant compte des conditions de positivite
!!! Si egalite_pentes==.true., on a C+ = C-
!!! On ne regarde que la positivite du volume specifique et de p+pi
!!! En Stiffened Gas, cela implique aussi la positivite de l'energie interne
!!!
!!! TODO : ordre 2
   subroutine pentes_solveur_1d(b, direction, egalite_pentes, Cm, Cp)

     type (STR_BLOC), pointer :: b
     integer, intent(in) :: direction
     logical, intent(in) :: egalite_pentes
     real*8, dimension(:,:), pointer :: Cm, Cp

     integer :: l, m
     integer :: lg, ld, mg, md, hl, hm

     real*8, dimension(:,:,:), pointer :: n

     real*8  :: r, correction_pente, Cm_ini
     logical :: positivite_EI, positivite_EIg, positivite_EId

     ! Il faut avoir fait au prealable la mise a jour des mailles fictives EULER
     !r = 1  ! (=  C_p / C_m )

     select case(direction)
     case (DIR_X)
        n => b%n_dir_l
        hl=1; hm=0
     case (DIR_Y)
        n => b%n_dir_m
        hl=0; hm=1
     end select

!!! On met a 0 les pentes car on fait un max
     Cm(:,:) = 0.d0; Cp(:,:) = 0.d0

     ! cdp : contraintes de positivite
     call critere_positivite_volume(b, direction, egalite_pentes, Cm)

     select case (modele_diphasique)
     case(MONOPHASIQUE, ALLAIRE)
        call critere_positivite_c2(b, direction, egalite_pentes, Cm)
        do m = b%md-hm, b%mf
           do l = b%ld-hl, b%lf
              lg = l; ld = l+hl
              mg = m; md = m+hm
              r = b%rho(ld,md) * b%c(ld,md) / ( b%rho(lg,mg) * b%c(lg,mg) )
              if ( egalite_pentes ) r = 1.d0
              Cm(l,m) = b%K_Riemann * Cm(l,m)
              Cp(l,m) = r*Cm(l,m)
           end do
        end do
     case default
        print*, "pentes_solveur_1d : modele_diphasique inconnu"
        call arret_code
     end select

   end subroutine pentes_solveur_1d

  subroutine critere_positivite_volume(bloc, direction, egalite_pentes, Cm)
    type (STR_BLOC)       , pointer    :: bloc
    integer               , intent(in) :: direction
    logical               , intent(in) :: egalite_pentes
    real*8, dimension(:,:), pointer    :: Cm

    integer :: l, m
    integer :: lg, ld, mg, md, hl, hm
    real*8  :: a, b, c, discriminant, r
    real*8  :: cdp_rho_g, cdp_rho_d
    real*8  :: delta_U, delta_P
    real*8  :: p_g, c_g, un_g, vol_g, U_g(2)
    real*8  :: p_d, c_d, un_d, vol_d, U_d(2)
    real*8, dimension(:,:,:), pointer :: n

    select case(direction)
     case (DIR_X)
        n => bloc%n_dir_l
        hl=1; hm=0
     case (DIR_Y)
        n => bloc%n_dir_m
        hl=0; hm=1
     end select

     do m = bloc%md-hm, bloc%mf
        do l = bloc%ld-hl, bloc%lf

           lg = l; ld = l+hl
           mg = m; md = m+hm

           ! on l'initialise a 0, au cas ou on les calcule pas
           cdp_rho_g=0.d0; cdp_rho_d=0.d0

           ! grandeurs des deux cotes de l'interface
           p_g = bloc%p(lg,mg); p_d = bloc%p(ld,md)
           c_g = bloc%c(lg,mg); c_d = bloc%c(ld,md)
           vol_g = 1.d0/bloc%rho(lg,mg)
           vol_d = 1.d0/bloc%rho(ld,md)  ! volume specifique

           U_g=(/bloc%u_x(lg,mg), bloc%u_y(lg,mg)/)
           U_d=(/bloc%u_x(ld,md), bloc%u_y(ld,md)/)
           Un_g = dot_product(U_g, n(:,l,m))
           Un_d = dot_product(U_d, n(:,l,m))

           r = c_d/vol_d / ( c_g/vol_g )
           if ( egalite_pentes ) r = 1.d0

           delta_u = un_d - un_g
           delta_p = p_d - p_g

!!! Condition pour que le volume specifique de l'EI* gauche soit positif
           a = (1.d0+r) * vol_g
           b = r * delta_U
           c = -delta_p
           discriminant = b**2.d0 - 4.d0*a*c
           if (discriminant >= 0.d0) cdp_rho_g = (-b + dsqrt(discriminant))/(2.d0*a)

!!! Condition pour que le volume specifique de l'EI* droit soit positif
           a = r * (1.d0+r) * vol_d
           b = delta_U
           c = delta_p
           discriminant = b**2.d0 - 4.d0*a*c
           if (discriminant >= 0.d0) cdp_rho_d = (-b + dsqrt(discriminant))/(2.d0*a)

!!! Remplissage des tableaux avec les pentes qui respectent tous les criteres de positivite
           Cm(l,m) = max( c_g/vol_g, cdp_rho_g, cdp_rho_d )
           if ( egalite_pentes ) Cm(l,m) = max( Cm(l,m), c_d/vol_d )

        end do
     end do

  end subroutine critere_positivite_volume

  subroutine critere_positivite_c2(bloc, direction, egalite_pentes, Cm)
    type (STR_BLOC)       , pointer    :: bloc
    integer               , intent(in) :: direction
    logical               , intent(in) :: egalite_pentes
    real*8, dimension(:,:), pointer    :: Cm

    integer :: l, m
    integer :: lg, ld, mg, md, hl, hm
    real*8  :: a, b, c, discriminant, r
    real*8  :: cdp_c2_g, cdp_c2_d
    real*8  :: delta_U, delta_P
    real*8  :: p_g, c_g, un_g, vol_g, U_g(2), eps_g, pi_g, epshat_g, GPi_g
    real*8  :: p_d, c_d, un_d, vol_d, U_d(2), eps_d, pi_d, epshat_d, GPi_d
    real*8, dimension(:,:,:), pointer :: n

    select case(direction)
     case (DIR_X)
        n => bloc%n_dir_l
        hl=1; hm=0
     case (DIR_Y)
        n => bloc%n_dir_m
        hl=0; hm=1
     end select

    do m = bloc%md-hm, bloc%mf
       do l = bloc%ld-hl, bloc%lf
          lg = l; ld = l+hl
          mg = m; md = m+hm

          ! on l'initialise a 0, au cas ou on les calcule pas
          cdp_c2_g=0.d0; cdp_c2_d=0.d0

          ! grandeurs des deux cotes de l'interface
          p_g = bloc%p(lg,mg); p_d = bloc%p(ld,md)
          c_g = bloc%c(lg,mg); c_d = bloc%c(ld,md)
          vol_g = 1.d0/bloc%rho(lg,mg)
          vol_d = 1.d0/bloc%rho(ld,md)  ! volume specifique

          U_g=(/bloc%u_x(lg,mg), bloc%u_y(lg,mg)/)
          U_d=(/bloc%u_x(ld,md), bloc%u_y(ld,md)/)
          Un_g = dot_product(U_g, n(:,l,m))
          Un_d = dot_product(U_d, n(:,l,m))
          delta_u = un_d - un_g
          delta_p = p_d - p_g

          eps_g = bloc%E(lg,mg) - 0.5d0 * (U_g(1)**2 + U_g(2)**2)
          eps_d = bloc%E(ld,md) - 0.5d0 * (U_d(1)**2 + U_d(2)**2)
          pi_g = bloc%pi_eq(lg,mg)
          pi_d = bloc%pi_eq(ld,md)
          epshat_g = eps_g - pi_g*vol_g
          epshat_d = eps_d - pi_d*vol_d

          r = (c_d/vol_d) / ( c_g/vol_g )
          if ( egalite_pentes ) r = 1.d0

          GPi_g = 2.d0 * ( r*P_g + P_d ) - delta_p
          GPi_d = 2.d0 * ( r*P_g + P_d ) + r * delta_p

!!! Condition pour l'hyperbolicité de l'etat intermediaire * gauche
          a = 2.d0 * (1.d0+r)**2.d0 * epshat_g + r**2.d0 * delta_u**2.d0
          b = - r * delta_U * ( delta_p + GPi_g + 2.d0 * (1.d0+r) * pi_g )
          c = delta_p * ( GPi_g + 2.d0 * (1.d0+r) * pi_g )
          discriminant = b**2.0d0 - 4.d0*a*c
          if (discriminant >= 0.d0) cdp_c2_g = (-b + dsqrt(discriminant))/(2.d0*a)

!!! Condition pour l'hyperbolicité de l'etat intermediaire * droite
          a = 2.d0 * r * (1.d0+r)**2.d0 * epshat_d + r * delta_u**2.d0
          b = - delta_U * ( GPi_d - r * delta_p  + 2.d0 * (1.d0+r) * pi_d )
          c = - delta_p * ( GPi_d + 2.d0 * (1.d0+r) * pi_d )
          discriminant = b**2.0d0 - 4.d0 * a * c
          if (discriminant >= 0.d0) cdp_c2_d = (-b + dsqrt(discriminant))/(2.d0*a)

          !!! Remplissage des tableaux avec les pentes qui respectent tous les criteres de positivite
           Cm(l,m) = max( Cm(l,m), cdp_c2_g, cdp_c2_d )
       end do
    end do

  end subroutine critere_positivite_c2

  subroutine critere_z(bloc, direction, egalite_pentes, Cm)
    type (STR_BLOC)       , pointer    :: bloc
    integer               , intent(in) :: direction
    logical               , intent(in) :: egalite_pentes
    real*8, dimension(:,:), pointer    :: Cm

    integer :: l, m
    integer :: lg, ld, mg, md, hl, hm
    real*8  :: a, b, c, discriminant, r
    real*8  :: cdp_zg, cdp_1mzg, cdp_zd, cdp_1mzd
    real*8  :: delta_U, delta_P
    real*8  :: z_g, c_g, vol_g, U_g(2), rho1c1_2
    real*8  :: z_d, c_d, vol_d, U_d(2), rho2c2_2
    real*8, dimension(:,:,:), pointer :: n

    select case(direction)
     case (DIR_X)
        n => bloc%n_dir_l
        hl=1; hm=0
     case (DIR_Y)
        n => bloc%n_dir_m
        hl=0; hm=1
     end select

    do m = bloc%md-hm, bloc%mf
       do l = bloc%ld-hl, bloc%lf

          lg = l; ld = l+hl
          mg = m; md = m+hm

          ! on l'initialise a 0, au cas ou on les calcule pas
          cdp_zg=0.d0; cdp_1mzg=0.d0
          cdp_zd=0.d0; cdp_1mzd=0.d0

          ! grandeurs des deux cotes de l'interface
          z_g = bloc%z(lg,mg); z_d = bloc%z(ld,md)
          c_g = bloc%c(lg,mg); c_d = bloc%c(ld,md)
          vol_g = 1.d0/bloc%rho(lg,mg)
          vol_d = 1.d0/bloc%rho(ld,md)  ! volume specifique

          U_g=(/bloc%u_x(lg,mg), bloc%u_y(lg,mg)/)
          U_d=(/bloc%u_x(ld,md), bloc%u_y(ld,md)/)

          delta_u = dot_product(U_d, n(:,l,m)) - dot_product(U_g, n(:,l,m))
          delta_p = bloc%p(ld,md) - bloc%p(lg,mg)

          r = c_d/vol_d / ( c_g/vol_g )
          if ( egalite_pentes ) r = 1.d0

          rho1c1_2 = fluide1%EOS%gamma*(bloc%p(lg,mg) + fluide1%EOS%pi)
          rho2c2_2 = fluide2%EOS%gamma*(bloc%p(lg,mg) + fluide2%EOS%pi)

!!! Condition pour que la fraction volumique de l'EI* gauche soit entre 0 et 1
          a = (1.d0+r) * vol_g * (z_g + (1.d0-z_g) * rho1c1_2/rho2c2_2 )
          b = r * delta_U
          c = -delta_p
          discriminant = b**2.0d0 - 4.d0*a*c
          if (discriminant >= 0.d0) cdp_zg = (-b + dsqrt(discriminant))/(2.d0*a)

          a = (1.d0+r) * vol_g * (z_g * rho2c2_2/rho1c1_2 + 1.d0 - z_g)
          discriminant = b**2.0d0 - 4.d0*a*c
          if (discriminant >= 0.d0) cdp_1mzg = (-b + dsqrt(discriminant))/(2.d0*a)  !!! 1-z>0

!!! Condition pour que la fraction volumique de l'EI* droit soit entre 0 et 1
          a = r*(1.d0+r) * vol_d * (z_d + (1.d0-z_d)*rho1c1_2/rho2c2_2)
          b = delta_U
          c = delta_p
          discriminant = b**2.0d0 - 4.d0*a*c
          if (discriminant >= 0.d0) cdp_zd = (-b + dsqrt(discriminant))/(2.d0*a)

          a = r*(1.d0+r) * vol_d * (z_d * rho2c2_2/rho1c1_2 + 1.d0-z_d)
          discriminant = b**2.0d0 - 4.d0*a*c
          if (discriminant >= 0.d0) cdp_1mzd = (-b + dsqrt(discriminant))/(2.d0*a)  !!! 1-z>0

!!! Remplissage des tableaux avec les pentes qui respectent tous les criteres de positivite
          Cm(l,m) = max( Cm(l,m), cdp_zg, cdp_1mzg, cdp_zd, cdp_1mzd )
       end do
    end do
  end subroutine critere_z

end module outils_pas_acoustique
