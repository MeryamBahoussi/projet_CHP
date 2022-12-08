module decodage
  use m_MPI
  use m_IO
  use EquationOfState
  use Multiespeces
  use parametres_fluide
  use outils_fluide

  use decodage_diphasique
  use decodage_chimie
  use m_transport_diffusion
  use thermo_fluide

  use OMP_LIB

  implicit none

contains

  subroutine encode_U(cell, U)
    type (STR_cell), intent(in)  :: cell
    real*8         , intent(out) :: U(:)

    integer :: i, ne
    real*8  :: rho, z

    ne = fluide1%ne  ! soit 1, nb_esp ou nb_ele

    z = max( min(cell%z, 1.d0-epsilon_z), 0.d0+epsilon_z )
    rho = z*cell%rho1 + (1.d0-z)*cell%rho2

    do i = 1, ne
       U(irho1z1+i-1) = cell%rho1 * z * cell%Ci(i)           ! rho1*z1*Ci
    end do
    U(irho2z2) = cell%rho2 * (1.d0-z)                        ! rho2*(1-z)
    U(irhou) = rho * cell%u                                  ! rho*u
    U(irhov) = rho * cell%v                                  ! rho*v
    U(irhoE) = rho * cell%h - cell%p + .5d0*rho*(cell%u**2+cell%v**2) ! rho*e
    U(iz   ) = z
  end subroutine encode_U

  subroutine encode_Ub(b, ldeb, lfin, mdeb, mfin)
    type (STR_bloc), pointer     :: b
    integer        , intent(in)  :: ldeb, lfin, mdeb, mfin

    integer         :: i, ne, l, m
    real*8, pointer :: ci(:,:,:)

    ne = fluide1%ne

    select case (fluide1%type)
    case(EOS, CHIMIE_FIGEE)
       ci => b%c_esp
    end select

    do m = mdeb, mfin
       do l = ldeb, lfin
          do i = 1, ne
             b%U(irho1z1+i-1,l,m) = b%rho1(l,m) * b%z(l,m) * ci(i,l,m)
          end do
          b%U(irho2z2,l,m) = b%rho2(l,m) * (1.d0 - b%z(l,m))    ! rho2*(1-z)
          b%U(irhou,l,m) = b%rho(l,m) * b%u_x(l,m)              ! rho*u
          b%U(irhov,l,m) = b%rho(l,m) * b%u_y(l,m)              ! rho*v
          b%U(irhoE,l,m) = b%rho(l,m) * b%E(l,m)                ! rho*E
          b%U(iz,l,m) = b%z(l,m)
       end do
    end do
  end subroutine encode_Ub

  subroutine decode(sb, decode_MF)

    type (STR_SUPER_BLOC), pointer    :: sb
    logical              , intent(in) :: decode_MF

    type (STR_BLOC), pointer :: b
    integer :: i, ib
    integer :: nMF_loc ! nombre de rangees de mailles fictives

    integer :: err_mpi
    integer :: ierreur, ierreur_MPI
!!! si ierreur=1, c'est que le decodage s'est mal passe
!!! on fait les sorties visit et on arrete le programme

    call debut_watchTime(wT_decodage)

    nMF_loc = 0
    if (decode_MF) nMF_loc = nMF

    ierreur = 0

    do i = 1, sb%nb_membre  ! boucle sur les blocs
       ib = sb%membre(i)
       b => acces_bloc(ib)%pointeur
       if (bloc_est_local (ib)) then
          call decode_U(b%ld-nMF_loc, b%lf+nMF_loc, b%md-nMF_loc, b%mf+nMF_loc, b, ierreur)
       end if
    end do

    ierreur_MPI = 0
    call MPI_ALLREDUCE(ierreur,ierreur_mpi,1,mpi_int,mpi_MAX,sb%mpi_comm,err_mpi)

    ! si il y a eu une erreur dans le decodage
    if (ierreur_MPI==1) then
       do i = 1, sb%nb_membre  ! boucle sur les blocs
          ib = sb%membre(i)
          b => acces_bloc(ib)%pointeur
          if (bloc_est_local (ib)) then
             call sorties_visit_nblocs(b, -2, it_out)
          end if
       end do
       call MPI_BARRIER(sb%mpi_comm, err_mpi)
       call arret_code
    end if

  call fin_watchtime(wT_decodage)

  end subroutine decode

  subroutine decode_U2cell(U, cell)
    use variables_globales
    use chimie

    real*8, dimension(:), intent(in)  :: U
    type (STR_cell)     , intent(out) :: cell

    integer :: ierreur, i, ne
    real*8  :: somme, eps
    real*8  :: entropie, hkapa, hk1, cp, cv

    ne = fluide1%ne

    call rho_from_U(U, cell%rho) ! densite

    if (cell%rho < 0.d0) then
       print*, "densite negative", cell%rho
       call arret_code
    end if

    cell%u = U(irhou)/cell%rho
    cell%v = U(irhov)/cell%rho

    ! fraction massique des especes du melange gazeux
    cell%ci(:) = 0.d0
    if (U(irho1z1) > tol_z) then
       somme = 0.d0
       do i = 1, ne
          cell%ci(i) = U(irho1z1+i-1)   ! ci = rho*y*ci
          cell%ci(i) = max(cell%ci(i),0.d0)
          somme = somme + cell%ci(i)
       end do
       cell%ci(:) = cell%ci(:)/somme  ! ci = rho*y*ci/(rho*y)
       cell%y = somme / cell%rho      ! y = rho*y/rho
    else
       print*, "rho1z1 negatif", U(irho1z1)
       ierreur = 1
    end if

    select case(MODELE_DIPHASIQUE)
    case (MONOPHASIQUE)
       cell%z = 1.d0
       cell%y = 1.d0
       cell%rho1 = cell%rho
       cell%rho2=0.d0

       select case(fluide1%type)
       case(EOS)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! calcul des Cv, Cp et des coefficients de transport K, mu, D
!!! puis calcul des parametres de la loi de melange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! on determine gamma, pi et q, mu, k, cv du melange
          call grandeurs_melange(fluide1, fluide2, cell%z, cell%y, cell%fluide)
!!! decodage de la pression
          ! epsilon = (E-1/2*u**2)
          eps = U(irhoE)/cell%rho - 0.5d0*(cell%u**2+cell%v**2)
          cell%p = EOS_p_from_rho_eps(cell%fluide%EOS, cell%rho, eps)
!!! vitesse du son
          call vitesse_du_son_melange_Allaire(cell%fluide%EOS%gamma, cell%fluide%EOS%pi, cell%p, cell%rho, cell%z, cell%c, ierreur)
!!! temperature
          cell%T = EOS_T_from_p_rho(fluide1%EOS, cell%p, cell%rho)
!!! enthalpie
          cell%h = eps + cell%p/cell%rho
       case(CHIMIE_FIGEE)
          call decode_chimie_loc(U, cell, ierreur)
       case default
          print*, "decode_u2cell : type de fluide inconnu"
          call arret_code
       end select
     case default
       call decode_diphasique_local(U, cell)
    end select

  end subroutine decode_U2cell

  subroutine decode_U(ldeb, lfin, mdeb, mfin, b, ierreur)
    use variables_globales
    use chimie

    integer        , intent(in)  :: ldeb, lfin, mdeb, mfin
    type (STR_BLOC), pointer     :: b
    integer        , intent(out) :: ierreur

    integer :: l, m
    real*8  :: eps
    real*8  :: entropie, hkapa, hk1, cp, cv, h
    real*8, dimension(:,:,:), pointer :: U
    real*8, dimension(:,:)  , pointer :: rho
    integer :: tid

    call debut_watchTime(wT_decodage_U)

    U => b%U
    rho => b%rho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! decodage des quantites thermo
!!! rho, ci, u_x, u_y, E
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--------------- DO PROJET ----------------!
    print *, " m = [", mdeb, "," , mfin  , "]"
    print *, " l = [", ldeb, "," , lfin  , "]"

  !$OMP PARALLEL
    !! tid = OMP_GET_THREAD_NUM()
  !$OMP DO
    do m = mdeb, mfin
       do l = ldeb, lfin

          call rho_from_U(U(:,l,m), rho(l,m))
          if (rho(l,m) < 0.d0) then
             print*, "densite negative", rho(l,m)
             ierreur = 1
          end if

          b%u_x(l,m) = U(irhou,l,m)/rho(l,m)
          b%u_y(l,m) = U(irhov,l,m)/rho(l,m)
          b%E(l,m)   = U(irhoE,l,m)/rho(l,m)

       end do
    end do
    !$OMP END DO

    print *, "I am thread", OMP_GET_THREAD_NUM(), " m = [", mdeb, "," , mfin  , "]"
    print *, "I am thread", OMP_GET_THREAD_NUM(), " l = [", ldeb, "," , lfin  , "]"

    !$OMP END PARALLEL

    !! Verification de la parallélisation
    ! print *, "somme b%u_x", sum(b%u_x)
    ! print *, "somme b%u_y", sum(b%u_y)
    ! print *, "somme b%E", sum(b%E)

    ! Open(15, File='parallele_parametres.dat')
    !
    ! write(15,*) b%u_x
    ! close(15)

    select case (MODELE_DIPHASIQUE)
    case (MONOPHASIQUE)
!!! on est en monophasique (utile ?)
       b%z(ldeb:lfin,mdeb:mfin) = 1.d0
       b%y(ldeb:lfin,mdeb:mfin) = 1.d0
       b%rho1(ldeb:lfin,mdeb:mfin) = rho(ldeb:lfin,mdeb:mfin)
       b%rho2(ldeb:lfin,mdeb:mfin) = 0.d0

       select case(fluide1%type)
       case(EOS)
!!! on stocke les cv1, gamma1, et cv, gamma, pi a partir du fluide1
!!! ca rempli les b%gamma_eq, b%pi_eq...
          b%cv1(ldeb:lfin,mdeb:mfin) = fluide1%EOS%cv
          b%gamma1(ldeb:lfin,mdeb:mfin) = fluide1%EOS%gamma

          do m = mdeb, mfin
             do l = ldeb, lfin
                call stocke_grandeurs_melange(fluide1, b, l, m)
             end do
          end do




!!! decodage de la pression
          select case (b%eq_energie)
          case (0)

             do m = mdeb, mfin
                do l = ldeb, lfin
                   fluide1%EOS%gamma = b%gamma_eq(l,m)
                   fluide1%EOS%pi = b%pi_eq(l,m)
                   fluide1%EOS%q = b%q_eq(l,m)

                   eps = U(irhoE,l,m)/rho(l,m) - 0.5d0*(b%u_x(l,m)**2+b%u_y(l,m)**2)
                   b%p(l,m) = EOS_p_from_rho_eps(fluide1%EOS, rho(l,m), eps)
                end do
             end do

          case (1) ! isentropique
 
             do m = mdeb, mfin
                do l = ldeb, lfin
                   b%p(l,m) = b%entropie_ref*rho(l,m)**fluide1%EOS%gamma - fluide1%EOS%pi
                end do
             end do

          case (2) ! isotherme
             print*, "isotherme non gere"
             call arret_code
          end select

!!! vitesse du son

          do m = mdeb, mfin
             do l = ldeb, lfin
                call vitesse_du_son_melange_Allaire(b%gamma_eq(l,m), b%pi_eq(l,m), b%p(l,m), rho(l,m), 1.d0, b%c(l,m), ierreur)
             end do
          end do


!!! Temperature
          if (b%visqueux) then

             do m = mdeb, mfin
                do l = ldeb, lfin
                   b%T(l,m) = EOS_T_from_p_rho(fluide1%EOS, b%p(l,m), rho(l,m))
                end do
             end do

          end if
       case (CHIMIE_FIGEE)
          call decode_chimie(ldeb, lfin, mdeb, mfin, b, ierreur)
       case default
          print*, "decode_U : type de fluide non defini en monophasique"
          call arret_code
       end select

       if (b%visqueux) then
          b%T1(ldeb:lfin,mdeb:mfin) = b%T(ldeb:lfin,mdeb:mfin)
          b%T2(ldeb:lfin,mdeb:mfin) = 0.d0
       end if

       ! calcul de la viscosite, de la conductivite et de la diffusion
       if (b%visqueux) call calcul_coeff_transport_b(ldeb, lfin, mdeb, mfin, fluide1, b)

    case (ALLAIRE) ! diphasique
       call decode_diphasique(ldeb, lfin, mdeb, mfin, b, ierreur)
       ! calcul de la viscosite, de la conductivite et de la diffusion
       if (b%visqueux) call calcul_coeff_transport_diphasique(ldeb, lfin, mdeb, mfin, fluide1, fluide2, b)
    case default
       print*, "decode U : cas non gere"
       call arret_code
    end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! enthalpie totale : E + p/rho
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do m = mdeb, mfin
       do l = ldeb, lfin
          b%Htot(l,m) = U(irhoE,l,m)/rho(l,m) + b%p(l,m)/rho(l,m)
       end do
    end do


  call fin_watchtime(wT_decodage_U)

  end subroutine decode_U

end module decodage
