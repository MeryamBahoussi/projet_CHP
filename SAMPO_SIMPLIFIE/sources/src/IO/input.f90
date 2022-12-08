module m_input
  use m_struct
  use parametres_globaux, only : FLUIDE, PETSc
  use parametres_fluide
  use decodage
  use variables_globales
  use parametres_IO
  USE MOD_Readintools
  use input_chimie

  implicit none

contains

  subroutine read_input(tache)
    use m_temps
    use m_verification
    use m_gestion_interbloc

    type (STR_TACHE), intent(inout) :: tache
    type (STR_SUPER_BLOC), pointer :: sb
    integer :: nb_super_bloc, nb_bloc
    integer :: i, isb
    integer :: nb_t_out

    !open(iunitInput,file=trim(path_input)//"input.dat",status="unknown")

    iter_ini   = GETINT('iter_deb',0)                ! Iteration de départ
    itermax    = GETINT('iter_max')                  ! Nombre d'iteration maximal
    tache%tmax = GETREAL("temps_Max")                ! Temps maximal de simulation
    film       = GETLOGICAL("sortie_instationnaire") ! sorties film visit ?
    output_1d  = GETLOGICAL("sortie_1D","F")         ! sorties 1d gnuplot ? (False par defaut)
    debug_MF   = GETLOGICAL("sortie_MF","F")         ! sorties des mailles fictives ? (False par defaut)
    it_restart = GETINT('prot_frequence',itermax+1)  ! taux d'impression des restarts
    it_debug   = GETINT('debug_frequence',itermax+1) ! taux d'impression des fichiers de debug

    tache%type_freqout = GETSTR("sortie_type_freq") ! type d'impression des resultats
    select case (trim(adjustl(tache%type_freqout)))
    case("iter")   ! sortie toutes les x iterations
       it_out = GETINT('sortie_frequence')
    case("temps")  ! sortie toutes les x instants
       it_out = itermax
       tache%freq_out = GETREAL("sortie_frequence")

!!!!  sortie en frequence temporelle
       nb_t_out = int(tache%tmax/tache%freq_out)
       if (nb_t_out > 500) then
          print*, "nombre de sorties temporelles trop grand"
          call arret_code
       end if
       if (tache%tmax /= nb_t_out*tache%freq_out) then
          ! le temps final n'est pas un multiple de la frequence
          allocate(tache%t_imposes(1:nb_t_out+1))
       else
          allocate(tache%t_imposes(1:nb_t_out))
       end if
       do i = 1, nb_t_out
          tache%t_imposes(i) = i*tache%freq_out
       end do
       if (tache%tmax /= nb_t_out*tache%freq_out) then
          tache%t_imposes(nb_t_out+1) = tache%tmax
       end if
       tache%ind_timp = 1
    case default
       print*, "il faut indique le type de sorties : iter ou temps"
       call arret_code
    end select

    nb_super_bloc = GETINT("nb_super_bloc",1) ! nombre de super blocs
    nb_bloc       = GETINT("nb_bloc",1)       ! nombre total de blocs

    tache%nb_membre = nb_super_bloc
    tache%fluide=.false.
    ! tache
    allocate(tache%membre(tache%nb_membre))
    do i = 1, tache%nb_membre
       tache%membre(i)=i
    end do
    ! super blocs :
    allocate(acces_super_bloc(nb_super_bloc))
    do i = 1, nb_super_bloc
       nullify(acces_super_bloc(i)%pointeur)
    end do
    ! blocs :
    allocate(acces_bloc(nb_bloc))
    do i = 1, nb_bloc
       nullify(acces_bloc(i)%pointeur)
    end do

    do i = 1, tache%nb_membre
       isb = tache%membre(i)
       allocate( acces_super_bloc(isb)%pointeur )
       sb => acces_super_bloc(isb)%pointeur

       call read_input_super_bloc(sb, isb)

       if (sb%nature==FLUIDE) tache%fluide=.true.
    end do

    if (nbprocs/=1) call decoupage
    
    !close(iunitInput)

    call verification_input(tache)

  end subroutine read_input

  subroutine read_input_super_bloc(sb, isb)

    type (STR_SUPER_BLOC), intent(inout) :: sb
    integer              , intent(in)    :: isb

    character(len=4) :: sb_num

    write(sb_num, FMT='(I4)') isb

    sb%nature = GETINT("sb_"//trim(adjustl(sb_num))//"_nature")

    select case(sb%nature)
    case(FLUIDE)
       call read_input_super_bloc_fluide(sb, "sb_"//trim(adjustl(sb_num)))
    case default
       if (numproc==0) print*, "Nature du super bloc inconnue"
       call arret_code
    end select

  end subroutine read_input_super_bloc

  subroutine read_input_super_bloc_fluide(sb, sb_num)
    type (STR_SUPER_BLOC), intent(inout) :: sb
    character(len=*)     , intent(in)    :: sb_num


    integer :: i, ib, ne
    type (STR_BLOC), pointer :: b
    integer :: ordre_Eulerdirect
    integer :: calcul_gradient, calcul_gradient_methode
    real*8  :: calcul_gradient_poids
    integer :: correction_son
    real*8  :: sigma, angle

    sigma = 0.d0
    angle = 0.d0

    calcul_gradient = -1
    calcul_gradient_methode = 0
    calcul_gradient_poids = 0.d0

    MODELE_DIPHASIQUE     = GETINT(trim(sb_num)//"_modele_diphasique")  ! MONOPHASIQUE ou modele ALLAIRE
    sb%visqueux           = GETLOGICAL(trim(sb_num)//"_visqueux","F")   ! visqueux ?
    sb%eq_energie         = GETINT(trim(sb_num)//"_eq_energie", "0")    ! 0:equation sur l'energie, 1:isentropique, 2:isotherme
    sb%gravite            = GETLOGICAL(trim(sb_num)//"_gravite","F")   ! gravite ?)

    if(sb%gravite) call lecture_gravite(sb)


    select case(sb%eq_energie)
    case (0)
       ! ok rien a lire de plus
    case (1)  ! isentropique
       if (MODELE_DIPHASIQUE /= MONOPHASIQUE) then
          print*, "isentropique seulement disponible en monophasique"
          call arret_code
       end if
       sb%entropie_ref = GETREAL(trim(sb_num)//"_entropie_ref")       ! entropie de référence
    case(2) ! isotherme
       if (MODELE_DIPHASIQUE /= MONOPHASIQUE) then
          print*, "isotherme seulement disponible en monophasique"
          call arret_code
       end if
       sb%temperature_ref = GETREAL(trim(sb_num)//"_temperature_ref") ! temperature de référence
    case default
       if (numproc==0) print*, "methode de résolution de l'énergie inconnue"
       call arret_code
    end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! gestion pas de temps !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    sb%sol_stat%actif = GETLOGICAL(trim(sb_num)//"_stationnaire", "F") ! calcul stationnaire ?

!!! cas ou on recherche un etat stationnaire
    if (sb%sol_stat%actif .eqv. .true.) then
       sb%sol_stat%kmax = GETINT(trim(sb_num)//"_kmax")           ! nombre d'iterations max pour converger
       sb%sol_stat%tol  = GETREAL(trim(sb_num)//"_tol")           ! critere d'arret pour converger

       sb%preconvergence%actif = GETLOGICAL(trim(sb_num)//"_preconvergence") ! preconvergence
       if (sb%preconvergence%actif) then
          sb%preconvergence%kmax = GETINT(trim(sb_num)//"_preconv_kmax")
          sb%preconvergence%tol  = GETREAL(trim(sb_num)//"_preconv_tol")

          sb%preconvergence%montee_cfl%actif = GETLOGICAL(trim(sb_num)//"_preconv_montee_cfl") ! montee en cfl dans la preconvergence ?
          if (sb%preconvergence%montee_cfl%actif) then
             sb%preconvergence%montee_cfl%iter_ini      = GETINT(trim(sb_num)//"_preconv_montee_cfl_it_ini")
             sb%preconvergence%montee_cfl%iter_trans    = GETINT(trim(sb_num)//"_preconv_montee_cfl_it_trans")
             sb%preconvergence%montee_cfl%cfl_ini       = GETREAL(trim(sb_num)//"_preconv_montee_cfl_cfl_ini")
             sb%preconvergence%montee_cfl%cfl_end       = GETREAL(trim(sb_num)//"_preconv_montee_cfl_cfl_end")
             sb%preconvergence%montee_cfl%cfl_croisiere = GETREAL(trim(sb_num)//"_preconv_montee_cfl_cfl_max")
          end if
       end if
    end if

    sb%montee_cfl%actif = GETLOGICAL(trim(sb_num)//"_montee_cfl", "F") ! montee en cfl ?
    if (sb%montee_cfl%actif) then
       sb%montee_cfl%iter_ini      = GETINT(trim(sb_num)//"_montee_cfl_it_ini")
       sb%montee_cfl%iter_trans    = GETINT(trim(sb_num)//"_montee_cfl_it_trans")
       sb%montee_cfl%cfl_ini       = GETREAL(trim(sb_num)//"_montee_cfl_cfl_ini")
       sb%montee_cfl%cfl_end       = GETREAL(trim(sb_num)//"_montee_cfl_cfl_end")
       sb%montee_cfl%cfl_croisiere = GETREAL(trim(sb_num)//"_montee_cfl_cfl_max")
    else
       sb%cfl = GETREAL(trim(sb_num)//"_cfl")
       sb%montee_cfl%cfl_croisiere = sb%cfl
    end if

!!! cfl_ini : correspond a la cfl initiale ou a celle du dernier pas de temps avant la reprise
    sb%cfl_ini = 0.d0
    if (iter_ini==0) then ! on n'est pas en reprise
       if (sb%sol_stat%actif) then
          if (sb%preconvergence%actif) then
             if (sb%preconvergence%montee_cfl%actif) then
                sb%cfl_ini = sb%preconvergence%montee_cfl%cfl_ini
             end if
          else
             if (sb%montee_cfl%actif) then
                sb%cfl_ini = sb%montee_cfl%cfl_ini
             else
                sb%cfl_ini = sb%cfl
             end if
          end if
       else
          if (sb%montee_cfl%actif) then
             sb%cfl_ini = sb%montee_cfl%cfl_ini
          else
             sb%cfl_ini = sb%cfl
          end if
       end if
    else
       if (sb%montee_cfl%actif) then
          sb%cfl_ini = sb%montee_cfl%cfl_croisiere
       else
          sb%cfl_ini = sb%cfl
       end if
    end if
    if (sb%cfl_ini == 0.d0) then
       print*, 'bug input : cfl_ini = 0'
       call arret_code
    end if
!!!!! fin gestion pas de temps

    sb%nb_membre = GETINT(trim(sb_num)//"_nb_membre", "1")  ! nombre de blocs
    allocate(sb%membre(sb%nb_membre))
    if (sb%nb_membre == 1) then
       sb%membre =  GETINTARRAY(trim(sb_num)//"_membre", sb%nb_membre, "1") ! numero des blocs du superbloc considere
    else
       sb%membre =  GETINTARRAY(trim(sb_num)//"_membre", sb%nb_membre) !  numero des blocs du superbloc considere
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! parametres numeriques
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    sb%solveur_Riemann = GETINT(trim(sb_num)//"_solveur_Riemann")                      ! solveur de Riemann
    sb%egalite_pentes_solveur = GETLOGICAL(trim(sb_num)//"_solveur_Riemann_eq_pentes") ! egalite des pentes du solveur de Riemann ?
    sb%K_Riemann = GETREAL(trim(sb_num)//"_solveur_Riemann_coeff_pentes")              ! coefficient K tel que a=K * rho*c
    correction_son = GETINT(trim(sb_num)//"_correction_son",0)
    ! terme qui permet de rajouter de la diffusion dans la partie décentrée
    sb%delta_epsilon = GETREAL(trim(sb_num)//"_delta_epsilon")
    sb%delta_epsilon_couche_lim = GETREAL(trim(sb_num)//"_delta_epsilon_couche_lim")

!!! parametres numeriques pour l'implicitation
    sb%implicite = GETLOGICAL(trim(sb_num)//"_implicite")    ! implicite ?
    if (sb%implicite .eqv. .true.) then
       ! solveur pour l'implicite: MKL=0, BALAYAGE=1, PETSc=2
       sb%solveur_implicite = GETINT(trim(sb_num)//"_implicite_solveur",PETSC)
       
       sb%inv_up%kmax = GETINT(trim(sb_num)//"_implicite_kmax") ! nb ite max pour l'inversion a la Chalons en U-p
       if (sb%inv_up%kmax /= 1) then
          sb%inv_up%tol = GETREAL(trim(sb_num)//"_implicite_tol") ! tolerance pour l'inversion a la Chalons en U-p
       end if
       if (sb%solveur_implicite == PETSC) then
          allocate(sb%inv_up%petsc) ! parametres petsc
          sb%inv_up%petsc%kmax = GETINT(trim(sb_num)//"_implicite_petsc_kmax")
          sb%inv_up%petsc%tol = GETREAL(trim(sb_num)//"_implicite_petsc_tol")
       end if

    end if
 
!!! ordre du schéma
    ordre_Eulerdirect = GETINT(trim(sb_num)//"_ordre",1)
    
    sb%nb_membre_input = sb%nb_membre

    sb%schema_temps = GETINT(trim(sb_num)//"_schema_temps",EULER_EXPLICITE)

    call read_parametres_fluides(sb%visqueux, fluide1, fluide2)
    ! ne : nombre d'elements ou d 'especes
    select case (fluide1%type)
    case(EOS)
       fluide1%ne = 1
    case(CHIMIE_FIGEE)
       fluide1%ne = fluide1%nb_esp
    case default
       print*, "type de fluide non defini"
       call arret_code
    end select
    ne = fluide1%ne
    call def_indices_fluide(ne)

    do i = 1, sb%nb_membre
       allocate(b)
       call read_input_bloc(b, sb%membre(i))
       call calc_stencil_bloc(sb, calcul_gradient, b%st_coin, b%st_grad)

       ! variables communes a tous les blocs d'un superbloc
       b%nature = sb%nature
       b%ne = ne
       b%neq = 6 + ne - 1  ! FLUIDE + QUANTITES TRANSPORTEES
       b%implicite = sb%implicite
       b%visqueux = sb%visqueux
       b%solveur_implicite = sb%solveur_implicite
       b%K_Riemann = sb%K_Riemann
       b%delta_epsilon = sb%delta_epsilon
       b%delta_epsilon_couche_lim = sb%delta_epsilon_couche_lim
       b%calcul_gradient         = calcul_gradient
       b%calcul_gradient_poids   = calcul_gradient_poids
       b%calcul_gradient_methode = calcul_gradient_methode
       b%eq_energie = sb%eq_energie
       b%entropie_ref = sb%entropie_ref
       b%temperature_ref = sb%temperature_ref

       ! on initialise pour ne pas tester des variables non initialisées
       b%ordre_Eulerdirect=1
       b%correction_son = correction_son
       b%ordre_Eulerdirect = ordre_Eulerdirect
       
       ib = b%numero
       if (.not. associated(acces_bloc(ib)%pointeur)) then
          allocate( acces_bloc(ib)%pointeur )
       else
          print*, " input.f90 : deux blocs ont le meme numero"
          call arret_code
       end if
       acces_bloc(ib)%pointeur = b

       deallocate(b)
    end do

  end subroutine read_input_super_bloc_fluide

  subroutine read_input_bloc(b, ib)

    type (STR_BLOC), pointer    :: b
    integer        , intent(in) :: ib

    integer :: iface, nb_patch(4)
    real*8  :: rtab_tmp(2)
    character(len=6) :: b_num
    type (STR_FACE), pointer :: face

    write(b_num, FMT='(I5)') ib
    b_num = "b_"//trim(adjustl(b_num))

    b%numero = GETINT(trim(b_num)//'_numero') ! numero du bloc
    if (b%numero /= ib) then
       print*, "erreur de numerotation dans la lecture du bloc ", b_num
       call arret_code
    end if
    b%num_input = b%numero

    b%type_initialisation = GETSTR(trim(b_num)//"_type_init") ! type d'initialisation

    b%maillage%type = GETSTR(trim(b_num)//"_maillage_type") ! type de maillage
    b%maillage%resserrer_maillage = GETLOGICAL(trim(b_num)//"_maillage_resserrer", "F") ! resserrer le maillage ?
    if (b%maillage%resserrer_maillage) then
       b%maillage%face = GETINT(trim(b_num)//"_maillage_resserrer_face")
       select case (b%maillage%face)
       case(1,2,12)  ! faces 1, 2 ou les deux
          b%maillage%beta_l = GETREAL(trim(b_num)//"_maillage_resserrer_coeff_l")
       case(3,4,34)  ! faces 3, 4 ou les deux
          b%maillage%beta_m = GETREAL(trim(b_num)//"_maillage_resserrer_coeff_m")
       case default
          b%maillage%beta_l = GETREAL(trim(b_num)//"_maillage_resserrer_coeff_l")
          b%maillage%beta_m = GETREAL(trim(b_num)//"_maillage_resserrer_coeff_m")
       end select
    end if
    rtab_tmp = GETREALARRAY(trim(b_num)//"_xin_xout", 2)
    b%xin  = rtab_tmp(1); b%xout = rtab_tmp(2)
    rtab_tmp = GETREALARRAY(trim(b_num)//"_yin_yout", 2)
    b%yin  = rtab_tmp(1); b%yout = rtab_tmp(2)

    b%ld_tot = GETINT(trim(b_num)//"_ld", 1)
    b%lf_tot = GETINT(trim(b_num)//"_lf")
    b%md_tot = GETINT(trim(b_num)//"_md", 1)
    b%mf_tot = GETINT(trim(b_num)//"_mf")

    ! nombre de patch par face
    nb_patch = GETINTARRAY(trim(b_num)//"_nb_patch",4,"1,1,1,1")

!!! de base il n'y a pas de deplacement du maillage
!!! ces variables peuvent etre modifiees par read_input_face
    b%maillage%move = .false.
    b%maillage%move_l = .false.
    b%maillage%move_m = .false.
    do iface = 1, 4
       allocate(face)
       call read_input_face(b, face, nb_patch(iface), iface)
!!! on peut donner les faces dans n'importe quel ordre
       b%face(face%patch(1)%bc%position) = face
       deallocate(face)
    end do

    b%numproc= numproc
    b%ld = b%ld_tot; b%lf = b%lf_tot
    b%md = b%md_tot; b%mf = b%mf_tot
    b%proc_ini = 0
    b%proc_fin = 0
    call MPI_charge_bloc(b)
    allocate(b%Nbcoupes(size(acces_bloc, 1), 2))
    b%NbCoupes(b%num_input,:)=1

  end subroutine read_input_bloc


  subroutine read_input_face(b, face, nb_patch, iface)

    type (STR_BLOC), pointer :: b
    type (STR_FACE), intent(inout) :: face
    integer, intent(in) :: nb_patch
    integer, intent(in) :: iface

    integer                               :: ip
    character(len=5)                      :: position
    character(len=2)                      :: p_num
    character(len=6)                      :: b_num
    character(len=20)                     :: p_name
    integer, dimension(4)                 :: tmp
    type (STR_PATCH)            , pointer :: patch
    type (STR_BoundaryCondition), pointer :: bc

    face%nb_patch = nb_patch
    allocate(face%patch(nb_patch))

    write(b_num, FMT='(I4)') b%numero
    b_num = "b_"//trim(adjustl(b_num))

    select case (iface)
    case(1)
       position = "West"
    case(2)
       position = "East"
    case(3)
       position = "South"
    case(4)
       position = "North"
    case default
       print*, "numero de face inconnu"
       call arret_code
    end select

    do ip = 1, face%nb_patch
       patch => face%patch(ip)
       bc => patch%bc

       if (face%nb_patch == 1) then
          p_name = trim(b_num)//"_CL_"//trim(position)
       else
          write(p_num, FMT='(I2)') ip
          p_name = trim(b_num)//"_CL_"//trim(position)//"_"//trim(adjustl(p_num))
       end if

       bc%condition = GETSTR(trim(p_name))
       ! si on a qu'un seul patch, pas besoin de donner les indices de debut et de fin du patch
       if (face%nb_patch == 1) then
          select case(position)
          case("West", "East")
             patch%ideb = b%md_tot
             patch%ifin = b%mf_tot
          case("North", "South")
             patch%ideb = b%ld_tot
             patch%ifin = b%lf_tot
          end select
       else
          patch%ideb = GETINT(trim(p_name)//"_ideb")
          patch%ifin = GETINT(trim(p_name)//"_ifin")
       end if

       select case(position)
       case("West")
          bc%position = 1
       case("East")
          bc%position = 2
       case("South")
          bc%position = 3
       case("North")
          bc%position = 4
       case default
          print*, "position de face inconnue", position
          print*, "West, East, South ou North"
       end select

!!! condition interbloc ou interface de couplage
       if ( bc%condition == "interbloc" .or.&
            (index(bc%condition,"itf")/=0 .and. bc%condition /="itf_fusion_sans_solide") ) then

          tmp  = GETINTARRAY(trim(p_name)//"_interbloc", 4)
          patch%ib_adj = tmp(1)
          patch%if_adj = tmp(2)
          patch%sens = tmp(3)
          patch%ideb_adj = tmp(4)
       else if (bc%condition =="itf_fusion_sans_solide") then
          patch%ib_adj = b%numero
          patch%if_adj = bc%position
          patch%sens=1
          patch%ideb_adj = patch%ideb
       else
          patch%ib_adj = -1  ! bord du domaine
       end if

!!! on regarde si on le maillage va bouger au cours du calcul
       if (index(bc%condition,"itf")/=0 .and. bc%condition /="itf_couplage") then
          b%maillage%move = .true.
          if (bc%position == 1 .or. bc%position == 2) b%maillage%move_l =.true.
          if (bc%position == 3 .or. bc%position == 4) b%maillage%move_m =.true.
          ! les variables move_l et move_m ne servent qu'a connaitre la direction
          ! de maillage (donc en (l,m)) suivant laquelle le maillage bouge.
          ! en curviligne, on peut avoir move_l = .false. mais X_t non nul
       end if

       ! lecture des donnees supplementaires pour certaines conditons aux limites
       ! entree supersonique ou subsonique, T imposee
       call read_donnees_condition_limites(patch, bc, trim(p_name))
    end do

  end subroutine read_input_face


  subroutine read_donnees_condition_limites(patch, bc, nom_patch)

    type (STR_PATCH), pointer :: patch
    type (STR_BoundaryCondition), pointer ::  bc
    character(len=*), intent(in) :: nom_patch

    ! entree ou sortie
    if ( index(bc%condition, "entree") /= 0 .or. &
         index(bc%condition, "sortie") /= 0 ) then
       call read_donnees_cell_infi(bc)
    end if

    if ( index(bc%condition, "paroi") /=0 ) then
       call read_donnees_cl_paroi(patch, bc, nom_patch)
    end if

    if (trim(bc%condition) == "injection_carbone") then
       call read_donnees_cl_injection_carbone(patch, bc, nom_patch)
    end if

  end subroutine read_donnees_condition_limites

  subroutine read_donnees_cl_injection_carbone(patch, bc, nom_patch)
    type (STR_PATCH)            , intent(in)    :: patch
    type (STR_BoundaryCondition), intent(inout) :: bc
    character(len=*)            , intent(in)    :: nom_patch

    integer           :: ne
    real*8            :: T_imposee

    if (fluide1%type /= CHIMIE_FIGEE) then
       print*, "injection_carbone seulement possible avec un melange de Gaz Parfaits"
       call arret_code
    end if

    ne = fluide1%ne

    allocate(bc%z_injection(1:patch%ifin-patch%ideb+1))   ! fraction volumique
    allocate(bc%u_injection(1:patch%ifin-patch%ideb+1))   ! vitesse de la paroi
    allocate(bc%rho_injection(1:patch%ifin-patch%ideb+1)) ! densite a la paroi
    allocate(bc%Tw(1:patch%ifin-patch%ideb+1))            ! temperature de paroi
    allocate(bc%Ci_paroi(ne,1:patch%ifin-patch%ideb+1))   ! concentration d'especes a la paroi
    allocate(bc%D_paroi(1:patch%ifin-patch%ideb+1))       ! coefficient de Diffusion a la paroi
    allocate(bc%K_paroi(1:patch%ifin-patch%ideb+1))       ! conductivite a la paroi
    allocate(bc%mu_paroi(1:patch%ifin-patch%ideb+1))      ! viscosite a la paroi
    allocate(bc%Hi_paroi(ne,1:patch%ifin-patch%ideb+1))   ! enthalpie a la paroi

    bc%z_injection(:) = 1.d0 ! on injecte du fluide 1


    T_imposee = GETREAL(nom_patch//"_T_imposee")
    bc%Tw(:) = T_imposee

  end subroutine read_donnees_cl_injection_carbone


  subroutine read_donnees_cl_paroi(patch, bc, nom_patch)
    type (STR_PATCH)            , intent(in)    :: patch
    type (STR_BoundaryCondition), intent(inout) :: bc
    character(len=*)            , intent(in)    :: nom_patch

    integer           :: iunit_fichier
    character(len=50) :: nom_fichier
    real*8            :: V_imposee, T_imposee, flux_impose

    nullify(bc%uw)
    nullify(bc%Tw)
    nullify(bc%flux_impose)

!!! attention la numerotation des tableaux de bc commencent à 1
!!! car la structure itf pointe vers ces tableaux
!!! (et on est oblige de commencer a 1 pour l'itf car ce n est pas obligatoirement
!!! la meme numerotation de part et d'autre de l'interface.)
    select case(trim(bc%condition))
    case("paroi", "paroi_entrainee", "paroi_entrainee_fichier")
       allocate(bc%uw(1:patch%ifin-patch%ideb+1)) ! vitesse de la paroi
    case("paroi_T_imposee", "paroi_entrainee_T_imposee", "paroi_T_imposee_fichier", &
         "paroi_entrainee_T_imposee_fichier")
       allocate(bc%Tw(1:patch%ifin-patch%ideb+1)) ! temperature de paroi
       allocate(bc%uw(1:patch%ifin-patch%ideb+1)) ! vitesse de la paroi
    case("paroi_flux_impose")
       allocate(bc%uw(1:patch%ifin-patch%ideb+1)) ! vitesse de la paroi
       allocate(bc%flux_impose(1:patch%ifin-patch%ideb+1))
    end select

    iunit_fichier = iunitInput+10
    select case(trim(bc%condition))
    case("paroi")
       bc%uw(:) = 0.d0
    case("paroi_entrainee")
       V_imposee = GETREAL(nom_patch//"_V_imposee")
       bc%uw(:) = V_imposee
    case("paroi_T_imposee")
       T_imposee = GETREAL(nom_patch//"_T_imposee")
       bc%Tw(:) = T_imposee
       bc%uw(:) = 0.d0
    case("paroi_entrainee_T_imposee")
       V_imposee = GETREAL(nom_patch//"_V_imposee")
       T_imposee = GETREAL(nom_patch//"_T_imposee")
       bc%uw(:) = V_imposee
       bc%Tw(:) = T_imposee
    case("paroi_entrainee_fichier")
       nom_fichier = GETSTR(nom_patch//"_V_imposee_fichier")
       open(unit = iunit_fichier, file=trim(path_input)//nom_fichier)
       read(iunit_fichier,*) bc%uw
       close(iunit_fichier)
    case("paroi_T_imposee_fichier")
       nom_fichier = GETSTR(nom_patch//"_T_imposee_fichier")
       open(unit = iunit_fichier, file=trim(path_input)//nom_fichier)
       read(iunit_fichier,*) bc%Tw
       close(iunit_fichier)
       bc%uw(:) = 0.d0
    case("paroi_entrainee_T_imposee_fichier")
       nom_fichier = GETSTR(nom_patch//"_T_u_imposee_fichier")
       open(unit = iunit_fichier, file=trim(path_input)//nom_fichier)
       read(iunit_fichier,*) bc%Tw
       read(iunit_fichier,*) !""""""""Vitesse paroi"""""""
       read(iunit_fichier,*) bc%uw
       close(iunit_fichier)

!!! flux_imposee
    case( "paroi_flux_impose" )
       flux_impose = GETREAL(nom_patch//"_flux_impose")
       bc%flux_impose(:) = flux_impose
       bc%uw(:) = 0.d0

    case default
       print*, "type de condition aux limites de paroi inconnue", trim(bc%condition)
       call arret_code
    end select

    select case(trim(bc%condition))
    case("paroi", "paroi_entrainee", "paroi_entrainee_fichier")
       bc%condition = "paroi_global_adiabatique"
    case("paroi_T_imposee", "paroi_entrainee_T_imposee", "paroi_T_imposee_fichier", &
         "paroi_entrainee_T_imposee_fichier")
       bc%condition = "paroi_global"
    end select


  end subroutine read_donnees_cl_paroi

!!! lecture des donnees d'initialisation
  subroutine read_donnees_initialisation(b)

    type (STR_BLOC), intent(inout) ::  b
    integer :: i, nb_zones_ini

    open(iunitIniFld,file=trim(path_input)//"ini_fld.dat",status="unknown")
    read(iunitIniFld,*) !fluide1
    read(iunitIniFld,*) !fluide2
    read(iunitIniFld,*) nb_zones_ini  ! nombre de zones pour l'initialisation
    allocate(b%tab_cl_ini(nb_zones_ini))

    do i = 1, nb_zones_ini
       call read_donnees_fluide(b%tab_cl_ini(i), iunitIniFld)
    end do

    close(iunitIniFld)

!!!! verif si on a le bon nombre de zones pour le type d'initialisation
    if ( index(b%type_initialisation,"restart") == 0 ) then
       select case(b%type_initialisation)
       case("vortex", "vortex_diphasique")
          print*, "L'initialisation est définie en dur dans le code pour ce cas"
       case("zone1", "melange", "vitesse_lineaire")
          if (nb_zones_ini < 1) then
             print*, "read_donnees_initialisation : pas assez de zones pour l'initialisation"
             call arret_code
          end if
       case("zone2", "sod_1d", "sod_1d_trizone", "bulle_1d", "bulle_2d", "double_bulle_2d", "carre_2d", "bulle_murrone", "couette",&
            "steady_discontinuity_advection","zalesak_disk_advection","etoile","circle","etoile2",&
            "creneau","Kothe_Rider","kh_instability","rabbit","ellipse_2d")
          if (nb_zones_ini < 2) then
             print*, "read_donnees_initialisation : pas assez de zones pour l'initialisation"
             call arret_code
          end if
       case("bulle_choc_1d", "bulle_choc")
          if (nb_zones_ini < 3) then
             print*, "read_donnees_initialisation : pas assez de zones pour l'initialisation"
             call arret_code
          end if
       case("sod_2d")
          if (nb_zones_ini < 4) then
             print*, "read_donnees_initialisation : pas assez de zones pour l'initialisation"
             call arret_code
          end if
       case default
          print*,""
          print*," Type d'initialisation non definie", b%type_initialisation
          print*,""
          print*,"ERREUR : read_donnees_initialisation"
          call arret_code
    end select
 end if
  end subroutine read_donnees_initialisation

  subroutine read_donnees_fluide(cell, iunit)
    use outils_fluide, only : recup_ind_esp, recup_ind_ele
    use chimie, only : xele_from_cele

    type (STR_cell), intent(inout) :: cell
    integer, intent(in) :: iunit

    integer :: i, j, ne
    real*8  :: Xi_read, Cv, Cp
    real*8  :: eps, masse_totale
    real*8, dimension(:), allocatable :: Xi
    character(len=8) :: nomEsp

    read(iunit,*) !################
    read(iunit,*) !#    ZONE i    #
    read(iunit,*) !################
    read(iunit,*) cell%z
    ! Dans le cas de l'équilibre chimique, on rentre p et T
    ! sinon on donne rho et p
    if (cell%z==1.d0) then
       read(iunit,*) cell%rho1
       cell%rho2 = 0.d0
    else if (cell%z==0.d0) then
       read(iunit,*) cell%rho2
       cell%rho1 = 0.d0
    else
       read(iunit,*) cell%rho1
       read(iunit,*) cell%rho2
    end if
    
!!! calcul des fractions massiques Ci
    if (cell%z /= 0.d0) then ! si on a du fluide 1
       ne = fluide1%ne
       allocate(cell%Ci(1:ne))
       select case(fluide1%type)
       case(EOS, GAZ_PARFAIT)
          cell%Ci(1) = 1.d0
       case(CHIMIE_FIGEE)
          allocate(Xi(1:ne))
          Xi(:) = -1.d30
          masse_totale = 0.d0  ! masse molaire totale
          do i = 1, ne
             read(iunit,*) nomEsp, Xi_read

             ! on recupere l'indice de l'espece qu'on vient de lire
             call recup_ind_esp(nomEsp, j)
             Xi(j) = Xi_read

             masse_totale = masse_totale + Xi(j) * fluide1%especes(j)%masse
             cell%Ci(j) = Xi(j) * fluide1%especes(j)%masse
          end do
          cell%Ci(:) = cell%Ci(:) / masse_totale

          if (minval(Xi) == -1.d30) then
             print*, "Probleme dans les conditions initiales. Une espece n'a pas ete lue."
             call arret_code
          end if
          if (sum(Xi(:)) /= 1.d0) then
             print*, "Probleme dans les conditions initiales. La somme des fractions molaires n'est pas egale a 1."
             call arret_code
          end if
          deallocate(Xi)

          call calcul_Cv_Cp_fige(fluide1, fluide1%especes, cell%ci(:), Cv, Cp)

       case default
          print*, "type de fermeture pour le fluide 1 non defini"
          call arret_code
       end select
    end if

    read(iunit,*) cell%u
    read(iunit,*) cell%v
    read(iunit,*) cell%p

    select case(fluide1%type)
    case( EOS, CHIMIE_FIGEE)
       ! Loi d'etat
       cell%rho = cell%rho1*cell%z+cell%rho2*(1.d0-cell%z)
       cell%y = cell%rho1*cell%z/cell%rho
       ! on determine gamma, pi, et q
       call grandeurs_melange(fluide1, fluide2, cell%z, cell%y, cell%fluide)
       eps = EOS_eps_from_p_rho(cell%fluide%EOS, cell%p, cell%rho)
       cell%h = eps + cell%p/cell%rho
    case default
       call arret_code()
    end select

  end subroutine read_donnees_fluide

  subroutine read_donnees_cell_infi(bc)
    use variables_globales, only : Htot_infini_amont

    type (STR_BoundaryCondition), intent(inout) ::  bc
    integer :: ne

    if ( index(bc%condition, "subsonique") /= 0 ) then
       open(iunitCLInfi,file=trim(path_input)//"cell_infi_subsonique.dat",status="unknown")
    else if ( index(bc%condition, "supersonique") /= 0 ) then
       open(iunitCLInfi,file=trim(path_input)//"cell_infi.dat",status="unknown")
    else
       print*, "condition limite d'entree ni subsonique ni supersonique"
       call arret_code
    end if
    allocate(bc%cell_infi)
    call read_donnees_fluide(bc%cell_infi, iunitCLInfi)
    close(iunitCLInfi)

    ne=fluide1%ne

    allocate(bc%U_infi(1:6+ne-1)) ! 6 = b%neq
    ! encodage puis decodage
    call encode_U(bc%cell_infi, bc%U_infi)
    call decode_U2cell(bc%U_infi, bc%cell_infi)

!!! on stocke l'enthapie totale infini amont pour annuler la correction d'entropie dans la couche limite (schema de Roe)
    Htot_infini_amont = (bc%U_infi(irhoE) + bc%cell_infi%p)/bc%cell_infi%rho
  end subroutine read_donnees_cell_infi

  subroutine read_parametres_fluides(visqueux, fluide1, fluide2)

    logical, intent(in) :: visqueux
    type (STR_fluide), intent(out) :: fluide1, fluide2

    character(len=90) :: fichier_fluide1, fichier_fluide2

    open(iunitIniFld,file=trim(path_input)//"ini_fld.dat",status="unknown")
    read(iunitIniFld,*) fichier_fluide1
    read(iunitIniFld,*) fichier_fluide2
    close(iunitIniFld)

!!! fluide 1
    if ( index(fichier_fluide1 ,"parfait") /= 0 ) then
       call read_parametres_gaz_parfait(fichier_fluide1, visqueux, fluide1)
    else if ( index(fichier_fluide1 ,"stiffened") /= 0 ) then
       call read_parametres_stiffened_gas(fichier_fluide1, visqueux, fluide1)
    else if ( index(fichier_fluide1 ,"mie_gruneisen") /= 0 ) then
       call read_parametres_mie_gruneisen(fichier_fluide1, visqueux, fluide1)
    else if ( index(fichier_fluide1 ,"melange_GP") /= 0 ) then
       call read_parametres_melange_GP(fichier_fluide1, fluide1)
    else
       print*,""
       print*,"Loi d'etat pour le fluide 1 inconnue"
       print*,""
       call arret_code
    end if

 !!! fluide 2
    if ( index(fichier_fluide2 ,"parfait") /= 0 ) then
       call read_parametres_gaz_parfait(fichier_fluide2, visqueux, fluide2)
    else if ( index(fichier_fluide2 ,"stiffened") /= 0 ) then
       call read_parametres_stiffened_gas(fichier_fluide2, visqueux, fluide2)
    else if ( index(fichier_fluide2 ,"mie_gruneisen") /= 0 ) then
       call read_parametres_mie_gruneisen(fichier_fluide2, visqueux, fluide2)
    else
       print*,""
       print*,"Loi d'etat pour le fluide 2 inconnue"
       print*,""
       call arret_code
    end if

  end subroutine read_parametres_fluides

  subroutine read_parametres_gaz_parfait(fichier, visqueux, fluidek)

    character(len=90), intent(in)  :: fichier
    logical          , intent(in)  :: visqueux
    type(STR_fluide) , intent(out) :: fluidek

    fluidek%type = EOS
    fluidek%EOS%nom = gaz_parfait
    fluidek%nb_esp = 1
    fluidek%nb_ele = 1

    open(iunitDataFld,file=trim(path_input)//trim(fichier)//".dat",status="unknown")
    read(iunitDataFld,*) fluidek%EOS%gamma
    if (visqueux) then
       call read_proprietes_transport(iunitDataFld, fluidek)
       read(iunitDataFld,*) fluidek%EOS%Cv
    else
       fluidek%loi_viscosite    = -1
       fluidek%loi_conductivite = -1
       fluidek%loi_diffusion    = -1
       fluidek%EOS%Cv = 1.d80
    end if
    close(iunitDataFld)

    fluidek%EOS%pi=0.d0
    fluidek%EOS%q=0.d0
    fluidek%EOS%T0=0.d0
    !fluidek%EOS%rho0=0.d0
    fluidek%EOS%p0=0.d0
    fluidek%EOS%eps0=0.d0

    fluidek%EOS%Cp = fluidek%EOS%gamma*fluidek%EOS%Cv

  end subroutine read_parametres_gaz_parfait

  subroutine read_parametres_stiffened_gas(fichier, visqueux, fluidek)

    character(len=90), intent(in) :: fichier
    logical, intent(in) :: visqueux
    type(STR_fluide), target, intent(out) :: fluidek

    fluidek%type = EOS
    fluidek%EOS%nom = stiffened_gas
    fluidek%nb_esp = 1
    fluidek%nb_ele = 1

    open(iunitDataFld,file=trim(path_input)//trim(fichier)//".dat",status="unknown")
    read(iunitDataFld,*) fluidek%EOS%gamma, fluidek%EOS%pi
    if (visqueux) then
       read(iunitDataFld,*) fluidek%enth_formation, fluidek%EOS%T0, fluidek%EOS%p0, fluidek%EOS%rho0

       ! viscosite dynamique, conductivite thermique, capacite thermique
       call read_proprietes_transport(iunitDataFld, fluidek)
       read(iunitDataFld,*) fluidek%EOS%Cv
    else
       fluidek%loi_viscosite    = -1
       fluidek%loi_conductivite = -1
       fluidek%loi_diffusion    = -1
       fluidek%EOS%Cv = 1.d80
    end if
    close(iunitDataFld)

    fluidek%EOS%q=0.d0
    !fluidek%EOS%T0=0.d0
    !fluidek%EOS%rho0=0.d0
    !fluidek%EOS%p0=0.d0
    fluidek%EOS%eps0= 0.d0

    fluidek%EOS%Cp = fluidek%EOS%gamma*fluidek%EOS%Cv

!!$!!! modifs
!!$    fluidek%eps0 = fluidek%enth_formation-fluidek%p0/fluidek%rho0
!!$    fluidek%q = fluidek%enth_formation - fluidek%Cv*fluidek%T0 - (fluidek%p0+fluidek%pi)/fluidek%rho0

    !!! les parametres rho0, T0 et p0 sont ils vraiment utiles en stiffened gaz ??

  end subroutine read_parametres_stiffened_gas

  subroutine read_parametres_mie_gruneisen(fichier, visqueux, fluidek)

    character(len=90), intent(in) :: fichier
    logical, intent(in) :: visqueux
    type(STR_fluide), target, intent(out) :: fluidek

    real*8 :: alpha_p, beta_T
    real*8 :: gGamma, gPi
    real*8, pointer :: gamma, pi, q, Cv, rho0, T0, p0, eps0, enth_formation

    gamma => fluidek%EOS%gamma
    pi    => fluidek%EOS%pi
    q     => fluidek%EOS%q
    rho0  => fluidek%EOS%rho0
    T0    => fluidek%EOS%T0
    p0    => fluidek%EOS%p0
    Cv    => fluidek%EOS%Cv
    eps0  => fluidek%EOS%eps0
    enth_formation => fluidek%enth_formation

    fluidek%type = EOS
    fluidek%EOS%nom = mie_gruneisen
    fluidek%nb_esp = 1
    fluidek%nb_ele = 1

    open(iunitDataFld,file=trim(path_input)//trim(fichier)//".dat",status="unknown")
    read(iunitDataFld,*) alpha_p, beta_T
    read(iunitDataFld,*) enth_formation, T0, p0, rho0
    if (visqueux) then
       ! viscosite dynamique, conductivite thermique, capacite thermique
       call read_proprietes_transport(iunitDataFld, fluidek)
       read(iunitDataFld,*) fluidek%EOS%Cv
    else
       fluidek%loi_viscosite    = -1
       fluidek%loi_conductivite = -1
       fluidek%loi_diffusion    = -1
       fluidek%EOS%Cv = 1.d80
    end if
    close(iunitDataFld)

    ! grand gamma et grand pi
    gGamma = alpha_p/(Cv*rho0*beta_T)
    gPi = 1.d0/beta_T - gGamma*p0 + gGamma**2.d0*Cv*T0*rho0

    gamma = gGamma+1.d0
    pi = (gPi-P0)/gamma
    eps0 = enth_formation-p0/rho0
    q = eps0 - gPi/(rho0*(gamma-1.d0))

    ! pas sur, a verifier !! (formule de Manu)
    fluidek%EOS%Cp = Cv * ( 1.d0 + (gamma-1.d0)**2.d0 * Cv*T0*rho0 / gPi )

  end subroutine read_parametres_mie_gruneisen

  subroutine read_parametres_melange_GP(fichier, fluidek)
    use parametres_fluide, only : R_gaz

    character(len=90), intent(in)  :: fichier
    type(STR_fluide) , intent(out) :: fluidek

    integer :: i, ne

    fluidek%type = CHIMIE_FIGEE
    fluidek%EOS%nom = -1
    fluidek%nb_ele = 1

    open(iunitDataFld,file=trim(path_input)//trim(fichier)//".dat",status="unknown")
    read(iunitDataFld,*) fluidek%nb_esp         ! nombre d'especes

    call read_proprietes_transport(iunitDataFld, fluidek)

    ne = fluidek%nb_esp
    allocate(fluidek%especes(1:ne))

    do i = 1, ne
       read(iunitDataFld,*) fluidek%especes(i)%nom, fluidek%especes(i)%Mu, fluidek%especes(i)%masse! , fluidek%especes(i)%Cv
    end do
    close(iunitDataFld)

!!! Calcul du R par especes R = R_gaz/masse
    do i = 1, ne
       fluidek%especes(i)%R = R_gaz/fluidek%especes(i)%masse
    end do

    fluidek%EOS%pi=0.d0
    fluidek%EOS%q=0.d0
    fluidek%EOS%T0=0.d0
    !fluidek%EOS%rho0=0.d0
    fluidek%EOS%p0=0.d0
    fluidek%EOS%eps0=0.d0

    ! pour savoir si on est diatomique ou non
    ! degueu et pas bon (pour le sodium Na par exemple)
    do i = 1, ne
       if (len(trim(fluidek%especes(i)%nom)) == 1) then
          fluidek%especes(i)%diatomique = 0
       else
          fluidek%especes(i)%diatomique = 1
       end if
    end do

  end subroutine read_parametres_melange_GP


!!! lecture des proprietes de transport
  ! lois de viscosite dynamique, de conductivite thermique et de diffusion
  ! ainsi que les paramètres associes
  subroutine read_proprietes_transport(iunit, fluidek)
    integer         , intent(in)    :: iunit
    type(STR_fluide), intent(inout) :: fluidek

    character(len=90) :: loi_viscosite, loi_conductivite, loi_diffusion

!!! lecture de la loi de viscosite
    fluidek%mu_cst = 1.d80
    read(iunit,*) loi_viscosite
    if ( index(loi_viscosite ,"constante") /= 0 ) then
       fluidek%loi_viscosite = CONSTANTE
       backspace(iunit)
       read(iunit,*) loi_viscosite, fluidek%mu_cst
    else if ( index(loi_viscosite ,"sutherland") /= 0 ) then
       fluidek%loi_viscosite = SUTHERLAND
    else if ( index(loi_viscosite ,"wilke") /= 0 ) then
       if (fluidek%type /= CHIMIE_FIGEE) then
          print*, "Loi de WILKE pour la viscosite seulement disponible en chimie figee"
          call arret_code
       end if
       fluidek%loi_viscosite = WILKE
    else
       print*, "Loi de viscosite inconnue", loi_viscosite
       call arret_code
    end if

!!! lecture de la loi de conductivite
    fluidek%K_cst = 1.d80
    read(iunit,*) loi_conductivite
    if ( index(loi_conductivite ,"constante") /= 0 ) then
       fluidek%loi_conductivite = CONSTANTE
       backspace(iunit)
       read(iunit,*) loi_conductivite, fluidek%K_cst
    else if ( index(loi_conductivite ,"prandtl") /= 0 ) then
       fluidek%loi_conductivite = PRANDTL
       backspace(iunit)
       read(iunit,*) loi_conductivite, fluidek%Pr
    else
       print*, "Loi de conductivite inconnue", loi_conductivite
       call arret_code
    end if

!!! lecture de la loi de diffusion
    fluidek%D_cst = 1.d80
    select case (fluidek%type)
    case (EOS)
       ! rien a lire
    case (CHIMIE_FIGEE)
       read(iunit,*) loi_diffusion
       if ( index(loi_diffusion ,"constante") /= 0 ) then
          fluidek%loi_diffusion = CONSTANTE
          backspace(iunit)
          read(iunit,*) loi_diffusion, fluidek%D_cst
       else if ( index(loi_diffusion ,"lewis") /= 0 ) then
          fluidek%loi_diffusion = LEWIS
          backspace(iunit)
          read(iunit,*) loi_diffusion, fluidek%Le
       else
          print*, "Loi de diffusion inconnue", loi_diffusion
          call arret_code
       end if
    case default
       call arret_code
    end select
  end subroutine read_proprietes_transport

  subroutine calc_stencil_bloc(sb, calcul_gradient, st_coin, st_grad)
    type (STR_SUPER_BLOC), intent(inout) :: sb
    integer              , intent(in)    :: calcul_gradient
    integer              , intent(out)   :: st_coin, st_grad

    st_coin = 0
    st_grad = -1

    ! stencil pour les calculs des gradients
    select case (calcul_gradient)
       case (WLSQ5,GG)
          st_grad = 5
       case (WLSQ9)
          st_grad = 9
    end select

  end subroutine calc_stencil_bloc

  subroutine lecture_gravite(sb)
    use variables_globales
    type (STR_SUPER_BLOC), intent(inout) :: sb
      real*8 :: g(2)

      open(unit=21,file=trim(path_input)//"parametres_gravite.dat")
      read(21,*) g(:)
      close(21)

      acc_pesanteur = g
  end subroutine lecture_gravite

end module m_input
