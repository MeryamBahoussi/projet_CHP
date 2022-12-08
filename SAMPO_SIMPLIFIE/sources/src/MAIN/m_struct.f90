module m_struct
  use m_struct_solveur, only : STR_PETSC, STR_CSR
  implicit none

  type STR_montee_cfl
     logical :: actif
     integer :: iter_ini, iter_trans
     real*8  :: cfl_ini, cfl_end, cfl_croisiere
  end type STR_montee_cfl

  type derivees_partielles
     real*8, dimension(:), allocatable :: d_drho1z1Ci
     real*8                            :: d_drho2z2, d_drhoeps, d_dz
  end type derivees_partielles

  type variable_reconstruite  ! reconstruction pour l'ordre 2
!!! pour la partie acoustique du splitting, on reconstruit : ux, uy, p
     real*8, dimension(:,:)  , pointer :: uxL, uyL, pL
     real*8, dimension(:,:)  , pointer :: uxR, uyR, pR

!!! dans le cas Euler direct, on reconstruit : (rhoyci)_1,ne, rhoy2, ux, uy, p et z
     real*8, dimension(:,:,:), pointer :: rhoyciL, rhoyciR
     real*8, dimension(:,:)  , pointer :: rhoy2L, zL
     real*8, dimension(:,:)  , pointer :: rhoy2R, zR
     ! on en deduit E et c pour le solveur de Riemann
     real*8, dimension(:,:)  , pointer :: EL, cL
     real*8, dimension(:,:)  , pointer :: ER, cR
  end type variable_reconstruite

  ! variables utilisees par un processus iteratif (Newton ou recherche etat stat par ex)
  type process_iteratif
     logical               :: actif
     integer               :: kmax       ! nombre d'iteration maximal
     real*8                :: tol        ! tolerance
     type (STR_PETSC), pointer :: petsc      ! pour le solveur Petsc si il y a une inversion
     type (STR_montee_cfl) :: montee_cfl ! montee en cfl si besoin
  end type process_iteratif

  type matrice_tridiagonale
     integer :: dim_bloc, dim_mat
     real*8, dimension(:, :, :), pointer :: A, B, C
  end type matrice_tridiagonale

  type STR_MATRICE
     integer :: stencil   ! 5 ou 9 points
     integer :: dim_bloc
     real*8, dimension(:,:,:,:), pointer :: A,B,C,D,E
     real*8, dimension(:,:,:,:), pointer :: AD,CD,AE,CE  ! coins
  end type STR_MATRICE

  type STR_EOS
     integer :: nom      ! Gaz Parfait, Stiffened Gas ou Mie Gruneisen
     real*8  :: gamma    ! coefficient adiabatique (Loi d'etat)
     real*8  :: pi       ! pression de reference   (Loi d'etat)
     real*8  :: q        ! energie de formation    (Loi d'etat)
     real*8  :: Cv       ! capacite thermique a volume constant
     real*8  :: Cp       ! capacite thermique a pression constante
     real*8  :: rho0     ! densite de reference
     real*8  :: T0       ! temperature de reference
     real*8  :: p0       ! pression de reference
     real*8  :: eps0     ! energie interne de reference
  end type STR_EOS

  !****************************************************************
  !*                                                              *
  !*       Structure informatique pour la chimie                  *
  !*           especes,                                           *
  !*           elements,                                          *
  !*           reactions,                                         *
  !*                                                              *
  !****************************************************************

  type STR_REACTION
     integer :: nesp                            ! nombre d'especes dans la réaction chimique
     integer :: D_nu                            ! somme des coefficients stoechio de la réaction
     integer, dimension(:), pointer :: nu       ! coefficients stochio des espèces dans la réaction
     integer, dimension(:), pointer :: ist      ! numéro de l'espece correspondante dans le tableau esp global
     integer, dimension(:), pointer :: ist_inv  ! l'inverse de ist
  end type STR_REACTION

  type STR_ELEMENT
     character(2)     :: nom               ! nom
     real*8           :: masse             ! masse molaire
     integer, pointer :: liste_especes(:)  ! liste des indices des especes qui contiennent cet element
     integer, pointer :: multiplicite(:)   ! nombre d'especes i dans cet element
  end type STR_ELEMENT

  type STR_ESPECE
     character(8)     :: nom               ! nom
     integer          :: diatomique        ! 1 si l'espece est diatomique, 0 sinon
     real*8           :: masse             ! masse molaire
     real*8           :: R                 ! R universel divise par la masse molaire
     ! real*8           :: Cv                ! capacite thermique a volume constant
     real*8           :: Mu                ! viscosite
     real*8           :: h0                ! enthalpie molaire de formation
     real*8           :: s0                ! entropie molaire de formation

     integer, pointer :: liste_elements(:) ! liste des indices des elements qui contenus cette espece
     integer, pointer :: multiplicite(:)   ! nombre d'element k dans cette espece

     !integer          :: nb_poly_glenn     ! nombre de polynomes de Glenn
     !integer          :: nb_coeff_glenn    ! nombre de coefficients par polynomes de Glenn
     real*8 , pointer :: coeff_glenn(:,:)  ! coefficients des polynomes de Glenn (notés Ai)
     real*8 , pointer :: coeff_C0(:)       ! coefficients C0 pour chaque polynome
     real*8 , pointer :: coeff_C1(:)       ! coefficients C1 pour chaque polynome
  end type STR_ESPECE

  type STR_fluide
     integer :: type    ! fermeture thermodynamique : EOS, chimie figée ou equilibre chimique
     integer :: ne      ! nombre d'especes ou d'elements (correspond au vecteur conservatif)
     integer :: nb_esp  ! nombre d'especes chimique
     integer :: nb_ele  ! nombre d'elements

     type (STR_EOS)               :: EOS
     type (STR_ESPECE)  , pointer :: especes(:)
     type (STR_ELEMENT) , pointer :: elements(:)
     type (STR_REACTION), pointer :: reactions(:)
     integer            , pointer :: mat_multiplicite(:,:)  ! Matrice des multiplicites (M_ij : nombre d'espece i dans l'element j)

     ! proprietes de transport
     integer :: loi_viscosite     ! loi de viscosite dynamique
     integer :: loi_conductivite  ! loi de conductivite thermique
     integer :: loi_diffusion     ! loi de diffusion inter especes
     real*8 :: mu_cst             ! viscosite dynamique
     real*8 :: k_cst              ! conductivite thermique
     real*8 :: D_cst              ! coefficient de diffusion
     real*8 :: Pr                 ! nombre de Prandtl
     real*8 :: Le                 ! nombre de Lewis

     real*8 :: enth_formation ! enthalpie de formation

  end type STR_fluide

  type STR_cell  ! peut on le supprimer ??
     type (STR_fluide) :: fluide  ! donnees du fluide
     real*8 :: rho      ! densite du melange
     real*8 :: rho1     ! densite du fluide 1
     real*8 :: rho2     ! densite du fluide 2
     real*8 :: z        ! fraction volumique du fluide 1
     real*8 :: y        ! fraction massique du fluide 1
     real*8 :: u,v      ! composantes de la vitesse
     real*8 :: p        ! pression
     real*8 :: h        ! enthalpie ( h=epsilon + p/rho ) avec epsilon l'energie interne
     real*8 :: c        ! vitesse du son
     real*8 :: T        ! temperature
     real*8 :: T1, T2   ! temperature de chaque phase
     real*8, dimension(:), pointer :: Ci   ! fraction massique  (Multi-especes)
     real*8, dimension(:), pointer :: x_esp! fraction molaire   (Chimie)
  end type STR_cell

  type STR_gradient     ! gradient associé a un sommet
     real*8 :: x, y !, t
  end type STR_gradient

  type STR_BoundaryCondition
     integer                         :: position  ! 1:West, 2:East, 3:South, 4:North
     character(len=50)               :: condition ! type de CL
     type (STR_cell)       , pointer :: cell_infi
     real*8, dimension(:)  , pointer :: U_infi
     real*8, dimension(:)  , pointer :: V_infi
     real*8, dimension(:)  , pointer :: flux_impose
     real*8, dimension(:)  , pointer :: u_injection
     real*8, dimension(:)  , pointer :: rho_injection
     real*8, dimension(:)  , pointer :: z_injection
     real*8, dimension(:,:), pointer :: Ci_paroi, Hi_paroi
     real*8, dimension(:)  , pointer :: Tw, uw    ! vitesse et temperature de paroi
     real*8, dimension(:)  , pointer :: D_paroi
     real*8, dimension(:)  , pointer :: K_paroi
     real*8, dimension(:)  , pointer :: mu_paroi
  end type STR_BoundaryCondition

  type STR_PATCH
     integer :: ib_adj   ! numero du bloc adjacent (-1 si c'est le bord)
     integer :: if_adj   ! numero de la face adjacente
     integer :: sens     ! sens de parcours (+- 1)  (pas implemente)
     integer :: ideb, ifin
     integer :: ideb_adj
     type (STR_BoundaryCondition)    :: BC
  end type STR_PATCH

  type STR_FACE
     integer :: nb_patch
     type (STR_PATCH), dimension(:), pointer :: patch
  end type STR_FACE

  type STR_COIN
     integer :: icoin                            ! numero du coin
     integer :: condition                        ! type de CL
     integer :: indices(8)                       ! ld, lf, md, mf du coin
     integer :: ib_opp                           ! numero du bloc oppose (vis a vis du coin)
     type (STR_PATCH)    , pointer :: patch1, patch2 ! patchs contenant le coin
     real*8, dimension(:), pointer :: n1, n2         ! normales (face1 et face2)
  end type STR_COIN

!!! Infos sur le maillage
  type STR_MAILLAGE
     character(len=90) :: type       ! type de maillage
     logical :: resserrer_maillage   ! maillage raffine sur les bords?
     integer :: face                 ! sur quelle face on raffine (1,2,3,4,12,34,1234,...)
     real*8  :: beta_l               ! coefficient en l
     real*8  :: beta_m               ! coefficient en m
     logical :: move                 ! le maillage bouge-t-il au cours du calcul ?
     logical :: move_l               ! en l?
     logical :: move_m               ! en m?
  end type STR_MAILLAGE

  type STR_BLOC
     integer :: num_input                         ! numero du bloc donne en entree
     integer :: numero                            ! numero du bloc dans le code
     integer :: nature                            ! nature du bloc (0:Fluide, 1:Solide)
     integer :: neq                               ! nombre d'equations (FLUIDE = 6)
     integer :: ne                                ! nombre d'especes ou d'elements (Multi-especes)
     integer :: ld_tot, lf_tot                    ! numerotation en l donnee en entree
     integer :: md_tot, mf_tot                    ! numerotation en m donnee en entree
     integer :: ld, lf, md, mf                    ! numerotation locale
     integer :: numproc                           ! numero du proc qui travaille sur ce bloc (//)
     integer :: numproc_sb                        ! numero du proc dans le communicateur du superbloc (//)
     integer :: numproc_fluide                    ! numero du proc dans le communicateur fluide (//)
     integer :: numproc_solide                    ! numero du proc dans le communicateur solide (//)
     integer :: proc_ini                          ! numero du premier proc du bloc pere      (//)
     integer :: proc_fin                          ! numero du dernier proc du bloc pere      (//)
     integer, dimension(:,:), pointer :: l_charge ! indices (ld,lf) initiaux pour tous les procs (//)
     integer, dimension(:,:), pointer :: m_charge ! indices (md,mf) initiaux pour tous les procs (//)
     integer, dimension(:,:), pointer :: NbCoupes ! nombre de decoupes suivant chaque direction  (//)
     real*8 :: xin, xout                          ! taille du domaine en x
     real*8 :: yin, yout                          ! taille du domaine en y
     real*8 :: K_Riemann                          ! facteur de sécurité pour les pentes du solveur de Riemann (critere de positivite)
     real*8 :: delta_epsilon                      ! terme qui permet de rajouter de la diffusion dans la partie décentrée
     real*8 :: delta_epsilon_couche_lim           ! correction du epsilon pour l'annuler dans la couche limite visqueuse
     logical :: visqueux                          ! termes visqueux ?
     logical :: implicite                         ! schema implicite ?
     integer :: ordre_Eulerdirect                 ! ordre spatial pour le direct
     integer :: calcul_gradient                   ! méthode de calcul du gradient (GG, WLSQ)
     integer :: calcul_gradient_methode           ! inverse direct ou decomposition QR
     integer :: st_coin                           ! stencil pour décoder les coins
     integer :: st_grad                           ! stencil du calcul du gradient
     real*8  :: calcul_gradient_poids                    ! poids dans le calcul du gradient (=q de ||dx||^-q)
     integer :: eq_energie                        ! 0:equation sur l'energie, 1:isentropique, 2:isotherme
     real*8 :: entropie_ref                       ! entropie de reference (si isentropique)
     real*8 :: temperature_ref                    ! temperature de reference (si isotherme)
     integer :: correction_son
     
     type (STR_MAILLAGE) :: maillage  ! maillage, si on raffine, quelle face?, coef

!!! pour l'initialisation
     type (STR_cell), dimension(:), pointer :: tab_cl_ini
     character(len=90) :: type_initialisation

!!! Faces
     type (STR_FACE), dimension(4) :: face
!!! Coins
     type (STR_COIN), dimension(4) :: coin

!!! tableaux
     !! valeurs aux centres des mailles
     real*8, dimension(:,:,:), pointer :: U, U_n   ! variables Euler           (Fluide)
     real*8, dimension(:,:)  , pointer :: rho      ! densite                   (F,S)
     real*8, dimension(:,:)  , pointer :: u_x, u_y ! composantes de la vitesse (F)
     real*8, dimension(:,:)  , pointer :: E        ! Energie totale            (F)
     real*8, dimension(:,:)  , pointer :: y        ! fraction massique         (F)
     real*8, dimension(:,:)  , pointer :: p, p_n   ! pression                  (F)
     real*8, dimension(:,:)  , pointer :: c        ! vitesse du son            (F)
     real*8, dimension(:,:)  , pointer :: z        ! fraction volumique        (F)
     real*8, dimension(:,:)  , pointer :: rho1     ! densite fluide 1          (F)
     real*8, dimension(:,:)  , pointer :: rho2     ! densite fluide 2          (F)
     real*8, dimension(:,:)  , pointer :: zeros    ! zeros pour unifier le codage
     real*8, dimension(:,:)  , pointer :: mu       ! viscosite                 (F:visqueux)
     real*8, dimension(:,:)  , pointer :: k        ! conductivite thermique    (F:visqueux,S)
     real*8, dimension(:,:)  , pointer :: cv       ! capacite thermique a volume constant (F:visqueux, S)
     real*8, dimension(:,:)  , pointer :: cv1      ! capacite thermique a volume constant (F:visqueux multiespeces)
     real*8, dimension(:,:)  , pointer :: T        ! temperature               (F:visqueux,S)
     real*8, dimension(:,:)  , pointer :: T1, T2   ! temperature phases 1 et 2 (F:visqueux)
     real*8, dimension(:,:)  , pointer :: Htot     ! enthalpie totale          (F)
     real*8, dimension(:,:)  , pointer :: h        ! enthalpie                 (F)

     real*8, dimension(:,:)  , pointer :: dt_vv    ! dt lie au visqueux    (debug)
     real*8, dimension(:,:)  , pointer :: dt_vT    ! dt lie au visqueux    (debug)

     real*8, dimension(:,:)  , pointer :: gamma_eq ! coefficient adiabatique   (F:EOS)
     real*8, dimension(:,:)  , pointer :: gamma1   ! coefficient adiabatique du fluide1  (F:EOS multiespeces)
     real*8, dimension(:,:)  , pointer :: pi_eq    ! pression de reference     (F:EOS)
     real*8, dimension(:,:)  , pointer :: q_eq     ! energie interne de ref    (F:EOS)
     real*8, dimension(:,:,:), pointer :: c_esp    ! fraction massique par especes  (F)
     real*8, dimension(:,:,:), pointer :: c_ele    ! fraction massique par elements (F)
     real*8, dimension(:,:,:), pointer :: Hi       ! enthalpie                 (F:Multi-especes)
     real*8, dimension(:,:)  , pointer :: D        ! coeff diffusion           (F:Multi-especes)

     real*8, dimension(:,:)  , pointer :: aire_old ! aire au temps n
     real*8, dimension(:,:)  , pointer :: aire     ! aire de la cellule
     real*8, dimension(:,:)  , pointer :: aire_new ! aire au temps n+1
     real*8, dimension(:,:)  , pointer :: dm       ! aire*rho                  (F)
     real*8, dimension(:,:,:), pointer :: centre_cl ! centre de la cellule
     real*8, dimension(:,:,:), pointer :: Mpinv   ! pseudo-inverse pour les calculs de gradient

     type (STR_cell), dimension(:,:), pointer :: cell ! peut on s'en passer ?

     !! valeurs aux noeuds
     real*8, dimension(:,:,:), pointer :: coord_nd        ! coordonnees des noeuds
     real*8, dimension(:,:,:), pointer :: coord_nd_old    ! noeuds au temps t^n
     real*8, dimension(:,:,:), pointer :: coord_nd_01x01  ! noeuds sur le maillage cartesien unitaire

     !! valeurs aux centres des faces
     real*8, dimension(:,:)  , pointer :: a_l             ! vitesse du son Lagrangienne
     real*8, dimension(:,:)  , pointer :: a_m             !
     real*8, dimension(:,:)  , pointer :: Cm_l            ! pente C- du solveur de Riemann
     real*8, dimension(:,:)  , pointer :: Cp_l            ! pente C+ du solveur de Riemann
     real*8, dimension(:,:)  , pointer :: Cm_m            !
     real*8, dimension(:,:)  , pointer :: Cp_m            !
     real*8, dimension(:,:)  , pointer :: dx_dt           ! vitesse du maillage mobile
     real*8, dimension(:,:)  , pointer :: dy_dt           !
     real*8, dimension(:,:)  , pointer :: dl_l            ! longueur de la face
     real*8, dimension(:,:)  , pointer :: dl_m            !
     real*8, dimension(:,:,:), pointer :: n_dir_l         ! normale (unitaire) a la face
     real*8, dimension(:,:,:), pointer :: n_dir_m         !
     real*8, dimension(:,:,:), pointer :: F               ! flux numerique (pour le direct ou la partie transport)
     real*8, dimension(:,:,:), pointer :: Fg_acoustique   ! flux numerique pour la partie acoustique
     real*8, dimension(:,:,:), pointer :: Fd_acoustique   ! on en a deux a cause d'EUCCLHYD
     real*8, dimension(:,:,:,:), pointer :: dV_dU         ! matrice de passage pour l'implicite direct
     real*8, dimension(:,:,:), pointer :: metrique_grad_l ! metrique pour calculer les gradients en visqeux
     real*8, dimension(:,:,:), pointer :: metrique_grad_m ! metrique pour calculer les gradients en visqeux

     !! recontructions des variables (rho1,rho2,rho1eps1,rho2eps2,ux,uy,z)
     real*8, dimension(:,:,:), pointer :: Uprim
     real*8, dimension(:,:,:), pointer :: Uprim_W
     real*8, dimension(:,:,:), pointer :: Uprim_E
     real*8, dimension(:,:,:), pointer :: Uprim_S
     real*8, dimension(:,:,:), pointer :: Uprim_N

     !! gradient pour U-muscl version Nishikawa
     real*8, dimension(:,:,:), pointer :: gradvar

     !! debug flag c<0
     integer, dimension(:,:), pointer :: debug_c

     !! est ce que TVB ?
     integer :: TVB

     real*8, dimension(:,:,:), pointer :: FGg, FGd               ! flux Gallice pour le schema Cnimp e
     real*8, dimension(:,:,:), pointer :: FEgx, FEdx, FEgy, FEdy ! flux EUCCLHYD pour le schema Cnimp e
     real*8, dimension(:,:,:), pointer :: sm_CN, sm_CN_ein

     !! etat intermediaire pour les schémas runge kutta
     real*8, dimension(:,:,:), pointer :: U_rk

     type (variable_reconstruite)      :: ordre2          ! variables pour ordre2

     type(derivees_partielles), dimension(:,:), pointer :: deriv_p,deriv_T ! pour implicite Euler direct

     !! second membre
     real*8, dimension(:,:,:), pointer :: sm_acoustique   ! pour la partie acoustique (la taille change si on est en explicite ou en implicite)
     real*8, dimension(:,:,:), pointer :: sm              ! pour le transport ou pour le direct


!!! Pour le solveur MKL, structure CSR
     type (STR_CSR), pointer     :: mat_csr
     type (matrice_tridiagonale) :: matrice

!!! matrices, vecteur d'inconnues, second membre pour l'implicite
     integer :: solveur_implicite
     type (STR_MATRICE) :: mat, mat_scalaire
     real*8, dimension(:,:,:), pointer :: sm_scalaire

     real*8, dimension(:,:,:), pointer :: X_imp, X_old, X_0

     real*8, dimension(:,:,:,:), pointer :: grad_e, grad_f

     real*8, dimension(:,:,:), pointer :: n_dir_l_dual
     real*8, dimension(:,:,:), pointer :: n_dir_m_dual
     real*8, dimension(:,:), pointer :: dl_l_dual
     real*8, dimension(:,:), pointer :: dl_m_dual
     real*8, dimension(:,:), pointer :: aire_dual

  end type STR_BLOC

  type STR_ACCES_BLOC ! structure d'acces aux blocs
     type (STR_BLOC), pointer :: pointeur
  end type STR_ACCES_BLOC

  type STR_SUPER_BLOC
     integer :: nature           ! fluide ou solide
     integer :: nb_membre_input  ! nombre de bloc du fichier d entrees
     integer :: nb_membre        ! nombre de bloc au cours du calcul
     integer :: numproc    
     integer :: mpi_comm         ! communicateur MPI 
     integer :: nbprocs
     integer, dimension(:), pointer :: membre  ! numeros des bloc du super bloc
!     integer, dimension(:), pointer :: nbprocs ! nombre de procs par bloc (les blocs donnes en entree)
     logical :: periodique
     
     real*8  :: dt                       ! pas de temps
     integer :: solveur_implicite        ! type de solveur lineaire pour implicite
     integer :: solveur_Riemann          ! type de schema de Riemann pour l'acoustique
     logical :: egalite_pentes_solveur   ! egalité des pentes du solveur ?
     integer :: schema_temps             ! schema en temps pour le transport (Euler explicite ou RK2)
     type (STR_montee_cfl) :: montee_cfl ! parametres de la montee en cfl
     real*8  :: cfl                      ! coefficient cfl (peut etre modifie par la partie transport)
     real*8  :: cfl_ini                  ! coeff cfl initial (n'est pas modifie)
     real*8  :: K_riemann                ! facteur de sécurité pour les pentes du solveur de Riemann (critere de positivite)
     real*8 :: delta_epsilon             ! terme qui permet de rajouter de la diffusion dans la partie décentrée
     real*8 :: delta_epsilon_couche_lim  ! correction du epsilon pour l'annuler dans la couche limite visqueuse
     logical :: implicite, visqueux, capillaire
     logical :: gravite
     integer :: eq_energie               ! 0:equation sur l'energie, 1:isentropique, 2:isotherme
     real*8  :: entropie_ref=0.d0        ! entropie de reference (si isentropique)
     real*8  :: temperature_ref=0.d0     ! temperature de reference (si isotherme)

     type (process_iteratif) :: inv_up   ! implicite en (u,p) pour la partie acoustique
     type (process_iteratif) :: inv_E    ! implicite equation energie (termes visqueux)
                                         ! implicite T pour un bloc solide
     type (process_iteratif) :: sol_stat ! calcul stationnaire
     type (process_iteratif) :: preconvergence ! preconvergence du fluide

  end type STR_SUPER_BLOC

  type STR_ACCES_SUPER ! structure d'acces aux super-blocs
     type (STR_SUPER_BLOC), pointer :: pointeur
  end type STR_ACCES_SUPER

  type STR_TACHE ! definition du probleme
     integer :: nb_membre
     integer, dimension(:), pointer :: membre
     logical :: fluide=.false.   ! fluide ?

     real*8 :: dt                ! pas de temps
     real*8 :: tmax              ! temps maximal

     character(len=5) :: type_freqout ! type de frequence pour les sorties
                                      ! ("iter" ou "temp")
     real*8 :: freq_out               ! frenquence (temporelle) des sorties
     real*8, pointer :: t_imposes(:)  ! liste des temps imposes pour les sorties
     integer :: ind_timp              ! indice dans la liste des temps imposes

  end type STR_TACHE

!!! variables globales
  type (STR_ACCES_SUPER), dimension(:), allocatable, save :: acces_super_bloc
  type (STR_ACCES_BLOC), dimension(:), allocatable, save :: acces_bloc

  logical :: delta_forme

  type (STR_fluide), target :: fluide1, fluide2 ! parametres des deux fluides
  
end module m_struct
