module m_initialisation_fluide
  use m_struct
  use m_MPI
  use m_maillage
  use variables_globales
  use decodage

  implicit none

  public initialisation_fluide

  private initialisation_reprise, initialisation_par_melange
  private initialisation_sod_1d, initialisation_sod_2d
  private initialisation_sod_1d_trizone
  private initialisation_bulle_1d, initialisation_bulle_2d
  private initialisation_ellipse_2d
  private initialisation_bulle_murrone, initialisation_carre_2d
  private initialisation_bulle_choc_1d, initialisation_bulle_choc
  private initialisation_vitesse_lineaire, initialisation_vortex
  private initialisation_steady_discontinuity,initialisation_zalesak_disk,initialisation_circle
  private initialisation_etoile,initialisation_etoile2,initialisation_creneau
  private initialisation_Kothe_Rider
  private initialisation_Kelvin_Helmholtz_instability
  private initialisation_rabbit
contains

  subroutine initialisation_fluide(b, type_maillage, type_initialisation)
    use m_input
    use decodage
    use outils_pas_acoustique
    use outils_fluide

    implicit none

    type (STR_BLOC), pointer :: b
    character(len=90), intent(in) :: type_maillage, type_initialisation

    integer :: l, m
    integer :: ierreur
    real*8, dimension(:,:,:), pointer     :: U
    type (STR_cell), dimension(:,:), pointer :: cell

    call creation_maillage(b, type_maillage)
    if (b%maillage%move) b%aire_new = b%aire

    ! On lit le fichier meme en repise pour avoir les gamma et pi
    call read_donnees_initialisation(b)

    if ( index(type_initialisation,"restart") /= 0 ) then
       call initialisation_reprise(b, type_initialisation)
     !  call init_reprise_plus_decoupe(b, type_initialisation)
    else
       select case(type_initialisation)
       case("zone1")
          do m = b%md-nMF, b%mf+nMF
             do l = b%ld-nMF, b%lf+nMF
                b%cell(l,m) = b%tab_cl_ini(1)
             end do
          end do
       case("zone2")
          do m = b%md-nMF, b%mf+nMF
             do l = b%ld-nMF, b%lf+nMF
                b%cell(l,m) = b%tab_cl_ini(2)
             end do
          end do
       case("melange")
          call initialisation_par_melange(b, b%tab_cl_ini)
       case("sod_1d")
          call initialisation_sod_1d(b, b%tab_cl_ini)
       case("sod_1d_trizone")
          call initialisation_sod_1d_trizone(b, b%tab_cl_ini)
       case("sod_2d")
          call initialisation_sod_2d(b, b%tab_cl_ini, .5d0, .5d0)
       case("bulle_1d")
          call initialisation_bulle_1d(b, b%tab_cl_ini, .4d0, .6d0)
       case("bulle_2d")
          call initialisation_bulle_2d(b, b%tab_cl_ini)
       case("ellipse_2d")
          call initialisation_ellipse_2d(b, b%tab_cl_ini)
       case("double_bulle_2d")
          call initialisation_double_bulle_2d(b, b%tab_cl_ini)
       case("carre_2d")
          call initialisation_carre_2d(b, b%tab_cl_ini)
       case("bulle_murrone")
          call initialisation_bulle_murrone(b, b%tab_cl_ini)
       case("bulle_choc_1d")
          call initialisation_bulle_choc_1d(b, b%tab_cl_ini)
       case("bulle_choc")
          call initialisation_bulle_choc(b, b%tab_cl_ini)
       case("vitesse_lineaire")
          call initialisation_vitesse_lineaire(b, b%tab_cl_ini)
       case("vortex")
          call initialisation_vortex(b)
       case("vortex_diphasique")
          call initialisation_vortex_diphasique(b)
       case("couette")
          call initialisation_couette(b, b%tab_cl_ini)
       case("steady_discontinuity_advection")
          call initialisation_steady_discontinuity(b,b%tab_cl_ini)
       case("zalesak_disk_advection")
          call initialisation_zalesak_disk(b,b%tab_cl_ini)
       case("etoile")
          call initialisation_etoile(b,b%tab_cl_ini)
       case("etoile2")
          call initialisation_etoile2(b,b%tab_cl_ini)
       case("circle")
          call initialisation_circle(b,b%tab_cl_ini)
       case("creneau")
          call initialisation_creneau(b,b%tab_cl_ini)
       case("Kothe_Rider")
          call initialisation_Kothe_Rider(b,b%tab_cl_ini)
       case("kh_instability")
          call initialisation_Kelvin_Helmholtz_instability(b,b%tab_cl_ini)
       case("rabbit")
          call initialisation_rabbit(b,b%tab_cl_ini)
       case default
          print*,""
          print*," Type d'initialisation non definie", type_initialisation
          print*,""
          print*,"ERREUR :initialisation"
          call arret_code
       end select

!!! Remplissage des variables conservatives
       cell => b%cell
       U=>b%U

       select case (fluide1%type)
       case(EOS)
          ! pour avoir une programmation unifiee
          do m = b%md-nMF, b%mf+nMF
             do l = b%ld-nMF, b%lf+nMF
                allocate(cell(l,m)%Ci(1))
                cell(l,m)%Ci(:) = 1.d0
             end do
          end do
       end select

       do m = b%md-nMF, b%mf+nMF
          do l = b%ld-nMF, b%lf+nMF
             call encode_U(cell(l,m), U(:,l,m))
          end do
       end do
       
    end if


!!! on decode le vecteur des variables conservatives car on peut avoir besoin de
!!!  - rho pour calculer la vitesse de fusion
!!!  - p pour les conditions aux limites de type subsonique
!!!  - ...
    call decode_U(b%ld-nMF, b%lf+nMF, b%md-nMF, b%mf+nMF, b, ierreur)

    call save_Un_bloc(b)
    deallocate(b%tab_cl_ini)

    if ( index(type_initialisation,"restart") /= 0 ) then
    else
       deallocate(b%cell)
    end if

  end subroutine initialisation_fluide

  subroutine initialisation_reprise(b, type_ini)
    use m_temps
    use m_MPI
    implicit none

    type (STR_BLOC), pointer :: b
    character(len=90), intent(in) :: type_ini
    integer :: iunit=47
    character(len=90) :: fichier !, iter_c
    integer :: iter_read
    integer :: neq_read, ld_read, lf_read, md_read, mf_read
    integer :: ld, lf, md, mf
    real*8, dimension(:,:,:), allocatable :: U_global

    fichier = trim(path_output)//'restart/'//trim(type_ini)//'_cellules_fld.dat'
    open(unit=iunit, file=fichier, status='old', form='unformatted' )
    read(iunit) iter_read, temps
    read(iunit) neq_read, ld_read, lf_read, md_read, mf_read

    if ( iter_read /= iter_ini .and. iter_read /= -1) then
       print*, ""
       print*, "L'ITERATION DE DEPART N'EST PAS CELLE DU FICHIER DE REPRISE"
       print*, "iter_read", iter_read, "iter_ini", iter_ini
       print*, ""
       ! call arret_code
    end if
    if ( neq_read /= b%neq .or. ld_read /= b%ld_tot .or. lf_read /= b%lf_tot &
         & .or. md_read /= b%md_tot .or. mf_read /= b%mf_tot ) then
       print*, ""
       print*, "IL Y A UNE INCOMPATIBILITE ENTRE LES INDICES DU RESTART ET CEUX DES DONNEES"
       print*, neq_read, b%neq
       print*, ld_read, b%ld_tot, lf_read, b%lf_tot
       print*, md_read, b%md_tot, mf_read, b%mf_tot
       print*, ""
       call arret_code
    end if

    allocate( U_global(b%neq, b%ld_tot-nMF:b%lf_tot+nMF, b%md_tot-nMF:b%mf_tot+nMF) )

    ! Todo : enlever les mailles fictives des reprises
    read(iunit) U_global(1:b%neq, b%ld_tot-nMF:b%lf_tot+nMF, b%md_tot-nMF:b%mf_tot+nMF)
    close(iunit)

    ld = b%l_charge(1, numproc); lf = b%l_charge(2, numproc)
    md = b%m_charge(1, numproc); mf = b%m_charge(2, numproc)
    b%U(:, b%ld-nMF:b%lf+nMF, b%md-nMF:b%mf+nMF) = U_global(:, ld-nMF:lf+nMF, md-nMF:mf+nMF)

    deallocate(U_global)
  end subroutine initialisation_reprise

  subroutine init_reprise_plus_decoupe(b, type_ini)
    use m_temps
    use m_MPI
    implicit none

    type (STR_BLOC), pointer :: b
    character(len=90), intent(in) :: type_ini
    integer :: iunit=47
    character(len=90) :: fichier !, iter_c
    integer :: iter_read
    integer :: neq_read, ld_read, lf_read, md_read, mf_read
    integer :: ld, lf, md, mf
    real*8, dimension(:,:,:), allocatable :: U_global

    fichier = trim(path_output)//'restart/'//trim(type_ini)//'_cellules_fld.dat'
    open(unit=iunit, file=fichier, status='old', form='unformatted' )
    read(iunit) iter_read, temps
    read(iunit) neq_read, ld_read, lf_read, md_read, mf_read

    if ( iter_read /= iter_ini .and. iter_read /= -1) then
       print*, ""
       print*, "L'ITERATION DE DEPART N'EST PAS CELLE DU FICHIER DE REPRISE"
       print*, "iter_read", iter_read, "iter_ini", iter_ini
       print*, ""
       ! call arret_code
    end if
    if ( neq_read /= b%neq .or. ld_read /= b%ld_tot .or. lf_read /= b%lf_tot &
         & .or. md_read /= b%md_tot .or. mf_read /= b%mf_tot ) then
       print*, ""
       print*, "IL Y A UNE INCOMPATIBILITE ENTRE LES INDICES DU RESTART ET CEUX DES DONNEES"
       print*, neq_read, b%neq
       print*, ld_read, b%ld_tot, lf_read, b%lf_tot
       print*, md_read, b%md_tot, mf_read, b%mf_tot
       print*, ""
       !call arret_code
    end if

    allocate( U_global(b%neq, ld_read-nMF:lf_read+nMF, md_read-nMF:mf_read+nMF) )

    read(iunit) U_global(1:b%neq, ld_read-nMF:lf_read+nMF, md_read-nMF:mf_read+nMF)
    close(iunit)

    ld = b%l_charge(1, numproc); lf = b%l_charge(2, numproc)
    md = b%m_charge(1, numproc); mf = b%m_charge(2, numproc)

    b%U(:, b%ld-nMF:b%lf+nMF, b%md-nMF:b%mf+nMF) = U_global(:, ld-nMF:lf+nMF, md-nMF:mf+nMF)

    deallocate(U_global)

!!! on initialise rho car on peut en avoir besoin pour calculer la vitesse de fusion
    b%rho(:,:) = b%U(1,:,:)

    call redresse_condition_ini(b)
  end subroutine init_reprise_plus_decoupe

  subroutine redresse_condition_ini(b)
    type (STR_BLOC), pointer :: b

    integer :: l,m

    ! Inclusion des mailles fictives via nMF ?


!!$    do m = b%md-2, b%mf+2
!!$       do l = b%ld-2, b%lf+2
!!$
!!$          !! on met z=0 ou z=1
!!$          if (abs(b%U(iz,l,m)-1.d0)<1.d-3) then
!!$             b%U(iz,l,m)=1.d0
!!$          end if
!!$          if (abs(b%U(iz,l,m))<1.d-3) then
!!$             b%U(iz,l,m)=0.d0
!!$          end if
!!$
!!$          !! la maille de melange
!!$          if (b%U(iz,l,m)/=1.d0 .and. b%U(iz,l,m)/=0.d0) then
!!$             b%U(iz,l,m) = 1.d0
!!$             b%U(:,l,m) = b%U(:,l-1,m)
!!$          end if
!!$
!!$       end do
!!$    end do

    if (b%lf /= b%ld) call arret_code
    do m = b%md-nMF, b%mf+nMF
       do l = b%ld, b%lf

          !! on met z=0 ou z=1
          if (abs(b%U(iz,l,m)-1.d0)<1.d-3) then
             b%U(iz,l,m)=1.d0
          end if
          if (abs(b%U(iz,l,m))<1.d-3) then
             b%U(iz,l,m)=0.d0
          end if

          !! la maille de melange
          if (b%U(iz,l,m)/=1.d0 .and. b%U(iz,l,m)/=0.d0) then
             if (b%U(iz,l,m)<.5d0) then
                b%U(iz,l,m) = 0.d0
                b%U(:,l,m) = b%U(:,l,m+1)
             else
                b%U(iz,l,m) = 1.d0
                b%U(:,l,m) = b%U(:,l,m-1)
             end if
          end if

       end do
    end do
    b%U(:,b%ld-2,:) = b%U(:,b%ld,:)
    b%U(:,b%ld-1,:) = b%U(:,b%ld,:)
    b%U(:,b%lf+1,:) = b%U(:,b%lf,:)
    b%U(:,b%lf+2,:) = b%U(:,b%lf,:)


  end subroutine redresse_condition_ini

!!! Routines des differents types d'initialisation pour le fluide

  subroutine initialisation_par_melange(b, tab_cl_ini)
    use decodage
    use decodage

    type (STR_BLOC), pointer :: b
    type (STR_cell), dimension(3), target, intent(in) :: tab_cl_ini

    type (STR_cell), pointer :: cell1, cell2
    type (STR_cell) :: melange
    real*8 :: z1, z2, eps

    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)

    z1 = cell1%z
    z2 = cell2%z

    if (z1==0.5d0 .and. z2==0.5d0) then
       melange%z = 0.5d0
    else
       print*, "erreur initialisation melange"
       melange%z = 0.5d0
       z1 = .5d0; z2 = .5d0
       cell1%z = .5d0
       cell2%z = .5d0
       call arret_code
    end if

    melange%rho1 = cell1%rho
    melange%rho2 = cell2%rho
    melange%y = (cell1%rho*z1)/(cell1%rho*z1 + cell2%rho*z2)
    melange%rho = cell1%rho * z1 + cell2%rho * z2
    melange%u =  cell1%u * z1 + cell2%u * z2
    melange%v =  cell1%v * z1 + cell2%v * z2
    melange%p =  cell1%p * z1 + cell2%p * z2

    ! on determine gamma, pi, q, mu, cv, K
    call grandeurs_melange(fluide1, fluide2, melange%z, melange%y, melange%fluide)

    ! energie interne
    eps = EOS_eps_from_p_rho(melange%fluide%EOS, melange%p, melange%rho)
    melange%h = eps + melange%p/melange%rho

    melange%c=sqrt(melange%fluide%EOS%gamma*(melange%P+melange%fluide%EOS%pi)/melange%rho)

    b%cell = melange

  end subroutine initialisation_par_melange

  subroutine initialisation_sod_1d(b, tab_cl_ini)
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(2), target, intent(in) :: tab_cl_ini

    integer :: l, m, ld, lf, md, mf
    real*8 :: x_int
    real*8, dimension(:,:,:), pointer :: centre_cl
    character :: sens_interface
    type (STR_cell), pointer :: cell1, cell2
    type (STR_cell), dimension(:,:), pointer :: cell

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)

    ! lecture des donnees
    open(unit=20,file=trim(path_input)//"parametres_initialisation.dat")
    read(20,*) sens_interface, x_int   !position de l'interface
    close(20)

    if (sens_interface == "y")  then
       ! fluide 2
       ! --------
       ! fluide 1
       do m = md-nMF, mf+nMF
          do l = ld-nMF, lf+nMF
             if (centre_cl(2,l,m)<=x_int) then ! fluide 1
                cell(l,m) = cell1
             else
                cell(l,m) = cell2
             end if
          end do
       end do
    else if (sens_interface == "x")  then
       ! fluide 1  |  fluide 2
       do m = md-nMF, mf+nMF
          do l = ld-nMF, lf+nMF
             if (centre_cl(1,l,m)<=x_int) then ! fluide 1
                cell(l,m) = cell1
             else
                cell(l,m) = cell2
             end if
          end do
       end do
    else
       print*, "orientation de l'interface entre les deux fluides inconnue"
       print*, sens_interface
       call arret_code
    end if

  end subroutine initialisation_sod_1d

    subroutine initialisation_sod_1d_trizone(b, tab_cl_ini)
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(3), target, intent(in) :: tab_cl_ini

    integer :: l, m, ld, lf, md, mf
    real*8 :: x_1, x_2
    real*8, dimension(:,:,:), pointer :: centre_cl
    character :: sens_interface
    type (STR_cell), pointer :: cell1, cell2, cell3
    type (STR_cell), dimension(:,:), pointer :: cell

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)
    cell3 => tab_cl_ini(3)

    ! lecture des donnees
    open(unit=20,file=trim(path_input)//"parametres_initialisation.dat")
    read(20,*) sens_interface, x_1 , x_2   !position de l'interface
    close(20)

    if (sens_interface == "y")  then
       ! fluide 3
       ! --------
       ! fluide 2
       ! --------
       ! fluide 1
       do m = md-nMF, mf+nMF
          do l = ld-nMF, lf+nMF
             if (centre_cl(2,l,m)<=x_1) then ! fluide 1
                cell(l,m) = cell1
             elseif(centre_cl(2,l,m)>x_1 .and.centre_cl(2,l,m)<=x_2)then
                cell(l,m) = cell2
             else
                cell(l,m) = cell3
             end if
          end do
       end do
    else if (sens_interface == "x")  then
       ! fluide 1  |  fluide 2 | fluide 3
       do m = md-nMF, mf+nMF
          do l = ld-nMF, lf+nMF
             if (centre_cl(1,l,m)<=x_1) then ! fluide 1
                cell(l,m) = cell1
             elseif(centre_cl(1,l,m)>x_1 .and.centre_cl(1,l,m)<=x_2)then
                cell(l,m) = cell2
             else
                cell(l,m) = cell3
             end if
          end do
       end do
    else
       print*, "orientation de l'interface entre les deux fluides inconnue"
       print*, sens_interface
       call arret_code
    end if

  end subroutine initialisation_sod_1d_trizone


  subroutine initialisation_sod_2d(b, tab_cl_ini, x0, y0)
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(4), target, intent(in) :: tab_cl_ini
    real*8, intent(in) :: x0, y0

    integer :: l, m, ld, lf, md, mf
    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), pointer :: cell1, cell2, cell3, cell4
    type (STR_cell), dimension(:,:), pointer :: cell

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)
    cell3 => tab_cl_ini(3)
    cell4 => tab_cl_ini(4)

    ! repartition des zones
    ! | 3 4 |
    ! | 1 2 |

    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          if (centre_cl(1,l,m)<=x0 .and. centre_cl(2,l,m)<=y0) then ! fluide 1
             cell(l,m) = cell1
          else if (centre_cl(1,l,m)>x0 .and. centre_cl(2,l,m)<=y0) then ! fluide 2
             cell(l,m) = cell2
          else if (centre_cl(1,l,m)<=x0 .and. centre_cl(2,l,m)>y0) then ! fluide 3
             cell(l,m) = cell3
          else if (centre_cl(1,l,m)>x0 .and. centre_cl(2,l,m)>y0) then ! fluide 4
             cell(l,m) = cell4
          end if
       end do
    end do

  end subroutine initialisation_sod_2d

  ! initialisation d'une bulle de fluide
  ! tab_cl_ini(2) pour x1<x<x2
  ! tab_cl_ini(1) autour
  subroutine initialisation_bulle_1d(b, tab_cl_ini, x1, x2)
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(2), target, intent(in) :: tab_cl_ini
    real*8, intent(in) :: x1, x2

    integer :: l, m, ld, lf, md, mf
    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), pointer :: cell1, cell2
    type (STR_cell), dimension(:,:), pointer :: cell

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)

    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          if (centre_cl(1,l,m)>=x1 .and. centre_cl(1,l,m)<=x2) then ! bulle = fluide2
             cell(l,m) = cell2
          else
             cell(l,m) = cell1
          end if
       end do
    end do

  end subroutine initialisation_bulle_1d

   ! initialition d'une bulle en 2 dimensions
  ! fluide 1 : fluide entourant la bulle
  ! fluide 2 : bulle
  subroutine initialisation_bulle_2d(b, tab_cl_ini)
    use m_maillage
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(2), target, intent(in) :: tab_cl_ini

    integer :: l, m, ld, lf, md, mf
    real*8 :: rayon, centre(2)
    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), pointer :: cell1, cell2
    type (STR_cell), dimension(:,:), pointer :: cell

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)

    open(unit=20,file=trim(path_input)//"parametres_initialisation.dat")
    read(20,*) centre
    read(20,*) rayon
    close(20)

    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          if (distance(centre_cl(:,l,m), centre(:))<=rayon) then
             cell(l,m) = cell1
          else
             cell(l,m) = cell2
          endif
       enddo
    enddo

  end subroutine initialisation_bulle_2d

   ! initialition d'une bulle en 2 dimensions
  ! fluide 1 : fluide entourant la bulle
  ! fluide 2 : bulle
  subroutine initialisation_ellipse_2d(b, tab_cl_ini)
    use m_maillage
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(2), target, intent(in) :: tab_cl_ini

    integer :: l, m, ld, lf, md, mf
    real*8 :: pa(2), centre(2), e
    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), pointer :: cell1, cell2
    type (STR_cell), dimension(:,:), pointer :: cell

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)

    open(unit=20,file=trim(path_input)//"parametres_initialisation.dat")
    read(20,*) centre
    read(20,*) pa
    close(20)

    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          e = ((centre_cl(1,l,m) - centre(1))/pa(1))**2 + ((centre_cl(2,l,m) - centre(2))/pa(2))**2
          if (e <= 1.d0) then
             cell(l,m) = cell1
          else
             cell(l,m) = cell2
          endif
       enddo
    enddo

  end subroutine initialisation_ellipse_2d


   ! initialition d'une bulle en 2 dimensions
  ! fluide 1 : fluide entourant la bulle
  ! fluide 2 : bulle 1
  ! fluide 3 : bulle 2
  subroutine initialisation_double_bulle_2d(b, tab_cl_ini)
    use m_maillage
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(3), target, intent(in) :: tab_cl_ini

    integer :: l, m, ld, lf, md, mf
    real*8 :: rayon, centre(4)
    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), pointer :: cell1, cell2, cell3
    type (STR_cell), dimension(:,:), pointer :: cell

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)
    cell3 => tab_cl_ini(3)

    open(unit=20,file=trim(path_input)//"parametres_initialisation.dat")
    read(20,*) centre
    read(20,*) rayon
    close(20)

    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          if (distance(centre_cl(:,l,m), centre(1:2))<=rayon) then
             cell(l,m) = cell2
          else if(distance(centre_cl(:,l,m), centre(3:4))<=rayon)then
             cell(l,m) = cell3
          else
             cell(l,m) = cell1
          endif
       enddo
    enddo

  end subroutine initialisation_double_bulle_2d


   ! initialition d'une bulle en 2 dimensions
  ! fluide 1 : fluide entourant la bulle
  ! fluide 2 : bulle
  subroutine initialisation_carre_2d(b, tab_cl_ini)
    use m_maillage
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(2), target, intent(in) :: tab_cl_ini

    integer :: l, m, ld, lf, md, mf
    real*8 :: rayon, centre(2), norme
    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), pointer :: cell1, cell2
    type (STR_cell), dimension(:,:), pointer :: cell

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)

    open(unit=20,file=trim(path_input)//"parametres_initialisation.dat")
    read(20,*) centre
    read(20,*) rayon
    close(20)

    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          norme = max(abs(centre_cl(1,l,m)- centre(1)),abs(centre_cl(2,l,m)- centre(2)))
          if (norme<=rayon) then
             cell(l,m) = cell1
          else
             cell(l,m) = cell2
          endif
       enddo
    enddo

  end subroutine initialisation_carre_2d

  subroutine initialisation_double_carre_2d(b, tab_cl_ini)
    use m_maillage
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(2), target, intent(in) :: tab_cl_ini

    integer :: l, m, ld, lf, md, mf
    real*8 :: rayon, centre(4), norme1, norme2
    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), pointer :: cell1, cell2
    type (STR_cell), dimension(:,:), pointer :: cell

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)

    open(unit=20,file=trim(path_input)//"parametres_initialisation.dat")
    read(20,*) centre
    read(20,*) rayon
    close(20)

    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          norme1 = max(abs(centre_cl(1,l,m)- centre(1)),abs(centre_cl(2,l,m)- centre(2)))
          norme2 = max(abs(centre_cl(1,l,m)- centre(3)),abs(centre_cl(2,l,m)- centre(4)))

          if (distance(centre_cl(:,l,m), centre(1:2))<=rayon .or. &
            distance(centre_cl(:,l,m), centre(3:4))<=rayon) then
             cell(l,m) = cell1
          else
             cell(l,m) = cell2
          endif
       enddo
    enddo

  end subroutine initialisation_double_carre_2d


  subroutine initialisation_Kelvin_Helmholtz_instability(b, tab_cl_ini)
    use m_maillage
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(2), target, intent(in) :: tab_cl_ini

    integer :: l, m, ld, lf, md, mf
    real*8 :: alpha, beta, k, f, pos_interface
    character :: sens_interface
    real*8 :: pi

    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), pointer :: cell1, cell2
    type (STR_cell), dimension(:,:), pointer :: cell

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)

    open(unit=20,file=trim(path_input)//"parametres_initialisation.dat")
    read(20,*) sens_interface, pos_interface
    read(20,*) alpha, beta, k
    close(20)

    pi = 4.d0*atan(1.d0)
    if(sens_interface == 'x')then
      do m = md-nMF, mf+nMF
         do l = ld-nMF, lf+nMF
           if(centre_cl(2,l,m)>=alpha .and. centre_cl(2,l,m) <= beta)then
             f = pos_interface - k*sin(pi*((centre_cl(2,l,m)-alpha)/(beta-alpha)))
           else 
             f = pos_interface
           end if

           if (centre_cl(1,l,m)<=f) then
             cell(l,m) = cell2
           else
             cell(l,m) = cell1
           endif
         enddo
      enddo
    else if(sens_interface == 'y')then
      do m = md-nMF, mf+nMF
         do l = ld-nMF, lf+nMF
           if(centre_cl(1,l,m)>=alpha .and. centre_cl(1,l,m) <= beta)then
             f = centre_cl(2,l,m) - pos_interface + k*sin(pi*((centre_cl(1,l,m)-alpha)/(beta-alpha)))
           else 
             f = centre_cl(2,l,m) - pos_interface
           end if

           if (f > 0.d0) then
             cell(l,m) = cell1
           else
             cell(l,m) = cell2
           endif
         enddo
      enddo
    end if

  end subroutine initialisation_Kelvin_Helmholtz_instability

 subroutine initialisation_bulle_murrone(b, tab_cl_ini)
   use m_maillage

   type (STR_BLOC), intent(inout) :: b
   type (STR_cell), dimension(2), target, intent(in) :: tab_cl_ini

   integer :: l, m, ld, lf, md, mf
   real*8 :: rayon, centre(2), dh
   real*8, dimension(:,:,:), pointer :: centre_cl
   type (STR_cell), pointer :: cell1, cell2
   type (STR_cell), dimension(:,:), pointer :: cell

!!!!!!!
    real*8 :: gamma1, pi1
    real*8 :: gamma2, pi2

    gamma1 = fluide1%EOS%gamma; pi1 = fluide1%EOS%pi
    gamma2 = fluide2%EOS%gamma; pi2 = fluide2%EOS%pi
!!!!!!!

   ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
   centre_cl => b%centre_cl
   cell => b%cell
   cell1 => tab_cl_ini(1)
   cell2 => tab_cl_ini(2)

   open(unit=20,file=trim(path_input)//"parametres_initialisation.dat")
   read(20,*) centre
   read(20,*) rayon
   close(20)

   do m = md-nMF, mf+nMF
      do l = ld-nMF, lf+nMF
         if (distance(centre_cl(:,l,m), centre(:))<=rayon) then
            cell(l,m) = cell2

            dh = b%centre_cl(2,l,m) - centre(2)
            cell(l,m)%p = cell(l,m)%p + cell(l,m)%rho * acc_pesanteur(2) * dh
            cell(l,m)%h = gamma2 * (cell(l,m)%p +  pi2) / (cell(l,m)%rho2*(gamma2-1.d0))
         else
            cell(l,m) = cell1

            dh = b%centre_cl(2,l,m) - centre(2)
            cell(l,m)%p = cell(l,m)%p + cell(l,m)%rho * acc_pesanteur(2) * dh
            cell(l,m)%h = gamma1 * (cell(l,m)%p +  pi1) / (cell(l,m)%rho1*(gamma1-1.d0))
         endif
      enddo
   enddo

  end subroutine initialisation_bulle_murrone

  subroutine initialisation_bulle_choc_1d(b, tab_cl_ini)
    type (STR_BLOC), pointer :: b

    type (STR_cell), dimension(3), target, intent(in) :: tab_cl_ini

    real*8 :: x12, x23, x32
    ! repartition des zones
    !       1       |    2    |  3  |    2
    ! Liquide(choc) | Liquide | gaz | Liquide

    integer :: l, m, ld, lf, md, mf
    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), pointer :: cell1, cell2, cell3
    type (STR_cell), dimension(:,:), pointer :: cell

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)
    cell3 => tab_cl_ini(3)

    open(unit=20,file=trim(path_input)//"parametres_initialisation.dat")
    read(20,*) x12
    read(20,*) x23
    read(20,*) x32
    close(20)

    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          if (centre_cl(1,l,m)<=x12) then
             cell(l,m) = cell1
          else if (x12<centre_cl(1,l,m) .and. centre_cl(1,l,m)<x23) then
             cell(l,m) = cell2
          else if (x23<=centre_cl(1,l,m) .and. centre_cl(1,l,m)<=x32) then
             cell(l,m) = cell3
          else if (x32<centre_cl(1,l,m)) then
             cell(l,m) = cell2
          else
             print*, "erreur dans l'initialisation"
             call arret_code
          end if
       end do
    end do
  end subroutine initialisation_bulle_choc_1d

  subroutine initialisation_bulle_choc(b, tab_cl_ini)
    use m_maillage
    implicit none

    type (STR_BLOC), pointer :: b
    type (STR_cell), dimension(3), target, intent(in) :: tab_cl_ini

    real*8 :: centre(2), rayon, position_choc
    character(len=4) :: sens_choc

    integer :: l, m, ld, lf, md, mf
    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), pointer :: cell1, cell2, cell3
    type (STR_cell), dimension(:,:), pointer :: cell


    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)
    cell3 => tab_cl_ini(3)

    open(unit=20,file=trim(path_input)//"parametres_initialisation.dat")
    read(20,*) centre
    read(20,*) rayon
    read(20,*) position_choc
    read(20,*) sens_choc
    close(20)

    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          if (distance(centre_cl(:,l,m), centre(:))<=rayon) then
             cell(l,m) = cell3
          else if (centre_cl(1,l,m)<=position_choc .and. sens_choc=="g->d" .or. &
               ( centre_cl(1,l,m)>=position_choc .and. sens_choc=="g<-d") ) then
             ! else if (centre_cl(2,l,m)<=position_choc .and. sens_choc=="g->d" .or. &
             !      ( centre_cl(2,l,m)>=position_choc .and. sens_choc=="g<-d") ) then
             cell(l,m) = cell1
          else
             cell(l,m) = cell2
          endif
       enddo
    enddo

    if (sens_choc/="g->d" .and. sens_choc/="g<-d") then
       print*, "Erreur initialisation_bulle_choc : sens de propagation du choc non defini", sens_choc
       call arret_code
    end if

  end subroutine initialisation_bulle_choc

  subroutine initialisation_vitesse_lineaire(b, tab_cl_ini)
    use m_maillage

    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(1), target, intent(in) :: tab_cl_ini

    integer :: l, m, ld, lf, md, mf
    real*8 :: u_bas, u_haut, pt(2)
    type (STR_cell), dimension(:,:), pointer :: cell

    cell => b%cell
    cell = tab_cl_ini(1)

    u_bas = 50.d0
    u_haut = 0.d0
    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf

    pt = (/b%xout, 0.5d0*(b%yin+b%yout)/)
    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          cell(l,m)%v = u_bas + (u_haut-u_bas)/(b%xout - b%xin) * distance(b%centre_cl(:,l,m), pt)
          cell(l,m)%u=0.d0
       end do
    end do

  end subroutine initialisation_vitesse_lineaire

  subroutine initialisation_vortex(b)
    implicit none

    type (STR_BLOC), intent(inout) :: b
    integer :: l, m, ld, lf, md, mf
    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), dimension(:,:), pointer :: cell
    real*8 :: pi_

!!!!!!!
    real*8 :: gamma1, pi1
    real*8 :: gamma2, pi2

    gamma1 = fluide1%EOS%gamma; pi1 = fluide1%EOS%pi
    gamma2 = fluide2%EOS%gamma; pi2 = fluide2%EOS%pi
!!!!!!!

    pi_ = acos(-1.d0)
    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell

    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          cell(l,m)%rho1 = 1.d0-.5d0*tanh(centre_cl(2,l,m)-.5d0)
          cell(l,m)%rho2 = 0.d0
          cell(l,m)%rho  = cell(l,m)%rho1
          cell(l,m)%z = 1.d0
          cell(l,m)%u = 2.d0*sin(pi_*centre_cl(1,l,m))**2 &
               & * sin(pi_*centre_cl(2,l,m))*cos(pi_*centre_cl(2,l,m))
          cell(l,m)%v = -2.d0*sin(pi_*centre_cl(1,l,m))*cos(pi_*centre_cl(1,l,m)) &
               & * sin(pi_*centre_cl(2,l,m))**2
          cell(l,m)%p = 1000.d0

          cell(l,m)%h = gamma1 * (cell(l,m)%p +  pi1) / ( cell(l,m)%rho*(gamma1-1.d0) )
       end do
    end do

  end subroutine initialisation_vortex

  subroutine initialisation_vortex_diphasique(b)
    use decodage
    implicit none

    type (STR_BLOC), intent(inout) :: b
    integer :: l, m, ld, lf, md, mf
    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), dimension(:,:), pointer :: cell
    real*8 :: pi_, centre(2), rayon, eps
    type(STR_fluide) :: melange

!!!!!!!
    real*8 :: gamma1, pi1
    real*8 :: gamma2, pi2

    gamma1 = fluide1%EOS%gamma; pi1 = fluide1%EOS%pi
    gamma2 = fluide2%EOS%gamma; pi2 = fluide2%EOS%pi
!!!!!!!

    pi_ = acos(-1.d0)
    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell

    centre=(/0.5d0, 0.25d0/)
    rayon=0.1d0

    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          if (distance(centre_cl(:,l,m), centre(:))<=rayon) then
             cell(l,m)%rho1 = 1.d0
             cell(l,m)%rho2 = 0.d0
             cell(l,m)%z = 1.d0
          else
             cell(l,m)%rho1 = 0.d0
             cell(l,m)%rho2 = 10.d0
             cell(l,m)%z = 0.d0
          end if
          cell(l,m)%u = 2.d0*sin(pi_*centre_cl(1,l,m))**2 &
               & * sin(pi_*centre_cl(2,l,m))*cos(pi_*centre_cl(2,l,m))
          cell(l,m)%v = -2.d0*sin(pi_*centre_cl(1,l,m))*cos(pi_*centre_cl(1,l,m)) &
               & * sin(pi_*centre_cl(2,l,m))**2
          cell(l,m)%p = 1000.d0

          cell(l,m)%rho = cell(l,m)%z*cell(l,m)%rho1 + (1.d0-cell(l,m)%z)*cell(l,m)%rho2
          cell(l,m)%y = cell(l,m)%rho1*cell(l,m)%z/cell(l,m)%rho
          call grandeurs_melange(fluide1, fluide2, cell(l,m)%z, cell(l,m)%y, melange)
          eps = EOS_eps_from_p_rho(melange%EOS, cell(l,m)%p, cell(l,m)%rho)
          cell(l,m)%h = eps + cell(l,m)%p/cell(l,m)%rho
       end do
    end do
  end subroutine initialisation_vortex_diphasique


  subroutine initialisation_couette(b, tab_cl_ini)
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(2), target, intent(in) :: tab_cl_ini

    integer :: l, m, ld, lf, md, mf
    type (STR_cell), pointer :: cell1, cell2
    type (STR_cell), dimension(:,:), pointer :: cell
    real*8 :: xd, Td, Tf, T

!!!!!!!
    real*8 :: gamma1, pi1
    real*8 :: gamma2, pi2

    gamma1 = fluide1%EOS%gamma; pi1 = fluide1%EOS%pi
    gamma2 = fluide2%EOS%gamma; pi2 = fluide2%EOS%pi
!!!!!!!

    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)

    xd = b%yout - b%yin
    Td = 4000.d0
    Tf = 933.d0

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          cell(l,m)%p = fluide2%EOS%Cv*(gamma2-1.d0)*cell2%rho*xd*(Td-Tf)/(xd*log(Td/Tf))
          T = (Td-Tf)/xd * b%centre_cl(2,l,m) + Tf

          cell(l,m)%rho1 = 0.d0
          cell(l,m)%rho2 = cell(l,m)%p/(fluide2%EOS%Cv*(gamma2-1.d0)*T)
          cell(l,m)%z = 0.d0
          cell(l,m)%u = 0.d0
          cell(l,m)%v = 0.d0

          cell(l,m)%h = gamma2 * (cell(l,m)%p +  pi2) /(cell(l,m)%rho2*(gamma2-1.d0))
       end do
    end do
  end subroutine initialisation_couette

  subroutine initialisation_etoile(b, tab_cl_ini)
    use m_maillage
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(2), target, intent(in) :: tab_cl_ini

    integer :: l, m, ld, lf, md, mf
    real*8 :: rayon, centre(2),angle,frequence,r,taille
    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), pointer :: cell1, cell2
    type (STR_cell), dimension(:,:), pointer :: cell

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)

    open(unit=20,file=trim(path_input)//"parametres_initialisation.dat")
    read(20,*) centre
    read(20,*) rayon
    read(20,*) frequence
    read(20,*) taille
    close(20)

    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          r=distance(centre_cl(:,l,m), centre(:))
          if(centre_cl(2,l,m)>=0)then
             angle = dacos(centre_cl(1,l,m)/r)
          else
             angle = -dacos(centre_cl(1,l,m)/r)
          end if
          if (r<=rayon+taille*cos(frequence*angle)) then
             cell(l,m) = cell2
          else
             cell(l,m) = cell1
          endif
       enddo
    enddo

  end subroutine initialisation_etoile

  subroutine initialisation_etoile2(b, tab_cl_ini)
    use m_maillage
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(2), target, intent(in) :: tab_cl_ini

    integer :: l, m, ld, lf, md, mf
    real*8 :: absolue
    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), pointer :: cell1, cell2
    type (STR_cell), dimension(:,:), pointer :: cell

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)


    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          absolue = abs(centre_cl(1,l,m)-0.5d0)
          if ((centre_cl(2,l,m)>(1.d0/3.d0-absolue)) .and. (centre_cl(2,l,m)<(2.d0/3.d0-4.d0*absolue) )) then
             cell(l,m) = cell2
          else if ((centre_cl(2,l,m)>(1.d0/3.d0+absolue)) .and. (centre_cl(2,l,m)<0.5d0)) then
             cell(l,m) = cell2
          else
             cell(l,m) = cell1
          end if
       end do
    end do
  end subroutine initialisation_etoile2

  subroutine initialisation_steady_discontinuity(b, tab_cl_ini)
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(2), target, intent(in) :: tab_cl_ini

    integer :: l, m, ld, lf, md, mf
    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), pointer :: cell1, cell2
    type (STR_cell), dimension(:,:), pointer :: cell

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)

    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          if (centre_cl(2,l,m)<=centre_cl(1,l,m)*0.5d0) then
             cell(l,m) = cell1
          else
             cell(l,m) = cell2
          endif
       enddo
    enddo

  end subroutine initialisation_steady_discontinuity

  subroutine initialisation_zalesak_disk(b, tab_cl_ini)
    use m_maillage
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(2), target, intent(in) :: tab_cl_ini

    integer :: l, m, ld, lf, md, mf
    real*8 :: rayon, centre(2),rayon_in,epaisseur
    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), pointer :: cell1, cell2
    type (STR_cell), dimension(:,:), pointer :: cell

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)

    open(unit=20,file=trim(path_input)//"parametres_initialisation.dat")
    read(20,*) centre
    read(20,*) rayon
    read(20,*) epaisseur
    read(20,*) rayon_in
    close(20)

    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          if (distance(centre_cl(:,l,m), centre(:))<=rayon .and. (abs(centre(1)-centre_cl(1,l,m))>=0.5d0*epaisseur .or. centre_cl(2,l,m)-centre(2) >=rayon_in )) then
             cell(l,m) = cell1
          else
             cell(l,m) = cell2
          endif
       enddo
    enddo
    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          cell(l,m)%u =0.5d0-centre_cl(2,l,m)
          cell(l,m)%v =centre_cl(1,l,m)-0.5d0
       end do
    end do

  end subroutine initialisation_zalesak_disk

  subroutine initialisation_rabbit(b,tab_cl_ini)
    use m_maillage
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(2), target, intent(in) :: tab_cl_ini

    integer :: l, m, ld, lf, md, mf, i, j, n, l1, m1
    type (STR_cell), pointer :: cell1, cell2
    type (STR_cell), dimension(:,:), pointer :: cell
    character(1),dimension(1536) :: ligne
    integer, dimension(768,768) :: tab_in

    n = 768
    open(unit=20,file=trim(path_input)//"Rabbit.dat")
    do j=1,n
      read(20,'(1536a1)') ligne
      do i=1,n
          if(ligne(2*i-1)=='0') tab_in(i,j) = 0
          if(ligne(2*i-1)=='1') tab_in(i,j) = 1
       enddo
    enddo
    close(20)

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)

    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
        if(m<md.or.m>mf.or.l<ld.or.l>lf)then
             cell(l,m) = cell2
        else
          l1 = l - ld + b%l_charge(1,numproc)
          m1 = m - md + b%m_charge(1,numproc)
          if (tab_in(l1,m1)==1) then
             cell(l,m) = cell1
          else
             cell(l,m) = cell2
          endif
        end if
       enddo
    enddo

  end subroutine initialisation_rabbit


  subroutine initialisation_Kothe_Rider(b,tab_cl_ini)
    use m_maillage
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(2), target, intent(in) :: tab_cl_ini

    integer :: l, m, ld, lf, md, mf
    real*8 :: rayon, centre(2)
    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), pointer :: cell1, cell2
    type (STR_cell), dimension(:,:), pointer :: cell

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)

    open(unit=20,file=trim(path_input)//"parametres_initialisation.dat")
    read(20,*) centre
    read(20,*) rayon

    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          if (distance(centre_cl(:,l,m), centre(:))<=rayon) then
             cell(l,m) = cell1
          else
             cell(l,m) = cell2
          endif
       end do
    end do
  end subroutine initialisation_Kothe_Rider

  subroutine initialisation_creneau(b,tab_cl_ini)
    use m_maillage
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(2), target, intent(in) :: tab_cl_ini

    integer :: l, m, ld, lf, md, mf
    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), pointer :: cell1, cell2
    type (STR_cell), dimension(:,:), pointer :: cell

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)


    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          if (centre_cl(1,l,m)<=0.6 .and. centre_cl(1,l,m)>=0.4 .and. centre_cl(2,l,m)>= 0.05 .and. centre_cl(2,l,m)<=0.15) then
             cell(l,m) = cell1
          else
             cell(l,m) = cell2
          endif
       enddo
    enddo
  end subroutine initialisation_creneau

  subroutine initialisation_circle(b, tab_cl_ini)
    use m_maillage
    implicit none

    type (STR_BLOC), intent(inout) :: b
    type (STR_cell), dimension(2), target, intent(in) :: tab_cl_ini

    integer :: l, m, ld, lf, md, mf
    real*8 :: rayon, centre(2)
    real*8, dimension(:,:,:), pointer :: centre_cl
    type (STR_cell), pointer :: cell1, cell2
    type (STR_cell), dimension(:,:), pointer :: cell

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    centre_cl => b%centre_cl
    cell => b%cell
    cell1 => tab_cl_ini(1)
    cell2 => tab_cl_ini(2)

    open(unit=20,file=trim(path_input)//"parametres_initialisation.dat")
    read(20,*) centre
    read(20,*) rayon
    close(20)

    do m = md-nMF, mf+nMF
       do l = ld-nMF, lf+nMF
          if (distance(centre_cl(:,l,m), centre(:))<=rayon) then
             cell(l,m) = cell1
          else
             cell(l,m) = cell2
          endif
       enddo
    enddo

  end subroutine initialisation_circle

end module m_initialisation_fluide
