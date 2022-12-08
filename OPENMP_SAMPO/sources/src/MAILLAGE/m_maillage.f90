module m_maillage
  use m_raffinement_maillage
  use metrique_maillage
  use generation_maillage
  use lecture_maillage_cgns
  use synchronisation_interbloc
  use parametres_fluide
  implicit none

contains

  subroutine creation_maillage(b, type_maillage)

    type (STR_BLOC), pointer ::  b
    character(len=90), intent(in) :: type_maillage
    integer :: i

    if ( index(type_maillage,"restart") /= 0 ) then
       call maillage_restart(b, type_maillage)
       ! call maillage_restart_plus_decoupe(b, type_maillage)
    else if ( index(type_maillage,".cgns") /= 0 ) then
       call maillage_cgns(b, type_maillage)
    else
       select case(type_maillage)
       case('cartesien')
          call maillage_cartesien(b)
       case('quad')
          call maillage_a_partir_de_4pts(b)
       case('ellipse')
          call maillage_ellipse(b)
       case("sphere")
          call maillage_sphere(b)
       case("ellipse_cartesien")
          call maillage_ellipse_cartesien(b)
       case("bosse_dellacherie")
          call maillage_avec_bosse_dellacherie(b)
       case("plaque_plane_avec_creux")
          call maillage_plaque_plane_avec_creux(b)
       case("random")
          call maillage_random(b)
       case("kershaw")
          call maillage_kershaw(b)
       case default
          print*,""
          print*," Type de maillage non defini", type_maillage
          print*," ERREUR :Creation_maillage "
          print*,""
          call arret_code
       end select

       !! On resserre le maillage que l'on a genere
       if ( b%maillage%resserrer_maillage ) call resserrement_de_maillage(b, b%coord_nd)
    end if

    do i=1, size(acces_super_bloc,1)
       if (acces_super_bloc(i)%pointeur%nature == b%nature) then
          call synchro_interbloc(acces_super_bloc(i)%pointeur, sync_coord_nd)
       end if
    end do

    ! Calcul de la métrique
    call calcul_centre_cellule(b, b%coord_nd, b%centre_cl)
    call calcul_aire(b, b%coord_nd, b%aire)
    call calcul_normales(b, b%coord_nd)

    if (b%nature == FLUIDE) then
       select case(b%calcul_gradient)
       case(WLSQ5,WLSQ9)
          select case(b%calcul_gradient_methode)
          case(Psinv)
             ! A * gradz = dz -> pseudo inverse : gradz = (tAA)-1 * tA * dz
             call calc_Mpinv_direct(b,b%Mpinv,b%st_grad,b%ld-(nMF-1), b%lf+(nMF-1),b%md-(nMF-1), b%mf+(nMF-1))
          case(QR)
             ! on passe par la decomposition A = QR : x = R^-1 Q^T b
             call calc_Mpinv_from_QR(b,b%Mpinv,b%st_grad,b%ld-(nMF-1), b%lf+(nMF-1),b%md-(nMF-1), b%mf+(nMF-1))
          case default
             print*, ' methode non reconnue pour les moindres carrés'
             call arret_code
          end select
       case(GG)
          ! rien à faire
       end select
       
       if(b%visqueux) then
          call calc_metrique_gradient(b)
       end if
    end if

    call coordonees_maillage_01x01(b, b%coord_nd, b%coord_nd_01x01)

    !call symetrisation(b)
  end subroutine creation_maillage

  subroutine maillage_restart(b, type_ini)
    use m_temps
    use variables_globales

    type (STR_BLOC), pointer :: b
    character(len=90), intent(in) :: type_ini
    integer :: iunit=47
    integer :: ld, lf, md, mf
    character(len=90) :: fichier
    character(len=5) :: iter_c_read
    integer :: iter_read, it_restart_read
    integer :: ld_read, lf_read, md_read, mf_read
    real*8, dimension(:,:,:), allocatable :: coord_nd_global

    if (b%nature==FLUIDE) then
       fichier = trim(path_output)//'restart/'//trim(type_ini)//'_maillage_fld.dat'
    end if

    open(unit=iunit, file=fichier, status='old', form='unformatted' )
    read(iunit) iter_read, temps, iter_c_read, it_restart_read
    read(iunit) ld_read, lf_read, md_read, mf_read

    if ( ld_read /= b%ld_tot .or. lf_read /= b%lf_tot &
         .or. md_read /= b%md_tot .or. mf_read /= b%mf_tot ) then
       print*, ""
       print*, "IL Y A UNE INCOMPATIBILITE ENTRE LES INDICES DU RESTART ET CEUX DES DONNEES"
       print*, ld_read, b%ld_tot, lf_read, b%lf_tot
       print*, md_read, b%md_tot, mf_read, b%mf_tot
       print*, ""
       call arret_code
    end if

    allocate(coord_nd_global(1:2, b%ld_tot-(nMF+1):b%lf_tot+nMF, b%md_tot-(nMF+1):b%mf_tot+nMF))

    read(iunit) b%xin, b%xout, b%yin, b%yout
    read(iunit)
    read(iunit)
    read(iunit) coord_nd_global(1:2, b%ld_tot-(nMF+1):b%lf_tot+nMF, b%md_tot-(nMF+1):b%mf_tot+nMF)
    close(iunit)

    ld = b%l_charge(1, numproc)
    lf = b%l_charge(2, numproc)
    md = b%m_charge(1, numproc)
    mf = b%m_charge(2, numproc)

    b%coord_nd(:, b%ld-(nMF+1):b%lf+nMF, b%md-(nMF+1):b%mf+nMF) = coord_nd_global(:, ld-(nMF+1):lf+nMF, md-(nMF+1):mf+nMF)

    deallocate(coord_nd_global)

  end subroutine maillage_restart

  subroutine maillage_restart_plus_decoupe(b, type_ini)
    use m_temps
    use variables_globales

    type (STR_BLOC), pointer :: b
    character(len=90), intent(in) :: type_ini
    integer :: iunit=47
    integer :: ld, lf, md, mf
    character(len=90) :: fichier
    character(len=5) :: iter_c_read
    integer :: iter_read, it_restart_read
    integer :: ld_read, lf_read, md_read, mf_read
    real*8, dimension(:,:,:), allocatable :: coord_nd_global

    if (b%nature==FLUIDE) then
       fichier = trim(path_output)//'restart/'//trim(type_ini)//'_maillage_fld.dat'
    end if

    open(unit=iunit, file=fichier, status='old', form='unformatted' )
    read(iunit) iter_read, temps, iter_c_read, it_restart_read
    read(iunit) ld_read, lf_read, md_read, mf_read

    if ( ld_read /= b%ld_tot .or. lf_read /= b%lf_tot &
         .or. md_read /= b%md_tot .or. mf_read /= b%mf_tot ) then
       print*, ""
       print*, "IL Y A UNE INCOMPATIBILITE ENTRE LES INDICES DU RESTART ET CEUX DES DONNEES"
       print*, ld_read, b%ld_tot, lf_read, b%lf_tot
       print*, md_read, b%md_tot, mf_read, b%mf_tot
       print*, ""
       !call arret_code
    end if

    allocate(coord_nd_global(1:2, ld_read-(nMF+1):lf_read+nMF, md_read-(nMF+1):mf_read+nMF))

    read(iunit) b%xin, b%xout, b%yin, b%yout
    read(iunit)
    read(iunit)
    read(iunit) coord_nd_global(1:2, ld_read-(nMF+1):lf_read+nMF, md_read-(nMF+1):mf_read+nMF)
    close(iunit)

    ld = b%l_charge(1, numproc)
    lf = b%l_charge(2, numproc)
    md = b%m_charge(1, numproc)
    mf = b%m_charge(2, numproc)

    b%coord_nd(:,b%ld-(nMF+1):b%lf+nMF, b%md-(nMF+1):b%mf+nMF) = coord_nd_global(:,ld-(nMF+1):lf+nMF, md-(nMF+1):mf+nMF)

    deallocate(coord_nd_global)

  end subroutine maillage_restart_plus_decoupe

  subroutine symetrisation(b)

    type (STR_BLOC), pointer :: b
    integer :: l, m, lsym


    ! Inclusion des mailles fictives via nMF ?
    do m = b%md-3, b%mf+2
       do l = b%ld-3, b%lf/2-1
          lsym = b%lf-l

          b%coord_nd(1,l,m) =  b%coord_nd(1,lsym,m)
          b%coord_nd(2,l,m) = -b%coord_nd(2,lsym,m)
       end do
       b%coord_nd(2,b%lf/2,m) = 0.d0
    end do
    do m = b%md-2, b%mf+2
       do l = b%ld-2, b%lf/2
          lsym = b%lf-l+1

          b%centre_cl(1,l,m) =  b%centre_cl(1,lsym,m)
          b%centre_cl(2,l,m) = -b%centre_cl(2,lsym,m)

          b%aire(l,m) = b%aire(lsym,m)
       end do
    end do

    do m = b%md-1, b%mf
       do l = b%ld-1, b%lf/2-1
          lsym = b%lf-l

          b%dl_l(l,m) = b%dl_l(lsym,m)
          b%n_dir_l(1,l,m) = - b%n_dir_l(1,lsym,m)
          b%n_dir_l(2,l,m) =  b%n_dir_l(2,lsym,m)
       end do
    end do

    do m = b%md-1, b%mf
       do l = b%ld, b%lf/2
          lsym = b%lf-l+1

          b%dl_m(l,m) = b%dl_m(lsym,m)

          b%n_dir_m(1,l,m) =  b%n_dir_m(1,lsym,m)
          b%n_dir_m(2,l,m) = - b%n_dir_m(2,lsym,m)
       end do
    end do

  end subroutine symetrisation

end module m_maillage
