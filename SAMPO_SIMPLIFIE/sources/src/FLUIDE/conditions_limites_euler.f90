module CL_euler

  use m_struct
  use decodage
  use synchronisation_interbloc
  use outils_conditions_limites
  use parametres_fluide
  use outils_fluide
  use parametres_globaux, only : nMF

  implicit none

  private cond_limite_injection_temps

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     Conditions limites pour la partie Euler      !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine condition_limites_Euler(sb)
    type (STR_SUPER_BLOC), pointer :: sb
    type (STR_BLOC), pointer :: b
    type (STR_PATCH), pointer :: patch
    integer :: i, ib, iface, ip

    if (sb%nb_membre /=1 .or. sb%periodique) then
       call synchro_interbloc(sb, sync_U)
    end if

    do i = 1, sb%nb_membre
       ib = sb%membre(i)
       b=>acces_bloc(ib)%pointeur

       if ( bloc_est_local (ib)) then
          do iface = 1, 4
             do ip = 1, b%face(iface)%nb_patch
                patch => b%face(iface)%patch(ip)
                call CL_patch_Euler(b, sb%dt, patch)
             end do
          end do
       end if
    end do

    !! pour la gestion des coins (EUCCLHYD, UMUSCL, ...)
    do i = 1, sb%nb_membre
       ib = sb%membre(i)
       b=>acces_bloc(ib)%pointeur
       if ( bloc_est_local (ib).and. b%st_coin>0) call CL_coins_Euler(b, sb)
    end do
  end subroutine condition_limites_Euler


  subroutine CL_patch_Euler(b, dt, patch)
    use m_temps

    type (STR_BLOC), pointer :: b
    real*8, intent(in) ::dt
    type (STR_PATCH), pointer :: patch

    integer :: l, m, ld, lf, md, mf, ll, mm
    integer :: l_deb, l_fin, m_deb, m_fin
    integer :: l_face, m_face
    integer :: n_ref, n_face, h_l, h_m
    integer :: ideb, ifin
    integer :: l_int, m_int
    integer :: xd !! Dimension de la direction ou l (respectivement m) est constante pour assurer la condition periodique
    integer :: ierreur
    real*8 :: matrice_sym(1:2,1:2), V_int(1:2), V_ext(1:2)
    real*8, dimension(:,:), allocatable :: normale
    real*8 :: p, rho_eps, eps ! pour l injection subsonique
    real*8 :: u_inj, v_inj
    type (STR_BoundaryCondition), pointer :: bc
    type (STR_cell), pointer :: infi
    real*8 :: fac, u_paroi, v_paroi, rho, rho_int

    ld = b%ld; lf=b%lf; md=b%md; mf=b%mf
    ideb = patch%ideb; ifin = patch%ifin
    bc => patch%bc

    select case (bc%position)
    case(4) !North
       l_deb = ideb; l_fin = ifin
       m_deb = mf+1; m_fin = mf+nMF
       n_ref = mf
       n_face = mf
       h_l = 1
       h_m = 0

       allocate( normale(2,l_deb:l_fin) )
       normale(1:2,l_deb:l_fin) = b%n_dir_m(1:2,l_deb:l_fin,mf)

    case(3) ! South
       l_deb = ideb; l_fin = ifin
       m_deb = md-nMF; m_fin = md-1
       n_ref = md
       n_face = md-1
       h_l = 1
       h_m = 0

       allocate( normale(2,l_deb:l_fin) )
       normale(1:2,l_deb:l_fin) = b%n_dir_m(1:2,l_deb:l_fin,md-1)

    case(2) ! East
       l_deb = lf+1; l_fin = lf+nMF
       m_deb = ideb; m_fin = ifin
       n_ref = lf
       n_face = lf
       h_l = 0
       h_m = 1

       allocate( normale(2,m_deb:m_fin) )
       normale(1:2,m_deb:m_fin) = b%n_dir_l(1:2,lf,m_deb:m_fin)

    case(1) ! West
       l_deb = ld-nMF; l_fin = ld-1
       m_deb = ideb; m_fin = ifin
       n_ref = ld
       n_face = ld-1
       h_l = 0
       h_m = 1

       allocate( normale(2,m_deb:m_fin) )
       normale(1:2,m_deb:m_fin) = b%n_dir_l(1:2,ld-1,m_deb:m_fin)

    case default
       print*, "Position de la condition aux limites inconnue ", bc%position
       call arret_code
    end select

    xd = (mf-md+1) * h_l + (lf-ld+1) * h_m

    select case(trim(bc%condition))
    case("interbloc", "interproc", "periodique")
       ! deja fait dans synchro_interbloc

    case("flux_nul")  ! condition aux limites transparentes
       do m = m_deb, m_fin
          do l = l_deb, l_fin
             b%U(:,l,m) = b%U( :, n_ref*(1-h_l)+l*h_l, n_ref*(1-h_m)+m*h_m )
          end do
       end do
       
    case("paroi_global", "symetrie", "paroi_global_adiabatique", "paroi_flux_impose")
       do m = m_deb, m_fin
          do l = l_deb, l_fin
             l_face = 2*n_face + 1 - l
             m_face = 2*n_face + 1 - m
             b%U(:,l,m) = b%U( :,l_face*(1-h_l)+l*h_l, m_face*(1-h_m)+m*h_m )

             matrice_sym(1,1) = 1.d0-2.d0*normale(1,l*h_l+m*h_m)**2
             matrice_sym(2,1) = -2.d0*normale(1,l*h_l+m*h_m)*normale(2,l*h_l+m*h_m)
             matrice_sym(1,2) = matrice_sym(2,1)
             matrice_sym(2,2) = 1.d0-2.d0*normale(2,l*h_l+m*h_m)**2

             ! on travaille sur rho*u et rho*v mais c'est bon car on a rho_ext=rho_int
             V_int = (/ b%U(irhou, l_face*(1-h_l)+l*h_l, m_face*(1-h_m)+m*h_m), &
                  & b%U(irhov, l_face*(1-h_l)+l*h_l, m_face*(1-h_m)+m*h_m) /)
             V_ext = matmul( matrice_sym, V_int )
             b%U(irhou,l,m) = V_ext(1)
             b%U(irhov,l,m) = V_ext(2)
          end do
       end do

    case("injection_carbone")  ! comme symmetrie, ... mais on rempli bc%Ci
       do m = m_deb, m_fin
          do l = l_deb, l_fin
             l_face = 2*n_face + 1 - l
             m_face = 2*n_face + 1 - m
             b%U(:,l,m) = b%U( :,l_face*(1-h_l)+l*h_l, m_face*(1-h_m)+m*h_m )

             matrice_sym(1,1) = 1.d0-2.d0*normale(1,l*h_l+m*h_m)**2
             matrice_sym(2,1) = -2.d0*normale(1,l*h_l+m*h_m)*normale(2,l*h_l+m*h_m)
             matrice_sym(1,2) = matrice_sym(2,1)
             matrice_sym(2,2) = 1.d0-2.d0*normale(2,l*h_l+m*h_m)**2

             ! on travaille sur rho*u et rho*v mais c'est bon car on a rho_ext=rho_int
             V_int = (/ b%U(irhou, l_face*(1-h_l)+l*h_l, m_face*(1-h_m)+m*h_m), &
                  & b%U(irhov, l_face*(1-h_l)+l*h_l, m_face*(1-h_m)+m*h_m) /)
             V_ext = matmul( matrice_sym, V_int )
             b%U(irhou,l,m) = V_ext(1)
             b%U(irhov,l,m) = V_ext(2)

             ! ToDo : utile ?
             !bc%Ci_paroi(:,(l-l_deb+1)*h_l+(m-m_deb+1)*h_m) = b%Ci(:,l_face*(1-h_l)+l*h_l, m_face*(1-h_m)+m*h_m )
          end do
       end do

    case("entree_supersonique")
       do m = m_deb, m_fin
          do l = l_deb, l_fin
             b%U(:,l,m) = bc%U_infi(:)
          end do
       end do

    case("injection")
       do m = m_deb, m_fin
          do l = l_deb, l_fin
             ! indice qui commence a 1 (pour les tableaux de l'itf)
             ll = l - l_deb + 1
             mm = m - m_deb + 1
             ! indice du point interieur
             l_int = n_ref*(1-h_l)+l*h_l
             m_int = n_ref*(1-h_m)+m*h_m

             ! fraction volumique
             b%U(iz,l,m) = bc%z_injection(ll*h_l+mm*h_m)
             if (b%U(iz,l,m)==1.d0) then
                b%U(irho1z1,l,m) = bc%rho_injection(ll*h_l+mm*h_m)
                b%U(irho2z2,l,m) = 0.d0
             else if (b%U(iz,l,m)==0.d0) then
                b%U(irho1z1,l,m) = 0.d0
                b%U(irho2z2,l,m) = bc%rho_injection(ll*h_l+mm*h_m)
             end if
!!! densite de melange
             call rho_from_U(b%U(:,l,m), rho)
             call rho_from_U(b%U(:,l_int,m_int), rho_int)

!!! on veut imposer la vitesse d'injection a la paroi
             u_paroi = bc%u_injection(l*h_l+m*h_m)*normale(1,l*h_l+m*h_m)
             v_paroi = bc%u_injection(l*h_l+m*h_m)*normale(2,l*h_l+m*h_m)
             ! u et v
             b%U(irhou,l,m) = 2.d0*u_paroi - b%U(irhou,l_int,m_int)/rho_int
             b%U(irhov,l,m) = 2.d0*v_paroi - b%U(irhov,l_int,m_int)/rho_int
             ! rho u et rho v
             b%U(irhou,l,m) = rho * b%U(irhou,l,m) ! rho u
             b%U(irhov,l,m) = rho * b%U(irhov,l,m) ! rho v

             p = b%p(l,m)
             fac = b%c(l,m)*dt/b%aire(l,m)
             u_inj = bc%u_injection(ll*h_l+mm*h_m)*normale(1,l*h_l+m*h_m) ! bizarre (diff avec u_paroi ??)
             v_inj = bc%u_injection(ll*h_l+mm*h_m)*normale(2,l*h_l+m*h_m)

             call pression_entree_subsonique_2d(p, &
                  (/b%u_x(l_int,m_int), b%u_y(l_int,m_int)/), b%p(l_int,m_int), &
                  (/u_inj, v_inj/), b%p_n(l,m), normale(:,l*h_l+m*h_m), &
                  b%rho(l,m)*b%c(l,m), fac)

             if (b%U(iz,l,m)==1.d0) then
                eps = EOS_eps_from_p_rho(fluide1%EOS, p, rho)
             else if (b%U(iz,l,m)==0.d0) then
                eps = EOS_eps_from_p_rho(fluide2%EOS, p, rho) !b%U(irho,l,m)
             else
                print*, "pb injection, z ne vaut ni 0 ni 1"
                call arret_code
             end if
             rho_eps = rho*eps

             b%U(irhoE,l,m) = rho_eps + &
                  .5d0*( b%U(irhou,l,m)**2+b%U(irhov,l,m)**2 )/rho
          end do
       end do


    case("entree_subsonique", "entree_subsonique_t") ! on impose u, v, rho, z de l'exterieur
       do m = m_deb, m_fin
          do l = l_deb, l_fin

             infi => bc%cell_infi

             !! cas de la compression de bulle
             if (trim(bc%condition) == "entree_subsonique_t") then
                infi%u = cond_limite_injection_temps(b%centre_cl(2,l,m), temps_n)
                call rho_from_U(bc%U_infi, rho)
                bc%U_infi(irhou) = rho*infi%u
             end if

             b%U(irho1z1:irhov,l,m) = bc%U_infi(irho1z1:irhov)

             l_int = n_ref*(1-h_l)+l*h_l
             m_int = n_ref*(1-h_m)+m*h_m

             p = b%p(l,m)
             fac = b%c(l,m)*dt/b%aire(l,m)
             call rho_from_U(b%U(:,l,m), rho)
             call pression_entree_subsonique_2d(p, &
                  (/b%u_x(l_int,m_int), b%u_y(l_int,m_int)/), b%p(l_int,m_int), &
                  (/infi%u, infi%v/), b%p_n(l,m), normale(:,l*h_l+m*h_m), &
                  b%rho(l,m)*b%c(l,m), fac)
             eps = EOS_eps_from_p_rho(infi%fluide%EOS, p, rho)
             rho_eps = rho*eps

             b%U(irhoE,l,m) = rho_eps + &
                  .5d0*( b%U(irhou,l,m)**2+b%U(irhov,l,m)**2 )/rho

          end do
       end do

    case("sortie_subsonique")
       ! rho, rho y, z sont ceux de l'interieur du domaine
       ! la pression est donnee par l'exterieur du domaine
       ! les vitesses u_x et u_y sont calculees grace a la caracteristique
       ! (rho e est mis a jour en consequence)
       do m = m_deb, m_fin
          do l = l_deb, l_fin
             l_int = n_ref*(1-h_l)+l*h_l
             m_int = n_ref*(1-h_m)+m*h_m

             b%U(irho1z1:irho2z2,l,m) = b%U(irho1z1:irho2z2,l_int,m_int)
             b%U(irhov,l,m) = b%U(irhov,l_int,m_int)
             b%U(iz,l,m)    = b%U(iz,l_int, m_int)

             infi => bc%cell_infi
             fac = b%c(l,m)*dt/b%aire(l,m)

             call rho_from_U(b%U(:,l,m), rho)
             call rho_from_U(b%U_n(:,l,m), rho_int) !rho_n
             call vitesse_sortie_subsonique(b%U(irhou:irhov,l,m), &
                  b%U(irhou:irhov,l_int,m_int)/rho, &
                  b%U_n(irhou:irhov,l,m)/rho_int, &
                  b%p(l_int, m_int), infi%p, b%rho(l,m)*b%c(l,m), fac)

             b%U(irhou:irhov,l,m) = rho*b%U(irhou:irhov,l,m)

             eps = EOS_eps_from_p_rho(infi%fluide%EOS, infi%p, rho)
             rho_eps = rho * eps
             b%U(irhoE,l,m) = rho_eps + &
                  .5d0*( b%U(irhou,l,m)**2+b%U(irhov,l,m)**2)/rho
          end do
       end do

    case default
       print*, "Type de condition aux limites inconnu ", bc%condition
       call arret_code
    end select


!!! on decode U pour avoir les variables  rho, u_x
    ierreur = 0
    call decode_U(l_deb, l_fin, m_deb, m_fin, b, ierreur)
    if (ierreur /= 0) then
       print*, "erreur: pb dans le decodage : CL_patch_Euler"
!       call arret_code
    end if

!!$    select case(trim(bc%condition))
!!$    case("paroi_global", "symetrie", "paroi_global_adiabatique", "paroi_flux_impose")
!!$       do m = m_deb, m_fin
!!$          do l = l_deb, l_fin
!!$             l_face = 2*n_face + 1 - l
!!$             m_face = 2*n_face + 1 - m
!!$             b%p(l,m) = b%p(l_face*(1-h_l)+l*h_l, m_face*(1-h_m)+m*h_m )
!!$             b%rho_eps(l,m) = b%rho_eps(l_face*(1-h_l)+l*h_l, m_face*(1-h_m)+m*h_m )
!!$          end do
!!$       end do
!!$    end select

    deallocate (normale)

  end subroutine CL_patch_Euler

!!! utile seulement avec le schema a 9 points
subroutine CL_coins_Euler(b, sb)
  type (STR_BLOC), pointer :: b
  type (STR_SUPER_BLOC), pointer :: sb

  integer :: ld, lf, md, mf

  ld = b%ld; lf=b%lf; md=b%md; mf=b%mf

  !!! coins lies au decoupage MPI
  if (sb%nb_membre /=1 .or. sb%periodique) call synchro_interbloc(sb, sync_U)

  !!! coin entre la face 1 et 3
  call coin_Euler(b, b%coin(1),ld  ,md  ,ld-1,md-1)
  call coin_Euler(b, b%coin(1),ld+1,md  ,ld-2,md-1)
  call coin_Euler(b, b%coin(1),ld  ,md+1,ld-1,md-2)
  call coin_Euler(b, b%coin(1),ld+1,md+1,ld-2,md-2)

  if(nMF>2)then
    call coin_Euler(b, b%coin(1),ld+2,md  ,ld-3,md-1)
    call coin_Euler(b, b%coin(1),ld+2,md+1,ld-3,md-2)
    call coin_Euler(b, b%coin(1),ld+2,md+2,ld-3,md-3)
    call coin_Euler(b, b%coin(1),ld+1,md+2,ld-2,md-3)
    call coin_Euler(b, b%coin(1),ld  ,md+2,ld-1,md-3)
  end if

  !!! coin entre la face 3 et 2
  call coin_Euler(b, b%coin(2),lf  ,md  ,lf+1,md-1)
  call coin_Euler(b, b%coin(2),lf-1,md  ,lf+2,md-1)
  call coin_Euler(b, b%coin(2),lf  ,md+1,lf+1,md-2)
  call coin_Euler(b, b%coin(2),lf-1,md+1,lf+2,md-2)

  if(nMF>2)then
    call coin_Euler(b, b%coin(2),lf-2,md  ,lf+3,md-1)
    call coin_Euler(b, b%coin(2),lf-2,md+1,lf+3,md-2)
    call coin_Euler(b, b%coin(2),lf-2,md+2,lf+3,md-3)
    call coin_Euler(b, b%coin(2),lf-1,md+2,lf+2,md-3)
    call coin_Euler(b, b%coin(2),lf  ,md+2,lf+1,md-3)
  end if

  !!! coin entre la face 2 et 4
  call coin_Euler(b, b%coin(3),lf  ,mf  ,lf+1,mf+1)
  call coin_Euler(b, b%coin(3),lf-1,mf  ,lf+2,mf+1)
  call coin_Euler(b, b%coin(3),lf  ,mf-1,lf+1,mf+2)
  call coin_Euler(b, b%coin(3),lf-1,mf-1,lf+2,mf+2)

  if(nMF>2)then
    call coin_Euler(b, b%coin(3),lf-2,mf  ,lf+3,mf+1)
    call coin_Euler(b, b%coin(3),lf-2,mf-1,lf+3,mf+2)
    call coin_Euler(b, b%coin(3),lf-2,mf-2,lf+3,mf+3)
    call coin_Euler(b, b%coin(3),lf-1,mf-2,lf+2,mf+3)
    call coin_Euler(b, b%coin(3),lf  ,mf-2,lf+1,mf+3)
  end if

  !!! coin entre la face 4 et 1
  call coin_Euler(b, b%coin(4),ld  ,mf  ,ld-1,mf+1)
  call coin_Euler(b, b%coin(4),ld+1,mf  ,ld-2,mf+1)
  call coin_Euler(b, b%coin(4),ld  ,mf-1,ld-1,mf+2)
  call coin_Euler(b, b%coin(4),ld+1,mf-1,ld-2,mf+2)

  if(nMF>2)then
    call coin_Euler(b, b%coin(4),ld+2,mf  ,ld-3,mf+1)
    call coin_Euler(b, b%coin(4),ld+2,mf-1,ld-3,mf+2)
    call coin_Euler(b, b%coin(4),ld+2,mf-2,ld-3,mf+3)
    call coin_Euler(b, b%coin(4),ld+1,mf-2,ld-2,mf+3)
    call coin_Euler(b, b%coin(4),ld  ,mf-2,ld-1,mf+3)
  end if

end subroutine CL_coins_Euler


subroutine coin_Euler(b, coin, l_int,m_int, l_ext,m_ext)
  use decodage

  type (STR_BLOC), pointer     :: b
  type (STR_COIN), intent(in)  :: coin
  integer        , intent(in)  :: l_int,m_int, l_ext,m_ext

  integer :: ierreur, i
  integer :: lp, mp
  integer :: ld, lf, md, mf
  real*8, dimension(b%neq,b%neq)    :: mat
  real*8, dimension(:,:,:), pointer :: U
  U => b%U

  ld = b%ld; lf=b%lf; md=b%md; mf=b%mf

  select case (coin%condition)
  case(CL_INTERPROCvs, CL_INTERBLOCvsINTERBLOC,&
    CL_INTERBLOCvsFLUX_NUL   , CL_FLUX_NULvsINTERBLOC, &
    CL_INTERBLOCvsSYMETRIE   , CL_SYMETRIEvsINTERBLOC, &
    CL_INTERBLOCvsPAROI      , CL_PAROIvsINTERBLOC, &
    CL_INTERBLOCvsPAROI_ADIAB, CL_PAROI_ADIABvsINTERBLOC, &
    CL_INTERBLOCvsINJECTION  , CL_INJECTIONvsINTERBLOC, &
    CL_INTERBLOCvsENTREE_SUPER, CL_ENTREE_SUPERvsINTERBLOC, &
    CL_INTERBLOCvsPERIODIQUE , CL_PERIODIQUEvsINTERBLOC)
    ! U deja synchronise entre les blocs
  case(CL_FLUX_NULvsFLUX_NUL)
    U(:,l_ext,m_ext) = U(:,l_int,m_int)
  case(CL_PAROIvsPAROI, CL_PAROIvsPAROI_ADIAB, &
    CL_PAROI_ADIABvsPAROI, CL_PAROI_ADIABvsPAROI_ADIAB, &
    CL_PAROIvsSYMETRIE, CL_SYMETRIEvsPAROI, &
    CL_PAROI_ADIABvsSYMETRIE, CL_SYMETRIEvsPAROI_ADIAB, &
    CL_SYMETRIEvsSYMETRIE)
    U(:    ,l_ext,m_ext) = U(:     ,l_int,m_int)
    U(irhou,l_ext,m_ext) =-U(irhou,l_int,m_int)
    U(irhov,l_ext,m_ext) =-U(irhov,l_int,m_int)

  case (CL_PAROIvsFLUX_NUL, CL_PAROI_ADIABvsFLUX_NUL, CL_SYMETRIEvsFLUX_NUL)
    mat = 0.d0
    mat(irhou:irhov,irhou:irhov) = - 2.d0*produit_tensoriel(2,coin%n1(:),coin%n1(:))
    do i = 1, b%neq
      mat(i,i) = mat(i,i) + 1.d0
    end do
    U(:,l_ext,m_ext) = matmul(mat,U(:,l_int,m_int))

  case(CL_FLUX_NULvsPAROI, CL_FLUX_NULvsPAROI_ADIAB, CL_FLUX_NULvsSYMETRIE)
    mat = 0.d0
    mat(irhou:irhov,irhou:irhov) = - 2.d0*produit_tensoriel(2,coin%n2(:),coin%n2(:))
    do i = 1, b%neq
      mat(i,i) = mat(i,i) + 1.d0
    end do
    U(:,l_ext,m_ext) = matmul(mat,U(:,l_int,m_int))

  case(CL_ENTREE_SUPERvsFLUX_NUL, CL_ENTREE_SUPERvsSYMETRIE)
    U(:,l_ext,m_ext) = coin%patch1%bc%U_infi(:)

  case(CL_FLUX_NULvsENTREE_SUPER, CL_SYMETRIEvsENTREE_SUPER)
    U(:,l_ext,m_ext) = coin%patch2%bc%U_infi(:)

  case (CL_PAROI_ADIABvsENTREE_SUPER,CL_PAROIvsENTREE_SUPER)
    mat = 0.d0
    mat(irhou:irhov,irhou:irhov) = - 2.d0*produit_tensoriel(2,coin%n1(:),coin%n1(:))
    do i = 1, b%neq
      mat(i,i) = mat(i,i) + 1.d0
    end do
    U(:,l_ext,m_ext) = matmul(mat,coin%patch2%bc%U_infi(:))

  case(CL_ENTREE_SUPERvsPAROI_ADIAB,CL_ENTREE_SUPERvsPAROI)
    mat = 0.d0
    mat(irhou:irhov,irhou:irhov) = - 2.d0*produit_tensoriel(2,coin%n2(:),coin%n2(:))
    do i = 1, b%neq
      mat(i,i) = mat(i,i) + 1.d0
    end do
    U(:,l_ext,m_ext) = matmul(mat,coin%patch1%bc%U_infi(:))

  case(CL_PERIODIQUEvsPERIODIQUE)
    lp = mod(l_ext+lf-ld,lf-ld+1)+1
    mp = mod(m_ext+mf-md,mf-md+1)+1
    U(:,l_ext,m_ext) = U(:,lp,mp)

  case (CL_PAROIvsPERIODIQUE,CL_PAROI_ADIABvsPERIODIQUE)
    lp = mod(l_ext+lf-ld,lf-ld+1)+1
    mp = m_ext
    U(:,l_ext,m_ext) = U(:,lp,mp)

  case(CL_PERIODIQUEvsPAROI,CL_PERIODIQUEvsPAROI_ADIAB)
    lp = l_ext
    mp = mod(m_ext+mf-md,mf-md+1)+1
    U(:,l_ext,m_ext) = U(:,lp,mp)

  case default
    print*, "cas non gere pour les coins en Euler"
    print*, coin%patch1%bc%condition, coin%patch2%bc%condition
    call arret_code
  end select

  ! decodage
  call decode_U(l_ext, l_ext, m_ext, m_ext, b, ierreur)

end subroutine coin_Euler

end module CL_euler
