module CL_visqueux

  use m_struct
  use synchronisation_interbloc

  implicit none
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   Conditions limites pour la partie visqueuse    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine condition_limites_visqueux(b, sb, no_synchro)
    
    type (STR_BLOC)      , pointer  :: b
    type (STR_SUPER_BLOC), pointer  :: sb
    logical              , optional :: no_synchro
    
    integer :: ld, lf, md, mf
    integer, dimension(6) :: indices
    real*8, dimension(:,:), pointer :: normale
    type (STR_FACE), pointer  :: face
    
    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf

!!!--------------------------------------------------!!!
!!!             Traitement des 4 faces               !!!
!!!--------------------------------------------------!!!
    ! indices = premier indice en n a remplir    (n = l ou m suivant la face)
    !           dernier indice en n a remplir
    !           indice en n du point a l'interieur du domaine 
    !           indice constant de la face du bord
    !           h_l : 1 si la face est telle que m=constante, 0 sinon
    !           h_m : 1 si la face est telle que l=constante, 0 sinon
    
!!! En visqueux il ne faut remplir qu'une seule rangee de mailles fictives
!!! la deuxieme rangee est pour la reconstruction (ordre 2) qui n'est faite que pour les termes Euler
    
    ! face 1
    face => b%face(1)
    indices = (/ld-1, ld-1, ld, ld-1, 0, 1/) 
    normale => b%n_dir_l(1:2,ld-1,md:mf)
    call CL_face_visqueux(b, face, indices, normale)
   
    ! face 2
    face => b%face(2)
    indices = (/lf+1, lf+1, lf, lf, 0, 1/) 
    normale => b%n_dir_l(1:2,lf,md:mf)
    call CL_face_visqueux(b, face, indices, normale)
    
    ! face 3
    face => b%face(3)
    indices = (/md-1, md-1, md, md-1, 1, 0/) 
    normale => b%n_dir_m(1:2,ld:lf,md-1) 
    call CL_face_visqueux(b, face, indices, normale)

    ! face 4
    face => b%face(4)
    indices = (/mf+1, mf+1, mf, mf, 1, 0/)
    normale => b%n_dir_m(1:2,ld:lf,mf)
    call CL_face_visqueux(b, face, indices, normale)
    
!!!--------------------------------------------------!!!
!!!               Traitement des coins               !!!
!!!--------------------------------------------------!!!
   ! if ((b%lf/=1 .and. b%mf/=1) .and. b%maillage%type/="cartesien") then
    if (b%lf/=1 .and. b%mf/=1) then
    ! indices = premier indice en l a remplir
       !           dernier indice en l a remplir
       !           premier indice en m a remplir
       !           dernier indice en m a remplir
       !           indice en l du point a l'interieur du domaine    
       !           indice en m du point a l'interieur du domaine    
       
!!$       ! coin entre la face 1 et 3
!!$       coin => b%coin(1)
!!$       if (coin%ib_opp == -1) then
!!$          indices = (/ld-2, ld-1, md-2, md-1, ld, md/)
!!$          call CL_coin_visqueux(b, indices, coin%patch1, coin%patch2)
!!$       end if
!!$       
!!$       ! coin entre la face 3 et 2
!!$       coin => b%coin(2)
!!$       if (coin%ib_opp == -1) then
!!$          indices = (/lf+1, lf+2, md-2, md-1, lf, md/)
!!$          call CL_coin_visqueux(b, indices, coin%patch1, coin%patch2)
!!$       end if
!!$       
!!$       ! coin entre la face 2 et 4
!!$       coin => b%coin(3)
!!$       if (coin%ib_opp == -1) then
!!$          indices = (/lf+1, lf+2, mf+1, mf+2, lf, mf/)
!!$          call CL_coin_visqueux(b, indices, coin%patch1, coin%patch2)
!!$       end if
!!$             
!!$       ! coin entre la face 4 et 1
!!$       coin => b%coin(4)
!!$       if (coin%ib_opp == -1) then
!!$          indices = (/ld-2, ld-1, mf+1, mf+2, ld, mf/)
!!$          call CL_coin_visqueux(b, indices, coin%patch1, coin%patch2)
!!$       end if
       
!!! coin entre la face 1 et 3
       call coin_visqueux(b, b%coin(1), ld,md, ld-1,md-1)
!!! coin entre la face 3 et 2
       call coin_visqueux(b, b%coin(2), lf,md, lf+1,md-1)
!!! coin entre la face 2 et 4
       call coin_visqueux(b, b%coin(3), lf,mf, lf+1,mf+1)
!!! coin entre la face 4 et 1
       call coin_visqueux(b, b%coin(4), ld,mf, ld-1,mf+1)
       
    end if

    if (.not. present(no_synchro)) then
       if (sb%nb_membre /=1 .or. sb%periodique) call synchro_interbloc(sb, sync_uvT)
    end if

  end subroutine condition_limites_visqueux

  ! on doit remplir la vitesse (u,v) et la temperature des mailles fictives
  subroutine CL_face_visqueux(b, face, indices, normale)
    use outils_conditions_limites
    use parametres_fluide, only : EOS, CHIMIE_FIGEE

    type (STR_BLOC), pointer :: b
    type (STR_FACE), pointer :: face
    integer, dimension(6), intent(in) :: indices
    real*8, dimension(:,:), intent(in) :: normale
    
    integer :: ip, l, m, k
    integer :: l_deb, l_fin, m_deb, m_fin
    integer :: n_ref, n_face, h_l, h_m
    integer :: lint, mint
    type (STR_PATCH), pointer :: patch
    type (STR_BoundaryCondition), pointer :: bc
    type (STR_cell), pointer :: infi

    real*8 :: uw, Tw, u_paroi, v_paroi, T_paroi
    real*8, dimension(:,:)  , pointer :: u, v, T
    real*8, dimension(:,:)  , pointer :: dl
    real*8, dimension(:,:,:), pointer :: ci


    u => b%u_x; v => b%u_y; T => b%T
    !!! TODO : il faut aussi mettre les CL sur les Ci (chimie figee)
    select case (fluide1%type)
    case(EOS, CHIMIE_FIGEE)
       ci => b%c_esp
    end select

    ! on recupere les indices utiles
    n_ref = indices(3); n_face = indices(4)
    h_l = indices(5); h_m = indices(6)
    
    do ip = 1, face%nb_patch
       patch => face%patch(ip)
       bc => patch%bc
       
       if (h_l==1) then ! face suivant m=constante
          l_deb = patch%ideb; l_fin = patch%ifin
          m_deb = indices(1); m_fin = indices(2)
  
       else ! face suivant l=constante
         l_deb = indices(1); l_fin = indices(2)
         m_deb = patch%ideb; m_fin = patch%ifin
       end if
       
       select case(trim(bc%condition))
       case("interbloc", "interproc")
          ! deja fait dans synchro_interbloc

       case("entree_supersonique", "flux_nul", "symetrie", &
            "entree_subsonique_t", "sortie_subsonique")
          ! deja fait par les conditions aux limites de la partie Euler
          ! faut il le refaire pour l'equation sur l'energie ?

       case("entree_subsonique")
          ! il faut etre coherent avec le choix fait dans outils_conditions_limites
          do m = m_deb, m_fin
             do l = l_deb, l_fin
                lint = n_ref*(1-h_l)+l*h_l
                mint = n_ref*(1-h_m)+m*h_m
                
                infi => bc%cell_infi

                u(l,m) = infi%u
                v(l,m) = infi%v
                
                call temperature_entree_subsonique(T(l,m), b%T(lint, mint), &
                     b%u_x(lint, mint), infi, b%rho(l,m)*b%c(l,m))
             end do
          end do
       
       case("injection_carbone")
          do m = m_deb, m_fin
             do l = l_deb, l_fin
                ! indice du point interieur
                lint = n_ref*(1-h_l)+l*h_l ! l_face*(1-h_l)+l*h_l
                mint = n_ref*(1-h_m)+m*h_m ! m_face*(1-h_m)+m*h_m

                k = (l-l_deb+1)*h_l+(m-m_deb+1)*h_m ! l'indice commence a 1
                
                ! T_ext = bc%T_imposee((l-l_deb+1)*h_l+(m-m_deb+1)*h_m)
!!! on veut imposer la vitesse d'injection et la temperature de fusion a la paroi
                u_paroi = bc%u_injection(k)*normale(1,l*h_l+m*h_m)
                v_paroi = bc%u_injection(k)*normale(2,l*h_l+m*h_m)
                T_paroi = bc%Tw(k)
                
                u(l,m) = 2.d0*u_paroi - u(lint, mint)
                v(l,m) = 2.d0*v_paroi - v(lint, mint)
                T(l,m) = 2.d0*T_paroi - T(lint, mint)
                Ci(:,l,m)   = 2.d0 * bc%Ci_paroi(:,k)    - Ci(:,lint, mint)
                b%rho(l,m)  = 2.d0 * bc%rho_injection(k) - b%rho(lint, mint)
                b%D(l,m)    = 2.d0 * bc%D_paroi(k)       - b%D(lint, mint)
                b%K(l,m)    = 2.d0 * bc%K_paroi(k)       - b%K(lint, mint)
                b%mu(l,m)   = 2.d0 * bc%mu_paroi(k)      - b%mu(lint, mint)
                b%Hi(:,l,m) = 2.d0 * bc%Hi_paroi(:,k)    - b%Hi(:,lint, mint)
                b%y(l,m)    = 1.d0  ! b%rho1(l,m) * b%z(l,m) / b%rho(l,m)
             end do
          end do
          
          
       case("injection") 
          do m = m_deb, m_fin
             do l = l_deb, l_fin
                ! indice du point interieur
                lint = n_ref*(1-h_l)+l*h_l ! l_face*(1-h_l)+l*h_l
                mint = n_ref*(1-h_m)+m*h_m ! m_face*(1-h_m)+m*h_m
                
                ! T_ext = bc%T_imposee((l-l_deb+1)*h_l+(m-m_deb+1)*h_m)
!!! on veut imposer la vitesse d'injection et la temperature de fusion a la paroi
                u_paroi = bc%u_injection(l*h_l+m*h_m)*normale(1,l*h_l+m*h_m)
                v_paroi = bc%u_injection(l*h_l+m*h_m)*normale(2,l*h_l+m*h_m)
                T_paroi = bc%Tw(l*h_l+m*h_m)

                !u(l,m) = 2.d0*u_paroi - u(lint, mint)
                !v(l,m) = 2.d0*v_paroi - v(lint, mint)
                u(l,m) = u_paroi
                v(l,m) = v_paroi
                T(l,m) = 2.d0*T_paroi - T(lint, mint)
             end do
          end do
          
       case("paroi_global")
          do m= m_deb, m_fin
             do l= l_deb, l_fin
                lint = n_ref*(1-h_l)+l*h_l
                mint = n_ref*(1-h_m)+m*h_m

                uw = bc%uw((l-l_deb+1)*h_l+(m-m_deb+1)*h_m)
                Tw = bc%Tw((l-l_deb+1)*h_l+(m-m_deb+1)*h_m) 

                ! vecteur tangent : t = (ny, -nx) | n ^ t = (0,0,-1)
                u(l,m) = 2.d0*uw*normale(2,l*h_l+m*h_m) - u(lint,mint)
                v(l,m) =-2.d0*uw*normale(1,l*h_l+m*h_m) - v(lint,mint)
                T(l,m) = 2.d0*Tw - T(lint,mint)
             end do
          end do
          
       case("paroi_global_adiabatique")
          do m= m_deb, m_fin
             do l= l_deb, l_fin
                lint = n_ref*(1-h_l)+l*h_l
                mint = n_ref*(1-h_m)+m*h_m
                
                uw = bc%uw((l-l_deb+1)*h_l+(m-m_deb+1)*h_m)
                
                ! vecteur tangent : t = (-ny, nx) | n ^ t = (0,0,1)
                u(l,m) = 2.d0*uw*normale(2,l*h_l+m*h_m) - u(lint,mint)
                v(l,m) =-2.d0*uw*normale(1,l*h_l+m*h_m) - v(lint,mint)
                T(l,m) = T(lint,mint)
             end do
          end do
          
       case("paroi_flux_impose") 
          if (h_l==1) then
             dl => b%dl_m
          else
             dl => b%dl_l
          end if
          do m = m_deb, m_fin
             do l = l_deb, l_fin
                lint = n_ref*(1-h_l)+l*h_l
                mint = n_ref*(1-h_m)+m*h_m
                
                uw = bc%uw((l-l_deb+1)*h_l+(m-m_deb+1)*h_m)
                
                ! vecteur tangent : t = (-ny, nx) | n ^ t = (0,0,1)
                u(l,m) = 2.d0*uw*normale(2,l*h_l+m*h_m) - u(lint,mint)
                v(l,m) =-2.d0*uw*normale(1,l*h_l+m*h_m) - v(lint,mint)
!!! Attention : ne fonctionne qu'en cartesien !!!
                T(l,m) = dl(l,m)/b%K(l,m)*bc%flux_impose((l-l_deb+1)*h_l+(m-m_deb+1)*h_m) + T(lint, mint)
             end do
          end do

       case default
          print*, "viscous_boundary_face : Type de condition aux limites inconnu ", bc%condition
          call arret_code
       end select
    end do
  end subroutine CL_face_visqueux
  
  subroutine coin_visqueux(b, coin, l_int,m_int, l_ext,m_ext)
    use parametres_fluide 
    
    type (STR_BLOC), pointer    :: b
    type (STR_COIN), intent(in) :: coin
    integer        , intent(in) :: l_int,m_int, l_ext,m_ext

    real*8                      :: u_int, v_int, T_int
    real*8         , pointer    :: u_ext, v_ext, T_ext
    type (STR_cell), pointer    :: cell_infi1, cell_infi2
    real*8, dimension(:,:,:), pointer :: ci
    
    integer :: k
    real*8 :: u_paroi, v_paroi

    cell_infi1 => coin%patch1%bc%cell_infi
    cell_infi2 => coin%patch2%bc%cell_infi
    
    u_ext => b%u_x(l_ext,m_ext); u_int = b%u_x(l_int,m_int)
    v_ext => b%u_y(l_ext,m_ext); v_int = b%u_y(l_int,m_int)
    T_ext => b%T(l_ext,m_ext)  ; T_int = b%T(l_int,m_int)
     
    select case (fluide1%type)
    case(EOS, CHIMIE_FIGEE)
       ci => b%c_esp
    end select
    
    select case (coin%condition)
    case(CL_INTERPROCvs, CL_INTERBLOCvsINTERBLOC,&
         CL_INTERBLOCvsFLUX_NUL   , CL_FLUX_NULvsINTERBLOC, &
         CL_INTERBLOCvsSYMETRIE   , CL_SYMETRIEvsINTERBLOC, &
         CL_INTERBLOCvsPAROI      , CL_PAROIvsINTERBLOC, &
         CL_INTERBLOCvsPAROI_ADIAB, CL_PAROI_ADIABvsINTERBLOC, &
         CL_INTERBLOCvsINJECTION  , CL_INJECTIONvsINTERBLOC  ) 
       ! (u,v,T) deja synchronise entre les blocs 
       ! il y a peut etre un probleme (a verifier)
    case(CL_PAROI_ADIABvsPAROI_ADIAB, CL_PAROIvsPAROI, &
         CL_PAROI_ADIABvsPAROI, CL_PAROIvsPAROI_ADIAB)
       u_ext =-u_int
       v_ext =-v_int
       T_ext = T_int
    case(CL_FLUX_NULvsFLUX_NUL)
       u_ext = u_int
       v_ext = v_int
       T_ext = T_int

    case(CL_PAROIvsFLUX_NUL, CL_PAROIvsSYMETRIE)
!       u_ext = 2.d0*cell_infi1%u - u_int
!       v_ext = 2.d0*cell_infi1%v - v_int
       u_ext =-u_int
       v_ext =-v_int
       T_ext = 2.d0*coin%patch1%bc%Tw(1) - T_int
       
    case(CL_FLUX_NULvsPAROI, CL_SYMETRIEvsPAROI) 
       u_ext = 2.d0*cell_infi2%u - u_int
       v_ext = 2.d0*cell_infi2%v - v_int
       T_ext = 2.d0*cell_infi2%T - T_int
       
    case(CL_PAROI_ADIABvsFLUX_NUL, CL_PAROI_ADIABvsSYMETRIE)
      ! u_ext = 2.d0*cell_infi1%u - u_int
      ! v_ext = 2.d0*cell_infi1%v - v_int
       u_ext =-u_int
       v_ext =-v_int
       T_ext = T_int
       
    case(CL_FLUX_NULvsPAROI_ADIAB, CL_SYMETRIEvsPAROI_ADIAB) 
       u_ext = 2.d0*cell_infi2%u - u_int
       v_ext = 2.d0*cell_infi2%v - v_int
       T_ext = T_int

    case(CL_ENTREE_SUPERvsFLUX_NUL, CL_ENTREE_SUPERvsSYMETRIE, CL_ENTREE_SUPERvsINTERBLOC)
       u_ext = cell_infi1%u
       v_ext = cell_infi1%v
       T_ext = cell_infi1%T
    case(CL_FLUX_NULvsENTREE_SUPER, CL_SYMETRIEvsENTREE_SUPER, CL_INTERBLOCvsENTREE_SUPER) 
       u_ext = cell_infi2%u
       v_ext = cell_infi2%v
       T_ext = cell_infi2%T
    case(CL_ENTREE_SUPERvsPAROI_ADIAB)
       u_ext =-cell_infi1%u
       v_ext =-cell_infi1%v
       T_ext = cell_infi1%T
    case(CL_PAROI_ADIABvsENTREE_SUPER)
       u_ext =-cell_infi2%u
       v_ext =-cell_infi2%v
       T_ext = cell_infi2%T
    case(CL_ENTREE_SUPERvsPAROI)
       u_ext =-cell_infi1%u
       v_ext =-cell_infi1%v
       T_ext = 2.d0*cell_infi2%T - T_int
    case(CL_PAROIvsENTREE_SUPER)
       u_ext =-cell_infi2%u
       v_ext =-cell_infi2%v
       !T_ext = 2.d0*cell_infi1%T - T_int
       T_ext = 2.d0*coin%patch1%bc%Tw(1) - T_int

    case(CL_PAROIvsINJECTION, CL_INJECTIONvsPAROI, &
         CL_PAROI_ADIABvsINJECTION, CL_INJECTIONvsPAROI_ADIAB)
       u_ext =-u_int
       v_ext =-v_int
       T_ext = 2.d0*933.d0 - T_int  ! pas bon
       
    case default
!!!    !!! ToDo a remettre !!!
       if (coin%patch1%bc%condition == "injection_carbone") then
          k = l_int ! a changer 
          u_paroi = coin%patch1%bc%u_injection(k)*coin%n1(1)
          v_paroi = coin%patch1%bc%u_injection(k)*coin%n1(2)
                   
          u_ext = 2.d0*u_paroi - u_int
          v_ext = 2.d0*v_paroi - v_int
          T_ext = 2.d0*coin%patch1%bc%Tw(k) - T_int
          Ci(:,l_ext,m_ext)   = 2.d0*coin%patch1%bc%Ci_paroi(:,k)   - Ci(:,l_int,m_int)
          b%rho(l_ext,m_ext)  = 2.d0*coin%patch1%bc%rho_injection(k)- b%rho(l_int,m_int)
          b%D(l_ext,m_ext)    = 2.d0*coin%patch1%bc%D_paroi(k)      - b%D(l_int,m_int)
          b%K(l_ext,m_ext)    = 2.d0*coin%patch1%bc%K_paroi(k)      - b%K(l_int,m_int)
          b%mu(l_ext,m_ext)   = 2.d0*coin%patch1%bc%mu_paroi(k)     - b%mu(l_int,m_int)
          b%Hi(:,l_ext,m_ext) = 2.d0*coin%patch1%bc%Hi_paroi(:,k)   - b%Hi(:,l_int,m_int)
          b%y(l_ext,m_ext)    = 1.d0 ! b%rho1(l,m) * b%z(l,m) / b%rho(l,m)
       end if
       
       !print*, "cas non gere pour les coins en visqueux"
       !print*, coin%patch1%bc%condition, coin%patch2%bc%condition
       !call arret_code
    end select
    
  end subroutine coin_visqueux

  subroutine CL_coin_visqueux(b, indices, patch1, patch2)
    type (STR_BLOC), pointer :: b
    integer, dimension(6) :: indices
    type (STR_PATCH), pointer :: patch1, patch2

    integer :: ldeb, lfin, mdeb, mfin, lref, mref
    real*8, dimension(:,:), pointer :: u, v, T
    character(len=50) :: bc1, bc2
        
    ldeb = indices(1); lfin = indices(2)
    mdeb = indices(3); mfin = indices(4)
    lref = indices(5); mref = indices(6)

    u => b%u_x; v => b%u_y; T => b%T
    
    bc1 = patch1%bc%condition
    bc2 = patch2%bc%condition

    if ( trim(bc1) == "interproc" .or. &
         trim(bc2) == "interproc" ) then
       ! (u,v,T) deja synchronise entre les blocs 
       ! il y a peut etre un probleme (a verifier)
    else if ( trim(bc1) == "interbloc" .or. &
         trim(bc2) == "interbloc" ) then
       ! (u,v,T) deja synchronise entre les blocs 
       ! il y a peut etre un probleme (a verifier)         
!!$    else if ( (trim(bc1) == "paroi" .and. trim(bc2) == "paroi" ) &
!!$         .or. (trim(bc1) == "flux_nul" .and. trim(bc2) == "paroi") &
!!$         .or. (trim(bc1) == "paroi" .and. trim(bc2) == "flux_nul")) then
!!$       u(ldeb:lfin,mdeb:mfin) = -u(lref,mref)
!!$       v(ldeb:lfin,mdeb:mfin) = -v(lref,mref)
!!$       T(ldeb:lfin,mdeb:mfin) =  T(lref,mref)
       
!!$    else if ( trim(bc1) == "flux_nul" .and. trim(bc2) == "flux_nul") then
!!$       u(ldeb:lfin,mdeb:mfin) = u(lref,mref)
!!$       v(ldeb:lfin,mdeb:mfin) = v(lref,mref)
!!$       T(ldeb:lfin,mdeb:mfin) = T(lref,mref)
       
!    else if ( (trim(bc1) == "paroi_T_imposee" .and. trim(bc2) == "flux_nul") &
!         .or. (trim(bc1) == "paroi_T_imposee" .and. trim(bc2) == "symetrie") &
!         .or. (trim(bc1) == "paroi_T_imposee" .and. trim(bc2) == "paroi_T_imposee")) then
!       u(ldeb:lfin,mdeb:mfin) = -u(lref,mref)
!       v(ldeb:lfin,mdeb:mfin) = -v(lref,mref)
!       T(ldeb:lfin,mdeb:mfin) = 2.d0*patch1%bc%cell_infi%T - T(lref,mref)
       
!!$    else if ( (trim(bc1) == "flux_nul" .and. trim(bc2) == "paroi_T_imposee") &
!!$         .or. (trim(bc1) == "symetrie" .and. trim(bc2) == "paroi_T_imposee")) then
!!$       u(ldeb:lfin,mdeb:mfin) = -u(lref,mref)
!!$       v(ldeb:lfin,mdeb:mfin) = -v(lref,mref)
!!$       T(ldeb:lfin,mdeb:mfin) = 2.d0*patch2%bc%cell_infi%T - T(lref,mref)
   
!!$    else if (  trim(bc1) == "entree_supersonique" &
!!$         .and. trim(bc2) == "entree_supersonique" ) then    
!!$       u(ldeb:lfin,mdeb:mfin) = patch1%bc%cell_infi%u
!!$       v(ldeb:lfin,mdeb:mfin) = patch1%bc%cell_infi%v
!!$       T(ldeb:lfin,mdeb:mfin) = patch1%bc%cell_infi%T  
!!$       
!!$    else if ( (trim(bc1) == "entree_supersonique" .and. trim(bc2) == "flux_nul") &
!!$         .or. (trim(bc1) == "entree_supersonique" .and. trim(bc2) == "symetrie")) then
!!$       u(ldeb:lfin,mdeb:mfin) = patch1%bc%cell_infi%u
!!$       v(ldeb:lfin,mdeb:mfin) = patch1%bc%cell_infi%v
!!$       T(ldeb:lfin,mdeb:mfin) = patch1%bc%cell_infi%T
!!$       
!!$    else if ( (trim(bc1) == "flux_nul" .and. trim(bc2) == "entree_supersonique") &
!!$         .or. (trim(bc1) == "symetrie" .and. trim(bc2) == "entree_supersonique")) then
!!$       u(ldeb:lfin,mdeb:mfin) = patch2%bc%cell_infi%u
!!$       v(ldeb:lfin,mdeb:mfin) = patch2%bc%cell_infi%v
!!$       T(ldeb:lfin,mdeb:mfin) = patch2%bc%cell_infi%T
!!$       
!     else if ( (trim(bc1) == "paroi" .and. trim(bc2) == "entree_supersonique") &
!          .or. (trim(bc1) == "entree_supersonique" .and. trim(bc2) == "paroi")) then
!        u(ldeb:lfin,mdeb:mfin) = -u(lref,mref)
!        v(ldeb:lfin,mdeb:mfin) = -v(lref,mref)
!        T(ldeb:lfin,mdeb:mfin) =  T(lref,mref)
       
! !!! pas bon
!     else if ( (trim(bc1) == "paroi_global" .and. trim(bc2) == "entree_supersonique") &
!          .or. (trim(bc2) == "paroi_global" .and. trim(bc1) == "entree_supersonique") &
!          .or. (trim(bc1) == "paroi_global" .and. trim(bc2) == "flux_nul") &
!          .or. (trim(bc2) == "paroi_global" .and. trim(bc1) == "flux_nul")) then
!        u(ldeb:lfin,mdeb:mfin) = -u(lref,mref)
!        v(ldeb:lfin,mdeb:mfin) = -v(lref,mref)
!        T(ldeb:lfin,mdeb:mfin) =  T(lref,mref)
       
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! pas sur !!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else if ( (bc1=="paroi_T_imposee" .and. bc2=="paroi") &
         .or. (bc2=="paroi_T_imposee" .and. bc1=="paroi")) then
       u(ldeb:lfin,mdeb:mfin) = -u(lref,mref)
       v(ldeb:lfin,mdeb:mfin) = -v(lref,mref)
       T(ldeb:lfin,mdeb:mfin) =  T(lref,mref)

       !! a voir !!
    else if ( (trim(bc1) == "paroi_T_imposee" .and. trim(bc2) == "entree_supersonique") &
         .or. (trim(bc1) == "entree_supersonique" .and. trim(bc2) == "paroi_T_imposee")) then 
       u(ldeb:lfin,mdeb:mfin) = -u(lref,mref)
       v(ldeb:lfin,mdeb:mfin) = -v(lref,mref)
       T(ldeb:lfin,mdeb:mfin) =  T(lref,mref)
       

       ! pas sur 
       ! pas bon, il faut peut etre mettre la vitesse de la paroi
       ! pour ca, il faut passer la normale en argument
    else if ( (bc1=="paroi_entrainee_T_imposee" .and. bc2=="paroi") &
         .or. (bc2=="paroi_entrainee_T_imposee" .and. bc1=="paroi") &
         .or. (bc1=="paroi_entrainee" .and. bc2=="paroi") &
         .or. (bc2=="paroi_entrainee" .and. bc1=="paroi") ) then
       u(ldeb:lfin,mdeb:mfin) = -u(lref,mref)
       v(ldeb:lfin,mdeb:mfin) = -v(lref,mref)
       T(ldeb:lfin,mdeb:mfin) =  T(lref,mref)
    else if ( (bc1=="paroi_global_adiabatique" .and. bc2=="paroi_global_adiabatique") ) then
       u(ldeb:lfin,mdeb:mfin) = -u(lref,mref)
       v(ldeb:lfin,mdeb:mfin) = -v(lref,mref)
       T(ldeb:lfin,mdeb:mfin) =  T(lref,mref)
    else if ( (bc1=="paroi_global" .and. bc2=="paroi_global_adiabatique") &
         .or. (bc2=="paroi_global" .and. bc1=="paroi_global_adiabatique") &
         .or. (bc1=="paroi_flux_impose" .and. bc2=="paroi_global_adiabatique") & 
         .or. (bc2=="paroi_flux_impose" .and. bc1=="paroi_global_adiabatique") ) then
       u(ldeb:lfin,mdeb:mfin) = -u(lref,mref)
       v(ldeb:lfin,mdeb:mfin) = -v(lref,mref)
       T(ldeb:lfin,mdeb:mfin) =  T(lref,mref)
       
!!! pas bon !!!
    else if ( (bc1=="injection" .and. bc2=="paroi_T_imposee") &
         .or. (bc2=="injection" .and. bc1=="paroi_T_imposee") ) then
       u(ldeb:lfin,mdeb:mfin) = -u(lref,mref)
       v(ldeb:lfin,mdeb:mfin) = -v(lref,mref)
       T(ldeb:lfin,mdeb:mfin) =  2.d0*933.d0 - T(lref,mref)
    else if ( (bc1=="injection" .and. bc2=="paroi") &
         .or. (bc2=="injection" .and. bc1=="paroi") ) then
       u(ldeb:lfin,mdeb:mfin) = -u(lref,mref)
       v(ldeb:lfin,mdeb:mfin) = -v(lref,mref)
       T(ldeb:lfin,mdeb:mfin) = 2.d0*933.d0 - T(lref,mref)
    else if ( (bc1=="injection" .and. bc2=="paroi_global_adiabatique") &
         .or. (bc2=="injection" .and. bc1=="paroi_global_adiabatique") &
         .or. (bc1=="injection" .and. bc2=="paroi_global") &
         .or. (bc2=="injection" .and. bc1=="paroi_global") &
         .or. (bc1=="injection" .and. bc2=="flux_nul") &
         .or. (bc2=="injection" .and. bc1=="flux_nul")) then
       u(ldeb:lfin,mdeb:mfin) = -u(lref,mref)
       v(ldeb:lfin,mdeb:mfin) = -v(lref,mref)
       T(ldeb:lfin,mdeb:mfin) = 2.d0*933.d0 - T(lref,mref)  
       
!!! que faire entre une paroi entrainee et une paroi fixe (singularite sur la vitesse)
!!$    else if ( (bc1=="paroi_entrainee" .and. bc2=="paroi") &
!!$         .or. (bc2=="paroi_entrainee" .and. bc1=="paroi") ) then
!!$       u(ldeb:lfin,mdeb:mfin) = -u(lref,mref)
!!$       v(ldeb:lfin,mdeb:mfin) = -v(lref,mref)
!!$       T(ldeb:lfin,mdeb:mfin) =  T(lref,mref)

      
    else
       print*, "cas non gere pour les coins en visqueux"
       print*, patch1%bc%condition, patch2%bc%condition
       call arret_code
    end if
  end subroutine CL_coin_visqueux
  
end module CL_visqueux
