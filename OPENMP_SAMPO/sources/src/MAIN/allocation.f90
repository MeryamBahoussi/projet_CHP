module m_allocation
  use m_struct
  use M_MPI
  use parametres_globaux
  implicit none

contains

  subroutine allocation(tache)
    type (STR_TACHE), intent(in) ::  tache

    type (STR_SUPER_BLOC), pointer :: sb
    type (STR_BLOC)      , pointer :: b
    integer :: i, j, isb, ib

    do i = 1, tache%nb_membre
       isb = tache%membre(i)
       sb => acces_super_bloc(isb)%pointeur

       do j = 1, sb%nb_membre
          ib = sb%membre(j)
          if (bloc_est_local (ib)) then
             b=>acces_bloc(ib)%pointeur
             call allocation_bloc(sb, b)
          end if
       end do
    end do

  end subroutine allocation

  subroutine allocation_bloc(sb, b)
    use a_mettre_en_donnees

    type (STR_SUPER_BLOC), pointer :: sb
    type (STR_BLOC)      , pointer    :: b

    integer :: ld, lf, md, mf, l, m
    integer :: ne              ! nombre d'especes ou d elements
    integer :: neq             ! nombre d'equations
    integer :: neq_acoustique  ! nombre d'equations pour la partie acoustique

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    neq = b%neq
    neq_acoustique = b%neq-1


    if (b%nature==FLUIDE) then
       ne = b%ne

!!! valeurs aux centres des mailles (+ mailles fictives)
       ! var Euler
       allocate(b%U(  neq, ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%U_n(neq, ld-nMF:lf+nMF, md-nMF:mf+nMF))
       ! variables
       allocate(b%u_x(ld-nMF:lf+nMF, md-nMF:mf+nMF), b%u_y(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%E(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%h(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%y(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%p(ld-nMF:lf+nMF, md-nMF:mf+nMF), b%p_n(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%c(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%z(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%rho1(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%rho2(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%Htot(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%mu(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%T1(ld-nMF:lf+nMF, md-nMF:mf+nMF), b%T2(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       !debug
       allocate(b%dt_vv(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%dt_vT(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       ! Loi d'etat
       allocate(b%gamma_eq(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%pi_eq(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%q_eq(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%Cv1(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%gamma1(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       ! si on a plusieurs especes et ou elements
       allocate(b%c_esp(1:fluide1%nb_esp, ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%c_ele(1:fluide1%nb_ele, ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%Hi(1:fluide1%nb_esp, ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%D( ld-nMF:lf+nMF, md-nMF:mf+nMF))
       ! initialisation pour gerer tous les cas
       b%c_esp(:,:,:) = 1.d0
       b%c_ele(:,:,:) = 1.d0
       b%Hi(:,:,:)=0.d0; b%D(:,:) = 0.d0

       allocate(b%dm(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       allocate(b%cell(ld-nMF:lf+nMF, md-nMF:mf+nMF))

!!! valeurs aux centres des faces
       allocate(b%Fg_acoustique(neq_acoustique, ld-1:lf, md-1:mf)) ! flux suivant l
       allocate(b%Fd_acoustique(neq_acoustique, ld-1:lf, md-1:mf)) ! flux suivant m
       ! vitesse du son Lagrangienne
       allocate(b%a_l(ld-1:lf, md-1:mf), b%a_m(ld-1:lf, md-1:mf))
       ! pentes du solveur de Riemann
       allocate(b%Cm_l(ld-1:lf, md:mf), b%Cp_l(ld-1:lf, md:mf))
       allocate(b%Cm_m(ld:lf, md-1:mf), b%Cp_m(ld:lf, md-1:mf))

       !metrique en visqueux
       allocate(b%metrique_grad_l(4,ld-1:lf,md:mf))
       allocate(b%metrique_grad_m(4,ld:lf,md-1:mf))


       ! tableau de la pseudoinverse lors du calcul du gradient par moindre carrés
       ! ne dépend que de la géométrie
       !ATTENTION : la matrice change si le maillage se déforme...
       allocate(b%Mpinv(2*(b%st_grad-1),b%ld-(nMF-1):b%lf+(nMF-1),b%md-(nMF-1):b%mf+(nMF-1)))

       !allocation debuggage c<0
       allocate(b%debug_c(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       b%debug_c = 0

       ! creation d'un maillage dual
       allocate(b%aire_dual(ld-1:lf, md-1:mf))
       allocate(b%n_dir_l_dual(1:2,ld-2:lf, md-1:mf))
       allocate(b%n_dir_m_dual(1:2,ld-1:lf, md-2:mf))
       ! longueurs
       allocate(b%dl_l_dual(ld-2:lf, md-1:mf))
       allocate(b%dl_m_dual(ld-1:lf, md-2:mf))



       !!!! todo : on alloue dans tous les cas pour avoir rhoyCi bien
       call allocation_ordre2(b, sb%solveur_Riemann)

!!! second membre
       if (.not. b%implicite) allocate(b%sm_acoustique(neq_acoustique,ld:lf,md:mf))
    end if ! nature == FLUIDE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Variables communes aux blocs fluide et solide
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! valeurs aux centres des mailles (+ mailles fictives)
    allocate(b%rho(ld-nMF:lf+nMF, md-nMF:mf+nMF))
    allocate(b%T(ld-nMF:lf+nMF, md-nMF:mf+nMF))
    allocate(b%K(ld-nMF:lf+nMF, md-nMF:mf+nMF))
    allocate(b%Cv(ld-nMF:lf+nMF, md-nMF:mf+nMF))
    allocate(b%centre_cl(1:2, ld-nMF:lf+nMF, md-nMF:mf+nMF))
    allocate(b%aire(ld-nMF:lf+nMF, md-nMF:mf+nMF))
    allocate(b%zeros(ld-nMF:lf+nMF, md-nMF:mf+nMF))
    b%zeros(:,:) = 0.d0 ! pour unifier le codage
!!! valeurs au noeuds (+ mailles fictives)
    allocate(b%coord_nd(1:2, ld-(nMF+1):lf+nMF, md-(nMF+1):mf+nMF))
    allocate(b%coord_nd_old(1:2, ld-(nMF+1):lf+nMF, md-(nMF+1):mf+nMF))
    allocate(b%coord_nd_01x01(1:2, ld-(nMF+1):lf+nMF, md-(nMF+1):mf+nMF))
!!! valeurs aux centres des faces
    allocate(b%F(neq, ld-1:lf, md-1:mf)) ! flux
!!! valeurs aux centres des faces (+1 maille fictive)
    ! normales
    allocate(b%n_dir_l(1:2,ld-nMF:lf+(nMF-1), md-(nMF-1):mf+(nMF-1)))
    allocate(b%n_dir_m(1:2,ld-(nMF-1):lf+(nMF-1), md-nMF:mf+(nMF-1)))
    ! longueurs
    allocate(b%dl_l(ld-nMF:lf+(nMF-1), md-(nMF-1):mf+(nMF-1)))
    allocate(b%dl_m(ld-(nMF-1):lf+(nMF-1), md-nMF:mf+(nMF-1)))

!!! second membre
    allocate(b%sm(neq,ld:lf,md:mf))

!!! derivees de la pression (pas seulement en implicite car necessaire le gaz reel + roe explicite)
    if (b%nature == FLUIDE) then
       ! on alloue plus grand les derivees de p parce qu'on les calcule dans le decodage (et on decode 2 rangees de mailles fictives)
       allocate(b%deriv_p(ld-nMF:lf+nMF, md-nMF:mf+nMF))
       do m = md-nMF, mf+nMF
          do l = ld-nMF, lf+nMF
             allocate(b%deriv_p(l,m)%d_drho1z1Ci(1:b%ne))
          end do
       end do
    end if

    if (b%implicite) then
       call allocation_implicite(sb, b)
    end if

    if (b%maillage%move) then
!!! valeurs aux centres des mailles (+ mailles fictives)
       allocate(b%aire_new(ld-nMF:lf+nMF, md-nMF:mf+nMF))
!!! valeurs aux centres des faces (+1 maille fictive)
       allocate(b%dx_dt(ld-2:lf+1, md-1:mf  ))
       allocate(b%dy_dt(ld-1:lf  , md-2:mf+1))
       b%dx_dt = 0.d0
       b%dy_dt = 0.d0
    end if

  end subroutine allocation_bloc

!!$ld-1:lf, md-1:mf  ! face
!!$
!!$ld-2:lf+2, md-2:mf+2  ! centre + MF
!!$ld-3:lf+2, md-3:mf+2  ! noeuds + MF


  subroutine allocation_implicite(sb, b)
    use a_mettre_en_donnees

    type (STR_SUPER_BLOC), pointer :: sb
    type (STR_BLOC)      , pointer    :: b

    integer :: l, m, ld, lf, md, mf, dim_bloc

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf


    ! choix du stencil de la matrice
    b%mat%stencil = 5
    
    
    dim_bloc = b%neq  ! on inverse tout
    allocate(b%deriv_T(ld-1:lf+1, md-1:mf+1))
    do m = md-1, mf+1
       do l = ld-1, lf+1
          allocate(b%deriv_T(l,m)%d_drho1z1Ci(1:b%ne))
       end do
    end do
    allocate(b%dV_dU(1:dim_bloc,1:dim_bloc,ld-1:lf+1, md-1:mf+1))

    allocate(b%grad_e(b%neq,b%neq,ld-1:lf+1,md-1:mf+1))
    allocate(b%grad_f(b%neq,b%neq,ld-1:lf+1,md-1:mf+1))

    ! pour aller plus vite en stationnaire on ne tient pas compte des coins
    if (b%visqueux .and. sb%sol_stat%actif .eqv. .false.) b%mat%stencil = 9



    if (b%solveur_implicite == GS) then
       call allocation_matrice1d(b, DIR_X, 1, b%matrice)
    else
       call allocation_matrice(b, dim_bloc, b%mat)
    end if

    if (sb%solveur_implicite==MKL .and. numproc==b%proc_ini) call allocation_matrice_CSR(b,dim_bloc,b%mat_csr)

    if (implicite==implicite_Up) then
       allocate(b%sm_acoustique(dim_bloc,ld:lf,md:mf))
    else if (implicite==implicite_total) then
       ! pour pouvoir utiliser la routine sm_euler
       allocate(b%sm_acoustique(dim_bloc+1,ld:lf,md:mf))
    end if

    allocate(b%X_imp(dim_bloc,ld-nMF:lf+nMF,md-nMF:mf+nMF))
    allocate(b%X_old(dim_bloc,ld-nMF:lf+nMF,md-nMF:mf+nMF))
    allocate(b%X_0(dim_bloc,ld-nMF:lf+nMF,md-nMF:mf+nMF))

    ! pour l'equation sur l energie
    if (b%eq_energie == 0) then
       b%mat_scalaire%stencil=5
       if (b%solveur_implicite == GS) then
          call allocation_matrice1d(b, DIR_X, 1, b%matrice)
       else if (b%solveur_implicite == PETSc) then
          call allocation_matrice(b, 1, b%mat_scalaire)
       end if
       allocate(b%sm_scalaire(1,ld:lf,md:mf))
    end if


  end subroutine allocation_implicite

  subroutine allocation_matrice1d(b, direction, dim, mat)
    type (STR_BLOC), pointer ::  b
    integer, intent(in) :: direction, dim
    type (matrice_tridiagonale), intent(inout) :: mat

    integer :: ld, lf, md, mf

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    mat%dim_bloc=dim

    if (direction==DIR_X) then  ! balayage dans le sens des X
       b%matrice%dim_mat = mf
       allocate(mat%a(dim,dim,md:mf))
       allocate(mat%b(dim,dim,md:mf))
       allocate(mat%c(dim,dim,md:mf))
    else if (direction==DIR_Y) then
       b%matrice%dim_mat = lf
       allocate(mat%a(dim,dim,ld:lf))
       allocate(mat%b(dim,dim,ld:lf))
       allocate(mat%c(dim,dim,ld:lf))
    end if

  end subroutine allocation_matrice1d

  subroutine allocation_matrice(b, dim, mat)
    type (STR_BLOC), pointer ::  b
    integer, intent(in) :: dim
    type (STR_MATRICE), intent(inout) :: mat

    integer :: ld, lf, md, mf

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf

    mat%dim_bloc=dim
    allocate(mat%A(dim,dim,ld:lf,md:mf))
    allocate(mat%B(dim,dim,ld:lf,md:mf))
    allocate(mat%C(dim,dim,ld:lf,md:mf))
    allocate(mat%D(dim,dim,ld:lf,md:mf))
    allocate(mat%E(dim,dim,ld:lf,md:mf))

    allocate(mat%AD(dim,dim,ld:lf,md:mf))
    allocate(mat%CD(dim,dim,ld:lf,md:mf))
    allocate(mat%AE(dim,dim,ld:lf,md:mf))
    allocate(mat%CE(dim,dim,ld:lf,md:mf))

  end subroutine allocation_matrice

  subroutine allocation_matrice_CSR(b, dim, mat_csr)
    use m_outils_CSR
    implicit none

    type (STR_BLOC), pointer    ::  b
    integer        , intent(in) :: dim  ! dimension de l'inconnue
    type (STR_CSR) , pointer    :: mat_csr
    ! variables locales
    integer              :: ld, lf, md, mf
    integer, allocatable :: graphe(:,:), Nombre_non_zero(:)

    ld = b%ld_tot; lf = b%lf_tot
    md = b%md_tot; mf = b%mf_tot

    allocate(mat_csr)

    mat_csr%nb_rows = ((mf-md+1)*(lf-ld+1)*dim)

    allocate( graphe(1:dim*9, 1: mat_csr%nb_rows ), Nombre_non_zero(1: mat_csr%nb_rows ) )

    call Structure_graphe_connectivites( lf-ld+1, mf-md+1, graphe,Nombre_non_zero, mat_csr%nb_rows, dim, .false.  )

    mat_csr%nNonZeros = sum(Nombre_non_zero(:))

    allocate(mat_csr%rowIndex(1: mat_csr%nb_rows +1))
    allocate(mat_csr%columns(mat_csr%nNonZeros))
    allocate(mat_csr%Values(mat_csr%nNonZeros))

    call structure_CSR(graphe,Nombre_non_zero,mat_csr%columns, mat_csr%rowIndex,mat_csr%nNonZeros,mat_csr%nb_rows, dim*9 )

    deallocate(graphe,Nombre_non_zero)

  end subroutine allocation_matrice_CSR

  ! allocation des tableaux pour l'ordre 2
  subroutine allocation_ordre2(b, solveur_Riemann)
    use parametres_fluide

    type (STR_BLOC), pointer    :: b
    integer        , intent(in) :: solveur_Riemann

    integer :: ld, lf, md, mf, ne

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf

    ! on alloue plus grand pour pouvoir utiliser le meme tableau peut importe
    ! la direction (X ou Y)
    ! on reconstruit (rhoyci)_1,ne, rhoy2, ux, uy, p et z
    ne = b%ne
    allocate(b%ordre2%rhoyciL(ne, ld-1:lf+1, md-1:mf+1))
    allocate(b%ordre2%rhoyciR(ne, ld-1:lf+1, md-1:mf+1))
    allocate(b%ordre2%rhoy2L(ld-1:lf+1, md-1:mf+1))
    allocate(b%ordre2%rhoy2R(ld-1:lf+1, md-1:mf+1))
    allocate(b%ordre2%uxL(ld-1:lf+1, md-1:mf+1))
    allocate(b%ordre2%uxR(ld-1:lf+1, md-1:mf+1))
    allocate(b%ordre2%uyL(ld-1:lf+1, md-1:mf+1))
    allocate(b%ordre2%uyR(ld-1:lf+1, md-1:mf+1))
    allocate(b%ordre2%pL(ld-1:lf+1, md-1:mf+1))
    allocate(b%ordre2%pR(ld-1:lf+1, md-1:mf+1))
    allocate(b%ordre2%zL(ld-1:lf+1, md-1:mf+1))
    allocate(b%ordre2%zR(ld-1:lf+1, md-1:mf+1))
    ! on en deduit E et c pour le solveur de Riemann
    allocate(b%ordre2%EL(ld-1:lf+1, md-1:mf+1))
    allocate(b%ordre2%ER(ld-1:lf+1, md-1:mf+1))
    allocate(b%ordre2%cL(ld-1:lf+1, md-1:mf+1))
    allocate(b%ordre2%cR(ld-1:lf+1, md-1:mf+1))
  end subroutine allocation_ordre2

  subroutine deallocation(tache)
    type (STR_TACHE), intent(in) ::  tache

    integer                        :: i, j, isb, ib
    type (STR_SUPER_BLOC), pointer :: sb
    type (STR_BLOC)      , pointer :: b

    do i = 1, tache%nb_membre
       isb = tache%membre(i)
       sb => acces_super_bloc(isb)%pointeur

       do j = 1, sb%nb_membre
          ib = sb%membre(j)
          if (bloc_est_local (ib)) then
             b=>acces_bloc(ib)%pointeur
             deallocate(b)
          end if
       end do
    end do

  end subroutine deallocation

end module m_allocation
