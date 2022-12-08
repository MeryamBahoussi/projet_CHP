module outils_maillage
  use m_struct
  use m_MPI
  use parametres_globaux, only : nMF

  implicit none 
  
contains
  
  function distance(A,B)
    implicit none
    real*8, dimension(2) :: A, B
    real*8 :: distance
    distance = dsqrt( (A(1)-B(1))**2 + (A(2)-B(2))**2 )
  end function distance

  ! produit mixture entre deux vecteurs et ez
  ! (u^v).ez
  function produitMixteAvecEz(u, v)
    real*8, dimension(2) :: u, v
    real*8 :: produitMixteAvecEz
    
    produitMixteAvecEz = u(1)*v(2)-u(2)*v(1)
  end function produitMixteAvecEz

  function dist_pt_droite(a,b,c,X)
    implicit none 
    real*8 :: X(2)
    real*8 :: a,b,c
    real*8 :: dist_pt_droite
    if ( dsqrt( a**2 + b**2 ) /= 0.d0) then
       dist_pt_droite = abs(a*X(1)+b*X(2)+c)/( dsqrt( a**2 + b**2 ) )
    else
       print*, "pb dist_pt_droite"
    end if
  end function dist_pt_droite

  subroutine calc_coeff_droite_w_2pts(a,b,c,M,N)
    implicit none 
    real*8,intent(in) :: M(1:2),N(1:2)
    real*8,intent(out) :: a,b,c

    ! Calcul des coefs a, b et c des droites ax+by+c=0 
    a=M(2)-N(2)
    b=N(1)-M(1)
    c=-(a*M(1)+b*M(2))
  end subroutine calc_coeff_droite_w_2pts

  subroutine Coordonnees_pt_sur_droite(O,vect_dir,distance_,M)
    implicit none   
    real*8,intent(in) ::  vect_dir(2),O(2),distance_
    real*8,intent(out) :: M(2)
    
    ! Calcul des coordonnées du pt M tq vect(OM) = distance * vect(vect_dir) 
    ! M = point distant de la longeueur =distance du pt O et qui appartient 
    ! a la droite passant par O ddont la vecteur direceteur est vect(vect_dir) 
    
    M(1) = O(1) + distance_ * vect_dir(1)
    M(2) = O(2) + distance_ * vect_dir(2)
  
  end subroutine Coordonnees_pt_sur_droite
  
  subroutine symetrie_axiale_d_un_pt(A, B, pt, pt_sym)
    implicit none 
    real*8, intent(in) :: A(1:2), B(1:2), pt(1:2)
    real*8, intent(inout) :: pt_sym(1:2)
    real*8 :: H(1:2) ! projete de pt sur la droite AB
    real*8 :: dx, dy

    if (A(1) .eq. B(1) .and.  A(2) .eq. B(2)) then
       print*, "erreur dans le calcul du symetrique par rapport a la droite (AB)"
       print*, "A = B"
       print*, "ERREUR : outils_maillage : symetrie_axiale_d_un_pt"
    else

       dx = (B(1)-A(1)); dy = (B(2)-A(2))

       H(1) = pt(1)*dx**2 + (pt(2)+A(2))*dx*dy + A(1)*dy**2
       H(1) = H(1) / (dx**2+dy**2)

       H(2) = A(2)*dx**2 - (A(1)-pt(1))*dx*dy + pt(2)*dy**2
       H(2) = H(2) / (dx**2+dy**2)

       pt_sym(1) = 2*H(1)-pt(1)
       pt_sym(2) = 2*H(2)-pt(2)

    end if

  end subroutine symetrie_axiale_d_un_pt
  
  ! trouve les noeuds d'une face (celle definie dans le fichier de donnees)
  subroutine find_noeud_bord(b, iface, noeud, bd, bf, nd_bord)
    type (STR_BLOC)           , pointer     :: b
    integer                   , intent(in)  :: iface
    real*8, dimension(2,b%ld-(nMF+1):b%lf+nMF,b%md-(nMF+1):b%mf+nMF), intent(in) :: noeud
    integer                   , intent(in)  :: bd, bf
    real*8, dimension(2,bd-(nMF+1):bf+nMF), intent(out) :: nd_bord

    integer :: iface_opp, me_bord, me_bord_opp
    integer :: ib, k, p_adj, taille, ibord
    integer :: STATINFO, MPI_FLAG=75
    integer :: MPI_status(MPI_STATUS_SIZE)
    integer, allocatable :: liste_procs(:)

    ib = b%numero
    call find_procs_jusquau_bord(ib, iface, me_bord, liste_procs)

    taille = 2*(bf+nMF-(bd-(nMF+1))+1)

    ! le proc me_bord envoie l'info a tous ceux qui sont entre lui et l'autre bord
    if (numproc==me_bord) then

       ! on trouve l'autre bord et l'indice du bord
       select case(iface)
       case(1)
          iface_opp = 2
          ibord = b%ld-1
          nd_bord(1:2,bd-(nMF+1):bf+nMF) = noeud(:,ibord,:)
       case(2)
          ibord = b%lf
          iface_opp = 1
          nd_bord(1:2,bd-(nMF+1):bf+nMF) = noeud(:,ibord,:)
       case(3)
          ibord = b%md-1
          iface_opp = 4
          nd_bord(1:2,bd-(nMF+1):bf+nMF) = noeud(:,:,ibord)
       case(4)
          ibord = b%mf
          iface_opp = 3
          nd_bord(1:2,bd-(nMF+1):bf+nMF) = noeud(:,:,ibord)
       end select

       call find_procs_jusquau_bord(ib, iface_opp, me_bord_opp, liste_procs)
       k = 1
       if (me_bord_opp/=me_bord) then
          do 
             p_adj = liste_procs(k)
             call MPI_SEND( nd_bord(1:2,bd-(nMF+1):bf+nMF), taille, MPI_REAL8, &
                  p_adj, MPI_FLAG, MPI_COMM_WORLD, STATINFO )
             k = k+1
             if (p_adj == me_bord_opp) exit
          end do
       end if
    else
       ! reception des points du bord
       call MPI_RECV( nd_bord(1:2,bd-(nMF+1):bf+nMF), taille, MPI_REAL8, &
            me_bord, MPI_FLAG, MPI_COMM_WORLD, MPI_status, STATINFO )     
    end if

  end subroutine find_noeud_bord

!!! on trouve le proc du bloc sur le bord
!!! comment definir le bord ??
  subroutine find_procs_jusquau_bord(ib_start, iface, numproc_bord, liste_p)
    integer, intent(in)  :: ib_start ! bloc de depart
    integer, intent(in)  :: iface
    integer, intent(out) :: numproc_bord
    integer, allocatable, intent(out) :: liste_p(:)

    integer :: ib, ip, ib_adj, k
    type (STR_BLOC) , pointer :: b
    type (STR_PATCH), pointer :: patch

    allocate(liste_p(size(acces_bloc,1)))
    liste_p = -1
    k = 1
    ib = ib_start
    boucle_blocs : do  ! on boucle tant qu on n a pas trouve le bord 
       b => acces_bloc(ib)%pointeur
       do ip = 1, b%face(iface)%nb_patch
          patch => b%face(iface)%patch(ip)
          ib_adj = patch%ib_adj

          if ( ib_adj == -1 .or. &
               (index(patch%bc%condition, "itf") /=0) .or. &
               (index(patch%bc%condition, "periodique") /=0) ) then
             numproc_bord = b%numproc
             exit boucle_blocs
          else
             if ( b%nature /= acces_bloc(ib_adj)%pointeur%nature .or. &
                  b%num_input /= acces_bloc(ib_adj)%pointeur%num_input ) then
                numproc_bord = b%numproc
                exit boucle_blocs
             else
                ib = ib_adj
                liste_p(k) = acces_bloc(ib)%pointeur%numproc
                k = k+1
             end if
          end if
       end do
    end do boucle_blocs

  end subroutine find_procs_jusquau_bord

end module outils_maillage
