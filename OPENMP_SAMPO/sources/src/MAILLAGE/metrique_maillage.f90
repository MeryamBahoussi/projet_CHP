module metrique_maillage
  use m_struct
  use m_MPI             , only : numproc, arret_code
  use outils_solveur    , only : inverse_mat22, QR_decomposition
  use outils_maillage   , only : distance, produitMixteAvecEz, find_noeud_bord
  use parametres_globaux, only : DIR_X, DIR_Y, nMF
 
  implicit none 

contains

  subroutine calcul_centre_cellule(b, noeud, centre)

    type (STR_BLOC), pointer :: b
    real*8, dimension(:,:,:), pointer :: noeud, centre

    integer :: l, m
    real*8, dimension(2) :: pI, pM
    real*8 :: a1, a2, b1, b2

    !Coordonnées centre des mailles (isobarycentre des 4 sommets)
    ! do m = b%md-nMF, b%mf+nMF
    !    do l = b%ld-nMF, b%lf+nMF
    !       centre(:,l,m)=0.25d0*(noeud(:,l-1,m) + noeud(:,l,m) +&
    !            noeud(:,l,m-1) + noeud(:,l-1,m-1) )
    !    end do
    ! end do


    ! centre de masse
    do m = b%md-nMF, b%mf+nMF
       do l = b%ld-nMF, b%lf+nMF
          ! pI = intersection des diagonales
          if(noeud(1,l-1,m-1) == noeud(1,l,m))then
             pI(1) = noeud(1,l-1,m-1)
             a2 = (noeud(2,l-1,m) - noeud(2,l,m-1))/(noeud(1,l-1,m) - noeud(1,l,m-1))
             b2 = noeud(2,l,m-1) - a2*noeud(1,l,m-1)
             pi(2) = a2*pI(1) + b2
          else if(noeud(1,l,m-1) == noeud(1,l-1,m))then
             pI(1) = noeud(1,l,m-1)
             a1 = (noeud(2,l-1,m-1) - noeud(2,l,m))/(noeud(1,l-1,m-1) - noeud(1,l,m))
             b1 = noeud(2,l-1,m-1) - a1*noeud(1,l-1,m-1)
             pI(2) = a1*pI(1) + b1
          else
             a1 = (noeud(2,l-1,m-1) - noeud(2,l,m))/(noeud(1,l-1,m-1) - noeud(1,l,m))
             b1 = noeud(2,l-1,m-1) - a1*noeud(1,l-1,m-1)
             a2 = (noeud(2,l-1,m) - noeud(2,l,m-1))/(noeud(1,l-1,m) - noeud(1,l,m-1))
             b2 = noeud(2,l,m-1) - a2*noeud(1,l,m-1)

             pI(1) = (b2-b1)/(a1-a2)
             pI(2) = a1*pI(1) + b1
          end if

          ! pM = isobarycentre des 4 sommets
          pM(:) = 0.25d0*(noeud(:,l-1,m) + noeud(:,l,m) +&
               noeud(:,l,m-1) + noeud(:,l-1,m-1) )

          ! centre de masse (IG = 4/3 IM) th de Wittenbauer
          centre(:,l,m) = (4.d0*pM(:) - pI(:))/3.d0
       end do
    end do

  end subroutine calcul_centre_cellule

  subroutine calcul_aire(b, nd, aire)
    !use variables_debug, only : cmp_manu

    type (STR_BLOC), pointer :: b
    real*8, dimension(:,:,:), pointer :: nd  ! coordonnee des sommets
    real*8, dimension(:,:), pointer :: aire

    integer :: l, m
    !    real*8 :: a_, b_, c_, high_sup, high_inf, base
    real*8 :: x0(2)

!!! On utilise la formule Aire = 1/2 * somme_i [x_i, x_i+1, ez]
!!! ou x_i sont les sommets du polygone, ez=(0,0,1)
!!! Ici, on se base pas par rapport à l'origine mais par rapport à nd(:,l-1,m-1)
    do m = b%md-nMF, b%mf+nMF
       do l = b%ld-nMF, b%lf+nMF
          x0 = nd(:,l-1,m-1)
          aire(l,m) = 0.5d0 * ( &
               produitMixteAvecEz(nd(:,l,m-1)-x0, nd(:,l,m)-x0) + &
               produitMixteAvecEz(nd(:,l,m)-x0, nd(:,l-1,m)-x0) )
       end do
    end do

!!$    if (cmp_manu) then
!!$       !Calcul d'une aire quelconque, décomposition en deux triangles
!!$       do m = b%md-2, b%mf+2
!!$          do l = b%ld-2, b%lf+2
!!$             call calc_coeff_droite_w_2pts(a_,b_,c_,nd(:,l,m),nd(:,l-1,m-1))
!!$             high_sup = dist_pt_droite(a_,b_,c_,nd(:,l-1,m))
!!$             high_inf = dist_pt_droite(a_,b_,c_,nd(:,l,m-1))
!!$             base = distance(nd(:,l,m), nd(:,l-1,m-1))
!!$             aire(l,m) = 0.5d0 * base * ( high_sup + high_inf )
!!$          end do
!!$       end do
!!$    end if

  end subroutine calcul_aire

  subroutine calcul_aire_dual(b, nd, aire)
    type (STR_BLOC), pointer :: b
    real*8, dimension(:,:,:), pointer :: nd  ! coordonnee des centres
    real*8, dimension(:,:), pointer :: aire

    integer :: l, m
    real*8 :: x0(2)

!!! On utilise la formule Aire = 1/2 * somme_i [x_i, x_i+1, ez]
!!! ou x_i sont les sommets du polygone, ez=(0,0,1)
!!! Ici, on se base pas par rapport à l'origine mais par rapport à nd(:,l-1,m-1)
    do m = b%md-1, b%mf
       do l = b%ld-1, b%lf
          x0 = nd(:,l,m)
          aire(l,m) = 0.5d0 * ( &
               produitMixteAvecEz(nd(:,l+1,m)-x0, nd(:,l+1,m+1)-x0) + &
               produitMixteAvecEz(nd(:,l+1,m+1)-x0, nd(:,l,m+1)-x0) )
       end do
    end do

  end subroutine calcul_aire_dual

  subroutine calcul_normales(b, noeud)

    type (STR_BLOC), pointer :: b
    real*8, dimension(:,:,:), pointer :: noeud

    integer :: l,m
    real*8, dimension(:,:,:), pointer :: n_dir_l, n_dir_m
    real*8, dimension(:,:), pointer :: dl_l, dl_m

    n_dir_l => b%n_dir_l; n_dir_m => b%n_dir_m
    dl_l => b%dl_l; dl_m => b%dl_m

    do m = b%md-(nMF-1), b%mf+(nMF-1)
       do l = b%ld-nMF, b%lf+(nMF-1)
          ! longueurs
          dl_l(l,m) = distance(noeud(:,l,m), noeud(:,l,m-1))

          !Normale direction l
          n_dir_l(1,l,m) = ( noeud(2,l,m)-noeud(2,l,m-1) ) / dl_l(l,m)
          n_dir_l(2,l,m) = - ( noeud(1,l,m)-noeud(1,l,m-1) ) / dl_l(l,m)
       end do
    end do

    do m = b%md-nMF, b%mf+(nMF-1)
       do l = b%ld-(nMF-1), b%lf+(nMF-1)
          ! longueurs
          dl_m(l,m) = distance(noeud(:,l,m), noeud(:,l-1,m))

          !Normale direction m
          n_dir_m(1,l,m) = - ( noeud(2,l,m) - noeud(2,l-1,m) ) / dl_m(l,m)
          n_dir_m(2,l,m) = ( noeud(1,l,m) - noeud(1,l-1,m) ) / dl_m(l,m)
       end do
    end do

!!$    do m = b%md-1, b%mf+1
!!$       do l = b%ld-1, b%lf+1
!!$          ! longueurs
!!$          dl_l(l,m) = distance(noeud(:,l,m), noeud(:,l,m-1))
!!$          dl_m(l,m) = distance(noeud(:,l,m), noeud(:,l-1,m))
!!$
!!$          !Normale direction l
!!$          n_dir_l(1,l,m) = ( noeud(2,l,m)-noeud(2,l,m-1) ) / dl_l(l,m)
!!$          n_dir_l(2,l,m) = - ( noeud(1,l,m)-noeud(1,l,m-1) ) / dl_l(l,m)
!!$
!!$          !Normale direction m
!!$          n_dir_m(1,l,m) = - ( noeud(2,l,m) - noeud(2,l-1,m) ) / dl_m(l,m)
!!$          n_dir_m(2,l,m) = ( noeud(1,l,m) - noeud(1,l-1,m) ) / dl_m(l,m)
!!$       end do
!!$    end do


  end subroutine calcul_normales

  subroutine calcul_normales_dual(b, noeud)

    type (STR_BLOC), pointer :: b
    real*8, dimension(:,:,:), pointer :: noeud

    integer :: l,m
    real*8, dimension(:,:,:), pointer :: n_dir_l, n_dir_m
    real*8, dimension(:,:), pointer :: dl_l, dl_m

    n_dir_l => b%n_dir_l_dual; n_dir_m => b%n_dir_m_dual
    dl_l => b%dl_l_dual; dl_m => b%dl_m_dual

    do m = b%md-1, b%mf
       do l = b%ld-2, b%lf
          ! longueurs
          dl_l(l,m) = distance(noeud(:,l+1,m+1), noeud(:,l+1,m))

          !Normale direction l
          n_dir_l(1,l,m) = ( noeud(2,l+1,m+1)-noeud(2,l+1,m) ) / dl_l(l,m)
          n_dir_l(2,l,m) = - ( noeud(1,l+1,m+1)-noeud(1,l+1,m) ) / dl_l(l,m)
       end do
    end do

    do m = b%md-2, b%mf
       do l = b%ld-1, b%lf
          ! longueurs
          dl_m(l,m) = distance(noeud(:,l+1,m+1), noeud(:,l,m+1))

          !Normale direction m
          n_dir_m(1,l,m) = - ( noeud(2,l+1,m+1) - noeud(2,l,m+1) ) / dl_m(l,m)
          n_dir_m(2,l,m) = ( noeud(1,l+1,m+1) - noeud(1,l,m+1) ) / dl_m(l,m)
       end do
    end do

  end subroutine calcul_normales_dual

  subroutine calc_metrique_gradient(b)
    type (STR_BLOC), pointer :: b

    integer :: l, m

    ! metrique pour le calcul du gradient aux faces
    ! on range dans l'ordre suivant : X_x, X_y, Y_x, Y_y

    !selon x
    do m=b%md, b%mf
       do l=b%ld-1, b%lf
          b%metrique_grad_l(1,l,m) = b%n_dir_l(1,l,m)*b%dl_l(l,m)
          b%metrique_grad_l(2,l,m) = b%n_dir_l(2,l,m)*b%dl_l(l,m)
          b%metrique_grad_l(3,l,m) = .25d0 * ( b%n_dir_m(1,l  ,m  )*b%dl_m(l  ,m  ) + &
               b%n_dir_m(1,l+1,m  )*b%dl_m(l+1,m  ) + &
               b%n_dir_m(1,l  ,m-1)*b%dl_m(l  ,m-1) + &
               b%n_dir_m(1,l+1,m-1)*b%dl_m(l+1,m-1))
          b%metrique_grad_l(4,l,m) = .25d0 * ( b%n_dir_m(2,l  ,m  )*b%dl_m(l  ,m  ) + &
               b%n_dir_m(2,l+1,m  )*b%dl_m(l+1,m  ) + &
               b%n_dir_m(2,l  ,m-1)*b%dl_m(l  ,m-1) + &
               b%n_dir_m(2,l+1,m-1)*b%dl_m(l+1,m-1))
       enddo
    enddo

    ! selon y
    do m=b%md-1, b%mf
       do l=b%ld, b%lf
          b%metrique_grad_m(1,l,m) = 0.25d0 * ( b%n_dir_l(1,l  ,m  )*b%dl_l(l  ,m  ) + &
               b%n_dir_l(1,l  ,m+1)*b%dl_l(l  ,m+1) + &
               b%n_dir_l(1,l-1,m  )*b%dl_l(l-1,m  ) + &
               b%n_dir_l(1,l-1,m+1)*b%dl_l(l-1,m+1))
          b%metrique_grad_m(2,l,m) = 0.25d0 * ( b%n_dir_l(2,l  ,m  )*b%dl_l(l  ,m  ) + &
               b%n_dir_l(2,l  ,m+1)*b%dl_l(l  ,m+1) + &
               b%n_dir_l(2,l-1,m  )*b%dl_l(l-1,m  ) + &
               b%n_dir_l(2,l-1,m+1)*b%dl_l(l-1,m+1))
          b%metrique_grad_m(3,l,m) = b%n_dir_m(1,l,m)*b%dl_m(l,m)
          b%metrique_grad_m(4,l,m) = b%n_dir_m(2,l,m)*b%dl_m(l,m)
       enddo
    enddo

  end subroutine calc_metrique_gradient

  ! On trouve la position de chaque noeud sur le carre unitaire
  ! utile pour deplacer le maillage (reconstruction des points interieurs) en //
  subroutine coordonees_maillage_01x01(b, noeud, nd_01x01)
    type (STR_BLOC), pointer :: b
    real*8, dimension(:,:,:), pointer :: noeud
    real*8, dimension(:,:,:), pointer :: nd_01x01

    integer :: l, m
    real*8, dimension(:,:), allocatable :: nd_, nf_
    real*8, dimension(2)                :: nd(2), nf(2)

!!! en l
    allocate(nd_(2,b%md-(nMF+1):b%mf+nMF))
    allocate(nf_(2,b%md-(nMF+1):b%mf+nMF))
    call find_noeud_bord(b, 1, noeud, b%md, b%mf, nd_(:,:))
    call find_noeud_bord(b, 2, noeud, b%md, b%mf, nf_(:,:))
    do l = b%ld-(nMF+1), b%lf+nMF
       do m = b%md-(nMF+1), b%mf+nMF
          nd = nd_(:,m); nf = nf_(:,m)
          nd_01x01(1,l,m) = distance(nd, noeud(:,l,m))/distance(nd, nf)
          if (l + b%l_charge(1,numproc)-1<b%ld_tot) nd_01x01(1,l,m) = -nd_01x01(1,l,m) ! mailles fictives a gauche
       end do
    end do
    deallocate(nd_, nf_)
!!! en m
    allocate(nd_(2,b%ld-(nMF+1):b%lf+nMF))
    allocate(nf_(2,b%ld-(nMF+1):b%lf+nMF))
    call find_noeud_bord(b, 3, noeud, b%ld, b%lf, nd_(:,:))
    call find_noeud_bord(b, 4, noeud, b%ld, b%lf, nf_(:,:))
    do l = b%ld-(nMF+1), b%lf+nMF
       nd = nd_(:,l); nf = nf_(:,l)
       do m = b%md-(nMF+1), b%mf+nMF
          nd_01x01(2,l,m) = distance(nd, noeud(:,l,m))/distance(nd, nf)
          if (m + b%m_charge(1,numproc)-1<b%md_tot) nd_01x01(2,l,m) = -nd_01x01(2,l,m)
       end do
    end do
    deallocate(nd_, nf_)

  end subroutine coordonees_maillage_01x01

  ! Calcul direct de Mpinv en passant par la pseudo-inverse
  subroutine calc_Mpinv_direct(b,Mpinv,stencil,ld,lf,md,mf)
    type (STR_BLOC), pointer :: b
    real*8, dimension(:,:,:), pointer :: Mpinv
    integer,intent(in) :: ld,lf,md,mf,stencil

    real*8, dimension(stencil-1,2)  :: A
    Real*8, dimension(2,stencil-1)  :: pseudoinv
    real*8, dimension(2*(stencil-1)):: d
    real*8 :: p

    integer :: l, m, i

    p = b%calcul_gradient_poids

    select case(stencil)
    case(5)
       do m = md, mf
          do l = ld, lf
             ! poids sur d pour retrouver la formule d'interpolztion exacte pour les fct lineaire
             d(:) = 0.d0

             A(1,1) = b%centre_cl(1,l-1,m)  -b%centre_cl(1,l,m)
             A(2,1) = b%centre_cl(1,l+1,m)  -b%centre_cl(1,l,m)
             A(3,1) = b%centre_cl(1,l,m-1)  -b%centre_cl(1,l,m)
             A(4,1) = b%centre_cl(1,l,m+1)  -b%centre_cl(1,l,m)

             A(1,2) = b%centre_cl(2,l-1,m)  -b%centre_cl(2,l,m)
             A(2,2) = b%centre_cl(2,l+1,m)  -b%centre_cl(2,l,m)
             A(3,2) = b%centre_cl(2,l,m-1)  -b%centre_cl(2,l,m)
             A(4,2) = b%centre_cl(2,l,m+1)  -b%centre_cl(2,l,m)

             !if poids (A -> WA)
             do i=1,4
                d(i) = sqrt(A(i,1)**2+A(i,2)**2)**(-p)
                A(i,:) = A(i,:)*d(i)
             end do

             pseudoinv = matmul(inverse_mat22(matmul(transpose(A),A)),transpose(A))

             !if poids
             do i=1,4
                pseudoinv(:,i) = pseudoinv(:,i)*d(i)
             end do

             Mpinv(1:4,l,m) = pseudoinv(1,:)
             Mpinv(5:8,l,m) = pseudoinv(2,:)

          end do
       end do

    case(9)
       do m = md, mf
          do l = ld, lf
             ! poids sur d pour retrouver la formule d'interpolztion exacte pour les fct lineaire
             d(:) = 0.d0

             A(1,1) = b%centre_cl(1,l-1,m)  -b%centre_cl(1,l,m)
             A(2,1) = b%centre_cl(1,l+1,m)  -b%centre_cl(1,l,m)
             A(3,1) = b%centre_cl(1,l,m-1)  -b%centre_cl(1,l,m)
             A(4,1) = b%centre_cl(1,l,m+1)  -b%centre_cl(1,l,m)
             A(5,1) = b%centre_cl(1,l-1,m-1)-b%centre_cl(1,l,m)
             A(6,1) = b%centre_cl(1,l+1,m-1)-b%centre_cl(1,l,m)
             A(7,1) = b%centre_cl(1,l-1,m+1)-b%centre_cl(1,l,m)
             A(8,1) = b%centre_cl(1,l+1,m+1)-b%centre_cl(1,l,m)

             A(1,2) = b%centre_cl(2,l-1,m)  -b%centre_cl(2,l,m)
             A(2,2) = b%centre_cl(2,l+1,m)  -b%centre_cl(2,l,m)
             A(3,2) = b%centre_cl(2,l,m-1)  -b%centre_cl(2,l,m)
             A(4,2) = b%centre_cl(2,l,m+1)  -b%centre_cl(2,l,m)
             A(5,2) = b%centre_cl(2,l-1,m-1)-b%centre_cl(2,l,m)
             A(6,2) = b%centre_cl(2,l+1,m-1)-b%centre_cl(2,l,m)
             A(7,2) = b%centre_cl(2,l-1,m+1)-b%centre_cl(2,l,m)
             A(8,2) = b%centre_cl(2,l+1,m+1)-b%centre_cl(2,l,m)

             !if poids (A -> WA)
             do i=1,8
                d(i) = sqrt(A(i,1)**2+A(i,2)**2)**(-p)
                A(i,:) = A(i,:)*d(i)
             end do

             pseudoinv = matmul(inverse_mat22(matmul(transpose(A),A)),transpose(A))

             !if poids
             do i=1,8
                pseudoinv(:,i) = pseudoinv(:,i)*d(i)
             end do

             Mpinv(1:8,l,m)  = pseudoinv(1,:)
             Mpinv(9:16,l,m) = pseudoinv(2,:)

          end do
       end do
    end select

  end subroutine calc_Mpinv_direct

  subroutine calc_Mpinv_from_QR(b,Mpinv,stencil,ld,lf,md,mf)
    type (STR_BLOC), pointer :: b
    real*8, dimension(:,:,:), pointer :: Mpinv
    integer,intent(in) :: ld,lf,md,mf,stencil

    real*8, dimension(stencil-1,2) :: A, Q
    Real*8, dimension(2,stencil-1) :: pseudoinv
    real*8, dimension(2,2)         :: R
    real*8, dimension(2*(stencil-1)):: d
    real*8 :: p

    integer :: l, m, i

    p = b%calcul_gradient_poids

    select case(stencil)
    case(5)
       do m = md, mf
          do l = ld, lf
             ! poids sur d pour retrouver la formule d'interpolztion exacte pour les fct lineaire
             d(:) = 0.d0

             A(1,1) = b%centre_cl(1,l-1,m)  -b%centre_cl(1,l,m)
             A(2,1) = b%centre_cl(1,l+1,m)  -b%centre_cl(1,l,m)
             A(3,1) = b%centre_cl(1,l,m-1)  -b%centre_cl(1,l,m)
             A(4,1) = b%centre_cl(1,l,m+1)  -b%centre_cl(1,l,m)

             A(1,2) = b%centre_cl(2,l-1,m)  -b%centre_cl(2,l,m)
             A(2,2) = b%centre_cl(2,l+1,m)  -b%centre_cl(2,l,m)
             A(3,2) = b%centre_cl(2,l,m-1)  -b%centre_cl(2,l,m)
             A(4,2) = b%centre_cl(2,l,m+1)  -b%centre_cl(2,l,m)

             !if poids (A -> WA)
             do i=1,4
                d(i) = sqrt(A(i,1)**2+A(i,2)**2)**(-p)
                A(i,:) = A(i,:)*d(i)
             end do

             call QR_decomposition(4,2,A,Q,R)

             pseudoinv = matmul(inverse_mat22(R),transpose(Q))

             !if poids
             do i=1,4
                pseudoinv(:,i) = pseudoinv(:,i)*d(i)
             end do

             Mpinv(1:4,l,m) = pseudoinv(1,:)
             Mpinv(5:8,l,m) = pseudoinv(2,:)

          end do
       end do

    case(9)
       do m = md, mf
          do l = ld, lf
             ! poids sur d pour retrouver la formule d'interpolztion exacte pour les fct lineaire
             d(:) = 0.d0

             A(1,1) = b%centre_cl(1,l-1,m)  -b%centre_cl(1,l,m)
             A(2,1) = b%centre_cl(1,l+1,m)  -b%centre_cl(1,l,m)
             A(3,1) = b%centre_cl(1,l,m-1)  -b%centre_cl(1,l,m)
             A(4,1) = b%centre_cl(1,l,m+1)  -b%centre_cl(1,l,m)
             A(5,1) = b%centre_cl(1,l-1,m-1)-b%centre_cl(1,l,m)
             A(6,1) = b%centre_cl(1,l+1,m-1)-b%centre_cl(1,l,m)
             A(7,1) = b%centre_cl(1,l-1,m+1)-b%centre_cl(1,l,m)
             A(8,1) = b%centre_cl(1,l+1,m+1)-b%centre_cl(1,l,m)

             A(1,2) = b%centre_cl(2,l-1,m)  -b%centre_cl(2,l,m)
             A(2,2) = b%centre_cl(2,l+1,m)  -b%centre_cl(2,l,m)
             A(3,2) = b%centre_cl(2,l,m-1)  -b%centre_cl(2,l,m)
             A(4,2) = b%centre_cl(2,l,m+1)  -b%centre_cl(2,l,m)
             A(5,2) = b%centre_cl(2,l-1,m-1)-b%centre_cl(2,l,m)
             A(6,2) = b%centre_cl(2,l+1,m-1)-b%centre_cl(2,l,m)
             A(7,2) = b%centre_cl(2,l-1,m+1)-b%centre_cl(2,l,m)
             A(8,2) = b%centre_cl(2,l+1,m+1)-b%centre_cl(2,l,m)

             !if poids (A -> WA)
             do i=1,8
                d(i) = sqrt(A(i,1)**2+A(i,2)**2)**(-p)
                A(i,:) = A(i,:)*d(i)
             end do

             call QR_decomposition(8,2,A,Q,R)

             pseudoinv = matmul(inverse_mat22(R),transpose(Q))

             !if poids
             do i=1,8
                pseudoinv(:,i) = pseudoinv(:,i)*d(i)
             end do

             Mpinv(1:8,l,m)  = pseudoinv(1,:)
             Mpinv(9:16,l,m) = pseudoinv(2,:)

          end do
       end do
    end select

  end subroutine calc_Mpinv_from_QR

end module metrique_maillage
