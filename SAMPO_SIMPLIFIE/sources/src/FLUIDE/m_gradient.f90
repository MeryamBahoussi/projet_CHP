module m_reconstruction_gradient
  use m_struct
  use m_MPI
  use parametres_globaux

  implicit none

  private norme

  contains

  ! methode des moindres carrés à 5 pts
  subroutine reconstruction_grad_WLSQ5(b, tab, grad, l, m)
    type (STR_BLOC), pointer     :: b
    integer        , intent(in)  :: l, m
    real*8, dimension(l-1:l+1,m-1:m+1), intent(in)  :: tab
    real*8         , dimension(2),intent(out) :: grad

    ! variables locales
    real*8, dimension(2,4) :: Mpinv_loc
    real*8, dimension(4)   :: dz

    !A * gradz = dz -> pseudo inverse : gradz = (tAA)-1 * tA * dz

    Mpinv_loc(1,:) = b%Mpinv(1:4,l,m)
    Mpinv_loc(2,:) = b%Mpinv(5:8,l,m)

    dz(1) = tab(l-1,m) - tab(l,m)
    dz(2) = tab(l+1,m) - tab(l,m)
    dz(3) = tab(l,m-1) - tab(l,m)
    dz(4) = tab(l,m+1) - tab(l,m)

    grad = matmul(Mpinv_loc,dz)
  end subroutine reconstruction_grad_WLSQ5

  ! methode des moindres carrés à 9 pts
  subroutine reconstruction_grad_WLSQ9(b, tab, grad, l, m)
    type (STR_BLOC), pointer     :: b
    integer        , intent(in)  :: l, m
    real*8, dimension(l-1:l+1,m-1:m+1), intent(in)  :: tab
    real*8         , dimension(2),intent(out) :: grad

    ! variables locales
    real*8, dimension(2,8) :: Mpinv_loc
    real*8, dimension(8)   :: dz

    !A * gradz = dz -> pseudo inverse : gradz = (tAA)-1 * tA * dz

    Mpinv_loc(1,:) = b%Mpinv(1:8,l,m)
    Mpinv_loc(2,:) = b%Mpinv(9:16,l,m)

    dz(1) = tab(l-1,m)   - tab(l,m)
    dz(2) = tab(l+1,m)   - tab(l,m)
    dz(3) = tab(l,m-1)   - tab(l,m)
    dz(4) = tab(l,m+1)   - tab(l,m)
    dz(5) = tab(l-1,m-1) - tab(l,m)
    dz(6) = tab(l+1,m-1) - tab(l,m)
    dz(7) = tab(l-1,m+1) - tab(l,m)
    dz(8) = tab(l+1,m+1) - tab(l,m)

    grad = matmul(Mpinv_loc,dz)
  end subroutine reconstruction_grad_WLSQ9

  ! methode de green-gauss centré sur les cellules
  subroutine reconstruction_grad_GG(b, tab, grad, l, m)
    type (STR_BLOC), pointer     :: b
    real*8, dimension(l-1:l+1,m-1:m+1), intent(in)  :: tab
    real*8         , dimension(2),intent(out) :: grad
    integer        , intent(in)  :: l, m

    real*8, pointer :: A
    real*8, pointer :: nxd(:), nxg(:), nyd(:), nyg(:)
    real*8, pointer :: dlxd, dlxg, dlyd, dlyg

    nxd => b%n_dir_l(:,l,m)
    nxg => b%n_dir_l(:,l-1,m)
    nyd => b%n_dir_m(:,l,m)
    nyg => b%n_dir_m(:,l,m-1)

    dlxd => b%dl_l(l,m)
    dlxg => b%dl_l(l-1,m)
    dlyd => b%dl_m(l,m)
    dlyg => b%dl_m(l,m-1)

    A => b%aire(l,m)

    grad(:) = (tab(l,m)+tab(l+1,m))*nxd*dlxd - &
              (tab(l,m)+tab(l-1,m))*nxg*dlxg + &
              (tab(l,m)+tab(l,m+1))*nyd*dlyd - &
              (tab(l,m)+tab(l,m-1))*nyg*dlyg
    grad(:) = 0.5d0*grad(:)/A
  end subroutine reconstruction_grad_GG


  ! methode de green-gauss centré sur les sommets
  subroutine reconstruction_grad_GG_vertex(b, tab, grad, l, m)
    type (STR_BLOC), pointer     :: b
    real*8, dimension(b%ld-nMF:,b%md-nMF:), intent(in)  :: tab
    real*8         , dimension(2),intent(out) :: grad
    integer        , intent(in)  :: l, m

    real*8, pointer :: A
    real*8, pointer :: nxd(:), nxg(:), nyd(:), nyg(:)
    real*8, pointer :: dlxd, dlxg, dlyd, dlyg

    nxd => b%n_dir_l_dual(:,l,m)
    nxg => b%n_dir_l_dual(:,l-1,m)
    nyd => b%n_dir_m_dual(:,l,m)
    nyg => b%n_dir_m_dual(:,l,m-1)

    dlxd => b%dl_l_dual(l,m)
    dlxg => b%dl_l_dual(l-1,m)
    dlyd => b%dl_m_dual(l,m)
    dlyg => b%dl_m_dual(l,m-1)

    A => b%aire_dual(l,m)

    grad(:) = -(tab(l,m)    + tab(l+1,m)  )*nyg*dlyg &
              +(tab(l+1,m)  + tab(l+1,m+1))*nxd*dlxd &
              +(tab(l+1,m+1)+ tab(l,m+1)  )*nyd*dlyd &
              -(tab(l,m+1)  + tab(l,m)    )*nxg*dlxg
    grad(:) = 0.5d0*grad(:)/A
  end subroutine reconstruction_grad_GG_vertex


  ! gradients X sur maillage curviligne, pour le visqueux et la diffusion
  subroutine calc_grad_face_x(metrique_grad, grad, tab, l, m)
    real*8, dimension(4), intent(in) :: metrique_grad
    type( STR_gradient ),intent(out) :: grad
    integer             , intent(in) :: l, m
    real*8, dimension(l-1:l+1,m-1:m+1), intent(in)  :: tab

    real*8 :: X_x, Y_y, X_y ,Y_x
    type( STR_gradient ) :: grad_XY

    grad_XY%X = tab(l+1,m)- tab(l,m)
    grad_XY%Y = .25d0 * (tab(l+1,m+1) - tab(l+1,m-1) + tab(l,m+1) - tab(l,m-1))

    ! propriétées metriques
    X_x = metrique_grad(1)
    X_y = metrique_grad(2)
    Y_x = metrique_grad(3)
    Y_y = metrique_grad(4)

    ! calcul du gradient sur la face (l,m)
    grad%x = X_x * grad_XY%X + Y_x * grad_XY%Y
    grad%y = X_y * grad_XY%X + Y_y * grad_XY%Y
  end subroutine calc_grad_face_x

  ! gradients Y sur maillage curviligne, pour le visqueux et la diffusion
  subroutine calc_grad_face_y(metrique_grad, grad, tab, l, m)
    real*8, dimension(4), intent(in) :: metrique_grad
    type( STR_gradient ),intent(out) :: grad
    integer             , intent(in) :: l, m
    real*8, dimension(l-1:l+1,m-1:m+1), intent(in)  :: tab

    real*8 :: X_x, Y_y, X_y ,Y_x
    type( STR_gradient ) :: grad_XY

    grad_XY%X = .25d0 * (tab(l+1,m+1) - tab(l+1,m) + tab(l-1,m+1) - tab(l-1,m))
    grad_XY%Y = tab(l,m+1) - tab(l,m)

    ! propriétées metriques
    X_x = metrique_grad(1)
    X_y = metrique_grad(2)
    Y_x = metrique_grad(3)
    Y_y = metrique_grad(4)

    ! calcul du gradient sur la face (l,m)
    grad%x = X_x * grad_XY%X + Y_x * grad_XY%Y
    grad%y = X_y * grad_XY%X + Y_y * grad_XY%Y
  end subroutine calc_grad_face_y

  function norme(u)
    real*8 :: u(2), norme

    norme = sqrt(u(1)*u(1) + u(2)*u(2))
  end function norme

end module
