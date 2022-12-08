module outils_solveur

  implicit none 

contains

  !***************************************************************
  !*                                                             *
  !*             Inversion d'une matrice 2x2                     *
  !*                                                             *
  !***************************************************************
  function inverse_mat22(mat) result (inv)
    real*8, dimension(2,2) :: mat, inv
    real*8                 :: det

    det = mat(1,1)*mat(2,2)-mat(1,2)*mat(1,2)
    if ( abs(det) < 1.d-15 ) then
       print*, ' Erreur inverse_mat22 : matrice non inversible', det
       inv = 0.d0
    else
       inv(1,1) =  mat(2,2)
       inv(1,2) = -mat(1,2)
       inv(2,1) = -mat(2,1)
       inv(2,2) =  mat(1,1)
       inv(:,:) = inv(:,:)/det
    end if
  end function inverse_mat22

  !***************************************************************
  !*                                                             *
  !*             Inversion d'une matrice 3x3                     *
  !*                                                             *
  !***************************************************************
  function inverse_mat33(M) result(invM)
    real*8, dimension(3,3) :: M, invM

    real*8 :: det

    det =   M(1,1)*M(2,2)*M(3,3) &
         +  M(1,2)*M(2,3)*M(3,1) &
         +  M(2,1)*M(3,2)*M(1,3) &
         -  M(3,1)*M(2,2)*M(1,3) &
         -  M(2,1)*M(1,2)*M(3,3) &
         -  M(3,2)*M(2,3)*M(1,1)

    if ( abs(det) < 1.d-15 ) then
       print*, ' Erreur inverse_mat33 : matrice non inversible', det
       invM = 0.d0
    else
       invM(1,1) = M(2,2)*M(3,3) - M(3,2)*M(2,3)
       invM(1,2) = M(1,3)*M(3,2) - M(1,2)*M(3,3)
       invM(1,3) = M(1,2)*M(2,3) - M(1,3)*M(2,2)
       invM(2,1) = M(2,3)*M(3,1) - M(2,1)*M(3,3)
       invM(2,2) = M(1,1)*M(3,3) - M(1,3)*M(3,1)
       invM(2,3) = M(1,3)*M(2,1) - M(1,1)*M(2,3)
       invM(3,1) = M(2,1)*M(3,2) - M(2,2)*M(3,1)
       invM(3,2) = M(1,2)*M(3,1) - M(1,1)*M(3,2)
       invM(3,3) = M(1,1)*M(2,2) - M(1,2)*M(2,1)
       invM      = invM / det
    end if
  end function inverse_mat33

  !Gram schmidt modifié
  subroutine QR_decomposition(m,n,A,Q,R)
    integer               , intent(in)    :: m,n
    real*8, dimension(m,n), intent(inout) :: A
    real*8, dimension(m,n), intent(out)   :: Q
    real*8, dimension(n,n), intent(out)   :: R

    integer :: i, j, k
    real*8  :: s

    Q = 0.d0
    R = 0.d0

    do k=1,n
       s = 0.d0
       do j=1,m
          s = s + A(j,k)**2
       enddo
       R(k,k) = sqrt(s)
       do j=1,m
          Q(j,k) = A(j,k)/R(k,k)
       enddo
       do i=k+1,n
          s = 0.d0
          do j=1,m
             s = s + A(j,i)*Q(j,k)
          enddo
          R(k,i) = s
          do j=1,m
             A(j,i) = A(j,i) - R(k,i)*Q(j,k)
          enddo
       enddo
    enddo
  end subroutine QR_decomposition
end module outils_solveur
