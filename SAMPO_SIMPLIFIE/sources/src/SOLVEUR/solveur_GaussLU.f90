module Gauss_LU
  
  implicit none

!!! Methode d elimination de Gauss
!!! GAUSSLU : factorisation de la matrice A
!!! GAUSSRES : resolution de Ax=b
!!! cf livre Analyse numerique matricielle appliquee a l'art de l'ingenieur; 1) methodes directes, pages 195-198

contains
  
  subroutine GAUSSLU(A, PIV)
    real*8 ,dimension(:,:),intent(inout):: A
    integer,dimension(:)  ,intent(out)  :: PIV
    
    integer :: i, j, k, n, JMAX
    real*8  :: AMAX, PIVINV, C
    
    n = size(PIV)
    do i = 1, n
       PIV(i) = i
    end do
    
    do k = 1, n-1
       AMAX = abs(A(PIV(k),k))
       JMAX = k

       do i = k+1, n
          if (abs(A(PIV(i),k)) > AMAX) then
             AMAX = abs(A(PIV(i),k))
             JMAX = i
          end if
       end do
       i = PIV(k)
       PIV(k) = PIV(JMAX)
       PIV(JMAX) = i
       if (A(PIV(k),k) /= 0.d0 ) then      !!!!!!!!!!!!!!!!!ATTENTION Modif
          PIVINV = 1.d0/A(PIV(k),k)
       else
          PIVINV=0.d0
       end if
       
       do i = k+1, n
          C = A(PIV(i),k) * PIVINV
          A(PIV(i),k) = C
          do j = k+1,n
             A(PIV(i),j) = A(PIV(i),j)-C*A(PIV(k),j)
          end do
       end do
    end do
  end subroutine GAUSSLU

  subroutine GAUSSRES(A, PIV, b, x)
    real*8 , dimension(:,:),intent(in)  :: A
    integer, dimension(:)  ,intent(in)  :: PIV
    real*8 , dimension(:)  ,intent(in)  :: b
    real*8 , dimension(:)  ,intent(out) :: x
    
    integer :: i, j, n
    real*8  :: C
    
    n = size(PIV)
    x(1) = b(PIV(1))
    do i = 2, n
       C = 0.d0
       do j = 1, i-1
          C = C + A(PIV(i),j)*x(j)
       end do
       x(i) = b(PIV(i))-C
    end do
    
    do i = n, 1, -1
       C = 0.d0
       do j = i+1, n
          C = C + A(PIV(i),j)*x(j)
       end do
       
       if (A(PIV(i),i) /= 0.d0) then     !!!!!!!!!!!!!!!!!ATTENTION modif
          x(i) = (x(i)-C)/A(PIV(i),i)
       elseif (A(PIV(i),i) /=0.d0 .and. x(i)-C == 0.d0) then
          x(i) = 0.d0
       else
          print*, "division par 0 GaussLU"
          x(i)=0.d0
          !call arret_code
       end if
    end do
  end subroutine GAUSSRES

end module Gauss_LU
