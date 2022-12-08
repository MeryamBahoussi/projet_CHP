module solveur_1d
  use m_MPI
  use m_struct

  implicit none

  public RL_JACOBI

  private initialisation
contains


!!!***************************************************************
!!!*                                                             *
!!!*    Inversion matricielle par blocs par méthode de Jacobi    *
!!!*                                                             *
!!!***************************************************************
  subroutine RL_JACOBI(uc, sm, matrice)

    type (matrice_tridiagonale)    , intent(inout)  :: matrice
    real*8        , dimension(:,-1:)      , intent(inout) :: uc
    real*8, dimension(:, :) , intent(inout) :: sm
    real*8, dimension(1:matrice%dim_bloc, 1:matrice%dim_mat)  :: sm_eq
    real*8, dimension(1:matrice%dim_bloc, 1:matrice%dim_mat)  :: d_Uc
    !Variables locales
    integer :: i, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call LU_decomposition(matrice)

    do l = 1 , matrice%dim_mat
       do i = 1, matrice%dim_bloc
          sm_eq(i, l) = sm(i, l)
       end do
    end do



    call LU_resolution(matrice, sm_eq) 



    d_Uc = sm_eq



    do l = 1, matrice%dim_mat
       do i = 1, matrice%dim_bloc
          uc(i, l) = uc(i, l) + d_Uc(i, l)
       end do
    end do

  end subroutine RL_JACOBI



  !***************************************************************
  !*                                                             *
  !*                         DECOMPOSITION LU                    *
  !*            on suppose : matrice = (a,b,c) = L * U           *
  !*            L est donnée par les couples (x_i,l_i)           *
  !*      avec x_i sous la diagonale et l_i le terme diagonal    *
  !*            U est donnée par les couples (y_i,u_i)           *
  !*      avec y_i sur la diagonale et u_i le terme diagonal     *
  !*      La décomposition est faite en supposant u_i = Id,      *
  !*      ce qui implique que :                                  *
  !*      x_i = a_i, l_i = b_i, x_i * yi, yi = l_i^-1 *c_i       *
  !*      En réalité, on stocke x_i dans aa_i, y_i dans cc_i,    *
  !*      et l_i^-1 dans bb_i                                    *
  !***************************************************************

  subroutine LU_decomposition(matrice)


    type (matrice_tridiagonale),   intent(inout)     :: matrice

    real*8 ,dimension(1:matrice%dim_bloc, 1:matrice%dim_bloc) :: bloc, blocm1
    integer                                                      :: l, i, j, k, dim_bloc, dim_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real*8,        pointer,    dimension(:,:,:)   :: aa, bb, cc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    aa => matrice%a; bb => matrice%b; cc => matrice%c

    dim_bloc = matrice%dim_bloc; dim_mat = matrice%dim_mat

    !   Initialisations (premiers termes du systeme):

    bloc(:,:) = bb(:,:,1)

    call matrice_inversion(bloc, blocm1, 0, 0, 1)
    !   cc(:,:,1) = matmul( blocm1(:,:), cc(:,:,1) )

    bloc = 0.d0

    do k = 1, dim_bloc
       do j = 1, dim_bloc
          do i = 1, dim_bloc
             bloc(i, j) = bloc(i, j) + blocm1(i, k) * cc(k, j, 1)
          end do
       end do
    end do

    cc(:,:,1) = bloc(:,:)
    bb(:,:,1) = blocm1(:,:)        

    !   Construction des matrices  L et U :
    do l = 2, dim_mat
       !   bb(:,:,l) = bb(:,:,l) - matmul( aa(:,:,l), cc(:,:,l-1) )

       bloc(:,:) = bb(:,:,l)

       do k = 1, dim_bloc
          do j = 1, dim_bloc
             do i = 1, dim_bloc
                bloc(i, j) = bloc(i, j) - aa(i, k, l) * cc(k, j, l-1)
             end do
          end do
       end do

       call matrice_inversion(bloc, blocm1 , 0, 0, 1)

       !   cc(:,:,l) = matmul( blocm1(:,:), cc(:,:,l) )

       bloc = 0.d0

       do k = 1, dim_bloc
          do j = 1, dim_bloc
             do i = 1, dim_bloc
                bloc(i, j) = bloc(i, j) + blocm1(i, k) * cc(k, j, l)
             end do
          end do
       end do


       cc(:,:,l) = bloc(:,:)
       bb(:,:,l) = blocm1(:,:)  

    end do

  end subroutine LU_decomposition

  !***************************************************************
  !*                                                             *
  !*      solveur du systeme tridiagonal PAR DECOMPOSITION LU    *
  !*            on suppose : matrice = (a,b,c) = L * U           *
  !*            L est donnée par les couples (x_i,l_i)           *
  !*      avec x_i sous la diagonale et l_i le terme diagonal    *
  !*            U est donnée par les couples (y_i,u_i)           *
  !*      avec y_i sur la diagonale et u_i le terme diagonal     *
  !*      La décomposition est faite en supposant u_i = Id,      *
  !*      ce qui implique que :                                  *
  !*      x_i = a_i, l_i = b_i, x_i * yi, yi = l_i^-1 *c_i       *
  !*      En réalité, on stocke x_i dans aa_i, y_i dans cc_i,    *
  !*      et l_i^-1 dans bb_i                                    *
  !***************************************************************


  subroutine LU_resolution(matrice, sm)

    type (matrice_tridiagonale), intent(inout) :: matrice
    real*8, dimension(matrice%dim_bloc, matrice%dim_mat), intent(inout) :: sm 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real*8, dimension(1:matrice%dim_bloc) :: sms, smp
    integer                               :: l, i, k
    integer                               :: dim_bloc, dim_mat
    real*8, pointer, dimension(:,:,:)     :: aa, bb, cc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    aa => matrice%a; bb => matrice%b; cc => matrice%c
    dim_bloc = matrice%dim_bloc; dim_mat = matrice%dim_mat


    if(ubound(sm,1) /= dim_bloc ) then
       print *,'subroutine LU_resolution : les blocs du second membre n"ont pas la même taille que ceux de la matrice'
       call arret_code
    end if

    !   Initialisation en l = 1
    !   sms(:) = matmul( bb(:,:,1), sm(:,1) )

    sms = 0.d0

    do k = 1, dim_bloc
       do i = 1, dim_bloc
          sms(i) = sms(i) + bb(i, k, 1) * sm(k, 1)
       end do
    end do

    sm(:, 1) = sms(:)

    !   Construction des vecteurs-solutions
    !   Descente :
    do l = 2, dim_mat

       sms(:) = sm(:,l)

       !   sms(:) = sms(:) - matmul( aa(:,:,l), sm(:,l-1) )
       do k = 1, dim_bloc
          do i = 1, dim_bloc
             sms(i) = sms(i) - aa(i, k, l) * sm(k, l-1)
          end do
       end do

       !   smp(:) = matmul( bb(:,:,l), sms(:) )

       smp = 0.d0

       do k = 1, dim_bloc
          do i = 1, dim_bloc
             smp(i) = smp(i) + bb(i, k, l) * sms(k)
          end do
       end do

       sm(:,l) = smp(:)

    end do

    !   RemontÃ©e :
    do l = dim_mat-1 , 1, -1

       sms(:) = sm(:,l)

       !   sm(:,l) = sm(:,l) - matmul( cc(:,:,l), sm(:,l+1) )
       do k = 1, dim_bloc
          do i = 1, dim_bloc
             sms(i) = sms(i) - cc(i, k, l) * sm(k, l+1)
          end do
       end do

       sm(:,l) = sms(:)


    end do

  end subroutine LU_resolution



  !***************************************************************
  !*                                                             *
  !*             Inversion d'une matrice                         *
  !*                                                             *
  !***************************************************************

  subroutine matrice_inversion(a1, am1, pivotage, nb_secmb, inversion)

    real*8        ,dimension(:, :)     , intent(inout)      :: a1     
    real*8        ,dimension(:, :)     , intent(out)        :: am1     
    integer                               , intent(in)         :: pivotage
    integer                               , intent(in)         :: nb_secmb, inversion

    real*8        ,dimension(:)        , allocatable        :: val     
    integer                               , save               :: nf, nf_col, k, i, j
    integer                               , save               :: pivot

    !   nb_secmb = nombre de second-membres stockes dans la matrice am1
    !   nf       = nombre de lignes et colonnes de la matrice a1
    !   nf_col   = nombre total de colonnes de la matrice am1

    !   nf_col = nb_secmb si on n'a pas besoin de l'inverse de a1
    !   am1 ne sert alors qu a stocker les solutions correspondant 
    !   aux second-membres donnes en entree (cas inversion = NON)

    nf = ubound(a1,1)
    nf_col = nf + nb_secmb

    if (inversion == 1) then
       nf_col = nf + nb_secmb
    else
       nf_col = nb_secmb
    end if

    allocate(val(nf + nb_secmb))

    if (inversion == 1) then
       do j =1, nf
          do i =1, nf
             am1(i, j) = 0.d0
          end do
          am1(j,j) = 1.d0
       end do
    end if

    !  Triangularisation sup (-> a1) et inf (-> am1) :
    do k =1, nf-1
       if (pivotage == 1) then
          pivot = k
          val(1) = abs(a1(k,k))
          do i = k+1, nf
             if (abs(a1(i,k)).gt.val(1)) then
                pivot = i
                val(1) = abs(a1(i,k))
             end if
          end do
          do j = k, nf
             val(j) = a1(k,j)
             a1(k,j) = a1(pivot,j)
          end do
          do j = k, nf
             a1(pivot,j) = val(j)
          end do
          do j = 1, nf_col
             val(j) = am1(k,j)
             am1(k,j) = am1(pivot,j)
          end do
          do j = 1, nf_col
             am1(pivot,j) = val(j)
          end do
       end if

       do j = k+1, nf 
          do i = k+1, nf  
             a1(i,j) = a1(i,j) - a1(i,k) * a1(k,j) / a1(k,k)  
          end do
       end do
       do i = k+1, nf  
          do j = 1, nf_col 
             am1(i,j) = am1(i,j) - a1(i,k) * am1(k,j) / a1(k,k)
          end do
       end do
    end do

    !   Remontée :

    do j = 1, nf_col
       am1(nf,j) = am1(nf,j) / a1(nf,nf) 
    end do

    do k = nf-1, 1, -1    
       do i = k+1, nf 
          do j = 1, nf_col    
             am1(k,j) = am1(k,j) - a1(k,i) * am1(i,j)
          end do
       end do
       do j = 1, nf_col    
          am1(k,j) = am1(k,j) / a1(k,k)
       end do
    end do
    deallocate(val)

  end subroutine matrice_inversion


  !***************************************************************
  !*                                                             *
  !*      solveur du systeme tridiagonal PAR DECOMPOSITION LU    *
  !*                                                             *
  !***************************************************************

  subroutine initialisation(matrice)

    type (matrice_tridiagonale),   intent(inout)     :: matrice
    real*8,        pointer,    dimension(:,:,:)   :: a, b, c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                 :: i, l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    a => matrice%a; b => matrice%b; c => matrice%c


    do l = 1, matrice%dim_mat
       a(:,:,l) = 0.d0
       b(:,:,l) = 0.d0
       c(:,:,l) = 0.d0
    end do

    do l = 1, matrice%dim_mat
       do i = 1, matrice%dim_bloc
          a(i,i,l) = -1.d0
          b(i,i,l) = 8.d0
          c(i,i,l) = -1.d0
       end do
    end do

  end subroutine initialisation



  !***************************************************************
  !*                                                             *
  !*             multiplication de 2 matrices                    *
  !*                                                             *
  !***************************************************************

  subroutine multiplication_matrices(a, b, c, dim_mat)

    real*8        ,dimension(:, :)     , intent(inout)      :: a, b, c   
    integer           , intent(in)                  :: dim_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                 :: i, j, k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    do j = 1, dim_mat
       do i = 1, dim_mat
          c(i, j)= 0.d0
          do k = 1, dim_mat
             c(i, j) = c(i, j) + a(i, k) * b(k, j)
          end do
       end do
    end do

  end subroutine multiplication_matrices

  !***************************************************************
  !*                                                             *
  !*    multiplication matrice tridiagonale par blocs vecteur    *
  !*                                                             *
  !***************************************************************

  subroutine multiplication_matrice_tridia_bloc_vecteur(matrice, x, y)

    type (matrice_tridiagonale),   intent(in)     :: matrice
    real*8,  dimension(matrice%dim_bloc,matrice%dim_mat), intent(in) :: x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real*8,  dimension(matrice%dim_bloc,matrice%dim_mat), intent(out) :: y

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                :: dim_bloc, dim_mat
    integer                                :: l
    real*8,  dimension(matrice%dim_bloc) :: y1, y2, y3
    real*8,        pointer,    dimension(:,:,:)   :: a, b, c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    a => matrice%a; b => matrice%b; c => matrice%c

    dim_bloc = matrice%dim_bloc; dim_mat = matrice%dim_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call A_fois_X(b(:,:,1), x(:,1), y2, dim_bloc) 
    call A_fois_X(c(:,:,1), x(:,2), y3, dim_bloc) 

    y(:,1) = y2 + y3

    do l = 2, dim_mat - 1

       call A_fois_X(a(:,:,l), x(:,l-1), y1, dim_bloc) 
       call A_fois_X(b(:,:,l), x(:,l), y2, dim_bloc) 
       call A_fois_X(c(:,:,l), x(:,l+1), y3, dim_bloc) 

       y(:,l) = y1 + y2 + y3

    end do

    call A_fois_X(a(:,:,dim_mat), x(:,dim_mat-1), y1, dim_bloc) 
    call A_fois_X(b(:,:,dim_mat), x(:,dim_mat), y2, dim_bloc) 

    y(:,dim_mat) = y1 + y2

  end subroutine multiplication_matrice_tridia_bloc_vecteur

  subroutine A_fois_X(A, x, y, n)

    real*8,  dimension(n,n), intent(in) :: a
    real*8,  dimension(n), intent(in)   :: x
    integer  ,  intent(in)                 :: n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real*8,  dimension(n), intent(out) :: y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                :: i, j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    do i = 1 , n
       y(i) = 0.d0
       do j = 1, n
          y(i) = y(i) + a(i,j) * x(j)
       end do
    end do

  end subroutine A_fois_X



  !***************************************************************
  !*                                                             *
  !*    Inversion matricielle scalaire par méthode de Jacobi    *
  !*                                                             *
  !***************************************************************
  subroutine RL_JACOBI_scalaire(uc, sm, matrice)

    type (matrice_tridiagonale)    , intent(inout)  :: matrice
    real*8        , dimension(-1:)      , intent(inout) :: uc
    real*8, dimension(:) , intent(inout) :: sm
    real*8, dimension(1:matrice%dim_mat)  :: sm_eq
    real*8, dimension(1:matrice%dim_mat)  :: d_Uc
    !Variables locales
    integer :: l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call LU_decomposition(matrice)

    do l = 1 , matrice%dim_mat
       sm_eq(l) = sm(l)
    end do



    call LU_resolution(matrice, sm_eq) 



    do l = 1 , matrice%dim_mat
       d_Uc(l) = sm_eq(l)
    end do



    do l = 1, matrice%dim_mat
       uc(l) = uc(l) + d_Uc(l)
    end do

  end subroutine RL_JACOBI_scalaire

end module solveur_1d
