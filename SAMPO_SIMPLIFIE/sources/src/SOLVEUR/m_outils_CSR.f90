  include 'mkl_dss.f90'

module m_outils_CSR
  use m_MPI
  use mkl_dss
  implicit none 

contains
  !________________________________________________________________
  !                                                              
  !!                          
  ! 	SUBROUTINE de remplissage des tableaux Columns et rowIndex du stockage CSR avec les données issue	
  !   du tableaux des connectivités (graphe) et de celui du nombre-de-non-zero par ligne	
  !     			
  !   nb_rows : nombre de ligne de la matrice	
  !   connect_max = 9 pour des inconnues scalaires ou 18 pour des vecteurs à 2 composantes
  !________________________________________________________________
  subroutine structure_CSR(graphe,Nombre_non_zero,columns, rowIndex,nNonZeros, nb_rows, connect_max)
    implicit none
    integer, intent(in) :: nNonZeros, nb_rows,  connect_max
    integer, intent(out) :: columns(1:nNonZeros), rowIndex(1:nb_rows+1) ! Struture CSR
    integer, intent(in) :: graphe(1:connect_max, 1:nb_rows), Nombre_non_zero(1:nb_rows) ! graphe de connectivités
    integer :: k

    rowIndex(1) = 1
    do k = 1, nb_rows
       rowIndex(k+1) = rowIndex(k) + Nombre_non_zero(k)

       columns(rowIndex(k):rowIndex(k+1)-1) = graphe(1:Nombre_non_zero(k),k)
    end do

    if (rowIndex(nb_rows+1)/=nNonZeros+1) then
       print*, "ERREUR :INCOHERENCE ENTRE NNONZEROS ET ROWINDEX"
       call arret_code
    end if
  end subroutine structure_CSR

  subroutine Structure_graphe_connectivites(Imax,Jmax, graphe, Nombre_non_zero, nb_rows,dim_Uc,periodique)
    implicit none 
    integer, intent(in) :: Imax,Jmax
    integer, intent(in) :: nb_rows,dim_Uc
    integer, intent(inout)  :: graphe(9*dim_Uc,1:nb_rows), Nombre_non_zero(1:nb_rows)
    logical, intent(in) :: periodique
    integer :: k_pt,k_elt, ligne_pt
    integer :: i,j,r


    ! dim_Uc donne la taille des vecteurs solutions par cellule : si dim_Uc = 1 il s'agit d'un scalaire

    ! la matrice sera donnée par des tableaux correspondant à chaque diagonale {f,d,g,a,b,c,h,e,i} 
    ! ( il s'agit de diagonales composées par des blocs dim_Uc*dim_Uc )
    ! 
    ! Connectivitées (schéma à 9 points) pour le point ij avec la position des termes des  diagonales
    ! ____________________
    !|      |      |      |
    !|  h   |  e   |  i   |
    !|      |      |      |
    !|------|------|------|
    !|  a   |  b   |  c   |
    !|      | Pt_ij|      |UNIT_DEBUG_MATRICE
    !|------|------|------|
    !|  f   |  d   |  g   |
    !| 	    |      |      |
    !|______|______|______|

    Nombre_non_zero = 0 

    !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
    !  Dans un premier temps on rempli la premiere ligne pour chaque point ( indice ligne_pt ) avec les connectivités 
    !  Il s'agit de remplir le tableau des connectivités avec la premiere ligne de chaque bloc {f,d,g,a,b,c,h,e,i} pour le pt ij
    !  Ensuite, on reparcourt le graphe pour dupliquer cette premiere ligne de chaque pt sur les (dim_Uc-1) lignes suivantes 
    !  qui appartiennent au meme point
    !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

    ! Termes de f  
    do j = 2, Jmax  
       do i = 2, Imax

          k_pt = i + (j-1) * Imax
          ligne_pt  = dim_Uc * (k_pt - 1) + 1

          do r = 1, dim_Uc
             k_elt = dim_Uc * ( (k_pt - Imax - 1)  - 1) + r
             Nombre_non_zero( ligne_pt ) = Nombre_non_zero( ligne_pt ) + 1
             graphe( Nombre_non_zero( ligne_pt ), ligne_pt ) = k_elt
          end do

       end do
    end do

    ! Termes de d  
    do j = 2, Jmax  
       do i = 1, Imax

          k_pt = i + (j-1) * Imax
          ligne_pt  = dim_Uc * (k_pt - 1) + 1

          do r = 1, dim_Uc
             k_elt = dim_Uc * ( (k_pt - Imax)  - 1) + r
             Nombre_non_zero( ligne_pt ) = Nombre_non_zero( ligne_pt ) + 1
             graphe( Nombre_non_zero( ligne_pt ), ligne_pt ) = k_elt
          end do

       end do
    end do

    ! Termes de g  
    do j = 2, Jmax  
       do i = 1, Imax-1

          k_pt = i + (j-1) * Imax
          ligne_pt  = dim_Uc * (k_pt - 1) + 1

          do r = 1, dim_Uc
             k_elt = dim_Uc * ( (k_pt - Imax + 1)  - 1) + r
             Nombre_non_zero( ligne_pt ) = Nombre_non_zero( ligne_pt ) + 1
             graphe( Nombre_non_zero( ligne_pt ), ligne_pt ) = k_elt
          end do

       end do
    end do

    ! Termes de a  
    do j = 1, Jmax
       do i = 2, Imax

          k_pt = i + (j-1) * Imax
          ligne_pt  = dim_Uc * (k_pt - 1) + 1

          do r = 1, dim_Uc
             k_elt = dim_Uc * ( (k_pt - 1)  - 1) + r
             Nombre_non_zero( ligne_pt ) = Nombre_non_zero( ligne_pt ) + 1
             graphe( Nombre_non_zero( ligne_pt ), ligne_pt ) = k_elt
          end do

       end do
    end do

    ! Termes de b  
    do j = 1, Jmax
       do i = 1, Imax

          k_pt = i + (j-1) * Imax
          ligne_pt  = dim_Uc * (k_pt - 1) + 1

          do r = 1, dim_Uc
             k_elt = dim_Uc * ( (k_pt)  - 1) + r
             Nombre_non_zero( ligne_pt ) = Nombre_non_zero( ligne_pt ) + 1
             graphe( Nombre_non_zero( ligne_pt ), ligne_pt ) = k_elt
          end do

       end do
    end do

    ! Termes de c  
    do j = 1, Jmax
       do i = 1, Imax-1

          k_pt = i + (j-1) * Imax
          ligne_pt  = dim_Uc * (k_pt - 1) + 1

          do r = 1, dim_Uc
             k_elt = dim_Uc * ( (k_pt + 1)  - 1) + r
             Nombre_non_zero( ligne_pt ) = Nombre_non_zero( ligne_pt ) + 1
             graphe( Nombre_non_zero( ligne_pt ), ligne_pt ) = k_elt
          end do

       end do
    end do

    ! Termes de h  
    do j = 1, Jmax-1
       do i = 2, Imax

          k_pt = i + (j-1) * Imax
          ligne_pt  = dim_Uc * (k_pt - 1) + 1

          do r = 1, dim_Uc
             k_elt = dim_Uc * ( (k_pt + Imax - 1)  - 1) + r
             Nombre_non_zero( ligne_pt ) = Nombre_non_zero( ligne_pt ) + 1
             graphe( Nombre_non_zero( ligne_pt ), ligne_pt ) = k_elt
          end do

       end do
    end do

    ! Termes de e  
    do j = 1, Jmax-1
       do i = 1, Imax

          k_pt = i + (j-1) * Imax
          ligne_pt  = dim_Uc * (k_pt - 1) + 1

          do r = 1, dim_Uc
             k_elt = dim_Uc * ( (k_pt + Imax)  - 1) + r
             Nombre_non_zero( ligne_pt ) = Nombre_non_zero( ligne_pt ) + 1
             graphe( Nombre_non_zero( ligne_pt ), ligne_pt ) = k_elt
          end do
       end do
    end do

    ! Termes de i  
    do j = 1, Jmax-1
       do i = 1, Imax -1

          k_pt = i + (j-1) * Imax
          ligne_pt  = dim_Uc * (k_pt - 1) + 1

          do r = 1, dim_Uc
             k_elt = dim_Uc * ( (k_pt + Imax + 1)  - 1) + r
             Nombre_non_zero( ligne_pt ) = Nombre_non_zero( ligne_pt ) + 1
             graphe( Nombre_non_zero( ligne_pt ), ligne_pt ) = k_elt
          end do

       end do
    end do

    !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤ 
    ! Applications des conditions limites : seule la condition "periodique" impose une 
    ! modification de la matrice
    !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

    ! ¤¤¤¤¤ East
    if ( periodique ) then

       ! Termes de g
       do j= 2, Jmax
          k_pt = i + (j-1) * Imax
          ligne_pt  = dim_Uc * (k_pt - 1) + 1
          do r = 1, dim_Uc
             k_elt = dim_Uc * ( (k_pt - (2*Imax-1))  - 1) + r
             Nombre_non_zero( ligne_pt ) = Nombre_non_zero( ligne_pt ) + 1
             graphe( Nombre_non_zero( ligne_pt ), ligne_pt ) = k_elt
          end do
       end do

       ! Termes de c 
       do j= 1, Jmax
          k_pt = i + (j-1) * Imax
          ligne_pt  = dim_Uc * (k_pt - 1) + 1
          do r = 1, dim_Uc
             k_elt = dim_Uc * ( (k_pt - Imax - 1)  - 1) + r
             Nombre_non_zero( ligne_pt ) = Nombre_non_zero( ligne_pt ) + 1
             graphe( Nombre_non_zero( ligne_pt ), ligne_pt ) = k_elt
          end do
       end do

       ! Termes de i
       do j= 1, Jmax-1
          k_pt = i + (j-1) * Imax
          ligne_pt  = dim_Uc * (k_pt - 1) + 1
          do r = 1, dim_Uc
             k_elt = dim_Uc * ( (k_pt + 1)  - 1) + r
             Nombre_non_zero( ligne_pt ) = Nombre_non_zero( ligne_pt ) + 1
             graphe( Nombre_non_zero( ligne_pt ), ligne_pt ) = k_elt
          end do
       end do

    end if


    ! ¤¤¤¤¤ west : 
    if ( periodique ) then

       ! Termes de f 
       do j= 2, Jmax
          k_pt = i + (j-1) * Imax
          ligne_pt  = dim_Uc * (k_pt - 1) + 1
          do r = 1, dim_Uc
             k_elt = dim_Uc * ( (k_pt - 1)  - 1) + r
             Nombre_non_zero( ligne_pt ) = Nombre_non_zero( ligne_pt ) + 1
             graphe( Nombre_non_zero( ligne_pt ), ligne_pt ) = k_elt
          end do
       end do

       ! Termes de a 
       do j= 1, Jmax
          k_pt = i + (j-1) * Imax
          ligne_pt  = dim_Uc * (k_pt - Imax) + 1
          do r = 1, dim_Uc
             k_elt = dim_Uc * ( (k_pt + Imax - 1)  - 1) + r
             Nombre_non_zero( ligne_pt ) = Nombre_non_zero( ligne_pt ) + 1
             graphe( Nombre_non_zero( ligne_pt ), ligne_pt ) = k_elt
          end do
       end do

       ! Termes de h 
       do j= 1, Jmax
          k_pt = i + (j-1) * Imax
          ligne_pt  = dim_Uc * (k_pt - 1) + 1
          do r = 1, dim_Uc
             k_elt = dim_Uc * ( (k_pt + 2 * Imax - 1)  - 1) + r
             Nombre_non_zero( ligne_pt ) = Nombre_non_zero( ligne_pt ) + 1
             graphe( Nombre_non_zero( ligne_pt ), ligne_pt ) = k_elt
          end do
       end do

    end if


    !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤ 
    ! Recopie des lignes car les blocs  de chaque diagonales sont des matrices  dim_Uc*dim_Uc
    !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

    if (dim_Uc > 1) then 
       do j = 1, Jmax  
          do i = 1, Imax

             k_pt = i + (j-1) * Imax
             ligne_pt  = dim_Uc * (k_pt - 1) + 1

             do r = 2, dim_Uc
                k_elt = dim_Uc * ( (k_pt )  - 1) + r
                Nombre_non_zero( k_elt ) = Nombre_non_zero( ligne_pt ) 
                graphe( :, k_elt ) = graphe( :, ligne_pt ) 

             end do

          end do
       end do
    end if

  end subroutine Structure_graphe_connectivites


  subroutine Remplissage_matrice_CSR(Imax,Jmax, values, nNonZeros,  &
       & Mat_a, Mat_b, Mat_c, Mat_d, Mat_e, Mat_f, Mat_g, Mat_h, Mat_i, &
       & dim_Uc,periodique)
    implicit none 

    ! LE STOCKAGE DE LA MATRICE EST EN CSR : on stock ligne par ligne 
    ! rowIndex(*) : Nombre d'élément non nul pour chaque ligne
    ! columns(*)  : indice de colonnes des valeurs non nuls de chaque ligne par ordre croissant
    ! values (*)  : valeurs des coefficients de la matrice correspondant au indice de colonne.

    ! Connectivitées (schéma à 9 points) pour le point ij avec la position des termes des  diagonales
    ! ____________________
    !|      |      |      |
    !|  h   |  e   |  i   |
    !|      |      |      |
    !|------|------|------|
    !|  a   |  b   |  c   |
    !|      | Pt_ij|      |UNIT_DEBUG_MATRICE
    !|------|------|------|
    !|  f   |  d   |  g   |
    !| 	    |      |      |
    !|______|______|______|	    

    ! Variables     
    integer, intent(in) :: Imax,Jmax
    integer, intent(in) :: dim_Uc
    real*8, dimension(1:,1:,1:,1:), intent(in):: Mat_a, Mat_b, Mat_c, Mat_d, Mat_e, Mat_f, Mat_g, Mat_h, Mat_i
    logical, intent(in) :: periodique

    integer, intent(in) :: nNonZeros ! Nombre de non zeros
    real*8, dimension(:), pointer, intent(out) :: values   ! coefficients de la matrice

    integer :: i,j
    integer :: r, nnz_ind

    ! Remplissage de la matrice ligne par ligne
    values = 0.d0
    nnz_ind = 1

    do j = 1, Jmax
       do i = 1, Imax

          do r = 1, dim_Uc

             if ( i == Imax .and. j /= 1 .and. periodique ) then 
                values(nnz_ind : nnz_ind + dim_Uc -1) = Mat_g( r, 1:dim_Uc, i, j)
                nnz_ind = nnz_ind + dim_Uc
             end if

             if ( i /= 1 .and. j /= 1 )  then 
                values(nnz_ind : nnz_ind + dim_Uc -1) = Mat_f( r, 1:dim_Uc, i, j)
                nnz_ind = nnz_ind + dim_Uc
             end if

             if ( j /= 1 )  then 
                values(nnz_ind : nnz_ind + dim_Uc -1) = Mat_d( r, 1:dim_Uc, i, j)
                nnz_ind = nnz_ind + dim_Uc
             end if

             if ( j /= 1 .and. i /= Imax )  then 
                values(nnz_ind : nnz_ind + dim_Uc -1) = Mat_g( r, 1:dim_Uc, i, j)
                nnz_ind = nnz_ind + dim_Uc 
             end if

             if ( i == 1 .and. j /= 1 .and. periodique ) then 
                values(nnz_ind : nnz_ind + dim_Uc -1) = Mat_f( r, 1:dim_Uc, i, j)
                nnz_ind = nnz_ind + dim_Uc
             end if

             if ( i == Imax .and. periodique ) then 
                values(nnz_ind : nnz_ind + dim_Uc -1) = Mat_c( r, 1:dim_Uc, i, j)
                nnz_ind = nnz_ind + dim_Uc
             end if

             if ( i /= 1 )  then 
                values(nnz_ind : nnz_ind + dim_Uc -1) = Mat_a( r, 1:dim_Uc, i, j)
                nnz_ind = nnz_ind + dim_Uc
             end if

             values(nnz_ind : nnz_ind + dim_Uc -1) = Mat_b( r, 1:dim_Uc, i, j)
             nnz_ind = nnz_ind + dim_Uc

             if ( i /= Imax )  then 
                values(nnz_ind : nnz_ind + dim_Uc -1) = Mat_c( r, 1:dim_Uc, i, j)
                nnz_ind = nnz_ind + dim_Uc
             end if

             if ( i == 1 .and. periodique ) then 
                values(nnz_ind : nnz_ind + dim_Uc -1) = Mat_a( r, 1:dim_Uc, i, j)
                nnz_ind = nnz_ind + dim_Uc
             end if

             if ( i == Imax .and. j /= Jmax .and. periodique ) then 
                values(nnz_ind : nnz_ind + dim_Uc -1) = Mat_i( r, 1:dim_Uc, i, j)
                nnz_ind = nnz_ind + dim_Uc
             end if

             if ( i /= 1 .and. j /= Jmax )  then
                values(nnz_ind : nnz_ind + dim_Uc -1) = Mat_h( r, 1:dim_Uc, i, j)
                nnz_ind = nnz_ind + dim_Uc
             end if

             if ( j /= Jmax )  then 
                values(nnz_ind : nnz_ind + dim_Uc -1) = Mat_e( r, 1:dim_Uc, i, j)
                nnz_ind = nnz_ind + dim_Uc
             end if

             if ( i /= Imax .and. j /= Jmax )  then 
                values(nnz_ind : nnz_ind + dim_Uc -1) = Mat_i( r, 1:dim_Uc, i, j)
                nnz_ind = nnz_ind + dim_Uc
             end if

             if ( i == 1 .and. j /= Jmax .and. periodique ) then 
                values(nnz_ind : nnz_ind + dim_Uc -1) = Mat_f( r, 1:dim_Uc, i, j)
                nnz_ind = nnz_ind + dim_Uc
             end if

          end do
       end do
    end do


    if ( nnz_ind-1 /= nNonZeros ) then 
       print*,"PROBLEME DE REMPLISSAGE DE LA MATRICE CSR DANS 'solveur_MKL_for_CSR'", nnz_ind-1, nNonZeros
       print*,"ERREUR :Remplissage_matrice_CSR"
       call arret_code
    end if

    ! call verification_values_CSR(Imax, Jmax, values, nNonZeros,dim_Uc,rowIndex)

  end subroutine Remplissage_matrice_CSR

  subroutine solveur_MKL_for_CSR(Imax,Jmax,rowIndex, columns, values, Sm, nNonZeros, nRows, dim_Uc)
    implicit none

    ! LE STOCKAGE DE LA MATRICE EST EN CSR : on stocke ligne par ligne 
    ! rowIndex(*) : Nombre d'élément non nul pour chaque ligne
    ! columns(*)  : indice de colonnes des valeurs non nuls de chaque ligne par ordre croissant
    ! values (*)  : valeurs des coefficients de la matrice correspondant au indice de colonne.

    ! Variables   
    integer, intent(in) :: Imax,Jmax  
    integer, intent(in) :: dim_Uc
    integer, intent(in) :: nRows     ! Nombre de lignes
    integer :: nCols     ! Nombre de colonnes 
    integer, intent(in) :: nNonZeros ! Nombre de non zeros
    real*8, dimension(1:,1:,1:), intent(inout) :: Sm
    integer :: nRhs    ! Nombre de seconds membres
    real*8, dimension(:), pointer :: rhs      ! second membre
    integer, dimension(:), pointer   :: columns  ! indices de colonne des non zeros 
    integer, dimension(:), pointer   :: rowIndex ! indices de debut de ligne dans columns(:)
    real*8, dimension(:), pointer :: values   ! coefficients de la matrice
    real*8, dimension(:), pointer :: solution ! vecteur solution de la matrice
#ifdef MKL_ACTIF
    TYPE(MKL_DSS_HANDLE) :: handle
    integer :: error ! code d'erreur de retour des fonctions du solveur DSS
    integer :: perm(1),i,j

    integer :: k_pt, r, k_elt

    perm(1)=0
    nCols = nRows
    nRhs = 1

    ! Allocation de la memoire  
    ALLOCATE(solution(nRows))

    !Allocation et remplissage du second menbre 
    allocate(rhs(1:nNonZeros))
    rhs = 0.d0
    do j = 1, Jmax
       do i = 1, Imax
          k_pt = i + (j-1) * Imax
          do r = 1, dim_Uc
             k_elt = dim_Uc * ( (k_pt )  - 1) + r
             rhs(k_elt) = Sm(r,i,j)
          end do
       end do
    end do


    ! Initialisation du solveur
    error = dss_create(handle, MKL_DSS_DEFAULTS)

    ! Definition de la structure des non-zeros de la matrice
    error = dss_define_structure(handle, MKL_DSS_SYMMETRIC_STRUCTURE, rowIndex, nRows, nCols, columns, nNonZeros)
    !error = dss_define_structure(handle, MKL_DSS_NON_SYMMETRIC, rowIndex, nRows, nCols, columns, nNonZeros)

    ! Renumerotation de la matrice
    error = dss_reorder(handle, MKL_DSS_DEFAULTS, perm(1))

    ! Factorisation de la matrice
    error = dss_factor_real(handle, MKL_DSS_DEFAULTS, values)

    ! Resolution du probleme
    error = dss_solve_real(handle, MKL_DSS_DEFAULTS, rhs, nRhs, solution)

    ! Recopie du vecteur solution 
    do j = 1, Jmax
       do i = 1, Imax
          k_pt = i + (j-1) * Imax
          do r = 1, dim_Uc
             k_elt = dim_Uc * ( (k_pt )  - 1) + r
             Sm(r,i,j) = solution(k_elt)
          end do
       end do
    end do

    ! Liberation de la memoire
    error = dss_delete(handle, MKL_DSS_DEFAULTS)     
    DEALLOCATE(solution,rhs)
#endif
  end subroutine solveur_MKL_for_CSR

!!$  subroutine evaluation_conditionnement_matrice_CSR(Imax,Jmax,rowIndex, columns, values, nNonZeros, nRows, dim_Uc, MG, MP, cond)
!!$    implicit none
!!$
!!$    ! LE STOCKAGE DE LA MATRICE EST EN CSR : on stock ligne par ligne 
!!$    ! rowIndex(*) : Nombre d'élément non nul pour chaque ligne
!!$    ! columns(*)  : indice de colonnes des valeurs non nuls de chaque ligne par ordre croissant
!!$    ! values (*)  : valeurs des coefficients de la matrice correspondant au indice de colonne.
!!$
!!$    ! Variables     
!!$    integer, intent(in) :: Imax,Jmax
!!$    integer, intent(in) :: dim_Uc
!!$    integer, intent(in) :: nRows     ! Nombre de lignes
!!$
!!$    integer, intent(in) :: nNonZeros ! Nombre de non zeros
!!$    real*8, intent(out) :: MG, MP, cond
!!$    real*8, dimension(1:nRows) ::  v_0, v_1, v
!!$    real*8, dimension(1:dim_Uc,1:Imax,1:Jmax) :: u_0,u_1, u
!!$    integer, dimension(:), pointer   :: columns  ! indices de colonne des non zeros 
!!$    integer, dimension(:), pointer   :: rowIndex ! indices de debut de ligne dans columns(:)
!!$    real*8, dimension(:), pointer :: values   ! coefficients de la matrice
!!$
!!$    integer :: max_iter, i, iter, j
!!$    real*8 :: eps, n_0, n_1, MG_0, MG_1, MP_0, MP_1
!!$    integer :: r
!!$
!!$    max_iter = ( nRows )!**3 !2000
!!$    eps      = 1.d-6
!!$    MG_0 = 0.d0
!!$    MP_0 = 0.d0
!!$
!!$    !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
!!$    !¤¤¤¤¤¤  METHODE D'ITERATION DIREC ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
!!$    !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
!!$    ! La 
!!$
!!$
!!$    !$==========================================
!!$    !$... CALCUL DE LA PLUS GRANDE VALEUR PROPRE
!!$    !$==========================================
!!$    !$... On se donne un vecteur aleatoire initial. que l'on normalise
!!$    call random_number( v )
!!$
!!$    v_0 = v 
!!$
!!$    ! Calcul de la norme de u_0
!!$    n_0 = sqrt( dot_product(v_0,v_0) )
!!$
!!$    ! Normalisation du vecteur v_0
!!$    v_0 = v_0 / n_0
!!$
!!$    !$... On initialise le processus itératif
!!$    ! v_1 = matrice * v_0
!!$    v_1 = 0.d0
!!$    do i = 1, nRows
!!$       do r = rowIndex(i), rowIndex(i+1)-1
!!$          v_1(i) = v_1(i) + values( r ) * v_0( columns(r) )
!!$       end do
!!$    end do
!!$
!!$    ! Calcul de la norme de v_1
!!$    n_1 = sqrt( dot_product(v_1,v_1) )
!!$
!!$    MG_1 = n_1 / n_0
!!$
!!$    iter = 0
!!$
!!$    do while ( abs(MG_0-MG_1)>eps .and. iter<=max_iter ) 
!!$
!!$       !$... Le suivant devient le precedent :
!!$       MG_0 = MG_1
!!$       v_0  = v_1
!!$       n_0  = n_1
!!$
!!$       ! v_1 = matrice * v_0
!!$       v_1 = 0.d0
!!$       do i = 1, nRows
!!$          do r = rowIndex(i), rowIndex(i+1)-1
!!$             v_1(i) = v_1(i) + values( r ) * v_0( columns(r) )
!!$          end do
!!$       end do
!!$
!!$       MG_1 = abs( dot_product( v_0, v_1 ) ) / n_0**2
!!$
!!$       v_1 = v_1 / n_0
!!$
!!$       n_1 = sqrt( dot_product(v_1,v_1) )
!!$
!!$       iter = iter + 1
!!$
!!$    end do
!!$
!!$
!!$    !$==========================================
!!$    !$... CALCUL DE LA PLUS PETITE VALEUR PROPRE
!!$    !$==========================================
!!$    !$... On fait un processus identique : La plus petite valeur propre de A
!!$    !$... est la plus grande valeur propre de Inverse(A)
!!$    !$... On se donne un vecteur aleatoire initial. que l'on normalise
!!$
!!$    call random_number( u )
!!$
!!$    u_0 = u 
!!$
!!$    ! Calcul de la norme de u_0
!!$    n_0 = 0.d0
!!$    do j = 1, Jmax 
!!$       do i = 1, Imax 
!!$          do r = 1, dim_Uc
!!$             n_0 = n_0 + ( u_0(r,i,j)**2 )
!!$          end do
!!$       end do
!!$    end do
!!$    n_0 = sqrt( n_0 )
!!$
!!$    u_0 = u_0 / n_0
!!$
!!$    !$... On initialise le processus itératif
!!$    u_1 = u_0
!!$    call solveur_MKL_for_CSR(Imax,Jmax,rowIndex, columns, values, u_1, nNonZeros, nRows , dim_Uc)
!!$
!!$
!!$    ! Calcul de la norme de u_1
!!$    n_1 = 0.d0
!!$    do j = 1, Jmax 
!!$       do i = 1, Imax 
!!$          do r = 1, dim_Uc
!!$             n_1 = n_1 + ( u_1(r,i,j)**2 )
!!$          end do
!!$       end do
!!$    end do
!!$    n_1 = sqrt( n_1 )
!!$
!!$    MP_1 = n_1 / n_0
!!$
!!$    iter = 0
!!$
!!$    do while ( abs(MP_0-MP_1)>eps .and. iter<=max_iter ) 
!!$
!!$       !$... Le suivant devient le precedent :
!!$       MP_0 = MP_1
!!$       u_0  = u_1
!!$       n_0  = n_1
!!$
!!$       call solveur_MKL_for_CSR(Imax,Jmax,rowIndex, columns, values, u_1, nNonZeros, nRows , dim_Uc)
!!$
!!$       MP_1 = 0.d0
!!$       do j = 1, Jmax 
!!$          do i = 1, Imax 
!!$             do r = 1, dim_Uc
!!$                MP_1 = MP_1 + ( u_1(r,i,j) *  u_0(r,i,j) )
!!$             end do
!!$          end do
!!$       end do
!!$       MP_1 = abs( MP_1 ) / n_0**2
!!$
!!$       u_1 = u_1 / n_0
!!$
!!$       ! Calcul de la norme de u_1
!!$       n_1 = 0.d0
!!$       do j = 1, Jmax 
!!$          do i = 1, Imax 
!!$             do r = 1, dim_Uc
!!$                n_1 = n_1 + ( u_1(r,i,j)**2 )
!!$             end do
!!$          end do
!!$       end do
!!$       n_1 = sqrt( n_1 )
!!$
!!$       iter = iter + 1
!!$
!!$    end do
!!$
!!$    !$...
!!$    !$... Sortie des résultats dans MG, MP, cond
!!$    !$...
!!$    MG   = MG_1
!!$    MP   = 1.d0 / MP_1
!!$    !cond = sqrt( MG / MP ) 
!!$    cond =  MG / MP 
!!$
!!$  end subroutine evaluation_conditionnement_matrice_CSR


end module m_outils_CSR
