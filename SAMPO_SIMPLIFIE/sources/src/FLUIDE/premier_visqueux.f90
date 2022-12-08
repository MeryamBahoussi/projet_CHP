module premier_visqueux
  use m_struct
  use a_mettre_en_donnees
  use m_MPI
  use outils_fluide, only : moyenne
  use parametres_fluide
  implicit none 
  
  private grad_visqueux_x, grad_visqueux_y, implicitation_M_visqueux

contains

  subroutine premier_membre_visqueux(sb, b, matrice)
    type (STR_SUPER_BLOC), pointer ::  sb
    type(STR_BLOC), pointer :: b
    type (STR_MATRICE), pointer :: matrice
    
    integer :: ld, lf, md, mf
    real*8, dimension(:,:,:,:), pointer :: aa, bb, cc, dd, ee
    real*8, dimension(:,:,:,:), pointer :: ad, ae, cd, ce
    
    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    
    aa => matrice%A; bb => matrice%B
    cc => matrice%C; dd => matrice%D
    ee => matrice%E 

    AD => b%mat%AD; CD => b%mat%CD
    AE => b%mat%AE; CE => b%mat%CE
!!$    ! on remet a zero les matrices des coins
!!$    AD(:,:,:,:) = 0.d0; CD(:,:,:,:) = 0.d0
!!$    AE(:,:,:,:) = 0.d0; CE(:,:,:,:) = 0.d0
    

!!!   matrices
!!!   D_lm X_l,m-1 + A_lm X_l-1,m + B_lm X_lm + C_lm X_l+1,m + E_lm X_l,m+1 = X^n   
!!! 
!!!   _______________________
!!!  |       |       |       |
!!!  | (A E) |   E   | (C E) |   ^ m
!!!  |_______|_______|_______|   |
!!!  |       |       |       |   |---> l
!!!  |   A   |   B   |   C   |
!!!  |_______|_______|_______|
!!!  |       |       |       |
!!!  | (A D) |   D   | (C D) |
!!!  |_______|_______|_______|
!!!   
   

!!! pour avoir la bonne implicitation
!!! on ne veut que la partie visqueuse a ces endroits
   ! call efface_mat_euler(b, aa, cc, dd, ee)
    call efface_mat_euler(b, matrice)
    
    call grad_visqueux_x(b, aa, bb, cc, dd, ee, AD, CD, AE, CE)

    call grad_visqueux_y(b, aa, bb, cc, dd, ee, AD, CD, AE, CE)

!!!! implicitation de la matrice
    call implicitation_M_visqueux(b,sb%dt,b%face(1),aa,bb,b%n_dir_l(:,ld-1,md:mf),md,mf)
    call implicitation_M_visqueux(b,sb%dt,b%face(2),cc,bb,b%n_dir_l(:,lf,md:mf),md,mf)
    call implicitation_M_visqueux(b,sb%dt,b%face(3),dd,bb,b%n_dir_m(:,ld:lf,md-1),ld,lf)
    call implicitation_M_visqueux(b,sb%dt,b%face(4),ee,bb,b%n_dir_m(:,ld:lf,mf),ld,lf) 
    
  end subroutine premier_membre_visqueux
 
  ! Pour l'implicitation de la partie visqueuse, 
  ! on ne veut que la contribution de la matrice visqeuse dans la matrice totale
  ! on efface donc la partie euler aux endroits ou on doit impliciter
  subroutine efface_mat_euler(b, matrice)
    type (STR_BLOC)   , pointer :: b
    type (STR_MATRICE), pointer :: matrice
    
    integer                   :: ipatch
    type (STR_PATCH), pointer :: patch
    type(STR_FACE)  , pointer :: face

    real*8, dimension(:,:,:,:), pointer :: AA, BB, CC, DD, EE
   real*8, dimension(:,:,:,:), pointer :: ad, ae, cd, ce

    AE => matrice%AE; EE => matrice%E; CE => matrice%CE
    AA => matrice%A ; BB => matrice%B; CC => matrice%C 
    AD => matrice%AD; DD => matrice%D; CD => matrice%CD
    
    ! face 1 : matrice A
    face => b%face(1)
    do ipatch = 1, face%nb_patch
       patch => face%patch(ipatch)
       
       select case (trim(patch%bc%condition))
       case("interproc", 'interbloc')
          ! on ne met pas a 0
          ! on veut garder la matrice complete pour assembler avec PETSc
       case default
          AE(:,:,b%ld,patch%ideb:patch%ifin)=0.d0
          AA(:,:,b%ld,patch%ideb:patch%ifin)=0.d0
          AD(:,:,b%ld,patch%ideb:patch%ifin)=0.d0
       end select
    end do
    
    ! face 2 : matrice C
    face => b%face(2)
    do ipatch = 1, face%nb_patch
       patch => face%patch(ipatch)
       
       select case (trim(patch%bc%condition))
       case("interproc", 'interbloc')
       case default
          CE(:,:,b%lf,patch%ideb:patch%ifin)=0.d0
          CC(:,:,b%lf,patch%ideb:patch%ifin)=0.d0
          CD(:,:,b%lf,patch%ideb:patch%ifin)=0.d0
       end select
    end do
    
    ! face 3 : matrice D
    face => b%face(3)
    do ipatch = 1, face%nb_patch
       patch => face%patch(ipatch)
       
       select case (trim(patch%bc%condition))
       case("interproc", 'interbloc')
       case default
          AD(:,:,patch%ideb:patch%ifin, b%md)=0.d0
          DD(:,:,patch%ideb:patch%ifin, b%md)=0.d0
          CD(:,:,patch%ideb:patch%ifin, b%md)=0.d0
       end select
    end do

    ! face 4 : matrice E
    face => b%face(4)
    do ipatch = 1, face%nb_patch
       patch => face%patch(ipatch)
       
       select case (trim(patch%bc%condition))
       case("interproc", 'interbloc')
       case default
          AE(:,:,patch%ideb:patch%ifin, b%mf)=0.d0
          EE(:,:,patch%ideb:patch%ifin, b%mf)=0.d0
          CE(:,:,patch%ideb:patch%ifin, b%mf)=0.d0
       end select
    end do

  end subroutine efface_mat_euler
  
  subroutine grad_visqueux_x(bloc, A, B, C, D, E, AD, CD, AE, CE)
    
    type (STR_BLOC), pointer :: bloc
    real*8, dimension(3,3,bloc%ld:bloc%lf,bloc%md:bloc%mf), intent(inout) :: A,B,C,D,E
    real*8, dimension(3,3,bloc%ld:bloc%lf,bloc%md:bloc%mf), intent(inout) :: AD,CD,AE,CE

    integer :: l, m, lg, ld, mg, md
    real*8 :: mu_moy, Jacobien_moy, Jac_g, Jac_d
    real*8 :: grad_XY_X(3,2), grad_XY_Y(3,2)
    real*8 :: d_x(3,2), d_y(3,2)
    real*8 :: X_x, X_y, Y_x, Y_y
    real*8 :: flux(3,3,3,2)
    real*8, dimension(:,:,:), pointer :: n_l, n_m
    real*8, dimension(:,:), pointer :: dl_l, dl_m
    
    n_l => bloc%n_dir_l
    n_m => bloc%n_dir_m
    dl_l => bloc%dl_l
    dl_m => bloc%dl_m

    ! on initialise a zero car il y a des composantes qui ne sont jamais calculees
    flux = 0.d0
    grad_XY_X = 0.d0; grad_XY_Y = 0.d0  
    
    do m = bloc%md, bloc%mf
       do l = bloc%ld-1, bloc%lf
          !   ___________________
          !  |  l,m+1  | l+1,m+1 |
          !  |   3,1   |   3,2   |
          !  |_________|_________|
          !  |   l,m   |  l+1,m  |
          !  |   2,1   x   2,2   |
          !  |_________|_________|
          !  |  l,m-1  | l+1,m-1 |
          !  |   1,1   |   1,2   |
          !  |_________|_________|
          !

          lg = l  ; mg = m  ! cell_g = cell(l,m)
          ld = l+1; md = m  ! cell_d = cell(l+1,m)

!!! Grandeurs moyennes a l'interface
          mu_moy = moyenne(bloc%mu(lg,mg), bloc%mu(ld,md))

!!! Jacobiens
          Jac_g = 1.d0 / bloc%aire(lg,mg)
          Jac_d = 1.d0 / bloc%aire(ld,md)
          Jacobien_moy = 2.d0*Jac_g*Jac_d/(Jac_g + Jac_d)

!!! Gradient a l'interface dans le maillage cartesien unitaire (X,Y)
          ! grad_._XY%X =  .(ld,md) - .(lg,mg)
          grad_XY_X(2,1) =-1.d0; grad_XY_X(2,2) = 1.d0
          
          ! grad_._XY%Y = .25d0 * (.(l+1,m+1) - .(l+1,m-1) + .(l,m+1) - .(l,m-1))
          grad_XY_Y(3,1) = .25d0; grad_XY_Y(3,2) = .25d0
          grad_XY_Y(1,1) =-.25d0; grad_XY_Y(1,2) =-.25d0
          
!!! Calcul des proprietes metriques 
          X_x = n_l(1,l,m)*dl_l(l,m)
          X_y = n_l(2,l,m)*dl_l(l,m)  
          Y_x = .25d0 * ( n_m(1,l,m)*dl_m(l,m) + n_m(1,l+1,m)*dl_m(l+1,m) + &
               n_m(1,l,m-1)*dl_m(l,m-1) + n_m(1,l+1,m-1)*dl_m(l+1,m-1) ) 
          Y_y = .25d0 * ( n_m(2,l,m)*dl_m(l,m) + n_m(2,l+1,m)*dl_m(l+1,m) + &
               n_m(2,l,m-1)*dl_m(l,m-1) + n_m(2,l+1,m-1)*dl_m(l+1,m-1) ) 

!!! Gradient a l'interface sur le maillage curviligne (repere (x,y))
          !._x(:,:) = X_x * grad_._XY_X(:,:) + Y_x * grad_._XY_Y(:,:)
          !._y(:,:) = X_y * grad_._XY_X(:,:) + Y_y * grad_._XY_Y(:,:)
          d_x(:,:) = X_x * grad_XY_X(:,:) + Y_x * grad_XY_Y(:,:)
          d_y(:,:) = X_y * grad_XY_X(:,:) + Y_y * grad_XY_Y(:,:)
         
          ! tenseur des contraintes
          !tau_xx = mu_moy * ( 4.d0/3.d0 * u_x - 2.d0/3.d0 * v_y )
          !tau_xy = mu_moy * ( u_y + v_x )
          !tau_yy = mu_moy * ( 4.d0/3.d0 * v_y - 2.d0/3.d0 * u_x )
          !ev(iXu) = tau_xx ; fv(iXu) = tau_xy  
          !ev(iXv) = tau_xy ; fv(iXv) = tau_yy
          !flux(:,l,m) = Jacobien_moy*dl_l(l,m)*( n_l(1,l,m)*ev(:) + n_l(2,l,m)*fv(:) )

!!! flux(x,y,l,m) : Pour l'equation sur x, composante relative a y au point (l,m)  
          flux(iXu,iXu,:,:) = n_l(1,l,m) *  4.d0/3.d0 * d_x(:,:) + n_l(2,l,m) * d_y(:,:)
          flux(iXu,iXv,:,:) = n_l(1,l,m) *(-2.d0/3.d0)* d_y(:,:) + n_l(2,l,m) * d_x(:,:)
          flux(iXv,iXu,:,:) = n_l(1,l,m) * d_y(:,:) + n_l(2,l,m) *(-2.d0/3.d0)* d_x(:,:)
          flux(iXv,iXv,:,:) = n_l(1,l,m) * d_x(:,:) + n_l(2,l,m) *  4.d0/3.d0 * d_y(:,:) 
          !flux(iXp, :,:,:) = 0.d0
          flux(:,:,:,:) = Jacobien_moy * mu_moy * dl_l(l,m) * flux(:,:,:,:) 

!!! Assemblage des matrices
          ! sm(1,l,m) = sm(1,l,m) - (flux(l,m)-flux(l-1,m))
          if ( l /= bloc%ld-1 ) then ! sm(1,l,m) = sm(1,l,m) - flux(l,m)
             D( :,:,l,m) = D( :,:,l,m) - flux(:,:,1,1)
             CD(:,:,l,m) = CD(:,:,l,m) - flux(:,:,1,2) !! coins
             B( :,:,l,m) = B( :,:,l,m) - flux(:,:,2,1)
             C( :,:,l,m) = C( :,:,l,m) - flux(:,:,2,2)
             E( :,:,l,m) = E( :,:,l,m) - flux(:,:,3,1)
             CE(:,:,l,m) = CE(:,:,l,m) - flux(:,:,3,2) !!
          end if

          if ( l /= bloc%lf ) then ! sm(1,l+1,m) = sm(1,l+1,m) + flux(l,m)
             AD(:,:,l+1,m) = AD(:,:,l+1,m) + flux(:,:,1,1) !!
             D( :,:,l+1,m) = D( :,:,l+1,m) + flux(:,:,1,2)
             A( :,:,l+1,m) = A( :,:,l+1,m) + flux(:,:,2,1)
             B( :,:,l+1,m) = B( :,:,l+1,m) + flux(:,:,2,2)
             AE(:,:,l+1,m) = AE(:,:,l+1,m) + flux(:,:,3,1) !!
             E( :,:,l+1,m) = E( :,:,l+1,m) + flux(:,:,3,2)
          end if
       end do
    end do
 
  end subroutine grad_visqueux_x

  subroutine grad_visqueux_y(bloc, A, B, C, D, E, AD, CD, AE, CE) 
    
    type (STR_BLOC), pointer :: bloc
    real*8, dimension(3,3,bloc%ld:bloc%lf,bloc%md:bloc%mf), intent(inout) :: A,B,C,D,E
    real*8, dimension(3,3,bloc%ld:bloc%lf,bloc%md:bloc%mf), intent(inout) :: AD,CD,AE,CE

    integer :: l, m, lg, ld, mg, md
    real*8 :: mu_moy, Jacobien_moy, Jac_g, Jac_d
    real*8 :: grad_XY_X(2,3), grad_XY_Y(2,3)
    real*8 :: d_x(2,3), d_y(2,3)
    real*8 :: X_x, X_y, Y_x, Y_y
    real*8 :: flux(3,3,2,3)
    real*8, dimension(:,:,:), pointer :: n_l, n_m
    real*8, dimension(:,:), pointer :: dl_l, dl_m
    
    n_l => bloc%n_dir_l
    n_m => bloc%n_dir_m
    dl_l => bloc%dl_l
    dl_m => bloc%dl_m
  
    ! on initialise a zero car il y a des composantes qui ne sont jamais calculees
    flux = 0.d0
    grad_XY_X = 0.d0; grad_XY_Y = 0.d0  
    
    do m = bloc%md-1, bloc%mf
       do l = bloc%ld, bloc%lf
          !   _____________________________
          !  | l-1,m+1 |  l,m+1  | l+1,m+1 |
          !  |   2,1   |   2,2   |   2,3   |
          !  |_________|____x____|_________|
          !  |  l-1,m  |   l,m   |  l+1,m  |
          !  |   1,1   |   1,2   |   1,3   |
          !  |_________|_________|_________|
          ! 
          
          lg = l  ; mg = m    ! cell_g = cell(l,m)
          ld = l  ; md = m+1  ! cell_d = cell(l,m+1)

!!! Grandeurs moyennes a l'interface
          mu_moy = moyenne(bloc%mu(lg,mg), bloc%mu(ld,md))

!!! Jacobiens
          Jac_g = 1.d0 / bloc%aire(lg,mg)
          Jac_d = 1.d0 / bloc%aire(ld,md)
          Jacobien_moy = 2.d0*Jac_g*Jac_d/(Jac_g + Jac_d)

!!! Gradient a l'interface dans le maillage cartesien unitaire (X,Y)
          !grad_._XY%Y = .(ld,md) - .(lg,mg)
          grad_XY_Y(1,2) =-1.d0; grad_XY_Y(2,2) = 1.d0
         
          !grad_._XY%X = .25d0 * ( .(l+1,m+1) + .(l+1,m) - .(l-1,m+1) - .(l-1,m) )
          grad_XY_X(1,3) = .25d0; grad_XY_X(2,3) = .25d0
          grad_XY_X(1,1) =-.25d0; grad_XY_X(2,1) =-.25d0
          

!!! Calculs des proprietes metriques 
          Y_x = n_m(1,l,m)*dl_m(l,m)  
          Y_y = n_m(2,l,m)*dl_m(l,m)  
          X_x = .25d0 * ( n_l(1,l,m)*dl_l(l,m) + n_l(1,l,m+1)*dl_l(l,m+1) + &
               n_l(1,l-1,m)*dl_l(l-1,m) + n_l(1,l-1,m+1)*dl_l(l-1,m+1) ) 
          X_y = .25d0 * ( n_l(2,l,m)*dl_l(l,m) + n_l(2,l,m+1)*dl_l(l,m+1) + &
               n_l(2,l-1,m)*dl_l(l-1,m) + n_l(2,l-1,m+1)*dl_l(l-1,m+1) )
          
!!! Gradient a l'interface sur le maillage curviligne (repere (x,y))
          !._x = X_x * grad_._XY%X  + Y_x * grad_._XY%Y
          !._y = X_y * grad_._XY%X  + Y_y * grad_._XY%Y
          d_x(:,:) = X_x * grad_XY_X(:,:) + Y_x * grad_XY_Y(:,:)
          d_y(:,:) = X_y * grad_XY_X(:,:) + Y_y * grad_XY_Y(:,:)
             
          ! tenseur des contraintes
          ! tau_xx = mu_moy * ( 4.d0/3.d0 * u_x  - 2.d0/3.d0 * v_y )
          ! tau_yy = mu_moy * ( 4.d0/3.d0 * v_y  - 2.d0/3.d0 * u_x )
          ! tau_xy = mu_moy * ( u_y + v_x )
          ! ev(iXu) = tau_xx ; fv(iXu) = tau_xy
          ! ev(iXv) = tau_xy ; fv(iXv) = tau_yy 
          ! flux(:,l,m) = Jacobien_moy*dl_m(l,m)*( n_m(1,l,m)*ev(:) + n_m(2,l,m)*fv(:) )
          
!!! flux(x,y,l,m) : Pour l'equation sur x, composante relative a y au point (l,m)  
          flux(iXu,iXu,:,:) = n_m(1,l,m) *  4.d0/3.d0 *d_x(:,:) + n_m(2,l,m) * d_y(:,:)
          flux(iXu,iXv,:,:) = n_m(1,l,m) *(-2.d0/3.d0)*d_y(:,:) + n_m(2,l,m) * d_x(:,:)
          flux(iXv,iXu,:,:) = n_m(1,l,m) * d_y(:,:) + n_m(2,l,m) *(-2.d0/3.d0)*d_x(:,:)
          flux(iXv,iXv,:,:) = n_m(1,l,m) * d_x(:,:) + n_m(2,l,m) *  4.d0/3.d0*d_y(:,:)
          !flux(iXp, :,:,:) = 0.d0
          flux(:,:,:,:) = Jacobien_moy * mu_moy * dl_m(l,m) * flux(:,:,:,:) 
          
          
!!! Assemblage des matrices
          ! sm(1,l,m) = sm(1,l,m) - (flux(l,m)-flux(l-1,m))
          if ( m /= bloc%md-1 ) then ! sm(1,l,m) = sm(1,l,m) - flux(l,m)
             A( :,:,l,m) = A( :,:,l,m) - flux(:,:,1,1)
             B( :,:,l,m) = B( :,:,l,m) - flux(:,:,1,2)
             C( :,:,l,m) = C( :,:,l,m) - flux(:,:,1,3)
             AE(:,:,l,m) = AE(:,:,l,m) - flux(:,:,2,1) !! coins
             E( :,:,l,m) = E( :,:,l,m) - flux(:,:,2,2)
             CE(:,:,l,m) = CE(:,:,l,m) - flux(:,:,2,3) !!
          end if

          if ( m /= bloc%mf ) then ! sm(1,l+1,m) = sm(1,l+1,m) + flux(l,m)
             AD(:,:,l,m+1) = AD(:,:,l,m+1) + flux(:,:,1,1) !!
             D( :,:,l,m+1) = D( :,:,l,m+1) + flux(:,:,1,2)
             CD(:,:,l,m+1) = CD(:,:,l,m+1) + flux(:,:,1,3) !!
             A( :,:,l,m+1) = A( :,:,l,m+1) + flux(:,:,2,1)
             B( :,:,l,m+1) = B( :,:,l,m+1) + flux(:,:,2,2)
             C( :,:,l,m+1) = C( :,:,l,m+1) + flux(:,:,2,3)
          end if
          
       end do
    end do
    
  end subroutine grad_visqueux_y

  subroutine implicitation_M_visqueux(b,dt,face,aa,bb,n,taille_inf,taille_sup)
    use m_MPI
    use outils_conditions_limites

    type (STR_BLOC), intent(in) :: b
    real*8, intent(in) :: dt
    type(STR_FACE), intent(in) :: face
    integer, intent(in) :: taille_inf,taille_sup
    
    real*8, dimension(1:3,1:3, b%ld:b%lf, b%md:b%mf), intent(in)    :: aa ! matrice extradiagonale
    real*8, dimension(1:3,1:3, b%ld:b%lf, b%md:b%mf), intent(inout) :: bb ! matrice diagonale
    real*8, intent(in) :: n(2,taille_inf:taille_sup)

    real*8 :: normale(2)
    real*8, dimension(1:2,1:2) :: matrice_sym
    real*8, dimension(1:3,1:3) :: mat
    integer :: id, if, jd, jf
    integer :: l, m
    integer :: h_l,h_m
    integer :: ipatch
    type (STR_PATCH), pointer :: patch
    real*8 :: a, fac

    do ipatch = 1, face%nb_patch
       patch => face%patch(ipatch)

       select case (patch%bc%position)
       case(4) ! North
          id = patch%ideb; jd = b%mf
          if = patch%ifin; jf = b%mf
          h_l = 1; h_m = 0
       case(3) ! South
          id = patch%ideb; jd = b%md
          if = patch%ifin; jf = b%md
          h_l = 1; h_m = 0     
       case(2) ! East
          id = b%lf; jd = patch%ideb
          if = b%lf; jf = patch%ifin
          h_l = 0; h_m = 1     
       case(1) ! West
          id = b%ld; jd = patch%ideb
          if = b%ld; jf = patch%ifin  
          h_l = 0; h_m = 1       
       end select

       select case (trim(patch%bc%condition))
       case("interproc", 'interbloc')
          ! on fait la resolution avec MKL sur le proc 0 uniquement
          ! inversion en parallele avec Petsc 
       case("entree_supersonique")
          ! tout est au second membre 
          ! ou on ne fait rien si on resoud en delta forme
       case("symetrie")
          do m = jd, jf
             do l = id, if
                matrice_sym = 0.d0
                normale(1:2) = n(1:2,l*h_l+m*h_m)
                
                matrice_sym(1,1) = 1.0d0-2.0d0*normale(1)**2
                matrice_sym(2,1) = -2.0d0*normale(1)*normale(2)
                matrice_sym(1,2) = matrice_sym(2,1)
                matrice_sym(2,2) = 1.0d0-2.0d0*normale(2)**2
                
                bb(iXu:iXv,iXu:iXv,l,m) = bb(iXu:iXv,iXu:iXv,l,m) + &
                     matmul(aa(iXu:iXv,iXu:iXv,l,m),matrice_sym(1:2,1:2))
             end do
          end do
       case("paroi_global", "paroi_global_adiabatique", "paroi_flux_impose")
          !, "injection") 
          ! u_MF = - u_int
          ! v_MF = - v_int
          ! rien sur la pression dans la partie visqueuse
          do m = jd, jf
             do l = id, if
                matrice_sym = 0.d0
                matrice_sym(1,1) = -1.d0
                matrice_sym(2,2) = -1.d0
                
                bb(iXu:iXv,iXu:iXv,l,m) = bb(iXu:iXv,iXu:iXv,l,m) + &
                     matmul(aa(iXu:iXv,iXu:iXv,l,m),matrice_sym(1:2,1:2))
                
             end do
          end do
          
       case("injection") 
          ! tout est au second membre

       case("entree_subsonique", "entree_subsonique_t")
          if (patch%bc%position /= 1) then
             print*, "la condition d'entree subsonique n'est implicitée que pour l=ld"
             call arret_code
          end if
          ! X_exterieur = Mat X_interieur
          do m = jd, jf
             do l = id, if
                a = b%rho(0,m)*b%c(0, m)
                fac = a*dt/b%dm(0,m) 

                call mat_pression_entree_subsonique(a, fac, mat)
                
                bb(:,:,l,m) = bb(:,:,l,m) + matmul(aa(:,:,l,m),mat(:,:))
             end do
          end do
          
       case("flux_nul")
          do m = jd, jf
             do l = id, if
                bb(:,iXu:iXv,l,m) = bb(:,iXu:iXv,l,m) + aa(:,iXu:iXv,l,m)
             end do
          end do
          
       case default
          print*, "Condition inconnue pour l'implicitation de la matrice visqueuse ", patch%bc%condition
          call arret_code
       end select
    end do
  end subroutine implicitation_M_visqueux

end module premier_visqueux
