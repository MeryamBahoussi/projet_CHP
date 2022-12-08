module SOLVEUR_balayage
  use solveur_1d
  use m_struct
  implicit none
  
contains
  
  subroutine RL_par_balayage(b, X_imp)
    use m_MPI
    implicit none
    
    type (STR_BLOC), intent(inout) ::  b
    real*8, dimension(3,b%ld-2:b%lf+2,b%md-2:b%mf+2), intent(inout) :: X_imp 
!!$    integer :: i, l, m
!!$    integer :: ld, lf, md, mf
!!$    integer :: k, kmax=40
!!$    integer :: STATINFO
!!$    real*8 :: residu, residu_MPI, tol=1.d-8
!!$    real*8, dimension(3) :: erreur, norme_X0 
!!$    real*8, dimension(:,:,:), pointer :: V
!!$    real*8, dimension(:,:,:), allocatable :: X_0, X_k
!!$    real*8, dimension(3,b%md:b%mf) :: Sm,DX,EX
!!$    
!!$    real*8, dimension(:,:,:,:), pointer :: aa, bb, cc, dd, ee
!!$    aa => b%mat_a; bb => b%mat_b
!!$    cc => b%mat_c; dd => b%mat_d
!!$    ee => b%mat_e
!!$    
    call arret_code

!!$       ! implicitation de la matrice
!!$       call implicitation_matrice(b,b%face(3),aa,bb,b%n_dir_m(:,ld:lf,md-1),ld,lf)
!!$       call implicitation_matrice(b,b%face(4),cc,bb,b%n_dir_m(:,ld:lf,mf),ld,lf)
!!$       call implicitation_cond_lim_mat_subsonique(b,bb(:,:,ld,:),dd(:,:,ld,:))




!!$    allocate( X_0(3,b%ld:b%lf, b%md:b%mf) )
!!$    allocate( X_k(3,b%ld-2:b%lf+2, b%md-2:b%mf+2) )
!!$    X_0(:,ld:lf,md:mf) = X_imp(:,ld:lf,md:mf)
!!$    X_k(:,:,:) = X_imp(:,:,:)
!!$
!!$!!! norme de X0 pour le calcul du residu
!!$    norme_X0 = 0.d0
!!$    do m = md, mf
!!$       do l = ld, lf
!!$          do i = 1, 3
!!$             norme_X0(i) = norme_X0(i) + X_0(i,l,m)**2.d0
!!$          end do
!!$       end do
!!$    end do
!!$    norme_X0(:) = sqrt(max(norme_X0(:), 1.d-12))
!!$
!!$    k = 0; residu=1.d70; residu_MPI=1.d70;
!!$    do while(k<kmax .and. residu_MPI>tol)
!!$       
!!$!!! conditions aux limites, remplissage des mailles fictives   
!!$       call MF_pour_implicite(b,X_k,2,4)   ! mailles fictives pour les faces 2 a 4
!!$       call MF_pour_implicite(b,X_imp,1,1) ! mailles fictives pour la face 1
!!$
!!$       do l = ld, lf
!!$          b%matrice%a = aa(:,:,l,:)
!!$          b%matrice%b = bb(:,:,l,:)
!!$          b%matrice%c = cc(:,:,l,:)
!!$
!!$!!! second membre
!!$          call mult_mat_vect(b,dd(:,:,l,:),X_imp(:,l-1,md:mf),DX)
!!$          call mult_mat_vect(b,ee(:,:,l,:),X_k(:,l+1,md:mf),EX)
!!$          Sm = X_0(:,l,md:mf) - DX - EX
!!$
!!$          ! modification du sm en fonction des conditions limites
!!$          call implicitation_cond_lim_sm(b, X_k(:,l,md-2:mf+2), Sm)
!!$          
!!$!!! resolution
!!$          X_imp(:,l,:) = 0.d0
!!$          call RL_JACOBI(X_imp(:,l,:), Sm, b%matrice)
!!$       end do
!!$
!!$!!! calcul du residu
!!$       erreur = 0.d0
!!$       do m = md, mf
!!$          do l = ld, lf
!!$             do i = 1, 3
!!$                erreur(i) = erreur(i) + (X_imp(i,l,m)-X_k(i,l,m))**2.d0
!!$             end do
!!$          end do
!!$       end do
!!$       do i = 1, 3
!!$          erreur(i) = sqrt(erreur(i)) / norme_X0(i)
!!$       end do
!!$       residu = maxval(erreur)
!!$       call MPI_ALLREDUCE( residu, residu_MPI, 1, MPI_REAL8, &
!!$            MPI_MAX, MPI_COMM_WORLD, STATINFO) 
!!$
!!$       !if (numproc==0) print*, k, residu
!!$
!!$!!! preparation de l iteration suivante
!!$       k = k+1
!!$       X_k(:,ld:lf,md:mf) = X_imp(:,ld:lf,md:mf)
!!$    end do
!!$
!!$    if (numproc==0 .and. k>=kmax .and. residu_MPI>1.d-6)  then 
!!$       print*, "Inversion non convergee : ", residu_MPI
!!$    end if
    
  end subroutine RL_par_balayage

!!$  
!!$  subroutine mult_mat_vect(b,A,x,y)
!!$    
!!$    type (STR_BLOC), intent(in) ::  b
!!$    real*8, dimension(3,3,b%md:b%mf), intent(in) :: A
!!$    real*8, dimension(1:3,b%md:b%mf), intent(in) :: x
!!$    real*8, dimension(1:3,b%md:b%mf), intent(out) :: y
!!$    integer :: m
!!$
!!$    do m = b%md, b%mf
!!$       call A_fois_X(a(:,:,m),x(:,m),y(:,m),3)
!!$    end do
!!$  end subroutine mult_mat_vect
!!$
!!$  subroutine implicitation_cond_lim_mat_subsonique(b,bb,dd)
!!$    use m_temps
!!$    use m_MPI
!!$    implicit none
!!$
!!$    type (STR_BLOC), intent(in) ::  b
!!$    real*8, dimension(3,3,b%md:b%mf), intent(in) :: dd
!!$    real*8, dimension(3,3,b%md:b%mf), intent(out) :: bb
!!$
!!$    type (STR_PATCH), pointer :: patch
!!$    integer :: ip
!!$    integer :: m  
!!$    real*8 :: fac, a
!!$
!!$    do ip = 1, b%face(1)%nb_patch
!!$       patch => b%face(1)%patch(ip)
!!$
!!$       select case (trim(patch%bc%condition))
!!$       case("entree_subsonique")
!!$
!!$          do m = patch%ideb, patch%ifin
!!$             a = b%rho(0,m)*b%c(0, m)
!!$             fac = a*dt/b%dm(0,m) 
!!$
!!$             bb(:,1,m) = bb(:,1,m) - a*fac/(1.d0+fac)*dd(:,3,m)
!!$             bb(:,3,m) = bb(:,3,m) + fac/(1.d0+fac)*dd(:,3,m)
!!$          end do
!!$       end select
!!$    end do
!!$  end subroutine implicitation_cond_lim_mat_subsonique
!!$
!!$  subroutine implicitation_cond_lim_sm(b, X, Sm)
!!$    use m_MPI
!!$    implicit none
!!$
!!$    type (STR_BLOC), intent(in) ::  b
!!$    real*8, dimension(3,b%md-2:b%mf+2), intent(in) :: X
!!$    real*8, dimension(3,b%md:b%mf), intent(inout) :: Sm
!!$
!!$    real*8, dimension(:,:,:), pointer :: mat_extradiag
!!$    type (STR_PATCH), pointer :: patch
!!$    integer :: iface, ip, md, mf
!!$    integer :: n_ref, n_ext
!!$
!!$    md = b%md; mf = b%mf
!!$
!!$    do iface = 3, 4
!!$       do ip = 1, b%face(iface)%nb_patch
!!$          patch => b%face(iface)%patch(ip)
!!$          select case (trim(patch%bc%position))
!!$          case("North")
!!$             n_ref = mf
!!$             n_ext = mf+1
!!$             mat_extradiag => b%matrice%c                                        
!!$
!!$          case("South")
!!$             n_ref = md
!!$             n_ext = md-1
!!$             mat_extradiag => b%matrice%a
!!$
!!$          case default
!!$             print*, "Face inconnue pour l'implicitation du second membre ", patch%bc%position
!!$             call arret_code
!!$          end select
!!$
!!$
!!$          select case (trim(patch%bc%condition))
!!$          case("flux_nul", "paroi", "symetrie")
!!$             ! on ne fait rien ici
!!$             ! deja fait dans implicitation_cond_lim_mat
!!$
!!$          case("entree_supersonique")
!!$             Sm(:,n_ref) = Sm(:,n_ref) - matmul(mat_extradiag(:,:,n_ref), X(:,n_ext)) 
!!$
!!$          case default
!!$             print*, "Condition inconnue pour l'implicitation du second membre ", patch%bc%condition
!!$             call arret_code
!!$          end select
!!$
!!$       end do
!!$    end do
!!$
!!$  end subroutine implicitation_cond_lim_sm
end module SOLVEUR_balayage
