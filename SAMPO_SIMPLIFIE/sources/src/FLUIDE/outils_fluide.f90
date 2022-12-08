module outils_fluide
  use m_struct
  use m_MPI
  use parametres_globaux,only : nMF

  implicit none
  
  public save_Un, save_Un_bloc
  public save_Urk, save_Urk_bloc, combine_SSPRK, combine_SSPRK_bloc
  public rho_from_U
  public recup_ind_esp, recup_ind_ele
  public val_abs
  public moyenne, produit_tensoriel
  private

contains
!!! sauvegarde le vecteur conservatif et la pression a l'instant tn
!!! utile pour les conditions subsonique, le calcul du residu
!!! et si la condition CFL n est pas verifiee pour le transport
  subroutine save_Un(sb)
    type (STR_SUPER_BLOC), pointer :: sb
    
    integer :: i, ib
    type (STR_BLOC), pointer :: b
    
    do i = 1, sb%nb_membre 
       ib = sb%membre(i)
       b => acces_bloc(ib)%pointeur
       if (bloc_est_local (ib)) then
          call save_Un_bloc(b)
       end if
    end do
  end subroutine save_Un
  
  subroutine save_Un_bloc(b)
    type (STR_BLOC), pointer :: b
    
    integer :: l, m
    
    do m = b%md-nMF, b%mf+nMF
       do l = b%ld-nMF, b%lf+nMF
          b%U_n(:,l,m) = b%U(:,l,m) 
          b%p_n(l,m) = b%p(l,m)
       end do
    end do
  end subroutine save_Un_bloc

    subroutine save_Urk(sb)
    type (STR_SUPER_BLOC), pointer :: sb
    
    integer :: i, ib
    type (STR_BLOC), pointer :: b
    
    do i = 1, sb%nb_membre 
       ib = sb%membre(i)
       b => acces_bloc(ib)%pointeur
       if (bloc_est_local (ib)) then
          call save_Urk_bloc(b)
       end if
    end do
  end subroutine save_Urk
  
  subroutine save_Urk_bloc(b)
    type (STR_BLOC), pointer :: b
    
    integer :: l, m
    
    do m = b%md-nMF, b%mf+nMF
       do l = b%ld-nMF, b%lf+nMF
          b%U_rk(:,l,m) = b%U(:,l,m) 
       end do
    end do
  end subroutine save_Urk_bloc

  subroutine combine_SSPRK(sb,a1,a2)
    type (STR_SUPER_BLOC), pointer :: sb
    real*8, intent(in) :: a1, a2

    integer :: i, ib
    type (STR_BLOC), pointer :: b
    
    do i = 1, sb%nb_membre 
       ib = sb%membre(i)
       b => acces_bloc(ib)%pointeur
       if (bloc_est_local (ib)) then
          call combine_SSPRK_bloc(b,a1,a2)
       end if
    end do
  end subroutine combine_SSPRK

  subroutine combine_SSPRK_bloc(b,a1,a2)
    type (STR_BLOC), pointer :: b
    real*8, intent(in) :: a1, a2

    integer :: l, m
    
    do m = b%md-nMF, b%mf+nMF
       do l = b%ld-nMF, b%lf+nMF
          b%U(:,l,m) = a1*b%U_rk(:,l,m) + a2*b%U(:,l,m)
       end do
    end do
  end subroutine combine_SSPRK_bloc
  
  ! calcul de la densite a partir des variables conservatives
  subroutine rho_from_U(U,rho)
    use parametres_fluide

    real*8, dimension(:), intent(in)  :: U
    real*8              , intent(out) :: rho

    rho = sum(U(irho1z1:irho2z2))
  end subroutine rho_from_U
  
!!! Trouve l'indice correspondant a une espece
  subroutine recup_ind_esp(nom_esp, ind)
    use m_MPI
    use m_struct, only : fluide1

    character(*), intent(in)  :: nom_esp
    integer     , intent(out) :: ind
    
    integer :: i
    
    ind = -1
    do i = 1, fluide1%nb_esp
       if (trim(fluide1%especes(i)%nom) == trim(nom_esp)) then
          ind = i
          exit
       end if
    end do
    if (ind == -1) then
       print*, "Erreur dans recup_ind_esp : l'espece suivante n'existe pas ", nom_esp
       call arret_code
    end if
  end subroutine recup_ind_esp

 !!! Trouve l'indice correspondant a un element
  subroutine recup_ind_ele(nom_ele, ind)
    use m_MPI
    use m_struct, only : fluide1

    character(*), intent(in)  :: nom_ele
    integer     , intent(out) :: ind
    
    integer :: k
    
    ind = -1
    do k = 1, fluide1%nb_ele
       if (trim(fluide1%elements(k)%nom) == trim(nom_ele)) then
          ind = k
          exit
       end if
    end do
    if (ind == -1) then
       print*, "Erreur dans recup_ind_ele : l'element suivant n'existe pas ", nom_ele
       call arret_code
    end if
  end subroutine recup_ind_ele
  
  function val_abs(epsilon,val)
    ! Correspond a la fonction Q_visc dans second_euler.f90 du code Aero du CESTA
    implicit none
    real*8 :: val_abs, epsilon, val

    val_abs = abs(val) + 0.5d0 * ( dim(epsilon, abs(val)))**2.d0 / max(epsilon, 1.d-60)
  end function val_abs

  
  ! moyenne entre deux reels
  function moyenne(a, b)
    implicit none

    real*8 :: a, b
    real*8 :: moyenne

    ! moyenne harmonique
    if (a==0.d0 .and. b==0.d0) then
       print*, "moyenne a=0 et b=0"
       call arret_code
       moyenne = 0.d0
    else
       moyenne = 2.d0*a*b/(a+b)
    end if

!!$    moyenne = 0.5d0*(a+b)

  end function moyenne
  
  function produit_tensoriel(dim,t1,t2) result (t1t2)
    integer                    :: dim
    real*8 ,dimension(dim)     :: t1, t2
    real*8 ,dimension(dim,dim) :: t1t2

    integer :: i,j

    do j = 1,dim
       do i = 1,dim
          t1t2(i,j) = t1(i)*t2(j)
       end do
    end do
  end function produit_tensoriel
  
end module outils_fluide
