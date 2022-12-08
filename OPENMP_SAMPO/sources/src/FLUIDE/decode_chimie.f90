module decodage_chimie

  use chimie
  use decodage_diphasique
  use m_transport_diffusion
  use Multiespeces

  implicit none

contains
  
  subroutine decode_chimie(ldeb, lfin, mdeb, mfin, b, ierreur)
    integer        , intent(in)  :: ldeb, lfin, mdeb, mfin
    type (STR_BLOC), pointer     :: b
    integer        , intent(out) :: ierreur
    
    integer :: l, m, i
    real*8  :: Cp, Cvi, somme, eps

    do m = mdeb, mfin
       do l = ldeb, lfin
          
          b%c_esp(:,l,m) = 0.d0
          somme = 0.d0
          do i = 1, b%ne
             b%c_esp(i,l,m) = b%U(irho1z1+i-1,l,m)  ! ci = rho*y*ci
             b%c_esp(i,l,m) = max(b%c_esp(i,l,m),0.d0)
             somme = somme + b%c_esp(i,l,m)
          end do
          b%c_esp(:,l,m) = b%c_esp(:,l,m)/somme    ! ci = rho*y*ci/(rho*y)
          
          call calcul_Cv_Cp_fige(fluide1, fluide1%especes, b%c_esp(:,l,m), b%cv1(l,m), Cp)
                         
          b%gamma1(l,m)   = fluide1%EOS%gamma
          b%gamma_eq(l,m) = fluide1%EOS%gamma
          b%Cv(l,m) = b%cv1(l,m)
          
!!! decodage de la pression
          eps = b%U(irhoE,l,m)/b%rho(l,m) - 0.5d0*(b%u_x(l,m)**2+b%u_y(l,m)**2) 
          b%p(l,m) = MESP_p_from_rho_eps(b%rho(l,m), eps, b%gamma1(l,m))
        
!!! vitesse du son (on est en monophasique : pi = 0 et z=1)
          call vitesse_du_son_melange_Allaire(b%gamma1(l,m), 0.d0, b%p(l,m), b%rho(l,m), 1.d0, b%c(l,m), ierreur)
       
!!! Temperature
          b%T(l,m) = MESP_T_from_p_rho(b%p(l,m),b%rho(l,m),b%cv1(l,m),b%gamma1(l,m))
          
!!! Enthalpies spécifiques
!!! Ici on suppose que l'enthalpie de formation de toutes les especes est nulle
          b%Hi(:,l,m) = 0.d0
          do i = 1, b%ne
             Cvi = Cv_esp(fluide1%especes(i), b%T(l,m))
             b%Hi(i,l,m) = ( Cvi + fluide1%especes(i)%R ) * b%T(l,m)
          end do

       end do
    end do    
  end subroutine decode_chimie
  
  subroutine decode_chimie_loc(U, cell, ierreur)
    real*8, dimension(:), intent(in)    :: U
    type (STR_cell)     , intent(inout) :: cell
    integer             , intent(out)   :: ierreur
    
    integer :: i
    real*8  :: Cv, Cp, Cvi, somme, eps
    real*8  :: mu, K, D
    
    cell%ci(:) = 0.d0
    somme = 0.d0
    do i = 1, fluide1%ne
       cell%ci(i) = U(irho1z1+i-1)  ! ci = rho*y*ci
       cell%ci(i) = max(cell%ci(i),0.d0)
       somme = somme + cell%ci(i)
    end do
    cell%ci(:) = cell%ci(:)/somme    ! ci = rho*y*ci/(rho*y)

    call calcul_Cv_Cp_fige(fluide1, fluide1%especes, cell%ci(:), Cv, Cp)

!!! decodage de la pression
    eps = U(irhoE)/cell%rho - 0.5d0*(cell%u**2+cell%v**2) 
    cell%p = MESP_p_from_rho_eps(cell%rho, eps, Cp/Cv)
        
!!! vitesse du son (on est en monophasique : pi = 0 et z=1)
    call vitesse_du_son_melange_Allaire(fluide1%EOS%gamma, 0.d0, cell%p, cell%rho, 1.d0, cell%c, ierreur)
    
!!! Temperature
    cell%T = MESP_T_from_p_rho(cell%p, cell%rho, Cv, fluide1%EOS%gamma)
    
!!! Enthalpies spécifiques
!!! Ici on suppose que l'enthalpie de formation de toutes les especes est nulle
    do i = 1, fluide1%ne
       Cvi = Cv_esp(fluide1%especes(i), cell%T)
       !cell%Hi(i) = ( Cvi + fluide1%especes(i)%R ) * cell%T
    end do

  end subroutine decode_chimie_loc

  subroutine calcul_Cv_Cp_fige(fld, especes, c_esp, Cv, Cp)
    use chimie, only : Cv_esp
    type(STR_fluide), intent(inout) :: fld
    type(STR_ESPECE), intent(in)    :: especes(:)
    real*8          , intent(in)    :: c_esp(1:size(especes))
    real*8          , intent(out)   :: Cv, Cp

    integer :: i
    real*8  :: R, Cvi 

    R = 0.d0
    Cv = 0.d0
    do i = 1, size(especes)
       R = R + c_esp(i) * especes(i)%R
       Cvi = Cv_esp(especes(i), -1.d70)  ! il faut passer T en argument ??
       Cv = Cv + c_esp(i) * Cvi
    end do
    Cp = R + Cv

!!! on stocke dans fld pour pouvoir utilise grandeurs_melange
    fld%EOS%gamma = Cp/Cv
    fld%EOS%pi = 0.d0  ! melange de gaz parfait
    fld%EOS%q = 0.d0
    fld%EOS%Cv = Cv
    fld%EOS%Cp = Cp
  end subroutine calcul_Cv_Cp_fige

end module decodage_chimie
