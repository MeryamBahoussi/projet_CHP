module chimie

  use m_struct
  use m_MPI
  implicit none

  logical :: Glenn = .false.
  
contains
  
!!!-----------------------------------------------!!!
!!!  Grandeurs utiles pour l'equilibre chimique   !!!
!!!    Constante de réaction                      !!!
!!!    Enthalpie libre standard de reaction       !!!
!!!    et leur derivees                           !!!
!!!-----------------------------------------------!!!
!!! calcul de la constante de reaction
  function Kp(reaction, T)
    use parametres_fluide, only : ONEATM, R_gaz
    
    type(STR_REACTION) :: reaction
    real*8             :: T
    real*8             :: A
    real*8             :: Kp
  
    A = delta_G0(reaction,T)
    Kp = ONEATM**reaction%D_nu * exp(-A/(R_gaz*T))
  end function Kp

!!! calcul de la derivee de la constante de reaction 
  Function dKp(reaction, T)
    use parametres_fluide, only : ONEATM, R_gaz
    
    type(STR_REACTION) :: reaction
    real*8             :: T
    real*8             :: dKp
    
    real*8 :: A, B
    
    A = delta_G0(reaction, T)
    B = ddelta_G0(reaction, T)

    dKp = ONEATM**reaction%D_nu * exp(-A/(R_gaz*T))*(R_gaz*(A-T*B))/(R_gaz*T)**2
  end function dKp
 
!!! calcul de l'enthalpie libre standard de reaction  
  function delta_G0(reaction, T)
    type(STR_REACTION) :: reaction
    real*8             :: T
    real*8             :: delta_G0
    
    integer :: i
    
    delta_G0 = 0.d0
    do i = 1, reaction%nesp
       delta_G0 = delta_G0 + reaction%nu(i) * G0(fluide1%especes(reaction%ist(i)), T)
    end do
  end function delta_G0

!!! calcul de la dérivée de l'enthalpie libre standard de reaction
  function ddelta_G0(reaction, T)
    type(STR_REACTION) :: reaction
    real*8             :: T
    real*8             :: ddelta_G0
    
    integer :: i
   
    ddelta_G0 = 0.d0
    do i = 1, reaction%nesp
       ddelta_G0 = ddelta_G0 + reaction%nu(i) * dG0(fluide1%especes(reaction%ist(i)), T) 
    end do
  end function ddelta_G0
  
!!!------------------------------!!!
!!!    Modele thermodynamique    !!!
!!!------------------------------!!!
!!! calcul de l'energie interne pour toutes les especes
  function calc_Eint_esp(especes, T) result (Eint)
    use parametres_fluide, only : R_gaz
   
    type(STR_ESPECE) :: especes(:)
    real*8           :: T
    real*8           :: Eint(size(especes))

    integer :: i
    real*8  :: Cvi

    if (.not. Glenn) then
       do i = 1, size(especes)
          Cvi = Cv_esp(especes(i), T)
          Eint(i) = especes(i)%h0 + Cvi * T
       end do
    else
       do i = 1, size(especes)
          Eint(i) = (Glenn_enthalpie(especes(i), T) - R_gaz * T)/especes(i)%masse
       end do
    end if
  end function calc_Eint_esp

!!! calcul de la capacite calorifique massique a volume constant pour une espece donnee
  function Cv_esp(espece, T) result (Cv)
    type(STR_ESPECE) :: espece
    real*8           :: T
    real*8           :: Cv
  
    if (.not. Glenn) then
       ! Cv = (2.5 ou 1.5)* R_gaz/masse
       Cv = (1.5d0 + 1.d0*espece%diatomique) * espece%R
    else    
       Cv = Glenn_cp(espece, T) - espece%R
    end if

  end function Cv_esp

!!! retrouve le numero du polynome de Glenn pour une temperature donnee
  function Glenn_find_num_polynome(espece, T) result (p)
    type(STR_ESPECE) :: espece
    real*8           :: T
    integer          :: p
    
    integer :: k
    logical :: ok
    
    if( T < 1000.d0) then
       p = 1
    else if (T >= 1000.d0 .and. T < 60000.d0) then
       p = 2
    else if (T >= 60000.d0) then
       p = 3
    end if
    
!!$    ok = .false.
!!$    ! boucle sur les polynomes de Glenn : nombre de plages de temperature
!!$    do k = 1, ubound(espece%coeff_glenn,2)
!!$       if( tmin(j,ie) <= T .and. T <= tmax(j,ie)) then
!!$          p = k
!!$          ok = .true.
!!$          exit
!!$       end if
!!$    end do
!!$    if (ok .eqv. .false.) then
!!$       print*, "T hors des plages de donnees pour Glenn"
!!$       call arret_code
!!$    end if
  end function Glenn_find_num_polynome

!!! Calcul de Cp avec les polynomes de Glenn (NASA-9)
  function Glenn_Cp(espece, T) result (Cp)
    type(STR_ESPECE) :: espece
    real*8           :: T
    real*8           :: Cp

    integer         :: p
    real*8          :: Tm2, Tm1, T2, T3, T4
    real*8, pointer :: A(:,:)
    
    ! on retrouve le polynome de Glenn pour un temperature donnee
    p = Glenn_find_num_polynome(espece, T)
    
    Tm1 = 1.d0/T
    Tm2 = Tm1*Tm1
    T2 = T*T
    T3 = T*T2
    T4 = T*T3
    
    A => espece%coeff_glenn

    Cp =   A(1,p)*Tm2 &
         + A(2,p)*Tm1 &
         + A(3,p)     &
         + A(4,p)*T   &
         + A(5,p)*T2  &
         + A(6,p)*T3  &
         + A(7,p)*T4 
    Cp = Cp * espece%R
  end function Glenn_Cp
  
!!! Calcul de l'enthalpie avec les polynomes de Glenn (NASA-9)
!!! h_i = \int Cp_i dT + b1*R_i
  function Glenn_enthalpie(espece, T) result (h)
    type(STR_ESPECE) :: espece
    real*8           :: T
    real*8           :: h
       
    integer         :: p
    real*8          :: Tm1, T2, T3, T4, T5
    real*8, pointer :: A(:,:)
    real*8, pointer :: b1(:)
        
    ! on retrouve le polynome de Glenn pour un temperature donnee
    p = Glenn_find_num_polynome(espece, T)

    Tm1 = 1.d0/T
    T2 = T*T
    T3 = T*T2
    T4 = T*T3
    T5 = T*T4
    
    A => espece%coeff_glenn
    b1 => espece%coeff_C0
    
    h =  - A(1,p) * Tm1       &
         - A(2,p) * log(T)    &
         + A(3,p) * T         &
         + A(4,p) * T2*0.5d0  &
         + A(5,p) * T3/3.d0   &
         + A(6,p) * T4*0.25d0 &
         + A(7,p) * T5*0.2d0  &
         + b1(p)
    h = h * espece%R
  end function Glenn_enthalpie

!!! Calcul de l'entropie avec les polynomes de Glenn (NASA-9)
!!! s_i = \int Cp_i/T dT + b2*R_i
  function Glenn_entropie(espece, T) result (s)  
    type(STR_ESPECE) :: espece
    real*8           :: T
    real*8           :: s
    
    integer         :: p
    real*8          :: Tm2, Tm1, T2, T3, T4
    real*8, pointer :: A(:,:)
    real*8, pointer :: b2(:)
    
    A => espece%coeff_glenn
    b2 => espece%coeff_C1
    
    ! on retrouve le polynome de Glenn pour un temperature donnee
    p = Glenn_find_num_polynome(espece, T)
    
    Tm1 = 1.d0/T
    Tm2 = Tm1*Tm1
    T2 = T*T
    T3 = T*T2
    T4 = T*T3
    
    s =  - A(1,p) * Tm2*0.5d0 &
         - A(2,p) * Tm1       &
         + A(3,p) * log(T)    &
         + A(4,p) * T         &
         + A(5,p) * T2*0.5d0  &
         + A(6,p) * T3/3.d0   &
         + A(7,p) * T4*0.25d0 &
         + b2(p)
    s = s * espece%R 
  end function Glenn_entropie
  
!!! calcul de l'enthalpie libre standard d'une espece
  function G0(espece, T)
    type(STR_ESPECE) :: espece
    real*8           :: T
    real*8           :: G0

    real*8 :: S0, H0, Cv
    character(len=16) :: toto

    if (.not. glenn) then
       Cv = Cv_esp(espece, T)
       H0 = (espece%h0 + (Cv + espece%R) * T) * espece%masse   !!! ??? *masse ???
       S0 = (espece%s0 + (Cv + espece%R) * log(T)) * espece%masse
    else
       H0 = Glenn_enthalpie(espece, T)
       S0 = Glenn_entropie(espece, T)
    end if
    G0 = H0 - T * S0
    
    ! on oblige le compilateur a vraiment calculer H0
    ! sans ca, on a pas les memes resultats avec et sans optimisation du compilateur ...
    ! write(toto,'(e16.7)') H0 
  end function G0

!!! calcul de la derivee de l'energie libre standard
  function dG0(espece, T)
    type(STR_ESPECE) :: espece
    real*8           :: T
    real*8           :: dG0

    real*8 :: S0, dS0, dH0, Cv

    if (.not. glenn) then
       Cv = Cv_esp(espece, T)
       dH0 = (Cv + espece%R) * espece%masse
       S0  = (espece%s0 + (Cv + espece%R) * log(T)) * espece%masse
       dS0 = ((Cv + espece%R)/T) * espece%masse
       dG0 = dH0 - T * dS0 - S0
    else
!!! dG0 = dH0 - T * dS0 - S0 = -S0
       S0 = Glenn_entropie(espece, T)
       dG0 = - S0
    end if

  end function dG0

!!!------------------------------------------!!!
!!! Routines de conversion                   !!!
!!! c_ele : fractions massiques en elements  !!!
!!! c_esp : fractions massiques en especes   !!!
!!! x_ele : fractions molaires en especes    !!!
!!! x_esp : fractions molaires en especes    !!!
!!!------------------------------------------!!!
  subroutine cesp_from_xesp(especes,x_esp,c_esp)
    type(STR_ESPECE) , dimension(:), intent(in)  :: especes
    real*8           , dimension(:), intent(in)  :: x_esp
    real*8           , dimension(:), intent(out) :: c_esp
    
    integer   :: i
    real*8    :: somme
    
    somme = 0.d0
    do i = 1, size(especes)
       c_esp(i) = especes(i)%masse * x_esp(i)
       somme = somme + c_esp(i)
    end do
    c_esp(:) = c_esp(:)/somme
  end subroutine cesp_from_xesp

  subroutine cele_from_xele(elements,x_ele,c_ele)
    type(STR_ELEMENT), dimension(:), intent(in)  :: elements
    real*8           , dimension(:), intent(in)  :: x_ele
    real*8           , dimension(:), intent(out) :: c_ele

    integer   :: k
    real*8    :: somme

    somme = 0.d0
    do k = 1, size(elements)
       c_ele(k) = elements(k)%masse * x_ele(k)
       somme = somme + c_ele(k)
    end do
    c_ele(:) = c_ele(:)/somme
  end subroutine cele_from_xele

  subroutine xesp_from_cesp(especes,c_esp,x_esp)
    type(STR_ESPECE) , dimension(:), intent(in)  :: especes
    real*8           , dimension(:), intent(in)  :: c_esp
    real*8           , dimension(:), intent(out) :: x_esp

    integer   :: i
    real*8    :: somme

    somme = 0.d0
    do i = 1, size(especes)
       x_esp(i) = c_esp(i)/especes(i)%masse
       somme = somme + x_esp(i)
    end do
    x_esp(:) = x_esp(:) / somme
  end subroutine xesp_from_cesp

  subroutine xele_from_cele(elements,c_ele,x_ele)
    type(STR_ELEMENT), dimension(:), intent(in)  :: elements
    real*8           , dimension(:), intent(in)  :: c_ele
    real*8           , dimension(:), intent(out) :: x_ele

    integer   :: k
    real*8    :: somme

    somme = 0.d0
    do k = 1, size(elements)
       x_ele(k) = c_ele(k)/elements(k)%masse
       somme = somme + x_ele(k)
    end do
    x_ele(:) = x_ele(:)/somme
  end subroutine xele_from_cele

  subroutine cele_from_cesp(elements,especes,c_esp,c_ele)
    type(STR_ELEMENT), intent(in)  :: elements(:)
    type(STR_ESPECE) , intent(in)  :: especes(:)
    real*8           , intent(in)  :: c_esp(1:size(especes))
    real*8           , intent(out) :: c_ele(1:size(elements))

    integer :: i,k,iesp

    do k = 1, size(elements) 
       c_ele(k) = 0.d0
       do i = 1, size(elements(k)%liste_especes(:))
          iesp = elements(k)%liste_especes(i)
          c_ele(k) = c_ele(k) + elements(k)%multiplicite(i)*elements(k)%masse/especes(iesp)%masse * c_esp(iesp)
       end do
    end do

    do k = 1, size(elements)
       if (c_ele(k) < 0.d0) c_ele(k) = 0.d0
    end do
    if (sum(c_ele(:))/=0.d0) c_ele(:) = c_ele(:)/sum(c_ele(:))

  end subroutine cele_from_cesp

  subroutine xele_from_xesp(elements,especes,x_esp,x_ele)
    type(STR_ELEMENT), intent(in)  :: elements(:)
    type(STR_ESPECE) , intent(in)  :: especes(:)
    real*8           , intent(in)  :: x_esp(1:size(especes))
    real*8           , intent(out) :: x_ele(1:size(elements))

    integer   :: i,k,iesp
    real*8    :: somme

    somme = 0.d0
    do k = 1, size(elements) 
       x_ele(k) = 0.d0
       do i = 1, size(elements(k)%liste_especes(:))
          iesp = elements(k)%liste_especes(i)
          x_ele(k) = x_ele(k) + elements(k)%multiplicite(i) * x_esp(iesp)
       end do
       somme = somme + x_ele(k)
    end do
    x_ele(:) = x_ele(:)/somme

    if (abs(sum(x_ele(:)) -1.d0) > 1.d-15) then
       print*, "Erreur xele_from_xesp"
       print*, sum(x_ele(:)), abs(sum(x_ele(:)) -1.d0)
       call arret_code
    end if
  end subroutine xele_from_xesp

end module chimie
