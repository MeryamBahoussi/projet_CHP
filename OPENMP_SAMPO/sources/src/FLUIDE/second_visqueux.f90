module second_Visqueux
  use m_struct
  use a_mettre_en_donnees
  use parametres_fluide
  use parametres_globaux
  use outils_fluide, only : moyenne
  use m_reconstruction_gradient

  implicit none
   
  real*8, parameter :: F2S3 = 2.d0 / 3.d0
  real*8, parameter :: F4S3 = 4.d0 / 3.d0

contains

  subroutine second_membre_visqueux(b, sm)
    implicit none

    type (STR_BLOC)                   , pointer       :: b
    real*8, dimension(1:,b%ld:, b%md:), intent(inout) :: sm

    integer :: l, m
    real*8, dimension(:,:,:), pointer :: flux

    flux => b%Fg_acoustique ! la composante en z n'est pas modifiee
    flux = 0.d0

    call flux_visqueux_x(b, flux) 

    do m = b%md, b%mf
       do l = b%ld, b%lf
          sm(:iz-1,l,m) = sm(:iz-1,l,m) - (flux(:iz-1,l,m)-flux(:iz-1,l-1,m))
       end do
    end do
    
    flux = 0.d0
    call flux_visqueux_y(b, flux)

    do m = b%md, b%mf
       do l = b%ld, b%lf
          sm(:iz-1,l,m) = sm(:iz-1,l,m) - (flux(:iz-1,l,m)-flux(:iz-1,l,m-1))
       end do
    end do
        
  end subroutine second_membre_visqueux

!!!-----------------------------------------------
!!! calcul du flux lié aux termes de dissipation
!!!-----------------------------------------------
!!! visqueux_bi_temperature : q = -\sum K_k \nabla T_k
!!! sinon                   : q = -K \nabla T avec K(Kk) et T(Cvk,yk,Tk) (cf decodage)
!!!
!!! le visqueux_bi_temperature ne marche pas (et ne doit pas etre utilise)
!!!
  subroutine flux_visqueux_x(b, flux)
  ! use variables_debug, only : cmp_manu
   implicit none

    type(STR_BLOC), pointer :: b
    real*8, dimension(1:b%neq-1,b%ld-1:b%lf, b%md-1:b%mf), intent(inout) :: flux

    integer :: i, j, l, m, lg, mg, ld, md
    real*8 :: X_x, Y_y, X_y ,Y_x
    real*8 :: tau_xx, tau_xy, tau_yy
    real*8 :: Jac_g, Jac_d, Jacobien_moy
    real*8 :: u_moy, v_moy, mu_moy, K_moy
    real*8 :: D_moy, rhoy_moy, Hi_moy                 ! (chimie figee)

    real*8, dimension(1:b%neq-1) :: ev, fv            ! Termes visqueux
    type( STR_gradient ) :: grad_u, grad_v, grad_T
    type( STR_gradient ) :: grad_Ci(1:b%ne)           ! (chimie figee)
    type( STR_gradient ) :: grad_Cj(1:fluide1%nb_esp) ! (Equilibre Chimique)

    real*8, dimension(:,:)  , pointer :: u, v, T
    real*8, dimension(:,:,:), pointer :: ci
    real*8, dimension(:,:,:), pointer :: cj           ! (Equilibre chimique)
    real*8, dimension(:,:,:), pointer :: metrique_grad

    u => b%u_x
    v => b%u_y
    T => b%T
    metrique_grad => b%metrique_grad_l

    select case (fluide1%type)
    case(EOS, CHIMIE_FIGEE)
       ci => b%c_esp
    end select

    do m = b%md, b%mf
       do l = b%ld-1, b%lf

          lg = l  ; mg = m  ! cell_g = cell(l,m)
          ld = l+1; md = m  ! cell_d = cell(l+1,m)

          ev = 0.d0 ; fv = 0.d0

          ! Grandeurs moyennes a l'interface
          u_moy = .5d0 * ( u(lg,mg) + u(ld,md) )
          v_moy = .5d0 * ( v(lg,mg) + v(ld,md) )

          K_moy = moyenne(b%K(lg,mg), b%K(ld,md))
          mu_moy = moyenne(b%mu(lg,mg), b%mu(ld,md))

          ! Jacobiens
          Jac_g = 1.d0 / b%aire(lg,mg)
          Jac_d = 1.d0 / b%aire(ld,md)
          Jacobien_moy = moyenne(Jac_g, Jac_d)

          ! Gradient a l'interface sur le maillage curviligne (repere (x,y))
          call calc_grad_face_x(metrique_grad(:,l,m),grad_u,u(l-1:l+1,m-1:m+1),l,m)
          call calc_grad_face_x(metrique_grad(:,l,m),grad_v,v(l-1:l+1,m-1:m+1),l,m)
          call calc_grad_face_x(metrique_grad(:,l,m),grad_T,T(l-1:l+1,m-1:m+1),l,m)

          ! tenseur des contraintes
          tau_xx = mu_moy * ( 4.d0/3.d0 * grad_u%x - 2.d0/3.d0 * grad_v%y )
          tau_yy = mu_moy * ( 4.d0/3.d0 * grad_v%y - 2.d0/3.d0 * grad_u%x )
          tau_xy = mu_moy * ( grad_u%y + grad_v%x )

          ev(iu) = tau_xx ; fv(iu) = tau_xy
          ev(iv) = tau_xy ; fv(iv) = tau_yy
          ev(iE) = u_moy * tau_xx + v_moy * tau_xy + K_moy * grad_T%x
          fv(iE) = u_moy * tau_xy + v_moy * tau_yy + K_moy * grad_T%y

          ! chimie figee
          ! div(rho y D grad(c_i)) pour l'equation en c_i
          ! div(sum_1^ne rho y D h_i grad c_i)
          if ( fluide1%type == CHIMIE_FIGEE ) then

             D_moy = moyenne(b%D(lg,mg), b%D(ld,md))
             rhoy_moy = .5d0 * ( b%y(lg,mg)*b%rho(lg,mg) + b%y(ld,md)*b%rho(ld,md) )

             do i = 1, b%ne
                ! TODO : !!  mettre les CL de coins pour les ci
                Hi_moy = .5d0 * ( b%Hi(i,lg,mg) + b%Hi(i,ld,md) )
                call calc_grad_face_x(metrique_grad(:,l,m),grad_Ci(i),ci(i,l-1:l+1,m-1:m+1),l,m)

                ev(irho1z1+i-1) = D_moy * rhoy_moy * grad_Ci(i)%x
                fv(irho1z1+i-1) = D_moy * rhoy_moy * grad_Ci(i)%y

                ev(iE) = ev(iE) + Hi_moy * ev(irho1z1+i-1)
                fv(iE) = fv(iE) + Hi_moy * fv(irho1z1+i-1)
             end do
          end if

          flux(:,l,m) = flux(:,l,m) + Jacobien_moy * ( metrique_grad(1,l,m) * ev(:) &
                                                     + metrique_grad(2,l,m) * fv(:) )

       end do
    end do
  end subroutine flux_visqueux_x

!!!-----------------------------------------------
!!! calcul du flux lié aux termes de dissipation
!!!-----------------------------------------------
!!! visqueux_bi_temperature : q = -\sum K_k \nabla T_k
!!! sinon                   : q = -K \nabla T avec K(Kk) et T(Cvk,yk,Tk) (cf decodage)
!!!
!!! le visqueux_bi_temperature ne doit pas etre utilise
!!!
  subroutine flux_visqueux_y(b, flux)
    !use variables_debug, only : cmp_manu
    implicit none

    type(STR_BLOC), pointer :: b
    real*8, dimension(1:b%neq-1,b%ld-1:b%lf, b%md-1:b%mf), intent(inout) :: flux

    integer :: i, j, l, m, lg, mg, ld, md
    real*8 :: tau_xx, tau_xy, tau_yy
    real*8 :: Jac_g, Jac_d, Jacobien_moy
    real*8 :: u_moy, v_moy, mu_moy, K_moy
    real*8 :: Ci_x(1:b%ne), Ci_y(1:b%ne) ! (chimie figee)
    real*8 :: Cj_x(1:fluide1%nb_esp), Cj_y(1:fluide1%nb_esp) ! (Equilibre chimique)
    real*8 :: D_moy, rhoy_moy, Hi_moy                        ! (chimie figee)

    real*8, dimension(1:b%neq-1) :: ev, fv ! Termes visqueux
    type( STR_gradient ) :: grad_u, grad_v, grad_T
    type( STR_gradient ) :: grad_Ci(1:b%ne) ! (chimie figee)
    type( STR_gradient ) :: grad_Cj(1:fluide1%nb_esp) ! (Equilibre Chimique)

    real*8, dimension(:,:)  , pointer :: u, v, T
    real*8, dimension(:,:,:), pointer :: ci
    real*8, dimension(:,:,:), pointer :: cj           ! (Equilibre chimique)
    real*8, dimension(:,:,:), pointer :: metrique_grad

    u => b%u_x
    v => b%u_y
    T => b%T
    metrique_grad => b%metrique_grad_m


    select case (fluide1%type)
    case(EOS, CHIMIE_FIGEE)
       ci => b%c_esp
    end select

    do m = b%md-1, b%mf
       do l = b%ld, b%lf

          lg = l  ; mg = m    ! cell_g = cell(l,m)
          ld = l  ; md = m+1  ! cell_d = cell(l,m+1)

          ev = 0.d0 ; fv = 0.d0

          ! Grandeurs moyennes a l'interface
          u_moy = .5d0 * ( u(lg,mg) + u(ld,md) )
          v_moy = .5d0 * ( v(lg,mg) + v(ld,md) )

          K_moy = moyenne(b%K(lg,mg), b%K(ld,md))
          mu_moy = moyenne(b%mu(lg,mg), b%mu(ld,md))

          Jac_g = 1.d0 / b%aire(lg,mg)
          Jac_d = 1.d0 / b%aire(ld,md)
          Jacobien_moy = moyenne(Jac_g, Jac_d)

          ! Gradient a l'interface sur le maillage curviligne (repere (x,y))
          call calc_grad_face_y(metrique_grad(:,l,m),grad_u,u(l-1:l+1,m-1:m+1),l,m)
          call calc_grad_face_y(metrique_grad(:,l,m),grad_v,v(l-1:l+1,m-1:m+1),l,m)
          call calc_grad_face_y(metrique_grad(:,l,m),grad_T,T(l-1:l+1,m-1:m+1),l,m)


          tau_xx = mu_moy * ( 4.d0/3.d0 * grad_u%x  - 2.d0/3.d0 * grad_v%y )
          tau_yy = mu_moy * ( 4.d0/3.d0 * grad_v%y  - 2.d0/3.d0 * grad_u%x )
          tau_xy = mu_moy * ( grad_u%y + grad_v%x )

          ev(iu) = tau_xx ; fv(iu) = tau_xy
          ev(iv) = tau_xy ; fv(iv) = tau_yy

          ev(iE) = u_moy * tau_xx + v_moy * tau_xy + K_moy * grad_T%x
          fv(iE) = u_moy * tau_xy + v_moy * tau_yy + K_moy * grad_T%y

          ! chimie figee
          if ( fluide1%type == CHIMIE_FIGEE ) then

             D_moy = moyenne(b%D(lg,mg), b%D(ld,md))
             rhoy_moy = .5d0 * ( b%y(lg,mg)*b%rho(lg,mg) + b%y(ld,md)*b%rho(ld,md) )

             do i = 1, b%ne
                Hi_moy = .5d0 * ( b%Hi(i,lg,mg) + b%Hi(i,ld,md) )

                call calc_grad_face_y(metrique_grad(:,l,m),grad_Ci(i),ci(i,l-1:l+1,m-1:m+1),l,m)

                ev(irho1z1+i-1) = D_moy * rhoy_moy * grad_Ci(i)%x
                fv(irho1z1+i-1) = D_moy * rhoy_moy * grad_Ci(i)%y

                ev(iE) = ev(iE) + Hi_moy * ev(irho1z1+i-1)
                fv(iE) = fv(iE) + Hi_moy * fv(irho1z1+i-1)
             end do
          end if

          flux(:,l,m) = flux(:,l,m) + Jacobien_moy * ( metrique_grad(3,l,m) * ev(:) &
                                                     + metrique_grad(4,l,m) * fv(:) )

       end do
    end do

  end subroutine flux_visqueux_y
  
end module second_Visqueux
