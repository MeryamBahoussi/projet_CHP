!!! Routine récupérée du code de Manuel Latige
module premier_euler_direct 
  use m_struct
  use parametres_fluide
  use outils_fluide
  use parametres_globaux, only : DIR_X, DIR_Y
  use premier_visqueux
  !use mod_tapenade

  implicit none 


contains

!!!***************************************************************
!!!*                                                             *
!!!*       calcul du premier-membre non-visqueux du systeme      *
!!!*                                                             *
!!!***************************************************************

  subroutine assemblage_matrice_Euler(sb, b, matrice)

    type (STR_SUPER_BLOC), pointer :: sb
    type (STR_BLOC)      , pointer :: b
    type (STR_MATRICE)   , pointer :: matrice

    integer :: l, m, i
    integer :: ld, lf, md, mf
    real*8 :: theta(b%neq)!, theta_z

    real*8, dimension(:,:,:,:), pointer :: AA, BB, CC, DD, EE

    theta(:) = 1.d0  !! TODO : a mettre en donnees
    
    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    
    AA => matrice%A ; BB => matrice%B; CC => matrice%C 
    DD => matrice%D ; EE => matrice%E

    call calc_derivees_pression(b)

    call efface_mat_euler(b, matrice)
    call gradient_flux_centres(b, AA, BB, CC, DD, EE)

    if (b%maillage%move .and. .not. sb%sol_stat%actif) then 
       !move = .true.
       print*, " euler direct + priseEnCompteInterfaceMobileDansFluide non gere"
       call arret_code
       !    call termes_en_temps_matrice_euler(b, AA, BB, CC, DD, EE) 
    else
       !move = .false.
    end if
    
    call calc_dissipation_implicite(b, DIR_X, sb%solveur_Riemann, AA, BB, CC)
    call calc_dissipation_implicite(b, DIR_Y, sb%solveur_Riemann, DD, BB, EE)
    
!!!!!!! Equation sur Z
    !if ( implicitation_total_advection_z ) then 
    !   call total_advection_z_direction_x(a,b,c)
    !else			
    call imp_advection_z(b, DIR_X, sb%solveur_Riemann, AA, BB, CC)
    call imp_advection_z(b, DIR_Y, sb%solveur_Riemann, DD, BB, EE)
    !end if

!!!!!!! IMPLICITATION condition limites
    call implicitation_CL_direct(b,sb%dt,b%face(1),AA,BB,b%n_dir_l(:,ld-1,md:mf),md,mf)
    call implicitation_CL_direct(b,sb%dt,b%face(2),CC,BB,b%n_dir_l(:,lf,md:mf),md,mf)
    call implicitation_CL_direct(b,sb%dt,b%face(3),DD,BB,b%n_dir_m(:,ld:lf,md-1),ld,lf)
    call implicitation_CL_direct(b,sb%dt,b%face(4),EE,BB,b%n_dir_m(:,ld:lf,mf),ld,lf)
!!!!!!!

!!!!!!! Assemblage
    do m = b%md, b%mf
       do l = b%ld, b%lf 
          do i = 1, b%neq 
             AA(i,:,l,m) = theta(i) * AA(i,:,l,m)
             BB(i,:,l,m) = theta(i) * BB(i,:,l,m)
             CC(i,:,l,m) = theta(i) * CC(i,:,l,m)
             DD(i,:,l,m) = theta(i) * DD(i,:,l,m)
             EE(i,:,l,m) = theta(i) * EE(i,:,l,m)
          end do
       end do
    end do
    
  end subroutine assemblage_matrice_Euler

  subroutine calc_derivees_pression(b)
    use decodage
    use chimie, only : Cv_esp

    type (STR_BLOC), pointer :: b

    integer :: i,l,m
    real*8  :: xi1, xi2, xi, delta1, delta2, rho1c1_2, rho2c2_2, coef_M
    real*8  :: R1, Cv1, dGamma1_drho1z1Ci, d1sXi_drho1z1Ci, rhoein
    real*8  :: gamma1, gamma2, pi1, pi2, Cvi
    type(derivees_partielles), pointer :: dp

    type(STR_fluide), pointer :: f1, f2

    f1 => fluide1; f2 => fluide2
    
    ! en gaz reel, les derivees de p sont calculées dans le decode
    !il faudrait faire ça pour tous les types de fluide
    do m = b%md-1, b%mf+1
       do l = b%ld-1, b%lf+1

          dp => b%deriv_p(l,m)

          gamma1 = b%gamma1(l,m)
          gamma2 = f2%EOS%gamma
          pi1 = f1%EOS%pi
          pi2 = f2%EOS%pi

!!! Calcul des coeffs (voir prop 1 dans la these ou l'article d'Allaire et al 2002) 
          ! xi : inverse du coeff de Gruneisen = 1/(gamma-1)
          xi1 = 1.d0/(gamma1-1.d0)
          xi2 = 1.d0/(gamma2-1.d0)
          xi  = b%z(l,m)*xi1 + (1.d0-b%z(l,m))*xi2
          ! delta_k : d(rho_k eps_k)/d(rho_k) a p constant 
          delta1 = 0.d0
          delta2 = 0.d0          ! a voir si Mie Gruneisen
          ! M = rho1*(delta1-eps1) - rho2*(delta2-eps2)
          ! 
          ! or rho_k*(delta_k-eps_k) = p - xi_k*rho_k*c_k^2
          rho1c1_2 = gamma1*(b%p(l,m) + pi1)
          rho2c2_2 = gamma2*(b%p(l,m) + pi2)
          coef_M = xi2*rho2c2_2 - xi1*rho1c1_2 

!!! derivees partielles de la pression
          dp%d_drhoeps      = 1.d0/xi
          dp%d_drho1z1Ci(1) = -delta1/xi
          dp%d_drho2z2      = -delta2/xi
          dp%d_dz           = coef_M/xi

          ! R_melange = R_melange + cel%C_i(i) * R_gaz / Masse_mol_i(i)
          ! Cv_melange = Cv_melange + cel%C_i(i) * Cv_i(i)
          if (f1%type == CHIMIE_FIGEE) then
             R1 = 0.d0; Cv1 = 0.d0
             do i = 1, b%ne
                R1 = R1 + b%c_esp(i,l,m) * f1%especes(i)%R
                Cvi = Cv_esp(f1%especes(i), -1.d70) !! T ?? 
                Cv1 = Cv1 + b%c_esp(i,l,m) * Cvi
             end do
             rhoein = b%U(irhoE,l,m)-0.5d0*b%rho(l,m)*(b%u_x(l,m)**2.0d0+b%u_y(l,m)**2.d0)
             do i = 1, b%ne
                Cvi = Cv_esp(f1%especes(i), -1.d70) !! T ?? 
                dGamma1_drho1z1Ci = (f1%especes(i)%R*Cv1 - R1*Cvi) / Cv1**2.d0
                d1sXi_drho1z1Ci = dGamma1_drho1z1Ci * b%z(l,m) / (xi*(gamma1-1.d0))**2.d0
                dp%d_drho1z1Ci(i) = ( rhoein - (1.d0-b%z(l,m))*pi2*gamma2*xi2) * d1sXi_drho1z1Ci
             end do
          end if

       end do
    end do
  end subroutine calc_derivees_pression


  !***************************************************************
  !*                                                             *
  !*     calcul des gradients des flux physiques E, F            *
  !*                                                             *
  !*     ON SUPPOSE grade, gradf        INITIALISE a 0.          *
  !*                                                             *
  !***************************************************************
  subroutine gradient_flux_centres(b, AA, BB, CC, DD, EE)
    type (STR_BLOC)                                   , pointer       :: b
    real*8, dimension(b%neq,b%neq,b%ld:b%lf,b%md:b%mf), intent(inout) :: AA,BB,CC,DD,EE

    integer :: l, m, i, j, nb_especes
    real*8  :: u, v, Gv2, htot, y1, y2
    real*8  :: dp_rho1z1Ci(1:b%ne)
    real*8  :: dp_rho2z2, dp_rhou, dp_rhov, dp_rhoe, dp_z
    real*8, dimension(:,:,:)  , pointer :: nl, nm
    real*8, dimension(:,:,:)  , pointer :: ci
    real*8, dimension(:,:,:,:), pointer :: grad_e, grad_f

    type(derivees_partielles), pointer :: dp

    nb_especes = b%ne

    nl => b%n_dir_l
    nm => b%n_dir_m

    grad_e => b%grad_e
    grad_f => b%grad_f

    select case (fluide1%type)
    case(EOS, CHIMIE_FIGEE)
       ci => b%c_esp
    end select
    
    do m = b%md-1, b%mf+1
       do l = b%ld-1, b%lf+1

          u = b%u_x(l,m);  v = b%u_y(l,m) 
          Gv2 = u**2 + v**2
          htot= b%Htot(l,m)

          dp => b%deriv_p(l,m)

          dp_rho1z1Ci(:) = dp%d_drho1z1Ci(:) + dp%d_drhoeps * Gv2 / 2.0d0
          dp_rho2z2      = dp%d_drho2z2      + dp%d_drhoeps * Gv2 / 2.0d0
          dp_rhou        =-dp%d_drhoeps * u
          dp_rhov        =-dp%d_drhoeps * v
          dp_rhoe        = dp%d_drhoeps
          dp_z           = dp%d_dz

          y1 = b%y(l,m); y2 = 1.d0-b%y(l,m)

!!!!!!!!! grad_e
          grad_e(:,:,l,m) = 0.d0
          do i = 1, nb_especes
             grad_e(irho1z1+i-1,:,l,m) = 0.d0 
             do j = 1, nb_especes
                grad_e(irho1z1+i-1,irho1z1+j-1,l,m) =-y1 * u * ci(i,l,m)
             end do
             ! on écrase pour le (i,i)
             grad_e(irho1z1+i-1,irho1z1+i-1,l,m) = u*(y2+y1*(1.d0-ci(i,l,m)))
             grad_e(irho1z1+i-1,irho2z2,l,m) =-y1 * ci(i,l,m) * u
             grad_e(irho1z1+i-1,irhou  ,l,m) = y1 * ci(i,l,m)
          end do

          grad_e(irho2z2,:,l,m) = 0.d0
          do i = 1, nb_especes
             grad_e(irho2z2,irho1z1+i-1,l,m) =-y2 * u
          end do
          grad_e(irho2z2,irho2z2,l,m) = y1 * u
          grad_e(irho2z2,irhou  ,l,m) = y2

          do i = 1, nb_especes
             grad_e(irhou,irho1z1+i-1,l,m) = - u**2 + dp_rho1z1Ci(i)  
          end do
          grad_e(irhou,irho2z2,l,m) =-u**2 + dp_rho2z2 
          grad_e(irhou,irhou  ,l,m) = 2.d0 * u + dp_rhou
          grad_e(irhou,irhov  ,l,m) = dp_rhov
          grad_e(irhou,irhoE  ,l,m) = dp_rhoe
          grad_e(irhou,iz     ,l,m) = dp_z

          do i = 1, nb_especes
             grad_e(irhov,irho1z1+i-1,l,m) =-u * v
          end do
          grad_e(irhov,irho2z2,l,m) =-u * v
          grad_e(irhov,irhou  ,l,m) = v
          grad_e(irhov,irhov  ,l,m) = u
          grad_e(irhov,irhoE  ,l,m) = 0.d0
          grad_e(irhov,iz     ,l,m) = 0.d0

          do i = 1, nb_especes
             grad_e(irhoE,irho1z1+i-1,l,m) =-u * htot + u * dp_rho1z1Ci(i)
          end do
          grad_e(irhoE,irho2z2,l,m) =-u * htot + u * dp_rho2z2
          grad_e(irhoE,irhou  ,l,m) = htot + u * dp_rhou
          grad_e(irhoE,irhov  ,l,m) = u * dp_rhov
          grad_e(irhoE,irhoE  ,l,m) = u * (1.d0 + dp_rhoe)
          grad_e(irhoE,iz     ,l,m) = u * dp_z
    
          grad_e(iz,:,l,m) = 0.d0

!!!!!!!!! grad_f
          grad_f(:,:,l,m) = 0.d0
          do i = 1, nb_especes
             grad_f(irho1z1+i-1,:,l,m) = 0.d0 
             do j = 1, nb_especes
                grad_f(irho1z1+i-1,irho1z1+j-1,l,m) =-y1 * v * ci(i,l,m)
             end do
             ! on écrase pour le (i,i)
             grad_f(irho1z1+i-1,irho1z1+i-1,l,m) = v * (y2+y1*(1.d0-ci(i,l,m)))
             grad_f(irho1z1+i-1,irho2z2,l,m) =-y1 * ci(i,l,m) * v
             grad_f(irho1z1+i-1,irhov  ,l,m) = y1 * ci(i,l,m)
          end do

          grad_f(irho2z2,:,l,m) = 0.d0 
          do i = 1, nb_especes
             grad_f(irho2z2,irho1z1+i-1,l,m) =-y2 * v
          end do
          grad_f(irho2z2,irho2z2,l,m) = y1 * v
          grad_f(irho2z2,irhov  ,l,m) = y2 

          grad_f(irhou,:,l,m) = 0.d0
          do i = 1, nb_especes
             grad_f(irhou,irho1z1+i-1,l,m) =-u * v
          end do
          grad_f(irhou,irho2z2,l,m) =-u * v
          grad_f(irhou,irhou  ,l,m) = v
          grad_f(irhou,irhov  ,l,m) = u        

          do i = 1, nb_especes
             grad_f(irhov,irho1z1+i-1,l,m) =-v**2 + dp_rho1z1Ci(i)
          end do
          grad_f(irhov,irho2z2,l,m) =-v**2 + dp_rho2z2 
          grad_f(irhov,irhou  ,l,m) = dp_rhou
          grad_f(irhov,irhov  ,l,m) = 2.d0 * v + dp_rhov
          grad_f(irhov,irhoE  ,l,m) = dp_rhoe
          grad_f(irhov,iz     ,l,m) = dp_z

          do i = 1, nb_especes
             grad_f(irhoE,irho1z1+i-1,l,m) =-v * htot + v * dp_rho1z1Ci(i)
          end do
          grad_f(irhoE,irho2z2,l,m) =-v * htot + v * dp_rho2z2
          grad_f(irhoE,irhou  ,l,m) = v * dp_rhou
          grad_f(irhoE,irhov  ,l,m) = htot + v * dp_rhov
          grad_f(irhoE,irhoE  ,l,m) = v * (1.d0 + dp_rhoe)
          grad_f(irhoE,iz     ,l,m) = v * dp_z

          grad_f(iz,:,l,m) = 0.d0
       end do
    end do
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Assemblage de la matrice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! direction x
    do m = b%md, b%mf
       do l = b%ld, b%lf
          AA(:,:,l,m) = AA(:,:,l,m) - 0.5d0 * b%dl_l(l-1,m) * (nl(1,l-1,m) * grad_e(:,:,l-1,m) &  
               &                                             + nl(2,l-1,m) * grad_f(:,:,l-1,m)) 
          
          BB(:,:,l,m) = BB(:,:,l,m) - 0.5d0 * b%dl_l(l-1,m) * (nl(1,l-1,m) * grad_e(:,:,l,m)   &  
               &                                             + nl(2,l-1,m) * grad_f(:,:,l,m))  &
               &                    + 0.5d0 * b%dl_l(l  ,m) * (nl(1,l  ,m) * grad_e(:,:,l,m)   &  
               &                                             + nl(2,l  ,m) * grad_f(:,:,l,m))
  
          CC(:,:,l,m) = CC(:,:,l,m) + 0.5d0 * b%dl_l(l  ,m) * (nl(1,l  ,m) * grad_e(:,:,l+1,m) &
               &                                             + nl(2,l  ,m) * grad_f(:,:,l+1,m))
       end do
    end do
    
!!!!! direction y
    do m = b%md, b%mf
       do l = b%ld, b%lf
          DD(:,:,l,m) = DD(:,:,l,m) - 0.5d0 * b%dl_m(l,m-1) * (nm(1,l,m-1) * grad_e(:,:,l,m-1) &
               &                                             + nm(2,l,m-1) * grad_f(:,:,l,m-1))
          
          BB(:,:,l,m) = BB(:,:,l,m) - 0.5d0 * b%dl_m(l,m-1) * (nm(1,l,m-1) * grad_e(:,:,l,m)   &
               &                                             + nm(2,l,m-1) * grad_f(:,:,l,m))  &      
               &                    + 0.5d0 * b%dl_m(l,m  ) * (nm(1,l,m  ) * grad_e(:,:,l,m)   &
               &                                             + nm(2,l,m  ) * grad_f(:,:,l,m)) 
          
          EE(:,:,l,m) = EE(:,:,l,m) + 0.5d0 * b%dl_m(l,m  ) * (nm(1,l,m  ) * grad_e(:,:,l,m+1) &  
               &                                             + nm(2,l,m  ) * grad_f(:,:,l,m+1))
       end do
    end do
  end subroutine gradient_flux_centres

  subroutine calc_dissipation_implicite(b, direction, solveurRiemann, AA, BB, CC)
    type (STR_BLOC)           , pointer       :: b
    integer                   , intent(in)    :: direction, solveurRiemann
    real*8, dimension(:,:,:,:), intent(inout) :: AA, BB, CC

    integer         :: i, l, m, hl, hm, lg, mg, ld, md
    real*8          :: Ug(2), Ud(2), Ung, Und, U_s, Un_roe, U_roe(2)
    real*8          :: xi_1, xi_2, xi_g, xi_d, xi_roe, moyenne_c2, c_roe, coef_g, coef_d
    real*8          :: Lambda1, Lambda2, Lambda3
    real*8          :: epsilon
    real*8, pointer :: disimp(:,:)
    real*8, pointer :: Cm(:,:), Cp(:,:), n(:,:,:), dl(:,:)

    select case(direction)
    case (DIR_X)
       Cm => b%Cm_l
       Cp => b%Cp_l
       n => b%n_dir_l
       dl => b%dl_l
       hl=1; hm=0
       !! AA, BB, CC sont vraiment AA, BB, CC
    case (DIR_Y)
       Cm => b%Cm_m
       Cp => b%Cp_m
       n => b%n_dir_m
       dl => b%dl_m
       hl=0; hm=1
       !! AA, BB, CC sont egaux a DD, BB, EE
    end select

    disimp => b%zeros
    
    select case (solveurRiemann)
    case(GALLICE)
       do m = b%md-hm, b%mf
          do l = b%ld-hl, b%lf

             lg = l   ; mg = m
             ld = l+hl; md = m+hm

             Ug(:) = (/ b%u_x(lg,mg), b%u_y(lg,mg) /)
             Ud(:) = (/ b%u_x(ld,md), b%u_y(ld,md) /)

             Ung = dot_product( Ug, n(:,l,m))
             Und = dot_product( Ud, n(:,l,m))

             U_s = (Cm(l,m)*Ung+Cp(l,m)*Und - (b%p(ld,md)-b%p(lg,mg)))/(Cm(l,m)+Cp(l,m))

             Lambda1 = Ung-Cm(l,m)/b%rho(lg,mg)
             Lambda2 = U_s
             Lambda3 = Und+Cp(l,m)/b%rho(ld,md)

             epsilon = max(abs(Lambda1),abs(Lambda2),abs(Lambda3)) * b%delta_epsilon
             disimp(l,m) = max(val_abs(epsilon,Lambda1),val_abs(epsilon,Lambda2),val_abs(epsilon,Lambda3))
             disimp(l,m) = 0.5d0 * disimp(l,m) * dl(l,m)
          end do
       end do
              
    case(ROE)
       do m = b%md-hm, b%mf
          do l = b%ld-hl, b%lf

             lg = l   ; mg = m
             ld = l+hl; md = m+hm

             Ug(:) = (/ b%u_x(lg,mg), b%u_y(lg,mg) /)
             Ud(:) = (/ b%u_x(ld,md), b%u_y(ld,md) /)

             Ung = dot_product( Ug, n(:,l,m))
             Und = dot_product( Ud, n(:,l,m))

             coef_g =  sqrt(b%rho(lg,mg)) / ( sqrt(b%rho(lg,mg)) + sqrt(b%rho(ld,md)) )
             coef_d =  sqrt(b%rho(ld,md)) / ( sqrt(b%rho(lg,mg)) + sqrt(b%rho(ld,md)) )

             U_roe(1) = Ug(1) * coef_g + Ud(1) * coef_d  
             U_roe(2) = Ug(2) * coef_g + Ud(2) * coef_d 
             Un_roe = Ung * coef_g + Und * coef_d

             xi_1 = 1.d0/(b%gamma1(l,m)-1.d0)
             xi_2 = 1.d0/(fluide2%EOS%gamma-1.d0)
             xi_g = b%z(lg,mg) * xi_1 + (1.d0-b%z(lg,mg)) * xi_2
             xi_d = b%z(ld,md) * xi_1 + (1.d0-b%z(ld,md)) * xi_2

             xi_roe = coef_g * xi_g + coef_d * xi_d
             moyenne_c2 = coef_g * xi_g * b%c(lg,mg)**2 + coef_d * xi_d * b%c(ld,md)**2 
             c_roe = sqrt( moyenne_c2 / xi_roe )

             Lambda1 = Un_roe - c_roe
             Lambda2 = Un_roe
             Lambda3 = Un_roe + c_roe

             epsilon = b%delta_epsilon * Max( abs(Lambda1), abs(Lambda2), abs(Lambda3) )
             disimp(l,m) = max(val_abs(epsilon,Lambda1),val_abs(epsilon,Lambda2),val_abs(epsilon,Lambda3))
             disimp(l,m) = 0.5d0 * disimp(l,m) * dl(l,m)
          end do
       end do

    case (HLL)
       do m = b%md-hm, b%mf
          do l = b%ld-hl, b%lf
             lg = l   ; mg = m
             ld = l+hl; md = m+hm

             Ug(:) = (/ b%u_x(lg,mg), b%u_y(lg,mg) /)
             Ud(:) = (/ b%u_x(ld,md), b%u_y(ld,md) /)
             
             Ung = dot_product( Ug, n(:,l,m))
             Und = dot_product( Ud, n(:,l,m))
             
             disimp(l,m) = max(0.d0, Ung+b%c(lg,mg), Und+b%c(ld,md))
             disimp(l,m) = 0.5d0 * disimp(l,m) * dl(l,m)
          end do
       end do

    case default
       print*, "ERREUR calc_dissipation_implicite : solveur de Riemann inconnu"
       call arret_code
    end select

    do m = b%md-hm, b%mf-hm
       do l = b%ld-hl, b%lf-hl
          do i = 1, b%neq-1 
             AA(i,i,l+hl,m+hm) = AA(i,i,l+hl,m+hm) - disimp(l,m)
             BB(i,i,l+hl,m+hm) = BB(i,i,l+hl,m+hm) + disimp(l,m)
          end do
       end do
    end do
    
    do m = b%md, b%mf
       do l = b%ld, b%lf
          do i = 1, b%neq-1 
             BB(i,i,l,m) = BB(i,i,l,m) + disimp(l,m) 
             CC(i,i,l,m) = CC(i,i,l,m) - disimp(l,m)
          end do
       end do
    end do
    
    b%zeros(:,:) = 0.d0
    
  end subroutine calc_dissipation_implicite

  !   subroutine termes_en_temps_matrice_euler(b, AA, BB, CC, DD, EE)
  !     type (STR_BLOC)                                   , pointer       :: b
  !     real*8, dimension(b%neq,b%neq,b%ld:b%lf,b%md:b%mf), intent(inout) :: AA,BB,CC,DD,EE

  !     integer :: l, m, i

  !     do m = b%md, b%mf
  !        do l = b%ld, b%lf
  !           do i = 1, b%neq !-1 (??)
  !              AA(i,i,l,m) = AA(i,i,l,m) - 0.5d0 * b%dl_l(l-1,m) * b%dX_dt(l-1,m) 
  !              CC(i,i,l,m) = CC(i,i,l,m) + 0.5d0 * b%dl_l(l  ,m) * b%dX_dt(l,m)
  !              BB(i,i,l,m) = BB(i,i,l,m) + 0.5d0 * b%dl_l(l-1,m) * b%dX_dt(l-1,m) &
  !                   &                    - 0.5d0 * b%dl_l(l  ,m) * b%dX_dt(l,m)   &
  !                   &                    + 0.5d0 * b%dl_m(l,m-1) * b%dY_dt(l,m-1) &      
  !                   &                    - 0.5d0 * b%dl_m(l,m)   * b%dY_dt(l,m)
  !              DD(i,i,l,m) = DD(i,i,l,m) - 0.5d0 * b%dl_m(l,m-1) * b%dY_dt(l,m-1)
  !              EE(i,i,l,m) = EE(i,i,l,m) + 0.5d0 * b%dl_m(l,m)   * b%dY_dt(l,m)
  !           end do
  !        end do
  !     end do

  !   end subroutine termes_en_temps_matrice_euler


  !***************************************************************
  !*                                                             *
  !*     Implicitation de l'equation sur z                       *
  !*                                                             *
  !***************************************************************	
  subroutine imp_advection_z(b, direction, solveurRiemann, AA, BB, CC)
    type (STR_BLOC)           , pointer       :: b
    integer                   , intent(in)    :: direction, solveurRiemann
    real*8, dimension(:,:,:,:), intent(inout) :: AA, BB, CC

    integer         :: l, m, hl, hm, lg, mg, ld, md
    real*8          :: U(2), Uim1, Ui, Uip1, Us_m, Us_p
    real*8          :: delta_p
    real*8, pointer :: Cm(:,:), Cp(:,:), n(:,:,:), dl(:,:)

    select case(direction)
    case (DIR_X)
       Cm => b%Cm_l
       Cp => b%Cp_l
       n => b%n_dir_l
       dl => b%dl_l
       hl=1; hm=0
       !! AA, BB, CC sont vraiment AA, BB, CC
    case (DIR_Y)
       Cm => b%Cm_m
       Cp => b%Cp_m
       n => b%n_dir_m
       dl => b%dl_m
       hl=0; hm=1
       !! AA, BB, CC sont egaux a DD, BB, EE
    end select

    select case(solveurRiemann)
    case( GALLICE )
       do m = b%md, b%mf
          do l = b%ld,b%lf

             lg = l-hl; mg = m-hm
             ld = l+hl; md = m+hm

             U(:) = (/ b%u_x(lg,mg), b%u_y(lg,mg) /)
             Uim1 = dot_product( n(:,lg,mg), U )
             U(:) = (/ b%u_x(l ,m ), b%u_y(l ,m ) /)
             Ui   = dot_product( n(:,lg,mg), U )
             delta_p = b%p(l,m) - b%p(lg,mg)

             Us_m = (Cm(lg,mg)*Uim1 + Cp(lg,mg)*Ui - delta_p) / (Cm(lg,mg)+Cp(lg,mg))

             U(:) = (/ b%u_x(l ,m ), b%u_y(l ,m ) /)
             Ui   = dot_product( n(:,l ,m ), U )
             U(:) = (/ b%u_x(ld,md), b%u_y(ld,md) /)
             Uip1 = dot_product( n(:,l ,m ), U )
             delta_p = b%p(ld,md) - b%p(l,m)

             Us_p = ( Cm(l,m)*Ui + Cp(l,m)*Uip1 - delta_p ) / (Cm(l,m) + Cp(l,m))

             AA(iz,iz,l,m) = AA(iz,iz,l,m) - 0.5d0*(Us_m+abs(Us_m)) * dl(lg,mg)
             BB(iz,iz,l,m) = BB(iz,iz,l,m) + 0.5d0*(Us_m+abs(Us_m)) * dl(lg,mg)
             BB(iz,iz,l,m) = BB(iz,iz,l,m) - 0.5d0*(Us_p-abs(Us_p)) * dl(l,m)
             CC(iz,iz,l,m) = CC(iz,iz,l,m) + 0.5d0*(Us_p-abs(Us_p)) * dl(l,m)
          end do
       end do

    case(ROE, HLL)

       do m = b%md, b%mf
          do l = b%ld,b%lf  
             lg = l-hl; mg = m-hm
             ld = l+hl; md = m+hm

             U(:) = (/ b%u_x(lg,mg), b%u_y(lg,mg) /)
             Uim1 = dot_product( n(:,lg,mg), U )
             U(:) = (/ b%u_x(l ,m ), b%u_y(l ,m ) /)
             Ui   = dot_product( n(:,lg,mg), U )
             ! moyenne de roe
             Us_m = (sqrt(b%rho(l,m))*Ui + sqrt(b%rho(lg,mg))*Uim1) / (sqrt(b%rho(l,m))+sqrt(b%rho(lg,mg)))

             U(:) = (/ b%u_x(l ,m ), b%u_y(l ,m ) /)
             Ui   = dot_product( n(:,l ,m ), U )
             U(:) = (/ b%u_x(ld,md), b%u_y(ld,md) /)
             Uip1 = dot_product( n(:,l ,m ), U )
             ! moyenne de roe
             Us_p = (sqrt(b%rho(l,m))*Ui + sqrt(b%rho(ld,md))*Uip1) / (sqrt(b%rho(l,m))+sqrt(b%rho(ld,md)))

             AA(iz,iz,l,m) = AA(iz,iz,l,m) - 0.5d0*(Us_m+abs(Us_m)) * dl(lg,mg)
             BB(iz,iz,l,m) = BB(iz,iz,l,m) + 0.5d0*(Us_m+abs(Us_m)) * dl(lg,mg)
             BB(iz,iz,l,m) = BB(iz,iz,l,m) - 0.5d0*(Us_p-abs(Us_p)) * dl(l,m)
             CC(iz,iz,l,m) = CC(iz,iz,l,m) + 0.5d0*(Us_p-abs(Us_p)) * dl(l,m)
          end do
       end do
    case default
       print*, "ERREUR imp_advection_z : solveur de Riemann inconnu"
       call arret_code
    end select

  end subroutine imp_advection_z

  subroutine implicitation_CL_direct(b,dt,face,aa,bb,n,taille_inf,taille_sup)
    use outils_conditions_limites

    type (STR_BLOC)           , intent(in)    :: b
    real*8                    , intent(in)    :: dt
    type(STR_FACE)            , intent(in)    :: face
    integer                   , intent(in)    :: taille_inf,taille_sup
    real*8, dimension(:,:,:,:), intent(in)    :: AA
    real*8, dimension(:,:,:,:), intent(inout) :: BB
    real*8                    , intent(in)    :: n(2,taille_inf:taille_sup)

    real*8  :: normale(2)
    real*8  :: mat(b%neq,b%neq)
    integer :: id, ifin, jd, jf
    integer :: i, l, m
    integer :: h_l,h_m
    integer :: ip
    type (STR_PATCH), pointer :: patch
   ! real*8  :: a, fac

    do ip = 1, face%nb_patch
       patch => face%patch(ip)

       select case (patch%bc%position)
       case(4) !North
          id = patch%ideb; jd = b%mf
          ifin = patch%ifin; jf = b%mf
          h_l = 1; h_m = 0
       case(3) !South
          id = patch%ideb; jd = b%md
          ifin = patch%ifin; jf = b%md
          h_l = 1; h_m = 0     
       case(2) !East
          id = b%lf; jd = patch%ideb
          ifin = b%lf; jf = patch%ifin
          h_l = 0; h_m = 1     
       case(1) !West
          id = b%ld; jd = patch%ideb
          ifin = b%ld; jf = patch%ifin  
          h_l = 0; h_m = 1       
       end select

       select case (trim(patch%bc%condition))
       case("interproc", 'interbloc')
          ! on fait la resolution avec MKL sur le proc 0 uniquement
          ! inversion en parallele avec Petsc 
       case("entree_supersonique")
          ! tout est au second membre 
          ! ou on ne fait rien si on resoud en delta forme

          ! Pour cette condition il n'y a strictement rien a faire dans la matrice

          ! Le newton s'ecrit (k indice d'itération du newton)
          ! delta_G/delta_U * (U_k+1-U_k) = - G(U_k)
          ! avec la matrice delta_G/delta_U = (AA,BB,CC)
          ! En detaillant les indices en espace, 
          ! la ligne de la matrice multipliee par le vecteur donne (U_k+1-U_k): 
          ! AA_i * (U_k+1,i-1 - U_k,i-1) + 
          ! BB_i * (U_k+1,i   - U_k,i  ) + 
          ! CC_i * (U_k+1,i+1 - U_k,i+1) = - G(U_k,i)

          ! La condition de dirichlet ne change pas au cours des itérations en k, donc
          ! on a  U_k+1,i-1 = U_k,i-1
          ! or pour les points du bord le terme AA_i n'est pas dans la matrice globale
          ! donc en ne faisant rien cela revient à faire 
          ! (U_k+1,i-1 - U_k,i-1) = 0

       case("paroi_global", "symetrie", "paroi_global_adiabatique", "paroi_flux_impose")
          do m = jd, jf
             do l = id, ifin
                mat = 0.d0

                normale(1:2) = n(1:2,l*h_l+m*h_m)

                mat(irhou,irhou) = -2.0d0*normale(1)**2
                mat(irhov,irhou) = -2.0d0*normale(1)*normale(2)
                mat(irhou,irhov) = mat(irhov,irhou)
                mat(irhov,irhov) = -2.0d0*normale(2)**2

                do i = 1, b%neq
                   mat(i,i) = mat(i,i) + 1.0d0
                end do

                BB(:,:,l,m) = BB(:,:,l,m) + matmul(AA(:,:,l,m),mat(:,:))

             end do
          end do

       case("injection_carbone") 
          !! rien n'est fait dans le code de manu
!!!   paroi, il y a injection donc il y a un dirichlet uniquement sur le z
!!$          do m = jd, jf
!!$             do l = id, ifin
!!$                mat = 0.d0
!!$                
!!$                normale(1:2) = n(1:2,l*h_l+m*h_m)
!!$                
!!$                mat(irhou,irhou) = -2.0d0*normale(1)**2
!!$                mat(irhov,irhou) = -2.0d0*normale(1)*normale(2)
!!$                mat(irhou,irhov) = mat(irhov,irhou)
!!$                mat(irhov,irhov) = -2.0d0*normale(2)**2
!!$                
!!$                do i = 1, b%neq
!!$                   mat(i,i) = mat(i,i) + 1.0d0
!!$                end do
!!$                mat(iz,iz) = 0.d0  ! dirichlet sur z
!!$                
!!$                BB(:,:,l,m) = BB(:,:,l,m) + matmul(AA(:,:,l,m),mat(:,:))
!!$             end do
!!$          end do
          
       case("flux_nul")
          do m = jd, jf
             do l = id, ifin
                BB(:,:,l,m) = BB(:,:,l,m) + AA(:,:,l,m)
             end do
          end do

       case default
          print*, "Condition inconnue pour l'implicitation de la matrice EULER", patch%bc%condition
          call arret_code
       end select
    end do
  end subroutine implicitation_CL_direct
  
end module premier_euler_direct
