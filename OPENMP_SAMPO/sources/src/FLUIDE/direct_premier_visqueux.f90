!!! Routine récupérée du code de Manuel Latige
module premier_visqueux_direct 
  use m_struct
  use parametres_globaux, only : DIR_X, DIR_Y
  use decodage
  use a_mettre_en_donnees
  use outils_fluide, only : moyenne
  use premier_visqueux
  
!  use variables_debug, only : cmp_manu
  implicit none 

contains

  subroutine assemblage_matrice_visqueux(sb, b, matrice)

    type (STR_SUPER_BLOC), pointer :: sb
    type (STR_BLOC)      , pointer :: b
    type (STR_MATRICE)   , pointer :: matrice

    integer :: l, m 
    real*8, dimension(:,:,:,:), pointer :: AA, BB, CC, DD, EE
    real*8, dimension(:,:,:,:), pointer :: AD, CD, AE, CE

    AE => matrice%AE; DD => matrice%D; CE => matrice%CE
    AA => matrice%A ; BB => matrice%B; CC => matrice%C 
    AD => matrice%AD; EE => matrice%E; CD => matrice%CD

    !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
    select case(fluide1%type)
    case(EOS)
       call calc_derivees_temperature(b)
    case( CHIMIE_FIGEE )
       call calc_derivee_temperature_melange(b)
    case default
       print*, "derivees temperature non implementees pour ce type de fluide en monophasique"
       call arret_code
    end select
    
    call calc_dv_du(b, b%U, b%dV_dU)

    !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
    call efface_mat_euler(b, matrice)
    
    call grad_visqueux_x(b, AA, BB, CC, DD, EE, AD, CD, AE, CE)
    call grad_visqueux_y(b, AA, BB, CC, DD, EE, AD, CD, AE, CE)
    
    call implicitation_CL_visqueux(b, matrice)

    do m = b%md, b%mf
       do l = b%ld, b%lf 
          AA(:,:,l,m) = matmul(AA(:,:,l,m), b%dV_dU(:,:,l-1,m  ))
          BB(:,:,l,m) = matmul(BB(:,:,l,m), b%dV_dU(:,:,l  ,m  ))
          CC(:,:,l,m) = matmul(CC(:,:,l,m), b%dV_dU(:,:,l+1,m  ))
          DD(:,:,l,m) = matmul(DD(:,:,l,m), b%dV_dU(:,:,l  ,m-1))
          EE(:,:,l,m) = matmul(EE(:,:,l,m), b%dV_dU(:,:,l  ,m+1))
          
          AD(:,:,l,m) = matmul(AD(:,:,l,m), b%dV_dU(:,:,l-1,m-1))
          CD(:,:,l,m) = matmul(CD(:,:,l,m), b%dV_dU(:,:,l+1,m-1))
          AE(:,:,l,m) = matmul(AE(:,:,l,m), b%dV_dU(:,:,l-1,m+1))
          CE(:,:,l,m) = matmul(CE(:,:,l,m), b%dV_dU(:,:,l+1,m+1))
       end do
    end do
    
  end subroutine assemblage_matrice_visqueux

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!       implicitation des conditions limites     !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine implicitation_CL_visqueux(b, matrice)
    type (STR_BLOC)   , pointer :: b
    type (STR_MATRICE), pointer :: matrice
    
    integer :: ld, lf, md, mf
    real*8, dimension(:,:,:,:), pointer :: AA, BB, CC, DD, EE
    real*8, dimension(:,:,:,:), pointer :: AD, CD, AE, CE
    
    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    AE => matrice%AE; DD => matrice%D; CE => matrice%CE
    AA => matrice%A ; BB => matrice%B; CC => matrice%C 
    AD => matrice%AD; EE => matrice%E; CD => matrice%CD
    
!!!!! implicitation de la matrice
!!! face 1
    call imp_CL_face_visqueux(b%face(1), b%neq, md, mf, b%n_dir_l(:,ld-1,:), DD(:,:,ld,:), AD(:,:,ld,:), BB(:,:,ld,:), AA(:,:,ld,:), EE(:,:,ld,:), AE(:,:,ld,:))
!!! face 2 
    call imp_CL_face_visqueux(b%face(2), b%neq, md, mf, b%n_dir_l(:,lf,:)  , DD(:,:,lf,:), CD(:,:,lf,:), BB(:,:,lf,:), CC(:,:,lf,:), EE(:,:,lf,:), CE(:,:,lf,:))
!!! face 3
    call imp_CL_face_visqueux(b%face(3), b%neq, ld, lf, b%n_dir_m(:,:,md-1), AA(:,:,:,md), AD(:,:,:,md), BB(:,:,:,md), DD(:,:,:,md), CC(:,:,:,md), CD(:,:,:,md))
!!! face 4
    call imp_CL_face_visqueux(b%face(4), b%neq, ld, lf, b%n_dir_m(:,:,mf)  , AA(:,:,:,mf), AE(:,:,:,mf), BB(:,:,:,mf), EE(:,:,:,mf), CC(:,:,:,mf), CE(:,:,:,mf))
    
    
    if (matrice%stencil == 9) then
!!!------------------
!!! coins 
!!!------------------
!!! coin entre la face 1 et 3
       call imp_CL_coin_visqueux(b%coin(1), b%neq, BB(:,:,ld,md), DD(:,:,ld,md), AA(:,:,ld,md), AD(:,:,ld,md))
!!! coin entre la face 3 et 2
       call imp_CL_coin_visqueux(b%coin(2), b%neq, BB(:,:,lf,md), DD(:,:,lf,md), CC(:,:,lf,md), CD(:,:,lf,md))
!!! coin entre la face 2 et 4
       call imp_CL_coin_visqueux(b%coin(3), b%neq, BB(:,:,lf,mf), EE(:,:,lf,mf), CC(:,:,lf,mf), CE(:,:,lf,mf))
!!! coin entre la face 4 et 1 
       call imp_CL_coin_visqueux(b%coin(4), b%neq, BB(:,:,ld,mf), EE(:,:,ld,mf), AA(:,:,ld,mf), AE(:,:,ld,mf))
    end if

  end subroutine implicitation_CL_visqueux
  
!!!!!  AE | EE   CE
!!!!!  AA | BB   CC
!!!!!  AD | DD   CD
  subroutine imp_CL_face_visqueux(face, dim_bloc, kd, kf, n, DD, AD, BB, AA, EE, AE)
    use outils_fluide

    type (STR_FACE)                           , intent(in)    :: face
    integer                                   , intent(in)    :: kd, kf, dim_bloc
    real*8, dimension(2,kd:kf)                , intent(in)    :: n  ! normale
    real*8, dimension(dim_bloc,dim_bloc,kd:kf), intent(inout) :: DD, AD, BB, AA, EE, AE
    
    integer                   :: ip, k, i
    real*8                    :: mat(dim_bloc,dim_bloc)
    type (STR_PATCH), pointer :: patch
    
    
    if (face%nb_patch>1) then
       print*, "implicitation_matrice_face_visqueux non gere pour nb_patch>1"
       call arret_code
       ! il faut changer les kd, kf en patch%ideb, patch%ifin
       ! pb (?) comment gerer les matrices de coins au bord des patchs ?
    end if
    do ip = 1, face%nb_patch
       patch => face%patch(ip)
       select case (trim(patch%bc%condition))
       case("interproc", 'interbloc')
          ! on fait la resolution avec MKL sur le proc 0 uniquement
          ! inversion en parallele avec Petsc
       case("entree_supersonique") 
          ! rien a faire, tout est au second membre  (? ne donne pas AX=Sm !?)
       case("flux_nul")
!!$          if (cmp_manu) then
!!$             do k = kd, kf
!!$                mat = 0.d0
!!$                do i = irhou, irhoE 
!!$                   mat(i,i) = 1.d0  !!! manu mais pourquoi ??
!!$                end do
!!$                if (k/=kd) DD(:,:,k) = DD(:,:,k) + matmul(AD(:,:,k),mat(:,:))
!!$                           BB(:,:,k) = BB(:,:,k) + matmul(AA(:,:,k),mat(:,:))
!!$                if (k/=kf) EE(:,:,k) = EE(:,:,k) + matmul(AE(:,:,k),mat(:,:))
!!$             end do
!!$          else
             do k = kd, kf
                if (k/=kd) DD(:,:,k) = DD(:,:,k) + AD(:,:,k)
                           BB(:,:,k) = BB(:,:,k) + AA(:,:,k)
                if (k/=kf) EE(:,:,k) = EE(:,:,k) + AE(:,:,k)
             end do
!!$          end if
       case("symetrie", "paroi_global_adiabatique", "paroi_flux_impose")
          do k = kd, kf
             mat = 0.d0
             mat(irhou:irhov,irhou:irhov) = - 2.d0*produit_tensoriel(2,n(:,k),n(:,k))
             do i = 1, dim_bloc 
                mat(i,i) = mat(i,i) + 1.d0
             end do
             
             if (k/=kd) DD(:,:,k) = DD(:,:,k) + matmul(AD(:,:,k),mat(:,:))
                        BB(:,:,k) = BB(:,:,k) + matmul(AA(:,:,k),mat(:,:))
             if (k/=kf) EE(:,:,k) = EE(:,:,k) + matmul(AE(:,:,k),mat(:,:))
          end do
          
       case("injection_carbone") 
          ! rien n'est fait dans le code de manu

       case("paroi_global")
          do k = kd, kf
             mat = 0.d0
             do i = 1, dim_bloc 
                mat(i,i) = mat(i,i) + 1.d0
             end do
             mat(irhou,irhou) = -1.d0
             mat(irhov,irhov) = -1.d0
           !  if (cmp_manu) then
           !     mat(irhoE,irhoE) = 0.d0
           !   else
                mat(irhoE,irhoE) = -1.d0
           !  end if

             if (k/=kd) DD(:,:,k) = DD(:,:,k) + matmul(AD(:,:,k),mat(:,:))
                        BB(:,:,k) = BB(:,:,k) + matmul(AA(:,:,k),mat(:,:))
             if (k/=kf) EE(:,:,k) = EE(:,:,k) + matmul(AE(:,:,k),mat(:,:))
          end do
          
       case default
          print*, "Condition inconnue pour l'implicitation de la matrice Visqueuse", patch%bc%condition
          call arret_code
       end select
    end do
  end subroutine imp_CL_face_visqueux

!!!!!  AE | EE   CE
!!!!!     | 
!!!!!  AA | BB   CC
!!!!!     |---------
!!!!!  AD   DD   CD
  subroutine imp_CL_coin_visqueux(coin, dim_bloc, BB, DD, AA, AD)
    use outils_fluide
    
    type (STR_COIN)                     , intent(in)    :: coin
    integer                             , intent(in)    :: dim_bloc
    real*8, dimension(dim_bloc,dim_bloc), intent(inout) :: BB, DD, AA, AD
    
    integer :: i
    real*8  :: n1(2), n2(2)
    real*8  :: mat(dim_bloc, dim_bloc)
  
    if (coin%ib_opp == -1) then
       n1 = coin%n1(:)
       n2 = coin%n2(:)
       
       select case (coin%condition)
       case(CL_INTERPROCvs, CL_INTERBLOCvsINTERBLOC)
          !
!!! coins lies au decoupage MPI
       case(CL_PAROIvsINTERBLOC)
          mat = 0.d0
          mat(irhou:irhov,irhou:irhov) = - 2.d0*produit_tensoriel(2,n1(:),n1(:))
          do i = irhou, irhov
             mat(i,i) = mat(i,i) + 1.d0
          end do
          mat(irhoE,irhoE) = -1.d0 
          AA(:,:) = AA(:,:) + matmul(AD(:,:),mat(:,:))
       
       case(CL_PAROI_ADIABvsINTERBLOC, CL_SYMETRIEvsINTERBLOC)
          mat = 0.d0
          mat(irhou:irhov,irhou:irhov) = - 2.d0*produit_tensoriel(2,n1(:),n1(:))
          do i = irhou, irhov
             mat(i,i) = mat(i,i) + 1.d0
          end do
          AA(:,:) = AA(:,:) + matmul(AD(:,:),mat(:,:))
       
       case(CL_INTERBLOCvsPAROI)
          mat = 0.d0
          mat(irhou:irhov,irhou:irhov) = - 2.d0*produit_tensoriel(2,n2(:),n2(:))
          do i = irhou, irhov
             mat(i,i) = mat(i,i) + 1.d0
          end do
          mat(irhoE,irhoE) = -1.d0 
          DD(:,:) = DD(:,:) + matmul(AD(:,:),mat(:,:))
       
       case(CL_INTERBLOCvsPAROI_ADIAB, CL_INTERBLOCvsSYMETRIE)
          mat = 0.d0
          mat(irhou:irhov,irhou:irhov) = - 2.d0*produit_tensoriel(2,n2(:),n2(:))
          do i = irhou, irhov
             mat(i,i) = mat(i,i) + 1.d0
          end do
          DD(:,:) = DD(:,:) + matmul(AD(:,:),mat(:,:))
          
       case(CL_FLUX_NULvsINTERBLOC)
          AA(:,:) = AA(:,:) + AD(:,:)
       case(CL_INTERBLOCvsFLUX_NUL)
          DD(:,:) = DD(:,:) + AD(:,:)
       
       case(CL_ENTREE_SUPERvsINTERBLOC  , CL_INTERBLOCvsENTREE_SUPER)
          ! rien a faire
          
!!! vrais coins du domaine de calcul 
       case(CL_PAROIvsPAROI, CL_PAROI_ADIABvsPAROI_ADIAB, &
            CL_PAROIvsPAROI_ADIAB, CL_PAROI_ADIABvsPAROI )
          BB(:,irhou:irhov) = BB(:,irhou:irhov) - AD(:,irhou:irhov)
          BB(:,irhoE)       = BB(:,irhoE)       + AD(:,irhoE)
      
       case(CL_FLUX_NULvsFLUX_NUL)
          BB(:,:) = BB(:,:) + AD(:,:)
          
       case(CL_PAROIvsFLUX_NUL, CL_PAROIvsSYMETRIE, &
            CL_FLUX_NULvsPAROI, CL_SYMETRIEvsPAROI)
          BB(:,irhou:irhov) = BB(:,irhou:irhov) - AD(:,irhou:irhov)
          BB(:,irhoE)       = BB(:,irhoE)       - AD(:,irhoE)
          
       case(CL_PAROI_ADIABvsFLUX_NUL, CL_PAROI_ADIABvsSYMETRIE, &
            CL_FLUX_NULvsPAROI_ADIAB, CL_SYMETRIEvsPAROI_ADIAB)
          BB(:,irhou:irhov) = BB(:,irhou:irhov) - AD(:,irhou:irhov)
          BB(:,irhoE)       = BB(:,irhoE)       + AD(:,irhoE)
          
       case(CL_ENTREE_SUPERvsFLUX_NUL   , CL_FLUX_NULvsENTREE_SUPER, &
            CL_ENTREE_SUPERvsSYMETRIE   , CL_SYMETRIEvsENTREE_SUPER, &
            CL_ENTREE_SUPERvsPAROI_ADIAB, CL_PAROI_ADIABvsENTREE_SUPER)
          ! rien a faire
          
       case(CL_ENTREE_SUPERvsPAROI, CL_PAROIvsENTREE_SUPER)
          BB(:,irhoE)       = BB(:,irhoE)       - AD(:,irhoE)
       
       case(CL_PAROIvsINJECTION, CL_INJECTIONvsPAROI, &
            CL_PAROI_ADIABvsINJECTION, CL_INJECTIONvsPAROI_ADIAB)
          BB(:,irhou:irhov) = BB(:,irhou:irhov) - AD(:,irhou:irhov)
          BB(:,irhoE)       = BB(:,irhoE)       - AD(:,irhoE)

       case default
          print*, "implicitation_matrice_coin_visqueux : CL non gerees ", coin%condition
         ! call arret_code
       end select
       
    end if
     
  end subroutine imp_CL_coin_visqueux

  
  subroutine calc_derivees_temperature(b)
    type (STR_BLOC), pointer :: b

    integer :: i,l,m
    real*8  :: rho1z1, rho1, z
    real*8  :: rho2z2, rho2, rhoCv
    real*8  :: G1, drho1_rho1z1, drho1_z, dG1_rho1
    real*8  :: G2, drho2_rho2z2, drho2_z, dG2_rho2
    type(derivees_partielles), pointer :: dT
    type(STR_fluide)         , pointer :: f1, f2
    
    f1 => fluide1; f2 => fluide2
        
    do m = b%md-1, b%mf+1
       do l = b%ld-1, b%lf+1
          
          dT => b%deriv_T(l,m)

          z = b%U(iz,l,m)
          rho1z1 = b%U(irho1z1,l,m)
          rho2z2 = b%U(irho2z2,l,m)

          G1 = 0.d0; drho1_rho1z1 = 0.d0; drho1_z = 0.d0; dG1_rho1 = 0.d0;
          G2 = 0.d0; drho2_rho2z2 = 0.d0; drho2_z = 0.d0; dG2_rho2 = 0.d0;
          if (abs(z) > tol_z .and. abs(rho1z1) > tol_z) then
             rho1 = rho1z1 / b%z(l,m)
             G1 = EOS_eps_from_rho_T(f1%EOS, rho1, 0.d0)
             
             drho1_rho1z1 = 1.d0 / z
             drho1_z      =-rho1 / z
             
             dG1_rho1 =-f1%EOS%pi / rho1**2.d0
          end if
          if (abs(1.d0-z) > tol_z .and. abs(rho2z2) > tol_z) then
             rho2 = rho2z2 / (1.d0-z)
             G2 = EOS_eps_from_rho_T(f2%EOS, rho2, 0.d0)

             drho2_rho2z2 = 1.d0 / (1.d0-z)
             drho2_z      = rho2 / (1.d0-z)

             dG2_rho2 =-f2%EOS%pi / rho2**2.d0  !! pas bon pour MG
          end if

          rhoCv = b%rho(l,m)*b%Cv(l,m)

          do i = 1, b%ne
             dT%d_drho1z1Ci(i) =-( G1+rho1z1*drho1_rho1z1*dG1_rho1+b%T(l,m)*f1%EOS%Cv ) / rhoCv ! ou rho*rho*cv ??
             ! G1+rho1z1*drho1_rho1z1*dG1_rho1 = 0
          end do
          dT%d_drho2z2 =-( G2+rho2z2*drho2_rho2z2*dG2_rho2+b%T(l,m)*f2%EOS%Cv ) / rhoCv 
          dT%d_drhoeps = 1.d0 / rhoCv 
          dT%d_dz      =-(rho1z1*drho1_z*dG1_rho1+rho2z2*drho2_z*dG2_rho2) / rhoCv

       end do
    end do
  
  end subroutine calc_derivees_temperature

  subroutine calc_derivee_temperature_melange(b)
    type (STR_BLOC), pointer :: b

    integer :: l,m
    real*8  :: z, rhoCv
    real*8  :: xi1, omega1
    real*8  :: xi2, omega2
    real*8  :: xi, gamma_tot, omega, dgamma_dz, dPi_dz, dCv_dz
    type(derivees_partielles), pointer :: dT
    type(STR_fluide)         , pointer :: f1, f2
     
    f1 => fluide1; f2 => fluide2
        
    do m = b%md-1, b%mf+1
       do l = b%ld-1, b%lf+1
          
          dT => b%deriv_T(l,m)
          z = b%U(iz,l,m)
          rhoCv = b%rho(l,m)*b%Cv(l,m)
          
          dT%d_drho1z1Ci(:) =-b%T(l,m) / b%rho(l,m)  ! rho et pas rhoCv ??
          dT%d_drho2z2      =-b%T(l,m) / b%rho(l,m) 
          dT%d_drhoeps      = 1.d0 / rhoCv 

          xi1 = 1.d0/(b%gamma1(l,m)-1.d0)
          xi2 = 1.d0/(f2%EOS%gamma-1.d0)
          xi  = z*xi1 + (1.d0-z)*xi2
          gamma_tot = 1.d0 / xi + 1.d0
          
          omega1 = f1%EOS%pi*b%gamma1(l,m)*xi1
          omega2 = f2%EOS%pi*f2%EOS%gamma*xi2
          omega  = z*omega1 + (1.d0-z)*omega2
          
          dgamma_dz = -(xi1-xi2) / xi**2

          dPi_dz = (omega1-omega2)*(1.d0-1.d0/gamma_tot)
          dPi_dz = dPi_dz + omega * dgamma_dz / gamma_tot**2

          if ( b%gamma1(l,m) /= f2%EOS%gamma ) then
             dCv_dz = -(f2%EOS%Cv - b%Cv1(l,m))*dgamma_dz*b%Cv(l,m) &
                  &/ (f2%EOS%Cv*(gamma_tot-f2%EOS%gamma) + b%Cv1(l,m)*(b%gamma1(l,m)-gamma_tot))
          else
             dCv_dz = 0.d0
          end if

          dT%d_dz = - dPi_dz/rhoCv  - b%T(l,m) / b%Cv(l,m) * dCv_dz
       end do
    end do

  end subroutine calc_derivee_temperature_melange

  
!!! En monophasique
!!! U = (rho ci, rho u, rho E) := (u_i, u_2, u_3)
!!! V = (ci, u, T) = (u_i/sum u_i, u_2/sum u_i, T_EOS(sum u_i, u_3 - 1/2 u_2^2/sum u_i))
!!!
!!! On exprime tout en fonction de d T_EOS/d rho_i et d T_EOS/d (rho*eps)
  subroutine calc_dV_dU(b, U, dV_dU)
    type (STR_BLOC)           , pointer :: b
    real*8, dimension(:,:,:)  , pointer :: U
    real*8, dimension(:,:,:,:), pointer :: dV_dU
    
    integer :: i, l,m, nb_especes
    real*8  :: rho, vol, u_x, u_y, u2, v2
    real*8  :: dT_rho2z2, dT_rhou, dT_rhov, dT_rhoE, dT_z
    real*8  :: dT_rho1z1(1:b%ne)
    type(derivees_partielles), pointer :: dT

    nb_especes = b%ne
    
    do m = b%md-1, b%mf+1
       do l = b%ld-1, b%lf+1
          
          dT => b%deriv_T(l,m)
          
          call rho_from_U(U(:,l,m), rho) ! densite
          vol = 1.d0/rho
          
          u_x = U(irhou,l,m)*vol
          u_y = U(irhov,l,m)*vol
          u2 = u_x**2.d0; v2 = u_y**2
          
          dT_rho1z1(1:nb_especes)= dT%d_drho1z1Ci(1:nb_especes) + .5d0*(u2+v2) * dT%d_drhoeps
          dT_rho2z2 = dT%d_drho2z2 + .5d0*(u2+v2) * dT%d_drhoeps
          dT_rhou   =-u_x * dT%d_drhoeps
          dT_rhov   =-u_y * dT%d_drhoeps
          dT_rhoE   = dT%d_drhoeps
          dT_z      = dT%d_dz
          
          dV_dU(:,:,l,m) = 0.d0
          
          !ci_ro = b_fld%cl(l,m)%y1 / b_fld%cl(l,m)%rho
          do i = 1, nb_especes
             dV_dU(i,irho1z1:irho2z2,l,m) = -U(irho1z1+i-1,l,m) * vol**2.d0
             dV_dU(i,i,l,m) = -U(irho1z1+i-1,l,m) * vol**2.d0 + vol
          end do
          ! pas coherent avec la chimie figee (erreur chez manu?? multiespeces)
          !dV_dU(irho1z1,irho1z1,l,m) =-U(irho1z1,l,m) * vol**2.d0
          !dV_dU(irho1z1,irho2z2,l,m) =-U(irho1z1,l,m) * vol**2.d0 + vol
          
          dV_dU(irho2z2,irho1z1:irho2z2,l,m) =-U(irho2z2,l,m) * vol**2.d0
          dV_dU(irho2z2,irho2z2        ,l,m) =-U(irho2z2,l,m) * vol**2.d0 + vol
          
          dV_dU(irhou,irhou          ,l,m) = vol
          dV_dU(irhou,irho1z1:irho2z2,l,m) =-u_x*vol

          dV_dU(irhov,irhov          ,l,m) = vol
          dV_dU(irhov,irho1z1:irho2z2,l,m) =-u_y*vol

          dV_dU(irhoe,irho1z1:irho2z2-1,l,m) = dT_rho1z1(1:nb_especes)
          dV_dU(irhoe,irho2z2,l,m) = dT_rho2z2
          dV_dU(irhoe,irhou  ,l,m) = dT_rhou
          dV_dU(irhoe,irhov  ,l,m) = dT_rhov
          dV_dU(irhoe,irhoe  ,l,m) = dT_rhoE
          dV_dU(irhoe,iz     ,l,m) = dT_z
       end do
    end do
          
  end subroutine calc_dV_dU
  
  subroutine grad_visqueux_x(b, AA, BB, CC, DD, EE, AD, CD, AE, CE) 
    type (STR_BLOC)                     , pointer       :: b
    real*8, dimension(1:,1:,b%ld:,b%md:), intent(inout) :: AA,BB,CC,DD,EE
    real*8, dimension(1:,1:,b%ld:,b%md:), intent(inout) :: AD,CD,AE,CE

    integer :: i, l, m, lg, ld, mg, md
    real*8 :: u_moy, v_moy, Hi_moy, rhoy_moy
    real*8 :: muJ_moy, KJ_moy, DJ_moy, Jac_g, Jac_d, J_moy
    real*8 :: grad_XY_X(3,2), grad_XY_Y(3,2)
    real*8 :: d_x(3,2), d_y(3,2)
    real*8 :: X_x, X_y, Y_x, Y_y
    real*8 :: flux(b%neq,b%neq,3,2)
    real*8, dimension(:,:,:), pointer :: n_l, n_m
    real*8, dimension(:,:)  , pointer :: dl_l, dl_m
    real*8, dimension(:,:)  , pointer :: u, v

    n_l => b%n_dir_l
    n_m => b%n_dir_m
    dl_l => b%dl_l
    dl_m => b%dl_m
    u => b%u_x
    v => b%u_y

    ! on initialise a zero car il y a des composantes qui ne sont jamais calculees
    flux = 0.d0
    grad_XY_X = 0.d0; grad_XY_Y = 0.d0  
    do m = b%md, b%mf
       do l = b%ld-1, b%lf
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
          u_moy  = .5d0 * ( u(lg,mg) + u(ld,md) )
          v_moy  = .5d0 * ( v(lg,mg) + v(ld,md) )

!!! Jacobiens
          Jac_g = 1.d0 / b%aire(lg,mg)
          Jac_d = 1.d0 / b%aire(ld,md)
          J_moy = moyenne(Jac_g, Jac_d)
          !if (cmp_manu) then
          !   KJ_moy  = J_moy * .5d0*(b%K(lg,mg)+b%K(ld,md))   ! manu
          !   muJ_moy = J_moy * .5d0*(b%mu(lg,mg)+b%mu(ld,md)) ! manu
          !else
             KJ_moy  = moyenne(b%K(lg,mg) *Jac_g, b%K(ld,md) *Jac_d)
             muJ_moy = moyenne(b%mu(lg,mg)*Jac_g, b%mu(ld,md)*Jac_d)
          !end if

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

!!! flux(x,y,l,m) : Pour l'equation sur x, composante relative a y au point (l,m)  
          !
!!$          if (.not. cmp_manu) then
             flux(irhou,irhou,:,:) = n_l(1,l,m)*  4.d0/3.d0 *d_x(:,:) + n_l(2,l,m)*d_y(:,:)
             flux(irhou,irhov,:,:) = n_l(1,l,m)*(-2.d0/3.d0)*d_y(:,:) + n_l(2,l,m)*d_x(:,:)
             !
             flux(irhov,irhou,:,:) = n_l(1,l,m)*d_y(:,:) + n_l(2,l,m)*(-2.d0/3.d0)*d_x(:,:)
             flux(irhov,irhov,:,:) = n_l(1,l,m)*d_x(:,:) + n_l(2,l,m)*  4.d0/3.d0 *d_y(:,:)
             !
             flux(irhoE,irhou,:,:) = flux(irhou,irhou,:,:)*u_moy+flux(irhou,irhov,:,:)*v_moy
             flux(irhoE,irhov,:,:) = flux(irhov,irhou,:,:)*u_moy+flux(irhov,irhov,:,:)*v_moy
             flux(irhoE,irhoE,:,:) = n_l(1,l,m)*d_x(:,:) + n_l(2,l,m)*d_y(:,:)
             !
!!$          else
!!$             ! !!!! simplification 
!!$             gradx2 = n_l(1,l,m)**2.d0 * dl_l(l,m)
!!$             grady2 = n_l(2,l,m)**2.d0 * dl_l(l,m)
!!$
!!$             flux(irhou,irhou,:,:) = 4.d0/3.d0 * gradx2 + grady2 
!!$             flux(irhou,irhov,:,:) = 1.d0/3.d0 * n_l(1,l,m) * n_l(2,l,m) * dl_l(l,m)
!!$             flux(irhov,irhou,:,:) = flux(irhou,irhov,:,:)
!!$             flux(irhov,irhov,:,:) = 4.d0/3.d0 * grady2 + gradx2
!!$             flux(irhoE,irhou,:,:) = flux(irhou,irhou,:,:)*u_moy+flux(irhou,irhov,:,:)*v_moy
!!$             flux(irhoE,irhov,:,:) = flux(irhov,irhou,:,:)*u_moy+flux(irhov,irhov,:,:)*v_moy
!!$             flux(irhoE,irhoE,:,:) = dl_l(l,m)
!!$             flux(:,:,2,1) = -flux(:,:,2,1)
!!$             ! !!!! fin simplification manu
!!$          end if
          flux(:,irhou,:,:) = muJ_moy * dl_l(l,m) * flux(:,irhou,:,:)
          flux(:,irhov,:,:) = muJ_moy * dl_l(l,m) * flux(:,irhov,:,:)
          flux(:,irhoE,:,:) =  KJ_moy * dl_l(l,m) * flux(:,irhoE,:,:)

          ! Partie chimie 
          if ( fluide1%type == CHIMIE_FIGEE ) then 
            ! if (cmp_manu) then
            !    DJ_moy = .5d0 * (b%D(lg,mg) + b%D(ld,md)) * J_moy ! manu
            ! else
                DJ_moy = moyenne(b%D(lg,mg)*Jac_g, b%D(ld,md)*Jac_d)
            ! end if
             rhoy_moy = .5d0 * ( b%y(lg,mg)*b%rho(lg,mg) + b%y(ld,md)*b%rho(ld,md) )

             do i = 1, b%ne
                Hi_moy = .5d0 * ( b%Hi(i,lg,mg) + b%Hi(i,ld,md) )
!!$                if (.not. cmp_manu) then
                   flux(irho1z1+i-1,irho1z1+i-1,:,:) = DJ_moy * rhoy_moy * (n_l(1,l,m)*d_x(:,:)+n_l(2,l,m)*d_y(:,:) ) * dl_l(l,m)
                   flux(irhoe,irho1z1+i-1,:,:) = Hi_moy*flux(irho1z1+i-1,irho1z1+i-1,:,:)
!!$                else
!!$                   ! simplication
!!$                   flux(irho1z1+i-1,irho1z1+i-1,:,:) = DJ_moy * rhoy_moy * dl_l(l,m) * dl_l(l,m)
!!$                   flux(irho1z1+i-1,irho1z1+i-1,2,1) =-flux(irho1z1+i-1,irho1z1+i-1,2,1) 
!!$                   flux(irhoe,irho1z1+i-1,:,:) = DJ_moy * rhoy_moy * Hi_moy * dl_l(l,m) * dl_l(l,m)
!!$                   flux(irhoe,irho1z1+i-1,2,1) =-flux(irhoe,irho1z1+i-1,2,1)
!!$                end if
             end do
          end if

!!! Assemblage des matrices
          ! sm(1,l,m) = sm(1,l,m) - (flux(l,m)-flux(l-1,m))
          if ( l /= b%ld-1 ) then ! sm(1,l,m) = sm(1,l,m) - flux(l,m)
             !DD(:,:,l,m) = DD(:,:,l,m) - flux(:,:,1,1)
             !CD(:,:,l,m) = CD(:,:,l,m) - flux(:,:,1,2) !! coins
             BB(:,:,l,m) = BB(:,:,l,m) - flux(:,:,2,1)
             CC(:,:,l,m) = CC(:,:,l,m) - flux(:,:,2,2)
             !EE(:,:,l,m) = EE(:,:,l,m) - flux(:,:,3,1)
             !CE(:,:,l,m) = CE(:,:,l,m) - flux(:,:,3,2) !!
          end if

          if ( l /= b%lf ) then ! sm(1,l+1,m) = sm(1,l+1,m) + flux(l,m)
             !AD(:,:,l+1,m) = AD(:,:,l+1,m) + flux(:,:,1,1) !!
             !DD(:,:,l+1,m) = DD(:,:,l+1,m) + flux(:,:,1,2)
             AA(:,:,l+1,m) = AA(:,:,l+1,m) + flux(:,:,2,1)
             BB(:,:,l+1,m) = BB(:,:,l+1,m) + flux(:,:,2,2)
             !AE(:,:,l+1,m) = AE(:,:,l+1,m) + flux(:,:,3,1) !!
             !EE(:,:,l+1,m) = EE(:,:,l+1,m) + flux(:,:,3,2)
          end if
       end do
    end do

  end subroutine grad_visqueux_x

  subroutine grad_visqueux_y(b, AA, BB, CC, DD, EE, AD, CD, AE, CE) 
    type (STR_BLOC)                     , pointer       :: b
    real*8, dimension(1:,1:,b%ld:,b%md:), intent(inout) :: AA,BB,CC,DD,EE
    real*8, dimension(1:,1:,b%ld:,b%md:), intent(inout) :: AD,CD,AE,CE

    integer :: i, l, m, lg, ld, mg, md
    real*8 :: u_moy, v_moy, Hi_moy, rhoy_moy
    real*8 :: muJ_moy, KJ_moy, DJ_moy, Jac_g, Jac_d, J_moy
    real*8 :: grad_XY_X(2,3), grad_XY_Y(2,3)
    real*8 :: d_x(2,3), d_y(2,3)
    real*8 :: X_x, X_y, Y_x, Y_y
    real*8 :: flux(b%neq,b%neq,2,3)
    real*8, dimension(:,:,:), pointer :: n_l, n_m
    real*8, dimension(:,:)  , pointer :: dl_l, dl_m
    real*8, dimension(:,:)  , pointer :: u, v

    n_l => b%n_dir_l
    n_m => b%n_dir_m
    dl_l => b%dl_l
    dl_m => b%dl_m
    u => b%u_x
    v => b%u_y

    ! on initialise a zero car il y a des composantes qui ne sont jamais calculees
    flux = 0.d0
    grad_XY_X = 0.d0; grad_XY_Y = 0.d0  
    do m = b%md-1, b%mf
       do l = b%ld, b%lf
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
          u_moy = .5d0 * ( u(lg,mg) + u(ld,md) )
          v_moy = .5d0 * ( v(lg,mg) + v(ld,md) )

!!! Jacobiens
          Jac_g = 1.d0 / b%aire(lg,mg)
          Jac_d = 1.d0 / b%aire(ld,md)
          J_moy = moyenne(Jac_g, Jac_d)
         ! if (cmp_manu) then
         !    KJ_moy  = J_moy * .5d0*(b%K(lg,mg)+b%K(ld,md))   ! manu
         !    muJ_moy = J_moy * .5d0*(b%mu(lg,mg)+b%mu(ld,md)) ! manu
         ! else
             KJ_moy  = moyenne(b%K(lg,mg) *Jac_g, b%K(ld,md) *Jac_d)
             muJ_moy = moyenne(b%mu(lg,mg)*Jac_g, b%mu(ld,md)*Jac_d)
         ! end if
          
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

!!! flux(x,y,l,m) : Pour l'equation sur x, composante relative a y au point (l,m)
          !
!!$          if (.not. cmp_manu) then
             flux(irhou,irhou,:,:) = n_m(1,l,m)*  4.d0/3.d0 *d_x(:,:) + n_m(2,l,m)*d_y(:,:)
             flux(irhou,irhov,:,:) = n_m(1,l,m)*(-2.d0/3.d0)*d_y(:,:) + n_m(2,l,m)*d_x(:,:)
             !
             flux(irhov,irhou,:,:) = n_m(1,l,m)*d_y(:,:) + n_m(2,l,m)*(-2.d0/3.d0)*d_x(:,:)
             flux(irhov,irhov,:,:) = n_m(1,l,m)*d_x(:,:) + n_m(2,l,m)*  4.d0/3.d0 *d_y(:,:)
             !
             flux(irhoE,irhou,:,:) = flux(irhou,irhou,:,:)*u_moy+flux(irhou,irhov,:,:)*v_moy
             flux(irhoE,irhov,:,:) = flux(irhov,irhou,:,:)*u_moy+flux(irhov,irhov,:,:)*v_moy
             flux(irhoE,irhoE,:,:) = n_m(1,l,m)*d_x(:,:) + n_m(2,l,m)*d_y(:,:)          
!!$          else
!!$             ! !!!! simplification
!!$             gradx2 = n_m(1,l,m)**2.d0 * dl_m(l,m) 
!!$             grady2 = n_m(2,l,m)**2.d0 * dl_m(l,m) 
!!$
!!$             flux(irhou,irhou,:,:) = 4.d0/3.d0 * gradx2 + grady2
!!$             flux(irhou,irhov,:,:) = n_m(1,l,m) * n_m(2,l,m) * dl_m(l,m)
!!$             flux(irhov,irhou,:,:) = flux(irhou,irhov,:,:)
!!$             flux(irhov,irhov,:,:) = 4.d0/3.d0 * grady2 + gradx2
!!$             flux(irhoE,irhou,:,:) = flux(irhou,irhou,:,:)*u_moy+flux(irhou,irhov,:,:)*v_moy
!!$             flux(irhoE,irhov,:,:) = flux(irhov,irhou,:,:)*u_moy+flux(irhov,irhov,:,:)*v_moy
!!$             flux(irhoE,irhoE,:,:) = dl_m(l,m)
!!$             flux(:,:,1,2) = -flux(:,:,1,2)
!!$             ! !!!! fin simplification manu 
!!$          end if

          flux(:,irhou,:,:) = muJ_moy * dl_m(l,m) * flux(:,irhou,:,:)
          flux(:,irhov,:,:) = muJ_moy * dl_m(l,m) * flux(:,irhov,:,:)
          flux(:,irhoE,:,:) =  KJ_moy * dl_m(l,m) * flux(:,irhoE,:,:)

          ! Partie chimie 
          if ( fluide1%type == CHIMIE_FIGEE ) then 
            ! if (cmp_manu) then
            !    DJ_moy = .5d0 * (b%D(lg,mg) + b%D(ld,md)) * J_moy ! manu
            ! else
                DJ_moy = moyenne(b%D(lg,mg)*Jac_g, b%D(ld,md)*Jac_d)
            ! end if
             rhoy_moy = .5d0 * ( b%y(lg,mg)*b%rho(lg,mg) + b%y(ld,md)*b%rho(ld,md) )

             do i = 1, b%ne
                Hi_moy = 0.5d0 * ( b%Hi(i,lg,mg) + b%Hi(i,ld,md) )
!!$                if (.not. cmp_manu) then
                   flux(irho1z1+i-1,irho1z1+i-1,:,:) = DJ_moy * rhoy_moy * (n_m(1,l,m)*d_x(:,:) + n_m(2,l,m)*d_y(:,:)) * dl_m(l,m)
                   flux(irhoe,irho1z1+i-1,:,:) = Hi_moy*flux(irho1z1+i-1,irho1z1+i-1,:,:)
!!$                else
!!$                   ! simplifications
!!$                   flux(irho1z1+i-1,irho1z1+i-1,:,:) = DJ_moy * rhoy_moy * dl_m(l,m) * dl_m(l,m)
!!$                   flux(irho1z1+i-1,irho1z1+i-1,1,2) =-flux(irho1z1+i-1,irho1z1+i-1,1,2)
!!$                   flux(irhoe,irho1z1+i-1,:,:) = DJ_moy * rhoy_moy * Hi_moy * dl_m(l,m) * dl_m(l,m)
!!$                   flux(irhoe,irho1z1+i-1,1,2) =-flux(irhoe,irho1z1+i-1,1,2)
!!$                end if

             end do
          end if

!!! Assemblage des matrices
          ! sm(1,l,m) = sm(1,l,m) - (flux(l,m)-flux(l-1,m))
          if ( m /= b%md-1 ) then ! sm(1,l,m) = sm(1,l,m) - flux(l,m)
             !AA(:,:,l,m) = AA(:,:,l,m) - flux(:,:,1,1)
             BB(:,:,l,m) = BB(:,:,l,m) - flux(:,:,1,2)
             !CC(:,:,l,m) = CC(:,:,l,m) - flux(:,:,1,3)
             !AE(:,:,l,m) = AE(:,:,l,m) - flux(:,:,2,1) !! coins
             EE(:,:,l,m) = EE(:,:,l,m) - flux(:,:,2,2)
             !CE(:,:,l,m) = CE(:,:,l,m) - flux(:,:,2,3) !!
          end if

          if ( m /= b%mf ) then ! sm(1,l+1,m) = sm(1,l+1,m) + flux(l,m)
             !AD(:,:,l,m+1) = AD(:,:,l,m+1) + flux(:,:,1,1) !!
             DD(:,:,l,m+1) = DD(:,:,l,m+1) + flux(:,:,1,2)
             !CD(:,:,l,m+1) = CD(:,:,l,m+1) + flux(:,:,1,3) !!
             !AA(:,:,l,m+1) = AA(:,:,l,m+1) + flux(:,:,2,1)
             BB(:,:,l,m+1) = BB(:,:,l,m+1) + flux(:,:,2,2)
             !CC(:,:,l,m+1) = CC(:,:,l,m+1) + flux(:,:,2,3)
          end if

       end do
    end do

  end subroutine grad_visqueux_y
  
  subroutine grad_visqueux_x_manu(b,AA,BB,CC)
    implicit none 
    type (STR_BLOC)                     , pointer       :: b
    double precision, dimension(1:, 1:, b%ld:, b%md:), intent(inout) :: AA,BB,CC
    double precision, dimension(b%ld-1:b%lf,b%md:b%mf) :: Jacobien_moy
    double precision, dimension(b%neq,b%neq, b%ld-1:b%lf) :: MM, flux_visq_g, flux_visq_d
    double precision, dimension(b%ld-1:b%lf) :: cuu, cuv, cvu, cvv
    double precision, dimension(b%ld-1:b%lf) :: u_moy, v_moy
    double precision, dimension(b%ld-1:b%lf) :: grad_jac, mu_moy, mu_jac, K_moy
    double precision, dimension(b%ld-1:b%lf) :: D_moy, rho1z1_moy
    double precision, dimension(:,:), allocatable :: H_i_moy
    double precision :: gradx2, grady2, Jac_1, Jac_2, inv_Jac
    real(kind=8), dimension(b%md:b%mf) :: Nx,Ny

    integer :: i, j, k, l, m
    integer :: nb_especes

    nb_especes = b%ne
    
    if ( fluide1%type == CHIMIE_FIGEE ) then 
       allocate( H_i_moy(1:nb_especes,b%ld-1:b%lf) )
       H_i_moy = 0.d0
    end if

    MM = 0.d0
    flux_visq_g = 0.d0
    flux_visq_d = 0.d0

    ! Calcul de la moyenne harmonique du Jacobien à chaque interface

    do m = b%md, b%mf
       do l = b%ld-1, b%lf
          Jac_1 = 1.d0 / b%Aire(l,m)
          Jac_2 = 1.d0 / b%Aire(l+1,m)
          inv_jac = 0.5d0 * ( 1.d0 / Jac_1 + 1.d0 / Jac_2 )
          Jacobien_moy(l,m) = 1.d0 / inv_jac
       end do
    end do

    Nx = b%n_dir_l(1,b%ld-1, b%md:b%mf)
    Ny = b%n_dir_l(2,b%ld-1, b%md:b%mf)
    
    ! a voir
    !call implicitation_cond_visq(b%West, size(b%West), b%dv_du(1:b%neq,1:b%neq,b%ld,b%md:b%mf), & 
    !     & b%dv_du(1:b%neq,1:b%neq,b%ld-1,b%md:b%mf),nx,ny,b%md,b%mf)	

    Nx = b%n_dir_l(1,b%lf, b%md:b%mf)
    Ny = b%n_dir_l(2,b%lf, b%md:b%mf)

    !call implicitation_cond_visq(b%East, size(b%East), b%dv_du(1:b%neq,1:b%neq,b%lf,b%md:b%mf), & 
    !     & b%dv_du(1:b%neq,1:b%neq,b%lf+1,b%md:b%mf),nx,ny,b%md,b%mf)	
    

    do m = b%md, b%mf

       MM = 0.d0
       flux_visq_g = 0.d0
       flux_visq_d = 0.d0

       do l = b%ld-1, b%lf
          gradx2 = ( b%n_dir_l(1,l,m)*b%dl_l(l,m) )**2.d0
          grady2 = ( b%n_dir_l(2,l,m)*b%dl_l(l,m) )**2.d0

          u_moy(l) = .5d0 * ( b%u_x(l,m) + b%u_x(l+1,m) )
          v_moy(l) = .5d0 * ( b%u_y(l,m) + b%u_y(l+1,m) )
          
          K_moy(l) =  0.5d0 * ( b%K(l,m) + b%K(l+1,m) )
          Mu_moy(l) =  0.5d0 * ( b%mu(l,m) + b%mu(l+1,m) )
          mu_jac(l) = mu_moy(l) * Jacobien_moy(l,m)
          grad_jac(l) = ( gradx2 + grady2 ) * Jacobien_moy(l,m)

          cuu(l) = mu_jac(l) * (4.d0/3.d0 * gradx2 + grady2 )
          cuv(l) = mu_jac(l) * 1.d0 / 3.d0 * b%n_dir_l(1,l,m) * b%n_dir_l(2,l,m) * b%dl_l(l,m)**2
          cvu(l) = cuv(l)
          cvv(l) =  mu_jac(l) * (4.d0/3.d0 * grady2 + gradx2 )

          ! Partie chimie 
          if ( fluide1%type == CHIMIE_FIGEE ) then 
             D_moy(l) = 0.5d0 * ( b%D(l,m) + b%D(l+1,m) )
             rho1z1_moy(l) = 0.5d0 * ( b%y(l,m) * b%rho(l,m) + b%y(l+1,m) * b%rho(l+1,m) )

             do i = 1, nb_especes
                H_i_moy(i,l) = 0.5d0 * ( b%Hi(i,l,m) + b%Hi(i,l+1,m) )
             end do
          end if

       end do

       do l = b%ld-1, b%lf
          MM(irhou,irhou,l) = cuu(l)
          MM(irhou,irhov,l) = cuv(l)

          MM(irhov,irhou,l) = cvu(l)
          MM(irhov,irhov,l) = cvv(l)

          MM(irhoe,irhou,l) = MM(irhou,irhou,l) * u_moy(l) + MM(irhou,irhov,l) * v_moy(l) 
          MM(irhoe,irhov,l) = MM(irhov,irhou,l) * u_moy(l) + MM(irhov,irhov,l) * v_moy(l)
          MM(irhoe,irhoe,l) = grad_jac(l) * K_moy(l)

          ! Partie chimie 
          if ( fluide1%type == CHIMIE_FIGEE ) then 
             do i = 1, nb_especes
                MM(irho1z1+i-1,irho1z1+i-1,l) = grad_jac(l) * D_moy(l) * rho1z1_moy(l)
                MM(irhoe,irho1z1+i-1,l) = grad_jac(l) * H_i_moy(i,l) * D_moy(l) * rho1z1_moy(l)
             end do
          end if

          do k = 1, b%neq
             do j = 1, b%neq
                do i = 1, b%neq
                   flux_visq_g(i,j,l) = flux_visq_g(i,j,l) + MM(i,k,l) * b%dv_du(k,j,l,m)
                   flux_visq_d(i,j,l) = flux_visq_d(i,j,l) + MM(i,k,l) * b%dv_du(k,j,l+1,m)
                end do
             end do
          end do
       end do

       do l = b%ld, b%lf-1
          do j = 1, b%neq
             do i = 1, b%neq
                AA(i,j,l+1,m) = AA(i,j,l+1,m) - flux_visq_g(i,j,l)
                CC(i,j,l,m) = CC(i,j,l,m) - flux_visq_d(i,j,l)
             end do
          end do
       end do

       do l = b%ld, b%lf
          do j = 1, b%neq
             do i = 1, b%neq
                BB(i,j,l,m) = BB(i,j,l,m) + flux_visq_d(i,j,l-1)
                BB(i,j,l,m) = BB(i,j,l,m) + flux_visq_g(i,j,l)
             end do
          end do
       end do

       do j = 1, b%neq
          do i = 1, b%neq
             BB(i,j,b%ld,m) = BB(i,j,b%ld,m) - flux_visq_g(i,j,b%ld-1)
             BB(i,j,b%lf,m) = BB(i,j,b%lf,m) - flux_visq_d(i,j,b%lf)
          end do
       end do

    end do

    if ( fluide1%type == CHIMIE_FIGEE ) deallocate( H_i_moy )

  end subroutine grad_visqueux_x_manu
  
  
  subroutine grad_visqueux_y_manu(b,AA,BB,CC)
    implicit none 
    type (STR_BLOC)                     , pointer       :: b
    double precision, dimension(1:, 1:, b%ld:, b%md:), intent(inout) :: AA,BB,CC
    double precision, dimension(b%ld:b%lf,b%md-1:b%mf) :: Jacobien_moy
    double precision, dimension(b%neq,b%neq, b%md-1:b%mf) :: MM, flux_visq_g, flux_visq_d
    double precision, dimension(b%md-1:b%mf) :: cuu, cuv, cvu, cvv
    double precision, dimension(b%md-1:b%mf) :: u_moy, v_moy
    double precision, dimension(b%md-1:b%mf) :: grad_jac, mu_moy, mu_jac, K_moy
    double precision, dimension(b%md-1:b%mf) :: D_moy, rho1z1_moy
    double precision, dimension(:,:), allocatable :: H_i_moy
    double precision :: gradx2, grady2, Jac_1, Jac_2, inv_jac
    real(kind=8), dimension(b%ld:b%lf) :: Nx,Ny

    integer :: i, j, k, l, m
    integer :: nb_especes

    nb_especes = b%ne
    
    if ( fluide1%type == CHIMIE_FIGEE ) then 
       allocate( H_i_moy(1:nb_especes,b%md-1:b%mf) )
       H_i_moy = 0.d0
    end if

    MM = 0.d0
    flux_visq_g = 0.d0
    flux_visq_d = 0.d0

    ! Calcul de la moyenne harmonique du Jacobien à chaque interface

    do m = b%md-1, b%mf
       do l = b%ld, b%lf
          Jac_1 = 1.d0 / b%Aire(l,m)
          Jac_2 = 1.d0 / b%Aire(l,m+1)
          inv_jac = 0.5d0 * ( 1.d0 / Jac_1 + 1.d0 / Jac_2 ) 
          Jacobien_moy(l,m) = 1.d0 / inv_jac
       end do
    end do

    Nx = b%n_dir_m(1,b%ld:b%lf,b%md-1)
    Ny = b%n_dir_m(2,b%ld:b%lf,b%md-1)

    !call implicitation_cond_visq(b%South, size(b%South), b%dv_du(1:b%neq,1:b%neq,b%ld:b%lf,b%md), & 
    !     & b%dv_du(1:b%neq,1:b%neq,b%ld:b%lf,b%md-1),nx,ny,b%ld,b%lf)
    b%dv_du(1:b%neq,1:b%neq,b%ld:b%lf,b%md-1) = 0.d0
    
    
    Nx = b%n_dir_m(1,b%ld:b%lf,b%mf)
    Ny = b%n_dir_m(2,b%ld:b%lf,b%mf)

    !call implicitation_cond_visq(b%North, size(b%North), b%dv_du(1:b%neq,1:b%neq,b%ld:b%lf,b%mf), & 
!         & b%dv_du(1:b%neq,1:b%neq,b%ld:b%lf,b%mf+1),nx,ny,b%ld,b%lf)
    b%dv_du(1:b%neq,1:b%neq,b%ld:b%lf,b%mf+1) = 0.d0
    
    do l = b%ld, b%lf

       MM = 0.d0
       flux_visq_g = 0.d0
       flux_visq_d = 0.d0

       do m = b%md-1, b%mf

          gradx2 = ( b%n_dir_m(1,l,m)*b%dl_m(l,m) )**2.d0
          grady2 = ( b%n_dir_m(2,l,m)*b%dl_m(l,m) )**2.d0
          
          u_moy(m) = .5d0 * ( b%u_x(l,m) + b%u_x(l,m+1) )
          v_moy(m) = .5d0 * ( b%u_y(l,m) + b%u_y(l,m+1) )
          
          K_moy(m) =  0.5d0 * ( b%K(l,m) + b%K(l,m+1) )
          Mu_moy(m) =  0.5d0 * ( b%mu(l,m) + b%mu(l,m+1) ) 
          mu_jac(m) = mu_moy(m) * Jacobien_moy(l,m)
          grad_jac(m) = ( gradx2 + grady2 ) * Jacobien_moy(l,m)


          cuu(m) = mu_jac(m) * ( 4.d0/3.d0 * gradx2 + grady2 )
          cuv(m) = mu_jac(m) *  b%n_dir_m(1,l,m) * b%n_dir_m(2,l,m) * b%dl_m(l,m)**2.d0
          cvu(m) = cuv(m)
          cvv(m) = mu_jac(m) * ( 4.d0/3.d0 * grady2 + gradx2 )

          ! Partie chimie 
          if ( fluide1%type == CHIMIE_FIGEE ) then 
             D_moy(m) = 0.5d0 * ( b%D(l,m) + b%D(l,m+1) )
             rho1z1_moy(m) = 0.5d0 * ( b%y(l,m) * b%rho(l,m) + b%y(l,m+1) * b%rho(l,m+1) )

             do i = 1, nb_especes
                H_i_moy(i,m) = 0.5d0 * ( b%Hi(i,l,m) + b%Hi(i,l,m+1) )
             end do
          end if
       end do

       do m = b%md-1, b%mf
          MM(irhou,irhou,m) = cuu(m)
          MM(irhou,irhov,m) = cuv(m)

          MM(irhov,irhou,m) = cvu(m)
          MM(irhov,irhov,m) = cvv(m)

          MM(irhoe,irhou,m) = MM(irhou,irhou,m) * U_moy(m) + MM(irhou,irhov,m) * v_moy(m)
          MM(irhoe,irhov,m) = MM(irhov,irhou,m) * U_moy(m) + MM(irhov,irhov,m) * v_moy(m)
          MM(irhoe,irhoe,m) = grad_jac(m) * K_moy(m)

          ! Partie chimie 
          if ( fluide1%type == CHIMIE_FIGEE ) then 
             do i = 1, nb_especes
                MM(irho1z1+i-1,irho1z1+i-1,m) = grad_jac(m) * D_moy(m) * rho1z1_moy(m)
                MM(irhoe,irho1z1+i-1,m) = grad_jac(m) * H_i_moy(i,m) * D_moy(m) * rho1z1_moy(m)
             end do
          end if

          do k = 1, b%neq
             do j = 1, b%neq
                do i = 1, b%neq
                   flux_visq_g(i,j,m) = flux_visq_g(i,j,m) + MM(i,k,m) * b%dv_du(k,j,l,m)
                   flux_visq_d(i,j,m) = flux_visq_d(i,j,m) + MM(i,k,m) * b%dv_du(k,j,l,m+1)
                end do
             end do
          end do
       end do

       do m = b%md+1, b%mf
          do j = 1, b%neq
             do i = 1, b%neq
                AA(i,j,l,m) = AA(i,j,l,m) - flux_visq_g(i,j,m-1)
             end do
          end do
       end do

       do m = b%md, b%mf
          do j = 1, b%neq
             do i = 1, b%neq
                BB(i,j,l,m) = BB(i,j,l,m) + flux_visq_d(i,j,m-1)
                BB(i,j,l,m) = BB(i,j,l,m) + flux_visq_g(i,j,m)
             end do
          end do
       end do

       do m = b%md, b%mf-1
          do j = 1, b%neq
             do i = 1, b%neq
                CC(i,j,l,m) = CC(i,j,l,m) - flux_visq_d(i,j,m)
             end do
          end do
       end do

       do j = 1, b%neq
          do i = 1, b%neq
             BB(i,j,l,b%md) = BB(i,j,l,b%md) - flux_visq_g(i,j,b%md-1)
             BB(i,j,l,b%mf) = BB(i,j,l,b%mf) - flux_visq_d(i,j,b%mf)
          end do
       end do

    end do

    if ( fluide1%type == CHIMIE_FIGEE ) deallocate( H_i_moy )

  end subroutine grad_visqueux_y_manu


end module premier_visqueux_direct
