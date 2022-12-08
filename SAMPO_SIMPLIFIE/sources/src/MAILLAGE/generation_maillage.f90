module generation_maillage
  use m_struct
  use m_MPI
  use outils_maillage
  use parametres_globaux

  implicit none

!!! generation de differents types de maillages possible
  ! maillage_cartesien_oblique
  ! maillage_ellipse_cartesien
  ! maillage_cartesien
  ! maillage_a_partir_de_4pts
  ! maillage_ellipse
  ! maillage_sphere
  ! maillage_avec_bosse_dellacherie
  ! maillage_plaque_plane_avec_creux
  ! maillage_random
  ! maillage_kershaw
contains

  subroutine maillage_ellipse_cartesien(b)
    use variables_globales, only : path_input

    type (STR_BLOC), pointer :: b

    integer :: l,m
    integer :: ld_sphere,lf_sphere,md_sphere,mf_sphere
    integer :: ld_cart,lf_cart,md_cart,mf_cart
    real*8 :: Ray, dR, Rayon_donne,dtheta
    real*8 :: r1, r2, theta1, theta2, excentricite2, rayon_ellipse, coef_b2_ellipse
    real*8 :: theta, centre_sphere(2)
    real*8 :: xa, ya, xc, yc
    real*8 :: dx, dy
    real*8 :: xb,yb,xd,yd
    integer :: ld, lf, md, mf
    real*8 :: xin, xout, yin, yout
    real*8 :: XXright,eps,xdeb,ydeb,u,v,x,y,D
    real*8 :: epsilon,Ratio,alpha
    real*8 :: dx1,dx2,dy1,dy2
    character(len=200) :: fichier_maillage

    ! Inclusion des mailles fictives via nMF ?

    ! ******************************************************
    ! Partie sphérique

    md_sphere=b%md; mf_sphere=b%mf;
    md_cart=b%md; mf_cart=b%mf;

    ! De base, 1/2 mailles sur la partie sphère et 1/2 sur la partie rectangle
    ld_sphere=b%ld; lf_sphere=b%lf/2;
    ld_cart=lf_sphere+1; lf_cart=b%lf;

    write(fichier_maillage, "('maillage_ellipse_cartesien_',i1.1,'.dat')") b%nature
    open(unit=20,file=trim(path_input)//trim(fichier_maillage), status="old") 
    read(20,*) centre_sphere
    read(20,*) rayon_donne
    read(20,*) r1   ! demi axe suivant x
    read(20,*) r2   ! demi axe suivant y
    read(20,*) theta1
    read(20,*) theta2
    read(20,*) D
    read(20,*) alpha
    read(20,*) epsilon
    close(20)

    ! if (r_max /= b%yout-b%yin) then
    !    print*, "Warning : Le rayon exterieur donne en parametre dans maillage.dat ne correspond au domaine decrit dans input.dat"
    !    r_max = b%yout-b%yin
    ! end if
    ! if (theta1 >= theta2) then
    !    print*, "Error : theta1 est superieur à theta2"
    !    call arret_code
    ! end if

    ! Conversion en radian
    theta1= acos(-1.d0)*theta1/180.d0
    theta2= acos(-1.d0)*theta2/180.d0

    XXright = D*sin(theta2)-rayon_donne*cos(theta2)

    ! ***********************************************************************

    if(epsilon < 1.0) then
       Ratio = int(1./epsilon)
       ld_sphere=b%ld; lf_sphere=b%lf/(1+Ratio);
       ld_cart=lf_sphere+1; lf_cart=b%lf
    else
       Ratio = int(epsilon)
       ld_sphere=b%ld; lf_sphere=b%lf/(1+Ratio)*Ratio;
       ld_cart=lf_sphere+1; lf_cart=b%lf
    endif

    ! ***********************************************************************

    ! Ellipse
    if (cos(theta1) == cos(theta2)) then
       excentricite2 = ( r1**2.d0 - r2**2.d0 ) / r1**2.d0
       coef_b2_ellipse = r1**2.0 * ( 1.d0 - excentricite2 )
    else
       excentricite2 = ( r1**2.d0 - r2**2.d0 ) / ( r1**2.d0 * cos( theta1 )**2.d0- r2**2.d0 * cos( theta2 )**2.d0)
       coef_b2_ellipse = r1**2.0 * ( 1.d0 - excentricite2 * cos( theta1 )**2.d0 )
    end if

    dtheta = abs(theta2-theta1)/dble(lf_sphere-ld_sphere+1)

    do l=b%ld-3,lf_sphere
       theta = theta1+dtheta*dble(l+b%l_charge(1,numproc)-1)
       rayon_ellipse =  sqrt( coef_b2_ellipse /( 1.d0  -  excentricite2 * cos(theta)**2.d0 ))
       dR = (rayon_ellipse - Rayon_donne)/dble(b%mf_tot-b%md_tot+1)

       xa=-r2*cos(theta)
       ya= r2*sin(theta)
       xc=-rayon_donne*cos(theta)
       yc= rayon_donne*sin(theta)

       do m=b%md-3,b%mf+2

          Ray= Rayon_donne + dR * ( m - b%md + 1 + b%m_charge(1,numproc)-1)
          b%coord_nd(1,l,m)= -Ray * cos(theta) + centre_sphere(1)
          b%coord_nd(2,l,m)= +Ray * sin(theta) + centre_sphere(2)

          ! b%coord_nd(1,l,m)=xc+dble(m)*(xa-xc)/dble(mf_sphere) + centre_sphere(1)
          ! b%coord_nd(2,l,m)=yc+dble(m)*(ya-yc)/dble(mf_sphere) + centre_sphere(2)
       end do

    end do

    ! ******************************************************
    ! Partie cartésienne

    eps = r2-rayon_donne

    xc = centre_sphere(1)-rayon_donne*cos(theta2)
    yc = centre_sphere(2)+rayon_donne*sin(theta2)

    xd = centre_sphere(1)+XXright
    yd = yc + (XXright+rayon_donne*cos(theta2))/sin(theta2)*cos(theta2)

    xa = centre_sphere(1)-r2*cos(theta2)
    ya = centre_sphere(2)+r2*sin(theta2)

    xb = xd - eps*cos(theta2)
    yb = yd + eps*sin(theta2)

    ! Propre à l'ellipse
    ! alpha = acos(-1.d0)*15./180.d0
    alpha = acos(-1.d0)*alpha/180.d0
    xb = xb - D/cos(alpha)*sin(alpha)*cos(theta2)
    yb = yb + D/cos(alpha)*sin(alpha)*sin(theta2)

    ! ici dx = xd-xc = xb-xa
    dx1 = (xd-xc)/dble(lf_cart - ld_cart+1)
    dx2 = (xb-xa)/dble(lf_cart - ld_cart+1)
    ! dx2 = dx1
    ! print*,dx1,dx2

    ! ici dy = yb-yd = ya-yc
    dy1 = (ya-yc)/dble(mf_cart - md_cart+1)
    dy2 = (yb-yd)/dble(mf_cart - md_cart+1)
    ! dy2 = dy1
    ! print*,dy1,dy2

    dy = (eps*sin(theta2))/dble(mf_cart - md_cart+1)
    do m=b%md-3,b%mf+2
       y = yc+dy*dble(m+b%m_charge(1,numproc)-1)
       ! y = yc+(1.-v)*dy1*dble(m+b%m_charge(1,numproc)-1)+v*dy2*dble(m+b%m_charge(1,numproc)-1)
       x = 1./(ya-yc)*(xc*(ya-y)+xa*(y-yc))
       do l=ld_cart,lf_cart+2
          u = (dble(l-lf_sphere+b%ld+b%l_charge(1,numproc)-1)-dble(ld_sphere+b%l_charge(1,numproc)-1))/dble(lf_cart - ld_cart+1)
          v = dble(m+b%m_charge(1,numproc)-1)/dble(mf_cart - md_cart+1)
          dx = (1.-v)*dx1+v*dx2
          b%coord_nd(1,l,m) = x + dx*(dble(l-lf_sphere+b%ld+b%l_charge(1,numproc)-1)-dble(ld_sphere+b%l_charge(1,numproc)-1))
       enddo
    end do

    dx = (XXright+rayon_donne*cos(theta2))/dble(lf_cart - ld_cart+1)
    do l=ld_cart,lf_cart+2
       x = xc+dx*(dble(l-lf_sphere+b%ld+b%l_charge(1,numproc)-1)-dble(ld_sphere+b%l_charge(1,numproc)-1))
       ! x = xc+(1.-u)*dx1*(dble(l-lf_sphere+b%ld+b%l_charge(1,numproc)-1)-dble(ld_sphere+b%l_charge(1,numproc)-1))+u*dx2*(dble(l-lf_sphere+b%ld+b%l_charge(1,numproc)-1)-dble(ld_sphere+b%l_charge(1,numproc)-1))
       y = 1./(xd-xc)*(yc*(xd-x)+yd*(x-xc))
       do m=b%md-3,b%mf+2
          u = (dble(l-lf_sphere+b%ld+b%l_charge(1,numproc)-1)-dble(ld_sphere+b%l_charge(1,numproc)-1))/dble(lf_cart - ld_cart+1)
          v = dble(m+b%m_charge(1,numproc)-1)/dble(mf_cart - md_cart+1)
          dy = (1.-u)*dy1+u*dy2
          b%coord_nd(2,l,m) = y + dy*dble(m+b%m_charge(1,numproc)-1)
       enddo
    end do


    ! ******************************************************



  end subroutine maillage_ellipse_cartesien








  subroutine maillage_cartesien(b)

    type (STR_BLOC), pointer :: b
    integer :: l,m
    integer :: ld, lf, md, mf
    real*8 :: xin, xout, yin, yout
    real*8 :: dx, dy

    xin = b%xin; xout = b%xout
    yin = b%yin; yout = b%yout

    ld=b%ld; lf=b%lf
    md=b%md; mf=b%mf

    dx = (xout - xin)/dble(b%lf_tot - b%ld_tot+1)
    dy = (yout - yin)/dble(b%mf_tot - b%md_tot+1)

    do m= md-(nMF+1), mf+nMF
       do l = ld-(nMF+1), lf+nMF
          b%coord_nd(1,l,m)= dx*dble(l+b%l_charge(1,numproc)-1) + xin
          b%coord_nd(2,l,m)= dy*dble(m+b%m_charge(1,numproc)-1) + yin
       end do
    end do

  end subroutine maillage_cartesien



  subroutine maillage_a_partir_de_4pts(b)
    use variables_globales, only : path_input

    type (STR_BLOC), pointer :: b
    integer :: l,m
    integer :: ld, lf, md, mf
    real*8 :: dx
    real*8 :: p1(2), p2(2), p3(2), p4(2), e1(2)
    real*8 :: dS, dN, dW, dE
    real*8 :: ndS(2), ndN(2), ndW(2), ndE(2)
    real*8, dimension(2,b%md-(nMF+1):b%mf+nMF) :: tabW, tabE
    real*8, dimension(2,b%ld-(nMF+1):b%lf+nMF) :: tabS, tabN

    ! lecture des 4 pts
    open(unit=20,file=trim(path_input)//"maillage.dat")
    read(20,*) p1
    read(20,*) p2
    read(20,*) p3
    read(20,*) p4
    close(20)

    b%xin  = p1(1) ; b%xout = p3(1)
    b%yin  = p1(2) ; b%yout = p3(2)

    ld=b%ld; lf=b%lf
    md=b%md; mf=b%mf

    ! pas d'espace des 4 aretes
    dS = norme(p1-p2)/dble(b%lf_tot - b%ld_tot+1)
    dN = norme(p3-p4)/dble(b%lf_tot - b%ld_tot+1)
    dW = norme(p1-p4)/dble(b%mf_tot - b%md_tot+1)
    dE = norme(p3-p2)/dble(b%mf_tot - b%md_tot+1)

    !les normales des 4 aretes
    ndS = (p2-p1)/norme(p1-p2)
    ndN = (p3-p4)/norme(p3-p4)
    ndW = (p4-p1)/norme(p4-p1)
    ndE = (p3-p2)/norme(p3-p2)

    do m = b%md-(nMF+1),b%mf+nMF
       tabE(:,m) = p2 + ndE*dE*dble(m+b%m_charge(1,numproc)-1)
       tabW(:,m) = p1 + ndW*dW*dble(m+b%m_charge(1,numproc)-1)
    end do

    do l = b%ld-(nMF+1),b%lf+nMF
       tabS(:,l) = p1 + ndS*dS*dble(l+b%l_charge(1,numproc)-1)
       tabN(:,l) = p4 + ndN*dN*dble(l+b%l_charge(1,numproc)-1)
    end do


    do m = md-(nMF+1), mf+nMF
       e1 = (tabE(:,m)-tabW(:,m))/norme(tabE(:,m)-tabW(:,m))
       dx = norme(tabE(:,m)-tabW(:,m))/dble(b%lf_tot - b%ld_tot+1)
       do l = ld-(nMF+1), lf+nMF
          b%coord_nd(1,l,m)= e1(1)*dx*dble(l+b%l_charge(1,numproc)-1) + tabW(1,m)
          b%coord_nd(2,l,m)= e1(2)*dx*dble(l+b%l_charge(1,numproc)-1) + tabW(2,m)
       end do
    end do

  contains
    function norme(u)
      real*8 :: u(2), norme

      norme = sqrt(u(1)*u(1) + u(2)*u(2))
    end function norme

  end subroutine maillage_a_partir_de_4pts

  subroutine maillage_ellipse(b)
    use variables_globales, only : path_input

    type (STR_BLOC), pointer :: b

    integer :: l,m
    real*8 :: Ray, dR, Rayon_donne, dtheta
    real*8 :: r1, r2, theta1, theta2, excentricite2, rayon_ellipse, coef_b2_ellipse, theta, centre_sphere(2)

    ! Inclusion des mailles fictives via nMF ?

    open(unit=20,file=trim(path_input)//"maillage.dat")
    read(20,*) centre_sphere
    read(20,*) rayon_donne
    read(20,*) r1   ! demi axe suivant x
    read(20,*) r2   ! demi axe suivant y
    read(20,*) theta1
    read(20,*) theta2
    close(20)

    if (r1 /= b%xout-b%xin .or. r2 /= b%yout-b%yin) then
       print*, "Warning : Les parametres donnes dans maillage.dat ne correspondent au domaine decrit dans input.dat"
       r1 = b%xout-b%xin
       r2 = b%yout-b%yin
    end if
    ! if (theta1 /= 0.d0 .or. theta2 /= 90.d0) then
    !    print*, "Error : cas non gere, on doit avoir theta1=0, theta2=90"
    !    call arret_code
    ! end if
    if (theta1 >= theta2) then
       print*, "Error : theta_min est superieur à theta_max", theta1, theta2
       call arret_code
    end if

    ! Conversion en radian
    theta1 = acos(-1.d0)*theta1/180.d0
    theta2 = acos(-1.d0)*theta2/180.d0

    ! Ellipse
    if (cos(theta1) == cos(theta2)) then
       excentricite2 = ( r1**2.d0 - r2**2.d0 ) / r1**2.d0
       coef_b2_ellipse = r1**2.0 * ( 1.d0 - excentricite2 )
    else
       excentricite2 = ( r1**2.d0 - r2**2.d0 ) / ( r1**2.d0 * cos( theta1 )**2.d0- r2**2.d0 * cos( theta2 )**2.d0)
       coef_b2_ellipse = r1**2.0 * ( 1.d0 - excentricite2 * cos( theta1 )**2.d0 )
    end if

    dtheta = (theta2-theta1)/dble(b%lf_tot-b%ld_tot+1)

    do l = b%ld-3, b%lf+2
       theta = theta1 + dble(l+b%l_charge(1,numproc)-1)*dtheta
       rayon_ellipse =  sqrt( coef_b2_ellipse /( 1.d0  -  excentricite2 * cos(theta)**2.d0 ))
       dR = (rayon_ellipse - Rayon_donne)/dble(b%mf_tot-b%md_tot+1)
       !Ray=rayon
       do m= b%md-3, b%mf+2
          Ray= Rayon_donne + dR * ( m - b%md + 1 + b%m_charge(1,numproc)-1)
          b%coord_nd(1,l,m)= -Ray * cos(theta) + centre_sphere(1)
          b%coord_nd(2,l,m)= +Ray * sin(theta) + centre_sphere(2)
       end do
    end do

  end subroutine maillage_ellipse

  subroutine maillage_sphere(b)
    use variables_globales, only : path_input

    type (STR_BLOC), pointer :: b

    integer :: l,m
    real*8 :: r_min, r_max, dtheta
    real*8 :: theta, theta_min, theta_max, centre_sphere(2)
    real*8 :: xa, ya, xc, yc

    ! Inclusion des mailles fictives via nMF ?

    open(unit=20,file=trim(path_input)//"maillage_.dat")
    read(20,*) centre_sphere
    read(20,*) r_min
    read(20,*) r_max
    read(20,*) theta_min
    read(20,*) theta_max
    close(20)

    if (r_max /= b%yout-b%yin) then
       print*, "Warning : Le rayon exterieur donne en parametre dans maillage.dat ne correspond au domaine decrit dans input.dat"
       r_max = b%yout-b%yin
    end if
    if (theta_min >= theta_max) then
       print*, "Error : theta_min est superieur à theta_max"
       call arret_code
    end if

    ! Conversion en radian
    theta_min = acos(-1.d0)*theta_min/180.d0
    theta_max = acos(-1.d0)*theta_max/180.d0

    dtheta = abs(theta_max-theta_min)/dble(b%lf_tot-b%ld_tot+1)

    do l=b%ld-3,b%lf+2
       theta = theta_min+dtheta*dble(l+b%l_charge(1,numproc)-1)

       xa=-r_max*cos(theta)
       ya= r_max*sin(theta)
       xc=-r_min*cos(theta)
       yc= r_min*sin(theta)
       do m=b%md-3,b%mf+2
          b%coord_nd(1,l,m)=xc+dble(m - b%md + 1 + b%m_charge(1,numproc)-1)*(xa-xc)/dble(b%mf_tot-b%md_tot+1) + centre_sphere(1)
          b%coord_nd(2,l,m)=yc+dble(m - b%md + 1 + b%m_charge(1,numproc)-1)*(ya-yc)/dble(b%mf_tot-b%md_tot+1) + centre_sphere(2)
       end do
    end do

  end subroutine maillage_sphere

  subroutine maillage_avec_bosse_dellacherie(b)
    implicit none

    type (STR_BLOC), pointer :: b
    real*8 :: hx, hy_1, hy_2
    integer :: l,m
    !   real*8 :: longueur, H , xc  ! H hauteur bosse ; l demi-largeur ; xc centre bosse
    real*8, parameter :: pi = acos(-1.d0)

    ! Inclusion des mailles fictives via nMF ?

    ! Cas test défini par dellacherie dans l'article
    ! Analysis of Godunov type schemes applied to compressible Euler systeme at low Mach number, JCP, 2009
    b%xout = 4.d0
    b%xin = 0.d0
    b%yout = 1.d0
    b%yin = 0.d0

    do l = b%ld-3, b%lf+2

       hx=(b%xout-b%xin)/dble(b%lf_tot)

       do m= b%md-3, b%mf+2

          b%coord_nd(1,l,m)=dble(l+b%l_charge(1,numproc)-1)*hx+b%xin

          if ((1.d0<b%coord_nd(1,l,m)) .and. (b%coord_nd(1,l,m)<3.d0)) then
             hy_2=b%yout
             hy_1= 0.1d0 * ( 1.d0 - cos((b%coord_nd(1,l,m)-1.d0) * pi) )
             b%coord_nd(2,l,m) = hy_1+dble(m+b%m_charge(1,numproc)-1)/dble(b%mf_tot)*(hy_2-hy_1)
          else
             hy_2=b%yout
             hy_1=b%yin
             b%coord_nd(2,l,m) = hy_1+dble(m+b%m_charge(1,numproc)-1)/dble(b%mf_tot)*(hy_2-hy_1)
          end if

       end do
    end do

  end subroutine maillage_avec_bosse_dellacherie


  subroutine maillage_plaque_plane_avec_creux(b)
    use variables_globales, only : path_input

    type (STR_BLOC), pointer :: b

    integer :: m,l
    real*8 :: y_plat1, y_plat2, y_creux, x_deb, dx_smooth1, dx_smooth2, dx_creux, dy, dx, x, y_deb, x_fin
    real*8 :: facteur_echelle, Lx, Ly
    logical :: zone_profil_couche_limite
    real*8, parameter :: pi = acos(-1.d0)

    open(unit=20,file=trim(path_input)//"maillage.dat")
    read(20,*) y_plat1, y_creux, y_plat2
    read(20,*) x_deb, dx_smooth1, dx_creux, dx_smooth2
    read(20,*) zone_profil_couche_limite
    close(20)

    facteur_echelle = 1.d-3 ! dessin de l'ordre du minimètre
    !b%xin = 0.d0 ;  b%xout = 2.d0 * facteur_echelle
    !b%yin = 0.d0 ;  b%yout = 2.d0 * facteur_echelle
    Lx = b%xout - b%xin
    Ly = b%yout - b%yin

    dx = Lx / dble(b%lf_tot-b%ld_tot+1)
    x_fin = x_deb + dx_smooth1 + dx_smooth2 + dx_creux

    do l = b%ld-(nMF+1), b%lf+nMF
       x = b%xin + dx * dble(l+b%l_charge(1,numproc)-1)
       if ( x <  x_deb ) then
          y_deb = y_plat1
       else if ( x >  x_fin ) then
          y_deb = y_plat2
       else
          if ( x < x_deb + dx_smooth1 ) then
             y_deb = y_plat1 + (y_plat1-y_creux) * (cos(pi * (x-x_deb)/(dx_smooth1))-1.d0)/2.d0
          else if ( x > x_fin - dx_smooth2 ) then
             y_deb = y_plat2 + (y_plat2-y_creux) * (cos(pi * (x_fin-x)/(dx_smooth2))-1.d0)/2.d0
          else
             y_deb = y_creux
          end if
       end if
       dy = (b%yout-y_deb) / dble(b%mf_tot-b%md_tot+1)
       do m = b%md-(nMF+1), b%mf+nMF
          b%coord_nd(1,l,m) =  x
          b%coord_nd(2,l,m) =  y_deb + dble(m+b%m_charge(1,numproc)-1) * dy
       end do
    end do

    ! on agrandit les premieres mailles pour permettre instauration d'une couche limite de type plaque_plane

    if ( zone_profil_couche_limite ) then
       do m = b%md-(nMF+1), b%mf+nMF
          do l = b%ld-(nMF+1),9
             if ( b%coord_nd(1,l,m) <= b%xin + dx * 9.d0 ) then
                b%coord_nd(1,l,m) = b%xin + dx * dble(l-9) * 20.d0 ! 20 ou 40.d0
             end if
          end do
       end do
    end if

  end subroutine maillage_plaque_plane_avec_creux

  subroutine maillage_random(b)
    use variables_globales, only : path_input

    type (STR_BLOC), pointer :: b

    ! variables locales
    integer :: l, m, seed
    integer :: ld, lf, md, mf
    integer :: ldeb, lfin, mdeb, mfin
    real*8  :: xin, xout, yin, yout
    real*8  :: dx, dy
    real*8  :: taux

    seed=0
    call srand(seed)

    xin = b%xin; xout = b%xout
    yin = b%yin; yout = b%yout

    ld=b%ld; lf=b%lf
    md=b%md; mf=b%mf

    dx = (xout - xin)/dble(b%lf_tot - b%ld_tot+1)
    dy = (yout - yin)/dble(b%mf_tot - b%md_tot+1)

    open(unit=20,file=trim(path_input)//"maillage.dat")
    read(20,*) taux
    close(20)

    do m= md-(nMF+1), mf+nMF
       do l = ld-(nMF+1), lf+nMF
          b%coord_nd(1,l,m)= dx*dble(l+b%l_charge(1,numproc)-1) + xin
          b%coord_nd(2,l,m)= dy*dble(m+b%m_charge(1,numproc)-1) + yin
       end do
    end do

    ldeb = ld+(nMF); lfin = lf-(nMF+1)
    mdeb = md+(nMF); mfin = mf-(nMF+1)
    if (b%m_charge(1,numproc) /= 1) mdeb = md-(nMF+1)
    if (b%l_charge(1,numproc) /= 1) ldeb = ld-(nMF+1)
    if (b%m_charge(2,numproc) /= b%mf_tot) mfin = mf+nMF
    if (b%l_charge(2,numproc) /= b%lf_tot) lfin = lf+nMF

    do m = mdeb, mfin
       do l = ldeb, lfin
          b%coord_nd(1,l,m)= dx*(dble(l+b%l_charge(1,numproc)-1)+taux*(rand()-.5d0)) + xin
          b%coord_nd(2,l,m)= dy*(dble(m+b%m_charge(1,numproc)-1)+taux*(rand()-.5d0)) + yin
       end do
    end do

  end subroutine maillage_random

  function rand()
    real*8 :: rand
    call random_number(rand)
  end function rand

!!$ Maillage carre (dx constant) (dy selon Kershaw 1/8)
  subroutine maillage_kershaw(b)
    type (STR_BLOC), pointer :: b

    ! variables locales
    integer :: l, m, lglob, mglob
    integer :: ld, lf, md, mf
    real*8  :: dx, dy, dy1, dy2, aa, bb
    real*8  :: x14, x34, x18, y14, y34, y18, y78

    ! Inclusion des mailles fictives via nMF ?

    ld=b%ld; lf=b%lf
    md=b%md; mf=b%mf

    dy2 = (b%xout - b%xin) / (0.5d0*dble(b%mf_tot - b%md_tot+1) + 1.d0)
    dy1 = 2.d0*dy2 / dble(b%mf_tot)

    dx = (b%xout - b%xin)/dble(b%lf_tot - b%ld_tot+1)
    dy = (b%yout - b%yin)/dble(b%mf_tot - b%md_tot+1)

    ! maillage cartesien de base pour les coordonnées en x
    do m= md-3, mf+2
       do l = ld-3, lf+2
          b%coord_nd(1,l,m)= dx*dble(l+b%l_charge(1,numproc)-1) + b%xin
          b%coord_nd(2,l,m)= 1.d60
       end do
    end do

    ! Kershaw 1/8 pour les coordonnées en y
    do m = md-3, mf+2
       mglob = m+b%m_charge(1,numproc)-1

       call calc_y14_y34_y18(mglob, b%mf_tot, dy1, dy2, y14, y34, y18, y78)
       x18 = dx*dble(b%lf_tot/8) + b%xin
       x14 = dx*dble(b%lf_tot/4) + b%xin
       x34 = dx*dble(3*b%lf_tot/4) + b%xin

       aa = (y14-y18) / (x14-x18)
       bb = (y34-y14) / (x34-x14)

       do l = ld-3, lf+2
          lglob = l+b%l_charge(1,numproc)-1
          if (lglob <= b%lf_tot/8) then
!!! zone 1
             b%coord_nd(2,l,m) = y18
          else if (lglob >= b%lf_tot/8+1 .and. lglob <= b%lf_tot/4-1) then
!!! zone 2
             b%coord_nd(2,l,m) = y18 + aa*(b%coord_nd(1,l,m)-x18)
          else if (lglob == b%lf_tot/4) then
!!! lf/4
             b%coord_nd(2,l,m) = y14
          else if (lglob >= b%lf_tot/4+1 .and. lglob <= 3*b%lf_tot/4-1) then
!!! zone 3
             b%coord_nd(2,l,m) = y14 + bb*(b%coord_nd(1,l,m)-x14)
          else if (lglob == 3*b%lf_tot/4) then
!!! 3/4 lf
             b%coord_nd(2,l,m) = y34
          else if (lglob >= 3*b%lf_tot/4+1 .and. lglob <= 7*b%lf_tot/8-1) then
!!! zone 4
             b%coord_nd(2,l,m) = y34 + aa*(b%coord_nd(1,l,m)-x34)
          else if (lglob >= 7*b%lf_tot/8) then
!!! zone 5
             b%coord_nd(2,l,m) = y78
          else
             print*, "bug maillage_kershaw", l, lglob, b%lf_tot
             call arret_code
             ! noeud deja calcule normalement
          end if
       end do
    end do

  contains
    subroutine calc_y14_y34_y18(m, mf, dy1, dy2, y14, y34, y18, y78)
      integer        , intent(in)  :: m, mf
      real*8         , intent(in)  :: dy1, dy2
      real*8         , intent(out) :: y14, y34, y18, y78

      if (m <= mf/2) then
         y14 = dble(m) * dy2
         y34 = dble(m) * dy1
         y18 = dble(m) * dy1
         y78 = dble(m) * dy2
      else
         y14 = 0.5d0 * mf*dy2 + dble(m-mf/2) * dy1
         y34 = (dble(m-mf/2)+1.d0) * dy2
         y18 = (dble(m-mf/2)+1.d0) * dy2
         y78 = 0.5d0 * mf*dy2 + dble(m-mf/2) * dy1
      end if
    end subroutine calc_y14_y34_y18

  end subroutine maillage_kershaw

end module generation_maillage
