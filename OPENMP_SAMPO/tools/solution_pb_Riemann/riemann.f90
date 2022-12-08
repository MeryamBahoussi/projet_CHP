program riemann
  implicit none
  double precision :: P_R,rho_R,gamma_R,P_inf_R,c_R, u_R
  double precision :: P_L,rho_L,gamma_L,P_inf_L,c_L, U_L 
  double precision :: P_star, u_star
  integer :: k,i
  double precision, dimension(:), allocatable :: rho, pression, vitesse, x_sortie, x
  double precision :: time, x_int
  integer :: Nb_pt
  character(len=40) :: test = "sod_LG"
  double precision :: long = 1.d0, xmin=0.d0, x_contact

  select case( test )

  case("sod_1d")
     P_L = 1.d0 
     rho_L = 1.d0
     gamma_L = 1.4d0
     P_inf_L = 0.d0
     c_L = sqrt( gamma_L * (P_L+P_inf_L) / rho_L )
     U_L = 0.d0

     P_R = 1.d-1
     rho_R =  1.25d-1
     gamma_R = 2.4d0
     P_inf_R = 0.d0
     c_R = sqrt( gamma_R * (P_R+P_inf_R)/rho_R )
     u_R = 0.d0

     time = 1.4d-1
     long = 1.d0
     x_int = 0.5d0
  case("sod_LG")
     P_L = 1.d9
     rho_L = 1.d3
     gamma_L = 4.4d0
     P_inf_L = 6.d8
     c_L = sqrt( gamma_L * (P_L+P_inf_L) / rho_L )
     U_L = 0.d0

     P_R = 1.d5
     rho_R =  5.d0
     gamma_R = 1.4d0
     P_inf_R = 0.d0
     c_R = sqrt( gamma_R * (P_R+P_inf_R)/rho_R )
     u_R = 0.d0

     time = 2.4d-4
     long = 1.d0
     x_int = 0.7d0
  case("sod_LeBlanc")
     P_L = 0.666666666d-1 
     rho_L = 1.d0
     gamma_L = 1.666666666d0
     P_inf_L = 0.d0
     c_L = sqrt( gamma_L * (P_L+P_inf_L) / rho_L )
     U_L = 0.d0

     P_R = 0.666666666d-10
     rho_R =  1.d-3
     gamma_R = 1.666666666d0
     P_inf_R = 0.d0
     c_R = sqrt( gamma_R * (P_R+P_inf_R)/rho_R )
     u_R = 0.d0

     time = 6.d0
     long = 9.d0
     x_int = 3.d0

  case("sod_123problem")
     P_L = 0.4d0 
     rho_L = 1.d0
     gamma_L = 1.4d0
     P_inf_L = 0.d0
     c_L = sqrt( gamma_L * (P_L+P_inf_L) / rho_L )
     U_L = -2.d0

     P_R = 0.4d0
     rho_R =  1.d0
     gamma_R = 1.4d0
     P_inf_R = 0.d0
     c_R = sqrt( gamma_R * (P_R+P_inf_R)/rho_R )
     u_R = 2.d0

     time = 1.d0
     long = 8.d0
     xmin =-4.d0
     x_int = 0.d0

  case("sod_2")
     P_L = 1.d5 
     rho_L = 1.d0
     gamma_L = 1.4d0
     P_inf_L = 0.d0
     c_L = sqrt( gamma_L * (P_L+P_inf_L) / rho_L )
     U_L = 0.d0

     P_R =  1.d4
     rho_R =  0.125d0
     gamma_R = 1.4d0
     P_inf_R = 0.d0
     c_R = sqrt( gamma_R * (P_R+P_inf_R)/rho_R )
     u_R = 0.d0

     time = 3.1d-4
     x_int = .5d0

  case("iollo_9")
     P_L = 1.d9 
     rho_L = 1000.d0
     gamma_L = 4.4d0
     P_inf_L = 6.d8
     c_L = sqrt( gamma_L * (P_L+P_inf_L) / rho_L )
     U_L = 0.0d0

     P_R =  1.0d5
     rho_R =  50.0d0
     gamma_R = 1.4d0
     P_inf_R = 0.d0
     c_R = sqrt( gamma_R * (P_R+P_inf_R)/rho_R )
     u_R = 0.d0

     time = 0.24d-3 !24.d-6
     x_int = .7d0
      
  case("injection")
     P_L = 1.0d5 
     rho_L = 10.d0
     gamma_L = 4.4d0
     P_inf_L = 6.d8
     c_L = sqrt( gamma_L * (P_L+P_inf_L) / rho_L )
     U_L = 10.0d0

     P_R =  1.0d5
     rho_R =  1.0d0
     gamma_R = 1.4d0
     P_inf_R = 0.d0
     c_R = sqrt( gamma_R * (P_R+P_inf_R)/rho_R )
     u_R = 100.d0

     time = 1.6d-3
     x_int = .001d0	

  case("manup156")
     P_L = 1.0d9 
     rho_L = 1000.d0
     gamma_L = 4.4d0
     P_inf_L = 6.d8
     c_L = sqrt( gamma_L * (P_L+P_inf_L) / rho_L )
     U_L = 0.0d0

     P_R =  1.0d5
     rho_R =  50.0d0
     gamma_R = 1.4d0
     P_inf_R = 0.d0
     c_R = sqrt( gamma_R * (P_R+P_inf_R)/rho_R )
     u_R = 0.d0

     time = 2.4d-4
     x_int = .7d0
  
  case("manup158")
     P_L = 1.0d12 
     rho_L = 1000.d0
     gamma_L = 4.4d0
     P_inf_L = 6.d8
     c_L = sqrt( gamma_L * (P_L+P_inf_L) / rho_L )
     U_L = 0.0d0

     P_R =  1.0d5
     rho_R =  10.0d0
     gamma_R = 1.4d0
     P_inf_R = 0.d0
     c_R = sqrt( gamma_R * (P_R+P_inf_R)/rho_R )
     u_R = 0.d0

     time = 8.3d-6
     x_int = .6d0

  end select

  print*,"quel est le nombre de points "
  read*, Nb_pt
  !print*,"le temps final"
  !read*,time


  allocate( rho(Nb_pt), pression(Nb_pt), vitesse(Nb_pt), x(Nb_pt), x_sortie(Nb_pt) )
  rho = 0.d0
  pression = 0.d0
  vitesse = 0.d0
  x = 0.d0
  x_sortie = 0.d0

  call recherche_p_star_u_star(P_R,rho_R,gamma_R,P_inf_R,c_R, u_R,P_L,rho_L,gamma_L,P_inf_L,c_L, U_L, U_star, P_star)

  print*,"fin", P_star, u_star

  do i = 1, Nb_pt
     x(i) = xmin + i * long / float(Nb_pt)
  end do

  x_sortie = x

  x = x - x_int


  call solution_riemann(P_R,rho_R,gamma_R,P_inf_R,c_R, u_R,P_L,rho_L,gamma_L,P_inf_L,c_L, U_L, U_star, P_star,rho, pression, vitesse, x,time ,Nb_pt )

  ! open(unit=1, file="solution_riemann.txt")
  ! write(1,*) 'x', ',' , 'rho' , ',' , 'vitesse', ',' , 'pression'
  ! do i = 1, Nb_pt
  !    write(1,'(e12.4,A,e12.4,A,e12.4,A,e12.4)') x_sortie(i),',',rho(i),',',vitesse(i),',',pression(i)
  ! end do
  ! close(1)

  open(unit=1, file="pression.dat")
  do i = 1, Nb_pt
     write(1,*) x_sortie(i), pression(i)
  end do
  close(1)

  open(unit=1, file="vitesse.dat")
  do i = 1, Nb_pt
     write(1,*) x_sortie(i), vitesse(i)
  end do
  close(1)	

  open(unit=1, file="rho.dat")
  do i = 1, Nb_pt
     write(1,*) x_sortie(i), rho(i)
  end do
  close(1)

  x_contact = x_int + time*u_star
  print*, 'position du contact', x_contact

  open(unit=1,file="z.dat")
  do i = 1, Nb_pt
     if(x_sortie(i)<=x_contact)then
        write(1,*) x_sortie(i), 1.d0
     else
        write(1,*) x_sortie(i), 0.d0
     end if
  end do
  print*,"FIN FIN FIN"

contains

  subroutine solution_riemann(P_R,rho_R,gamma_R,P_inf_R,c_R, u_R,P_L,rho_L,gamma_L,P_inf_L,c_L, U_L, U_star, P_star,rho, pression, vitesse, x,time ,Nb_pt )
    implicit none
    double precision, intent(in) :: P_R,rho_R,gamma_R,P_inf_R,c_R, u_R
    double precision, intent(in)  :: P_L,rho_L,gamma_L,P_inf_L,c_L, U_L
    double precision, intent(in) :: U_star, P_star	
    integer, intent(in) :: Nb_pt
    double precision, dimension(1:Nb_pt), intent(out) :: rho, pression, vitesse
    double precision, dimension(1:Nb_pt), intent(in) :: x
    double precision, intent(in) :: time
    double precision :: vitesse_choc, a, r, c_star, rho_star, x_t
    integer :: i

    do i = 1, Nb_pt

       ! Le point est du coté gauche
       !##################################
       if ( x(i) / time <= U_star ) then 
          !##################################

          !----------------------------------------
          if ( P_L < P_star ) then ! choc
             !----------------------------------------
             vitesse_choc = U_L - c_L * sqrt( (gamma_L+1.d0)/(2.d0*gamma_L) * (P_star+P_inf_L)/(P_L+P_inf_L) + (gamma_L-1.d0)/(2.d0*gamma_L) )

             if ( x(i) / time <= vitesse_choc ) then
                rho(i) = rho_L
                pression(i) = P_L
                vitesse(i) = U_L
             else 
                r = (P_star+P_inf_L)/(P_L+P_inf_L)
                a = (gamma_L+1.d0)/(gamma_L-1.d0)
                rho(i) = rho_L * (a * r + 1.d0) / (r+a) 
                pression(i) = P_star
                vitesse(i) = U_star				
             end if

             !----------------------------------------
          else ! détente
             !----------------------------------------

             r = (P_star+P_inf_L)/(P_L+P_inf_L)
             !rho_star = rho_L * 
             c_star = c_L + (gamma_L-1.d0)/2.d0 * (U_L-U_star)

             if ( x(i) / time <= u_L - c_L ) then 
                rho(i) = rho_L
                pression(i) = P_L
                vitesse(i) = U_L								
             else if ( x(i) / time >= u_star - c_star ) then
                rho(i) = rho_L * r**(1.d0/gamma_L)
                pression(i) = P_star
                vitesse(i) = U_star				
             else
                x_t = x(i) / time
                rho(i) = rho_L * (2.d0/(gamma_L+1)+(gamma_L-1.d0)/((gamma_L+1.d0)*c_L) * (u_L-x_t))**(2.d0/(gamma_L-1.d0))
                pression(i) = (P_L+P_inf_L) * (2.d0/(gamma_L+1)+(gamma_L-1)/((gamma_L+1)*c_L)*(u_L-x_t))**(2.d0*gamma_L/(gamma_L-1.d0)) - P_inf_L
                vitesse(i) = 2.d0/(gamma_L+1) * (c_L+(gamma_L-1.d0)/2.d0*U_L+x_t)			
             end if

          end if

          ! Le point est du coté droit
          !##################################
       else
          !##################################

          !----------------------------------------
          if ( P_R < P_star ) then ! choc
             !----------------------------------------
             vitesse_choc = U_R + c_R * sqrt( (gamma_R+1.d0)/(2.d0*gamma_R) * (P_star+P_inf_R)/(P_R+P_inf_R) + (gamma_R-1.d0)/(2.d0*gamma_R) )
             if ( x(i) / time >= vitesse_choc ) then
                rho(i) = rho_R
                pression(i) = P_R
                vitesse(i) = U_R
             else 
                r = (P_star+P_inf_R)/(P_R+P_inf_R)
                a = (gamma_R+1.d0)/(gamma_R-1.d0)
                rho(i) = rho_R * (a * r + 1.d0) / (r+a) 
                pression(i) = P_star
                vitesse(i) = U_star				
             end if

             !----------------------------------------
          else ! detente
             !----------------------------------------

             r = (P_star+P_inf_R)/(P_R+P_inf_R)
             !rho_star = rho_L * 
             c_star = c_R + (gamma_R-1.d0)/2.d0 * (U_star-U_R)
             
          
             if ( x(i) / time >= u_R + c_R ) then 
                rho(i) = rho_R
                pression(i) = P_R
                vitesse(i) = U_R								
             else if ( x(i) / time <= u_star + c_star ) then
                rho(i) = rho_R * r**(1.d0/gamma_R)
                pression(i) = P_star
                vitesse(i) = U_star				
             else
                x_t = x(i) / time
                rho(i) = rho_R * (2.d0/(gamma_R+1)-(gamma_R-1.d0)/((gamma_R+1.d0)*c_R) * (u_R-x_t))**(2.d0/(gamma_R-1.d0))
                pression(i) = (P_R+P_inf_R) * (2.d0/(gamma_R+1)-(gamma_R-1)/((gamma_R+1)*c_R)*(u_R-x_t))**(2.d0*gamma_R/(gamma_R-1.d0)) - P_inf_R
                vitesse(i) = 2.d0/(gamma_R+1) * (-c_R+(gamma_R-1.d0)/2.d0*U_R+x_t)			
             end if

          end if
       end if

    end do

  end subroutine solution_riemann




  subroutine recherche_p_star_u_star(P_R,rho_R,gamma_R,P_inf_R,c_R, u_R,P_L,rho_L,gamma_L,P_inf_L,c_L, U_L, U_star, P_star)
    implicit none
    double precision, intent(in) :: P_R,rho_R,gamma_R,P_inf_R,c_R, u_R
    double precision, intent(in)  :: P_L,rho_L,gamma_L,P_inf_L,c_L, U_L
    double precision, intent(out) :: U_star, P_star
    double precision :: f_L, f_R, f, P_star_min, P_star_max, U_star_R,U_star_L 
    integer :: k,i


    P_star_min = min(P_R,P_L)
    call f_de_P(P_star_min,u_L,P_L,rho_L,gamma_L,P_inf_L,c_L,U_R,P_R,rho_R,gamma_R,P_inf_R,c_R,f)
    do while ( f >= 0.d0 )
       P_star_min = P_star_min / 2.d0
       call f_de_P(P_star_min,u_L,P_L,rho_L,gamma_L,P_inf_L,c_L,U_R,P_R,rho_R,gamma_R,P_inf_R,c_R,f)
    end do

    P_star_max = max(P_R,P_L)
    call f_de_P(P_star_max,u_L,P_L,rho_L,gamma_L,P_inf_L,c_L,U_R,P_R,rho_R,gamma_R,P_inf_R,c_R,f)
    do while ( f <= 0.d0 )
       P_star_max = P_star_max * 2.d0
       call f_de_P(P_star_max,u_L,P_L,rho_L,gamma_L,P_inf_L,c_L,U_R,P_R,rho_R,gamma_R,P_inf_R,c_R,f)
    end do

    do k = 1, 22000

       P_star = 0.5d0 *( P_star_min+ P_star_max)

       call f_de_P(P_star,u_L,P_L,rho_L,gamma_L,P_inf_L,c_L,U_R,P_R,rho_R,gamma_R,P_inf_R,c_R,f)


       if ( f < 0.d0)  then 
          P_star_min  = P_star
       else if (f > 0.d0 ) then 
          P_star_max = P_star
       end if

       !print*,"p_star", k, P_star 

    end do

    call onde(P_star,P_R,rho_R,gamma_R,P_inf_R,c_R,f_R)
    U_star_R = U_R + f_R

    call onde(P_star,P_L,rho_L,gamma_L,P_inf_L,c_L,f_L)
    U_star_L = U_L - f_L

    print*,"**",  U_star_L, U_star_R

    U_star = U_star_L

  end subroutine recherche_p_star_u_star



  subroutine f_de_P(P_star,u_L,P_L,rho_L,gamma_L,P_inf_L,c_L,U_R,P_R,rho_R,gamma_R,P_inf_R,c_R,f)
    implicit none 
    double precision :: P_R,rho_R,gamma_R,P_inf_R,c_R, u_R
    double precision :: P_L,rho_L,gamma_L,P_inf_L,c_L, U_L
    double precision :: P_star, f_L, f_R, f	

    call onde(P_star,P_L,rho_L,gamma_L,P_inf_L,c_L,f_L)
    call onde(P_star,P_R,rho_R,gamma_R,P_inf_R,c_R,f_R)

    f = u_R-U_L + f_R + f_L

  end subroutine f_de_P

  subroutine  onde(P_star,P,rho,gamma,P_inf,c,f)
    implicit none
    double precision, intent(in) :: P_star,P,rho,gamma,P_inf,c
    double precision, intent(out) :: f

    double precision :: r, P_chap_star, P_chap


    r = (P_star+P_inf)/(P+P_inf)
    P_chap_star = P_star+P_inf
    P_chap = P +P_inf

    If ( P_star > P) then 

       f = sqrt( 2.d0 / (  rho * P_chap * ( (gamma+1.d0) * r + (gamma-1.d0) ) ) ) * (P_star-p)

    else
       f = 2.d0 * c /(gamma-1.d0) * (r**((gamma-1.d0)/2.d0/gamma) -1.d0 )

    end if

  end subroutine  onde

end program riemann
