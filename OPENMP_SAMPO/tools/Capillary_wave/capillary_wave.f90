program capillary_wave
	implicit none

	integer :: iunit_poly = 20
	real*8, parameter :: pi = 4.d0*atan(1.d0)

	real*8 :: xd    = -0.5d0
	real*8 :: xf    = 0.5d0
	real*8 :: rho0  = 1.d3
	real*8 :: rho1  = 1.d3
	real*8 :: nu    = 1.d-3
	real*8 :: sigma = 1.d-1
	real*8 :: L     = 2.d0*pi
	real*8 :: omega0, b, k, epsilon

	real*8 :: x
	integer :: i, n=100000

	k      = 2.d0*pi/L

	! ***** version avec nu et k *****
	omega0 = sqrt((sigma*k**3)/(rho0+rho1))
	b      = (rho0*rho1)/(rho0+rho1)**2

	print*, 'omega0^2 = ', omega0**2
	print*, 'b        = ', b
	print*, 'nu       = ', nu
	print*, 'k        = ', k

	open(unit=iunit_poly,file="polynome.dat")
    write(iunit_poly,*) "#  x        y"
	do i=0,n
		x = xd + i*(xf-xd)/dble(n)
		write(iunit_poly,*) x, polynome(b,k,nu,omega0,x)
	end do
    close(iunit_poly)



    ! ***** version avec epsilon *****
	! epsilon = 6.472d-2
	! omega0  = sqrt(6.778d0)
	! b      = (rho0*rho1)/(rho0+rho1)**2

	! open(unit=iunit_poly,file="polynome.dat")
 !    write(iunit_poly,*) "#x ", ", ", "y"
	! do i=0,n
	! 	x = xd + i*(xf-xd)/dble(n)
	! 	write(iunit_poly,*) x, ",", polynome2(b,epsilon,omega0,x)
	! end do
 !    close(iunit_poly)

contains


	function polynome(b,k,nu,omega0,x)
		real*8 :: b,k,nu,omega0,x
		real*8 :: polynome
		real*8 :: gamma

		gamma = nu*(k**2)

		polynome = x**4 &
		- 4.d0*b*sqrt(gamma)*x**3 &
		+ 2.d0*(1.d0-6.d0*b)*gamma*x**2 &
		+ 4.d0*(1.d0-3.d0*b)*(gamma)**(3.d0/2.d0)*x &
		+ (1.d0-4.d0*b)*(gamma)**2 + omega0**2
	end function

	function polynome2(b,epsilon,omega0,x)
		real*8 :: b,epsilon,omega0,x
		real*8 :: polynome2, gamma

		gamma = epsilon*omega0

		polynome2 = x**4 &
		- 4.d0*b*sqrt(gamma)*x**3 &
		+ 2.d0*(1.d0-6.d0*b)*gamma*x**2 &
		+ 4.d0*(1.d0-3.d0*b)*(gamma)**(3.d0/2.d0)*x &
		+ (1.d0-4.d0*b)*(gamma)**2 + omega0**2
	end function
end program
