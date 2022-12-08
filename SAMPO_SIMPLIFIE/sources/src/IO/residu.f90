module m_residu
  use m_struct
  use m_MPI
  use parametres_globaux, only : FLUIDE
  implicit none
  
  integer, parameter, private :: iunit_residu = 99
  integer, parameter, private :: iunit_volume = 100
  integer, parameter, private :: iunit_masse = 101
  integer, parameter, private :: iunit_conservation = 102
  integer, parameter, private :: iunit_dt = 103
  integer, parameter, private :: iunit_dilatation = 104
  integer, parameter, private :: iunit_energies = 105
  integer, parameter, private :: iunit_saut_pression = 106
  integer, parameter, private :: iunit_amplitude = 107
contains
  
  
  subroutine write_debug(tache, iter, temps)
    
    type (STR_TACHE), intent(in) ::  tache
    integer, intent(in) :: iter
    real*8, intent(in) :: temps

    integer :: i, isb
    type (STR_SUPER_BLOC), pointer :: sb

    do i = 1, tache%nb_membre 
       isb = tache%membre(i)
       sb => acces_super_bloc(isb)%pointeur
       if (sb%nature == FLUIDE)  call debug_fluide(sb, iter, temps)
    end do
    
    call check_dt_tache(iter, temps, tache%dt)

  end subroutine write_debug

  subroutine debug_fluide(sb, iter, temps)
    use variables_debug
    
    type (STR_SUPER_BLOC), pointer :: sb
    integer, intent(in) :: iter
    real*8, intent(in) :: temps
    
    if (debug_volume_bulle) call write_volume_bulle(sb, iter, temps)
    
    if (debug_conservation) call check_conservation(sb, iter)

  end subroutine debug_fluide

  function norme_Lp(n, X, p) result(norme_X)
    integer           :: n
    real*8            :: X(n)
    real*8            :: norme_X
    integer, optional :: p
    
    integer :: i

    if (.not. present(p)) p = 2    ! par defaut, on calcule une norme L2
           
    norme_X = 0.d0
    if (p /= 0) then
       do i = 1, n
          norme_X = norme_X + abs(X(i))**dble(p)
       end do
       norme_X = (norme_X)**(1.d0/dble(p))
    else ! norme L_infini
       do i = 1, n
          if (abs(X(i)) > norme_X) then
             norme_X = abs(X(i))
          end if
       end do
    end if

  end function norme_Lp
    
  subroutine norme(mpi_comm, n1, n2, n3, X, norme_X)
    integer                     , intent(in)  :: mpi_comm
    integer                     , intent(in)  :: n1, n2, n3
    real*8 , dimension(n1,n2,n3), intent(in)  :: X
    real*8 , dimension(n1)      , intent(out) :: norme_X
    
    real*8  :: p
    integer :: i, l, m, STATINFO
    real*8, dimension(n1) :: norme_proc

    p=2.d0

    norme_proc(:) = 0.d0
    do m = 1, n3
       do l = 1, n2
          do i = 1, n1
             norme_proc(i) = norme_proc(i) + abs(X(i,l,m))**p
          end do
       end do
    end do
    
    call MPI_ALLREDUCE(norme_proc, norme_X, n1, MPI_REAL8, &
         MPI_SUM, mpi_comm, STATINFO)
    
    norme_X(:) = (norme_X(:))**(1.d0/p)
    
  end subroutine norme

  subroutine calc_residu(mpi_comm,n1,n2,n3,X,X_old,norme_X0,residu)
    integer                    , intent(in)  :: mpi_comm
    integer                    , intent(in)  :: n1,n2,n3
    real*8, dimension(n1,n2,n3), intent(in)  :: X, X_old
    real*8, dimension(n1)      , intent(in)  :: norme_X0
    real*8, dimension(n1)      , intent(out) :: residu

    integer :: i
    real*8, dimension(n1) :: erreur
    real*8 :: denom

    erreur(:) = 0.d0; residu=1.d70

    call norme(mpi_comm,n1,n2,n3,X-X_old,erreur)

    do i = 1, n1
       if (norme_X0(i)<1.d-12) then
          denom=1.d0  ! pas terrible mais que faire
          erreur(i) = 0.d0
       else
          denom=norme_X0(i)
       end if
       erreur(i) = erreur(i) / denom
    end do
    
    ! pas besoin de mpi_allreduce car cela est fait dans norme
    residu(:) = erreur(:)
    
  end subroutine calc_residu

  
  subroutine calc_norme_sb(sb, neq, norme_sb)
    type (STR_SUPER_BLOC), pointer :: sb
    integer, intent(in) :: neq
    real*8, intent(out) :: norme_sb(neq)
    
    integer :: i, ib
    type (STR_BLOC), pointer :: b
    
    norme_sb(:) = 0.d0
    do i = 1, sb%nb_membre 
       ib = sb%membre(i)
       b => acces_bloc(ib)%pointeur
       if (bloc_est_local (ib)) then
          call norme(sb%mpi_comm,b%neq,b%lf,b%mf,b%U(:,b%ld:b%lf,b%md:b%mf),norme_sb)
          ! il y a deja un reduce dans le norme
       end if
    end do
    
  end subroutine calc_norme_sb

  subroutine calc_residu_sb(sb, neq, norme_U0, residu)
    use parametres_fluide, only : irhoE

    type (STR_SUPER_BLOC), pointer :: sb
    integer, intent(in)  :: neq
    real*8 , intent(in)  :: norme_U0(neq)
    real*8 , intent(out) :: residu(neq)
    
    integer :: i, ib
    type (STR_BLOC), pointer :: b
    
    ! calcul du residu |U_k+1 - U_k|/|U_0|
    do i = 1, sb%nb_membre 
       ib = sb%membre(i)
       b => acces_bloc(ib)%pointeur
       if (bloc_est_local (ib)) then
          call calc_residu(sb%mpi_comm,b%neq,b%lf,b%mf,b%U(:,b%ld:b%lf,b%md:b%mf), &
               & b%U_n(:,b%ld:b%lf,b%md:b%mf),norme_U0,residu)
          if (sb%eq_energie/=0) residu(irhoE)=0.d0
       end if
    end do
  end subroutine calc_residu_sb

  subroutine check_residu(iter, norme_U0, residu)
    use variables_globales, only : it_debug
    use parametres_fluide, only : EOS, CHIMIE_FIGEE
    integer, intent(in) :: iter
    real*8 , intent(in) :: norme_U0(:), residu(:)
    
    integer :: i
    
    if (numproc==0) then
       if (iter==0) then 
          ! pour l ecriture dans le fichier
          open(unit=iunit_residu,file="residu.dat")
          ! affichage a l ecran
          print*
          print*, "-------norme U0-------"
          ! do i = 1, ubound(norme_U0,1)
          !    write(*,'(i8,e13.5)') i, norme_U0(i)
          ! end do
          i = 1
          select case (fluide1%type)
          case(EOS)
             write(*,'(i8,e13.5,A)') i, norme_U0(i), "  rho y"; i=i+1
          case(CHIMIE_FIGEE)
             do i = 1, fluide1%ne
                write(*,'(i8,e13.5,A)') i, norme_U0(i), "  rho y c_"//fluide1%especes(i)%nom
             end do
             i = fluide1%ne + 1
          case default
             print*, 'type de fluide 1 non prévu'
          end select
          write(*,'(i8,e13.5,A)') i, norme_U0(i), "  rho (1-y)"; i=i+1         
          write(*,'(i8,e13.5,A)') i, norme_U0(i), "  rho u"; i=i+1         
          write(*,'(i8,e13.5,A)') i, norme_U0(i), "  rho v"; i=i+1
          write(*,'(i8,e13.5,A)') i, norme_U0(i), "  rho E"; i=i+1  
          write(*,'(i8,e13.5,A)') i, norme_U0(i), "  z"
          print*, "----------------------"
       else if (iter==1 .or. mod(iter, it_debug)==0) then
          write(iunit_residu,'(i8)', advance='no') iter
          do i = 1, ubound(norme_U0,1)
             write(iunit_residu,'(e13.5)',advance='no') real(residu(i))
          end do
          write(iunit_residu,*)
       end if
    end if
    
  end subroutine check_residu
  
  subroutine check_dt(sb, iter, temps)
    type (STR_SUPER_BLOC), pointer :: sb
    integer, intent(in) :: iter
    real*8, intent(in) :: temps

    if (numproc==0) then
       if (iter==0) then
          open(unit=iunit_dt,file="dt.dat", status='replace')
       else
          open(unit=iunit_dt,file="dt.dat", position='append', status='old')
          write(iunit_dt,*) iter, temps, sb%dt
          ! write(iunit_dt,*) iter-1, temps, sngl(sb%dt) ! pour comparer avec manu
       end if
       close(iunit_dt)
    end if
    
  end subroutine check_dt

  subroutine check_dt_tache(iter, temps, dt)
    integer, intent(in) :: iter
    real*8 , intent(in) :: temps, dt

    if (numproc==0) then
       if (iter==0) then
          open(unit=iunit_dt,file="dt.dat", status='replace')
       else
          open(unit=iunit_dt,file="dt.dat", position='append', status='old')
          write(iunit_dt,*) iter, temps, dt
       end if
       close(iunit_dt)
    end if
  end subroutine check_dt_tache
  

  subroutine check_conservation(sb, iter)
    type (STR_SUPER_BLOC), pointer :: sb
    integer, intent(in) :: iter
    
    integer :: iunit
    real*8, save :: conservation_ini(6)
    real*8, dimension(6) :: conservation
    
    iunit = iunit_conservation
    
    if (iter==0) then 
       call calc_conservation(sb, conservation_ini)
       if (numproc==0) then
          open(unit=iunit,file="conservation.dat", status='replace')
          write(iunit,*) "# iter ", "      rho  ", "       rho y  ", "      rho E"
          close(iunit)
       end if
    else
       
       call calc_conservation(sb, conservation)
       conservation(:) = abs(conservation(:) - conservation_ini(:)) / abs(conservation_ini(:))
       if (numproc==0) then
          open(unit=iunit,file="conservation.dat", position='append', status='old')
          write(iunit,1000) iter, conservation(1), conservation(2), conservation(5)
          close(iunit)
       end if
    end if

1000 format(i8,e13.5,e13.5,e13.5)
  end subroutine check_conservation

  subroutine calc_conservation(sb, conservation)
    
    type (STR_SUPER_BLOC), pointer :: sb
    real*8, dimension(6), intent(out) :: conservation

    integer :: i, ib, l, m, STATINFO
    real*8, dimension(6) :: conservation_proc
    type (STR_BLOC), pointer :: b

    if (mpi_comm_fluide == MPI_COMM_NULL) return
    
    conservation(:) = 0.d0
    conservation_proc(:) = 0.d0
    
    do i = 1, sb%nb_membre 
       ib = sb%membre(i)
       b => acces_bloc(ib)%pointeur
       if (bloc_est_local (ib)) then
          do m = b%md, b%mf
             do l = b%ld, b%lf
                conservation_proc(:) = conservation_proc(:)+b%U(:,l,m)
             end do
          end do
       end if
       call MPI_ALLREDUCE(conservation_proc, conservation, b%neq, MPI_REAL8, &
            MPI_SUM, sb%mpi_comm, STATINFO)
    end do
    
    ! pour eviter les divisions par zero
    do i = 1, 6
       conservation(i) = max(conservation(i), 1.d-15)
    end do
  end subroutine calc_conservation

  subroutine calc_masse(sb, masse, masse_new)
    
    type (STR_SUPER_BLOC), pointer :: sb
    real*8, intent(out) :: masse, masse_new

    integer :: i, ib, l, m, STATINFO
    real*8 :: masse_proc, masse_new_proc
    type (STR_BLOC), pointer :: b

    masse = 0.d0
    masse_proc = 0.d0
    masse_new = 0.d0
    masse_new_proc = 0.d0

    do i = 1, sb%nb_membre 
       ib = sb%membre(i)
       b => acces_bloc(ib)%pointeur
       if (bloc_est_local (ib)) then
          do m = b%md, b%mf
             do l = b%ld, b%lf
                masse_proc = masse_proc+b%U_n(1,l,m)*b%aire(l,m)
                masse_new_proc = masse_new_proc+b%U(1,l,m)*b%aire_new(l,m) 
             end do
          end do
       end if
       call MPI_ALLREDUCE(masse_proc, masse, 1, MPI_REAL8, &
            MPI_SUM, sb%mpi_comm, STATINFO)
       call MPI_ALLREDUCE(masse_new_proc, masse_new, 1, MPI_REAL8, &
            MPI_SUM, sb%mpi_comm, STATINFO)
    end do
     
  end subroutine calc_masse

!!! Verification de la dilatation du fluide dans le cas de la fusion 
!!! ( rho_0 V - \int_{\Omega} rho ) / (rho_0*V)
  subroutine check_dilatation(sb, iter)
    type (STR_SUPER_BLOC), pointer :: sb
    integer, intent(in) :: iter

    integer      :: iunit
    real*8, save :: masse_ini
    real*8       :: masse_new, masse
        
    iunit = iunit_dilatation
    
    if (iter==0) then
       call calc_masse(sb, masse, masse_ini)
       if (numproc==0) open(unit=iunit,file="dilatation.dat", status='replace')
    else
       call calc_masse(sb, masse, masse_new)
       
       ! ecriture
       if (numproc==0) then
          open(unit=iunit,file="dilatation.dat", position='append', status='old')
          write(iunit,*) iter, real((masse-masse_ini)/masse_ini)
          close(iunit)
       end if
       
    end if
    
  end subroutine check_dilatation
  
  subroutine write_volume_bulle(sb, iter, temps)
    type (STR_SUPER_BLOC), intent(in) :: sb
    integer, intent(in) :: iter
    real*8, intent(in) :: temps
    
    real*8 :: volume1, volume2, volume3
    type (STR_BLOC), pointer :: b
    integer :: i, ib

    do i = 1, sb%nb_membre 
       ib = sb%membre(i)
       b => acces_bloc(ib)%pointeur
       if (bloc_est_local (ib)) then
          call volume_bulle(b, volume1, 0.1d0)
          call volume_bulle(b, volume2, 0.5d0)
          call volume_bulle(b, volume3, 0.9d0)
       end if
    end do

    if (iter == 0 .and. numproc==0) then
       print*, "creation du fichier volume.dat" 
       open(unit=iunit_volume, file="volume.dat")
    end if

    if (numproc==0) then
       write(iunit_volume,*) real(temps), real(acos(-1.d0)*1.d-2-temps**2/120.d0), &
            & real(volume1), real(volume2), real(volume3)
    end if

  end subroutine write_volume_bulle

  subroutine volume_bulle(b, volume, seuil)
    use m_struct
    use m_MPI
    implicit none

    type (STR_BLOC), intent(in) :: b
    real*8, intent(in) :: seuil
    real*8, intent(out) :: volume

    integer :: l, m, STATINFO
    real*8 :: vol

    vol = 0.d0
    do m = b%md, b%mf
       do l = b%ld, b%lf
          if (b%z(l,m)<seuil) then
             vol = vol + b%aire(l,m)
          end if
       end do
    end do

    call MPI_ALLREDUCE(vol, volume, 1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,STATINFO)

    if ( trim(b%face(3)%patch(1)%bc%condition) == "symetrie" .or. &
         trim(b%face(4)%patch(1)%bc%condition) == "symetrie" ) then
       volume = 2.d0*volume
    end if

  end subroutine volume_bulle

end module m_residu
