module delta_t
  use m_struct
  use m_MPI
  use m_temps
  use parametres_globaux, only : FLUIDE
  use m_boundary_condition
  use outils_pas_acoustique

  implicit none

contains

  ! calcul du pas de temps
  subroutine deltat(tache, iter, arret, dt)
    use variables_globales, only : it_out

    type (STR_TACHE), intent(inout) :: tache
    integer         , intent(in)    :: iter
    logical         , intent(out)   :: arret
    real*8          , intent(inout) :: dt
    
    integer :: i, isb 
    real*8  :: dt_fld
    real*8  :: tnormal, t_objectif
    type (STR_SUPER_BLOC), pointer :: sb

    call debut_watchTime(wT_deltat)

    do i = 1, tache%nb_membre 
       isb = tache%membre(i)
       sb => acces_super_bloc(isb)%pointeur 
       
       !sb%cfl_ini = sb%cfl
       if (sb%montee_cfl%actif) call montee_en_cfl(iter, sb%montee_cfl, sb%cfl)

       call dt_fluide(sb, iter, dt_fld)
    end do

    dt = dt_fld
    
    do isb = 1, size(acces_super_bloc,1)
       acces_super_bloc(isb)%pointeur%dt=dt
    end do
    
!!!! Recherche du prochain temps que l'on veut atteindre
    if (tache%type_freqout == "iter") then
       t_objectif = tache%tmax
    else if (tache%type_freqout == "temps") then
       t_objectif = tache%t_imposes(tache%ind_timp)
    end if
    
    tnormal = temps + dt

    if(tnormal >= t_objectif) then
       dt = t_objectif - temps
       temps = t_objectif
       tache%ind_timp = tache%ind_timp+1
    else
       temps = tnormal 
    end if
    
    if (dt<=0.d0) then
       print*, "Erreur pas de temps negatif ou nul", dt
       call arret_code
    end if

    arret = .false.
    if (temps == tache%tmax) then
       arret = .true.
       tache%ind_timp = tache%ind_timp - 1
    end if

!!! debug : affichage des differents dt pendant le calcul
    if (numproc==0 .and. (mod(iter, it_out)==0 .or. iter==1 .or. temps==t_objectif)) then
       write(6,1002) iter, temps, dt
    end if

1002 format(' ITER=',i8,'   ',' Temps=',e13.5,'   ',' DT=',e13.5)

    call fin_watchTime(wT_deltat)
  end subroutine deltat
  
  subroutine dt_fluide(sb, iter, dt_fld)
    implicit none

    type (STR_SUPER_BLOC), pointer :: sb
    integer, intent(in) :: iter
    real*8, intent(out) :: dt_fld 

    real*8 :: dt_e, dt_v

    dt_e = 1.d70    ! dt pour le Euler direct
    dt_v = 1.d70    ! dt pour les termes visqueux
    dt_fld=1.d70    ! pas de temps du fluide

    if (mpi_comm_fluide == MPI_COMM_NULL) return

    ! decodage des points interieurs puis calcul des mailles fictives
    call debut_watchTime(wT_dt_decode)
    call decode(sb, decode_MF=.false.)
    call fin_watchTime(wT_dt_decode)

    call debut_watchTime(wT_dt_CL)
    call condition_limites_Euler(sb)
    call fin_watchTime(wT_dt_CL)

    ! calcul de la vitesse du son Lagrangienne
    call debut_watchTime(wT_dt_pente)
    call pentes_solveur_Riemann(sb)
    call fin_watchTime(wT_dt_pente)

    call debut_watchTime(wT_dt_Euler)
    call dt_euler(sb, dt_e)
    call fin_watchTime(wT_dt_Euler)

    if (sb%visqueux) call dt_visqueux(sb, dt_v)
       
    dt_fld = min(dt_e, dt_v)
       
!!! debug : affichage des differents dt pendant le calcul
    if (numproc==0) then
       if (.not. sb%sol_stat%actif .and. (mod(iter, it_out)==0 .or. iter==1)) then
          print*
          write(6,1001) dt_e, dt_v
       end if
    end if
        
1001 format(' dt euler =',e13.5,'   ',' dt visqueux =',e13.5)
  end subroutine dt_fluide
  
  subroutine dt_euler(sb, dt)
    type (STR_SUPER_BLOC), pointer     :: sb
    real*8               , intent(out) :: dt
    
    integer :: l, m, i, ib
    integer :: ld, lf, md, mf
    integer :: STATINFO
    real*8  :: invdt, invdt_MPI
    real*8  :: min_dl_l, min_dl_m, U(2), n_dir_l(2), n_dir_m(2)
    real*8  :: Un_x1, Un_x2
    real*8  :: Un_y1, Un_y2
    real*8  :: vx, vy
    real*8, dimension(:,:), pointer :: vmesh_x, vmesh_y
    type (STR_BLOC), pointer :: b

    invdt = -1.d70

    do i = 1, sb%nb_membre 
       ib = sb%membre(i)
       b => acces_bloc(ib)%pointeur
       
       if (bloc_est_local (ib)) then
          ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
          
          vmesh_x => b%zeros
          vmesh_y => b%zeros

          select case (sb%solveur_Riemann)
          case(GALLICE)
             do m = md, mf
                do l = ld, lf
                   min_dl_l = min(b%dl_l(l-1,m),b%dl_l(l,m))
                   min_dl_m = min(b%dl_m(l,m-1),b%dl_m(l,m))
                   U = (/ b%u_x(l,m), b%u_y(l,m) /)
                   
                   Un_x1 = dot_product(U,b%n_dir_l(:,l-1,m))
                   Un_x2 = dot_product(U,b%n_dir_l(:,l  ,m))
                   Un_y1 = dot_product(U,b%n_dir_m(:,l,m-1))
                   Un_y2 = dot_product(U,b%n_dir_m(:,l,m  ))
                   
                   vx= Max(abs(vmesh_x(l-1,m)) + abs(Un_x1) + b%Cp_l(l-1,m)/b%rho(l,m),&
                        &  abs(vmesh_x(l  ,m)) + abs(Un_x2) + b%Cm_l(l  ,m)/b%rho(l,m) )
                   vy= Max(abs(vmesh_y(l,m-1)) + abs(Un_y1) + b%Cp_m(l,m-1)/b%rho(l,m),&
                        &  abs(vmesh_y(l,m  )) + abs(Un_y2) + b%Cm_m(l,m  )/b%rho(l,m) )
                   
                   invdt = max( invdt, vx/min_dl_m ) ! min_dl_m = delta_x
                   invdt = max( invdt, vy/min_dl_l ) ! min_dl_l = delta_y
                end do
             end do
          case(HLL, ROE)
             do m = md, mf
                do l = ld, lf
                   min_dl_l = min(b%dl_l(l-1,m),b%dl_l(l,m))
                   min_dl_m = min(b%dl_m(l,m-1),b%dl_m(l,m))
                   U = (/ b%u_x(l,m), b%u_y(l,m) /)
                   
                   n_dir_l = .5d0 * (b%n_dir_l(:,l-1,m)+b%n_dir_l(:,l,m))
                   n_dir_m = .5d0 * (b%n_dir_m(:,l,m-1)+b%n_dir_m(:,l,m))
                   
                   vx = b%c(l,m) + abs(dot_product(U,n_dir_l)) + .5d0*(vmesh_x(l-1,m)+vmesh_x(l,m))
                   vy = b%c(l,m) + abs(dot_product(U,n_dir_m)) + .5d0*(vmesh_y(l,m-1)+vmesh_y(l,m))
                
                   invdt = max( invdt, vx/min_dl_m ) ! min_dl_m = delta_x
                   invdt = max( invdt, vy/min_dl_l ) ! min_dl_l = delta_y
                end do
             end do
          case default
             print*, "dt fluide non defini pour le solveur de Riemann choisi"
             call arret_code
          end select
       end if
    end do
    
    call MPI_ALLREDUCE(invdt,invdt_mpi,1,mpi_real8,mpi_MAX,sb%mpi_comm,statinfo)
    dt = sb%cfl/invdt_MPI
    
  end subroutine dt_euler

  subroutine dt_visqueux(sb, dt_v)
    use m_struct
    use m_MPI
    implicit none

    type (STR_SUPER_BLOC), intent(in) ::  sb
    real*8, intent(out) :: dt_v

    integer :: l, m, i, ib
    integer :: ld, lf, md, mf
    integer :: STATINFO
    real*8 :: invdt_mu, invdt_K, invdt_D, invdt, invdt_MPI
    
    real*8, dimension(:,:), pointer :: rho, dl_l, dl_m, Cv, mu, K, D
    type (STR_BLOC), pointer ::  b

    invdt = -1.d70
    
    do i = 1, sb%nb_membre 
       ib = sb%membre(i)
       b => acces_bloc(ib)%pointeur

       if (bloc_est_local (ib)) then
          ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
          rho => b%rho
          Cv => b%Cv
          mu => b%mu
          K => b%K
          D => b%D
          dl_l => b%dl_l; dl_m => b%dl_m
          
          if (fluide1%type /= CHIMIE_FIGEE) D(:,:) = 0.d0

          do m = md, mf
             do l = ld, lf
                ! contrainte liee a la viscosite
                invdt_mu = -1.d70
                invdt_mu = max(invdt_mu, .5d0*(mu(l  ,m  )+mu(l+1,m  ))/dl_l(l  ,m)**2)
                invdt_mu = max(invdt_mu, .5d0*(mu(l-1,m  )+mu(l  ,m  ))/dl_l(l-1,m)**2)
                invdt_mu = max(invdt_mu, .5d0*(mu(l  ,m  )+mu(l  ,m+1))/dl_m(l  ,m)**2)
                invdt_mu = max(invdt_mu, .5d0*(mu(l  ,m-1)+mu(l  ,m  ))/dl_m(l,m-1)**2)
                invdt_mu = 4.d0 * invdt_mu / rho(l,m) 
                
                ! contrainte liee a la conduction thermique
                invdt_K = -1.d70 
                invdt_K = max(invdt_K, .5d0*(K(l  ,m  )+K(l+1,m  ))/dl_l(l  ,m  )**2)
                invdt_K = max(invdt_K, .5d0*(K(l-1,m  )+K(l  ,m  ))/dl_l(l-1,m  )**2)
                invdt_K = max(invdt_K, .5d0*(K(l  ,m  )+K(l  ,m+1))/dl_m(l  ,m  )**2)
                invdt_K = max(invdt_K, .5d0*(K(l  ,m-1)+K(l  ,m  ))/dl_m(l  ,m-1)**2)
                invdt_K = 4.d0 * invdt_K / (rho(l,m)*Cv(l,m))
                
                ! contrainte liee a la diffusion de Fick
                invdt_D = -1.d70 
                invdt_D = max(invdt_D, .5d0*(D(l  ,m  )+D(l+1,m  ))/dl_l(l  ,m  )**2)
                invdt_D = max(invdt_D, .5d0*(D(l-1,m  )+D(l  ,m  ))/dl_l(l-1,m  )**2)
                invdt_D = max(invdt_D, .5d0*(D(l  ,m  )+D(l  ,m+1))/dl_m(l  ,m  )**2)
                invdt_D = max(invdt_D, .5d0*(D(l  ,m-1)+D(l  ,m  ))/dl_m(l  ,m-1)**2)
                invdt_D = 4.d0 * invdt_D
                
                invdt = max(invdt, invdt_mu, invdt_K, invdt_D)
                ! debug : stockage des differents pas de temps
                !b%dt_vv(l,m) = sb%cfl/invdt_mu
                !b%dt_vT(l,m) = sb%cfl/invdt_K
             end do
          end do
       end if
    end do
    
    call MPI_ALLREDUCE(invdt, invdt_mpi, 1,mpi_real8,mpi_MAX,sb%mpi_comm,statinfo)
    dt_v = sb%cfl/invdt_MPI

  end subroutine dt_visqueux

  subroutine montee_en_cfl(it, montee_cfl, cfl)
    implicit none

    integer              , intent(in)    :: it
    type (STR_montee_cfl), intent(in)    :: montee_cfl
    real*8               , intent(inout) :: cfl
    
    if ( it < montee_cfl%iter_ini  ) then
       cfl = montee_cfl%cfl_ini
    else if ( montee_cfl%iter_ini <= it .and. it < montee_cfl%iter_ini + montee_cfl%iter_trans  ) then
       cfl = (montee_cfl%cfl_end - montee_cfl%cfl_ini)/float(montee_cfl%iter_trans)*float(it-montee_cfl%iter_ini) + montee_cfl%cfl_ini
    else
       cfl = montee_cfl%cfl_croisiere
    endif
  
  end subroutine montee_en_cfl

end module delta_t
