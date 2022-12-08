module m_euler_direct
  use second_euler_direct
  use second_visqueux
  use m_boundary_condition
  use premier_euler_direct
  use premier_visqueux_direct
  use m_solveurs
  use m_residu
  use outils_pas_acoustique
  use parametres_globaux, only : nMF

  implicit none 

contains

  subroutine calcul_fluide_direct(sb)
    
    type (STR_SUPER_BLOC), pointer ::  sb

    integer :: i, ib
    type (STR_BLOC), pointer :: b

    integer :: l, m
    
    do i = 1, sb%nb_membre 
       ib = sb%membre(i)
       b => acces_bloc(ib)%pointeur

       if (bloc_est_local (ib)) then
          if (.not. b%implicite) then
             ! assemblage du second membre
             call debut_watchTime(wT_fluide_sm)
             call assemblage_sm(sb, b, sb%dt, b%U, b%U, b%sm)
             call fin_watchTime(wT_fluide_sm)

             ! mise à jour de la solution fluide
             call debut_watchTime(wT_fluide_maj)
             do m = b%md, b%mf
                do l = b%ld, b%lf
                   b%U(:,l,m) = b%U(:,l,m) + b%sm(:,l,m)
                end do
             end do
             call fin_watchTime(wT_fluide_maj)

             ! decodage de la solution
             call debut_watchTime(wT_fluide_decode)
             call decode(sb, decode_MF=.true.)
             call fin_watchTime(wT_fluide_decode)
          else 
             call Navier_stokes_NEWTON(b, sb)
          end if
       end if
    end do
  end subroutine calcul_fluide_direct
  
  subroutine Navier_stokes_NEWTON(b, sb)
    type (STR_BLOC)      , pointer :: b
    type (STR_SUPER_BLOC), pointer :: sb
    
    integer                    :: k
    integer                    :: ld, lf, md, mf
    real*8                     :: dt, norme_Xn(b%neq), residu(b%neq)
    real*8           , pointer :: Xkp1(:,:,:), Xk(:,:,:), Xn(:,:,:), sm(:,:,:)
    type(STR_MATRICE), pointer :: mat

    ld=b%ld; lf=b%lf; md=b%md; mf=b%mf
    dt = sb%dt
    mat => b%mat
    sm => b%sm
    
    Xkp1 => b%U
    Xk   => b%X_old; Xk(:,:,:) = Xkp1(:,:,:)
    Xn   => b%U_n
    call norme(sb%mpi_comm,b%neq,lf,mf,Xn(:,ld:lf,md:mf),norme_Xn)

    delta_forme =.true.  

    ! si theta schema : il faut stocker le second membre calcule avec X_0 (=U_n)
    
!!!!!!!!!!!!!!!!   DEBUT DE LA BOUCLE DE NEWTON 
    k = 1
    do
       !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤ Mise a jour des pentes et des conditions limites
       call condition_limites_Euler(sb)
       call pentes_solveur_Riemann(sb)
       
       !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤ ASSEMBLAGE DE LA MATRICE ET DU SECOND MEMBRE
       call assemblage_sm(sb, b, dt, Xkp1, Xn, sm)
       call assemblage_matrice_direct(sb, b, dt, mat)
       
       !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤ RESOLUTION DU SYSTEME LINEAIRE
       select case(b%solveur_implicite)
       case( MKL )
          call RL_MKL_sur_proc0(b, Xkp1(:,ld:lf,md:mf), mat, sm)
       case( PETSc )
          call RL_petsc_global(b, Xkp1(:,ld-2:lf+2,md-2:mf+2), sb, mat, sm, (/2,1/), sb%inv_up%petsc, 1)
       end select
       
       !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤ CALCUL DE LA CONVERGENCE DU NEWTON
       call calc_residu(sb%mpi_comm, b%neq,lf,mf,Xkp1(:,ld:lf,md:mf), &
            Xk(:,ld:lf,md:mf),norme_Xn,residu) 
       
       !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤ PREPARATION DE L'ITERATION SUIVANTE
       ! decodage de Xkp1 (les MF sont decodees dans condition_limites_Euler)
       call decode(sb, decode_MF=.false.) 
       k = k+1
       Xk(:,:,:) = Xkp1(:,:,:)
       
       if (k >= sb%inv_up%kmax .or. maxval(residu) < sb%inv_up%tol) exit
    end do
!!!!!!!!!!!!!!!!   FIN DE LA BOUCLE DE NEWTON 
    
    if (sb%inv_up%kmax>1) then
       if (k>=sb%inv_up%kmax .and. maxval(residu)>sb%inv_up%tol .and. numproc==b%proc_ini) then
          print '("Newton fluide non converge" :, e14.6,e14.6,e14.6,e14.6,e14.6,e14.6)', residu(irho1z1), residu(irho2z2), residu(irhou), residu(irhov), residu(irhoE), residu(iz)
          print*
       end if
    end if

  end subroutine Navier_stokes_NEWTON

  subroutine assemblage_sm(sb, b, dt, Xkp1, Xn, sm)
    use m_verification
    type (STR_SUPER_BLOC)   , pointer    :: sb
    type (STR_BLOC)         , pointer    :: b
    real*8                  , intent(in) :: dt
    real*8, dimension(:,:,:), pointer    :: Xkp1, Xn, sm

    integer :: l, m
    real*8, dimension(:,:), pointer :: An, Anp1

    sm(:,:,:) = 0.d0
    call second_membre_euler_direct(b, sb, sm)
    if (b%visqueux) then
       call condition_limites_visqueux(b, sb)
       call second_membre_visqueux(b, sm)
    end if
    
    if (b%maillage%move .and. .not. sb%sol_stat%actif) then 
       ! prise en compte du maillage mobile
       An => b%aire
       Anp1 => b%aire_new
       print*, 'maillage mobile + direct non gere'
       call arret_code
       ! call terme_maillage_mobile(b, sm)
    else
       An => b%aire
       Anp1 => b%aire
    end if
    
    do m = b%md, b%mf
       do l = b%ld, b%lf
          sm(:,l,m) = -(Xkp1(:,l,m) - (An(l,m)/Anp1(l,m))*Xn(:,l,m)) - dt/Anp1(l,m)*sm(:,l,m)
       end do
    end do

  end subroutine assemblage_sm
  
  subroutine assemblage_matrice_direct(sb, b, dt, mat)
    type (STR_SUPER_BLOC), pointer    :: sb
    type (STR_BLOC)      , pointer    :: b
    real*8               , intent(in) :: dt    
    type (STR_MATRICE)   , pointer    :: mat
    
    integer :: i, j, l, m
    real*8, dimension(:,:), pointer :: An, Anp1
    
    mat%AE = 0.d0; mat%E = 0.d0; mat%CE = 0.d0
    mat%A  = 0.d0; mat%B = 0.d0; mat%C  = 0.d0
    mat%AD = 0.d0; mat%D = 0.d0; mat%CD = 0.d0
    
!!! on commence par la matrice visqueuse car on doit la multiplier par dV/dU
!!! mat = mat_Eu + mat_Vs*dV/dU
!!! Ca permet de ne pas stoker les deux matrices
    if (b%visqueux) then
       call assemblage_matrice_visqueux(sb, b, mat)
    end if
    call condition_limites_Euler(sb) ! on recalcule les CL, car la matrice depend du vecteur des variables conservatives U
    call assemblage_matrice_Euler(sb, b, mat)
    
    if (b%maillage%move .and. .not. sb%sol_stat%actif) then 
       ! prise en compte du maillage mobile
       An => b%aire
       Anp1 => b%aire_new
    else
       An => b%aire
       Anp1 => b%aire
    end if

!!! mat = Id + dt/Aire * mat
    do m=b%md, b%mf
       do l = b%ld, b%lf
          do j = 1, b%neq
             do i = 1, b%neq
                mat%A(i,j,l,m) = dt/Anp1(l,m) * mat%A(i,j,l,m)
                mat%B(i,j,l,m) = dt/Anp1(l,m) * mat%B(i,j,l,m)
                mat%C(i,j,l,m) = dt/Anp1(l,m) * mat%C(i,j,l,m)
                mat%D(i,j,l,m) = dt/Anp1(l,m) * mat%D(i,j,l,m)
                mat%E(i,j,l,m) = dt/Anp1(l,m) * mat%E(i,j,l,m)
                
                mat%AD(i,j,l,m) = dt/Anp1(l,m) * mat%AD(i,j,l,m)
                mat%CD(i,j,l,m) = dt/Anp1(l,m) * mat%CD(i,j,l,m)
                mat%AE(i,j,l,m) = dt/Anp1(l,m) * mat%AE(i,j,l,m)
                mat%CE(i,j,l,m) = dt/Anp1(l,m) * mat%CE(i,j,l,m)
             end do
             
             mat%b(j,j,l,m) = 1.d0 + mat%b(j,j,l,m) 
          end do
       end do
    end do
    
    
    if (b%solveur_implicite == MKL)  then
       call MKL_assemblage(b, mat)
    end if
  end subroutine assemblage_matrice_direct

end module m_euler_direct
