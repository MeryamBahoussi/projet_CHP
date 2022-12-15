module m_calcul_fluide
  use m_euler_direct
  use m_boundary_condition
  use delta_t
  use outils_fluide
  implicit none

contains

  subroutine calcul_fluide(sb, iter)
    type (STR_SUPER_BLOC), pointer :: sb 
    integer           , intent(in) :: iter

    call debut_watchTime(wT_fluide_calcul)

    ! On garde l'etat a l'instant n en memoire
    call debut_watchTime(wT_fluide_Un)
    call save_Un(sb)
    call fin_watchTime(wT_fluide_Un)
    
    if (.not. sb%sol_stat%actif) then
       call calcul_fluide_direct(sb)
    else ! cas stationnaire
       call calcul_fluide_stationnaire(sb, iter)
    end if

    call fin_watchTime(wT_fluide_calcul)

  end subroutine calcul_fluide
  
  subroutine calcul_fluide_stationnaire(sb, iter, affichage)
    use m_residu

    type (STR_SUPER_BLOC), pointer :: sb
    integer, intent(in) :: iter
    logical, intent(in), optional :: affichage ! est ce qu on affiche les residus ?

    integer :: i, k, neq, taux_impression, kmax
    logical :: arret, affichage_residu
    real*8  :: tol
    real*8, dimension(:) , allocatable :: norme_Un, residu
    type (STR_montee_cfl), pointer     :: montee_cfl
    
    
    affichage_residu = .true.
    if (present(affichage)) affichage_residu = affichage
    
    k=0
    neq = acces_bloc(sb%membre(1))%pointeur%neq
    arret=.false.
    
    allocate(norme_Un(neq), residu(neq))
    residu(:) = 0.d0
    call calc_norme_sb(sb, neq, norme_Un)
    if (affichage_residu) call check_residu(k,norme_Un(:),residu(:))
    
    if (numproc==0 .and. affichage_residu) then
       print*
       print*, "------------ debut recherche etat stationnaire ------------"
    end if
    
!!! Modif preconvergence fluide
    if (iter==1 .and. sb%preconvergence%actif) then ! iter_ini+1 ??
       montee_cfl => sb%preconvergence%montee_cfl
       kmax = sb%preconvergence%kmax
       tol = sb%preconvergence%tol
    else
       montee_cfl => sb%montee_cfl
       kmax = sb%sol_stat%kmax
       tol = sb%sol_stat%tol
       sb%cfl = sb%montee_cfl%cfl_croisiere
    end if
    
    ! taux d'impression des residus
    taux_impression = max(1,kmax/12)
    
    do while (arret .eqv. .false.)
       
       ! On garde l'etat a l'instant n en memoire
       call save_Un(sb)
       
       ! calcul du pseudo pas de temps pour le fluide
       if (montee_cfl%actif) call montee_en_cfl(k, montee_cfl, sb%cfl)
       
       call dt_fluide(sb, k, sb%dt)
       call calcul_fluide_direct(sb)
              
       k = k+1
       
       ! Calcul de la convergence si on fait un calcul stationnaire
       call calc_residu_sb(sb, neq, norme_Un, residu)
       if (affichage_residu) then
          call check_residu(k,norme_Un(:),residu(:))
          ! debug : affichage des differents residu pendant le calcul
          if (numproc==0 .and. mod(k, taux_impression)==0) then
             write(6,1000) k, sb%cfl, sb%dt, maxval(residu)
          end if
       end if
       
       ! conditions d'arret
       if (maxval(residu) < tol) arret=.true.
       if (k>kmax) arret=.true.
       
    end do
    
    if (numproc==0) then
       print*
       if (maxval(residu)<tol) then 
          write(*,'("    ITER = ",i5, " solution stationnaire trouvee pour le fluide en",i5," iterations")') iter, k
       else
          write(*,'("    ITER = ",i5, " solution toujours pas stationnarisee apres",i5," iterations")') iter, k
       end if
       if (affichage_residu) then
          print*, "    residu final"
          write(*,'("    ")', advance='no')
          do i = 1, neq
             write(*,'(i12)', advance='no') i
          end do
          write(*,*)
          write(*,'("    ")', advance='no')
          do i = 1, neq
             write(*,'(e12.5)', advance='no') residu(i)
          end do
          print*, "------------- fin recherche etat stationnaire -------------"
          print*
       end if
    end if
1000 format("     k=",i5,"  cfl=",e13.5,"  dt=",e13.5,"  max(residu)=",e13.5)
    
    deallocate(norme_Un, residu)
  end subroutine calcul_fluide_stationnaire

end module m_calcul_fluide
