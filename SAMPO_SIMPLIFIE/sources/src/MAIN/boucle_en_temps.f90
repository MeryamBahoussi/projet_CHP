module m_boucle_en_temps
  use delta_t
  use m_calcul_fluide
  use variables_globales
  use m_MPI
  use m_gestion_IO
  implicit none
  
contains
  
  subroutine boucle_en_temps(tache)

    type (STR_TACHE), intent(inout) :: tache
    type (STR_SUPER_BLOC), pointer :: sb
    integer :: iter, i, isb
    logical :: arret=.false.
    integer :: ierr
    real*8  :: tnormal
    
    call init_petsc(tache)
    
    call creation_dossiers_output

    iter=iter_ini
    if (iter==0) then
       if (film) call sorties(tache, iter_ini, it_out)
    end if
    call write_debug(tache, iter, temps)

    if(tache%type_freqout == "temps") then
      i=1
      do while(tache%t_imposes(i)<=temps)
        i = i+1
      end do
      tache%ind_timp = i
    end if

    call debut_watchtime(wt_boucle_temps)

    do while (arret .eqv. .false.) ! boucle en temps
       
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       temps_n = temps  ! temps a l instant n
       iter = iter + 1

       call deltat(tache, iter, arret, tache%dt)

       do i = 1, tache%nb_membre 
          isb = tache%membre(i)
          sb => acces_super_bloc(isb)%pointeur
                  
          if (sb%mpi_comm == MPI_COMM_NULL) cycle

          call calcul_fluide(sb, iter)
          
          ! dans le cas des sorties en temps, il faut vérifier qu'on a pas fait de bétises sur le t_objectif
          tnormal = temps + tache%dt
          if (tache%type_freqout=="temps" .and. tnormal == tache%t_imposes(tache%ind_timp-1)) then
             tache%ind_timp = tache%ind_timp - 1
          end if
       end do
       
       call sorties_boucle_en_temps(tache, iter)
       if (iter >= itermax) arret=.true.   
       if (arret) exit

       ! todo : a remettre ici
       call sorties_boucle_en_temps(tache, iter)
    end do
    call fin_watchtime(wt_boucle_temps)

    call write_debug(tache, -1, temps)
    
#ifdef PETSC_ACTIF
    call PetscFinalize(ierr)
#endif
    
    if (numproc==0) write(6,1000) iter, temps, tache%dt
1000 format(' ITER=',i8,'   ',' Temps=',e13.5,'   ',' DT=',e13.5)
    
  end subroutine boucle_en_temps
  
end module m_boucle_en_temps
