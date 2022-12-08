module m_gestion_IO
  use m_boundary_condition
  use m_IO
  use m_residu
  use variables_sorties

contains

  subroutine creation_dossiers_output
    logical :: test_dir

    ! on teste si le dossier de sortie existe deja
    inquire(FILE=trim(path_output)//'visit', exist=test_dir)
    if (.not. test_dir) then
       call system("mkdir "//trim(path_output)//'visit')
    end if

    ! on teste si le dossier gnuplot existe deja
    inquire(FILE=trim(path_output)//'gnuplot', exist=test_dir)
    if (.not. test_dir) then
       call system("mkdir "//trim(path_output)//'gnuplot')
    end if

    ! on teste si le dossier de reprise existe deja
    inquire(FILE=trim(path_output)//'restart', exist=test_dir)
    if (.not. test_dir) then
       call system("mkdir "//trim(path_output)//'restart')
    end if

  end subroutine creation_dossiers_output

  subroutine sorties_boucle_en_temps(tache, iter)
    type (STR_TACHE), intent(in) :: tache
    integer         , intent(in) :: iter

    if (film) then
       if (mod(iter, it_restart) == 0) then
          call reprise(tache, iter)
       end if
       select case (tache%type_freqout)
       case("iter")
          if (mod(iter, it_out) == 0 .or. iter==1) then
             call sorties(tache, iter, it_out)
          end if
       case("temps")
          if (tache%ind_timp /= 1 ) then ! pour ne pas depasser le tableau
             if (temps==tache%t_imposes(tache%ind_timp-1)) then
                call sorties(tache, tache%ind_timp-1, 1)
             end if
          end if
       end select
    end if
    
!!!! Dans tous les cas, on fait les sorties paroi et de debuggage      
    if (mod(iter, it_debug) == 0) then 
       call write_debug(tache, iter, temps)
    end if

  end subroutine sorties_boucle_en_temps
  
  
  subroutine sorties(tache, iter, iter_out)

    type (STR_TACHE), intent(in) ::  tache
    integer, intent(in) :: iter, iter_out

    integer                        :: i, j, isb, ib
    character(len=6)               :: iter_c
    type (STR_BLOC)      , pointer :: b
    type (STR_SUPER_BLOC), pointer :: sb

    call debut_watchtime(wt_sorties)
    
    do i = 1, tache%nb_membre 
       isb = tache%membre(i)
       sb => acces_super_bloc(isb)%pointeur 
        
       if (sb%mpi_comm == MPI_COMM_NULL) cycle

       if (numproc == 0) then
          iter_c = write_iter_c(iter, iter_out)
          print *, 'Ecriture des fichiers ', trim(path_output)//'visit/output_xxxx_'//trim(iter_c)//'.vtk'
       end if

       call decode(sb, debug_MF)
       if (debug_MF) call condition_limites_Euler(sb)

       do j = 1, sb%nb_membre 
          ib = sb%membre(j)
          b => acces_bloc(ib)%pointeur
          if (bloc_est_local (ib)) then
             if (debug_MF.and.b%visqueux) call condition_limites_visqueux(b,sb) ! a faire ?
             call sorties_visit_nblocs(b, iter, iter_out)
          end if
       end do
       
       ! sortie 1D communes au fluide et au solide
       if (output_1d .and. iter==-1) then
          if (acces_bloc(1)%pointeur%mf==1) then
             call sorties_gnuplot(sb, DIR_X)
          else if (acces_bloc(1)%pointeur%lf==1) then 
             call sorties_gnuplot(sb, DIR_Y)
          else
             if (acces_bloc(1)%pointeur%lf>acces_bloc(1)%pointeur%mf) then
                print*, "On trace les sorties gnuplot le long de l'axe x"
                call sorties_gnuplot(sb, DIR_X)
             else
                print*, "On trace les sorties gnuplot le long de l'axe y"
                call sorties_gnuplot(sb, DIR_Y)
             end if
          end if
       end if
    end do

        call debut_watchtime(wt_sorties)
    
    do i = 1, tache%nb_membre 
       isb = tache%membre(i)
       sb => acces_super_bloc(isb)%pointeur 
        
       if (sb%mpi_comm == MPI_COMM_NULL) cycle

       if (numproc == 0) then
          iter_c = write_iter_c(iter, iter_out)
          print *, 'Ecriture des fichiers ', trim(path_output)//'visit/output_xxxx_'//trim(iter_c)//'.vtk'
       end if

       call decode(sb, debug_MF)
       if (debug_MF) call condition_limites_Euler(sb)

       do j = 1, sb%nb_membre 
          ib = sb%membre(j)
          b => acces_bloc(ib)%pointeur
          if (bloc_est_local (ib)) then
             if (debug_MF.and.b%visqueux) call condition_limites_visqueux(b,sb) ! a faire ?
             call sorties_visit_nblocs(b, iter, iter_out)
          end if
       end do
       
       ! sortie 1D communes au fluide et au solide
       if (output_1d .and. iter==-1) then
          if (acces_bloc(1)%pointeur%mf==1) then
             call sorties_gnuplot(sb, DIR_X)
          else if (acces_bloc(1)%pointeur%lf==1) then 
             call sorties_gnuplot(sb, DIR_Y)
          else
             if (acces_bloc(1)%pointeur%lf>acces_bloc(1)%pointeur%mf) then
                print*, "On trace les sorties gnuplot le long de l'axe x"
                call sorties_gnuplot(sb, DIR_X)
             else
                print*, "On trace les sorties gnuplot le long de l'axe y"
                call sorties_gnuplot(sb, DIR_Y)
             end if
          end if
       end if
    end do

    call fin_watchtime(wt_sorties)
  end subroutine sorties

  subroutine reprise(tache, iter)

    type (STR_TACHE), intent(in) ::  tache
    integer, intent(in) :: iter

    type (STR_SUPER_BLOC), pointer :: sb
    type (STR_BLOC), pointer :: b
    integer :: i, j, isb, ib
 
    do i = 1, tache%nb_membre 
       isb = tache%membre(i)
       sb => acces_super_bloc(isb)%pointeur 

       do j = 1, sb%nb_membre 
          ib = sb%membre(j)
          if (bloc_est_local (ib)) then
             b => acces_bloc(ib)%pointeur
             call protection_reprise(b, iter, sb%nb_membre_input)
          end if
       end do
    end do

  end subroutine reprise

end module m_gestion_IO
