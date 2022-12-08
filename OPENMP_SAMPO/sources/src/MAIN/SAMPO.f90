program main
  use m_MPI
  use m_input
  use m_allocation
  use m_initialisation
  use m_gestion_IO
  use m_boucle_en_temps
  use OMP_LIB
  implicit none

!#include "finclude/petsc.h"

  type (STR_TACHE) ::  tache
  integer :: STATINFO

  call MPI_initialisation
  call init_watchtime

  if ( numproc == 0 )  then
     print*,""
     print*,"DATE DE DEBUT D'EXECUTION"
     call system("date")
     print*,""
     print '("Calcul sur ", i4, " processeur(s) avec ", i2, " thread(s) OpenMP")', nbprocs, OMP_GET_MAX_THREADS()
     print*,""
  end if

  call read_input(tache)
  call allocation(tache)

  ! initialisation des differents communicateurs mpi
  call init_communicateurs(tache)

  call initialisation_calcul(tache)

  call boucle_en_temps(tache)

  call MPI_BARRIER(MPI_COMM_WORLD, STATINFO)
  !call reprise(tache, -1)
  call sorties(tache, -1, it_out)


  if ( numproc == 0 )  then
     print*,""
     print*,"DATE DE FIN D'EXECUTION"
     call system("date")

     print*,""
     print*, "---------------------------------------------------------------"
     print*, "- Analyse du temps de calcul des différentes portions du code -"
     print*, "---------------------------------------------------------------"
  end if
  call affichage_watchTime_total(wT_boucle_temps , "+ boucle en temps     ")
  call affichage_watchTime_total(wT_fluide_calcul, "+ calcul fluide       ", wT_boucle_temps%ttot)
  call affichage_watchTime_total(wT_fluide_Un    , "    - save Un           ", wT_fluide_calcul%ttot)
  call affichage_watchTime_total(wT_fluide_sm    , "    - assemblage sm     ", wT_fluide_calcul%ttot)
  call affichage_watchTime_total(wT_fluide_maj   , "    - maj fluide        ", wT_fluide_calcul%ttot)
  call affichage_watchTime_total(wT_fluide_decode, "    - decodage fluide   ", wT_fluide_calcul%ttot)
  call affichage_watchTime_total(wT_deltat       , "+ delta t             ", wT_boucle_temps%ttot)
  call affichage_watchTime_total(wT_dt_decode    , "    - dt decode         ", wT_deltat%ttot)
  call affichage_watchTime_total(wT_dt_CL        , "    - dt CL             ", wT_deltat%ttot)
  call affichage_watchTime_total(wT_dt_pente     , "    - dt pente          ", wT_deltat%ttot)
  call affichage_watchTime_total(wT_dt_Euler     , "    - dt Euler          ", wT_deltat%ttot)
  call affichage_watchTime_total(wT_decodage     , "+ decodage            ", wT_boucle_temps%ttot)
  call affichage_watchTime_total(wT_decodage_U     , "    - decode U         ", wT_boucle_temps%ttot)
  call affichage_watchTime_total(wT_sorties      , "+ sorties             ", wT_boucle_temps%ttot)
  call affichage_watchTime_total(wT_commMPI      , "+ communication MPI   ", wT_boucle_temps%ttot)
  call affichage_watchTime_procs(wT_commMPI      , "+ communication MPI   ")



  call deallocation(tache)
  call MPI_fin

end program main
