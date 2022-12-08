module m_temps
  use m_MPI

  implicit none

  real*8 :: temps=0.d0
  real*8 :: temps_n=0.d0

  ! pour calculer les temps cpu de certaines portions du code
  type STR_WATCHTIME
     real*8  :: tdeb, tfin, ttot
     integer :: nb
  end type STR_WATCHTIME

  type (STR_WATCHTIME) :: wT_boucle_temps
  type (STR_WATCHTIME) :: wT_fluide_calcul
  type (STR_WATCHTIME) :: wT_fluide_Un
  type (STR_WATCHTIME) :: wT_fluide_sm
  type (STR_WATCHTIME) :: wT_fluide_maj
  type (STR_WATCHTIME) :: wT_fluide_decode
  type (STR_WATCHTIME) :: wT_deltat
  type (STR_WATCHTIME) :: wT_commMPI
  type (STR_WATCHTIME) :: wT_sorties
  type (STR_WATCHTIME) :: wT_dt_decode
  type (STR_WATCHTIME) :: wT_dt_CL
  type (STR_WATCHTIME) :: wT_dt_pente
  type (STR_WATCHTIME) :: wT_dt_Euler
  type (STR_WATCHTIME) :: wT_decodage
  type (STR_WATCHTIME) :: wT_decodage_U

contains

  subroutine init_watchtime()
    wT_boucle_temps%ttot =0.d0; wT_boucle_temps%nb =0
    wT_fluide_calcul%ttot=0.d0; wT_fluide_calcul%nb=0
    wT_fluide_Un%ttot    =0.d0; wT_fluide_Un%nb    =0
    wT_fluide_maj%ttot   =0.d0; wT_fluide_maj%nb   =0
    wT_fluide_decode%ttot=0.d0; wT_fluide_decode%nb=0
    wT_deltat%ttot       =0.d0; wT_deltat%nb       =0
    wT_commMPI%ttot      =0.d0; wT_commMPI%nb      =0
    wT_sorties%ttot      =0.d0; wT_sorties%nb      =0

    wT_dt_decode%ttot    =0.d0; wT_dt_decode%nb    =0
    wT_dt_CL%ttot        =0.d0; wT_dt_CL%nb        =0
    wT_dt_pente%ttot     =0.d0; wT_dt_pente%nb     =0
    wT_dt_Euler%ttot     =0.d0; wT_dt_Euler%nb     =0

    wT_decodage%ttot     =0.d0; wT_decodage%nb     =0
    wT_decodage_U%ttot     =0.d0; wT_decodage_U%nb     =0

  end subroutine init_watchtime

  subroutine debut_watchtime(wT)
    type (STR_WATCHTIME), intent(out) :: wT

    !call cpu_time(wT%tdeb)
    wT%tdeb = MPI_WTIME()
  end subroutine debut_watchtime

  subroutine fin_watchtime(wT)
    type (STR_WATCHTIME), intent(inout) :: wT

    !call cpu_time(wT%tfin)
    wT%tfin = MPI_WTIME()
    wT%ttot = wT%ttot + wT%tfin-wT%tdeb
    wT%nb = wT%nb+1
  end subroutine fin_watchtime

  subroutine affichage_watchTime_procs(wT, msg)
    type (STR_WATCHTIME), intent(in) :: wT
    character(len=*)    , intent(in) :: msg

    integer :: STATINFO
    call MPI_BARRIER(MPI_COMM_WORLD, STATINFO)
    print '(" ", A, " (proc", i2,")", f18.10, " s")', msg, numproc, wT%ttot
    call MPI_BARRIER(MPI_COMM_WORLD, STATINFO)
  end subroutine affichage_watchTime_procs

  subroutine affichage_watchTime_total(wT, msg, total)
    type (STR_WATCHTIME), intent(in) :: wT
    character(len=*)    , intent(in) :: msg
    real*8              , intent(in), optional :: total

    ! variables locales
    integer :: STATINFO
    real*8  :: ttot_mpi, total_mpi
    character(len=len(msg)) :: void

    void=""

    if (present(total)) then
       call MPI_ALLREDUCE(total,total_mpi,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,STATINFO)
    end if
    call MPI_ALLREDUCE(wT%ttot,ttot_mpi,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,STATINFO)
    if ( numproc == 0 )  then
       if (present(total)) then
          print '(" ", A, " (total) ", f18.10, " s", f18.10, "%")', msg, ttot_mpi, ttot_mpi*100./total_mpi
       else
          print '(" ", A, " (total) ", f18.10, " s")', msg, ttot_mpi
       end if
       !print '(A,"  (nb appel)", i7)', void, wT%nb
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, STATINFO)
  end subroutine affichage_watchTime_total
end module m_temps
