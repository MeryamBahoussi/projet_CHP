module m_MPI
  use m_struct
  use parametres_globaux, only : FLUIDE, SOLIDE
  implicit none

  integer  :: nbprocs         ! Nombre de processeurs total
  integer  :: numproc         ! Rang du processeur courant dans le communicateur MPI_COMM_WORLD
  integer  :: numproc_fluide  ! Rang du processeur courant dans le communicateur fluide
  integer  :: numproc_solide  ! Rang du processeur courant dans le communicateur solide
  integer  :: mpi_comm_fluide ! Communicateur MPI du domaine fluide
  integer  :: mpi_comm_solide ! Communicateur MPI du domaine solide
 
  public MPI_initialisation, arret_code, MPI_fin
  
  include 'mpif.h'
contains

  ! Initialisation de la librairie MPI.
  subroutine MPI_initialisation()
    implicit none

    integer :: STATINFO

    call MPI_INIT(STATINFO)
    if (STATINFO /= MPI_SUCCESS) then
       print*, "erreur initialisation MPI"
       call arret_code
    end if

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nbprocs, STATINFO)
    if (STATINFO /= MPI_SUCCESS) then
       print*, "erreur initialisation MPI"
       call arret_code
    end if

    call MPI_COMM_RANK(MPI_COMM_WORLD, numproc, STATINFO)
    if (STATINFO /= MPI_SUCCESS) then
       print*, "erreur initialisation MPI"
       call arret_code
    end if

  end subroutine MPI_initialisation

  ! Fin de l'utilisation de la librairie MPI.
  subroutine MPI_fin()
    implicit none
    integer  :: STATINFO
    
    call MPI_FINALIZE(STATINFO)
  end subroutine MPI_fin

  ! Calcul de la charge de chaque proc
  subroutine MPI_charge(k1,kn,Np,me,i1,in)
    implicit none
    integer,intent(in)::k1,kn,Np,me
    integer,intent(in out)::i1,in
    integer::q,r,n

    n = kn - k1 + 1

    q=n/Np
    r=n-Np*q

    if(me<r) then 
       i1=me*(q+1)+1
       in=(me+1)*(q+1)
    else
       i1=1+r+me*q
       in=i1+q-1
    end if

    i1 = i1 + k1 - 1
    in = in + k1 - 1 
  end subroutine MPI_charge

  subroutine MPI_charge_bloc(b)
    use m_struct
    implicit none
    
    type (STR_BLOC), intent(inout) ::  b
    integer :: me 
    
    allocate( b%l_charge(1:2,0:nbprocs-1) )
    allocate( b%m_charge(1:2,0:nbprocs-1) )

    do me = 0, nbprocs-1
       call MPI_charge(b%ld_tot, b%lf_tot, nbprocs, me, b%l_charge(1,me), b%l_charge(2,me))
       call MPI_charge(b%md_tot, b%mf_tot, nbprocs, me, b%m_charge(1,me), b%m_charge(2,me))
    end do
    
  end subroutine MPI_charge_bloc


  ! renvoie le numero de processeur qui possede ce bloc
  integer function bloc_numproc (ib)
    use m_struct
    implicit none
    integer, intent (in) :: ib

    bloc_numproc = acces_bloc(ib)%pointeur%numproc

  end function bloc_numproc

  ! est-ce que le proc. possede le bloc
  logical function bloc_est_local(ib)
    use m_struct
    implicit none
    integer, intent(in)  :: ib

    bloc_est_local = bloc_numproc (ib) .EQ. numproc

  end function bloc_est_local
  
  subroutine arret_code()
    integer  :: ierr

    call MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
    stop
  end subroutine arret_code

  !****************************************************************
  !*                                                              *
  !*  initialisation de l ensemble des sous communicateurs        *
  !*                                                              *
  !****************************************************************
  subroutine init_communicateurs(tache)
    type (STR_TACHE), intent(in) :: tache 

    call init_comm_fluide(tache)
    call init_comm_solide(tache)
    call init_comm_super_bloc(tache)
    
  end subroutine init_communicateurs

  !****************************************************************
  !*                                                              *
  !*  initialisation des communicateurs intra super bloc          *
  !*                                                              *
  !****************************************************************
  subroutine init_comm_super_bloc(tache)
    type (STR_TACHE)  , intent(in) :: tache 

    ! variables locales
    type (STR_SUPER_BLOC), pointer :: s_b
    type (STR_BLOC)      , pointer :: b_i
    integer                        :: s, u, ibu, ib, ibi
    integer                        :: couleur, ierr, cle
    integer, dimension(nbprocs)    :: list_numproc_sb

    ! Parcours des super blocs
    do s = 1, tache%nb_membre
       if (.not.(associated(acces_super_bloc(tache%membre(s))%pointeur)))  cycle       
       s_b => acces_super_bloc(tache%membre(s))%pointeur

       couleur = MPI_UNDEFINED
       cle = -1

       ! Parcours des blocs informatiques
       do ib = 1, s_b%nb_membre
          ! Numéro du b_i
          ibi = s_b%membre(ib)
          b_i => acces_bloc(ibi)%pointeur

          if (bloc_est_local (b_i%numero)) then
             if (cle == -1) cle = ibi
             couleur = tache%membre(s)
          end if
       end do
       
       ! Creation des sous communicateurs intra super bloc
       call MPI_COMM_SPLIT(MPI_COMM_WORLD, couleur, cle, s_b%mpi_comm, ierr)

       if (s_b%mpi_comm /= MPI_COMM_NULL) then
          ! Calcul du rang du processeur courant dans le sous communicateur
          call MPI_COMM_RANK(s_b%mpi_comm, s_b%numproc, ierr)
          ! Calcul de la taille du sous communicateur
          call MPI_COMM_SIZE(s_b%mpi_comm, s_b%nbprocs, ierr)
       else
          s_b%numproc = -1
          s_b%nbprocs = 0
       end if
       
       list_numproc_sb(:) = -1

       ! Communications du rang entre tous les processeurs
       call MPI_ALLGATHER(s_b%numproc, 1, MPI_INTEGER4, list_numproc_sb, 1, MPI_INTEGER4, MPI_COMM_WORLD, ierr)

       ! Stockage du rang dans le bloc informatique         
       do ib = 1, s_b%nb_membre
          ibi = s_b%membre(ib)
          b_i => acces_bloc(ibi)%pointeur
          b_i%numproc_sb = list_numproc_sb(b_i%numproc + 1)
       end do
    end do
       
  end subroutine init_comm_super_bloc

  !****************************************************************
  !*                                                              *
  !*  initialisation du communicateur du domaine fluide           *
  !*                                                              *
  !****************************************************************
  subroutine init_comm_fluide(tache)
    type (STR_TACHE)  , intent(in) :: tache 

    ! variables locales
    type (STR_SUPER_BLOC), pointer :: s_b
    type (STR_BLOC)      , pointer :: b_i
    integer                        :: s, u, ibu, ib, ibi
    integer                        :: couleur, ierr
    integer, dimension(nbprocs)    :: list_numproc_fluide

    couleur = MPI_UNDEFINED

    ! Parcours des super blocs
    do s = 1, tache%nb_membre
       if (.not.(associated(acces_super_bloc(tache%membre(s))%pointeur)))  cycle       
       s_b => acces_super_bloc(tache%membre(s))%pointeur

       ! Parcours des blocs informatiques
       do ib = 1, s_b%nb_membre
          ! Numéro du b_i
          ibi = s_b%membre(ib)
          b_i => acces_bloc(ibi)%pointeur
          
          if (bloc_est_local (b_i%numero) .and. (b_i%nature == FLUIDE)) then
             couleur = 1
          end if
       end do
    end do

    ! Creation des sous communicateurs intra super bloc
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, couleur, numproc, mpi_comm_fluide, ierr)

    if (mpi_comm_fluide /= MPI_COMM_NULL) then
      ! call MPI_COMM_SIZE(mpi_comm_fluide, s, ierr)
       call MPI_COMM_RANK(mpi_comm_fluide, numproc_fluide, ierr)
    else
       numproc_fluide = -1
    end if

    ! Communications du rang entre tous les processeurs
    call MPI_ALLGATHER(numproc_fluide, 1, MPI_INTEGER4, list_numproc_fluide, 1, MPI_INTEGER4, MPI_COMM_WORLD, ierr)

    ! Stockage du rang dans le bloc informatique
    do s = 1, tache%nb_membre
       if (.not.(associated(acces_super_bloc(tache%membre(s))%pointeur)))  cycle       
       s_b => acces_super_bloc(tache%membre(s))%pointeur

       do ib = 1, s_b%nb_membre
          ibi = s_b%membre(ib)
          b_i => acces_bloc(ibi)%pointeur
          b_i%numproc_fluide = list_numproc_fluide(b_i%numproc + 1)
       end do
    end do

  end subroutine init_comm_fluide

  !****************************************************************
  !*                                                              *
  !*  initialisation du communicateur du domaine solide           *
  !*                                                              *
  !****************************************************************
  subroutine init_comm_solide(tache)
    type (STR_TACHE)  , intent(in) :: tache 

    ! variables locales
    type (STR_SUPER_BLOC), pointer :: s_b
    type (STR_BLOC)      , pointer :: b_i
    integer                        :: s, u, ibu, ib, ibi
    integer                        :: couleur, ierr
    integer, dimension(nbprocs)    :: list_numproc_solide

    couleur = MPI_UNDEFINED

    ! Parcours des super blocs
    do s = 1, tache%nb_membre
       if (.not.(associated(acces_super_bloc(tache%membre(s))%pointeur)))  cycle       
       s_b => acces_super_bloc(tache%membre(s))%pointeur

       ! Parcours des blocs informatiques
       do ib = 1, s_b%nb_membre
          ! Numéro du b_i
          ibi = s_b%membre(ib)
          b_i => acces_bloc(ibi)%pointeur
          
          if (bloc_est_local (b_i%numero) .and. (b_i%nature == SOLIDE)) then
             couleur = 1
          end if
       end do
    end do

    ! Creation des sous communicateurs intra super bloc
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, couleur, numproc, mpi_comm_solide, ierr)

    if (mpi_comm_solide /= MPI_COMM_NULL) then
       call MPI_COMM_RANK(mpi_comm_solide, numproc_solide, ierr)
    else
       numproc_solide = -1
    end if

    ! Communications du rang entre tous les processeurs
    call MPI_ALLGATHER(numproc_solide, 1, MPI_INTEGER4, list_numproc_solide, 1, MPI_INTEGER4, MPI_COMM_WORLD, ierr)

    ! Stockage du rang dans le bloc informatique
    do s = 1, tache%nb_membre
       if (.not.(associated(acces_super_bloc(tache%membre(s))%pointeur)))  cycle       
       s_b => acces_super_bloc(tache%membre(s))%pointeur

       do ib = 1, s_b%nb_membre
          ibi = s_b%membre(ib)
          b_i => acces_bloc(ibi)%pointeur
          b_i%numproc_solide = list_numproc_solide(b_i%numproc + 1)
       end do
    end do

  end subroutine init_comm_solide

end module m_MPI
