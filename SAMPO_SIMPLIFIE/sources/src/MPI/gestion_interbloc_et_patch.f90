module m_gestion_interbloc
  use m_struct
  use m_MPI
  use parametres_globaux
  implicit none

contains

!!! gestion du decoupage lie a MPI
  subroutine decoupage

    integer :: i, statinfo, choix_decoupe
    integer :: NbCoupes(size(acces_bloc,1),2)

    type (STR_BLOC)      , pointer :: b

    ! verification si on peut lancer le cas
    if (implicite_avec_petsc() .and. nbprocs < size(acces_bloc, 1)) then
       print*, "pas assez de processeurs pour lancer ce cas avec PETSc" 
       call arret_code
    end if

    ! lecture du nombre de decoupages suivant l et m
    call lecture_nb_decoupes(NbCoupes)

    ! determine le choix repartition des procs 
    call determine_type_decoupe(NbCoupes, choix_decoupe)

    call decoupe_blocs_2d(NbCoupes, choix_decoupe)

    ! on retrouve sb%membre et sb%nb_membre
    do i = 1, size(acces_super_bloc,1) 
       call retrouve_membre_sb_2d(acces_super_bloc(i)%pointeur, NbCoupes)
    end do

    call retrouve_num_bloc_adjacent_2d

    call MPI_BARRIER(MPI_COMM_WORLD, STATINFO) ! utile ?

    ! on change les conditions aux limites pour les interblocs
    do i=1, size(acces_bloc, 1)
       b=> acces_bloc(i)%pointeur
       call patch_interbloc_mpi_2d(b, Nbcoupes(b%num_input,:))
    end do

    ! call creation_interbloc
  end subroutine decoupage

  ! Si on utilise petsc, on doit avoir un bloc par proc
  function implicite_avec_petsc()

    logical :: implicite_avec_petsc

    integer :: i

    implicite_avec_petsc = .false.

    do i = 1, size(acces_super_bloc,1)
       if ( acces_super_bloc(i)%pointeur%implicite .and. &
            acces_super_bloc(i)%pointeur%solveur_implicite == PETSc) then
          implicite_avec_petsc =.true.
          return
       end if
    end do

  end function implicite_avec_petsc

  
  subroutine lecture_nb_decoupes(NbCoupes)
    use variables_globales, only : path_input

    integer, intent(out) :: NbCoupes(size(acces_bloc,1),2)

    integer :: i, ib, nbblocs
    integer :: iunit = 46
    type (STR_BLOC), pointer :: b
    
    open(iunit,file=trim(path_input)//"decoupage.dat",status="unknown")
    read(iunit,*) ! # numero du bloc | nb blocs suivant l | nb blocs suivant m
    nbblocs = 0
    do i =1, size(acces_bloc, 1)
       read(iunit,*) ib, NbCoupes(ib,1), NbCoupes(ib,2)
       
       nbblocs = nbblocs + NbCoupes(i,1)*NbCoupes(i,2)
       b => acces_bloc(ib)%pointeur
       b%Nbcoupes(i,:) = NbCoupes(i,:)
    end do

  end subroutine lecture_nb_decoupes

  ! deux choix pour le type de decoupe
  ! - choix_decoupe=1 : 1 proc <-> 1 bloc 
  ! - choix_decoupe=2 : 1 proc <-> 1 bloc (possible) dans chaque superbloc 
  subroutine determine_type_decoupe(NbCoupes, choix_decoupe)
    integer, intent(in)  :: NbCoupes(size(acces_bloc,1),2)
    integer, intent(out) :: choix_decoupe
    ! variables locales
    type (STR_SUPER_BLOC), pointer :: s_b
    integer                        :: isb, ib, ibi
    integer                        :: nb_procs_decoupe1, nb_procs_decoupe2, nb_proc_sb
    
    nb_procs_decoupe1 = 0
    do ib = 1, size(acces_bloc,1)
       nb_procs_decoupe1 = nb_procs_decoupe1 + product(NbCoupes(ib,:))
    end do

    nb_procs_decoupe2 = 0
    do isb = 1, size(acces_super_bloc,1)
       s_b => acces_super_bloc(isb)%pointeur
       nb_proc_sb = 0
       do ib = 1, s_b%nb_membre
          ibi = s_b%membre(ib)
          nb_proc_sb = nb_proc_sb + product(NbCoupes(ibi,:))
       end do
       nb_procs_decoupe2 = max(nb_procs_decoupe2, nb_proc_sb)
    end do
    
    if (nbprocs == nb_procs_decoupe1) then
       choix_decoupe = 1
    else if (nbprocs == nb_procs_decoupe2) then
       choix_decoupe = 2
    else
       print*, "le nombre de proc n'est pas compatible avec le fichier entrees/decoupage.dat"
       print*, "il faut utiliser ", nb_procs_decoupe1, "procs ou ", nb_procs_decoupe2, "procs"
       call arret_code
    end if
    
  end subroutine determine_type_decoupe

  !***********************************************!
  !         Routine de decoupe d'un bloc          !
  !***********************************************!
  subroutine decoupe_blocs_2d(NbCoupes, choix_decoupe)

    integer, intent(in) :: NbCoupes(size(acces_bloc,1),2)
    integer, intent(in) :: choix_decoupe

    ! variables locales
    integer                            :: isb, ibi
    integer                            :: i, ib, ib_new, proc_ini, me, nb_blocs, sb_ini
    integer                            :: mex, mey
    type (STR_BLOC)      , pointer     :: b_pere
    type (STR_BLOC)      , pointer     :: b_fils
    type (STR_SUPER_BLOC), pointer     :: s_b
    type (STR_ACCES_BLOC), allocatable :: acces_bloc_new(:)

    ! le nombre de blocs apres decoupe est independant du choix de la decoupe
    nb_blocs = 0
    do ib = 1, size(acces_bloc,1)
       nb_blocs = nb_blocs + product(NbCoupes(ib,:))
    end do
    
    allocate(acces_bloc_new(nb_blocs))
    do i = 1, nb_blocs
       nullify(acces_bloc_new(i)%pointeur)
    end do

    proc_ini = 0
    sb_ini = 0
    ! boucle sur les super-blocs
    do isb = 1, size(acces_super_bloc,1)
       s_b => acces_super_bloc(isb)%pointeur

       do ibi = 1, s_b%nb_membre
          ib = s_b%membre(ibi)

          ! bloc pere = bloc a decouper
          b_pere => acces_bloc(ib)%pointeur

          b_pere%proc_ini = proc_ini
          b_pere%proc_fin = proc_ini+product(NbCoupes(ib,:))-1

          deallocate(b_pere%l_charge)
          allocate( b_pere%l_charge(1:2,b_pere%proc_ini:b_pere%proc_fin) )
          deallocate(b_pere%m_charge)
          allocate( b_pere%m_charge(1:2,b_pere%proc_ini:b_pere%proc_fin) )
          b_pere%l_charge(:,:) = 0
          b_pere%m_charge(:,:) = 0

          do mey = 0, NbCoupes(ib,2)-1
             do mex = 0, NbCoupes(ib,1)-1

                allocate(b_fils)
                call copie_bloc(b_pere, b_fils)
                call MPI_charge(b_pere%ld_tot, b_pere%lf_tot, NbCoupes(ib,1), &
                     mex, b_fils%ld, b_fils%lf)
                call MPI_charge(b_pere%md_tot, b_pere%mf_tot, NbCoupes(ib,2), &
                     mey, b_fils%md, b_fils%mf)

                me = mex + mey*NbCoupes(ib,1) + b_pere%proc_ini
                b_fils%l_charge(1,me) = b_fils%ld
                b_fils%l_charge(2,me) = b_fils%lf
                b_fils%m_charge(1,me) = b_fils%md
                b_fils%m_charge(2,me) = b_fils%mf

                b_fils%num_input = ib
                b_fils%numero = sb_ini+me+1
                b_fils%numproc = me

                ib_new = sb_ini+me+1
                if (.not. associated(acces_bloc_new(ib_new)%pointeur)) then
                   allocate( acces_bloc_new(ib_new)%pointeur )
                else
                   print*, " input.f90 : deux blocs ont le meme numero apres la decoupe", ib_new
                   call arret_code
                end if
                acces_bloc_new(ib_new)%pointeur = b_fils
                deallocate(b_fils)
             end do
          end do

          proc_ini = proc_ini+product(NbCoupes(ib,:))
       end do
       if (choix_decoupe==2) then
          sb_ini = proc_ini
          proc_ini = 0
       end if
    end do
   
    ! La numerotation de chaque bloc commence a 1
    do ib = 1, nb_blocs
       b_fils => acces_bloc_new(ib)%pointeur
       b_fils%lf = b_fils%lf-b_fils%ld+1
       b_fils%ld = 1
       b_fils%mf = b_fils%mf-b_fils%md+1
       b_fils%md = 1
    end do

    ! creation du nouvel acces_bloc
    deallocate(acces_bloc)
    allocate(acces_bloc(1:nb_blocs)) 
    do i = 1, nb_blocs
       acces_bloc(i) = acces_bloc_new(i)
    end do

  end subroutine decoupe_blocs_2d
  
  subroutine copie_bloc(b_pere, b_fils)

    type (STR_BLOC), pointer :: b_fils, b_pere
    integer :: i

    b_fils = b_pere    ! pour toutes les variables communes
    do i = 1, 4
       allocate(b_fils%face(i)%patch(b_pere%face(i)%nb_patch))
       b_fils%face(i)%patch(:) = b_pere%face(i)%patch(:)
    end do

  end subroutine copie_bloc

  subroutine retrouve_membre_sb_2d(sb, NbCoupes)

    type (STR_SUPER_BLOC), pointer    :: sb
    integer              , intent(in) :: NbCoupes(:,:)

    integer                            :: i, j, k, ib, ib_old, ib_new 
    integer                            :: nb_membre_new
    integer, dimension(:), allocatable :: membre_new

    ! on retrouve le nombre de bloc dans un super-blocs
    nb_membre_new=0
    do i = 1, sb%nb_membre 
       ib = sb%membre(i)
       nb_membre_new = nb_membre_new + product(NbCoupes(ib,:))
    end do

    ! on retrouve le nouveau numero des blocs appartenant au super-bloc
    allocate(membre_new(1:nb_membre_new))
    k = 1 ! compteur sur les nouveaux membres
    ! boucle sur les nouveaux blocs
    do i = 1, size(acces_bloc,1)
       ib_new = acces_bloc(i)%pointeur%num_input
       ! boucle sur les membres du super bloc
       do j = 1, sb%nb_membre
          ib_old = sb%membre(j)
          if (ib_new == ib_old) then
             membre_new(k) = i
             k = k+1
          end if
       end do
    end do

    if (k/=nb_membre_new+1) then
       print*, "retrouve_membre_sb : erreur dans le decoupage"
       call arret_code
    end if

    sb%nb_membre = nb_membre_new
    deallocate(sb%membre)
    allocate(sb%membre(1:sb%nb_membre))
    sb%membre = membre_new
    deallocate(membre_new)

  end subroutine retrouve_membre_sb_2d

  subroutine retrouve_num_bloc_adjacent_2d()
    
    type (STR_BLOC), pointer :: b, b_adj
    type (STR_PATCH), pointer :: patch
    integer :: ib, iface, ip, ib_adj, proc_adj, if_adj

    integer :: ladj_old, madj_old, me

    do ib = 1, size(acces_bloc, 1)
       b => acces_bloc(ib)%pointeur
       !if (bloc_est_local(ib)) then
       me = b%numproc
          do iface = 1, 4
             do ip=1, b%face(iface)%nb_patch
                patch => b%face(iface)%patch(ip)
                ! test pour voir si on a un interbloc ou une iterface
                if ( patch%bc%condition=="interbloc" .or. &
                     index(patch%bc%condition, "itf") /=0 ) then
                   ! on recupere le numero de la face du bloc adjacent
                   if_adj = patch%if_adj

                   ! boucle sur les possibles blocs adjacents
                   b_b_adj : do ib_adj = 1, size(acces_bloc, 1)
                      b_adj => acces_bloc(ib_adj)%pointeur
                      ! test pour savoir si un bloc issu du bloc (pere) adjacent
                      if ( patch%ib_adj == b_adj%num_input) then
                         
                         if (if_adj==1) then
                            ladj_old = b_adj%ld_tot
                            madj_old = patch%ideb_adj+b%m_charge(1,me)-b%md_tot
                         else if (if_adj==2) then
                            ladj_old = b_adj%lf_tot
                            madj_old = patch%ideb_adj+b%m_charge(1,me)-b%md_tot
                         else if (if_adj==3) then
                            ladj_old = patch%ideb_adj+b%l_charge(1,me)-b%ld_tot
                            madj_old = b_adj%md_tot
                         else if (if_adj==4) then
                            ladj_old = patch%ideb_adj+b%l_charge(1,me)-b%ld_tot
                            madj_old = b_adj%mf_tot  
                         end if
                         
                         proc_adj = b_adj%numproc
                         if ( ladj_old >= b_adj%l_charge(1,proc_adj) .and. &
                              ladj_old <= b_adj%l_charge(2,proc_adj) .and. &
                              madj_old >= b_adj%m_charge(1,proc_adj) .and. &
                              madj_old <= b_adj%m_charge(2,proc_adj) ) then
                            patch%ib_adj = b_adj%numero
                            exit b_b_adj
                         end if
 
                      end if
                   end do b_b_adj
                end if
             end do
          end do
       !end if
    end do

  end subroutine retrouve_num_bloc_adjacent_2d

  subroutine patch_interbloc_mpi_2d(b, Nbcoupes)

    type (STR_BLOC)      , pointer    :: b
    integer, dimension(2), intent(in) :: Nbcoupes

    type (STR_BLOC) , pointer :: b_adj
    type (STR_FACE) , pointer :: face
    type (STR_PATCH), pointer :: patch
    integer                   :: ip
    integer         , pointer :: ideb, ifin, ideb_adj

    !face 1  west : (ld,m)
    if (b%l_charge(1,b%numproc) /= b%ld_tot) then
       ! ce n'est pas une frontiere du bloc pere
       face => b%face(1)
       face%nb_patch = 1
       deallocate(face%patch)
       allocate(face%patch(1))
       patch => face%patch(1)
       patch%bc%position = 1
       patch%bc%condition = "interbloc"
       patch%ideb = b%md
       patch%ifin = b%mf
       patch%ib_adj = b%numero-1
       patch%if_adj = 2
       patch%sens = 1
       patch%ideb_adj = b%md
    else
       do ip=1, b%face(1)%nb_patch
          patch => b%face(1)%patch(ip)
          ideb => patch%ideb
          ideb = max(ideb, b%md)
          ifin => patch%ifin
          ifin = min(ifin, b%mf)

          if (patch%bc%condition=="interbloc") then
             b_adj => acces_bloc(patch%ib_adj)%pointeur
             ideb_adj => patch%ideb_adj
             ideb_adj = max(ideb_adj, b_adj%ld)
          end if
       end do
    end if

    !face 2  east : (lf,m)
    if (b%l_charge(2,b%numproc) /= b%lf_tot) then
       face => b%face(2)
       face%nb_patch = 1
       deallocate(face%patch)
       allocate(face%patch(1))
       patch => face%patch(1)
       patch%bc%position = 2
       patch%bc%condition = "interbloc"
       patch%ideb = b%md
       patch%ifin = b%mf
       patch%ib_adj = b%numero+1
       patch%if_adj = 1
       patch%sens = 1
       patch%ideb_adj = b%md
    else
       do ip=1, b%face(2)%nb_patch
          patch => b%face(2)%patch(ip)
          ideb => patch%ideb
          ideb = max(ideb, b%md)
          ifin => patch%ifin
          ifin = min(ifin, b%mf)

          if (patch%bc%condition=="interbloc") then
             b_adj => acces_bloc(patch%ib_adj)%pointeur
             ideb_adj => patch%ideb_adj
             ideb_adj = max(ideb_adj, b_adj%ld)
          end if
       end do
    end if

    !face 3 south : (l,md)
    if (b%m_charge(1,b%numproc) /= b%md_tot) then
       face => b%face(3)
       face%nb_patch = 1
       deallocate(face%patch)
       allocate(face%patch(1))
       patch => face%patch(1)
       patch%bc%position = 3
       patch%bc%condition = "interbloc"
       patch%ideb = b%ld
       patch%ifin = b%lf
       patch%ib_adj = b%numero-NbCoupes(1)
       patch%if_adj = 4
       patch%sens = 1
       patch%ideb_adj = b%ld
    else
       do ip=1, b%face(3)%nb_patch
          patch => b%face(3)%patch(ip)
          ideb => patch%ideb
          ideb = max(ideb, b%ld)
          ifin => patch%ifin
          ifin = min(ifin, b%lf)

          if (patch%bc%condition=="interbloc") then
             b_adj => acces_bloc(patch%ib_adj)%pointeur
             ideb_adj => patch%ideb_adj
             ideb_adj = max(ideb_adj, b_adj%ld)
          end if
       end do
    end if

    !face 4 north : (l,mf)
    if (b%m_charge(2,b%numproc) /= b%mf_tot) then
       face => b%face(4)
       face%nb_patch = 1
       deallocate(face%patch)
       allocate(face%patch(1))
       patch => face%patch(1)
       patch%bc%position = 4
       patch%bc%condition = "interbloc"
       patch%ideb = b%ld
       patch%ifin = b%lf
       patch%ib_adj = b%numero+NbCoupes(1)
       patch%if_adj = 3
       patch%sens = 1
       patch%ideb_adj = b%ld
    else
       do ip=1, b%face(4)%nb_patch
          patch => b%face(4)%patch(ip)
          ideb => patch%ideb
          ideb = max(ideb, b%ld)
          ifin => patch%ifin
          ifin = min(ifin, b%lf)

          if (patch%bc%condition=="interbloc") then
             b_adj => acces_bloc(patch%ib_adj)%pointeur
             ideb_adj => patch%ideb_adj
             ideb_adj = max(ideb_adj, b_adj%ld)
          end if
       end do
    end if

  end subroutine patch_interbloc_mpi_2d

!!$!***********************************************!
!!$!       Routines de decoupage d'un bloc         !
!!$!                (decoupage 1D)                 !
!!$!***********************************************! 
!!$  subroutine decoupe_blocs(lf_total, nbprocs_b)
!!$    
!!$    integer, intent(in) :: lf_total
!!$    integer, intent(in) :: nbprocs_b(size(acces_bloc,1))
!!$
!!$    integer                  :: i, ib, proc_ini, me, nb_blocs, sb_ini
!!$    type (STR_BLOC), pointer :: b
!!$    type (STR_ACCES_BLOC), dimension(:), allocatable :: acces_bloc_new
!!$    
!!$    ! au final on doit avoir nbprocs blocs
!!$    if (choix_decoupe==1) then
!!$       nb_blocs = nbprocs
!!$    else
!!$       nb_blocs = nbprocs*size(acces_super_bloc,1)
!!$    end if
!!$    
!!$    allocate(acces_bloc_new(nb_blocs))
!!$    do i = 1, nb_blocs
!!$       allocate(acces_bloc_new(i)%pointeur)
!!$    end do
!!$    
!!$    proc_ini = 0
!!$    sb_ini = 0
!!$    do ib = 1, size(acces_bloc,1)
!!$       ! bloc a decouper
!!$       b => acces_bloc(ib)%pointeur
!!$       
!!$       b%proc_ini = proc_ini
!!$       b%proc_fin = proc_ini+nbprocs_b(ib)-1
!!$       
!!$       deallocate(b%l_charge)
!!$       allocate( b%l_charge(1:2,b%proc_ini:b%proc_fin) )
!!$       
!!$       do me = b%proc_ini, b%proc_fin
!!$          call MPI_charge(b%ld_tot, b%lf_tot, nbprocs_b(ib), &
!!$               me-proc_ini, b%ld, b%lf)
!!$          
!!$          b%l_charge(1,me) = b%ld
!!$          b%l_charge(2,me) = b%lf
!!$          
!!$          b%num_input = ib
!!$          b%numero = sb_ini+me+1
!!$          b%numproc = me
!!$          
!!$          acces_bloc_new(sb_ini+me+1)%pointeur = b
!!$       end do
!!$       ! La numerotation de chaque bloc commence a 1
!!$       do me = b%proc_ini, b%proc_fin
!!$          b => acces_bloc_new(sb_ini+me+1)%pointeur
!!$          b%lf = b%lf-b%ld+1
!!$          b%ld = 1
!!$       end do
!!$       
!!$       if (choix_decoupe==1) then
!!$          proc_ini = proc_ini+nbprocs_b(ib)
!!$          sb_ini = 0
!!$       else
!!$          proc_ini = mod(proc_ini+nbprocs_b(ib),nbprocs)
!!$          sb_ini = nbprocs ! * ?
!!$       end if
!!$    end do
!!$    
!!$    deallocate(acces_bloc)
!!$    allocate(acces_bloc(1:nb_blocs)) 
!!$    do i = 1, nb_blocs
!!$       acces_bloc(i) = acces_bloc_new(i)
!!$    end do
!!$    
!!$  end subroutine decoupe_blocs

!!$!!! calcul du nombre de processeurs par bloc (decoupage 1D)
!!$!!! On regarde le nombre de mailles suivant l
!!$!!! et on trouve le nombre de proc en fonction de la taille du bloc
!!$  subroutine calcul_nbprocs_par_bloc(lf_total, nbprocs_b)
!!$    
!!$    integer, intent(in)  :: lf_total 
!!$    integer, intent(out) :: nbprocs_b(size(acces_bloc,1)) ! nombre de procs par bloc
!!$    
!!$    integer                     :: i, nbprocs_tot
!!$    type (STR_BLOC), pointer    :: b
!!$        
!!$    do i =1, size(acces_bloc, 1)
!!$       b => acces_bloc(i)%pointeur
!!$       nbprocs_b(i) = int(nbprocs * dble((b%lf_tot-b%ld_tot+1))/lf_total)
!!$      
!!$       if (choix_decoupe==1) then
!!$          if ( nbprocs_b(i) == nbprocs ) then
!!$             nbprocs_b(i) = nbprocs_b(i) - 1   ! on enleve un proc si un bloc les a tous
!!$          else if ( nbprocs_b(i) == 0 ) then
!!$             nbprocs_b(i) = 1                  ! on en met au moins 1
!!$          end if
!!$       end if
!!$    end do
!!$    
!!$    ! verification 
!!$    if (choix_decoupe==1) then
!!$       nbprocs_tot = sum(nbprocs_b)
!!$       if (nbprocs_tot == nbprocs-1) then
!!$          nbprocs_b(1) = nbprocs_b(1)+1
!!$       else if (nbprocs_tot/=nbprocs) then
!!$          if (numproc==0) print*, "erreur dans le nombre de procs total", nbprocs_tot, nbprocs
!!$          call arret_code
!!$       end if
!!$    end if
!!$
!!$  end subroutine calcul_nbprocs_par_bloc
!!$
!!$  ! decoupage 1D
!!$  subroutine retrouve_num_bloc_adjacent
!!$    type (STR_BLOC), pointer :: b, b_adj
!!$    type (STR_PATCH), pointer :: patch
!!$    integer :: ib, iface, ip, ib_adj, proc_adj
!!$    
!!$    integer :: idadj_old
!!$    
!!$    do ib = 1, size(acces_bloc, 1)
!!$       b => acces_bloc(ib)%pointeur
!!$       if (bloc_est_local(ib)) then
!!$          do iface = 1, 4
!!$             do ip=1, b%face(iface)%nb_patch
!!$                patch => b%face(iface)%patch(ip)
!!$                if ( patch%bc%condition=="interbloc" .or. &
!!$                     index(patch%bc%condition, "itf") /=0 ) then
!!$                   
!!$                   do ib_adj = 1, size(acces_bloc, 1)
!!$                      b_adj => acces_bloc(ib_adj)%pointeur
!!$                      if ( patch%ib_adj == b_adj%num_input) then
!!$                         if (iface == 1) then
!!$                            patch%ib_adj =  b_adj%proc_fin+1
!!$                         else if  (iface == 2) then
!!$                            patch%ib_adj =  b_adj%proc_ini+1
!!$                         else if (iface == 3 .or. iface == 4) then
!!$                            ! boucle sur les procs du bloc (pere) adjacent
!!$                            do proc_adj = b_adj%proc_ini, b_adj%proc_fin
!!$                               
!!$                               idadj_old = b%l_charge(1,numproc) - patch%ideb + patch%ideb_adj
!!$                               if ( idadj_old >= b_adj%l_charge(1,proc_adj) .and. &
!!$                                    idadj_old <= b_adj%l_charge(2,proc_adj) ) then
!!$                                  if (choix_decoupe==1) then
!!$                                     patch%ib_adj = proc_adj+1
!!$                                  else
!!$                                     patch%ib_adj = b_adj%numero+proc_adj
!!$                                  end if
!!$                               end if
!!$                            end do
!!$                         end if
!!$                         exit
!!$                      end if
!!$                   end do
!!$                end if
!!$             end do
!!$          end do
!!$       end if
!!$    end do
!!$        
!!$  end subroutine retrouve_num_bloc_adjacent


!!! Numerotation des coins
!!!  4_____3
!!!  |     |
!!!  |  b  |
!!!  |_____|
!!!  1     2
!!!
!!! a chaque coin, il y a deux patchs associes
!!! le patch1 est celui qui varie en l
!!! le patch2 est celui qui varie en m
  subroutine initialisation_coins(b)
    type (STR_BLOC), pointer :: b
    
    type (STR_COIN), pointer :: coin

    ! coin entre la face 1 et 3
    coin => b%coin(1)
    coin%icoin = 1
    coin%indices(1:4) = (/b%ld  , b%ld+(nMF-1), b%md  , b%md+(nMF-1)/)
    coin%indices(5:8) = (/b%ld-nMF, b%ld-1, b%md-nMF, b%md-1/)
    coin%patch1 => b%face(3)%patch(1)
    coin%patch2 => b%face(1)%patch(1)
    coin%n1 => b%n_dir_m(1:2,b%ld,b%md-1)
    coin%n2 => b%n_dir_l(1:2,b%ld-1,b%md)
    call ib_opp_coin(b, coin, coin%ib_opp)
    if (b%nature == FLUIDE) call condition_coin(coin, coin%condition)
    

    ! coin entre la face 3 et 2
    coin => b%coin(2)
    coin%icoin = 2
    coin%indices(1:4) = (/b%lf-(nMF-1), b%lf  , b%md  , b%md+(nMF-1)/)
    coin%indices(5:8) = (/b%lf+1, b%lf+nMF, b%md-nMF, b%md-1/)
    coin%patch1 => b%face(3)%patch(b%face(3)%nb_patch)
    coin%patch2 => b%face(2)%patch(1)
    coin%n1 => b%n_dir_m(1:2,b%lf,b%md-1)
    coin%n2 => b%n_dir_l(1:2,b%lf,b%md)
    call ib_opp_coin(b, coin, coin%ib_opp)
    if (b%nature == FLUIDE) call condition_coin(coin, coin%condition)

    ! coin entre la face 2 et 4
    coin => b%coin(3)
    coin%icoin = 3
    coin%indices(1:4) = (/b%lf-(nMF-1), b%lf  , b%mf-(nMF-1), b%mf  /)
    coin%indices(5:8) = (/b%lf+1, b%lf+nMF, b%mf+1, b%mf+nMF/)
    coin%patch1 => b%face(4)%patch(b%face(4)%nb_patch)
    coin%patch2 => b%face(2)%patch(b%face(2)%nb_patch)
    coin%n1 => b%n_dir_m(1:2,b%lf,b%mf)
    coin%n2 => b%n_dir_l(1:2,b%lf,b%mf)
    call ib_opp_coin(b, coin, coin%ib_opp)
    if (b%nature == FLUIDE) call condition_coin(coin, coin%condition)
    
    ! coin entre la face 4 et 1
    coin => b%coin(4)
    coin%icoin = 4
    coin%indices(1:4) = (/b%ld  , b%ld+(nMF-1), b%mf-(nMF-1), b%mf  /)
    coin%indices(5:8) = (/b%ld-nMF, b%ld-1, b%mf+1, b%mf+nMF/)
    coin%patch1 => b%face(4)%patch(1)
    coin%patch2 => b%face(1)%patch(b%face(1)%nb_patch)
    coin%n1 => b%n_dir_m(1:2,b%ld,b%mf)
    coin%n2 => b%n_dir_l(1:2,b%ld-1,b%mf)
    call ib_opp_coin(b, coin, coin%ib_opp)
    if (b%nature == FLUIDE) call condition_coin(coin, coin%condition)

  end subroutine initialisation_coins

!!!!!!!!!!!
!!! Trouve le numero du bloc a l'oppose du coin
!!!!!!!!!!!
!!! pour chaque patch contenant le coin
!!! * on regarde le numero du bloc adjacent
!!!   si c'est le bord (=-1) ib_opp = -1
!!!   si non
!!! * on trouve le numero du patch_adj
!!! * contient il un coin ? 
!!!   si non => ib_opp = ib_adj
!!!   si oui
!!!     on trouve le patch qui contient aussi le coin dans ce bloc_adj
!!!     ib_opp = patch_opp%ib_adj
!!!
!!! il faut faire ca a partir des deux patchs qui contiennent le coin
!!! si ib_opp(patch1) /= ib_opp(patch2) (/=-1) c'est qu on est dans un cas non gere
!!! Exemple de cas non gere
!!!    b1  /
!!!   ____/  b2
!!!       \
!!!    b3  \
!!!
  subroutine ib_opp_coin(b, coin, ib_opp)
    type (STR_BLOC), pointer     :: b 
    type (STR_COIN), pointer     :: coin
    integer        , intent(out) :: ib_opp

    integer :: ib_opp1, ib_opp2
    
    call ib_opp_coin_patch(coin, 1, coin%patch1, ib_opp1)
    call ib_opp_coin_patch(coin, 2, coin%patch2, ib_opp2)

    if (ib_opp1 /= -1 .and. ib_opp2 /= -1) then
       if (ib_opp1 /= ib_opp2) then
          print*, "cas non gere pour les coins"
          call arret_code
       else
          ib_opp = ib_opp1
       end if
    else if (ib_opp1 == -1 .and. ib_opp2 /= -1) then
       ib_opp = ib_opp2
    else if (ib_opp1 /= -1 .and. ib_opp2 == -1) then
       ib_opp = ib_opp1
    else if (ib_opp1 == -1 .and. ib_opp2 == -1) then
       ib_opp = -1
    end if
    
    if (ib_opp /= -1) then
       if (acces_bloc(ib_opp)%pointeur%nature /= b%nature) then
          ib_opp = -1 ! coin entre bloc fluide et solide non gere ! a modifier
       end if
    end if

  end subroutine ib_opp_coin

  ! pos = 1 si le coin est au debut de la face
  ! pos = 2 si le coin est a la fin de la face
  subroutine ib_opp_coin_patch(coin, num_patch, patch, ib_opp)
    type (STR_COIN) , pointer     :: coin
    type (STR_PATCH), pointer     :: patch
    integer         , intent(in)  :: num_patch   !(1 ou 2) 
    integer         , intent(out) :: ib_opp
    
    type (STR_BLOC), pointer  :: b_adj
    type (STR_FACE), pointer  :: f_adj
    type (STR_PATCH), pointer :: patch_opp

    integer :: ip_adj
    integer :: pos, deb, fin
    deb = 1; fin = 2
    
    if (num_patch==1) then
       if (coin%icoin == 1 .or. coin%icoin == 4) then
          pos = deb
       else
          pos = fin
       end if
    else
       if (coin%icoin == 1 .or. coin%icoin == 2) then
          pos = deb
       else
          pos = fin
       end if
    end if

    if (patch%ib_adj /= -1) then
       b_adj => acces_bloc(patch%ib_adj)%pointeur
       f_adj => b_adj%face(patch%if_adj)
       
       call ipatch_adj(patch, ip_adj)

       if ((ip_adj==1 .and. pos==deb) .or. (ip_adj==f_adj%nb_patch .and. pos==fin)) then
          ! on trouve le patch qui contient aussi le coin dans ce bloc_adj
          call infos_patch_commun_coin(b_adj, patch%if_adj, pos, patch_opp)
          ib_opp = patch_opp%ib_adj
       else  ! le patch n'est pas sur un bord de la face
          ib_opp = patch%ib_adj
          ! dans ce cas, il faut faire quelque chose de particulier
          ! dans la communication des coins
          ! on va recevoir des infos mais ne pas en envoyer
       end if
    else
       ib_opp = -1
    end if

  end subroutine ib_opp_coin_patch
  
!!! Retrouve le numero du patch adjacent
  subroutine ipatch_adj(patch, ip_adj)
    type (STR_PATCH), pointer     :: patch
    integer         , intent(out) :: ip_adj

    type (STR_BLOC), pointer :: b_adj
    type (STR_FACE), pointer :: f_adj
    integer :: i

    b_adj => acces_bloc(patch%ib_adj)%pointeur
    f_adj => b_adj%face(patch%if_adj)

    ip_adj = -1

    do i = 1, f_adj%nb_patch
       if ( patch%ideb_adj >= f_adj%patch(i)%ideb .and. &
            patch%ideb_adj <= f_adj%patch(i)%ifin ) then

          ip_adj = i
       end if
    end do

    if (ip_adj == -1) then
       print*, "erreur dans ipatch_adj"
       call arret_code
    end if

  end subroutine ipatch_adj
  
  ! on trouve les infos du patch qui contient aussi le coin dans le (meme) bloc
  ! pos = 1, on regarde le debut de la face (i=ideb)
  ! pos = 2, on regarde la fin   de la face  (i=ifin)
  subroutine infos_patch_commun_coin(b, iface, pos, patch_commun)
    type (STR_BLOC) , pointer     :: b
    integer         , intent(in)  :: iface, pos
    type (STR_PATCH), pointer     :: patch_commun

    integer :: iface_commun, pos_commun

    select case(iface)
    case(1)
       if (pos==1) then
          iface_commun = 3; pos_commun = 1
       else
          iface_commun = 4; pos_commun = 1
       end if
    case(2)
       if (pos==1) then
          iface_commun = 3; pos_commun = 2
       else
          iface_commun = 4; pos_commun = 2
       end if
    case(3)
       if (pos==1) then
          iface_commun = 1; pos_commun = 1
       else
          iface_commun = 2; pos_commun = 1
       end if
    case(4)
       if (pos==1) then
          iface_commun = 1; pos_commun = 2
       else
          iface_commun = 2; pos_commun = 2
       end if
    end select
    
     if (pos_commun==1) patch_commun=>b%face(iface_commun)%patch(1)
  !  if (pos_commun==2) patch_commun=>b%face(iface_commun)%patch(b%face(iface)%nb_patch)
     if (pos_commun==2) patch_commun=>b%face(iface_commun)%patch(b%face(iface_commun)%nb_patch)     
   end subroutine infos_patch_commun_coin


   subroutine condition_coin(coin, condition)
     use parametres_fluide
     
     type (STR_COIN), pointer     :: coin
     integer        , intent(out) :: condition
     
     logical           :: erreur
     character(len=50) :: bc1, bc2
     bc1 = coin%patch1%bc%condition
     bc2 = coin%patch2%bc%condition
    
     erreur=.false.
     
     select case (trim(bc1))
     case("periodique")
        select case(trim(bc2))
        case("interbloc")
           condition = CL_PERIODIQUEvsINTERBLOC
        case("periodique")
           condition = CL_PERIODIQUEvsPERIODIQUE
        case("paroi_global")
          condition = CL_PERIODIQUEvsPAROI
        case("paroi_global_adiabatique")
          condition = CL_PERIODIQUEvsPAROI_ADIAB
        case default
           erreur=.true.
        end select
     case("flux_nul")
        select case (trim(bc2))
        case("flux_nul")
           condition = CL_FLUX_NULvsFLUX_NUL
        case("symetrie")
           condition = CL_FLUX_NULvsSYMETRIE
        case("interbloc")
           condition = CL_FLUX_NULvsINTERBLOC
        case("interproc")
           condition = CL_INTERPROCvs
        case("entree_supersonique")
           condition = CL_FLUX_NULvsENTREE_SUPER
        case("entree_subsonique")
           condition = CL_FLUX_NULvsENTREE_SUB
        case("entree_subsonique_t")  
           erreur=.true.
        case ("sortie_subsonique")
           condition = CL_FLUX_NULvsSORTIE_SUB
        case("paroi_global")
           condition = CL_FLUX_NULvsPAROI
        case("paroi_global_adiabatique")
           condition = CL_FLUX_NULvsPAROI_ADIAB
        case default
           erreur=.true.
        end select
     case ("symetrie")
        select case (trim(bc2))
        case("flux_nul")
           condition = CL_SYMETRIEvsFLUX_NUL
        case("symetrie")
           condition = CL_SYMETRIEvsSYMETRIE
        case("interbloc")
           condition = CL_SYMETRIEvsINTERBLOC
        case("interproc")
           condition = CL_INTERPROCvs
        case("entree_supersonique")
           condition = CL_SYMETRIEvsENTREE_SUPER
        case("entree_subsonique")
           condition = CL_SYMETRIEvsENTREE_SUB
        case("entree_subsonique_t")  
           erreur=.true.
        case ("sortie_subsonique")
           condition = CL_SYMETRIEvsSORTIE_SUB
        case("paroi_global")
           condition = CL_SYMETRIEvsPAROI
        case("paroi_global_adiabatique")
           condition = CL_SYMETRIEvsPAROI_ADIAB
        case default
           erreur=.true.
        end select
     case("interbloc")
        select case (trim(bc2))
        case("flux_nul")
           condition = CL_INTERBLOCvsFLUX_NUL
        case("symetrie")
           condition = CL_INTERBLOCvsSYMETRIE
        case("interbloc")
           condition = CL_INTERBLOCvsINTERBLOC
        case("interproc")
           condition = CL_INTERPROCvs
        case("entree_supersonique")
           condition = CL_INTERBLOCvsENTREE_SUPER
           ! case("entree_subsonique")
           ! case("entree_subsonique_t")  
           ! case ("sortie_subsonique")
        case("paroi_global")
           condition = CL_INTERBLOCvsPAROI
        case("paroi_global_adiabatique")
           condition = CL_INTERBLOCvsPAROI_ADIAB
        case("itf_fusion", "itf_fusion_sans_solide")
           condition = CL_INTERBLOCvsINJECTION
        case("periodique")
           condition = CL_INTERBLOCvsPERIODIQUE
        case default
           erreur=.true.
        end select
     case("interproc")
        condition = CL_INTERPROCvs
     case("entree_supersonique")
        select case (trim(bc2))
        case("flux_nul")
           condition = CL_ENTREE_SUPERvsFLUX_NUL
        case("symetrie")
           condition = CL_ENTREE_SUPERvsSYMETRIE
        case("interbloc")
           condition = CL_ENTREE_SUPERvsINTERBLOC
        case("interproc")
           condition = CL_INTERPROCvs
        case("entree_supersonique", "entree_subsonique","entree_subsonique_t")  
           erreur=.true.
        case ("sortie_subsonique")
           erreur=.true.
        case("paroi_global")
           condition = CL_ENTREE_SUPERvsPAROI
        case("paroi_global_adiabatique")
           condition = CL_ENTREE_SUPERvsPAROI_ADIAB
        case default
           erreur=.true.
        end select
     case("entree_subsonique")
        select case (trim(bc2))
        case("flux_nul")
           condition = CL_ENTREE_SUBvsFLUX_NUL
        case("symetrie")
           condition = CL_ENTREE_SUBvsSYMETRIE
        !case("interbloc")
        !   condition = CL_INTERBLOCvs
        case("interproc")
           condition = CL_INTERPROCvs
        case("entree_supersonique", "entree_subsonique","entree_subsonique_t")  
           erreur=.true.
        case ("sortie_subsonique")
           erreur=.true.
        case("paroi_global")
           condition = CL_ENTREE_SUBvsPAROI
        case("paroi_global_adiabatique")
           condition = CL_ENTREE_SUBvsPAROI_ADIAB
        case default
           erreur=.true.
        end select
     case("entree_subsonique_t")
        erreur=.true. ! doit etre possible de le faire
     case("sortie_subsonique")
        select case (trim(bc2))
        case("flux_nul")
           condition = CL_SORTIE_SUBvsFLUX_NUL
        case("symetrie")
           condition = CL_SORTIE_SUBvsSYMETRIE
        !case("interbloc")
        !   condition = CL_INTERBLOCvs
        case("interproc")
           condition = CL_INTERPROCvs
        case("entree_supersonique", "entree_subsonique","entree_subsonique_t")  
           erreur=.true.
        case ("sortie_subsonique")
           erreur=.true.
        case("paroi_global")
           condition = CL_SORTIE_SUBvsPAROI
        case("paroi_global_adiabatique")
           condition = CL_SORTIE_SUBvsPAROI_ADIAB
        case default
           erreur=.true.
        end select
     case("paroi_global")
        select case (trim(bc2))
        case("flux_nul")
           condition = CL_PAROIvsFLUX_NUL
        case("symetrie")
           condition = CL_PAROIvsSYMETRIE
        case("interbloc")
           condition = CL_PAROIvsINTERBLOC
        case("interproc")
           condition = CL_INTERPROCvs
        case("entree_supersonique")
           condition = CL_PAROIvsENTREE_SUPER
        case("entree_subsonique")
           condition = CL_PAROIvsENTREE_SUB
        case("entree_subsonique_t")  
           erreur=.true.
        case ("sortie_subsonique")
           condition = CL_PAROIvsSORTIE_SUB
        case("paroi_global")
           condition = CL_PAROIvsPAROI
        case("paroi_global_adiabatique")
           condition = CL_PAROIvsPAROI_ADIAB
        case("periodique")
           condition = CL_PAROIvsPERIODIQUE
        case default
           erreur=.true.
        end select
     case("paroi_global_adiabatique")
        select case (trim(bc2))
        case("flux_nul")
           condition = CL_PAROI_ADIABvsFLUX_NUL
        case("symetrie")
           condition = CL_PAROI_ADIABvsSYMETRIE
        case("interbloc")
           condition = CL_PAROI_ADIABvsINTERBLOC
        case("interproc")
           condition = CL_INTERPROCvs
        case("entree_supersonique")
           condition = CL_PAROI_ADIABvsENTREE_SUPER
        case("entree_subsonique")
           condition = CL_PAROI_ADIABvsENTREE_SUB
        case("entree_subsonique_t")  
           erreur=.true.
        case ("sortie_subsonique")
           condition = CL_PAROI_ADIABvsSORTIE_SUB
        case("paroi_global")
           condition = CL_PAROI_ADIABvsPAROI
        case("paroi_global_adiabatique")
           condition = CL_PAROI_ADIABvsPAROI_ADIAB
        case("itf_fusion", "itf_fusion_sans_solide")
           condition = CL_PAROI_ADIABvsINJECTION
        case("periodique")
           condition = CL_PAROI_ADIABvsPERIODIQUE
        case default
           erreur=.true.
        end select
        
     case("itf_fusion", "itf_fusion_sans_solide")
        select case (trim(bc2))
        case("interbloc")
           condition = CL_INJECTIONvsINTERBLOC
        case("interproc")
           condition = CL_INTERPROCvs
        case("paroi_global")
           condition = CL_INJECTIONvsPAROI
        case("paroi_global_adiabatique")
           condition = CL_INJECTIONvsPAROI_ADIAB
        case default
           erreur=.true.
        end select
        
     case("itf_couplage")
        !!! TODO : doit dependre de la nature 
        select case (trim(bc2))
        case("flux_nul")
           condition = CL_PAROIvsFLUX_NUL
        case("symetrie")
           condition = CL_PAROIvsSYMETRIE
        case("interbloc")
           condition = CL_PAROIvsINTERBLOC
        case("interproc")
           condition = CL_INTERPROCvs
        case("entree_supersonique")
           condition = CL_PAROIvsENTREE_SUPER
        case("entree_subsonique")
           condition = CL_PAROIvsENTREE_SUB
        case("entree_subsonique_t")  
           erreur=.true.
        case ("sortie_subsonique")
           condition = CL_PAROIvsSORTIE_SUB
        case("paroi_global")
           condition = CL_PAROIvsPAROI
        case("paroi_global_adiabatique")
           condition = CL_PAROIvsPAROI_ADIAB
        case default
           erreur=.true.
        end select
     case default
        erreur=.true.
     end select
     
     
     if (erreur) then
        print*, "condition_coin : cas non gere ", trim(bc1), " ",trim(bc2)
       ! call arret_code
     end if

   end subroutine condition_coin
   
end module m_gestion_interbloc




