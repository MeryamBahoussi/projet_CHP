module SOLVEUR_PETSC
  use m_struct
  use parametres_globaux
  use parametres_fluide
#ifdef PETSC_ACTIF
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
  use petscvec
  use petscmat
  use petscksp
  use m_MPI, only : numproc, bloc_est_local, arret_code
#else
  use m_MPI
#endif

  implicit none

  type STR_BORNES         !   definition des bornes de calcul
     integer                             :: ld
     integer                             :: lf
     integer                             :: md
     integer                             :: mf
  end type STR_BORNES

contains

  !****************************************************************
  !*                                                              *
  !*                        CONFIGURATION                         *
  !*                                                              *
  !****************************************************************
  
  ! Initialisation
  subroutine init_petsc(tache)

    type (STR_TACHE), intent(in) :: tache

#ifdef PETSC_ACTIF
    ! Variables locales
    PetscErrorCode                             ierr
    type (STR_SUPER_BLOC),  pointer           :: sb
    type (STR_BLOC)     ,   pointer           :: b
    integer                                   :: u, ib, mb, ibi, isb

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    PETSC_COMM_WORLD = MPI_COMM_WORLD
    
    ! TODO il faudrait mettre la structure petsc au niveau du bloc et non du super bloc
    do mb = 1, tache%nb_membre
       isb = tache%membre(mb)
       sb => acces_super_bloc(isb)%pointeur
       
       if (sb%solveur_implicite==PETSC) then
          do ib = 1, sb%nb_membre
             ibi = sb%membre(ib)
             b  => acces_bloc(ibi)%pointeur
             if (bloc_est_local(b%numero)) then
                if (b%implicite) then
                   call init_solveur_global(b, sb, (/2,1/), b%mat%dim_bloc, b%mat%stencil, sb%inv_up%petsc)
                end if
             end if
          end do
       end if
    end do
#endif
  end subroutine init_petsc
 

  !****************************************************************
  !*                                                              *
  !*                    SOLVEUR GLOBAL 2D                         *
  !*                                                              *
  !****************************************************************

  ! Initialisation du solveur
  subroutine init_solveur_global(b_i, s_b, directions, dim_Uc, stencil, petsc)

    ! Arguments
    type (STR_SUPER_BLOC)                , pointer       :: s_b
    type (STR_BLOC)                      , pointer       :: b_i
    integer       , dimension(2)         , intent(in)    :: directions
    integer                              , intent(in)    :: dim_Uc
    integer                              , intent(in)    :: stencil
    type (STR_PETSC)                     , pointer       :: petsc

#ifdef PETSC_ACTIF
    ! Variables locales
    type (STR_BORNES)         :: bn
    type (STR_PATCH), pointer :: patch
    PetscErrorCode               ierr
    PetscInt                     i, l, m, ind_nnz, d_nnz_tmp
    PetscInt                     f, p, ind_l, ind_m, nb_elts

    ! Matrice structure BAIJ
    PetscInt                block_size          ! taille des blocs
    PetscInt                n_loc               ! taille de la matrice locale
    PetscInt                n_glob              ! taille de la matrice locale
    PetscInt                d_nz                ! nombre de blocs nonzeros par bloc-rangee dans la matrice diagonale locale
    PetscInt, pointer ::    d_nnz(:)            ! nombre de blocs nonzeros pour chaque bloc-rangee dans la matrice diagonale locale
    PetscInt                o_nz                ! nombre de blocs nonzeros par bloc-rangee dans la matrice extra-diagonale locale
    PetscInt, pointer ::    o_nnz(:)            ! nombre de blocs nonzeros pour chaque bloc-rangee dans la matrice extra-diagonale locale
    !PetscInt                its                 ! iterations for convergence
    

    bn = STR_BORNES(b_i%ld, b_i%lf, b_i%md, b_i%mf)
    
    ! Initialisation des dimensions de la matrice
    n_loc = dim_Uc * bn%lf * bn%mf
    n_glob = 0
    do i=1, size(acces_bloc,1)
       if (associated(acces_bloc(i)%pointeur) .and. (acces_bloc(i)%pointeur%nature == s_b%nature)) then
          n_glob = n_glob + acces_bloc(i)%pointeur%lf * acces_bloc(i)%pointeur%mf
       end if
    end do
    n_glob = n_glob * dim_Uc
    block_size = dim_Uc
    d_nz = 0 ! cette variable est ignoree car la variable d_nnz est utilisee   
    ALLOCATE(d_nnz(bn%lf * bn%mf)) ! la taille est egale au nombre de blocs ligne

    ! Calcul du tableau d_nnz
    select case (directions(1)) ! direction principale
    case (DIR_X) ! DIMENSION 1 | 3 | 2
       do l = bn%ld, bn%lf
          do m = bn%md, bn%mf
             ind_nnz = m + (l-1) * bn%mf
             d_nnz(ind_nnz) = 5
             if (l == bn%ld) d_nnz(ind_nnz) = d_nnz(ind_nnz) - 1
             if (l == bn%lf) d_nnz(ind_nnz) = d_nnz(ind_nnz) - 1
             if (m == bn%md) d_nnz(ind_nnz) = d_nnz(ind_nnz) - 1
             if (m == bn%mf) d_nnz(ind_nnz) = d_nnz(ind_nnz) - 1
             if (stencil == 9) then
                if ((bn%lf /= 1) .and. (bn%mf /= 1)) then 
                   d_nnz_tmp = 4
                   if (l == bn%ld) d_nnz_tmp = d_nnz_tmp - 2
                   if (l == bn%lf) d_nnz_tmp = d_nnz_tmp - 2
                   if (m == bn%md) d_nnz_tmp = d_nnz_tmp / 2
                   if (m == bn%mf) d_nnz_tmp = d_nnz_tmp / 2
                   d_nnz(ind_nnz) = d_nnz(ind_nnz) + d_nnz_tmp
                end if
             end if
          end do
       end do
    case (DIR_Y) ! DIMENSION 2 | 3 | 1         
       do m = bn%md, bn%mf
          do l = bn%ld, bn%lf
             ind_nnz = l + (m-1) * bn%lf
             d_nnz(ind_nnz) = 5
             if (l == bn%ld) d_nnz(ind_nnz) = d_nnz(ind_nnz) - 1
             if (l == bn%lf) d_nnz(ind_nnz) = d_nnz(ind_nnz) - 1
             if (m == bn%md) d_nnz(ind_nnz) = d_nnz(ind_nnz) - 1
             if (m == bn%mf) d_nnz(ind_nnz) = d_nnz(ind_nnz) - 1
             if (stencil == 9) then
                if ((bn%lf /= 1) .and. (bn%mf /= 1)) then 
                   d_nnz_tmp = 4
                   if (l == bn%ld) d_nnz_tmp = d_nnz_tmp - 2
                   if (l == bn%lf) d_nnz_tmp = d_nnz_tmp - 2
                   if (m == bn%md) d_nnz_tmp = d_nnz_tmp / 2
                   if (m == bn%mf) d_nnz_tmp = d_nnz_tmp / 2
                   d_nnz(ind_nnz) = d_nnz(ind_nnz) + d_nnz_tmp
                end if
             end if
          end do
       end do
    end select

    o_nz = 0 ! cette variable est ignoree car la variable o_nnz est utilisee
    ALLOCATE(o_nnz(bn%lf * bn%mf)) ! la taille est egale au nombre de blocs ligne

    ! Calcul du tableau o_nnz
    o_nnz = 0
    do f = 1, 4
       do p = 1, b_i%face(f)%nb_patch
          patch => b_i%face(f)%patch(p)
          if (patch%bc%condition=="interbloc" .or. patch%bc%condition=="periodique") then
             nb_elts = patch%ifin - patch%ideb + 1

             do i = 1, nb_elts
                select case (f)
                case (1)
                   ind_l = 1
                   ind_m = patch%ideb + i-1
                case (2)
                   ind_l = bn%lf
                   ind_m = patch%ideb + i-1
                case (3)
                   ind_l = patch%ideb + i-1
                   ind_m = 1
                case (4)
                   ind_l = patch%ideb + i-1
                   ind_m = bn%mf
                end select

                select case (directions(1)) ! direction principale
                case (DIR_X) ! DIMENSION 1 | 3 | 2
                   ind_nnz = ind_m + (ind_l-1) * bn%mf
                case (DIR_Y) ! DIMENSION 2 | 3 | 1 
                   ind_nnz = ind_l + (ind_m-1) * bn%lf
                end select
                if (patch%ib_adj == b_i%numero) then
                   ! on ajoute a la matrice diagonale locale
                   d_nnz(ind_nnz) = d_nnz(ind_nnz) + 1
                   if (stencil == 9) then
                      if (nb_elts /= 1) d_nnz(ind_nnz) = d_nnz(ind_nnz) + 1
                      if (i /= 1 .and. i /= nb_elts) d_nnz(ind_nnz) = d_nnz(ind_nnz) + 1
                   end if
                else
                   ! on ajoute a la matrice extra-diagonale locale
                   o_nnz(ind_nnz) = o_nnz(ind_nnz) + 1
                   if (stencil == 9) then
                      if (nb_elts /= 1) o_nnz(ind_nnz) = o_nnz(ind_nnz) + 1
                      if (i /= 1 .and. i /= nb_elts) o_nnz(ind_nnz) = o_nnz(ind_nnz) + 1
                   end if
                end if               
             end do
          end if
       end do
    end do

    ! Creation de la matrice distribuee
    call MatCreateBAIJ(s_b%mpi_comm, block_size, n_loc, n_loc, n_glob, n_glob, d_nz, d_nnz, o_nz, o_nnz, petsc%mat_A, ierr)
    call MatSetFromOptions(petsc%mat_A, ierr)

    DEALLOCATE(d_nnz)
    DEALLOCATE(o_nnz)

    ! Creation des vecteurs distribues
    call VecCreateMPI(s_b%mpi_comm, n_loc, n_glob, petsc%vec_x, ierr)
    call VecSetFromOptions(petsc%vec_x, ierr)
    call VecDuplicate(petsc%vec_x, petsc%vec_b, ierr)

    ! Creation du contexte du solveur lineaire
    call KSPCreate(s_b%mpi_comm, petsc%ksp, ierr)
    call KSPSetType(petsc%ksp, KSPGMRES, ierr)

    ! Configuration du solveur et du preconditionneur  
    ! Ici, la matrice definissant le syteme lineaire
    ! sert aussi de preconditionneur
    call KSPGetPC(petsc%ksp, petsc%ksp_pc, ierr)
    call PCSetType(petsc%ksp_pc, PCBJACOBI, ierr)

    ! Options GMRES
    call KSPSetTolerances(petsc%ksp, petsc%tol, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, petsc%kmax, ierr)

    call PCFactorSetReuseOrdering(petsc%ksp_pc, PETSC_TRUE, ierr)

    call KSPSetFromOptions(petsc%ksp, ierr)
#endif
  end subroutine init_solveur_global

  ! Solveur global utilisant PETSc
  subroutine RL_petsc_global(b_i, uc, s_b, matrice, sm, directions, petsc, iter)

    ! Arguments
    type (STR_SUPER_BLOC)                , pointer       :: s_b
    type (STR_BLOC)                      , pointer       :: b_i
    type (STR_MATRICE)                   , pointer       :: matrice
    real*8        , dimension(:,-1:,-1:) , intent(inout) :: uc
    real*8        , dimension(:, :, :)   , intent(inout) :: sm
    integer       , dimension(2)         , intent(in)    :: directions
    type (STR_PETSC)                     , pointer       :: petsc
    integer                              , intent(in)    :: iter

#ifdef PETSC_ACTIF
    integer :: dim_Uc
    type (STR_BORNES) :: bn
    type (STR_PATCH), pointer :: patch
    type (STR_BLOC) , pointer :: b_i_adj
    PetscErrorCode               ierr
    PetscInt                     i, j, k, l, m
    PetscScalar, pointer      :: val_tab(:)
    PetscInt                     ind_i, ind_j, offset, offset_adj
    PetscInt                     f, p, ind_l, ind_m, ind_l_adj, ind_m_adj, nb_elts, signe
    PetscScalar, pointer :: b_array(:)
    PetscScalar, pointer :: x_array(:)
    ! Resolution
    PetscInt                its                 ! iterations for convergence
    PetscInt                its_fin             ! iterations for convergence
    PetscScalar             val
    real*8               :: tol, tol_MPI

    PetscInt                Istart, Iend
    KSPConvergedReason      reason


    dim_Uc = matrice%dim_bloc
    bn = STR_BORNES(b_i%ld, b_i%lf, b_i%md, b_i%mf)

    ! si c'est la première itération du Newton, on doit construire la matrice
    ! et la preconditionner, sinon on a pas besoin car PETSc la garde au chaud,
    ! il faut juste reconstruire le second membre
    if(iter==1)then

       offset = 0 ! Calcul du decalage du bloc diagonal
       do i=1, size(acces_bloc,1)
          if (associated(acces_bloc(i)%pointeur) .and. (acces_bloc(i)%pointeur%nature == s_b%nature)) then
             if (acces_bloc(i)%pointeur%numproc == numproc) then
                exit
             else
                offset = offset + acces_bloc(i)%pointeur%lf * acces_bloc(i)%pointeur%mf
             end if
          end if
       end do

       allocate(val_tab(dim_Uc * dim_Uc))

       ! Remplissage de la matrice
       select case (directions(1)) ! direction principale
       case (DIR_X) ! DIMENSION 1 | 3 | 2

          call addDiagMatPetscDirX(bn,0, 0,0, 0,dim_Uc,offset,matrice%b ,     0,petsc%mat_A)     !b
          call addDiagMatPetscDirX(bn,0, 0,1, 0,dim_Uc,offset,matrice%a ,    -1,petsc%mat_A)     !a
          call addDiagMatPetscDirX(bn,0, 0,0,-1,dim_Uc,offset,matrice%c ,     1,petsc%mat_A)     !c
          call addDiagMatPetscDirX(bn,1, 0,0, 0,dim_Uc,offset,matrice%d ,-bn%mf,petsc%mat_A)     !d
          call addDiagMatPetscDirX(bn,0,-1,0, 0,dim_Uc,offset,matrice%e , bn%mf,petsc%mat_A)     !e
          if (matrice%stencil == 9) then
             call addDiagMatPetscDirX(bn,0,-1,1, 0,dim_Uc,offset,matrice%ae, bn%mf-1,petsc%mat_A)   !ae
             call addDiagMatPetscDirX(bn,1, 0,1, 0,dim_Uc,offset,matrice%ad,-bn%mf-1,petsc%mat_A)   !ad
             call addDiagMatPetscDirX(bn,1, 0,0,-1,dim_Uc,offset,matrice%cd,-bn%mf+1,petsc%mat_A)   !cd
             call addDiagMatPetscDirX(bn,0,-1,0,-1,dim_Uc,offset,matrice%ce ,bn%mf+1,petsc%mat_A)   !ce
          end if

       case (DIR_Y) ! DIMENSION 2 | 3 | 1

          call addDiagMatPetscDirY(bn,0, 0,0, 0,dim_Uc,offset,matrice%b ,     0,petsc%mat_A)     !b
          call addDiagMatPetscDirY(bn,1, 0,0, 0,dim_Uc,offset,matrice%a ,    -1,petsc%mat_A)     !a
          call addDiagMatPetscDirY(bn,0,-1,0, 0,dim_Uc,offset,matrice%c ,     1,petsc%mat_A)     !c
          call addDiagMatPetscDirY(bn,0, 0,1, 0,dim_Uc,offset,matrice%d ,-bn%lf,petsc%mat_A)     !d
          call addDiagMatPetscDirY(bn,0, 0,0,-1,dim_Uc,offset,matrice%e , bn%lf,petsc%mat_A)     !e
          if (matrice%stencil == 9) then
             call addDiagMatPetscDirY(bn,1, 0,0,-1,dim_Uc,offset,matrice%ae, bn%lf-1,petsc%mat_A)     !ae
             call addDiagMatPetscDirY(bn,1, 0,1, 0,dim_Uc,offset,matrice%ad,-bn%lf-1,petsc%mat_A)     !ad
             call addDiagMatPetscDirY(bn,0,-1,0,-1,dim_Uc,offset,matrice%ce, bn%lf+1,petsc%mat_A)     !ce
             call addDiagMatPetscDirY(bn,0,-1,1, 0,dim_Uc,offset,matrice%cd,-bn%lf+1,petsc%mat_A)     !cd
          end if

       end select ! directions(1)

       ! Coefficients pour les liaisons interblocs
       do f = 1, 4
          do p = 1, b_i%face(f)%nb_patch
             patch => b_i%face(f)%patch(p)
             if (patch%bc%condition=="interbloc" .or. patch%bc%condition=="periodique") then
                nb_elts = patch%ifin - patch%ideb + 1
                offset_adj = 0
                do i=1, size(acces_bloc,1)
                   if (associated(acces_bloc(i)%pointeur) .and. (acces_bloc(i)%pointeur%nature == s_b%nature)) then
                      if (acces_bloc(i)%pointeur%numero == patch%ib_adj) then
                         b_i_adj => acces_bloc(i)%pointeur
                         exit
                      else
                         offset_adj = offset_adj + acces_bloc(i)%pointeur%lf * acces_bloc(i)%pointeur%mf
                      end if
                   end if
                end do

                do i = 1, nb_elts
                   ! Calcul des indices de la maille du bloc b
                   select case (f)
                   case (1)
                      ind_l = 1
                      ind_m = patch%ideb + i-1
                   case (2)
                      ind_l = bn%lf
                      ind_m = patch%ideb + i-1
                   case (3)
                      ind_l = patch%ideb + i-1
                      ind_m = 1
                   case (4)
                      ind_l = patch%ideb + i-1
                      ind_m = bn%mf
                   end select

                   if (patch%sens < 0) then
                      signe = -1
                   else
                      signe =  1
                   end if

                   ! Calcul des indices de la maille en vis a vis dans le bloc b_i_adj
                   select case (patch%if_adj)
                   case (1)
                      ind_l_adj = 1
                      ind_m_adj = patch%ideb_adj + signe * (i-1)
                   case (2)
                      ind_l_adj = b_i_adj%lf
                      ind_m_adj = patch%ideb_adj + signe * (i-1)
                   case (3)
                      ind_l_adj = patch%ideb_adj + signe * (i-1)
                      ind_m_adj = 1
                   case (4)
                      ind_l_adj = patch%ideb_adj + signe * (i-1)
                      ind_m_adj = b_i_adj%mf
                   end select

                   ! Calcul de l indice de ligne avec les indices du bloc b
                   select case (directions(1)) ! direction principale
                   case (DIR_X) ! DIMENSION 1 | 3 | 2
                      ind_i = offset + ind_m-1 + (ind_l-1) * bn%mf
                   case (DIR_Y) ! DIMENSION 2 | 3 | 1
                      ind_i = offset + ind_l-1 + (ind_m-1) * bn%lf
                   end select

                   ! Calcul de l indice de colonne avec les indices du bloc b_i_adj
                   ! select case (b_i_adj%principale) ! direction principale
                   select case (directions(1))
                   case (DIR_X) ! DIMENSION 1 | 3 | 2
                      ind_j = offset_adj + ind_m_adj-1 + (ind_l_adj-1) * b_i_adj%mf
                   case (DIR_Y) ! DIMENSION 2 | 3 | 1
                      ind_j = offset_adj + ind_l_adj-1 + (ind_m_adj-1) * b_i_adj%lf
                   end select

                   do j = 1, dim_Uc
                      do k = 1, dim_Uc
                         select case (f)
                         case (1)
                            if (directions(1) == DIR_X) then
                               val_tab(k + dim_Uc * (j-1)) = matrice%d(j,k,ind_l,ind_m)
                            else
                               val_tab(k + dim_Uc * (j-1)) = matrice%a(j,k,ind_l,ind_m)
                            end if
                         case (2)
                            if (directions(1) == DIR_X) then
                               val_tab(k + dim_Uc * (j-1)) = matrice%e(j,k,ind_l,ind_m)
                            else
                               val_tab(k + dim_Uc * (j-1)) = matrice%c(j,k,ind_l,ind_m)
                            end if
                         case (3)
                            if (directions(1) == DIR_X) then
                               val_tab(k + dim_Uc * (j-1)) = matrice%a(j,k,ind_l,ind_m)
                            else
                               val_tab(k + dim_Uc * (j-1)) = matrice%d(j,k,ind_l,ind_m)
                            end if
                         case (4)
                            if (directions(1) == DIR_X) then
                               val_tab(k + dim_Uc * (j-1)) = matrice%c(j,k,ind_l,ind_m)
                            else
                               val_tab(k + dim_Uc * (j-1)) = matrice%e(j,k,ind_l,ind_m)
                            end if
                         end select
                      end do
                   end do

                   call MatSetValuesBlocked(petsc%mat_A, 1, ind_i, 1, ind_j, val_tab, INSERT_VALUES, ierr)


                   if (matrice%stencil == 9) then
                      if (nb_elts /= 1) then
                         if (i /= 1) then 
                            select case (patch%if_adj)
                            case (1)
                               ind_l_adj = 1
                               ind_m_adj = patch%ideb_adj + signe * (i-2)
                            case (2)
                               ind_l_adj = b_i_adj%lf
                               ind_m_adj = patch%ideb_adj + signe * (i-2)
                            case (3)
                               ind_l_adj = patch%ideb_adj + signe * (i-2)
                               ind_m_adj = 1
                            case (4)
                               ind_l_adj = patch%ideb_adj + signe * (i-2)
                               ind_m_adj = b_i_adj%mf
                            end select

                            select case (directions(1)) ! direction principale
                            case (DIR_X) ! DIMENSION 1 | 3 | 2
                               ind_i = offset + ind_m-1 + (ind_l-1) * bn%mf
                            case (DIR_Y) ! DIMENSION 2 | 3 | 1
                               ind_i = offset + ind_l-1 + (ind_m-1) * bn%lf
                            end select

                            select case (directions(1))
                            case (DIR_X) ! DIMENSION 1 | 3 | 2
                               ind_j = offset_adj + ind_m_adj-1 + (ind_l_adj-1) * b_i_adj%mf
                            case (DIR_Y) ! DIMENSION 2 | 3 | 1
                               ind_j = offset_adj + ind_l_adj-1 + (ind_m_adj-1) * b_i_adj%lf
                            end select

                            ! a confirmer ?
                            do j = 1, dim_Uc
                               do k = 1, dim_Uc
                                  select case (f)
                                  case (1)
                                     if (directions(1) == DIR_X) then
                                        val_tab(k + dim_Uc * (j-1)) = matrice%ad(j,k,ind_l,ind_m)
                                     else
                                        val_tab(k + dim_Uc * (j-1)) = matrice%ad(j,k,ind_l,ind_m)
                                     end if
                                  case (2)
                                     if (directions(1) == DIR_X) then
                                        val_tab(k + dim_Uc * (j-1)) = matrice%ae(j,k,ind_l,ind_m)
                                     else
                                        val_tab(k + dim_Uc * (j-1)) = matrice%cd(j,k,ind_l,ind_m)
                                     end if
                                  case (3)
                                     if (directions(1) == DIR_X) then
                                        val_tab(k + dim_Uc * (j-1)) = matrice%ad(j,k,ind_l,ind_m)
                                     else
                                        val_tab(k + dim_Uc * (j-1)) = matrice%ad(j,k,ind_l,ind_m)
                                     end if
                                  case (4)
                                     if (directions(1) == DIR_X) then
                                        val_tab(k + dim_Uc * (j-1)) = matrice%cd(j,k,ind_l,ind_m)
                                     else
                                        val_tab(k + dim_Uc * (j-1)) = matrice%ae(j,k,ind_l,ind_m)
                                     end if
                                  end select
                               end do
                            end do
                            call MatSetValuesBlocked(petsc%mat_A, 1, ind_i, 1, ind_j, val_tab, INSERT_VALUES, ierr)
                         end if

                         if (i /= nb_elts) then
                            select case (patch%if_adj)
                            case (1)
                               ind_l_adj = 1
                               ind_m_adj = patch%ideb_adj + signe * (i)
                            case (2)
                               ind_l_adj = b_i_adj%lf
                               ind_m_adj = patch%ideb_adj + signe * (i)
                            case (3)
                               ind_l_adj = patch%ideb_adj + signe * (i)
                               ind_m_adj = 1
                            case (4)
                               ind_l_adj = patch%ideb_adj + signe * (i)
                               ind_m_adj = b_i_adj%mf
                            end select

                            select case (directions(1)) ! direction principale
                            case (DIR_X) ! DIMENSION 1 | 3 | 2
                               ind_i = offset + ind_m-1 + (ind_l-1) * bn%mf
                            case (DIR_Y) ! DIMENSION 2 | 3 | 1
                               ind_i = offset + ind_l-1 + (ind_m-1) * bn%lf
                            end select

                            select case (directions(1))
                            case (DIR_X) ! DIMENSION 1 | 3 | 2
                               ind_j = offset_adj + ind_m_adj-1 + (ind_l_adj-1) * b_i_adj%mf
                            case (DIR_Y) ! DIMENSION 2 | 3 | 1
                               ind_j = offset_adj + ind_l_adj-1 + (ind_m_adj-1) * b_i_adj%lf
                            end select

                            do j = 1, dim_Uc
                               do k = 1, dim_Uc
                                  select case (f)
                                  case (1)
                                     if (directions(1) == DIR_X) then
                                        val_tab(k + dim_Uc * (j-1)) = matrice%cd(j,k,ind_l,ind_m)
                                     else
                                        val_tab(k + dim_Uc * (j-1)) = matrice%ae(j,k,ind_l,ind_m)
                                     end if
                                  case (2)
                                     if (directions(1) == DIR_X) then
                                        val_tab(k + dim_Uc * (j-1)) = matrice%ce(j,k,ind_l,ind_m)
                                     else
                                        val_tab(k + dim_Uc * (j-1)) = matrice%ce(j,k,ind_l,ind_m)
                                     end if
                                  case (3)
                                     if (directions(1) == DIR_X) then
                                        val_tab(k + dim_Uc * (j-1)) = matrice%ae(j,k,ind_l,ind_m)
                                     else
                                        val_tab(k + dim_Uc * (j-1)) = matrice%cd(j,k,ind_l,ind_m)
                                     end if
                                  case (4)
                                     if (directions(1) == DIR_X) then
                                        val_tab(k + dim_Uc * (j-1)) = matrice%ce(j,k,ind_l,ind_m)
                                     else
                                        val_tab(k + dim_Uc * (j-1)) = matrice%ce(j,k,ind_l,ind_m)
                                     end if
                                  end select
                               end do
                            end do
                            call MatSetValuesBlocked(petsc%mat_A, 1, ind_i, 1, ind_j, val_tab, INSERT_VALUES, ierr)
                         end if
                      end if
                   end if
                end do
             end if
          end do
       end do


       ! Assemblage de la matrice
       call MatAssemblyBegin(petsc%mat_A, MAT_FINAL_ASSEMBLY, ierr)
       call MatAssemblyEnd(petsc%mat_A, MAT_FINAL_ASSEMBLY, ierr)

       !!! Je ne sais pas a quoi sert cette routine. Elle n'est plus appelee dans le code maison  
       ! call MatGetOwnershipRange(petsc%mat_A,Istart,Iend,ierr)

       ! Definition des Operators. Ici, la matrice definissant le syteme lineaire
       ! sert aussi de preconditionneur
       call KSPSetOperators(petsc%ksp, petsc%mat_A, petsc%mat_A, ierr)

    end if


    !call MatView(petsc%mat_A, PETSC_VIEWER_STDOUT_WORLD, ierr)

    ! Remplissage du vecteur second membre
    select case (directions(1)) ! direction principale
    case (DIR_X) ! DIMENSION 1 | 3 | 2
       call VecGetArrayF90(petsc%vec_b, b_array, ierr)
       do l = bn%ld, bn%lf
          do m = bn%md, bn%mf
             do i = 1, dim_Uc
                ind_i = i + (m-1) * dim_Uc + bn%mf * (l-1) * dim_Uc
                b_array(ind_i) = sm(i,l,m)
             end do
          end do
       end do
       call VecRestoreArrayF90(petsc%vec_b, b_array, ierr)
    case (DIR_Y) ! DIMENSION 2 | 3 | 1
       call VecGetArrayF90(petsc%vec_b, b_array, ierr)
       do m = bn%md, bn%mf
          do l = bn%ld, bn%lf
             do i = 1, dim_Uc
                ind_i = i + (l-1) * dim_Uc + bn%lf * (m-1) * dim_Uc
                b_array(ind_i) = sm(i,l,m)
             end do
          end do
       end do
       call VecRestoreArrayF90(petsc%vec_b, b_array, ierr)
    end select ! directions(1)



    ! Assemblage du vecteur
    call VecAssemblyBegin(petsc%vec_b, ierr)
    call VecAssemblyEnd(petsc%vec_b, ierr)

    !call VecView(vec_b, PETSC_VIEWER_STDOUT_WORLD, ierr)

    ! Resolution du systeme
    call KSPSolve(petsc%ksp, petsc%vec_b, petsc%vec_x, ierr)

!!$    ! Debug : on regarde si on a atteint la convergence
!!$    call KSPGetConvergedReason(ksp, reason);
!!$    if (reason==KSP_DIVERGED_INDEFINITE_PC) then
!!$       call PetscPrintf(sb%mpi_comm,"\nDivergence because of indefinite preconditioner\n", ierr)
!!$       call PetscPrintf(sb%mpi_comm,"\nRun the executable again but with '-pc_factor_shift_type POSITIVE_DEFINITE' option.\n", ierr)
!!$    else if (reason<0) then
!!$       call PetscPrintf(sb%mpi_comm,"\nOther kind of divergence: this should not happen.\n",ierr)
!!$       call KSPGetResidualNorm(ksp, tol, ierr)
!!$       call MPI_ALLREDUCE(tol, tol_MPI, 1,MPI_REAL8,MPI_MAX,sb%mpi_comm,ierr)
!!$       if (numproc==0) print*, "residu à la derniere iteration :", tol_MPI
!!$    else
!!$       call KSPGetIterationNumber(ksp,its_fin, ierr)
!!$       call KSPGetResidualNorm(ksp, tol, ierr)
!!$       call MPI_ALLREDUCE(tol, tol_MPI, 1,MPI_REAL8,MPI_MAX,sb%mpi_comm,ierr)
!!$       if (its_fin == its) then
!!$          if (numproc==0) then
!!$             print '("Toujours pas converge apres :", i4," iterations.")', its
!!$             print*, "residu à la derniere iteration :", tol_MPI
!!$          end if
!!$       else
!!$          !if (numproc==0 .and. dim_Uc == 3) print*, "convergence en", its_fin," iterations."
!!$          ! c'est converge
!!$       end if
!!$    end if

    ! call VecView(vec_x, PETSC_VIEWER_STDOUT_WORLD, ierr)

    ! Ecriture de vec_x dans uc
    select case (directions(1)) ! direction principale
    case (DIR_X) ! DIMENSION 1 | 3 | 2
       call VecGetArrayF90(petsc%vec_x, x_array, ierr)
       do l = bn%ld, bn%lf
          do m = bn%md, bn%mf
             do i = 1, dim_Uc
                ind_i = i + (m-1) * dim_Uc + bn%mf * (l-1) * dim_Uc
                uc(i,l,m) = uc(i,l,m) + x_array(ind_i)
             end do
          end do
       end do
       call VecRestoreArrayF90(petsc%vec_x, x_array, ierr)
    case (DIR_Y) ! DIMENSION 2 | 3 | 1
       call VecGetArrayF90(petsc%vec_x, x_array, ierr)
       do m = bn%md, bn%mf
          do l = bn%ld, bn%lf
             do i = 1, dim_Uc
                ind_i = i + (l-1) * dim_Uc + bn%lf * (m-1) * dim_Uc
                uc(i,l,m) = uc(i,l,m) + x_array(ind_i)
             end do
          end do
       end do
       call VecRestoreArrayF90(petsc%vec_x, x_array, ierr)
    end select ! directions(1)

    ! on DEALLOCATE seulement lorsqu'on construit la matrice (iter==1)
    if(iter==1) DEALLOCATE(val_tab)

#endif
  end subroutine RL_petsc_global

#ifdef PETSC_ACTIF
  ! Routine generique pour ajouter les coefficients dans la matrice globale de petsc
  subroutine addDiagMatPetscDirX(bn,dld,dlf,dmd,dmf,dim_Uc,offset,diag,decalagej,mat_A)
#include "petsc/finclude/petsc.h"
 
    type (STR_BORNES)         , intent(in) :: bn              ! dimension du bloc
    integer                   , intent(in) :: dld,dlf,dmd,dmf ! decalages d'indices de debut et de fin pour la diagonale consideree
    integer                   , intent(in) :: dim_Uc          ! dimension du bloc de la matrice
    integer                   , intent(in) :: offset          ! du decalage du bloc diagonal
    integer                   , intent(in) :: decalagej       ! decalage de la diagonale dans le stockage de petsc
    real*8, dimension(:,:,:,:), pointer    :: diag            ! diagonale de la matrice a ajouter
    Mat                                       mat_A           ! matrice petsc

    ! variables locales
    integer :: l, m, i, j
    integer :: ind_i, ind_j
    integer :: ierr
    PetscScalar, pointer :: val_tab(:)   
    
    allocate(val_tab(dim_Uc * dim_Uc))

    do l = bn%ld+dld, bn%lf+dlf
       do m = bn%md+dmd, bn%mf+dmf
          do i = 1, dim_Uc
             do j = 1, dim_Uc
                val_tab(j + dim_Uc * (i-1)) = diag(i,j,l,m) 
             end do
          end do
          ind_i = offset + m-1 + bn%mf * (l-1)
          ind_j = offset + m-1 + bn%mf * (l-1) + decalagej
          call MatSetValuesBlocked(mat_A, 1, ind_i, 1, ind_j, val_tab, INSERT_VALUES, ierr)
       end do
    end do
    
    deallocate(val_tab)

  end subroutine addDiagMatPetscDirX
#endif

#ifdef PETSC_ACTIF
  ! Routine generique pour ajouter les coefficients dans la matrice globale de petsc
  subroutine addDiagMatPetscDirY(bn,dld,dlf,dmd,dmf,dim_Uc,offset,diag,decalagej,mat_A)
#include "petsc/finclude/petsc.h"
    
    type (STR_BORNES)         , intent(in) :: bn              ! dimension du bloc
    integer                   , intent(in) :: dld,dlf,dmd,dmf ! decalages d'indices de debut et de fin pour la diagonale consideree
    integer                   , intent(in) :: dim_Uc          ! dimension du bloc de la matrice
    integer                   , intent(in) :: offset          ! du decalage du bloc diagonal
    integer                   , intent(in) :: decalagej       ! decalage de la diagonale dans le stockage de petsc
    real*8, dimension(:,:,:,:), pointer    :: diag            ! diagonale de la matrice a ajouter
    Mat                                       mat_A           ! matrice petsc

    ! variables locales
    integer :: l, m, i, j
    integer :: ind_i, ind_j
    integer :: ierr
    PetscScalar, pointer :: val_tab(:)
    
    allocate(val_tab(dim_Uc * dim_Uc))
    do m = bn%md+dmd, bn%mf+dmf
       do l = bn%ld+dld, bn%lf+dlf
          do i = 1, dim_Uc
             do j = 1, dim_Uc
                val_tab(j + dim_Uc * (i-1)) = diag(i,j,l,m) 
             end do
          end do
          ind_i = offset + l-1 + bn%lf * (m-1)
          ind_j = offset + l-1 + bn%lf * (m-1) + decalagej
          call MatSetValuesBlocked(mat_A, 1, ind_i, 1, ind_j, val_tab, INSERT_VALUES, ierr)
       end do
    end do
  end subroutine addDiagMatPetscDirY
#endif

end module SOLVEUR_PETSC
