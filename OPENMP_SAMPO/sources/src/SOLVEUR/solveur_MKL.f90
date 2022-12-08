module solveur_MKL
  use m_struct
  use m_outils_CSR
  !use m_verification
  use m_communication_MPI
  
  implicit none

contains
  
  subroutine MKL_assemblage(bloc, matrice)
    implicit none 

    type (STR_BLOC), pointer :: bloc
    type (STR_MATRICE), pointer :: matrice
    
    integer :: taille_l, taille_m 
    integer :: dim_X

    real*8, dimension(:,:,:,:), pointer :: A_global, B_global, C_global
    real*8, dimension(:,:,:,:), pointer :: D_global, E_global, nul
    real*8, dimension(:,:,:,:), pointer :: AD_global, CD_global
    real*8, dimension(:,:,:,:), pointer :: AE_global, CE_global

    real*8, dimension(:,:,:,:), pointer :: A, B, C, D, E
    real*8, dimension(:,:,:,:), pointer :: AD, CD, AE, CE

    A => matrice%A; B => matrice%B
    C => matrice%C; D => matrice%D
    E => matrice%E
    
    dim_X = matrice%dim_bloc

    if (matrice%stencil == 5) then
       if (numproc == bloc%proc_ini) then
          allocate( nul(dim_X,dim_X,bloc%ld_tot:bloc%lf_tot,bloc%md_tot:bloc%mf_tot) )
          nul(:,:,:,:) = 0.d0
       end if
    else
       AD => matrice%AD; CD => matrice%CD
       AE => matrice%AE; CE => matrice%CE
    end if
    
    ! Envoi au proc ini si on est en //
    if (bloc%proc_fin /= bloc%proc_ini) then
       call MPI_send_to_proc_ini(bloc, dim_X, A, A_global)
       call MPI_send_to_proc_ini(bloc, dim_X, B, B_global)
       call MPI_send_to_proc_ini(bloc, dim_X, C, C_global)
       call MPI_send_to_proc_ini(bloc, dim_X, D, D_global)
       call MPI_send_to_proc_ini(bloc, dim_X, E, E_global)
       if (matrice%stencil == 9) then
          call MPI_send_to_proc_ini(bloc, dim_X, AD, AD_global)
          call MPI_send_to_proc_ini(bloc, dim_X, CD, CD_global)
          call MPI_send_to_proc_ini(bloc, dim_X, AE, AE_global)
          call MPI_send_to_proc_ini(bloc, dim_X, CE, CE_global)
       end if
    else
       A_global => A
       B_global => B
       C_global => C
       D_global => D
       E_global => E
       
       if (matrice%stencil == 9) then
          AD_global => AD
          CD_global => CD
          AE_global => AE
          CE_global => CE
       end if
    end if

    ! Resolution du systeme par le processeur ini
    if ( numproc == bloc%proc_ini ) then
       taille_l = bloc%lf_tot - bloc%ld_tot + 1
       taille_m = bloc%mf_tot - bloc%md_tot + 1

       if (matrice%stencil == 9) then
          call Remplissage_matrice_CSR( taille_l, taille_m, &
               bloc%mat_csr%values, bloc%mat_csr%nNonZeros, &
               A_global, B_global, C_global, D_global, E_global, &
               AD_global, CD_global, AE_global, CE_global, dim_X, .false.)
       else
          call Remplissage_matrice_CSR( taille_l, taille_m, &
               bloc%mat_csr%values, bloc%mat_csr%nNonZeros, &
               A_global, B_global, C_global, D_global, E_global, &
               nul, nul, nul, nul, dim_X, .false.)
       end if
    end if
    
    !! liberation de la memoire
    if (numproc == bloc%proc_ini) then
       if (bloc%proc_fin /= bloc%proc_ini) then
          deallocate(A_global, B_global, C_global, D_global, E_global)
          if (matrice%stencil == 9) deallocate(AD_global,AE_global,CD_global,CE_global)
       end if
       if (matrice%stencil == 5) deallocate(nul)
    end if
    
  end subroutine MKL_assemblage

  subroutine RL_MKL_sur_proc0(b, X, matrice, sm)

    type (STR_BLOC), pointer ::  b
    type (STR_MATRICE), pointer :: matrice
    real*8, dimension(matrice%dim_bloc,b%ld:b%lf,b%md:b%mf), intent(out) :: X
    real*8, dimension(:,:,:), pointer :: Sm 
    
    integer :: me, taille
    integer :: dim_X
    integer :: STATINFO, MPI_FLAG_sm = 100
    integer :: MPI_status(MPI_STATUS_SIZE)

    real*8, dimension(:,:,:), allocatable :: sm_global
    integer, dimension(:,:), pointer :: l_pere, m_pere

    dim_X = matrice%dim_bloc

    l_pere => b%l_charge
    m_pere => b%m_charge
    
    ! Envoi au proc 0 le second membre
    call MPI_send_to_proc_ini(b, dim_X, sm, sm_global, .false.)
        
    ! Resolution du systeme par le processeur ini
    if ( numproc == b%proc_ini ) then
       
       call solveur_MKL_for_CSR( b%lf_tot-b%ld_tot+1,b%mf_tot-b%md_tot+1, &
            b%mat_csr%rowIndex, b%mat_csr%columns, &
            b%mat_csr%values, Sm_global, b%mat_csr%nNonZeros, b%mat_csr%nb_Rows, dim_X)
       
       if (.not. delta_forme) then
          X(:,b%ld:b%lf,b%md:b%mf) = sm_global(:,b%ld:b%lf,b%md:b%mf)
       else
          X(:,b%ld:b%lf,b%md:b%mf) = X(:,b%ld:b%lf,b%md:b%mf) + &
               sm_global(:,b%ld:b%lf,b%md:b%mf)
       end if
    end if
    
    ! Diffusion de la solution aux autres processeurs (utiliser un broadcast ??)
    nbprocs = b%proc_fin - b%proc_ini + 1

    if (nbprocs > 1) then
       if ( numproc == b%proc_ini ) then
          do me = 1, nbprocs-1
             taille = dim_X*(l_pere(2,me)-l_pere(1,me)+1)*(m_pere(2,me)-m_pere(1,me)+1)

             call MPI_SEND( &
                  sm_global(:,l_pere(1,me):l_pere(2,me),m_pere(1,me):m_pere(2,me)) ,&
                  taille, MPI_REAL8, &
                  me, MPI_FLAG_sm, MPI_COMM_WORLD, STATINFO )
          end do
       else
          taille = dim_X*( b%mf - b%md + 1 ) * ( b%lf - b%ld + 1 )
          call MPI_RECV( sm(:,b%ld:b%lf,b%md:b%mf),&
               taille, MPI_REAL8, &
               0, MPI_FLAG_sm, MPI_COMM_WORLD, MPI_STATUS, STATINFO )
       
          if (.not. delta_forme) then
             X(:,b%ld:b%lf,b%md:b%mf) = sm(:,b%ld:b%lf,b%md:b%mf)
          else
             X(:,b%ld:b%lf,b%md:b%mf) = X(:,b%ld:b%lf,b%md:b%mf) + &
                  sm(:,b%ld:b%lf,b%md:b%mf)
          end if
       end if
    end if

    if ( numproc == b%proc_ini ) deallocate(sm_global)

  end subroutine RL_MKL_sur_proc0

end module SOLVEUR_MKL
