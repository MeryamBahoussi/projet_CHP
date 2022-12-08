module m_verification
  use m_struct
  use m_MPI
  use m_maillage, only : distance
  use synchronisation_interbloc
  use parametres_fluide 
  use parametres_globaux
  
  implicit none
  
  public verification_input

  private
contains
  
  ! verification des donnees
  subroutine verification_input(tache)
    type (STR_TACHE), intent(inout) :: tache
    
    type (STR_SUPER_BLOC), pointer :: sb
    type (STR_BLOC), pointer :: b
    integer :: nb_bloc_verif
    integer :: i, j
    integer :: isb, ib
    integer :: cartesien
    logical :: ordre2
    logical :: periodique
    logical :: nonUniformBc

!!! verification du nombre de blocs
    nb_bloc_verif=0
    do i = 1, tache%nb_membre
       isb = tache%membre(i)
       sb => acces_super_bloc(isb)%pointeur
       nb_bloc_verif = nb_bloc_verif + sb%nb_membre
    end do
    if (nb_bloc_verif /= ubound(acces_bloc,1)) then
       if (numproc==0) print*, "verification.f90 : Incoherence sur le nombre de bloc"
       call arret_code
    end if

!!! implicite en //  
    do i = 1, tache%nb_membre
       isb = tache%membre(i)
       sb => acces_super_bloc(isb)%pointeur
       
       if (sb%implicite) then
          if (sb%solveur_implicite==PETSC .and. nbprocs < sb%nb_membre) then
             if (numproc==0) print*, "verification.f90 : Pas assez de processeurs pour lancer ce cas avec PETSc", nbprocs, sb%nb_membre
             call arret_code
          else if (sb%solveur_implicite/=PETSC .and. sb%nb_membre_input/=1) then
             if (numproc==0) print*, "verification.f90 : Le solveur utilisé ne gère pas les cas multiblocs"
             call arret_code
          end if
       end if
    end do

!!! verification de la topologie des faces (interbloc ...)
    do i = 1, tache%nb_membre 
       isb = tache%membre(i)
       sb => acces_super_bloc(isb)%pointeur 
       do j = 1, sb%nb_membre 
          ib = sb%membre(j)
          if (bloc_est_local(ib)) then
             b=>acces_bloc(ib)%pointeur
             call verification_faces_bloc(b)
            ! call print_infos_bloc(b)
          end if
       end do
    end do

!!!!!!!!! coherence des parametres
    do i = 1, tache%nb_membre 
       isb = tache%membre(i)
       sb => acces_super_bloc(isb)%pointeur 
              
       select case (sb%nature)
       case (FLUIDE)
          ! choix du solveur, implicite
          select case(sb%solveur_Riemann)
          case(GALLICE, HLL, ROE)
          case default
             print*, "Solveur inconnu pour la resolution fluide DIRECTE explicite : ", sb%solveur_Riemann
             call arret_code
          end select
             
          call is_ordre2_fluide(sb, ordre2)
          if (ordre2 .and. sb%cfl > 0.5) then
             print*, "!!!!!!! Warning : Calcul a l'ordre 2 et cfl > 0.5 !!!!!!"
          end if

          call is_cartesien_fluide(sb, cartesien)
         
       case default
          print*, "nature de super bloc inconnue :", sb%nature
          call arret_code
       end select

!!! pour les deux natures (fluide et solide)
       ! les cas avec des conditions limites non uniformes ne sont pas geres en parallele
       ! il y a un problème dans la gestion du decoupage des patchs et des tableaux de condition limite 
       call have_nonuniform_bc(sb, nonUniformBc)
       if (nonUniformBc .and. nbprocs > sb%nb_membre_input) then
          print*, "Erreur : les conditions conditions limites non uniformes ne sont pas gerees en parallele"
          call arret_code
       end if
    end do

  end subroutine verification_input

  subroutine is_ordre2_fluide(sb, ordre2)
    type (STR_SUPER_BLOC), intent(in)  :: sb
    logical              , intent(out) :: ordre2
    
    ! variables locale
    integer :: ordre_acoustique, ordre_transport, ordre_Eulerdirect
    
    if (sb%nature == SOLIDE) return
 
    ! tous les blocs ont les memes ordre
    ordre_Eulerdirect = acces_bloc(sb%membre(1))%pointeur%ordre_Eulerdirect

    ! on regarde si on est a l'ordre 2
    ordre2 = .false.
    if (ordre_Eulerdirect==2) ordre2 = .true.
    
  end subroutine is_ordre2_fluide
  
  subroutine is_cartesien_fluide(sb, cartesien)
    type (STR_SUPER_BLOC), intent(in)  :: sb
    integer              , intent(out) :: cartesien
    
    ! variables locales
    integer                  :: i, ib
    type (STR_BLOC), pointer :: b

    if (sb%nature == SOLIDE) return
    
    cartesien = -1000
    do i = 1, sb%nb_membre
       ib = sb%membre(i)
       b=>acces_bloc(ib)%pointeur
       
       if ( index(b%maillage%type,"restart") /= 0 ) then
          cartesien = 2 ! on est en reprise donc on est pas sur que le maillage soit cartesien
       else
          select case(b%maillage%type)
          case ("cartesien")
             cartesien = 1
          case default
             cartesien = 0 ! on a un bloc qui n'a pas de maillage cartesien
          end select
       end if
    end do
  end subroutine is_cartesien_fluide

  subroutine have_nonuniform_bc(sb, nonUniformBc)
    type (STR_SUPER_BLOC), intent(in)  :: sb
    logical              , intent(out) :: nonUniformBc
    
    ! variables locales
    integer                               :: i, ib, iface, ip
    type (STR_BLOC)             , pointer :: b
    type (STR_BoundaryCondition), pointer :: bc

    nonUniformBc = .false.
    do i = 1, sb%nb_membre
       ib = sb%membre(i)
       b => acces_bloc(ib)%pointeur
       
       do iface = 1, 4
          do ip = 1, b%face(iface)%nb_patch
             bc => b%face(iface)%patch(ip)%bc
             if (bc%condition == "paroi_global") then
                if ( associated(bc%uw)) then
                   if (abs(maxval(bc%uw(:)) - minval(bc%uw(:))) >= 1.d-10) nonUniformBc = .true.
                end if
                if ( associated(bc%Tw)) then
                   if (abs(maxval(bc%Tw(:)) - minval(bc%Tw(:))) >= 1.d-10) nonUniformBc = .true.
                end if
                if ( associated(bc%flux_impose)) then
                   if (abs(maxval(bc%flux_impose(:)) - minval(bc%flux_impose(:))) >= 1.d-10) nonUniformBc = .true.
                end if
             end if
          end do
       end do
    end do
  end subroutine have_nonuniform_bc

  ! verification des infos lues pour une face
  ! on regarde la coherence des
  !   indices du patch vis a vis de ceux du bloc
  !   conditions de type interbloc
  subroutine verification_faces_bloc(b)
    use m_MPI
    implicit none

    type (STR_BLOC), target, intent(in) :: b
    type (STR_BLOC), pointer :: b_adj
    type (STR_FACE), pointer :: face, face_adj
    type (STR_PATCH), pointer :: patch, patch_adj
    type (STR_BoundaryCondition), pointer :: bc, bc_adj
    integer :: iface, ip
    integer :: ib_adj, if_adj, ip_adj, ideb_adj
        
    do iface = 1, 4
       face => b%face(iface)
       do ip = 1, face%nb_patch
          patch => face%patch(ip)
          bc => patch%bc
          
          if ( (iface==1 .or. iface==2) .and. &
               & (patch%ideb < b%md_tot .or. patch%ifin > b%mf_tot)) then
             print*, "Erreur sur les indices des patchs en m"
             call arret_code
          end if
          
          if ( (iface==3 .or. iface==4) .and. &
               & (patch%ideb < b%ld_tot .or. patch%ifin > b%lf_tot)) then
             print*, "Erreur sur les indices des patchs en l"
             call arret_code
          end if

          if (bc%condition == "interbloc") then
             ib_adj = patch%ib_adj
             if_adj = patch%if_adj
             ideb_adj = patch%ideb_adj
             
             if (ib_adj==-1) then
                print*, "Erreur condition limite interbloc"
                print*, "Pas de bloc adjacent"
                call arret_code
             end if

             if (bloc_est_local(ib_adj)) then
                b_adj => acces_bloc(ib_adj)%pointeur
                face_adj => b_adj%face(if_adj)
                
                do ip_adj = 1, face_adj%nb_patch
                   patch_adj => face_adj%patch(ip_adj)
                   if (ideb_adj >= patch_adj%ideb .and. &
                        & ideb_adj <= patch_adj%ifin) then
                      bc_adj => patch_adj%bc
                      
                      if ( bc_adj%condition /= "interbloc") then
                         print*, "Erreur condition limite interbloc"
                         print*, "Le bloc adjacent n'a pas de condition limite interbloc"
                         print*, 'b', b%numero, 'b_adj', b_adj%numero
                         call arret_code
                      end if
                   end if
                end do
             end if
          end if

          if (index(bc%condition, "itf") /=0) then
             ib_adj = patch%ib_adj
             if_adj = patch%if_adj
             ideb_adj = patch%ideb_adj
             
             if (ib_adj==-1) then
                print*, "Erreur condition limite interface"
                print*, "Pas de bloc adjacent"
                call arret_code
             end if

             b_adj => acces_bloc(ib_adj)%pointeur
             face_adj => b_adj%face(if_adj)

             if (b%nature == b_adj%nature .and. bc%condition/="itf_fluide" .and. bc%condition/="itf_fusion_sans_solide") then
                print*, "Erreur condition limite iterface"
                print*, "Les deux blocs en vis à vis sont de même nature"
                print*, 'b', b%numero, 'b_adj', b_adj%numero
                call arret_code
             end if
             
             do ip_adj = 1, face_adj%nb_patch
                patch_adj => face_adj%patch(ip_adj)
                if (ideb_adj >= patch_adj%ideb .and. &
                     & ideb_adj <= patch_adj%ifin) then
                   bc_adj => patch_adj%bc
                   
                   if ( index(bc_adj%condition, "itf") ==0 ) then
                      print*, "Erreur condition limite interface"
                      print*, "Le bloc adjacent n'a pas de condition limite interbloc"
                      print*, 'b', b%numero, 'b_adj', b_adj%numero
                      call arret_code
                   end if
                end if
             end do
          end if
          
       end do
    end do

  end subroutine verification_faces_bloc
     
end module m_verification
