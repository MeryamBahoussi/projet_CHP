module m_initialisation
  use m_struct
  use m_MPI
  use m_initialisation_fluide
  use parametres_globaux, only : FLUIDE
  use m_gestion_interbloc

implicit none
  
contains
  
  subroutine initialisation_calcul(tache)
    implicit none
    
    type (STR_TACHE), intent(in) ::  tache
    
    integer                        :: i, j, isb, ib
    type (STR_SUPER_BLOC), pointer :: sb
    type (STR_BLOC)      , pointer :: b

!!! initialisation des infos des coins (utile pour les synchro_interbloc)
    do i=1, size(acces_bloc, 1)
       b=> acces_bloc(i)%pointeur
       call initialisation_coins(b)
    end do
    
!!! Initialisation des blocs fluides et solides
!!! pour le fluide, on decode aussi U et on stocke Un et pn
    do i = 1, tache%nb_membre 
       isb = tache%membre(i)
       sb => acces_super_bloc(isb)%pointeur 
       
       do j = 1, sb%nb_membre 
          ib = sb%membre(j)
          if (bloc_est_local (ib)) then
             b=>acces_bloc(ib)%pointeur
             if (b%nature == FLUIDE) then
                call initialisation_fluide(b, b%maillage%type, b%type_initialisation)
             end if
          end if
       end do
    end do
    
  end subroutine initialisation_calcul
    
end module m_initialisation
