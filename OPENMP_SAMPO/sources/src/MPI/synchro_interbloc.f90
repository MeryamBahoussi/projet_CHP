module synchronisation_interbloc
  use m_temps
  use m_struct
  use m_MPI
  use m_communication_MPI
  use parametres_fluide, only : iz

  implicit none

  integer, parameter :: sync_U = 1
  integer, parameter :: sync_uvT = 3
  integer, parameter :: sync_coord_nd = 4
  integer, parameter :: sync_T = 5 ! pour le solide (pas implemente)

contains
  
  subroutine synchro_interbloc(sb, send)
    
    type (STR_SUPER_BLOC), intent(inout) :: sb
    integer, intent(in) :: send ! 1 on synchronise U
    ! 2 on synchronise X_imp=(u,v,p) pour l'implicite
    ! 3 on synchronise (u,v,T) pour le visqueux
    ! 4 on synchronise le maillage

    integer :: i, ib, iface, ip, icoin
    type (STR_BLOC) , pointer :: b, b_adj
    type (STR_PATCH), pointer :: patch
    type (STR_COIN) , pointer :: coin

    call debut_watchTime(wT_commMPI)

    do i = 1, sb%nb_membre
       ib = sb%membre(i)
       b=>acces_bloc(ib)%pointeur

       if ( bloc_est_local (ib)) then
          
          ! communications interbloc (ceux definis dans les donnees)
          do iface = 1, 4
             do ip = 1, b%face(iface)%nb_patch
                patch => b%face(iface)%patch(ip)

                select case(trim(patch%bc%condition))
                case("interbloc")
                   b_adj => acces_bloc(patch%ib_adj)%pointeur
                   if (bloc_est_local (patch%ib_adj)) then
                      call comm_interbloc_local(iface, patch, b, b_adj,send)
                   else
                      call comm_interbloc_MPI(iface, patch, b, b_adj, send)
                      !call MPI_communication_interbloc(iface, patch, b, b_adj, send)
                   end if             
                end select
             end do
          end do
          ! communications periodiques
          do iface = 1, 4
             do ip = 1, b%face(iface)%nb_patch
                patch => b%face(iface)%patch(ip)

                select case(trim(patch%bc%condition))
                case("periodique")
                   if (send /= sync_coord_nd) then
                      b_adj => acces_bloc(patch%ib_adj)%pointeur
                      if (bloc_est_local (patch%ib_adj)) then
                         call comm_interbloc_local(iface, patch, b, b_adj,send)
                      else
                         call comm_interbloc_MPI(iface, patch, b, b_adj, send)
                      end if
                   end if
                end select
             end do
          end do
          ! communication des coins
          do icoin = 1, 4
             coin => b%coin(icoin)
             
             if (coin%ib_opp /= -1) then
                b_adj => acces_bloc(coin%ib_opp)%pointeur
                if (bloc_est_local (coin%ib_opp)) then
                   !call comm_interbloc_local(iface, patch, b, b_adj,send)
                   !call arret_code
                else
                   call comm_coin_MPI(coin, b, b_adj, send)
                end if
             end if
          end do
       end if
    end do

    call fin_watchTime(wT_commMPI)

  end subroutine synchro_interbloc

  subroutine comm_interbloc_local(iface, patch, b, b_adj,send) 
    integer, intent(in) :: iface
    type (STR_PATCH), intent(in) :: patch
    type (STR_BLOC), intent(inout) ::  b
    type (STR_BLOC), intent(in) ::  b_adj
    integer, intent(in) :: send

    integer :: ld, lf, md, mf
    integer :: ld_adj, md_adj
    integer :: l, m, ll, mm

    md = patch%ideb; mf = patch%ifin

    select case(iface)
    case(1)
       ld = b%ld-nMF; lf = b%ld-1 ! indices des mailles fictives à remplir
       md = patch%ideb; mf = patch%ifin
    case(2)
       ld = b%lf+1; lf = b%lf+nMF
       md = patch%ideb; mf = patch%ifin
    case(3)
       ld = patch%ideb; lf = patch%ifin
       md = b%md-nMF; mf=b%md-1
    case(4)
       ld = patch%ideb; lf = patch%ifin
       md = b%mf+1; mf=b%mf+nMF
    case default
       print*, "Erreur comm_interbloc_local"
       call arret_code
    end select

    select case(patch%if_adj)
    case(1) ! west : (ld,m)
       ld_adj = b_adj%ld 
       md_adj = patch%ideb_adj
    case(2) ! east : (lf,m)
       ld_adj = b_adj%lf-1
       md_adj = patch%ideb_adj
    case(3) ! south : (l,md)
       ld_adj = patch%ideb_adj
       md_adj = b_adj%md
    case(4) ! north : (l,mf)
       ld_adj = patch%ideb_adj
       md_adj = b_adj%mf-1
    end select

    select case(send)
    case(sync_U)
       do m = md, mf
          do l = ld, lf
             ll = ld_adj + (l-ld)
             mm = md_adj + (m-md)
             b%U(:,l,m) = b_adj%U(:, ll, mm)
          end do
       end do
    case(sync_uvT)
       do m = md, mf
          do l = ld, lf
             ll = ld_adj + (l-ld)
             mm = md_adj + (m-md)
             b%u_x(l,m) = b_adj%u_x(ll, mm)
             b%u_y(l,m) = b_adj%u_y(ll, mm)
             b%T(l,m) = b_adj%T(ll, mm)
          end do
       end do
    case(sync_coord_nd)
       print*, "communication en local pour le maillage non implementee"
       call arret_code
    end select

  end subroutine comm_interbloc_local

!!$  subroutine com_interbloc_local(iface, patch, b, b_adj,send) 
!!$    integer         , intent(in) :: iface
!!$    type (STR_BLOC) , pointer    :: b, b_adj
!!$    type (STR_PATCH), pointer    :: patch
!!$    integer         , intent(in) :: send
!!$    ! 1 on synchronise U
!!$    ! 2 on synchronise X_imp=(u,v,p) pour l'implicite
!!$    ! 3 on synchronise (u,v,T) pour le visqueux
!!$
!!$    integer :: bn(4), bn_adj(4)
!!$    
!!$    call indices_a_remplir(b, patch%ideb, patch%ifin, iface, indices)
!!$    call indices_a_envoyer(b_adj, patch%ideb_adj, patch%ideb_adj+patch%ifin-patch%ideb,iface,  indices_adj)
!!$
!!$    tab(bn(1):bn(2),bn(3):bn(4)) = tab_adj(bn_adj(1):bn_adj(2),bn_adj(3):bn_adj(4))
!!$    
!!$  end subroutine com_interbloc_local

  subroutine indices_a_envoyer(b, iface, ideb, ifin, indices)
    type (STR_BLOC)      , pointer       :: b
    integer              , intent(in)    :: iface, ideb, ifin
    integer, dimension(4), intent(inout) :: indices

    indices = 0

    select case(iface)
    case(1)
       indices(1:4) = (/b%ld  , b%ld+(nMF-1), ideb-nMF, ifin+nMF/)
    case(2)
       indices(1:4) = (/b%lf-(nMF-1), b%lf  , ideb-nMF, ifin+nMF/)
    case(3)
       indices(1:4) = (/ideb-nMF, ifin+nMF, b%md  , b%md+(nMF-1)/)
    case(4)
       indices(1:4) = (/ideb-nMF, ifin+nMF, b%mf-(nMF-1), b%mf  /)
    end select

  end subroutine indices_a_envoyer

  subroutine indices_a_remplir(b, iface, ideb, ifin, indices)
    type (STR_BLOC)      , pointer       :: b
    integer              , intent(in)    :: iface, ideb, ifin
    integer, dimension(4), intent(inout) :: indices

    indices = 0

    select case(iface)
    case(1)
       indices(1:4) = (/b%ld-nMF, b%ld-1, ideb-nMF, ifin+nMF/)
    case(2)
       indices(1:4) = (/b%lf+1, b%lf+nMF, ideb-nMF, ifin+nMF/)
    case(3)
       indices(1:4) = (/ideb-nMF, ifin+nMF, b%md-nMF, b%md-1/)
    case(4)
       indices(1:4) = (/ideb-nMF, ifin+nMF, b%mf+1, b%mf+nMF/)
    end select

  end subroutine indices_a_remplir
  
  subroutine indices_a_remplir_et_echanger(b, patch, iface, send, indices)

    type (STR_BLOC)      , pointer       :: b
    type (STR_PATCH)     , pointer       :: patch
    integer              , intent(in)    :: iface
    integer              , intent(in)    :: send
    integer, dimension(8), intent(inout) :: indices
    

    ! indices(1:4) : indices des mailles interieures a envoyer
    ! indices(5:8) : indices des mailles fictives a remplir

    indices(:) = 0
    
    if (send/=sync_coord_nd) then
       select case(iface)
       case(1)
          indices(5:8) = (/b%ld-nMF, b%ld-1      , patch%ideb-nMF, patch%ifin+nMF/)
          indices(1:4) = (/b%ld    , b%ld+(nMF-1), patch%ideb-nMF, patch%ifin+nMF/)
       case(2)
          indices(5:8) = (/b%lf+1      , b%lf+nMF, patch%ideb-nMF, patch%ifin+nMF/)
          indices(1:4) = (/b%lf-(nMF-1), b%lf    , patch%ideb-nMF, patch%ifin+nMF/)
       case(3)
          indices(5:8) = (/patch%ideb-nMF, patch%ifin+nMF, b%md-nMF, b%md-1/)
          indices(1:4) = (/patch%ideb-nMF, patch%ifin+nMF, b%md    , b%md+(nMF-1)/)
       case(4)
          indices(5:8) = (/patch%ideb-nMF, patch%ifin+nMF, b%mf+1      , b%mf+nMF/)
          indices(1:4) = (/patch%ideb-nMF, patch%ifin+nMF, b%mf-(nMF-1), b%mf  /)
       end select
    else
       select case(iface)
       case(1)
          indices(5:8) = (/b%ld-(nMF+1), b%ld-1  , patch%ideb-1-nmF, patch%ifin+nMF/)
          indices(1:4) = (/b%ld-1, b%ld+(nMF-1)  , patch%ideb-1-nmF, patch%ifin+nMF/)
       case(2)
          indices(5:8) = (/b%lf  , b%lf+nMF  , patch%ideb-1-nmF, patch%ifin+nMF/)
          indices(1:4) = (/b%lf-nMF, b%lf    , patch%ideb-1-nmF, patch%ifin+nMF/)
       case(3)
          indices(5:8) = (/patch%ideb-1-nmF, patch%ifin+nMF,   b%md-(nMF+1), b%md-1/)
          indices(1:4) = (/patch%ideb-1-nmF, patch%ifin+nMF,   b%md-1, b%md+(nMF-1)/)
       case(4)
          indices(5:8) = (/patch%ideb-1-nmF, patch%ifin+nMF,   b%mf  , b%mf+nMF/)
          indices(1:4) = (/patch%ideb-1-nmF, patch%ifin+nMF,   b%mf-nMF, b%mf  /)
       end select
    end if
    
  end subroutine indices_a_remplir_et_echanger

  subroutine comm_interbloc_MPI(iface, patch, b, b_adj, send) 
    
    integer         , intent(in) :: iface
    type (STR_BLOC) , pointer    :: b, b_adj
    type (STR_PATCH), pointer    :: patch
    integer         , intent(in) :: send
   
    integer :: indices(8)

    call indices_a_remplir_et_echanger(b, patch, iface, send, indices)
    
    select case(send)
    case(sync_U)
       call MPI_communication_vecteur(indices, b_adj%numproc, b%U)
    case(sync_uvT)
       call MPI_communication_scalaire(indices, b_adj%numproc, b%u_x)
       call MPI_communication_scalaire(indices, b_adj%numproc, b%u_y)
       call MPI_communication_scalaire(indices, b_adj%numproc, b%T)
    case(sync_T)
       call MPI_communication_scalaire(indices, b_adj%numproc, b%T)
    case(sync_coord_nd)
       call MPI_communication_vecteur(indices, b_adj%numproc, b%coord_nd)
    end select
  end subroutine comm_interbloc_MPI

  subroutine comm_coin_MPI(coin, b, b_adj, send) 
    type (STR_COIN), pointer    :: coin
    type (STR_BLOC), pointer    :: b, b_adj
    integer        , intent(in) :: send
    
    select case(send)
    case(sync_U)
       call MPI_communication_vecteur(coin%indices, b_adj%numproc, b%U) 
    case(sync_uvT)
       call MPI_communication_scalaire(coin%indices, b_adj%numproc, b%u_x)
       call MPI_communication_scalaire(coin%indices, b_adj%numproc, b%u_y)
       call MPI_communication_scalaire(coin%indices, b_adj%numproc, b%T)
    case(sync_T)
       call MPI_communication_scalaire(coin%indices, b_adj%numproc, b%T)
    case(sync_coord_nd)
     !  print*, "communication des coins pour le maillage non implementee"
     !  call arret_code
    end select
  end subroutine comm_coin_MPI

end module synchronisation_interbloc
