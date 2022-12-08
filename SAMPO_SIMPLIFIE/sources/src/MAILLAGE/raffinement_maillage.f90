module m_raffinement_maillage

  use m_struct
  use m_MPI
  use outils_maillage
  use parametres_globaux, only : nMF

  implicit none

contains

  subroutine resserrement_de_maillage( b, coord_nd )
    implicit none 

    type (STR_BLOC), pointer ::  b
    real*8, dimension(1:2,b%ld-(nMF+1):b%lf+nMF,b%md-(nMF+1):b%mf+nMF), intent(inout) :: coord_nd

    integer :: face

    face = b%maillage%face

    select case( face )
    case(1) ! resserrement sur le bord ouest
       call raffine_suivant_l( b, b%ld_tot-1, coord_nd )
    case(2) ! resserrement sur le bord est
       call raffine_suivant_l( b, b%lf_tot, coord_nd )
    case(3) ! resserrement sur le bord sud
       call raffine_suivant_m( b, b%md_tot-1, coord_nd )
    case(4) ! resserrement sur le bord nord
       call raffine_suivant_m( b, b%mf_tot, coord_nd )
    case(12)
       call raffine_double_suivant_l(b, coord_nd)
    case(34) ! resserrement sur les bords sud et nord
       call raffine_double_suivant_m(b, coord_nd)  
    case(43) ! a enlever
       call raffine_suivant_m( b, 10, coord_nd )
    case(13) ! resserrement sur les bords ouest et sud
       call raffine_suivant_l( b, b%ld_tot-1, coord_nd )
       call raffine_suivant_m( b, b%md_tot-1, coord_nd )
    case(23) ! resserrement sur les bords est et sud
       call raffine_suivant_l( b, b%lf_tot  , coord_nd )
       call raffine_suivant_m( b, b%md_tot-1, coord_nd )
    case(123) ! resserrement sur les bords est, ouest et sud
       call raffine_double_suivant_l(b, coord_nd )
       call raffine_suivant_m(b, b%md_tot-1, coord_nd)
    case(134) ! resserrement sur les bords ouest, sud et nord
       call raffine_suivant_l( b, b%ld_tot-1, coord_nd )
       call raffine_double_suivant_m(b, coord_nd)
    case(234) ! resserrement sur les bords est, sud et nord
       call raffine_suivant_l( b, b%lf_tot, coord_nd ) 
       call raffine_double_suivant_m(b, coord_nd)
    case(1234) ! resserrement sur tous les bords
       call raffine_double_suivant_l(b, coord_nd)
       call raffine_double_suivant_m(b, coord_nd)
       b%maillage%face = 34
    case default
       print*, "type de raffinement maillage non defini"
       call arret_code
    end select

  end subroutine resserrement_de_maillage
  

!!$  subroutine raffine_suivant_l( b, l0, coord_nd )
!!$    implicit none
!!$
!!$    type (STR_BLOC), pointer ::  b
!!$    integer :: l0
!!$    real*8, dimension(1:2, b%ld-3:b%lf+2, b%md-3:b%mf+2), intent(inout) :: coord_nd
!!$
!!$    integer :: l, m, l_pere
!!$    integer :: me_interf, statinfo, me_interf_MPI
!!$    real*8 :: beta
!!$    real*8 :: distance_max, vecteur(2), dx
!!$    real*8, dimension(:,:,:), allocatable :: new_nodes
!!$
!!$    real*8 :: dmax_g, dmax_d, vecteur_g(2), vecteur_d(2)
!!$    real*8 :: cd(2), cI(2), cf(2)
!!$
!!$    ! On verifie que la position est bien dans le domaine de calcul
!!$    if (l0<b%ld_tot-1 .or. l0>b%lf_tot) then
!!$       if ( numproc == 0 ) then 
!!$          print*," L'indice de resserrement de maillage n'est pas coherent"
!!$          print*," Raffinement autour de l=", l0
!!$          print*,"ERREUR : raffine_suivant_l"
!!$       end if
!!$       call arret_code  
!!$    end if
!!$
!!$    allocate( new_nodes(1:2, b%ld-3:b%lf+2, b%md-3:b%mf+2) ) 
!!$    new_nodes = coord_nd
!!$
!!$    beta = b%maillage%beta_l
!!$
!!$    ! On determine le proc qui connait coord_nd(:,l0, :)  
!!$    me_interf =-1
!!$    if (l0 >= b%l_charge(1,numproc)-1 .and. l0 <= b%l_charge(2,numproc)) me_interf=numproc
!!$    call MPI_ALLREDUCE(me_interf, me_interf_MPI, 1, MPI_INT, MPI_MAX, sb%mpi_comm,STATINFO)
!!$    me_interf = me_interf_MPI
!!$    
!!$    do m = b%md_tot-3, b%mf_tot+2
!!$       cd = 0.d0; cI = 0.d0; cf = 0.d0
!!$
!!$       ! tous les procs doivent connaitre :
!!$       ! cd = coord_nd(:,b%ld_tot-1,m), 
!!$       ! cf = coord_nd(:,b%lf_tot, m) et 
!!$       ! cI = coord_nd(:,l0, m)
!!$       if (numproc==b%proc_ini) cd = coord_nd(:,b%ld-1,m)
!!$       if (numproc==b%proc_fin) cf = coord_nd(:,b%lf, m)
!!$       if (numproc==me_interf) cI = coord_nd(:,l0-b%l_charge(1,numproc)+1, m)
!!$
!!$       if (b%proc_ini/=b%proc_fin) then
!!$          CALL MPI_BCAST(cd,2,MPI_REAL8,b%proc_ini,sb%mpi_comm,STATINFO)
!!$          CALL MPI_BCAST(cf,2,MPI_REAL8,b%proc_fin,sb%mpi_comm,STATINFO)
!!$          CALL MPI_BCAST(cI,2,MPI_REAL8,me_interf,sb%mpi_comm,STATINFO)
!!$       end if
!!$
!!$       dmax_g = 0.d0; dmax_d = 0.d0
!!$       vecteur_g = 0.d0; vecteur_d = 0.d0
!!$       if (l0 /= b%ld_tot-1) then
!!$          dmax_g = distance( cI,cd )
!!$          vecteur_g = ( cI-cd ) / dmax_g
!!$       end if
!!$       if (l0 /= b%lf_tot) then
!!$          dmax_d = distance( cf,cI )
!!$          vecteur_d = ( cf-cI ) / dmax_d
!!$       end if
!!$
!!$       do l = b%ld-1, b%lf
!!$          ! on se replace dans la numerotation initiale en l
!!$          l_pere = l + b%l_charge(1, numproc)-1
!!$          if (l_pere<l0) then
!!$             dx = float(l0 - l_pere) / float( l0 - b%ld_tot + 1 )
!!$             new_nodes(:,l,m) = cI - dmax_g * &
!!$                  transformation_de_resserrement_de_maillage(dx, beta) * vecteur_g(:)
!!$          else if (l_pere>l0) then
!!$             dx = float(l_pere - l0) / float( b%lf_tot - l0 )
!!$             new_nodes(:,l,m) = cI + dmax_d * &
!!$                  transformation_de_resserrement_de_maillage(dx, beta) * vecteur_d(:)
!!$          else if (l_pere==l0) then
!!$             new_nodes(:,l,m) = cI
!!$          end if
!!$       end do
!!$
!!$       ! cellules fictives
!!$       distance_max = distance( coord_nd(:,b%lf,m), coord_nd(:,b%ld-1,m) )
!!$       vecteur = ( coord_nd(:,b%lf,m) - coord_nd(:,b%ld-1,m) ) / distance_max
!!$       dx = distance( new_nodes(:,b%ld,m), new_nodes(:,b%ld-1,m) )
!!$       do l = b%ld-3, b%ld-2
!!$          new_nodes(:,l,m) = coord_nd(:,b%ld-1,m) + dx * ( l - (b%ld-1) ) * vecteur(:)
!!$       end do
!!$       dx = distance( new_nodes(:,b%lf,m), new_nodes(:,b%lf-1,m) )
!!$       do l = b%lf+1, b%lf+2
!!$          new_nodes(:,l,m) = coord_nd(:,b%lf,m) + dx * ( l - (b%lf) ) * vecteur(:)
!!$       end do
!!$
!!$    end do
!!$
!!$    coord_nd = new_nodes
!!$
!!$    deallocate(new_nodes)
!!$
!!$  end subroutine raffine_suivant_l

  subroutine raffine_suivant_l( b, l0, coord_nd )
    implicit none

    type (STR_BLOC), pointer ::  b
    integer :: l0
    real*8, dimension(1:2, b%ld-(nMF+1):b%lf+nMF,b%md-(nMF+1):b%mf+nMF), intent(inout) :: coord_nd

    integer :: l, m, l_pere
    integer :: ld, lf, ld_mfg, lf_mfg, ld_mfd, lf_mfd
    real*8 :: beta
    real*8 :: distance_max, vecteur(2), dx
    real*8, dimension(:,:,:), allocatable :: new_nodes
    real*8, dimension(:,:)  , allocatable :: cd_, cf_
    real*8, dimension(2)                  :: cd(2), cf(2)
    
    ! On verifie que la position est bien dans le domaine de calcul
    if (l0<b%ld_tot-1 .or. l0>b%lf_tot) then
       if ( numproc == 0 ) then 
          print*," L'indice de resserrement de maillage n'est pas coherent"
          print*," Raffinement autour de l=", l0
          print*," ERREUR : raffine_suivant_l"
       end if
       call arret_code  
    end if

    allocate( new_nodes(1:2,b%ld-(nMF+1):b%lf+nMF,b%md-(nMF+1):b%mf+nMF))
    new_nodes = 0.d0

    beta = b%maillage%beta_l

    allocate(cd_(2,b%md-(nMF+1):b%mf+nMF))
    allocate(cf_(2,b%md-(nMF+1):b%mf+nMF))
    call find_noeud_bord(b, 1, coord_nd, b%md, b%mf, cd_(:,:))
    call find_noeud_bord(b, 2, coord_nd, b%md, b%mf, cf_(:,:))
    
    do m = b%md-(nMF+1), b%mf+nMF
       cd = cd_(:,m); cf = cf_(:,m)

       distance_max = distance( cf, cd )
       vecteur = ( cf-cd ) / distance_max

       ! on ne traite que les mailles fictives du bord exterieur du domaine separement
       if (nbprocs == 1) then
          ld = b%ld-1; lf = b%lf   ! indices des points consideres comme interieur
          ld_mfg = b%ld-(nMF+1); lf_mfg = b%ld-2    ! mailles fictives a gauche (en l=ld)
          ld_mfd = b%lf+1; lf_mfd = b%lf+nMF    ! mailles fictives a droite (en l=lf)
       else
          if (b%l_charge(1,numproc) == b%ld_tot) then
             ld = b%ld-1
             ld_mfg = b%ld-(nMF+1); lf_mfg = b%ld-2
          else 
             ld = b%ld-(nMF-1)
             ld_mfg = 0 ; lf_mfg = -1         ! on ne fait pas les MF a gauche
          end if
          
          if (b%l_charge(2,numproc) == b%lf_tot) then
             lf = b%lf
             ld_mfd = b%lf+1; lf_mfd = b%lf+nMF
          else
             lf = b%lf+nMF
             ld_mfd = 0 ; lf_mfd = -1
          end if
       end if
       
       ! points interieurs
       if (l0 == b%ld_tot-1) then
          do l = ld, lf
             ! on se replace dans la numerotation initiale en l
             l_pere = l + b%l_charge(1, numproc)-1
             dx = float(l_pere) / float( b%lf_tot - b%ld_tot + 1 )
             new_nodes(:,l,m) = cd + distance_max * transformation_de_resserrement_de_maillage(dx, beta) * vecteur(:)
          end do
       else if (l0 == b%lf_tot) then
          do l = ld, lf
             ! on se replace dans la numerotation initiale en l
             l_pere = l + b%l_charge(1, numproc)-1
             dx = float(l_pere) / float( b%lf_tot - b%ld_tot + 1 )
             new_nodes(:,lf-l+ld,m) = cf - distance_max * transformation_de_resserrement_de_maillage(dx, beta) * vecteur(:)
          end do
       end if
       
       ! mailles fictives des bords
       dx = distance( new_nodes(:,b%ld,m), new_nodes(:,b%ld-1,m) )
       do l = ld_mfg, lf_mfg
          new_nodes(:,l,m) = coord_nd(:,b%ld-1,m) + dx*(l - (b%ld-1)) * vecteur(:)
       end do
       dx = distance( new_nodes(:,b%lf,m), new_nodes(:,b%lf-1,m) )
       do l = ld_mfd, lf_mfd
          new_nodes(:,l,m) = coord_nd(:,b%lf,m) + dx*(l - (b%lf)) * vecteur(:)
       end do
       
    end do

    coord_nd = new_nodes
    deallocate(new_nodes, cd_, cf_)

  end subroutine raffine_suivant_l

  subroutine raffine_double_suivant_l(b, coord_nd)
    implicit none

    type (STR_BLOC), pointer ::  b
    real*8, dimension(1:2, b%ld-(nMF+1):b%lf+nMF,b%md-(nMF+1):b%mf+nMF), intent(inout) :: coord_nd

    integer :: l, m, l_pere
    integer :: ld, lf, ld_mfg, lf_mfg, ld_mfd, lf_mfd
    real*8 :: beta
    real*8 :: distance_max, vecteur(2), dx
    real*8, dimension(:,:,:), allocatable :: new_nodes
    real*8, dimension(:,:)  , allocatable :: cd_, cf_
    real*8, dimension(2)                  :: cd(2), cf(2) 
    
    beta=b%maillage%beta_l

    allocate( new_nodes(1:2, b%ld-(nMF+1):b%lf+nMF,b%md-(nMF+1):b%mf+nMF) ) 
    new_nodes = 0.d0

    allocate(cd_(2,b%md-(nMF+1):b%mf+nMF))
    allocate(cf_(2,b%md-(nMF+1):b%mf+nMF))
    
    ! tous les procs doivent connaitre :
    ! cd = coord_nd(:,b%ld_tot-1,m) et 
    ! cf = coord_nd(:,b%lf_tot, m) 
   ! if (numproc==6) then
    call find_noeud_bord(b, 1, coord_nd, b%md, b%mf, cd_(:,:))
    call find_noeud_bord(b, 2, coord_nd, b%md, b%mf, cf_(:,:))
    !else
    !   cd_ = b%coord_nd(:,b%ld-1,b%md-3:b%mf+2)
    !   cf_ = b%coord_nd(:,b%lf, b%md-3:b%mf+2)
    !end if
    !beta = 1.1d0

    do m = b%md-(nMF+1), b%mf+nMF
       cd = cd_(:,m); cf = cf_(:,m)
!!$       cd = 0.d0; cf = 0.d0     
!!$       if (numproc==b%proc_ini) cd = coord_nd(:,b%ld-1,m)
!!$       if (numproc==b%proc_fin) cf = coord_nd(:,b%lf, m)
!!$
!!$       if (b%proc_ini/=b%proc_fin) then
!!$          CALL MPI_BCAST(cd,2,MPI_REAL8,b%proc_ini,sb%mpi_comm,STATINFO)
!!$          CALL MPI_BCAST(cf,2,MPI_REAL8,b%proc_fin,sb%mpi_comm,STATINFO)
!!$       end if
       
       distance_max = distance( cf, cd )
       vecteur = ( cf-cd ) / distance_max

       ! on ne traite que les mailles fictives du bord exterieur du domaine separement
       if (nbprocs == 1) then
          ld = b%ld-1; lf = b%lf   ! indices des points consideres comme interieur
          ld_mfg = b%ld-(nMF+1); lf_mfg = b%ld-2    ! mailles fictives a gauche (en l=ld)
          ld_mfd = b%lf+1; lf_mfd = b%lf+nMF    ! mailles fictives a droite (en l=lf)
       else
          if (b%l_charge(1,numproc) == b%ld_tot) then
             ld = b%ld-1
             ld_mfg = b%ld-(nMF+1); lf_mfg = b%ld-2
          else 
             ld = b%ld-(nMF+1)
             ld_mfg = 0 ; lf_mfg = -1         ! on ne fait pas les MF a gauche
          end if
          
          if (b%l_charge(2,numproc) == b%lf_tot) then
             lf = b%lf
             ld_mfd = b%lf+1; lf_mfd = b%lf+nMF
          else
             lf = b%lf+nMF
             ld_mfd = 0 ; lf_mfd = -1
          end if
       end if

       ! points interieurs
       do l = ld, lf
          ! on se replace dans la numerotation initiale en l
          l_pere = l + b%l_charge(1, numproc)-1
          dx = float(l_pere) / float( b%lf_tot - b%ld_tot + 1 )
          new_nodes(:,l,m) = cd + distance_max * transformation_de_bi_resserrement_de_maillage(dx, beta) * vecteur(:)
       end do

       ! mailles fictives des bords
       dx = distance( new_nodes(:,b%ld,m), new_nodes(:,b%ld-1,m) )
       do l = ld_mfg, lf_mfg
          new_nodes(:,l,m) = coord_nd(:,b%ld-1,m) + dx*(l - (b%ld-1)) * vecteur(:)
       end do
       dx = distance( new_nodes(:,b%lf,m), new_nodes(:,b%lf-1,m) )
       do l = ld_mfd, lf_mfd
          new_nodes(:,l,m) = coord_nd(:,b%lf,m) + dx*(l - (b%lf)) * vecteur(:)
       end do

    end do

    coord_nd = new_nodes
    deallocate(new_nodes, cd_, cf_)

  end subroutine raffine_double_suivant_l


  subroutine raffine_suivant_m( b, m0, coord_nd )
    implicit none

    type (STR_BLOC), pointer ::  b
    integer :: m0
    real*8, dimension(1:2, b%ld-(nMF+1):b%lf+nMF,b%md-(nMF+1):b%mf+nMF), intent(inout) :: coord_nd

    integer :: l, m, m_pere
    integer :: med, meI, mef
    real*8 :: beta
    real*8 :: distance_max, vecteur(2), dx
    real*8 :: dmax_g, dmax_d, vecteur_g(2), vecteur_d(2)
    real*8, dimension(:,:,:), allocatable :: new_nodes
    real*8, dimension(:,:)  , allocatable :: cd_, cI_, cf_
    real*8, dimension(2)                  :: cd(2), cI(2), cf(2) 

    ! On verifie que la position est bien dans le domaine de calcul
    if (m0<b%md_tot-1 .or. m0>b%mf_tot) then
       if ( numproc == 0 ) then 
          print*," L'indice de resserrement de maillage n'est pas coherent"
          print*," Raffinement autour de m=", m0
          print*," ERREUR : raffine_suivant_m"
       end if
       call arret_code  
    end if

    allocate( new_nodes(1:2, b%ld-(nMF+1):b%lf+nMF,b%md-(nMF+1):b%mf+nMF) )
    new_nodes = coord_nd

    beta = b%maillage%beta_m

    allocate(cd_(2,b%ld-(nMF+1):b%lf+nMF))
    allocate(cI_(2,b%ld-(nMF+1):b%lf+nMF))
    allocate(cf_(2,b%ld-(nMF+1):b%lf+nMF))
    call recuperation_info_resserrement(b, m0, b%ld, b%lf, cd_(:,:), cI_(:,:), cf_(:,:))

    do l = b%ld-(nMF+1), b%lf+nMF
       cd = cd_(:,l); cI = cI_(:,l); cf = cf_(:,l)
       med = -1 ; meI = -1 ; mef = -1

       ! tous les procs doivent connaitre :
       ! cd = coord_nd(:,l,b%md_tot-1), 
       ! cf = coord_nd(:,l,b%mf_tot) et 
       ! cI = coord_nd(:,l, m0)

       dmax_g = 0.d0; dmax_d = 0.d0
       vecteur_g = 0.d0; vecteur_d = 0.d0
       if (m0 /= b%md_tot-1) then
          dmax_g = distance( cI,cd )
          vecteur_g = ( cI-cd ) / dmax_g
       end if
       if (m0 /= b%mf_tot) then
          dmax_d = distance( cf,cI )
          vecteur_d = ( cf-cI ) / dmax_d
       end if

       do m = b%md-1, b%mf
          ! on se replace dans la numerotation initiale en m
          m_pere = m + b%m_charge(1, numproc)-1
          if (m_pere<m0) then
             dx = float(m0 - m_pere) / float( m0 - b%md_tot + 1 )
             new_nodes(:,l,m) = cI - dmax_g * &
                  transformation_de_resserrement_de_maillage(dx, beta) * vecteur_g(:)
          else if (m_pere>m0) then
             dx = float(m_pere - m0) / float( b%mf_tot - m0 )
             new_nodes(:,l,m) = cI + dmax_d * &
                  transformation_de_resserrement_de_maillage(dx, beta) * vecteur_d(:)
          else if (m_pere==m0) then
             new_nodes(:,l,m) = cI
          end if
       end do

       ! cellules fictives
       distance_max = distance( coord_nd(:,l,b%mf), coord_nd(:,l,b%md-1) )
       vecteur = ( coord_nd(:,l,b%mf) - coord_nd(:,l,b%md-1) ) / distance_max
       dx = distance( new_nodes(:,l,b%md), new_nodes(:,l,b%md-1) )
       do m = b%md-(nMF+1), b%md-2
          new_nodes(:,l,m) = coord_nd(:,l,b%md-1) + dx*(m - (b%md-1)) * vecteur(:)
       end do
       dx = distance( new_nodes(:,l,b%mf), new_nodes(:,l,b%mf-1) )
       do m = b%mf+1, b%mf+nMF
          new_nodes(:,l,m) = coord_nd(:,l,b%mf) + dx*(m - b%mf) * vecteur(:)
       end do

    end do

    coord_nd = new_nodes

    deallocate(new_nodes, cd_, cI_, cf_)

  end subroutine raffine_suivant_m

  subroutine raffine_suivant_m_old( b, m0, coord_nd )
    implicit none

    type (STR_BLOC), pointer ::  b
    integer :: m0
    real*8, dimension(1:2, b%ld-(nMF+1):b%lf+nMF,b%md-(nMF+1):b%mf+nMF), intent(inout) :: coord_nd

    integer :: l, m
    real*8 :: beta
    real*8 :: distance_max, vecteur(2), dx
    real*8, dimension(:,:,:), allocatable :: new_nodes

    real*8 :: dmax_g, dmax_d, vecteur_g(2), vecteur_d(2)
    real*8 :: cd(2), cI(2), cf(2)

    ! On verifie que la position est bien dans le domaine de calcul
    if (m0<b%md_tot-1 .or. m0>b%mf_tot) then
       if ( numproc == 0 ) then 
          print*," L'indice de resserrement de maillage n'est pas coherent"
          print*," Raffinement autour de m=", m0
          print*," ERREUR : raffine_suivant_m"
       end if
       call arret_code  
    end if

    allocate( new_nodes(1:2, b%ld-(nMF+1):b%lf+nMF,b%md-(nMF+1):b%mf+nMF) ) 
    new_nodes = coord_nd

    beta = b%maillage%beta_m

    do l = b%ld-(nMF+1), b%lf+nMF
       cd = coord_nd(:,l,b%md_tot-1)
       cf = coord_nd(:,l,b%mf_tot)
       cI = coord_nd(:,l,m0)

       dmax_g = 0.d0; dmax_d = 0.d0
       vecteur_g = 0.d0; vecteur_d = 0.d0
       if (m0 /= b%md_tot-1) then
          dmax_g = distance( cI,cd )
          vecteur_g = ( cI-cd ) / dmax_g
       end if
       if (m0 /= b%mf_tot) then
          dmax_d = distance( cf,cI )
          vecteur_d = ( cf-cI ) / dmax_d
       end if

       do m = b%md-1, b%mf
          if (m<m0) then
             dx = float(m0 - m) / float( m0 - b%md_tot + 1 )
             new_nodes(:,l,m) = cI - dmax_g * &
                  transformation_de_resserrement_de_maillage(dx, beta) * vecteur_g(:)
          else if (m>m0) then
             dx = float(m - m0) / float( b%mf_tot - m0 )
             new_nodes(:,l,m) = cI + dmax_d * &
                  transformation_de_resserrement_de_maillage(dx, beta) * vecteur_d(:)
          else if (m==m0) then
             new_nodes(:,l,m) = cI
          end if
       end do

       ! cellules fictives
       distance_max = distance( cf, cd )
       vecteur = ( cf-cd ) / distance_max

       dx = distance( new_nodes(:,l,b%md), new_nodes(:,l,b%md-1) )
       do m = b%md-(nMF+1), b%md-2
          new_nodes(:,l,m) = coord_nd(:,l,b%md-1) + dx*(m - (b%md-1)) * vecteur(:)
       end do

       dx = distance( new_nodes(:,l,b%mf), new_nodes(:,l,b%mf-1) )
       do m = b%mf+1, b%mf+nMF
          new_nodes(:,l,m) = coord_nd(:,l,b%mf) + dx*(m - (b%mf)) * vecteur(:)
       end do

    end do

    coord_nd = new_nodes
    deallocate(new_nodes)

  end subroutine raffine_suivant_m_old

  subroutine raffine_double_suivant_m_old(b, coord_nd)
    implicit none

    type (STR_BLOC), pointer ::  b
    real*8, dimension(1:2, b%ld-(nMF+1):b%lf+nMF,b%md-(nMF+1):b%mf+nMF), intent(inout) :: coord_nd

    integer :: l, m
    real*8 :: beta
    real*8 :: distance_max, vecteur(2), dx
    real*8, dimension(:,:,:), allocatable :: new_nodes

    beta=b%maillage%beta_m

    allocate( new_nodes(1:2, b%ld-(nMF+1):b%lf+nMF,b%md-(nMF+1):b%mf+nMF) ) 
    new_nodes = 0.d0

    do l = b%ld-(nMF+1), b%lf+nMF
       distance_max = distance( coord_nd(:,l,b%mf), coord_nd(:,l,b%md-1) )
       vecteur = ( coord_nd(:,l,b%mf) - coord_nd(:,l,b%md-1) ) / distance_max
       do m = b%md-1, b%mf
          dx = float(m) / float( b%mf - b%md + 1 )
          new_nodes(:,l,m) = coord_nd(:,l,b%md-1) + distance_max * transformation_de_bi_resserrement_de_maillage(dx, beta) * vecteur(:)
       end do

       dx = distance( new_nodes(:,l,b%md), new_nodes(:,l,b%md-1) )
       do m = b%md-(nMF+1), b%md-2
          new_nodes(:,l,m) = coord_nd(:,l,b%md-1) + dx*(m - (b%md-1)) * vecteur(:)
       end do

       dx = distance( new_nodes(:,l,b%mf), new_nodes(:,l,b%mf-1) )
       do m = b%mf+1, b%mf+nMF
          new_nodes(:,l,m) = coord_nd(:,l,b%mf) + dx*(m - (b%mf)) * vecteur(:)
       end do
    end do

    coord_nd = new_nodes
    deallocate(new_nodes)

  end subroutine raffine_double_suivant_m_old

  subroutine raffine_double_suivant_m(b, coord_nd)
    implicit none

    type (STR_BLOC), pointer ::  b
    real*8, dimension(1:2, b%ld-(nMF+1):b%lf+nMF,b%md-(nMF+1):b%mf+nMF), intent(inout) :: coord_nd

    integer :: l, m, m_pere
    integer :: md, mf, md_mfg, mf_mfg, md_mfd, mf_mfd
    real*8 :: beta
    real*8 :: dx, distance_max, vecteur(2)
    real*8, dimension(:,:,:), allocatable :: new_nodes
    real*8, dimension(:,:)  , allocatable :: cd_, cf_
    real*8, dimension(2)                  :: cd(2), cf(2) 

    beta = b%maillage%beta_m

    allocate( new_nodes(1:2, b%ld-(nMF+1):b%lf+nMF,b%md-(nMF+1):b%mf+nMF) ) 
    new_nodes = 0.d0
   
    allocate(cd_(2,b%ld-(nMF+1):b%lf+nMF))
    allocate(cf_(2,b%ld-(nMF+1):b%lf+nMF))
    
    call find_noeud_bord(b, 3, coord_nd, b%ld, b%lf, cd_(:,:))
    call find_noeud_bord(b, 4, coord_nd, b%ld, b%lf, cf_(:,:))

    do l = b%ld-(nMF+1), b%lf+nMF
       cd = cd_(:,l); cf = cf_(:,l)
       ! tous les procs doivent connaitre :
       ! cd = coord_nd(:,l,b%md_tot-1) et 
       ! cf = coord_nd(:,l,b%mf_tot)

       distance_max = distance( cf, cd )
       vecteur = ( cf-cd ) / distance_max
       
       ! on ne traite que les mailles fictives du bord exterieur du domaine separement
       if (nbprocs == 1) then
          md = b%md-1; mf = b%mf   ! indices des points consideres comme interieur
          md_mfg = b%md-(nMF+1); mf_mfg = b%md-2    ! mailles fictives a gauche (en l=ld)
          md_mfd = b%mf+1; mf_mfd = b%mf+nMF    ! mailles fictives a droite (en l=lf)
       else
          if (b%m_charge(1,numproc) == b%md_tot) then
             md = b%md-1
             md_mfg = b%md-(nMF+1); mf_mfg = b%md-2
          else 
             md = b%md-(nMF+1)
             md_mfg = 0 ; mf_mfg = -1         ! on ne fait pas les MF a gauche
          end if
          
          if (b%m_charge(2,numproc) == b%mf_tot) then
             mf = b%mf
             md_mfd = b%mf+1; mf_mfd = b%mf+nMF
          else
             mf = b%mf+nMF
             md_mfd = 0 ; mf_mfd = -1
          end if
       end if

       ! points interieurs
       do m = md, mf
          ! on se replace dans la numerotation initiale en m
          m_pere = m + b%m_charge(1, numproc)-1
          dx = float(m_pere) / float( b%mf_tot - b%md_tot + 1 )
          new_nodes(:,l,m) = cd + distance_max * transformation_de_bi_resserrement_de_maillage(dx, beta) * vecteur(:)
       end do
       
       ! mailles fictives des bords
       dx = distance( new_nodes(:,l,b%md), new_nodes(:,l,b%md-1) )
       do m = md_mfg, mf_mfg
          new_nodes(:,l,m) = coord_nd(:,l,b%md-1) + dx*(m - (b%md-1)) * vecteur(:)
       end do
       dx = distance( new_nodes(:,l,b%mf), new_nodes(:,l,b%mf-1) )
       do m = md_mfd, mf_mfd
          new_nodes(:,l,m) = coord_nd(:,l,b%mf) + dx*(m - b%mf) * vecteur(:)
       end do  

    end do

    coord_nd = new_nodes
    deallocate(new_nodes, cd_, cf_)

  end subroutine raffine_double_suivant_m

  function transformation_de_resserrement_de_maillage(x_regulier, beta)
    implicit none 
    real*8 :: x_regulier, beta
    real*8 :: transformation_de_resserrement_de_maillage

    real*8 ::  zlog, puiss, rapp

    ! La subroutine est inspirée d'un subroutine de Gerard Gallice:
    ! transformation de [0,1] sur [0,1] permettant de resserrer le maillage au voisinage du corps

    ! beta = coeff_resserrement_maillage_fld
    zlog = log( (beta+1) / (beta-1) )
    puiss = zlog * (1.d0 - x_regulier)
    puiss = exp(puiss)
    rapp = (1.d0 - puiss) / (1.d0 + puiss)

    transformation_de_resserrement_de_maillage = 1.d0 + beta * rapp

  end function transformation_de_resserrement_de_maillage

  function transformation_de_bi_resserrement_de_maillage(x_regulier, beta)
    implicit none 
    real*8 :: x_regulier, beta
    real*8 :: transformation_de_bi_resserrement_de_maillage

    transformation_de_bi_resserrement_de_maillage = transformation_de_resserrement_de_maillage(2.d0 * x_regulier, beta) / 2.d0

  end function transformation_de_bi_resserrement_de_maillage

  subroutine recuperation_info_resserrement(b, m0, ld, lf, cd, cI, cf)

    type (STR_BLOC)               , pointer     ::  b
    integer                       , intent(in)  :: m0, ld, lf
    real*8, dimension(2,ld-(nMF+1):lf+nMF), intent(out) :: cd, cI, cf

    integer :: med, meI, mef, j, mey
    integer :: taille, STATINFO, MPI_FLAG=75
    integer :: MPI_status(MPI_STATUS_SIZE)
    integer :: NbCoupes(2)

    NbCoupes = b%NbCoupes(b%num_input,:)

    ! Tout le monde connait le numero du processeur qui connait cd, cf
    med = b%proc_ini + mod(numproc-b%proc_ini,NbCoupes(1))
    mef = med + (NbCoupes(2)-1)*NbCoupes(1)

    taille = 2*(lf+nMF-(ld-(nMF+1))+1)

!!! points de la face 3
    ! le proc med envoie l'info a tous ceux au dessus de lui
    if (numproc==med) then
       cd(1:2,ld-(nMF+1):lf+nMF) = b%coord_nd(:,:,b%md-1)
       do j = 1, NbCoupes(2)
          call MPI_SEND( b%coord_nd(:,:,b%md-1), taille, MPI_REAL8, &
               med+(j-1)*NbCoupes(1), MPI_FLAG, MPI_COMM_WORLD, STATINFO )
       end do
    else
       ! reception des points de la face 3
       call MPI_RECV( cd(1:2,ld-(nMF+1):lf+nMF), taille, MPI_REAL8, &
            med, MPI_FLAG, MPI_COMM_WORLD, MPI_status, STATINFO )     
    end if

!!! points de la face 4
    ! le proc med envoie l'info a tous ceux en dessous de lui
    if (numproc==mef) then
       cf(1:2,ld-(nMF+1):lf+nMF) =  b%coord_nd(:,:,b%mf)
       do j = 1, NbCoupes(2)
          call MPI_SEND( b%coord_nd(:,:,b%mf), taille, MPI_REAL8, &
               med+(j-1)*NbCoupes(1), MPI_FLAG, MPI_COMM_WORLD, STATINFO )
       end do
    else
       ! reception des points de la face 4
       call MPI_RECV( cf(1:2,ld-(nMF+1):lf+nMF), taille, MPI_REAL8, &
            mef, MPI_FLAG, MPI_COMM_WORLD, MPI_status, STATINFO )     
    end if

    if (m0==b%md-1) then
       cI = cd
    else if (m0==b%mf) then
       cI = cf
    else
       do j = 1, NbCoupes(2)
          mey = med+(j-1)*NbCoupes(1)
          if (m0 >= b%m_charge(1,mey)-1 .and. m0 <= b%m_charge(2,mey)) then
             meI = mey
          end if
       end do
       if (numproc==meI) then
          cI(1:2,ld-(nMF+1):lf+nMF) = b%coord_nd(:,:,m0-b%m_charge(1,numproc)+1)
          do j = 1, NbCoupes(2)
             call MPI_SEND( b%coord_nd(:,:,m0-b%m_charge(1,numproc)+1), &
                  taille, MPI_REAL8, &
                  med+(j-1)*NbCoupes(1), MPI_FLAG, MPI_COMM_WORLD, STATINFO )
          end do
       else
          ! reception des points
          call MPI_RECV( cI(1:2,ld-(nMF+1):lf+nMF), taille, MPI_REAL8, &
               meI, MPI_FLAG, MPI_COMM_WORLD, MPI_status, STATINFO )     
       end if
    end if

  end subroutine recuperation_info_resserrement

end module m_raffinement_maillage
