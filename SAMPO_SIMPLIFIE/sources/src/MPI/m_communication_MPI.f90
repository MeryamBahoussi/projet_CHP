module m_communication_MPI
  use m_MPI
  use m_struct
  use parametres_globaux
  implicit none

  interface MPI_send_to_proc_ini
     module procedure MPI_send_scalaire_to_proc_ini  ! tableau de reel
     module procedure MPI_send_maillage_to_proc_ini  ! maillage
     module procedure MPI_send_matrice_to_proc_ini   ! matrice
     module procedure MPI_send_vecteur_to_proc_ini   ! tableau de vecteur
  end interface

contains

  subroutine MPI_send_scalaire_to_proc_ini(b, tab_in, tab_out, send_MF)
    implicit none
    
    type (STR_BLOC), intent(in) ::  b
    real*8, dimension(b%ld-nMF:b%lf+nMF, b%md-nMF:b%mf+nMF), intent(in) :: tab_in
    real*8, dimension(:,:), allocatable, intent(out) :: tab_out
    logical, intent(in) :: send_MF ! envoi des mailles fictives ?
    integer :: taille, me, proc_ini, proc_fin, nbprocs_b
    integer :: ld, lf, md, mf
    integer :: STATINFO, MPI_FLAG = 90
    integer :: nMF_loc
    integer, dimension(MPI_STATUS_SIZE) :: MPI_status
    integer, dimension(:,:), pointer :: l_pere, m_pere

    ld = b%ld; lf = b%lf
    md = b%md; mf = b%mf
    proc_ini = b%proc_ini
    proc_fin = b%proc_fin
    nbprocs_b = proc_fin-proc_ini+1
    l_pere => b%l_charge
    m_pere => b%m_charge

    nMF_loc = 0
    if (send_MF) nMF_loc = nMF
    
    if ( numproc == proc_ini ) then
       allocate(tab_out(b%ld_tot-nMF_loc:b%lf_tot+nMF_loc, b%md_tot-nMF_loc:b%mf_tot+nMF_loc))
       tab_out(l_pere(1,numproc)-nMF_loc:l_pere(2,numproc)+nMF_loc, m_pere(1,numproc)-nMF_loc:m_pere(2,numproc)+nMF_loc) = &
            tab_in(ld-nMF_loc:lf+nMF_loc, md-nMF_loc:mf+nMF_loc)
    end if
    
    nMF_loc = 0
    
    if (nbprocs_b /= 1) then
       if ( numproc /= proc_ini ) then 
          taille = (lf - ld+1+2*nMF_loc)*(mf - md+1+2*nMF_loc)
          call MPI_SEND( tab_in(ld-nMF_loc:lf+nMF_loc,md-nMF_loc:mf+nMF_loc), &  
               taille, MPI_REAL8, &
               proc_ini, MPI_FLAG, MPI_COMM_WORLD, STATINFO )
       else
          do me = proc_ini+1, proc_fin
             taille = (l_pere(2,me) - l_pere(1,me)+1+2*nMF_loc)*&
                  ( m_pere(2,me) - m_pere(1,me)+1+2*nMF_loc)
             call MPI_RECV( tab_out(l_pere(1,me)-nMF_loc:l_pere(2,me)+nMF_loc,&
                  m_pere(1,me)-nMF_loc:m_pere(2,me)+nMF_loc), &
                  taille, MPI_REAL8, &
                  me, MPI_FLAG, MPI_COMM_WORLD, MPI_status, STATINFO )
          end do
       end if
    end if
    
    if (send_MF) call MPI_send_scalaire_MF_to_proc_ini(b, tab_in, tab_out)
  end subroutine MPI_send_scalaire_to_proc_ini

  ! envoi des mailles fictives au proc ini
  subroutine MPI_send_scalaire_MF_to_proc_ini(b, tab_in, tab_out)
    type (STR_BLOC), intent(in) ::  b
    real*8, dimension(b%ld-nMF:b%lf+nMF, b%md-nMF:b%mf+nMF), intent(in) :: tab_in
    real*8, dimension(b%ld_tot-nMF:b%lf_tot+nMF, b%md_tot-nMF:b%mf_tot+nMF), intent(inout) :: tab_out
    integer :: taille, me, proc_ini, proc_fin, nbprocs_b
    integer :: ld, lf, md, mf
    integer :: STATINFO, MPI_FLAG = 91
    integer, dimension(MPI_STATUS_SIZE) :: MPI_status
    integer, dimension(:,:), pointer :: l_pere, m_pere

    ld = b%ld; lf = b%lf
    md = b%md; mf = b%mf
    proc_ini = b%proc_ini
    proc_fin = b%proc_fin
    nbprocs_b = proc_fin-proc_ini+1
    l_pere => b%l_charge
    m_pere => b%m_charge

    if ( numproc /= proc_ini ) then
       ! face 1 du bloc pere
       if (l_pere(1,b%numproc) == b%ld_tot) then
          taille = nMF*(mf-md+1)
          call MPI_SEND( tab_in(ld-nMF:ld-1,md:mf), &  
               taille, MPI_REAL8, &
               proc_ini, MPI_FLAG, MPI_COMM_WORLD, STATINFO )
       end if
       ! face 2 du bloc pere
       if (l_pere(2,b%numproc) == b%lf_tot) then
          taille = nMF*(mf-md+1)
          call MPI_SEND( tab_in(lf+1:lf+nMF,md:mf), &  
               taille, MPI_REAL8, &
               proc_ini, MPI_FLAG, MPI_COMM_WORLD, STATINFO )
       end if
       ! face 3 du bloc pere
       if (m_pere(1,b%numproc) == b%md_tot) then
          taille = nMF*(lf-ld+1)
          call MPI_SEND( tab_in(ld:lf,md-nMF:md-1), &  
               taille, MPI_REAL8, &
               proc_ini, MPI_FLAG, MPI_COMM_WORLD, STATINFO )
       end if
       ! face 4 du bloc pere
       if (m_pere(2,b%numproc) == b%mf_tot) then
          taille = nMF*(lf-ld+1)
          call MPI_SEND( tab_in(ld:lf,mf+1:mf+nMF), &  
               taille, MPI_REAL8, &
               proc_ini, MPI_FLAG, MPI_COMM_WORLD, STATINFO )
       end if
       ! coin 1-3 (en bas a gauche)
       if (l_pere(1,b%numproc) == b%ld_tot .and. m_pere(1,b%numproc) == b%md_tot) then
          taille = 4
          call MPI_SEND( tab_in(ld-nMF:ld-1,md-nMF:md-1), &  
               taille, MPI_REAL8, &
               proc_ini, MPI_FLAG, MPI_COMM_WORLD, STATINFO )
       end if
       ! coin 2-3 (en bas a droite)
       if (l_pere(2,b%numproc) == b%lf_tot .and. m_pere(1,b%numproc) == b%md_tot) then
          taille = 4
          call MPI_SEND( tab_in(lf+1:lf+nMF,md-nMF:md-1), &  
               taille, MPI_REAL8, &
               proc_ini, MPI_FLAG, MPI_COMM_WORLD, STATINFO )
       end if
       ! coin 1-4 (en haut a gauche)
       if (l_pere(1,b%numproc) == b%ld_tot .and. m_pere(2,b%numproc) == b%mf_tot) then
          taille = 4
          call MPI_SEND( tab_in(ld-nMF:ld-1,mf+1:mf+nMF), &  
               taille, MPI_REAL8, &
               proc_ini, MPI_FLAG, MPI_COMM_WORLD, STATINFO )
       end if
       ! coin 2-4 (en haut a droite)
       if (l_pere(2,b%numproc) == b%lf_tot .and. m_pere(2,b%numproc) == b%mf_tot) then
          taille = 4
          call MPI_SEND( tab_in(lf+1:lf+nMF,mf+1:mf+nMF), &  
               taille, MPI_REAL8, &
               proc_ini, MPI_FLAG, MPI_COMM_WORLD, STATINFO )
       end if
    else
       do me = proc_ini+1, proc_fin
          ! face 1 du bloc pere
          if (l_pere(1,me) == b%ld_tot) then
             taille = nMF*(m_pere(2,me)-m_pere(1,me)+1)
             call MPI_RECV( tab_out(l_pere(1,me)-nMF:l_pere(1,me)-1,&
                  m_pere(1,me):m_pere(2,me)), &
                  taille, MPI_REAL8, &
                  me, MPI_FLAG, MPI_COMM_WORLD, MPI_status, STATINFO )
          end if
          ! face 2 du bloc pere
          if (l_pere(2,me) == b%lf_tot) then
             taille = nMF*(m_pere(2,me)-m_pere(1,me)+1)
             call MPI_RECV( tab_out(l_pere(2,me)+1:l_pere(2,me)+nMF,&
                  m_pere(1,me):m_pere(2,me)), &
                  taille, MPI_REAL8, &
                  me, MPI_FLAG, MPI_COMM_WORLD, MPI_status, STATINFO )
          end if
          ! face 3 du bloc pere
          if (m_pere(1,me) == b%md_tot) then
             taille = nMF*(l_pere(2,me)-l_pere(1,me)+1)
             call MPI_RECV( tab_out(l_pere(1,me):l_pere(2,me),&
                  m_pere(1,me)-nMF:m_pere(1,me)-1), &
                  taille, MPI_REAL8, &
                  me, MPI_FLAG, MPI_COMM_WORLD, MPI_status, STATINFO )
          end if
          ! face 4 du bloc pere
          if (m_pere(2,me) == b%mf_tot) then
             taille = nMF*(l_pere(2,me)-l_pere(1,me)+1)
             call MPI_RECV( tab_out(l_pere(1,me):l_pere(2,me),&
                  m_pere(2,me)+1:m_pere(2,me)+nMF), &
                  taille, MPI_REAL8, &
                  me, MPI_FLAG, MPI_COMM_WORLD, MPI_status, STATINFO )
          end if
          ! coin 1-3 (en bas a gauche)
          if (l_pere(1,me) == b%ld_tot .and. m_pere(1,me) == b%md_tot) then
             taille = 4
             call MPI_RECV( tab_out(l_pere(1,me)-nMF:l_pere(1,me)-1,&
                  m_pere(1,me)-nMF:m_pere(1,me)-1), &
                  taille, MPI_REAL8, &
                  me, MPI_FLAG, MPI_COMM_WORLD, MPI_status, STATINFO )
          end if
          ! coin 2-3 (en bas a droite)
          if (l_pere(2,me) == b%lf_tot .and. m_pere(1,me) == b%md_tot) then
             taille = 4
             call MPI_RECV( tab_out(l_pere(2,me)+1:l_pere(2,me)+nMF,&
                  m_pere(1,me)-nMF:m_pere(1,me)-1), &
                  taille, MPI_REAL8, &
                  me, MPI_FLAG, MPI_COMM_WORLD, MPI_status, STATINFO )
          end if
          ! coin 1-4 (en haut a gauche)
          if (l_pere(1,me) == b%ld_tot .and. m_pere(2,me) == b%mf_tot) then
             taille = 4
             call MPI_RECV( tab_out(l_pere(1,me)-nMF:l_pere(1,me)-1,&
                  m_pere(2,me)+1:m_pere(2,me)+nMF), &
                  taille, MPI_REAL8, &
                  me, MPI_FLAG, MPI_COMM_WORLD, MPI_status, STATINFO )
          end if
          ! coin 2-4 (en haut a droite)
          if (l_pere(2,me) == b%lf_tot .and. m_pere(2,me) == b%mf_tot) then
             taille = 4
             call MPI_RECV( tab_out(l_pere(2,me)+1:l_pere(2,me)+nMF,&
                  m_pere(2,me)+1:m_pere(2,me)+nMF), &
                  taille, MPI_REAL8, &
                  me, MPI_FLAG, MPI_COMM_WORLD, MPI_status, STATINFO )
          end if
       end do
    end if
    
  end subroutine MPI_send_scalaire_MF_to_proc_ini

  subroutine MPI_send_maillage_to_proc_ini(b, n, tab_in, tab_out)
    implicit none
    
    type (STR_BLOC), intent(in) ::  b
    integer, intent(in) :: n
    real*8, dimension(n,b%ld-(nMF+1):b%lf+nMF, b%md-(nMF+1):b%mf+nMF), intent(in) :: tab_in
    real*8, dimension(:,:,:), allocatable, intent(out) :: tab_out
    integer :: taille, me, proc_ini, proc_fin, nbprocs_b
    integer :: ld, lf, md, mf
    integer :: STATINFO
    integer, dimension(MPI_STATUS_SIZE) :: MPI_status
    integer, dimension(:,:), pointer :: l_pere, m_pere

    ld = b%ld; lf = b%lf
    md = b%md; mf = b%mf
    proc_ini = b%proc_ini
    proc_fin = b%proc_fin
    nbprocs_b = proc_fin-proc_ini+1
    l_pere => b%l_charge
    m_pere => b%m_charge
    
    ! allocation du tableau de sortie tab_out pour le proc_ini
    ! On stocke la partie que connait le proc_ini dans tab_out
    if ( numproc == proc_ini ) then
       allocate(tab_out(1:n,b%ld_tot-(nMF+1):b%lf_tot+nMF, b%md_tot-(nMF+1):b%mf_tot+nMF))
       tab_out(:,l_pere(1,numproc)-(nMF+1):l_pere(2,numproc)+nMF, m_pere(1,numproc)-(nMF+1):m_pere(2,numproc)+nMF) = &
        tab_in(:,ld-(nMF+1):lf+nMF, md-(nMF+1):mf+nMF)
    end if

    ! si on a plusieurs proc sur le bloc (ceux donnes en entree)
    if ( nbprocs_b/=1 ) then
       if ( numproc /= proc_ini ) then 
          taille = n*(lf - ld+2*(nMF+1))*(mf - md+2*(nMF+1))
          call MPI_SEND( tab_in(1:n,ld-(nMF+1):lf+nMF,md-(nMF+1):mf+nMF), &  
               taille, MPI_REAL8, &
               proc_ini, 30, MPI_COMM_WORLD, STATINFO )
       else
          do me = proc_ini+1, proc_fin
             taille = n*(l_pere(2,me) - l_pere(1,me)+2*(nMF+1))*(m_pere(2,me) - m_pere(1,me)+2*(nMF+1))
             call MPI_RECV( tab_out(1:n,l_pere(1,me)-(nMF+1):l_pere(2,me)+nMF,m_pere(1,me)-(nMF+1):m_pere(2,me)+nMF), &
                  & taille, MPI_REAL8, &
                  & me, 30, MPI_COMM_WORLD, MPI_status, STATINFO )  
          end do
       end if
    end if
  end subroutine MPI_send_maillage_to_proc_ini

  subroutine MPI_send_matrice_to_proc_ini(b, n, mat_in, mat_out)
    implicit none
    
    type (STR_BLOC), intent(in) ::  b
    integer, intent(in) :: n
    real*8, dimension(n,n,b%ld:b%lf, b%md:b%mf), intent(in) :: mat_in
    real*8, dimension(:,:,:,:), pointer, intent(out) :: mat_out
    integer :: taille, me
    integer :: proc_ini, proc_fin, nbprocs_b
    integer :: ld, lf, md, mf
    integer :: STATINFO, MPI_FLAG = 100
    integer, dimension(MPI_STATUS_SIZE) :: MPI_status
    integer, dimension(:,:), pointer :: l_pere, m_pere

    ld = b%ld; lf = b%lf
    md = b%md; mf = b%mf
    proc_ini = b%proc_ini
    proc_fin = b%proc_fin
    nbprocs_b = proc_fin-proc_ini+1
    l_pere => b%l_charge
    m_pere => b%m_charge

    if ( numproc == proc_ini ) then
       allocate(mat_out(1:n,1:n,b%ld_tot:b%lf_tot, b%md_tot:b%mf_tot))
       mat_out(:,:,ld:lf, md:mf) = mat_in(:,:,ld:lf, md:mf)
    end if

    if ( nbprocs_b/=1 ) then
       if ( numproc /= proc_ini ) then 
          taille = n*n*(lf - ld+1)*(mf - md+1)
          call MPI_SEND( mat_in(:,:,ld:lf,md:mf), &  
               taille, MPI_REAL8, &
               proc_ini, MPI_FLAG, MPI_COMM_WORLD, STATINFO )
       else
          do me = proc_ini+1, proc_fin
             taille = n*n*(l_pere(2,me) - l_pere(1,me)+1)*(m_pere(2,me) - m_pere(1,me)+1)

             call MPI_RECV( mat_out(:,:,l_pere(1,me):l_pere(2,me),&
                  m_pere(1,me):m_pere(2,me)), &
                  taille, MPI_REAL8, &
                  me, MPI_FLAG, MPI_COMM_WORLD, MPI_status, STATINFO )  
          end do
       end if
    end if
  end subroutine MPI_send_matrice_to_proc_ini


  subroutine MPI_send_vecteur_to_proc_ini(b, n, vect_in, vect_out, send_MF)
    implicit none
    
    type (STR_BLOC), intent(in) ::  b
    integer, intent(in) :: n
   ! real*8, dimension(n,b%ld-2:b%lf+2, b%md-2:b%mf+2), intent(in) :: vect_in
    real*8, dimension(:,:,:), pointer :: vect_in

    real*8, dimension(:,:,:), allocatable, intent(out) :: vect_out
    logical, intent(in) :: send_MF ! envoi des mailles fictives ?
        
    integer :: taille, me, i
    integer :: proc_ini, proc_fin, nbprocs_b
    integer :: ld, lf, md, mf, nMF_loc
    integer :: STATINFO, MPI_FLAG = 110
    integer, dimension(MPI_STATUS_SIZE) :: MPI_status
    integer, dimension(:,:), pointer :: l_pere, m_pere
    
    ld = b%ld; lf = b%lf
    md = b%md; mf = b%mf
    proc_ini = b%proc_ini
    proc_fin = b%proc_fin
    nbprocs_b = proc_fin-proc_ini+1
    l_pere => b%l_charge
    m_pere => b%m_charge

    nMF_loc = 0
    if (send_MF) nMF_loc = nMF

    if ( numproc == b%proc_ini ) then
       allocate( vect_out(1:n,b%ld_tot-nMF_loc:b%lf_tot+nMF_loc, b%md_tot-nMF_loc:b%mf_tot+nMF_loc) ) 
       vect_out = -1.d70
       vect_out(:,l_pere(1,numproc)-nMF_loc:l_pere(2,numproc)+nMF_loc, m_pere(1,numproc)-nMF_loc:m_pere(2,numproc)+nMF_loc) = vect_in(:,ld-nMF_loc:lf+nMF_loc, md-nMF_loc:mf+nMF_loc)
    end if

    nMF_loc=0

    if ( nbprocs_b/=1 ) then
       if ( numproc /= b%proc_ini ) then 
          taille = n*(lf - ld+1+2*nMF_loc)*(mf - md+1+2*nMF_loc)
          call MPI_SEND( vect_in(1:n, ld-nMF_loc:lf+nMF_loc,md-nMF_loc:mf+nMF_loc), &  
               taille, MPI_REAL8, &
               b%proc_ini, MPI_FLAG, MPI_COMM_WORLD, STATINFO )
       else
          do me = b%proc_ini+1, b%proc_fin
             taille = n*(l_pere(2,me)-l_pere(1,me)+1+2*nMF_loc)*(m_pere(2,me)-m_pere(1,me)+1+2*nMF_loc)
             call MPI_RECV( vect_out(:,l_pere(1,me)-nMF_loc:l_pere(2,me)+nMF_loc,&
                  m_pere(1,me)-nMF_loc:m_pere(2,me)+nMF_loc), &
                  taille, MPI_REAL8, &
                  me, MPI_FLAG, MPI_COMM_WORLD, MPI_status, STATINFO )  
          end do
       end if
    end if
    
    if (send_MF) then
       do i = 1, n
          call MPI_send_scalaire_MF_to_proc_ini(b, vect_in(i,:,:), vect_out(i,:,:))
       end do
    end if

  end subroutine MPI_send_vecteur_to_proc_ini

  
  subroutine MPI_communication_scalaire(indices, numproc_adj, tab) 
    
    integer, dimension(1:8), intent(in) :: indices
    integer                , intent(in) :: numproc_adj
    real*8 , dimension(:,:), pointer    :: tab

    
    integer :: ldi, lfi, mdi, mfi
    integer :: lde, lfe, mde, mfe
    integer :: taille
    integer :: STATINFO
    integer :: MPI_status(MPI_STATUS_SIZE)
 
    ldi = indices(1); lfi = indices(2)
    mdi = indices(3); mfi = indices(4)
    lde = indices(5); lfe = indices(6)
    mde = indices(7); mfe = indices(8)
    
    taille = (lfi - ldi + 1) * (mfi - mdi + 1)
     
    if (numproc_adj > numproc) then
       call MPI_SEND( tab(ldi:lfi, mdi:mfi), &  
            taille, MPI_REAL8, &
            numproc_adj, 20+numproc, MPI_COMM_WORLD, STATINFO )
       call mpi_RECV( tab(lde:lfe,mde:mfe), & 
            taille, MPI_REAL8, &
            numproc_adj, 20+numproc_adj, MPI_COMM_WORLD, MPI_STATUS, STATINFO )
    else
       call mpi_RECV( tab(lde:lfe,mde:mfe), & 
            taille, MPI_REAL8, &
            numproc_adj, 20+numproc_adj, MPI_COMM_WORLD, MPI_STATUS, STATINFO )
       call MPI_SEND( tab(ldi:lfi, mdi:mfi), &  
            taille, MPI_REAL8, &
            numproc_adj, 20+numproc, MPI_COMM_WORLD, STATINFO )
    end if
  end subroutine MPI_communication_scalaire
    
  subroutine MPI_communication_vecteur(indices, numproc_adj, tab) 
    
    integer, dimension(1:8)  , intent(in) :: indices
    integer                  , intent(in) :: numproc_adj
    real*8 , dimension(:,:,:), pointer    :: tab

    integer :: i
    integer :: ldi, lfi, mdi, mfi
    integer :: lde, lfe, mde, mfe
    integer :: taille
    integer :: STATINFO
    integer :: MPI_status(MPI_STATUS_SIZE)
    
    ldi = indices(1); lfi = indices(2)
    mdi = indices(3); mfi = indices(4)
    lde = indices(5); lfe = indices(6)
    mde = indices(7); mfe = indices(8)
    
    taille = (lfi - ldi + 1) * (mfi - mdi + 1)
     
    if (numproc_adj > numproc) then
       do i = 1, size(tab,1)
          call MPI_SEND( tab(i,ldi:lfi, mdi:mfi), &  
               taille, MPI_REAL8, &
               numproc_adj, 20+numproc+i, MPI_COMM_WORLD, STATINFO )
          
          call mpi_RECV( tab(i,lde:lfe,mde:mfe), & 
               taille, MPI_REAL8, &
               numproc_adj, 20+numproc_adj+i, MPI_COMM_WORLD, MPI_STATUS, STATINFO )
       end do
    else
       do i = 1, size(tab,1)
          call mpi_RECV( tab(i,lde:lfe,mde:mfe), & 
               taille, MPI_REAL8, &
               numproc_adj, 20+numproc_adj+i, MPI_COMM_WORLD, MPI_STATUS, STATINFO )

          call MPI_SEND( tab(i,ldi:lfi, mdi:mfi), &  
               taille, MPI_REAL8, &
               numproc_adj, 20+numproc+i, MPI_COMM_WORLD, STATINFO )         
       end do
    end if
  end subroutine MPI_communication_vecteur
    
end module m_communication_MPI



