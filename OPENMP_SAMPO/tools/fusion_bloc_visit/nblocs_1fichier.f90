! source .../env_code.sh
! ifort -o nblocs_1fichier nblocs_1fichier.f90 $NETCDF_LDFLAGS -lnetcdff -lnetcdf -I$NETCDF_INC
! gcc -o nblocs_1fichier nblocs_1fichier.f90 $NETCDF_LDFLAGS -lnetcdff -lnetcdf -J$NETCDF_INC
! ./nblocs_1fichier

module m_fonction
  implicit none
  
  include 'netcdf.inc'
  
  logical, save :: temperature =.false.

  type donnees_fichier
     integer                            :: lf, mf        
     real*8                             :: temps
     real*8 , dimension(:,:,:), pointer :: centre_nd
     
     integer                                      :: nb_Rtabs, nb_Itabs
     real*8           , dimension(:,:,:), pointer :: Rtabs
     integer          , dimension(:,:,:), pointer :: Itabs
     character(len=20), dimension(:)    , pointer :: nom_Rtabs, nom_Itabs
     
  end type donnees_fichier

  type acces_donnees  ! structure d'acces aux donnees
     type (donnees_fichier), pointer :: pointeur
  end type acces_donnees
  
  interface ecrit_quantite
     module procedure    ecr_qutt_3D_rt, ecr_qutt_2D_rt, ecr_qutt_2D_it
     module procedure    ecr_qutt_0D_rt, ecr_qutt_0D_it
     module procedure    ecr_qutt_2D_rt_vtk, ecr_qutt_2D_it_vtk
  end interface

contains
  
  subroutine allocate_donnees(d)
    implicit none
    
    type(donnees_fichier), pointer :: d
    
    allocate(d%centre_nd(2,0:d%lf, 0:d%mf))
    allocate(d%Itabs(d%nb_Itabs,1:d%lf, 1:d%mf))
    allocate(d%Rtabs(d%nb_Rtabs,1:d%lf, 1:d%mf))
    allocate(d%nom_Rtabs(d%nb_Rtabs))
    allocate(d%nom_Itabs(d%nb_Itabs))
  end subroutine allocate_donnees

  subroutine deallocate_donnees(d)
    implicit none
    type(donnees_fichier), pointer :: d
    deallocate(d%centre_nd,d%Itabs,d%Rtabs,d%nom_Itabs,d%nom_Rtabs)
  end subroutine deallocate_donnees
  
  subroutine trouve_nb_variables(fichier, d)
    
    character(len=*), intent(in) :: fichier
    type(donnees_fichier), pointer :: d

    integer :: lf, mf
    integer :: l,m
    integer :: iunit = 20
    real*8 , pointer :: Rtoto(:,:)
    integer, pointer :: Itoto(:,:)
    character(len=10) :: toto
    character(len=20) :: rang_var, nom_var, type_var 
    integer :: err

    open(unit=20, file=trim(fichier)//'.vtk')
    read(iunit,'(1A26)') ! '# vtk DataFile Version 2.0'
    read(iunit,*) !'Flow'
    read(iunit,*) !'ASCII'
    read(iunit,*) !'DATASET STRUCTURED_GRID'
    read(iunit,*) !'FIELD FieldData 1'
    read(iunit,*) !'TIME 1 1 double'
    read(iunit,*) d%temps
    read(iunit,*) toto, d%lf, d%mf!, 1
    read(iunit,*) !'POINTS', (lf-1+2)*(mf-1+2), 'double'

    d%lf = d%lf-1; d%mf = d%mf-1
    lf = d%lf; mf = d%mf
    allocate(Rtoto(1:lf,1:mf))
    allocate(Itoto(1:lf,1:mf))
    
    do m=0,mf
       do l=0,lf
          read(iunit,*) Rtoto(1,1), Rtoto(lf,mf)
       end do
    end do
    read(iunit,*) !'CELL_DATA', (lf-1+1)*(mf-1+1)
    
    d%nb_Rtabs = 0
    d%nb_Itabs = 0
    do  
       read(iunit,*, iostat=err) rang_var, nom_var, type_var
       if (err==-1) exit
       if (trim(rang_var) == "VECTORS") then
          do m=1,mf
             do l=1,lf
                read(iunit,*) Rtoto(l,m), Rtoto(l,m), Rtoto(l,m)
             end do
          end do
       else if (trim(rang_var) == "SCALARS") then
          read(iunit,*) !'LOOKUP_TABLE default'
          if (trim(type_var) == "double") then
             read(iunit,*) Rtoto(1:lf,1:mf)
             d%nb_Rtabs = d%nb_Rtabs+1
          else if (trim(type_var) == "int") then
             read(iunit,*) Itoto(1:lf,1:mf)
             d%nb_Itabs = d%nb_Itabs+1
          else
             print*, 'fin anormale du fichier ou type de variables inconnu'
             stop
          end if
       else
          print*, 'fin anormale du fichier ou type de variables inconnu'
          stop
       end if
    end do
    
    deallocate(Rtoto, Itoto)
    close(20)
    
  end subroutine trouve_nb_variables
  
  subroutine lectureFichierVtk(fichier, d)
    
    character(len=*), intent(in) :: fichier
    type(donnees_fichier), pointer :: d

    integer :: lf, mf
    integer :: l,m, indI, indR
    integer :: iunit = 20
    character(len=10) :: toto
    character(len=20) :: rang_var, nom_var, type_var 
    integer :: err
    
    open(unit=20, file=trim(fichier)//'.vtk')
    read(iunit,'(1A26)') ! '# vtk DataFile Version 2.0'
    read(iunit,*) !'Flow'
    read(iunit,*) !'ASCII'
    read(iunit,*) !'DATASET STRUCTURED_GRID'
    read(iunit,*) !'FIELD FieldData 1'
    read(iunit,*) !'TIME 1 1 double'
    read(iunit,*) d%temps
    read(iunit,*) toto, d%lf, d%mf!, 1
    read(iunit,*) !'POINTS', (lf-1+2)*(mf-1+2), 'double'

    d%lf = d%lf-1; d%mf = d%mf-1
    
    call allocate_donnees(d)
    lf = d%lf; mf = d%mf
    do m=0,mf
       do l=0,lf
          read(iunit,*) d%centre_nd(1,l,m), d%centre_nd(2,l,m)!, 0.0d0
       end do
    end do
    
    read(iunit,*) !'CELL_DATA', (lf-1+1)*(mf-1+1)
    
    indR = 0
    indI = 0
    do  
       read(iunit,*, iostat=err) rang_var, nom_var, type_var
       if (err==-1) exit
       if (trim(rang_var) == "VECTORS") then
          ! on ne recupere pas les vecteurs
          do m=1,mf
             do l=1,lf
                read(iunit,*) toto, toto, toto
             end do
          end do
       else if (trim(rang_var) == "SCALARS") then
          read(iunit,*) !'LOOKUP_TABLE default'
          if (trim(type_var) == "double") then
             indR = indR + 1
             read(iunit,*) d%Rtabs(indR,1:lf,1:mf)
             d%nom_Rtabs(indR) = trim(nom_var)
          else if (trim(type_var) == "int") then
             indI = indI + 1
             read(iunit,*) d%Itabs(indI,1:lf,1:mf)
             d%nom_Itabs(indI) = trim(nom_var)
          else
             print*, 'fin anormale du fichier ou type de variables inconnu'
             stop
          end if
       end if
    end do
        
    close(20)
    
  end subroutine lectureFichierVtk
  
  subroutine assemble(d_tab,d_global, nb_proc_x, nb_proc_y)
    integer, intent(in) :: nb_proc_x, nb_proc_y
    type (acces_donnees), dimension(nb_proc_x*nb_proc_y) :: d_tab
    type(donnees_fichier), pointer :: d_global, d_loc
    
    integer :: me, i, lf, mf
    integer :: ll(1:2,1:nb_proc_x*nb_proc_y), mm(1:2,1:nb_proc_x*nb_proc_y)  ! indice du tableau global
    integer :: decalage_l, decalage_m
    integer :: nbproc
     
    print*, "debut assemble"

    nbproc = nb_proc_x*nb_proc_y

!!! on retrouve les lf et mf
    decalage_l = 0; decalage_m = 0
    do me = 1, nbproc
       d_loc => d_tab(me)%pointeur
       ll(1, me) = 1 + decalage_l
       ll(2, me) = d_loc%lf + decalage_l
       mm(1, me) = 1 + decalage_m
       mm(2, me) = d_loc%mf + decalage_m
       
       decalage_l = decalage_l + d_loc%lf
       if (mod(me,nb_proc_x)==0) then
          decalage_l = 0
          decalage_m = mm(2, me)
       end if
    end do
        
    d_global%lf = ll(2,nbproc)
    d_global%mf = mm(2,nbproc)
    d_global%nb_Itabs = d_tab(1)%pointeur%nb_Itabs
    d_global%nb_Rtabs = d_tab(1)%pointeur%nb_Rtabs
    call allocate_donnees(d_global)
   
    d_global%temps= d_tab(1)%pointeur%temps
    
    do me = 1, nbproc
       d_loc => d_tab(me)%pointeur
       lf = d_loc%lf; mf = d_loc%mf
       
       ! centres des mailles
       d_global%centre_nd(:,ll(1,me)-1:ll(2,me),mm(1,me)-1:mm(2,me)) = d_loc%centre_nd(:,0:lf,0:mf)
       
       ! tableaux de reels
       do i = 1, d_global%nb_Rtabs
          d_global%Rtabs(i,ll(1,me):ll(2,me),mm(1,me):mm(2,me))=d_loc%Rtabs(i,1:lf,1:mf)
          d_global%nom_Rtabs(i) = d_loc%nom_Rtabs(i)
       end do
       
       ! tableaux d'entiers
       do i = 1, d_global%nb_Itabs
          d_global%Itabs(i,ll(1,me):ll(2,me),mm(1,me):mm(2,me))=d_loc%Itabs(i,1:lf,1:mf)
          d_global%nom_Itabs(i) = d_loc%nom_Itabs(i)
       end do
    end do
    
    print*, "fin assemble"
  end subroutine assemble
  
  subroutine sorties_vtk(fichier, d)
    character(len=*)     , intent(in) :: fichier
    type(donnees_fichier), pointer    :: d
    
    integer :: lf, mf, l, m, dim(2), i
    integer :: iunit = 20
    
    lf = d%lf; mf = d%mf
    
    open(unit=iunit, file=trim(fichier)//'.vtk')
    write(iunit,'(1A26)') '# vtk DataFile Version 2.0'
    write(iunit,*) 'Flow'
    write(iunit,*) 'ASCII'
    write(iunit,*) 'DATASET STRUCTURED_GRID'
    write(iunit,*) 'FIELD FieldData 1'
    write(iunit,*) 'TIME 1 1 double'
    write(iunit,*)  d%temps
    write(iunit,*) 'DIMENSIONS', lf+1, mf+1, 1
    write(iunit,*) 'POINTS', (lf+1)*(mf+1), 'double'
    do m=0,mf
       do l=0,lf
          write(iunit,*) d%centre_nd(1,l,m), d%centre_nd(2,l,m), 0.0d0
       end do
    end do
    write(iunit,*) 'CELL_DATA', (lf-1+1)*(mf-1+1)
    dim = (/lf,mf/)
    do i = 1, d%nb_Itabs
       call ecrit_quantite(iunit,dim,d%Itabs(i,:,:),d%nom_Itabs(i))
    end do
    do i = 1, d%nb_Rtabs
       call ecrit_quantite(iunit,dim,d%Rtabs(i,:,:),d%nom_Rtabs(i))
    end do
    close(20)
  end subroutine sorties_vtk
  
  subroutine sortie_netcdf_definition_bloc(ncid, d, ch_num, XYZ_id, lf_id, mf_id, lfp1_id, mfp1_id)
    integer              , intent(in)  :: ncid
    type(donnees_fichier), pointer     :: d
    character(len=4)     , intent(in)  :: ch_num
    integer              , intent(in)  :: XYZ_id
    integer              , intent(out) :: lf_id, mf_id, lfp1_id, mfp1_id

    integer               :: i, iret, iret_glob, var_id
    integer, dimension(2) :: dim_ids
    character(len=20)     :: nom

    iret_glob = NF_NOERR
    iret = NF_DEF_DIM(ncid, 'LF_'//trim(adjustl(ch_num)), d%lf, lf_id)
    iret_glob = iret_glob + iret
    iret = NF_DEF_DIM(ncid, 'MF_'//trim(adjustl(ch_num)), d%mf, mf_id)
    iret_glob = iret_glob + iret
    iret = NF_DEF_DIM(ncid, 'LF_p_1_'//trim(adjustl(ch_num)), d%lf+1, lfp1_id)
    iret_glob = iret_glob + iret
    iret = NF_DEF_DIM(ncid, 'MF_p_1_'//trim(adjustl(ch_num)), d%mf+1, mfp1_id)
    iret_glob = iret_glob + iret
    call faute('definition dimension ', iret)

    nom = "maille_bloc_"//trim(adjustl(ch_num))
    iret = NF_DEF_VAR(ncid, nom, NF_DOUBLE, 3, (/lfp1_id,mfp1_id,XYZ_id/), var_id)
    iret = NF_PUT_ATT_TEXT(ncid, var_id, 'long_name', len(nom), nom)
    call faute('definition maillage ', iret)

    dim_ids =(/lf_id, mf_id/)
    do i = 1, d%nb_Rtabs
       nom = trim(d%nom_Rtabs(i))//"_bloc_"//trim(adjustl(ch_num))
       iret = NF_DEF_VAR(ncid, nom, NF_DOUBLE, 2, dim_ids, var_id)
       iret = NF_PUT_ATT_TEXT(ncid, var_id, 'long_name', len(nom), nom)   
       call faute('definition '//nom, iret)
    end do
    do i = 1, d%nb_Itabs
       nom = trim(d%nom_Rtabs(i))//"_bloc_"//trim(adjustl(ch_num))
       iret = NF_DEF_VAR(ncid, nom, NF_INT, 2, dim_ids, var_id)
       iret = NF_PUT_ATT_TEXT(ncid, var_id, 'long_name', len(nom), nom)
       call faute('definition '//nom, iret)
    end do

  end subroutine sortie_netcdf_definition_bloc
  
  subroutine sortie_netcdf_ecriture_bloc(ncid, d, ch_num, XYZ_id, lf_id, mf_id, lfp1_id, mfp1_id)
    integer              , intent(in)  :: ncid
    type(donnees_fichier), pointer     :: d
    character(len=4)     , intent(in)  :: ch_num
    integer              , intent(in)  :: XYZ_id
    integer              , intent(in)  :: lf_id, mf_id, lfp1_id, mfp1_id

    integer               :: i, iret, var_id
    integer, dimension(2) :: dim_ids
    character(len=20)     :: nom
    
    nom = "maille_bloc_"//trim(adjustl(ch_num))
    call ecrit_quantite(ncid, iret, (/lfp1_id,mfp1_id,XYZ_id/), d%centre_nd, nom)
    call faute('ecriture maillage ', iret)

    dim_ids =(/lf_id, mf_id/)
    do i = 1, d%nb_Rtabs
       nom = trim(d%nom_Rtabs(i))//"_bloc_"//trim(adjustl(ch_num))
       call ecrit_quantite(ncid, iret, dim_ids, d%Rtabs(i,:,:), nom)
       call faute('ecriture '//nom, iret)
    end do
    do i = 1, d%nb_Itabs
       nom = trim(d%nom_Itabs(i))//"_bloc_"//trim(adjustl(ch_num))
       call ecrit_quantite(ncid, iret, dim_ids, d%Itabs(i,:,:), nom) 
       call faute('ecriture '//nom, iret)
    end do

  end subroutine sortie_netcdf_ecriture_bloc

  subroutine sorties_netcdf(fichier, d)
    character(len=*)     , intent(in) :: fichier
    type(donnees_fichier), pointer    :: d

    integer :: iret
    integer :: ncid = 1
    integer :: lf_id, mf_id, lfp1_id, mfp1_id, i
    integer :: XYZ_id
    integer :: tempsid, altid, iterationid, var_id
    character(len=20) :: nom

    !   creation
    iret = NF_CREATE(trim(fichier)//'.nc', OR(NF_CLOBBER, NF_64BIT_OFFSET), ncid)
    call faute('ouverture fichier '//trim(fichier)//'.nc', iret)
   
    ! definition
    iret = NF_DEF_DIM(ncid, 'XYZ', 2, XYZ_id)
    
    call sortie_netcdf_definition_bloc(ncid, d, "1   ", XYZ_id, lf_id, mf_id, lfp1_id, mfp1_id)

    call ecrit_quantite(ncid, iret, tempsid, d%temps, "temps")
    call ecrit_quantite(ncid, iret, iterationid, 1, "iteration")
    call ecrit_quantite(ncid, iret, altid, 0.d0, "altitude")

    call sortie_netcdf_ecriture_bloc(ncid, d, "1   ", XYZ_id, lf_id, mf_id, lfp1_id, mfp1_id)
 
    iret = NF_CLOSE(ncid)
    call faute('fermeture fichier '//trim(fichier)//'.nc', iret)
    
  end subroutine sorties_netcdf


  subroutine sorties_netcdf_nblocs(fichier, d_tab, nb_proc_x, nb_proc_y)
    character(len=*)    , intent(in) :: fichier
    integer             , intent(in) :: nb_proc_x, nb_proc_y
    type (acces_donnees), intent(in) :: d_tab(nb_proc_x*nb_proc_y)

    integer :: iret
    integer :: ncid = 1
    integer :: i, me
    integer :: tempsid, altid, iterationid, var_id
    integer :: XYZ_id
    integer :: nbproc  
    integer, dimension(:), allocatable :: lf_id, mf_id, lfp1_id, mfp1_id
    character(len=20) :: nom
    character(len=4) :: ch_num

    nbproc = nb_proc_x*nb_proc_y

    !   creation
    iret = NF_CREATE(trim(fichier)//'.nc', OR(NF_CLOBBER, NF_64BIT_OFFSET), ncid)
    call faute('ouverture fichier '//trim(fichier)//'.nc', iret)
        
    ! definition
    iret = NF_DEF_DIM(ncid, 'XYZ', 2, XYZ_id)

    allocate(lf_id(nbproc))
    allocate(mf_id(nbproc))
    allocate(lfp1_id(nbproc))
    allocate(mfp1_id(nbproc))
    ! definition des dimensions et des variables
    do me = 1, nbproc
       write(ch_num, FMT='(I4)') me       
       call sortie_netcdf_definition_bloc(ncid, d_tab(me)%pointeur, ch_num, XYZ_id, lf_id(me), mf_id(me), lfp1_id(me), mfp1_id(me))
    end do

    call ecrit_quantite(ncid, iret, tempsid, d_tab(1)%pointeur%temps, "temps")
    call ecrit_quantite(ncid, iret, iterationid, 1, "iteration")
    call ecrit_quantite(ncid, iret, altid, 0.d0, "altitude")

    do me = 1, nbproc
       write(ch_num, FMT='(I4)') me
       call sortie_netcdf_ecriture_bloc(ncid, d_tab(me)%pointeur, ch_num, XYZ_id, lf_id(me), mf_id(me), lfp1_id(me), mfp1_id(me))      
    end do
    deallocate(lf_id, mf_id, lfp1_id, mfp1_id)

    iret = NF_CLOSE(ncid)
    call faute('fermeture fichier '//trim(fichier)//'.nc', iret)

  end subroutine sorties_netcdf_nblocs
  
  subroutine ecr_qutt_3D_rt(ncid, iret, dimids, var, var_text) 
    integer                 , intent(in)    :: ncid, dimids(3)
    integer                 , intent(inout) :: iret
    real*8, dimension(:,:,:), intent(in)    :: var
    character(len=*)        , intent(in)    :: var_text
    
    integer :: var_id
    
    var_id = -1
    iret = NF_INQ_VARID (ncid, var_text, var_id)
    call faute(var_text, iret)
    iret = NF_PUT_VARA_DOUBLE(ncid, var_id, (/1,1,1/), (/ubound(var,2),ubound(var,3),1/), var(1,:,:))  
    iret = NF_PUT_VARA_DOUBLE(ncid, var_id, (/1,1,2/), (/ubound(var,2),ubound(var,3),1/), var(2,:,:))
        
  end subroutine ecr_qutt_3D_rt
  
  subroutine ecr_qutt_2D_rt(ncid, iret, dimids, var, var_text) 
    integer               , intent(in)    :: ncid, dimids(2)
    integer               , intent(inout) :: iret
    real*8, dimension(:,:), intent(in)    :: var
    character(len=*)      , intent(in)    :: var_text
    
    integer :: var_id
        
    var_id = -1
    iret = NF_INQ_VARID (ncid, var_text, var_id)
    call faute(var_text, iret)
    iret = NF_PUT_VARA_DOUBLE(ncid, var_id, (/1,1/), (/ubound(var,1),ubound(var,2)/), var(:,:))
  end subroutine ecr_qutt_2D_rt
  
  subroutine  ecr_qutt_2D_it(ncid, iret, dimids, var, var_text) 
    integer                , intent(in)    :: ncid, dimids(2)
    integer                , intent(inout) :: iret
    integer, dimension(:,:), intent(in)    :: var
    character(len=*)       , intent(in)    :: var_text
    
    integer :: var_id
    
    var_id = -1
    iret = NF_INQ_VARID (ncid, var_text, var_id)
    call faute(var_text, iret)
    iret = NF_PUT_VARA_INT(ncid, var_id, (/1,1/), (/ubound(var,1),ubound(var,2)/), var(:,:))    
  end subroutine ecr_qutt_2D_it
  
  subroutine ecr_qutt_0D_rt(ncid, iret, dim_id, var, var_text) 
    integer         , intent(in)    :: ncid, dim_id
    integer         , intent(inout) :: iret
    real*8          , intent(in)    :: var
    character(len=*), intent(in)    :: var_text
    
    integer :: var_id
    
    iret = NF_REDEF(ncid)
    iret = NF_DEF_VAR(ncid, var_text, NF_DOUBLE, 0, dim_id, var_id)
    iret = NF_PUT_ATT_TEXT(ncid, var_id, 'long_name', len(var_text), var_text)
    iret = NF_ENDDEF(ncid)
    iret = NF_PUT_VAR(ncid, var_id, var)
  end subroutine ecr_qutt_0D_rt
  
  subroutine ecr_qutt_0D_it(ncid, iret, dim_id, var, var_text) 
    integer         , intent(in)    :: ncid, dim_id
    integer         , intent(inout) :: iret
    integer         , intent(in)    :: var
    character(len=*), intent(in)    :: var_text
    
    integer :: var_id
    
    iret = NF_REDEF(ncid)
    iret = NF_DEF_VAR(ncid, var_text, NF_INT, 0, dim_id, var_id)
    iret = NF_PUT_ATT_TEXT(ncid, var_id, 'long_name', len(var_text), var_text)
    iret = NF_ENDDEF(ncid)
    iret = NF_PUT_VAR(ncid, var_id, var)
  end subroutine ecr_qutt_0D_it

  subroutine ecr_qutt_2D_rt_vtk(file_id, dim, var, var_text) 
    integer                         , intent(in) :: file_id, dim(2)
    real*8          , dimension(:,:), intent(in) :: var
    character(len=*)                , intent(in) :: var_text
    
    integer :: l, m

    write(file_id,*) 'SCALARS '//trim(adjustl(var_text))//' double'
    write(file_id,*) 'LOOKUP_TABLE default'
    do m=1,dim(2)
       do l=1,dim(1)
          write(file_id,*) var(l,m)
       end do
    end do
  end subroutine ecr_qutt_2D_rt_vtk
  
  subroutine ecr_qutt_2D_it_vtk(file_id, dim, var, var_text) 
    integer                , intent(in) :: file_id, dim(2)
    integer, dimension(:,:), intent(in) :: var
    character(len=*)       , intent(in) :: var_text

    integer :: l, m

    write(file_id,*) 'SCALARS '//trim(adjustl(var_text))//' int'
    write(file_id,*) 'LOOKUP_TABLE int'
    do m=1,dim(2)
       do l=1,dim(1)
          write(file_id,*) var(l,m)
       end do
    end do
  end subroutine ecr_qutt_2D_it_vtk

  subroutine faute(name, iret)
    character(len=*)         , intent (in)  :: name
    integer                  , intent (in)  :: iret

    if (iret .NE. NF_NOERR) then
       print *, name, ': Attention : ', Trim(NF_STRERROR(iret))
    end if
  end subroutine faute

  subroutine deallocate_toutes_donnees(d_tab, d_global, ip_ini, nbproc)
    type (acces_donnees), dimension(:), allocatable :: d_tab
    type(donnees_fichier), pointer :: d_global
    integer, intent(in) :: ip_ini, nbproc
    
    integer :: ip,  ip_loc
    
    do ip = ip_ini, nbproc
       ip_loc = ip - ip_ini + 1
       call deallocate_donnees(d_tab(ip_loc)%pointeur)
    end do
    call deallocate_donnees(d_global)
    
  end subroutine deallocate_toutes_donnees
  
  subroutine nbblocs_to_unfichier(prefixe, suffixe, ib, ip_ini, nb_proc_x, nb_proc_y, type_output)
    character(len=*), intent(in) :: prefixe, suffixe
    integer         , intent(in) :: ib, ip_ini, nb_proc_x, nb_proc_y
    integer         , intent(in) :: type_output

    type(acces_donnees)  , allocatable :: d_tab(:)
    type(donnees_fichier), pointer     :: d_global
    integer            :: ip, ip_loc
    integer            :: nbproc
    character(len=300) :: fichier
    character(len=4)   :: ip_c
    character(len=2)   :: ib_c

    nbproc = nb_proc_x * nb_proc_y
    
    allocate(d_tab(nbproc))
    do ip = ip_ini, nbproc+ip_ini-1
       ip_loc = ip - ip_ini + 1
       allocate(d_tab(ip_loc)%pointeur)
       write(ip_c,'(i4.4)') ip
       fichier=prefixe//'_'//trim(ip_c)//'_'//suffixe
       
       if (ip_loc == 1) then
          call trouve_nb_variables(fichier, d_tab(ip_loc)%pointeur)
       else
          d_tab(ip_loc)%pointeur%nb_Rtabs = d_tab(1)%pointeur%nb_Rtabs
          d_tab(ip_loc)%pointeur%nb_Itabs = d_tab(1)%pointeur%nb_Itabs
       end if
       call lectureFichierVtk(fichier, d_tab(ip_loc)%pointeur)
    end do

    if (ib == 1) then
       fichier=prefixe//'_'//suffixe
    else
       write(ib_c,'(i2.2)') ib
       fichier=prefixe//'_'//trim(ib_c)//'_'//suffixe
    end if
    
    select case (type_output)
    case(1) ! .nc avec n bloc
       call sorties_netcdf_nblocs(fichier, d_tab,nb_proc_x,nb_proc_y)
    case(2) ! .nc avec un seul bloc
       allocate(d_global)
       call assemble(d_tab,d_global,nb_proc_x,nb_proc_y)
       call sorties_netcdf(fichier, d_global)
       call deallocate_donnees(d_global)
    case(3) ! .vtk
       allocate(d_global)
       call assemble(d_tab,d_global,nb_proc_x,nb_proc_y)
       call sorties_vtk(fichier, d_global)
       call deallocate_donnees(d_global)
    end select

    do ip = ip_ini, nbproc
       ip_loc = ip - ip_ini + 1
       call deallocate_donnees(d_tab(ip_loc)%pointeur)
    end do
  end subroutine nbblocs_to_unfichier
  
end module m_fonction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
  use m_fonction
  
  implicit none
  
  integer            :: it, i
  integer            :: itdeb=0
  integer            :: itfin=-1
  integer            :: nbproc=1 
  integer            :: nb_proc_x=1
  integer            :: nb_proc_y=1
  integer            :: ip_ini = 1
  integer            :: num_bloc=1
  integer            :: type_output=2 
  logical            :: ini=.false.
  logical            :: fin=.false.
  logical            :: error=.false.
  logical            :: fluide=.true. ! si false, c'est le solide 
  character(len=5)   :: suffixe
  character(len=400) :: path_donnees="."
  character(len=400) :: arg, arg2

!!!!!!!!!
!!  recuperation des arguments
!!!!!!!!!
  i = 1
  do while (i < iargc()+1)
     call getarg(i, arg)
     select case (trim(arg))
     case ("-n")
        i = i+1
        call getarg(i, arg2); read(arg2,'(i4)') nbproc
     case("-d")
        i = i+1
        call getarg(i, path_donnees)
     case("--iter_deb")
        i = i+1
        call getarg(i, arg2); read(arg2,'(i3)') itdeb
     case("--iter_fin")
        i = i+1
        call getarg(i, arg2); read(arg2,'(i3)') itfin
     case ("--ini")
        ini=.true.
     case("--fin")
        fin=.true.
     case("--error")
        error=.true.
     case("--num_bloc")
        i = i+1
        call getarg(i, arg2); read(arg2,'(i3)') num_bloc
     case("--fluide")
        fluide=.true.
     case("--solide")
        fluide=.false.
     case ("--iproc_ini")
        i = i+1
        call getarg(i, arg2); read(arg2,'(i3)') ip_ini
     case("--output")
        i = i+1
        call getarg(i, arg2); read(arg2,'(i3)') type_output
     case("-h", "--help")
        write(6,"(A)")
        write(6,"(A)") "Utilisation : nblocs_1fichier [OPTIONS]"
        write(6,"(A)") "Script permettant de fusionner plusieurs fichiers de sorties en un seul."
        write(6,"(A)") "Il faut se placer dans le dossier du cas test (pas dans le dossier visit)."
        write(6,"(A)")
        write(6,"(A)") "   -n          : nombre de procs                                              (défaut : 1)"
        write(6,"(A)") "   -d          : chemin relatif vers le repertoire du cas                     (défaut : .)"
        write(6,"(A)") "   --iter_deb  : iteration initiale                                           (défaut : 0)"
        write(6,"(A)") "   --iter_fin  : iteration finale                                             (défaut : -1)"
        write(6,"(A)") "   --ini       : gestion de la sortie initiale                                (défaut : désactivé)"
        write(6,"(A)") "   --fin       : gestion de la sortie finale                                  (défaut : désactivé)"
        write(6,"(A)") "   --error     : gestion de la sortie error                                   (défaut : désactivé)"
        write(6,"(A)") "   --num_bloc  : numero du bloc à traiter                                     (défaut : 1)"
        write(6,"(A)") "   --fluide    : pour des sorties de fluide                                   (défaut : activé)"
        write(6,"(A)") "   --solide    : pour des sorties de solide                                   (défaut : désactivé)"
        write(6,"(A)") "   --iproc_ini : indice du proc initial, utile pour les calculs d'aérothermie (défaut : 1)"
        write(6,"(A)") "   --output    : type d'extension du fichier de sortie" 
        write(6,"(A)") "               : 1 .nc avec nbprocs blocs                                     (défaut)"
        write(6,"(A)") "               : 2 .nc avec 1 bloc"
        write(6,"(A)") "               : 3 .vtk"
        write(6,"(A)")
        write(6,"(A)") "Exemple :"
     !   write(6,"(A)") "   ./nblocs_1fichier -n nbprocs -d cheminAbsolu --iter_deb 1 --iter_fin n [--ini] [--fin] [--fluide] [--solide] [--iproc_ini 1] [--num_bloc 1]"
        stop
     case default
        write(6,"(A)") "option inconnue : ", trim(arg)
        stop 
     end select
     i = i+1
  end do

!!! recuperation du nombre de decoupe en x et en y
  if (nbproc /= 1) then
     open(unit=47,file=trim(path_donnees)//"/entrees/decoupage.dat")
     read(47,*)
     i = -1
     do while (i /= num_bloc)
        read(47,*) i, nb_proc_x, nb_proc_y
     end do
     !if (.not. fluide)  read(47,*) i, nb_proc_x, nb_proc_y
     close(47)
     if (nb_proc_x*nb_proc_y /= nbproc) then
      !  write(6,"(A,A,2I8)"), "Erreur : le nombre de procs donnés en argument ne correspond pas au fichier : ", trim(path_donnees)//"/entrees/decoupage.dat", nb_proc_x*nb_proc_y, nbproc
        stop 
     end if
  end if
   
  if (fluide) then
     path_donnees=trim(path_donnees)//"/visit/output"
  else
     path_donnees=trim(path_donnees)//"/visit/solide"
  end if
  
  if (ini .eqv. .true.) then
     suffixe="ini"
     call nbblocs_to_unfichier(trim(path_donnees), trim(suffixe), num_bloc, ip_ini, nb_proc_x, nb_proc_y, type_output)
  end if
  
  do it=itdeb,itfin ! boucle sur les iterations
     write(suffixe,'(i5.5)') it
     call nbblocs_to_unfichier(trim(path_donnees), trim(suffixe), num_bloc, ip_ini, nb_proc_x, nb_proc_y, type_output)
  end do
     
  if (fin .eqv. .true.) then
     suffixe="fin"
     call nbblocs_to_unfichier(trim(path_donnees), trim(suffixe), num_bloc, ip_ini, nb_proc_x, nb_proc_y, type_output)
  end if
  
  if (error .eqv. .true.) then
     suffixe="error"
     call nbblocs_to_unfichier(trim(path_donnees), trim(suffixe), num_bloc, ip_ini, nb_proc_x, nb_proc_y, type_output)
  end if
end program main
