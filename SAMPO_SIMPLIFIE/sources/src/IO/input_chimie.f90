module input_chimie
  use m_struct
  use parametres_IO
  use variables_globales
  use m_MPI

  implicit none
  
  public read_chimie, read_thermo_database
  private
contains
  
  subroutine read_chimie(fluidek, espece, element, reaction)

    type(STR_fluide)  , intent(inout) :: fluidek
    type(STR_ESPECE)  , pointer       :: espece(:)
    type(STR_ELEMENT) , pointer       :: element(:)
    type(STR_REACTION), pointer       :: reaction(:)

    ! variables locales
    integer                        :: i, j, k, p
    integer                        :: nb_esp, nb_ele
    logical                        :: comparaison
    integer          , allocatable :: mat_mult(:,:)
    character(len=1) , allocatable :: tmp(:)
    type(STR_ELEMENT), allocatable :: element_decode(:)

    open(unit=iunitChimie,file=trim(path_input)//'chimie.dat')
    read(iunitChimie,*) 
    read(iunitChimie,*) 
    ! lecture du nombre d'element
    read(iunitChimie,*) ! Nombre d'éléments
    read(iunitChimie,*) nb_ele

    ! allocation
    allocate(element(nb_ele))    ! elements chimiques

    read(iunitChimie,*) ! Liste des elements
    do i = 1, nb_ele
       read(iunitChimie,*) element(i)%nom
    end do

    ! on lit les réactions chimiques et on en déduit la liste d'espèces
    call read_reactions_chimiques(reaction, espece)
    
    close(iunitChimie)

!!! calcul de la matrice de multiplicite
   ! allocate(fluide1%mat_multiplicite(nb_esp,nb_ele))
   ! mat_mult => fluide1%mat_multiplicite
    
    call decode_elements_from_especes(espece, element_decode, mat_mult)
    
    ! compararaion avec les elements lus
    call verif_elem(element_decode, element, comparaison)
    if (comparaison .eqv. .true.) then
       element = element_decode
    else
       print*, "incoherence entre les elements et les espèces du fichier chimie.dat"
       call arret_code
    end if

    fluidek%nb_esp = size(espece)
    nb_esp = size(espece)
    fluidek%nb_ele = size(element)

    ! creation et initialisation des structure liste_espece/element et multiplicite
    call decode_liste_elements(espece, mat_mult)
    call decode_liste_especes(element, mat_mult)
    
    allocate(fluidek%mat_multiplicite(nb_esp,nb_ele))
    fluidek%mat_multiplicite = mat_mult

    deallocate(mat_mult, element_decode)
  end subroutine read_chimie

  subroutine verif_elem(element_decode, element_lu, comparaison)
    ! pas propre du tout les allocatable et pointer ...
    type(STR_ELEMENT) , allocatable   :: element_decode(:)
    type(STR_ELEMENT) , pointer       ::  element_lu(:)
    logical           , intent(out)   :: comparaison
    
    integer :: k, nb_ele, l
    logical :: element_trouve

    comparaison = .false.
    
    if (size(element_decode) /= size(element_decode)) then
       comparaison = .false.
       return
    end if
    
    nb_ele = size(element_decode)
    
    do k = 1, nb_ele
       element_trouve = .false.
       do l = 1, nb_ele
          if (element_decode(l)%nom == element_lu(k)%nom) then
             element_trouve = .true.
          end if
       end do
       if ( element_trouve .eqv. .false. ) then
          return
       end if
    end do
    comparaison = .true.

  end subroutine verif_elem
    

  subroutine read_reactions_chimiques(reaction, espece)
    type(STR_REACTION), pointer       :: reaction(:)
    type(STR_ESPECE)  , pointer       :: espece(:)

    ! variables locales
    integer                       :: i, j, k, p, r, ir
    integer                       :: nb_reac, nb_espinreac, nb_esp
    integer                       :: sgn, compt
    character(len=8)              :: esp
    character(len=20)             :: ligne
    character(len=8), allocatable :: tabesp(:,:)
    character(len=8), allocatable :: test(:)
    character(len=1), allocatable :: tmp(:)
    logical                       :: new_esp
    
    p=0
    compt=0
    r=1
    
    read(iunitChimie,*) ! Nombre de réactions
    read(iunitChimie,*) nb_reac
     
    ! allocation
    allocate(reaction(nb_reac))  ! reactions chimiques
    allocate(tabesp(nb_reac,40)) ! 40 esp maximun par reaction (...)
    
    read(iunitChimie,*) ! Reactions
    allocate(tmp(20)) ! pour lire caractere par caractere la ligne (bof, bof, bof...)
    do i = 1, nb_reac
       read(iunitChimie,*) ligne
       do j=1,20
          read(ligne,'(20a1)') tmp(1), tmp(2), tmp(3), tmp(4), tmp(5), tmp(6), tmp(7), tmp(8), tmp(9), tmp(10), &
               &               tmp(11),tmp(12),tmp(13),tmp(14),tmp(15),tmp(16),tmp(17),tmp(18),tmp(19),tmp(20)
       end do
       j=1
       p=0
       r=1
       do while(tmp(j)/=' ')
          j=j+1
          if(tmp(j)=='+'.and.p==0)then
             r=r+1
          else if(tmp(j)=='=')then
             p=p+1
          else if(tmp(j)=='+'.and.p/=0)then
             p=p+1
          end if
       end do
       nb_espinreac = r+p
       reaction(i)%nesp = nb_espinreac

       compt=compt+r+p

       j=1
       k=1
       sgn=-1
       allocate(reaction(i)%nu(nb_espinreac)) ! coefficients stoechio des espèces dans la réaction

       do while(tmp(j)/=' ')
          read(tmp(j),*) r
          reaction(i)%nu(k)=sgn*r
          j=j+1
          esp=tmp(j)
          j=j+1
          do while(tmp(j)/='+'.and.tmp(j)/='='.and.tmp(j)/=' ')
             esp=trim(esp)//tmp(j)
             j=j+1
          end do
          if(tmp(j)=='=')then
             sgn=1
          end if
          tabesp(i,k)=esp
          j=j+1
          k=k+1
       end do
       reaction(i)%D_nu=0
       do j = 1, nb_espinreac
          reaction(i)%D_nu = reaction(i)%D_nu + reaction(i)%nu(j)
       end do
    end do
    deallocate(tmp)

    allocate(test(compt))
    test(1)=tabesp(1,1)
    nb_esp = 1
    do i = 1, nb_reac  ! boucle sur les reactions
       do j = 1, reaction(i)%nesp  ! boucle sur les especes presentes dans la reaction
          new_esp=.true.
          do k = 1, nb_esp
             if ( trim(test(k)) == trim(tabesp(i,j)) ) then
                new_esp=.false.
             end if
          end do
          if ( new_esp .eqv. .true. ) then
             nb_esp=nb_esp+1
             test(nb_esp) = tabesp(i,j)
          end if
       end do
    end do

    allocate(espece(nb_esp))
    !!! pour mettre le premier en dernier ???
    do i = 1, nb_esp-1
       espece(i)%nom = test(i+1)
    end do
    espece(nb_esp)%nom = test(1)
    
    do r = 1, nb_reac
       allocate(reaction(r)%ist(reaction(r)%nesp))
       allocate(reaction(r)%ist_inv(nb_esp))

       do ir = 1, reaction(r)%nesp
          do i = 1, nb_esp
             if(tabesp(r,ir) == espece(i)%nom)then
                reaction(r)%ist(ir)=i
                reaction(r)%ist_inv(i)=ir
             end if
          end do
       end do
    end do
    
    deallocate(tabesp)
  end subroutine read_reactions_chimiques


!!! Lecture des paramètres de chimie dans la base de données
  subroutine read_thermo_database(espece, element)
    use parametres_fluide, only : R_gaz
    use outils_fluide, only : recup_ind_esp

    type(STR_ESPECE) , intent(inout) :: espece(:)
    type(STR_ELEMENT), intent(inout) :: element(:)

    integer :: i, j, k, p
    integer :: nb_esp, nb_ele
    integer :: nb_poly_glenn, nb_coeff_glenn
    character(len=8) :: nom_lu

    nb_esp = size(espece)
    nb_ele = size(element)

    ! a mettre en arguments pour pouvoir choisir plusieurs bases thermo (UMUSCL 7 ici, mais UMUSCL 9 ou UMUSCL 11 possibles)
    nb_poly_glenn = 3
    nb_coeff_glenn = 7

    do i = 1, nb_esp
       allocate(espece(i)%coeff_glenn(nb_coeff_glenn, nb_poly_glenn))
       allocate(espece(i)%coeff_C0(nb_poly_glenn))
       allocate(espece(i)%coeff_C1(nb_poly_glenn))
    end do

    open(unit=9,file=trim(path_input)//'database.dat')
    read(9,*)
    read(9,*)
    read(9,*)
    read(9,*)
    read(9,*)
    read(9,*)
    read(9,*)
    i = 0
    do while ( i /= nb_esp )
       i = i+1

       read(9,*) nom_lu

       j = -1
       do k = 1, nb_esp
          if (trim(espece(k)%nom) == trim(nom_lu)) then
             j = k
             exit
          end if
       end do
       if (j == -1) then
          read(9,*)
          do p = 1, nb_poly_glenn
             read(9,*)
          end do
          cycle ! l'espece lue n'est pas dans la liste d'especes du cas
       end if

       read(9,*) espece(j)%masse, espece(j)%h0, espece(j)%s0, espece(j)%diatomique
       do p = 1, nb_poly_glenn
          read(9,*) espece(j)%coeff_glenn(:,p), espece(j)%coeff_C0(p), espece(j)%coeff_C1(p)
       end do
    end do
    close(9)

!!! Calcul et ajustage des masses
    call calcul_masses_elements(element, espece)
    call ajustage_masses_especes(espece, element)

!!! Calcul du R par especes R = R_gaz/masse
    do i = 1, nb_esp
       espece(i)%R = R_gaz/espece(i)%masse
    end do

    ! on initialise comme on peut pour que ca passe, mais il faudrait lire les valeurs
    do i = 1, nb_esp
       espece(i)%mu = 1.d0
    end do

  end subroutine read_thermo_database

  
!!! On retrouve les éléments présents dans une liste d'espèces donnée
  subroutine decode_elements_from_especes(especes, elements, mat_mult)
    type(STR_ESPECE)              , intent(in)  :: especes(:)
    type(STR_ELEMENT), allocatable, intent(out) :: elements(:)
    integer          , allocatable, intent(out) :: mat_mult(:,:)

    integer :: nb_elements, nb_especes
    integer :: i, k, ind_elmt, nbEM
    type(STR_ELEMENT), allocatable :: ele_pdoublons(:) ! elements plus les doublons
    integer          , allocatable :: mat_mult_pdoublons(:,:)

    nb_especes = size(especes)
    ! on fait une surestimation du nombre d elements pour pouvoir lire les doublons
    nb_elements = 2*nb_especes

    allocate(ele_pdoublons(nb_elements))
    allocate(mat_mult_pdoublons(nb_especes, nb_elements))
    mat_mult_pdoublons = 0

    ind_elmt=1
    ! boucle sur les especes
    do i = 1, nb_especes
       ! on lit l'espece en partant de la fin
       k = len_trim(adjustl(especes(i)%nom))
       lecture : do while (k>0)
          select case (especes(i)%nom(k:k))
          case("+", "-") ! cas des especes ionisees
             ele_pdoublons(ind_elmt)%nom = "e-"

             nbEM = 1 ! nombre d'electrons
             if (especes(i)%nom(k:k) == '+') nbEM = -1

             k = k-1

             mat_mult_pdoublons(i,ind_elmt) = mat_mult_pdoublons(i,ind_elmt) + nbEM
             ind_elmt = ind_elmt + 1
             if (k==1 .and. especes(i)%nom(k:k+1)=="e-") exit lecture ! cas de e-
          case("2", "3", "4", "5", "6", "7", "8", "9")
             read(especes(i)%nom(k:k),'(i1)') mat_mult_pdoublons(i,ind_elmt)
             k = k-1
          case default
             if (ichar(especes(i)%nom(k:k))>=97 .and. ichar(especes(i)%nom(k:k))<=122) then ! si on a une minuscule
                ele_pdoublons(ind_elmt)%nom = especes(i)%nom(k-1:k)
                k=k-1
             else
                ele_pdoublons(ind_elmt)%nom = especes(i)%nom(k:k)
             end if
             k = k-1
             if (mat_mult_pdoublons(i,ind_elmt) == 0) mat_mult_pdoublons(i,ind_elmt) = 1
             ind_elmt = ind_elmt + 1
          end select
       end do lecture
    end do

    nb_elements = ind_elmt-1
!!! gestion des doublons
    call enleve_doublons(ele_pdoublons(1:nb_elements), elements, mat_mult_pdoublons, mat_mult)

    deallocate(ele_pdoublons, mat_mult_pdoublons)

  end subroutine decode_elements_from_especes

!!! On supprime les doublons sur les elements que l'on a trouves
  subroutine enleve_doublons(elements_doublons, elements, mat_doublons, mat)
    implicit none

    type(STR_ELEMENT)             , intent(in)    :: elements_doublons(:)
    type(STR_ELEMENT), allocatable, intent(out)   :: elements(:)
    integer                       , intent(inout) :: mat_doublons(:,:)
    integer          , allocatable, intent(out)   :: mat(:,:)

    type(STR_ELEMENT) :: elements_tmp(size(elements_doublons))
    integer           :: mat_mult_tmp(ubound(mat_doublons,1),ubound(mat_doublons,2))

    integer :: k    ! indice pour la liste d elements avec les doublons
    integer :: nb_el
    integer :: l

    nb_el = 1

    elements_tmp(1) = elements_doublons(1)
    mat_mult_tmp(:,:) = 0
    mat_mult_tmp(:,1) = mat_doublons(:,1)
    outer : do k=2,size(elements_doublons) ! boucle sur la liste d'elements avec les doublons
       do l = 1, nb_el
          if (elements_tmp(l)%nom == elements_doublons(k)%nom) then
             ! on doit repercuter le changement de num de l'element dans la matrice de multiplicite
             ! on ajouter les deux colonnes correspondant au meme element (et on met a 0 l'ancienne colonne)
             mat_mult_tmp(:,l) = mat_mult_tmp(:,l) + mat_doublons(:,k)
             mat_doublons(:,k) = 0
             cycle outer
          end if
       end do
       !if (any(element_tmp%nom == elements_doublons(k)%nom )) cycle
       nb_el = nb_el + 1
       elements_tmp(nb_el) = elements_doublons(k)

       ! on doit repercuter le changement de num de l'element dans la matrice de multiplicite
       mat_mult_tmp(:,nb_el) = mat_mult_tmp(:,nb_el) + mat_doublons(:,k)
       mat_doublons(:,k) = 0
    end do outer

    allocate(elements(nb_el))
    do k=1,nb_el  ! boucle sur la liste d'elements sans les doublons
       elements(k)%nom = elements_tmp(k)%nom
    end do
    allocate(mat(1:ubound(mat_doublons,1),1:nb_el))
    mat(:,1:nb_el) = mat_mult_tmp(:,1:nb_el)

  end subroutine enleve_doublons


!!!--------------------------------------------------------------------------!!!
!!!  Routines permettant de retrouver                                        !!!
!!!    la liste des especes presentes dans un element (avec la multiplicite) !!!
!!!    la liste des elements contenus dans une espece (avec la multiplicite) !!!
!!!--------------------------------------------------------------------------!!!
  subroutine decode_liste_elements(especes, mat_multiplicite)
    implicit none

    type(STR_ESPECE) , dimension(:)  , intent(inout) :: especes
    integer          , dimension(:,:), intent(in)    :: mat_multiplicite

    integer :: i
    integer :: nb_elements_loc, k, k_loc

    ! creation des listes d'especes et d'elements (et les multiplicites associees)
    do i = 1, size(especes)
       ! nombre d'entrees non nulles sur la ligne
       nb_elements_loc =size(mat_multiplicite(i,:)) - count(mat_multiplicite(i,:)==0)
       allocate(especes(i)%multiplicite(nb_elements_loc))
       especes(i)%multiplicite(:)=pack(mat_multiplicite(i,:), mat_multiplicite(i,:) /=0)
       allocate(especes(i)%liste_elements(nb_elements_loc))
       k_loc = 1
       do k = 1, size(mat_multiplicite(i,:))
          if (mat_multiplicite(i,k) /= 0) then
             especes(i)%liste_elements(k_loc) = k
             k_loc = k_loc+1
          end if
       end do
    end do
  end subroutine decode_liste_elements

  subroutine decode_liste_especes(elements, mat_multiplicite)
    implicit none

    type(STR_ELEMENT), dimension(:)  , intent(inout) :: elements
    integer          , dimension(:,:), intent(in)    :: mat_multiplicite

    integer :: k
    integer :: nb_especes_loc, i, i_loc

    do k = 1, size(elements)
       nb_especes_loc =size(mat_multiplicite(:,k)) - count(mat_multiplicite(:,k)==0)
       allocate(elements(k)%multiplicite(nb_especes_loc))
       elements(k)%multiplicite(:)=pack(mat_multiplicite(:,k), mat_multiplicite(:,k) /=0)
       allocate(elements(k)%liste_especes(nb_especes_loc))
       i_loc = 1
       do i =1, size(mat_multiplicite(:,k))
          if (mat_multiplicite(i,k) /= 0) then
             elements(k)%liste_especes(i_loc) = i
             i_loc = i_loc+1
          end if
       end do
    end do
  end subroutine decode_liste_especes


  subroutine calcul_masses_elements(elements, especes)
    implicit none

    type(STR_ESPECE) , intent(in)  :: especes(:)
    type(STR_ELEMENT), intent(out) :: elements(:)

    integer :: i, k
    logical :: ok

    do k = 1, size(elements)
       ok = .false.
       ! boucle sur les especes possedant cet element
       do i = 1, size(elements(k)%liste_especes(:))
          ! si l'espece a le meme nom que l'element on recupere la masse
          if (elements(k)%nom == especes(elements(k)%liste_especes(i))%nom) then
             elements(k)%masse = especes(elements(k)%liste_especes(i))%masse
             ok = .true.
          end if
       end do

       if (.not. ok) then
          print*, "Erreur : on ne peut pas calculer la masse de l'element :", elements(k)%nom, "car il n'y pas l'espece correspondante"
          call arret_code()
       end if
    end do

  end subroutine calcul_masses_elements

  subroutine ajustage_masses_especes(especes, elements)
    implicit none

    type(STR_ELEMENT), intent(in)  :: elements(:)
    type(STR_ESPECE) , intent(out) :: especes(:)

    integer :: i,k

    do i = 1, size(especes)
       especes(i)%masse = 0.d0
       ! boucle sur les elements contenus dans l'espece consideree
       do k=1, size(especes(i)%liste_elements(:)) 
          especes(i)%masse = especes(i)%masse + especes(i)%multiplicite(k) * elements(especes(i)%liste_elements(k))%masse
       end do
    end do
  end subroutine ajustage_masses_especes

end module input_chimie
