module m_IO 
  use m_struct
  use variables_globales
  use parametres_globaux, only : FLUIDE, DIR_X, DIR_Y
  use m_communication_MPI
  use m_temps
  use EquationOfState
  use parametres_IO
  implicit none

  real*8, parameter :: tol_z_ = tol_z  ! 1.d-15

contains

  subroutine sorties_gnuplot(sb, dir)
    type (STR_SUPER_BLOC), pointer    :: sb
    integer              , intent(in) :: dir

    integer :: ld, lf, md, mf
    integer :: j, ib
    character(len=90) :: path 
    character(len=4) :: extension
    type (STR_BLOC), pointer :: b
    character(len=5) :: nb_bloc_c
    real*8, dimension(:,:), allocatable :: rho, p, u_x, z, y
    real*8, dimension(:,:), allocatable :: u_y, T, E, gamma, h, pi
    real*8, dimension(:,:,:), allocatable :: coord_nd

    integer :: pos

    path = trim(path_output)//"gnuplot/"
    extension = ".dat"

    do j = 1, sb%nb_membre 
       ib = sb%membre(j)
       b => acces_bloc(ib)%pointeur

       if (bloc_est_local (ib)) then
          if (sb%nb_membre_input /= 1) then
             write(nb_bloc_c,'(i2.2)') b%numero
             nb_bloc_c = '_'//trim(nb_bloc_c)
          else 
             nb_bloc_c=''
          end if

          call MPI_send_to_proc_ini(b, 2, b%coord_nd, coord_nd)
          call MPI_send_to_proc_ini(b, b%T  , T  , debug_MF)
          if (sb%nature == FLUIDE) then
             call MPI_send_to_proc_ini(b, 2, b%coord_nd, coord_nd)
             call MPI_send_to_proc_ini(b, b%rho, rho, debug_MF)
             call MPI_send_to_proc_ini(b, b%p  , p  , debug_MF)
             call MPI_send_to_proc_ini(b, b%u_x, u_x, debug_MF)
             call MPI_send_to_proc_ini(b, b%u_y, u_y, debug_MF)
             call MPI_send_to_proc_ini(b, b%z  , z  , debug_MF)
             call MPI_send_to_proc_ini(b, b%y  , y  , debug_MF)
             call MPI_send_to_proc_ini(b, b%h        , h      , debug_MF)
             call MPI_send_to_proc_ini(b, b%gamma_eq , gamma  , debug_MF)
             call MPI_send_to_proc_ini(b, b%pi_eq    , pi     , debug_MF)
          end if

          if (b%numproc == b%proc_ini) then

             if (debug_MF) then
                ld = b%ld_tot-nMF; lf = b%lf_tot+nMF
                md = b%md_tot-nMF; mf = b%mf_tot+nMF
             else
                ld = b%ld_tot; lf = b%lf_tot
                md = b%md_tot; mf = b%mf_tot
             end if
             
             select case (dir)
             case (DIR_X)
                pos = ceiling((b%mf_tot-b%md_tot+1)/2.d0)
                call write_gnuplot_scalar(ld, lf, coord_nd(1,ld-1:lf,pos), T(ld:lf,pos)    , path, nb_bloc_c, extension, "T")
                if (sb%nature == FLUIDE) then
                   call write_gnuplot_scalar(ld, lf, coord_nd(1,ld-1:lf,pos), rho(ld:lf,pos)  , path, nb_bloc_c, extension, "rho")
                   call write_gnuplot_scalar(ld, lf, coord_nd(1,ld-1:lf,pos), p(ld:lf,pos)    , path, nb_bloc_c, extension, "p")
                   call write_gnuplot_scalar(ld, lf, coord_nd(1,ld-1:lf,pos), u_x(ld:lf,pos)  , path, nb_bloc_c, extension, "u_x")
                   call write_gnuplot_scalar(ld, lf, coord_nd(1,ld-1:lf,pos), u_y(ld:lf,pos)  , path, nb_bloc_c, extension, "u_y")
                   call write_gnuplot_scalar(ld, lf, coord_nd(1,ld-1:lf,pos), z(ld:lf,pos)    , path, nb_bloc_c, extension, "z")
                   call write_gnuplot_scalar(ld, lf, coord_nd(1,ld-1:lf,pos), y(ld:lf,pos)    , path, nb_bloc_c, extension, "y")
                   call write_gnuplot_scalar(ld, lf, coord_nd(1,ld-1:lf,pos), h(ld:lf,pos)    , path, nb_bloc_c, extension, "enthalpie")
                   call write_gnuplot_scalar(ld, lf, coord_nd(1,ld-1:lf,pos), gamma(ld:lf,pos), path, nb_bloc_c, extension, "gamma")
                   call write_gnuplot_scalar(ld, lf, coord_nd(1,ld-1:lf,pos), pi(ld:lf,pos)   , path, nb_bloc_c, extension, "pi")
                end if
             case (DIR_Y)
                pos = ceiling((b%lf_tot-b%ld_tot+1)/2.d0)
                call write_gnuplot_scalar(md, mf, coord_nd(2,pos,md-1:mf), T(pos,md:mf)    , path, nb_bloc_c, extension, "T")
                if (sb%nature == FLUIDE) then
                   call write_gnuplot_scalar(md, mf, coord_nd(2,pos,md-1:mf), rho(pos,md:mf)  , path, nb_bloc_c, extension, "rho")
                   call write_gnuplot_scalar(md, mf, coord_nd(2,pos,md-1:mf), p(pos,md:mf)    , path, nb_bloc_c, extension, "p")
                   call write_gnuplot_scalar(md, mf, coord_nd(2,pos,md-1:mf), u_x(pos,md:mf)  , path, nb_bloc_c, extension, "u_x")
                   call write_gnuplot_scalar(md, mf, coord_nd(2,pos,md-1:mf), u_y(pos,md:mf)  , path, nb_bloc_c, extension, "u_y")
                   call write_gnuplot_scalar(md, mf, coord_nd(2,pos,md-1:mf), z(pos,md:mf)    , path, nb_bloc_c, extension, "z")
                   call write_gnuplot_scalar(md, mf, coord_nd(2,pos,md-1:mf), y(pos,md:mf)    , path, nb_bloc_c, extension, "y")
                   call write_gnuplot_scalar(md, mf, coord_nd(2,pos,md-1:mf), h(pos,md:mf)    , path, nb_bloc_c, extension, "enthalpie")
                   call write_gnuplot_scalar(md, mf, coord_nd(2,pos,md-1:mf), gamma(pos,md:mf), path, nb_bloc_c, extension, "gamma")
                   call write_gnuplot_scalar(md, mf, coord_nd(2,pos,md-1:mf), pi(pos,md:mf)   , path, nb_bloc_c, extension, "pi")
                end if
             end select
             
             deallocate(coord_nd, T)
             if (sb%nature == FLUIDE) then
                deallocate(rho, p, u_x, z, y)
                deallocate(u_y, gamma, h, pi)
             end if
          end if
       end if
    end do

  end subroutine sorties_gnuplot


  subroutine sortie_blasius(b)
    type (STR_BLOC), pointer :: b

    real*8, dimension(:), allocatable :: eta
    integer           :: iunit
    integer           :: l, m, ld, lf, md, mf
    character(len=90) :: fichier
    character(len=4)  :: nb_bloc_c
    real*8            :: x_pos, Re_x, u_inf, nu
    integer           :: l_pos
    logical           :: test_dir, test_pos

    test_pos = .false.


    ld = b%ld; lf = b%lf
    md = b%md; mf = b%mf

    write(nb_bloc_c,'(i4.4)') b%numero
    nb_bloc_c = trim(nb_bloc_c)

    x_pos = 0.9d0
    u_inf = 1.d2
    nu    = 2.d-4

    ! on trouve l'indice de la cellule qui correspond à x = x_pos
    do m = md,mf
      do l = ld,lf
        if(b%coord_nd(1,l-1,m)<= x_pos .and. x_pos <=b%coord_nd(1,l,m)) then
          l_pos = l
          test_pos = .true.
        end if
      enddo
    enddo

    if(test_pos) then
      ! on remplie le vecteur eta = y/x*(sqrt(Re_x))
      ! Re_x = u_inf * x / nu
      allocate(eta(md:mf))
      do m=md,mf
        Re_x = u_inf*b%centre_cl(1,l_pos,m)/nu
        eta(m) = b%centre_cl(2,l_pos,m)/b%centre_cl(1,l_pos,m) * sqrt(Re_x)
      enddo

      ! on teste si le dossier gnuplot existe deja
      inquire(FILE=trim(path_output)//'Blasius', exist=test_dir)
      if (.not. test_dir) then
         call system("mkdir "//trim(path_output)//'Blasius')
      end if



      iunit = iunitBlasius

      fichier = trim(path_output)//'/Blasius/blasius'//trim(nb_bloc_c)//'.dat'
      print *, 'Ecriture du fichier ', trim(fichier)
      open(unit=iunit,file=trim(fichier))
      do m = md, mf
         Re_x = u_inf*b%centre_cl(1,l_pos,m)/nu
         write(iunit,*) eta(m), b%u_x(l,m)/u_inf, b%u_y(l,m)* sqrt(Re_x)
      end do
      close(iunit)

      deallocate(eta)
    end if

  end subroutine sortie_blasius

!!!!!!!!!!!!!!!!!!!!!
!!!      Reprises
!!!!!!!!!!!!!!!!!!!!!

  subroutine protection_reprise(b, iter, nb_blocs)

    type (STR_BLOC), pointer    ::  b
    integer        , intent(in) :: iter, nb_blocs

    integer           :: iunit
    integer           :: ld, lf, md, mf
    character(len=5)  :: iter_c, nb_bloc_c 
    character(len=6)  :: type
    character(len=90) :: fichier
    
    real*8, dimension(:,:)  , allocatable :: T
    real*8, dimension(:,:,:), allocatable :: U, coord_nd

    ld = b%ld_tot; lf = b%lf_tot
    md = b%md_tot; mf = b%mf_tot

    iter_c = write_iter_c(iter, it_restart)
    if (nb_blocs /= 1) then
       write(nb_bloc_c,'(i2.2)') b%num_input
       nb_bloc_c = trim(nb_bloc_c)//'_'
    else 
       nb_bloc_c=''
    end if

    ! On envoie au proc les infos au bloc ini
    call MPI_send_to_proc_ini(b, 2, b%coord_nd, coord_nd)
    
    select case (b%nature)
    case (FLUIDE)
       iunit = iunitRepriseFld
       call MPI_send_to_proc_ini(b, b%neq, b%U, U, .true.)
       type = "fld"
    end select
    
    if (numproc == b%proc_ini) then
       fichier = trim(path_output)//'restart/restart_'//trim(nb_bloc_c)//trim(iter_c)
       open(iunit,file=trim(fichier)//'_maillage_'//trim(type)//'.dat', form="unformatted",status='unknown')
       write(iunit) iter, temps, iter_c, it_restart
       write(iunit) ld, lf, md, mf
       write(iunit) b%xin, b%xout, b%yin, b%yout
       write(iunit) ""  
       write(iunit) ""
       write(iunit) coord_nd( 1:2, ld-(nMF+1):lf+nMF, md-(nMF+1):mf+nMF)
       close(iunit)

       select case (b%nature)
       case (FLUIDE)
          open(iunit,file=trim(fichier)//'_cellules_fld.dat', form="unformatted",status='unknown')
          write(iunit) iter, temps
          write(iunit) b%neq, ld, lf, md, mf
          write(iunit) U(1:b%neq, ld-nMF:lf+nMF, md-nMF:mf+nMF)
          deallocate(U)
       end select
       close(iunit)

       deallocate(coord_nd)
    end if
  end subroutine protection_reprise

!!!!!!!!!!!!!!!!!!!!!
!!!      Sorties visit
!!!!!!!!!!!!!!!!!!!!!
  subroutine sorties_visit_nblocs(b, iter, iter_out)
    use m_maillage, only : distance

    type (STR_BLOC), pointer ::  b
    integer, intent(in) :: iter, iter_out

    integer :: iunit
    integer :: l, m, ld, lf, md, mf
!!$    real*8 :: lx, ly, l_moyen
    character(len=90) :: fichier
    character(len=6) :: iter_c, nb_bloc_c
    character(len=4) :: type
    real*8 :: div

    iter_c = write_iter_c(iter, iter_out)
    write(nb_bloc_c,'(i4.4)') b%numero
    nb_bloc_c = trim(nb_bloc_c)//'_'

    if (debug_MF) then
       ld = b%ld-nMF; lf = b%lf+nMF
       md = b%md-nMF; mf = b%mf+nMF
       !ld = b%ld_tot-1; lf = b%lf_tot+1
       !md = b%md_tot-1; mf = b%mf_tot+1
    else
       ld = b%ld; lf = b%lf
       md = b%md; mf = b%mf
    end if

    select case (b%nature)
    case (FLUIDE)
       iunit=iunitVisitFld
       fichier=trim(path_output)//'visit/output_'//trim(nb_bloc_c)//trim(iter_c)//'.vtk'
       type = 'Flow'
    end select
    
    if (numproc==0) print *, 'Ecriture du fichier ', trim(fichier)
    open(unit=iunit,file=fichier)
    write(iunit,'(1A26)') '# vtk DataFile Version 2.0'
    write(iunit,*) trim(type)
    write(iunit,*) 'ASCII'
    write(iunit,*) 'DATASET STRUCTURED_GRID'
    write(iunit,*) 'FIELD FieldData 1'
    write(iunit,*) 'TIME 1 1 double'
    write(iunit,*) temps
    write(iunit,*) 'DIMENSIONS', lf-ld+2,mf-md+2,1
    write(iunit,*) 'POINTS', (lf-ld+2)*(mf-md+2), 'double'
    do m=md-1,mf
       do l=ld-1,lf
          write(iunit,*) b%coord_nd(1,l,m), b%coord_nd(2,l,m), 0.0d0
       end do
    end do
    write(iunit,*) 'CELL_DATA', (lf-ld+1)*(mf-md+1)
    select case (b%nature)
    case (FLUIDE)
       call write_vtk_scalar(iunit, ld, lf, md, mf, b%u_x(ld:lf,md:mf), "u")
       call write_vtk_scalar(iunit, ld, lf, md, mf, b%u_y(ld:lf,md:mf), "v")
       call write_vtk_scalar(iunit, ld, lf, md, mf, b%rho(ld:lf,md:mf), "rho")
       call write_vtk_scalar(iunit, ld, lf, md, mf, b%z(ld:lf,md:mf)  , "z")
       call write_vtk_scalar(iunit, ld, lf, md, mf, b%p(ld:lf,md:mf)  , "p")
       call write_vtk_scalar(iunit, ld, lf, md, mf, b%c(ld:lf,md:mf)  , "c")

       write(iunit,*) 'SCALARS eps double'
       write(iunit,*) 'LOOKUP_TABLE default'
       do m=md,mf
          do l=ld,lf
             write(iunit,*) b%U(irhoE,l,m)/b%rho(l,m) - 0.5d0*(b%u_x(l,m)**2+b%u_y(l,m)**2)
          end do
       end do

       call write_vtk_scalar(iunit, ld, lf, md, mf, b%U(irho1z1,ld:lf,md:mf), "rho1_z")
      ! call write_vtk_scalar(iunit, ld, lf, md, mf, b%rho1(ld:lf,md:mf)     , "rho1")
      ! call write_vtk_scalar(iunit, ld, lf, md, mf, b%rho2(ld:lf,md:mf)     , "rho2")
      ! call write_vtk_scalar(iunit, ld, lf, md, mf, b%y(ld:lf,md:mf)        , "y")

       if(b%correction_son>0)then
           write(iunit,*) 'SCALARS flag_c double'
           write(iunit,*) 'LOOKUP_TABLE default'
           do m=md,mf
              do l=ld,lf
                 write(iunit,*) b%debug_c(l,m)
              end do
           end do
       end if

      ! write(iunit,*) 'SCALARS Mach double'
      ! write(iunit,*) 'LOOKUP_TABLE default'
      ! do m=md,mf
      !    do l=ld,lf
      !       write(iunit,*) sqrt(b%u_x(l,m)**2.d0 + b%u_y(l,m)**2.d0)/b%c(l,m)
      !    enddo
      ! enddo
      ! write(iunit,*) 'VECTORS U double'
      ! do m=md,mf
      !    do l=ld,lf
      !       write(iunit,*) b%u_x(l,m), b%u_y(l,m), 0.d0
      !    enddo
      ! enddo
    
      ! if (nbprocs /= 1) then 
      !    b%zeros(ld:lf,md:mf) = numproc
      !    call write_vtk_scalar(iunit, ld, lf, md, mf, b%zeros(ld:lf,md:mf) , "proc")
      !    b%zeros(ld:lf,md:mf) = 0.d0
      ! end if
       
       if (b%maillage%move) then
          call write_vtk_scalar(iunit, ld, lf, md, mf, b%dX_dt(ld:lf,md:mf), "Xt")
          call write_vtk_scalar(iunit, ld, lf, md, mf, b%dY_dt(ld:lf,md:mf), "Yt")
       end if     
       
       select case(fluide1%type)
       case (EOS)
          if (b%visqueux) then
             call write_vtk_scalar(iunit, ld, lf, md, mf, b%T(ld:lf,md:mf), "T")
          end if
       case (CHIMIE_FIGEE)
          call write_vtk_scalar(iunit, ld, lf, md, mf, b%T(ld:lf,md:mf), "T")
          do l = 1, fluide1%nb_esp
             call write_vtk_scalar(iunit, ld, lf, md, mf, b%c_esp(l,ld:lf,md:mf), "C_"//trim(fluide1%especes(l)%nom))
          end do
       end select

!!$    if (b%visqueux) then
!!$       write(iunit,*) 'SCALARS Re double'
!!$       write(iunit,*) 'LOOKUP_TABLE default'
!!$       do m=md,mf
!!$          do l=ld,lf
!!$             lx = .5d0 * ( distance(coord_nd(:,l,m),coord_nd(:,l,m-1)) + &
!!$                  distance(coord_nd(:,l-1,m),coord_nd(:,l-1,m-1)) )
!!$             ly = .5d0 * ( distance(coord_nd(:,l,m),coord_nd(:,l-1,m)) + &
!!$                  distance(coord_nd(:,l,m-1),coord_nd(:,l-1,m-1)) )
!!$             l_moyen = 2.d0 * ly * lx / (lx + ly) ! moyenne harmonique
!!$             write(iunit,*) l_moyen*rho(l,m)*sqrt(u_x(l,m)**2.d0 + u_y(l,m)**2.d0)/mu(l,m)
!!$          end do
!!$       end do
!!$    end if

       if (b%maillage%type == "bosse_dellacherie") then
          write(iunit,*) 'SCALARS p_normalized double'
          write(iunit,*) 'LOOKUP_TABLE default'
          do m=md,mf
             do l=ld,lf
                write(iunit,*) (b%p(l,m)-1.d5)/1.d5 
                !write(iunit,*) (b%p(l,m)-minval(b%p(ld:lf,md:mf))) / &
                !     (maxval(b%p(ld:lf,md:mf))-minval(b%p(ld:lf,md:mf)))
             end do
          end do
       end if

     !  write(iunit,*) 'SCALARS div_U double'
     !  write(iunit,*) 'LOOKUP_TABLE default'
     !  do m=md,mf
     !     do l=ld,lf
     !        div = 0.d0
     !        div = div + .5d0*(b%u_x(l,m)+b%u_x(l+1,m))*b%dl_l(l,m)
     !        div = div - .5d0*(b%u_x(l,m)+b%u_x(l-1,m))*b%dl_l(l-1,m)
     !        div = div + .5d0*(b%u_y(l,m)+b%u_y(l,m+1))*b%dl_m(l,m)
     !        div = div - .5d0*(b%u_y(l,m)+b%u_y(l,m-1))*b%dl_m(l,m-1)
     !        write(iunit,*) div
     !     end do
     !  end do

       !call write_vtk_scalar(iunit, ld, lf, md, mf, b%dt_a(ld:lf,md:mf) , "dt_a")
       !call write_vtk_scalar(iunit, ld, lf, md, mf, b%dt_vv(ld:lf,md:mf), "dt_vv")
       !call write_vtk_scalar(iunit, ld, lf, md, mf, b%dt_vT(ld:lf,md:mf), "dt_vT")
       !call write_vtk_scalar(iunit, ld, lf, md, mf, b%dt_t(ld:lf,md:mf) , "dt_t")

    end select
    close(iunit)

  end subroutine sorties_visit_nblocs

  subroutine write_vtk_scalar(iunit, ld, lf, md, mf, tab, nom)
    integer         , intent(in) :: iunit
    integer         , intent(in) :: ld, lf, md, mf
    real*8          , intent(in) :: tab(ld:lf,md:mf)
    character(len=*), intent(in) :: nom
    
    write(iunit,*) 'SCALARS '//trim(nom)//' double'
    write(iunit,*) 'LOOKUP_TABLE default'
    write(iunit,*) tab 

  end subroutine write_vtk_scalar

  subroutine write_gnuplot_scalar(ld, lf, coord_nd, tab, path, nb_bloc_c, extension, nom)
    integer         , intent(in) :: ld, lf
    real*8          , intent(in) :: coord_nd(ld-1:lf)
    real*8          , intent(in) :: tab(ld:lf)
    character(len=*), intent(in) :: path, nb_bloc_c, extension, nom
    
    integer           :: l, iunit
    character(len=90) :: fichier

    iunit = iunitGnuplot

    fichier = trim(path)//trim(nom)//trim(nb_bloc_c)//extension
    print *, 'Ecriture du fichier ', trim(fichier)
    open(unit=iunit,file=trim(fichier))
    do l = ld, lf
       write(iunit,*) .5d0*(coord_nd(l-1)+coord_nd(l)), tab(l)
    end do
    close(iunit)
    
  end subroutine write_gnuplot_scalar
  
  function write_iter_c(iter, freq)
    integer :: iter, freq
    character(len=5) :: write_iter_c

    if (iter == -1) then
       write_iter_c = "fin"
    else if (iter == -2) then
       write_iter_c = "error"
    else if (iter == 0) then
       write_iter_c = "ini"  
    else if (iter == 1 .and. freq /= 1) then
       write_iter_c = "00000" 
    else if (mod(iter, freq) == 0) then
       write(write_iter_c,'(i5.5)') iter/freq
    else
       write(write_iter_c,'(i5.5)') iter
    end if
  end function write_iter_c
end module m_IO
