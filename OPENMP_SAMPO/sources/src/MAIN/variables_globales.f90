module variables_globales
  implicit none
  
  integer :: iter_ini, itermax
  integer :: it_out, it_restart, it_debug
  logical :: film, debug_MF, output_1d

  logical :: CFLTransportOK

  ! g : acceleration de pesanteur
  real*8, dimension(2) :: acc_pesanteur = (/0.d0, -9.81d0/) 
  
  character(len=30) :: path_input = "entrees/"
  character(len=30) :: path_output= "./" !"output/"

  ! pour annuler la correction d'entropie dans la couche limite (schema de Roe)
  real*8 :: Htot_infini_amont
end module variables_globales


module parametres_globaux
  implicit none
  
  ! nature du super_bloc
  integer, parameter :: FLUIDE = 0
  integer, parameter :: SOLIDE = 1
  ! direction du maillage
  integer, parameter :: DIR_X = 1
  integer, parameter :: DIR_Y = 2
  ! methode de resolution du systeme lineaire
  integer, parameter :: MKL   = 0
  integer, parameter :: GS    = 1
  integer, parameter :: PETSc = 2
  ! geometrie
  integer, parameter :: PLAN2D = 1
  ! nombre de mailles fictives
  integer, parameter :: nMF = 2
end module parametres_globaux

module variables_debug
  implicit none
  
  logical :: debug_volume_bulle=.false.
  logical :: debug_residu=.false.
  logical :: debug_conservation=.true.
  logical :: energies=.true.
  logical :: saut_pression=.true.
  logical :: amplitude=.false.

 ! logical, parameter :: cmp_manu=.true. !.false.
end module variables_debug

module variables_sorties
  implicit none

  logical :: Blasius = .false.

end module variables_sorties


module a_mettre_en_donnees
  use parametres_fluide
  implicit none
  
  integer :: implicite = implicite_Up !implicite_total
  
end module a_mettre_en_donnees
