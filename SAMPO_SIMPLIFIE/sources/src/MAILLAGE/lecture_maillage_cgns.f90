module lecture_maillage_cgns
  use variables_globales
  use m_struct
  use m_MPI
  use parametres_globaux, only : nMF

  implicit none

contains

  subroutine maillage_cgns(b, filename)
    type (STR_BLOC), pointer :: b
    character(len=90), intent(in) :: filename
    
    print*
    print*, " ||---- : type de fichier non supporte"
    call arret_code
  end subroutine maillage_cgns

end module lecture_maillage_cgns

