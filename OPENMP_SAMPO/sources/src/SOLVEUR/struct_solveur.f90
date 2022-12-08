module m_struct_solveur
#ifdef PETSC_ACTIF
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
  use petscvec
  use petscmat
  use petscksp
#endif

  implicit none

  !* STRUCTURE POUR LES SOLVEURS PETSC
  type STR_PETSC
     real*8               :: tol                 ! tolerance
     integer              :: kmax                ! nb d'iter max pour l'inversion
#ifdef PETSC_ACTIF
     Vec                     vec_x               ! solution vector
     Vec                     vec_b               ! right-hand-side vector
     Mat                     mat_A               ! matrix that defines linear system  
     KSP                     ksp                 ! Krylov subspace method context (linear solver context)
     PC                      ksp_pc              ! preconditioner context
#endif
  end type STR_PETSC
  
  
  !* STRUCTURE POUR LE SOLVEUR MKL (stockage creux CSR)
  type STR_CSR
     integer :: nb_rows, nNonZeros
     integer, dimension(:), pointer :: columns, rowIndex
     real*8 , dimension(:), pointer :: Values
  end type STR_CSR
  
end module m_struct_solveur
