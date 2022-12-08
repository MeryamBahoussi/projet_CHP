module m_solveurs
  
  use solveur_MKL
  use solveur_balayage
  use solveur_petsc, only : init_petsc, init_solveur_global, RL_petsc_global, petsc, mkl
  
end module m_solveurs
