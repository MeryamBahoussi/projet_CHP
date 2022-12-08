module parametres_fluide
  implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    indices des differents tableaux   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! U = ( rho y, rho (1-y), rho u, rho v, rho E, z)
!!! U = ( rho y c_1, rho y c_2, ..., rho y c_ne, rho (1-y), rho u, rho v, rho E, z)
  integer, save :: irho1z1
  integer, save :: irho2z2
  integer, save :: irhou
  integer, save :: irhov
  integer, save :: irhoE
  integer, save :: iz

!!! Ureco = (rho1, rho2, rho1eps1, rho2eps2, u, v, z)
  integer, save :: irho1
  integer, save :: irho2
  integer, save :: irho1eps1
  integer, save :: irho2eps2
  integer, save :: iux
  integer, save :: iuy
  integer, save :: izz

!!! V = ( 1/rho, y, u, v, E ) pour la partie Lagrange
!!! V = ( 1/rho, y c1, y c_2, ..., y c_ne, u, c, E)
  integer, save :: ivol
  integer, save :: iy
  integer, save :: iu
  integer, save :: iv
  integer, save :: iE

!!! X = ( u, v, p ) pour l'implicite
  integer, parameter :: iXu = 1
  integer, parameter :: iXv = 2
  integer, parameter :: iXp = 3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! parametres & booleens pour le fluide !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! methodes de resolution
  integer, parameter :: EULER_DIRECT    = 2
  ! schemas pour la partie acoustique
  integer, parameter :: GALLICE  = 20
  integer, parameter :: HLL      = 22
  integer, parameter :: ROE      = 4
  ! schema en temps pour
  integer, parameter :: EULER_EXPLICITE = 1
  ! méthode de calcul du gradient
  integer, parameter :: WLSQ5 = 1 ! Weighted Least SQuare 5 points
  integer, parameter :: WLSQ9 = 2 ! Weighted Least SQuare 9 points
  integer, parameter :: GG    = 3 ! Green-Gauss cell based
  ! méthode pour résoudre le calcul ud gradient avec les moindres carrés (WLSQ)
  integer, parameter :: Psinv = 1 ! inverse direct de la pseudoinverse
  integer, parameter :: QR    = 2 ! decomposition QR


  ! correction EUCCLHYD
  integer, parameter :: hybrid       = 1
  integer, parameter :: conservative = 2

  integer, parameter :: implicite_Up=1, implicite_total=2
  logical :: cdp_Riemann=.true.
  logical :: cdp_pression=.false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! type de fermeture thermo pour le fluide !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: EOS = 1
  integer, parameter :: CHIMIE_FIGEE = 2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!             Lois d'Etat              !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: GAZ_PARFAIT = 0
  integer, parameter :: STIFFENED_GAS = 1
  integer, parameter :: MIE_GRUNEISEN = 2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Lois de viscosite et de conductivite !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: CONSTANTE  = 0
  integer, parameter :: SUTHERLAND = 1
  integer, parameter :: WILKE      = 2
  integer, parameter :: PRANDTL    = 3
  integer, parameter :: LEWIS      = 4
  integer, parameter :: MPP        = 5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!        Monophasique ou Allaire       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: MONOPHASIQUE = 0
  integer, parameter :: ALLAIRE      = 1
  integer, save :: modele_diphasique

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!       conditions aux limites         !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: CL_FLUX_NUL            = 1
  integer, parameter :: CL_SYMETRIE            = 2
  integer, parameter :: CL_INTERBLOC           = 3
  integer, parameter :: CL_INTERPROC           = 4
  integer, parameter :: CL_ENTREE_SUPERSONIQUE = 5
  integer, parameter :: CL_ENTREE_SUBSONIQUE   = 6
  integer, parameter :: CL_ENTREE_SUBSONIQUE_t = 7
  integer, parameter :: CL_SORTIE_SUBSONIQUE   = 8
  integer, parameter :: CL_PAROI               = 9
  integer, parameter :: CL_PAROI_ADIABATIQUE   = 10
  integer, parameter :: CL_INJECTION           = 11
  integer, parameter :: CL_PERIODIQUE          = 12

  !integer, parameter :: CL_INTERBLOCvs = 30

  integer, parameter :: CL_INTERBLOCvsFLUX_NUL       = 301
  integer, parameter :: CL_INTERBLOCvsSYMETRIE       = 302
  integer, parameter :: CL_INTERBLOCvsINTERBLOC      = 303
  integer, parameter :: CL_INTERBLOCvsENTREE_SUPER   = 305
  integer, parameter :: CL_INTERBLOCvsPAROI          = 309
  integer, parameter :: CL_INTERBLOCvsPAROI_ADIAB    = 3010
  integer, parameter :: CL_INTERBLOCvsINJECTION      = 3011
  integer, parameter :: CL_INTERBLOCvsPERIODIQUE     = 3012

  integer, parameter :: CL_INTERPROCvs = 40

  integer, parameter :: CL_FLUX_NULvsFLUX_NUL        = 101
  integer, parameter :: CL_FLUX_NULvsSYMETRIE        = 102
  integer, parameter :: CL_FLUX_NULvsINTERBLOC       = 103
  integer, parameter :: CL_FLUX_NULvsENTREE_SUPER    = 105
  integer, parameter :: CL_FLUX_NULvsENTREE_SUB      = 106
  integer, parameter :: CL_FLUX_NULvsSORTIE_SUB      = 108
  integer, parameter :: CL_FLUX_NULvsPAROI           = 109
  integer, parameter :: CL_FLUX_NULvsPAROI_ADIAB     = 1010

  integer, parameter :: CL_SYMETRIEvsFLUX_NUL        = 201
  integer, parameter :: CL_SYMETRIEvsSYMETRIE        = 202
  integer, parameter :: CL_SYMETRIEvsINTERBLOC       = 203
  integer, parameter :: CL_SYMETRIEvsENTREE_SUPER    = 205
  integer, parameter :: CL_SYMETRIEvsENTREE_SUB      = 206
  integer, parameter :: CL_SYMETRIEvsSORTIE_SUB      = 208
  integer, parameter :: CL_SYMETRIEvsPAROI           = 209
  integer, parameter :: CL_SYMETRIEvsPAROI_ADIAB     = 2010

  integer, parameter :: CL_ENTREE_SUPERvsFLUX_NUL    = 501
  integer, parameter :: CL_ENTREE_SUPERvsSYMETRIE    = 502
  integer, parameter :: CL_ENTREE_SUPERvsINTERBLOC   = 503
  integer, parameter :: CL_ENTREE_SUPERvsPAROI       = 509
  integer, parameter :: CL_ENTREE_SUPERvsPAROI_ADIAB = 5010

  integer, parameter :: CL_ENTREE_SUBvsFLUX_NUL      = 601
  integer, parameter :: CL_ENTREE_SUBvsSYMETRIE      = 602
  integer, parameter :: CL_ENTREE_SUBvsPAROI         = 609
  integer, parameter :: CL_ENTREE_SUBvsPAROI_ADIAB   = 6010

  integer, parameter :: CL_SORTIE_SUBvsFLUX_NUL      = 801
  integer, parameter :: CL_SORTIE_SUBvsSYMETRIE      = 802
  integer, parameter :: CL_SORTIE_SUBvsPAROI         = 809
  integer, parameter :: CL_SORTIE_SUBvsPAROI_ADIAB   = 8010

  integer, parameter :: CL_PAROIvsFLUX_NUL           = 901
  integer, parameter :: CL_PAROIvsSYMETRIE           = 902
  integer, parameter :: CL_PAROIvsINTERBLOC          = 903
  integer, parameter :: CL_PAROIvsENTREE_SUPER       = 905
  integer, parameter :: CL_PAROIvsENTREE_SUB         = 906
  integer, parameter :: CL_PAROIvsSORTIE_SUB         = 908
  integer, parameter :: CL_PAROIvsPAROI              = 909
  integer, parameter :: CL_PAROIvsPAROI_ADIAB        = 9010
  integer, parameter :: CL_PAROIvsINJECTION          = 9011
  integer, parameter :: CL_PAROIvsPERIODIQUE         = 9012

  integer, parameter :: CL_PAROI_ADIABvsFLUX_NUL     = 1001
  integer, parameter :: CL_PAROI_ADIABvsSYMETRIE     = 1002
  integer, parameter :: CL_PAROI_ADIABvsINTERBLOC    = 1003
  integer, parameter :: CL_PAROI_ADIABvsENTREE_SUPER = 1005
  integer, parameter :: CL_PAROI_ADIABvsENTREE_SUB   = 1006
  integer, parameter :: CL_PAROI_ADIABvsSORTIE_SUB   = 1008
  integer, parameter :: CL_PAROI_ADIABvsPAROI        = 1009
  integer, parameter :: CL_PAROI_ADIABvsPAROI_ADIAB  = 10010
  integer, parameter :: CL_PAROI_ADIABvsINJECTION    = 10011
  integer, parameter :: CL_PAROI_ADIABvsPERIODIQUE   = 10012

  integer, parameter :: CL_INJECTIONvsFLUX_NUL     = 1101
  integer, parameter :: CL_INJECTIONvsSYMETRIE     = 1102
  integer, parameter :: CL_INJECTIONvsINTERBLOC    = 1103
  integer, parameter :: CL_INJECTIONvsENTREE_SUPER = 1105
  integer, parameter :: CL_INJECTIONvsENTREE_SUB   = 1106
  integer, parameter :: CL_INJECTIONvsSORTIE_SUB   = 1108
  integer, parameter :: CL_INJECTIONvsPAROI        = 1109
  integer, parameter :: CL_INJECTIONvsPAROI_ADIAB  = 11010

  integer, parameter :: CL_PERIODIQUEvsINTERBLOC   = 1203
  integer, parameter :: CL_PERIODIQUEvsPAROI       = 1209
  integer, parameter :: CL_PERIODIQUEvsPAROI_ADIAB = 12010
  integer, parameter :: CL_PERIODIQUEvsPERIODIQUE  = 12012
 

!!! constantes physiques
  real*8 , parameter :: R_gaz = 8.314d0    ! constante universelle des gaz
  real*8 , parameter :: ONEATM = 101325.d0 ! 1 atmosphère en Pa

!!! epsilon et tolerance sur la fraction volumique
  real*8, parameter :: epsilon_z = 0.d0 !1.d-10
  real*8, parameter :: tol_z = 0.d0!1.d-9
contains

  subroutine def_indices_fluide(nb_especes)
    integer, intent(in) :: nb_especes

    irho1z1 = 1
    irho2z2 = irho1z1 + nb_especes
    irhou   = irho2z2 + 1
    irhov   = irhou + 1
    irhoE   = irhov + 1
    iz      = irhoE + 1

    ivol = 1
    iy   = 2
    iu   = iy + nb_especes
    iv   = iu + 1
    iE   = iv + 1

    irho1 = 1
    irho2 = irho1 + nb_especes
    irho1eps1 = irho2 + 1
    irho2eps2 = irho1eps1 + nb_especes
    iux = irho2eps2 + 1
    iuy = iux + 1
    izz = iuy + 1

  end subroutine def_indices_fluide

end module parametres_fluide
