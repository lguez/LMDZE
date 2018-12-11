module YOMCST

  implicit none

  ! From phylmd/YOMCST.h, version 1.2 2005/06/06 13:16:33

  ! A1.1.bis Constantes concernant l'orbite de la Terre
  ! Default values from AMIP II

  REAL:: R_ecc = 0.016715 ! excentricité
  real:: R_peri = 102.7 ! équinoxe
  real:: R_incl = 23.441 ! inclinaison, en degr\'es

contains

  subroutine read_YOMCST

    use unit_nml_m, only: unit_nml

    namelist /YOMCST_nml/ R_ecc, R_peri, R_incl

    !------------------------------------

    print *, "Enter namelist 'YOMCST_nml'."
    read(unit=*, nml=YOMCST_nml)
    write(unit_nml, nml=YOMCST_nml)

  end subroutine read_YOMCST

end module YOMCST
