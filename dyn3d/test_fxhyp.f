program test_fxhyp

  USE dimensions, ONLY: iim
  use dynetat0_m, only: read_serre
  use fxhyp_m, only: fxhyp
  use unit_nml_m, only: unit_nml, set_unit_nml

  implicit none

  REAL, dimension(iim + 1):: xprimm025, rlonv, xprimv, rlonu, xprimu, xprimp025

  !--------------------------------------------------------

  call set_unit_nml
  open(unit_nml, file="used_namelists.txt", status="replace", action="write")
  call read_serre
  call fxhyp(xprimm025, rlonv, xprimv, rlonu, xprimu, xprimp025)
  close(unit_nml)

  ! We can use the same unit number although we are not writing a namelist:
  open(unit_nml, file="test_fxhyp_out.txt", status="replace", action="write")
  write(unit_nml, fmt = *) "xprimm025 = ", xprimm025
  write(unit_nml, fmt = *) "rlonv = ", rlonv
  write(unit_nml, fmt = *) "xprimv = ", xprimv
  write(unit_nml, fmt = *) "rlonu = ", rlonu
  write(unit_nml, fmt = *) "xprimu = ", xprimu
  write(unit_nml, fmt = *) "xprimp025 = ", xprimp025
  close(unit_nml)

end program test_fxhyp
