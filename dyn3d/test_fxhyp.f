program test_fxhyp

  USE dimensions, ONLY: iim
  use fxhyp_m, only: fxhyp
  use jumble, only: new_unit
  use read_serre_m, only: read_serre
  use unit_nml_m, only: unit_nml

  implicit none

  REAL, dimension(iim + 1):: xprimm025, rlonv, xprimv, rlonu, xprimu, xprimp025
  integer unit

  !--------------------------------------------------------

  call new_unit(unit_nml)
  open(unit_nml, file="used_namelists.txt", status="replace", action="write")
  call read_serre
  call fxhyp(xprimm025, rlonv, xprimv, rlonu, xprimu, xprimp025)
  close(unit_nml)

  unit = unit_nml
  open(unit, file="test_fxhyp_out.txt", status="replace", action="write")
  write(unit, fmt = *) "xprimm025 = ", xprimm025
  write(unit, fmt = *) "rlonv = ", rlonv
  write(unit, fmt = *) "xprimv = ", xprimv
  write(unit, fmt = *) "rlonu = ", rlonu
  write(unit, fmt = *) "xprimu = ", xprimu
  write(unit, fmt = *) "xprimp025 = ", xprimp025
  close(unit)

end program test_fxhyp
