PROGRAM ce0l

  ! This program sets the initial and boundary values.

  use comdissnew, only: read_comdissnew
  use conf_gcm_m, only: conf_gcm
  use dimens_m, only: iim, jjm
  use etat0_mod, only: etat0
  use grilles_gcm_netcdf_sub_m, only: grilles_gcm_netcdf_sub
  use jumble, only: new_unit
  use limit_mod, only: limit
  use read_serre_m, only: read_serre
  use unit_nml_m, only: unit_nml

  implicit none

  REAL phis(iim + 1, jjm + 1) ! surface geopotential, in m2 s-2

  !-------------------------------------

  call new_unit(unit_nml)
  open(unit_nml, file="used_namelists.txt", status="replace", action="write")
  CALL conf_gcm
  call read_comdissnew
  call read_serre
  CALL etat0(phis)
  CALL limit
  CALL grilles_gcm_netcdf_sub(phis)
  close(unit_nml)

END PROGRAM ce0l
