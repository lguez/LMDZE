PROGRAM ce0l

  ! This program sets the initial and boundary values.

  use comdissnew, only: read_comdissnew
  use conf_gcm_m, only: conf_gcm
  use dimensions, only: iim, jjm
  use dimphy, only: klon, init_dimphy
  use dynetat0_chosen_m, only: read_serre
  use etat0_m, only: etat0
  use grilles_gcm_netcdf_sub_m, only: grilles_gcm_netcdf_sub
  use indicesol, only: nbsrf
  use limit_m, only: limit
  use paramet_m, only: paramet
  use unit_nml_m, only: unit_nml, set_unit_nml

  implicit none

  REAL, ALLOCATABLE:: phis(:, :) ! (iim + 1, jjm + 1)
  ! surface geopotential, in m2 s-2

  REAL, ALLOCATABLE:: pctsrf(:, :) ! (klon, nbsrf)
  ! ("pctsrf(i, :)" is the composition of the surface at horizontal
  ! position "i")

  !-------------------------------------

  ALLOCATE(phis(iim + 1, jjm + 1))
  ALLOCATE(pctsrf(klon, nbsrf))
  
  call set_unit_nml
  open(unit_nml, file="used_namelists.txt", status="replace", action="write")
  call paramet
  call init_dimphy
  CALL conf_gcm
  call read_comdissnew
  call read_serre
  CALL etat0(phis, pctsrf)
  CALL limit(pctsrf)
  CALL grilles_gcm_netcdf_sub(phis)
  close(unit_nml)
  print *, "ce0l: done"

END PROGRAM ce0l
