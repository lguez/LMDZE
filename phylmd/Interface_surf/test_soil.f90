program test_soil

  use comconst, only: iniconst
  use conf_gcm_m, only: conf_gcm
  USE dimsoil, only: nsoilmx
  USE indicesol, only: is_ter
  use soil_m, only: soil
  use unit_nml_m, only: unit_nml, set_unit_nml

  implicit none

  real tsoil(1, nsoilmx) ! (knon, nsoilmx)
  ! temperature inside the ground (K), layer 1 nearest to the surface

  REAL soilcap(1) ! (knon)
  ! specific heat per unit surface (in J m-2 K-1)

  REAL soilflux(1) ! (knon) 
  ! surface diffusive flux from ground (W m-2)

  !----------------------------------------------------------------

  call set_unit_nml
  open(unit_nml, file="used_namelists.txt", status="replace", action="write")
  CALL conf_gcm
  CALL iniconst
  tsoil = 270.
  call soil(is_ter, [0.], [288.], tsoil, soilcap, soilflux)
  print *, "tsoil = ", tsoil
  print *, "soilcap = ", soilcap
  print *, "soilflux = ", soilflux

end program test_soil
