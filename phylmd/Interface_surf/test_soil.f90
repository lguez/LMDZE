program test_soil

  use jumble, only: new_unit
  use mpi_f08, only: mpi_init, mpi_finalize, mpi_comm_size, mpi_comm_world, &
       mpi_abort
  use xios, only: xios_initialize, xios_finalize, xios_context_initialize, &
       xios_context_finalize, xios_close_context_definition

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

  integer unit, i, n_proc, return_comm

  !----------------------------------------------------------------

  call mpi_init
  call mpi_comm_size(mpi_comm_world, n_proc)
  if (n_proc /= 1) call mpi_abort(mpi_comm_world, errorcode = 1)
  call xios_initialize("LMDZE", return_comm = return_comm)
  CALL xios_context_initialize("LMDZE_context", return_comm)

  call set_unit_nml
  open(unit_nml, file="used_namelists.txt", status="replace", action="write")
  CALL conf_gcm
  call xios_close_context_definition
  tsoil = 270.
  call soil(is_ter, [0.], [288.], tsoil, soilcap, soilflux)
  call new_unit(unit)
  open(unit, file = "test_soil.csv", status = "replace", action = "write")
  write(unit, fmt = *) "K"
  write(unit, fmt = *) "T_soil"
  
  do i = 1, nsoilmx
     write(unit, fmt = *) tsoil(1, i)
  end do

  close(unit)
  print *, "soilcap = ", soilcap
  print *, "soilflux = ", soilflux
  CALL xios_context_finalize
  call xios_finalize
  call mpi_finalize

end program test_soil
