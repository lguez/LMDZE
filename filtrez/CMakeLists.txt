target_sources(LMDZE_objects PRIVATE filtreg_hemisph.f90
  filtreg_scal.f90 filtreg_v.f90 inifgn.f90 inifilr.f90
  inifilr_hemisph.f90)

add_executable(test_inifilr test_inifilr.f90)

target_link_libraries(test_inifilr PRIVATE LMDZE_objects
  NetCDF95::netcdf95 NetCDF_Fortran::netcdff NR_util::nr_util)
