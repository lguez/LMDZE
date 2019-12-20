program test_ozonecm

  ! This is a program in Fortran 2003.

  ! This program computes values of ozone abundance from Royer, on a
  ! latitude-pressure grid, and writes the values to a NetCDF file.
  ! The pressure grid is "presnivs" from "disvert".

  use nr_util, only: arth, assert
  use netcdf95, only: nf95_create, nf95_def_dim, nf95_def_var, nf95_put_att, &
       nf95_enddef, nf95_put_var, nf95_close
  use netcdf, only: nf90_clobber, nf90_float, nf90_global

  use dimensions, only: jjm, llm
  USE dimphy, ONLY : klon
  USE dimsoil, ONLY : nsoilmx
  use disvert_m, only: disvert, ap, bp, preff, presnivs
  use dynetat0_chosen_m, only: read_serre
  USE indicesol, ONLY : nbsrf
  use ozonecm_m, only: ozonecm
  use phyetat0_m, only: rlat, phyetat0
  use unit_nml_m, only: unit_nml, set_unit_nml

  implicit none

  real p(llm + 1) ! pressure at LMDZ layer interface, in Pa
  real wo(klon, llm, 360) ! column density of ozone in a cell, in kDU
  real tro3(klon, llm, 360) ! mole fraction of ozone
  integer julien, k
  real, parameter:: RG = 9.80665 ! acceleration of gravity, in m s-2
  real, parameter:: dobson_u = 2.1415e-05 ! Dobson unit, in kg m-2

  REAL pctsrf(klon, nbsrf)
  REAL ftsol(klon, nbsrf)
  REAL ftsoil(klon, nsoilmx, nbsrf)
  REAL qsurf(klon, nbsrf)
  REAL qsol(klon) ! column-density of water in soil, in kg m-2
  REAL snow(klon, nbsrf)
  REAL albe(klon, nbsrf)
  REAL rain_fall(klon)
  REAL snow_fall(klon)
  real solsw(klon)
  REAL sollw(klon)
  real fder(klon)
  REAL radsol(klon)
  REAL frugs(klon, nbsrf)
  REAL agesno(klon, nbsrf)
  REAL zmea(klon)
  REAL zstd(klon)
  REAL zsig(klon)
  REAL zgam(klon)
  REAL zthe(klon)
  REAL zpic(klon)
  REAL zval(klon)
  REAL t_ancien(klon, llm), q_ancien(klon, llm)
  LOGICAL ancien_ok
  real rnebcon(klon, llm), ratqs(klon, llm)
  REAL clwcon(klon, llm), run_off_lic_0(klon)
  real sig1(klon, llm) ! section adiabatic updraft

  real w01(klon, llm) 
  ! vertical velocity within adiabatic updraft

  integer ncid_startphy

  ! For NetCDF:
  integer ncid, dimid_time, dimid_plev, dimid_latitude
  integer varid_time, varid_plev, varid_latitude, varid_tro3

  !---------------------

  call assert(klon == jjm + 1, "test_ozonecm: iim should be 1")

  call set_unit_nml
  open(unit_nml, file="used_namelists.txt", status="replace", action="write")

  call read_serre
  call disvert
  p = ap + bp * preff
  call phyetat0(pctsrf, ftsol, ftsoil, qsurf, qsol, snow, albe, rain_fall, &
       snow_fall, solsw, sollw, fder, radsol, frugs, agesno, zmea, zstd, zsig, &
       zgam, zthe, zpic, zval, t_ancien, q_ancien, ancien_ok, rnebcon, ratqs, &
       clwcon, run_off_lic_0, sig1, w01, ncid_startphy)

  do julien = 1, 360
     wo(:, :, julien) = ozonecm(REAL(julien), spread(p, dim=1, ncopies=jjm+1))
  end do

  close(unit_nml)

  forall (k=1: llm)
     tro3(:, k, :) &
          = wo(:, k, :) * 1e3 * dobson_u / (p(k)-p(k+1)) * rg / 48. * 29.
  end forall

  call nf95_create("O3_Royer.nc", nf90_clobber, ncid)
  call nf95_def_dim(ncid, "time", 360, dimid_time)
  call nf95_def_dim(ncid, "plev", llm, dimid_plev)
  call nf95_def_dim(ncid, "latitude", jjm + 1, dimid_latitude)
  call nf95_def_var(ncid, "time", nf90_float, dimid_time, varid_time)
  call nf95_def_var(ncid, "plev", nf90_float, dimid_plev, varid_plev)
  call nf95_def_var(ncid, "latitude", nf90_float, dimid_latitude, &
       varid_latitude)
  call nf95_def_var(ncid, "tro3", nf90_float, &
       (/dimid_latitude, dimid_plev, dimid_time/), varid_tro3)

  call nf95_put_att(ncid, varid_time, "units", "days since 1976-1-1")
  call nf95_put_att(ncid, varid_time, "calendar", "360_day")
  call nf95_put_att(ncid, varid_time, "standard_name", "time")

  call nf95_put_att(ncid, varid_plev, "units", "mbar")
  call nf95_put_att(ncid, varid_plev, "long_name", "pressure")
  call nf95_put_att(ncid, varid_plev, "standard_name", "air_pressure")

  call nf95_put_att(ncid, varid_latitude, "units", "degrees_north")
  call nf95_put_att(ncid, varid_latitude, "standard_name", "latitude")

  call nf95_put_att(ncid, varid_tro3, "long_name", "O3 mole fraction")
  call nf95_put_att(ncid, varid_tro3, "standard_name", &
       "mole_fraction_of_ozone_in_air")

  call nf95_put_att(ncid, nf90_global, "title", &
       "ozone climatology from J.-F. Royer")

  call nf95_enddef(ncid)

  call nf95_put_var(ncid, varid_time, (/(julien - 0.5, julien = 1, 360)/))

  call nf95_put_var(ncid, varid_plev, presnivs(llm:1:-1) / 100.)
  ! sort in ascending order and convert to hPa

  call nf95_put_var(ncid, varid_latitude, rlat)

  call nf95_put_var(ncid, varid_tro3, tro3(:, llm:1:-1, :))
  ! "plev" is in ascending order whereas "p" is descending

  call nf95_close(ncid)

  print *, 'The file "O3_Royer.nc" has been created.'

end program test_ozonecm