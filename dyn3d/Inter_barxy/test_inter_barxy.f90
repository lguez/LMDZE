program test_inter_barxy

  use comdissnew, only: read_comdissnew
  use comgeom, only: inigeom
  use dimensions, only: iim, jjm, set_dimensions
  USE dynetat0_m, only: rlonu, rlatv, fyhyp, fxhyp
  use dynetat0_chosen_m, only: read_serre
  use inter_barxy_m, only: inter_barxy
  USE nr_util, ONLY : pi
  use paramet_m, only: paramet
  use unit_nml_m, only: set_unit_nml, unit_nml

  implicit none

  integer:: iml_dyn = 3, jml_dyn = 6
  real, allocatable:: lon_rad(:), lon_ini(:)
  real, allocatable:: lat_rad(:), lat_ini(:)
  integer i
  real, allocatable:: var_ana3d(:, :)
  REAL, allocatable:: var_tmp2d(:, :) ! (iim, jjm + 1)

  namelist /main/iml_dyn, jml_dyn

  !------------------------

  call set_unit_nml
  open(unit_nml, file="used_namelists.txt", status="replace", action="write")

  call set_dimensions
  allocate(var_tmp2d(iim, jjm + 1))
  call paramet
  
  print *, "Enter namelist 'main'."
  read (unit=*, nml=main)
  write(unit=unit_nml, nml=main)

  allocate(lon_rad(iml_dyn), lon_ini(iml_dyn))
  allocate(lat_rad(jml_dyn-1), lat_ini(jml_dyn))
  allocate(var_ana3d(iml_dyn, jml_dyn))

  call read_comdissnew
  call read_serre
  CALL fyhyp
  CALL fxhyp
  CALL inigeom

  lon_ini = - pi + 2 * pi / iml_dyn * (/(i, i = 0, iml_dyn - 1)/)

  forall (i = 1: iml_dyn - 1) lon_rad(i) = (lon_ini(i) + lon_ini(i + 1)) / 2
  lon_rad(iml_dyn) = (lon_ini(iml_dyn) + pi) / 2

  lat_ini = pi / 2 - pi / (jml_dyn - 1) * (/(i, i = 0, jml_dyn - 1)/)
  forall (i = 1: jml_dyn - 1) lat_rad(i) = (lat_ini(i) + lat_ini(i + 1)) / 2

  !!forall (i = 1: jml_dyn) var_ana3d(:, i) = i
  forall (i = 1: jml_dyn) var_ana3d(:, i) = 240. + real(i) / 2

  call inter_barxy(lon_rad, lat_rad, var_ana3d, rlonu(:iim), rlatv, var_tmp2d)

  !!print *, "lon_rad * 180. / pi = ", lon_rad * 180. / pi
  print *, "lat_rad * 180. / pi = ", lat_rad * 180. / pi
  !!print *, "var_tmp2d = ", var_tmp2d

  print *, "minval(var_tmp2d, dim=1) = ", minval(var_tmp2d, dim=1)
  print *, "maxval(var_tmp2d, dim=1) = ", maxval(var_tmp2d, dim=1)
  close(unit_nml)

end program test_inter_barxy
