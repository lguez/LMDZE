program test_inter_barxy

  use comconst, only: dtvr, daysec, iniconst
  use comgeom, only: inigeom
  use conf_gcm_m, only: conf_gcm, day_step
  use dimens_m, only: iim, jjm
  USE dynetat0_m, only: rlonu, rlatv
  use disvert_m, only: pa
  use inter_barxy_m, only: inter_barxy
  USE nr_util, ONLY : pi
  use read_serre_m, only: read_serre

  implicit none

  integer:: iml_dyn = 3, jml_dyn = 6
  real, allocatable:: lon_rad(:), lon_ini(:)
  real, allocatable:: lat_rad(:), lat_ini(:)
  integer i
  real, allocatable:: var_ana3d(:, :)
  REAL var_tmp2d(iim, jjm + 1)

  namelist /main/iml_dyn, jml_dyn

  !------------------------

  print *, "Enter namelist 'main'."
  read (unit=*, nml=main)
  write(unit=*, nml=main)

  allocate(lon_rad(iml_dyn), lon_ini(iml_dyn))
  allocate(lat_rad(jml_dyn-1), lat_ini(jml_dyn))
  allocate(var_ana3d(iml_dyn, jml_dyn))

  CALL conf_gcm
  dtvr = daysec / real(day_step)
  print *, 'dtvr = ', dtvr
  pa = 5e4
  CALL iniconst
  call read_serre
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

end program test_inter_barxy
