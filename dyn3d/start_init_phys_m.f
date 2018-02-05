MODULE start_init_phys_m

  ! From startvar.F, version 1.4
  ! 2006/01/27 15:14:22 Fairhead

  IMPLICIT NONE

CONTAINS

  SUBROUTINE start_init_phys(tsol_2d, qsol_2d)

    use conf_dat2d_m, only: conf_dat2d
    use dimens_m, only: iim, jjm
    use dynetat0_m, only: rlonu, rlatv
    use gr_int_dyn_m, only: gr_int_dyn
    use inter_barxy_m, only: inter_barxy
    use netcdf, only: nf90_nowrite
    use netcdf95, only: nf95_open, nf95_close, nf95_get_var, nf95_inq_varid, &
         nf95_gw_var, find_coord
    use nr_util, only: assert, pi

    REAL, intent(out):: tsol_2d(:, :), qsol_2d(:, :) ! (iim + 1, jjm + 1)

    ! Variables local to the procedure:

    INTEGER iml_phys, jml_phys, ncid, varid
    REAL, ALLOCATABLE:: lon_rad(:), lat_rad(:)
    REAL, ALLOCATABLE:: lon_ini(:), lat_ini(:) ! longitude and latitude in rad
    REAL, ALLOCATABLE:: var_ana(:, :)
    real tmp_var(iim, jjm + 1)

    !-----------------------------------

    print *, "Call sequence information: start_init_phys"

    call assert((/size(tsol_2d, 1), size(qsol_2d, 1)/) == iim + 1, &
         "start_init_phys 1")
    call assert((/size(tsol_2d, 2), size(qsol_2d, 2)/) == jjm + 1, &
         "start_init_phys 2")

    call nf95_open('ECPHY.nc', nf90_nowrite, ncid)

    call find_coord(ncid, varid=varid, std_name="longitude")
    call nf95_gw_var(ncid, varid, lon_ini)
    lon_ini = lon_ini * pi / 180. ! convert to rad
    iml_phys = size(lon_ini)

    call find_coord(ncid, varid=varid, std_name="latitude")
    call nf95_gw_var(ncid, varid, lat_ini)
    lat_ini = lat_ini * pi / 180. ! convert to rad
    jml_phys = size(lat_ini)

    ! Allocate the space we will need to get the data out of this file
    ALLOCATE(var_ana(iml_phys, jml_phys))

    ALLOCATE(lon_rad(iml_phys))
    ALLOCATE(lat_rad(jml_phys))

    ! Surface temperature:
    call nf95_inq_varid(ncid, 'ST', varid)
    call nf95_get_var(ncid, varid, var_ana)
    CALL conf_dat2d(lon_ini, lat_ini, lon_rad, lat_rad, var_ana)
    CALL inter_barxy(lon_rad, lat_rad(:jml_phys -1), var_ana, rlonu(:iim), &
         rlatv, tmp_var) 
    tsol_2d = gr_int_dyn(tmp_var)

    ! Soil moisture:
    call nf95_inq_varid(ncid, 'CDSW', varid)
    call nf95_get_var(ncid, varid, var_ana)
    CALL conf_dat2d(lon_ini, lat_ini, lon_rad, lat_rad, var_ana)
    CALL inter_barxy(lon_rad, lat_rad(:jml_phys -1), var_ana, rlonu(:iim), &
            rlatv, tmp_var)
    qsol_2d = gr_int_dyn(tmp_var)

    call nf95_close(ncid)

  END SUBROUTINE start_init_phys

END MODULE start_init_phys_m
