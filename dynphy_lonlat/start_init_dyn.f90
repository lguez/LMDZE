MODULE start_init_dyn_m

  ! From startvar.F, version 1.4
  ! January 27th, 2006 15:14:22

  INTEGER, protected:: iml_dyn, jml_dyn, llm_dyn

  REAL, protected, allocatable:: lon_ini(:), lat_ini(:)
  ! longitude and latitude from the input file, converted to rad

  real, protected, allocatable:: levdyn_ini(:)

CONTAINS

  SUBROUTINE start_init_dyn(tsol_2d, phis, ps)

    use netcdf, only: nf90_nowrite
    use netcdf95, only: nf95_open, nf95_close, nf95_get_var, nf95_inq_varid, &
         nf95_gw_var, find_coord
    use jumble, only: assert, pi

    use comgeom, only: aire_2d, apoln, apols
    use conf_dat2d_m, only: conf_dat2d
    use dimensions, only: iim, jjm
    use dynetat0_m, only: rlonu, rlatv
    use gr_int_dyn_m, only: gr_int_dyn
    use inter_barxy_m, only: inter_barxy

    REAL, intent(in):: tsol_2d(:, :) ! (iim + 1, jjm + 1)

    REAL, intent(in):: phis(:, :) ! (iim + 1, jjm + 1)
    ! surface geopotential, in m2 s-2

    REAL, intent(out):: ps(:, :) ! (iim + 1, jjm + 1) surface pressure, in Pa

    ! Local:

    INTEGER ncid, varid
    REAL, ALLOCATABLE:: lon_rad(:), lat_rad(:)
    REAL, ALLOCATABLE:: var_ana(:, :)
    real z(iim + 1, jjm + 1)
    real tmp_var(iim, jjm + 1)

    !--------------------------

    print *, "Call sequence information: start_init_dyn"
    call assert((/size(tsol_2d, 1), size(phis, 1), size(ps, 1)/) == iim + 1, &
         "start_init_dyn size 1")
    call assert((/size(tsol_2d, 2), size(phis, 2), size(ps, 2)/) == jjm + 1, &
         "start_init_dyn size 2")

    call nf95_open('ECDYN.nc', nf90_nowrite, ncid)

    call find_coord(ncid, varid=varid, std_name="longitude")
    call nf95_gw_var(ncid, varid, lon_ini)
    lon_ini = lon_ini * pi / 180. ! convert to rad
    iml_dyn = size(lon_ini)
    print *, "iml_dyn = ", iml_dyn

    call find_coord(ncid, varid=varid, std_name="latitude")
    call nf95_gw_var(ncid, varid, lat_ini)
    lat_ini = lat_ini * pi / 180. ! convert to rad
    jml_dyn = size(lat_ini)
    print *, "jml_dyn = ", jml_dyn

    call find_coord(ncid, varid=varid, std_name="plev")
    call nf95_gw_var(ncid, varid, levdyn_ini)
    llm_dyn = size(levdyn_ini)
    print *, "llm_dyn = ", llm_dyn

    ALLOCATE(var_ana(iml_dyn, jml_dyn), lon_rad(iml_dyn), lat_rad(jml_dyn))

    ! 'Z': Surface geopotential
    call nf95_inq_varid(ncid, 'Z', varid)
    call nf95_get_var(ncid, varid, var_ana)
    CALL conf_dat2d(lon_ini, lat_ini, lon_rad, lat_rad, var_ana)
    CALL inter_barxy(lon_rad, lat_rad(:jml_dyn -1), var_ana, rlonu(:iim), &
         rlatv, tmp_var)
    z = gr_int_dyn(tmp_var)

    ! 'SP': Surface pressure
    call nf95_inq_varid(ncid, 'SP', varid)
    call nf95_get_var(ncid, varid, var_ana)
    CALL conf_dat2d(lon_ini, lat_ini, lon_rad, lat_rad, var_ana)
    CALL inter_barxy(lon_rad, lat_rad(:jml_dyn -1), var_ana, rlonu(:iim), &
         rlatv, tmp_var) 
    ps = gr_int_dyn(tmp_var)

    call nf95_close(ncid)

    ! Adapt the surface pressure to "phis", with the hypsometric equation:
    ps(:iim, :) = ps(:iim, :) &
         * (1. + (z(:iim, :) - phis(:iim, :)) / 287. / tsol_2d(:iim, :))

    ps(iim + 1, :) = ps(1, :)
    ps(:, 1) = SUM(aire_2d(:iim, 1) * ps(:iim, 1)) / apoln
    ps(:, jjm + 1) = SUM(aire_2d(:iim, jjm + 1) * ps(:iim, jjm + 1)) / apols

  END SUBROUTINE start_init_dyn

END MODULE start_init_dyn_m
