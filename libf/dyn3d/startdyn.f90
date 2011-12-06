MODULE startdyn

  ! From startvar.F, version 1.4
  ! 2006/01/27 15:14:22 Fairhead

  INTEGER iml_dyn, jml_dyn, llm_dyn

  REAL, ALLOCATABLE:: lon_ini(:), lat_ini(:)
  ! (longitude and latitude from the input file, converted to rad)

  real, ALLOCATABLE:: levdyn_ini(:)

CONTAINS

  SUBROUTINE start_init_dyn(tsol_2d, psol)

    USE flincom, only: flininfo, flinopen_nozoom
    use comgeom, only: aire_2d, apoln, apols
    use conf_dat2d_m, only: conf_dat2d
    use inter_barxy_m, only: inter_barxy
    use comgeom, only: rlonu, rlatv
    use dimens_m, only: iim, jjm
    use gr_int_dyn_m, only: gr_int_dyn
    use start_init_orog_m, only: phis
    use nr_util, only: assert, pi
    use netcdf, only: nf90_nowrite
    use netcdf95, only: nf95_open, nf95_close, nf95_get_var, nf95_inq_varid

    REAL, intent(in):: tsol_2d(:, :) ! (iim + 1, jjm + 1)
    REAL, intent(out):: psol(:, :) ! (iim + 1, jjm + 1) surface pressure, in Pa

    ! Local:

    REAL date, dt
    INTEGER itau(1), ncid, varid, fid_dyn, ttm_dyn
    REAL, ALLOCATABLE:: lon_rad(:), lat_rad(:)

    REAL, ALLOCATABLE:: lon_dyn(:, :), lat_dyn(:, :)
    ! (longitude and latitude from the input file, in rad or degrees)

    REAL, ALLOCATABLE:: var_ana(:, :)
    real z(iim + 1, jjm + 1)
    real tmp_var(iim, jjm + 1)

    !--------------------------

    print *, "Call sequence information: start_init_dyn"
    call assert((/size(tsol_2d, 1), size(psol, 1)/) == iim + 1, &
         "start_init_dyn size 1")
    call assert((/size(tsol_2d, 2), size(psol, 2)/) == jjm + 1, &
         "start_init_dyn size 2")

    CALL flininfo('ECDYN.nc', iml_dyn, jml_dyn, llm_dyn, ttm_dyn, fid_dyn)
    print *, "iml_dyn = ", iml_dyn, ", jml_dyn = ", jml_dyn, &
         ", llm_dyn = ", llm_dyn, ", ttm_dyn = ", ttm_dyn

    ALLOCATE(lat_dyn(iml_dyn, jml_dyn))
    allocate(lon_dyn(iml_dyn, jml_dyn), levdyn_ini(llm_dyn))

    CALL flinopen_nozoom(iml_dyn, jml_dyn, llm_dyn, lon_dyn, lat_dyn, &
         levdyn_ini, ttm_dyn, itau, date, dt, fid_dyn)

    ALLOCATE(var_ana(iml_dyn, jml_dyn), lon_rad(iml_dyn), lon_ini(iml_dyn))

    IF (MAXVAL(lon_dyn) > pi) THEN
       ! Assume "lon_dyn" is in degrees 
       lon_ini = lon_dyn(:, 1) * pi / 180.
    ELSE
       lon_ini = lon_dyn(:, 1) 
    ENDIF

    ALLOCATE(lat_rad(jml_dyn), lat_ini(jml_dyn))

    IF (MAXVAL(lat_dyn) > pi) THEN
       lat_ini = lat_dyn(1, :) * pi / 180.
    ELSE
       lat_ini = lat_dyn(1, :) 
    ENDIF

    call nf95_open('ECDYN.nc', nf90_nowrite, ncid)

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
    psol = gr_int_dyn(tmp_var)

    call nf95_close(ncid)

    psol(:iim, :) = psol(:iim, :) &
         * (1. + (z(:iim, :) - phis(:iim, :)) / 287. / tsol_2d(:iim, :))
    psol(iim + 1, :) = psol(1, :)

    psol(:, 1) = SUM(aire_2d(:iim, 1) * psol(:iim, 1)) / apoln
    psol(:, jjm + 1) = SUM(aire_2d(:iim, jjm + 1) * psol(:iim, jjm + 1)) &
         / apols

  END SUBROUTINE start_init_dyn

END MODULE startdyn
