MODULE start_init_phys_m

  ! From startvar.F, version 1.4
  ! 2006/01/27 15:14:22 Fairhead

  IMPLICIT NONE

CONTAINS

  SUBROUTINE start_init_phys(tsol_2d, qsol_2d)

    USE flincom, only: flininfo, flinopen_nozoom, flinclo
    use conf_dat2d_m, only: conf_dat2d
    use inter_barxy_m, only: inter_barxy
    use gr_int_dyn_m, only: gr_int_dyn
    use comgeom, only: rlonu, rlatv
    use dimens_m, only: iim, jjm
    use nr_util, only: assert
    use netcdf, only: nf90_nowrite
    use netcdf95, only: nf95_open, nf95_close, nf95_get_var, nf95_inq_varid

    REAL, intent(out):: tsol_2d(:, :), qsol_2d(:, :) ! (iim + 1, jjm + 1)

    ! Variables local to the procedure:

    INTEGER fid_phys, iml_phys, jml_phys, ncid, varid
    REAL, ALLOCATABLE, DIMENSION(:, :):: lon_phys, lat_phys
    REAL date, dt
    REAL, ALLOCATABLE:: levphys_ini(:)

    INTEGER itau(1)
    INTEGER  llm_tmp, ttm_tmp

    REAL, ALLOCATABLE:: lon_rad(:), lat_rad(:)
    REAL, ALLOCATABLE:: lon_ini(:), lat_ini(:)
    REAL, ALLOCATABLE:: var_ana(:, :)
    real tmp_var(iim, jjm + 1)

    !-----------------------------------

    print *, "Call sequence information: start_init_phys"

    call assert((/size(tsol_2d, 1), size(qsol_2d, 1)/) == iim + 1, &
         "start_init_phys 1")
    call assert((/size(tsol_2d, 2), size(qsol_2d, 2)/) == jjm + 1, &
         "start_init_phys 2")

    CALL flininfo('ECPHY.nc', iml_phys, jml_phys, llm_tmp, ttm_tmp, fid_phys)

    ALLOCATE(lat_phys(iml_phys, jml_phys))
    ALLOCATE(lon_phys(iml_phys, jml_phys))
    ALLOCATE(levphys_ini(llm_tmp))

    CALL flinopen_nozoom(iml_phys, jml_phys, llm_tmp, lon_phys, lat_phys, &
         levphys_ini, ttm_tmp, itau, date, dt, fid_phys)
    CALL flinclo(fid_phys)

    DEALLOCATE(levphys_ini)

    ! Allocate the space we will need to get the data out of this file
    ALLOCATE(var_ana(iml_phys, jml_phys))

    !   In case we have a file which is in degrees we do the transformation
    ALLOCATE(lon_rad(iml_phys))
    ALLOCATE(lon_ini(iml_phys))

    IF ( MAXVAL(lon_phys) > 2.0 * ASIN(1.0) ) THEN
       lon_ini = lon_phys(:, 1) * 2.0 * ASIN(1.0) / 180.0
    ELSE
       lon_ini = lon_phys(:, 1) 
    ENDIF

    ALLOCATE(lat_rad(jml_phys))
    ALLOCATE(lat_ini(jml_phys))

    IF ( MAXVAL(lat_phys) > 2.0 * ASIN(1.0) ) THEN
       lat_ini = lat_phys(1, :) * 2.0 * ASIN(1.0) / 180.0
    ELSE
       lat_ini = lat_phys(1, :) 
    ENDIF

    call nf95_open('ECPHY.nc', nf90_nowrite, ncid)

    ! We get the two standard variables
    ! 'ST': surface temperature
    call nf95_inq_varid(ncid, 'ST', varid)
    call nf95_get_var(ncid, varid, var_ana)
    CALL conf_dat2d(lon_ini, lat_ini, lon_rad, lat_rad, var_ana)
    CALL inter_barxy(lon_rad, lat_rad(:jml_phys -1), var_ana, rlonu(:iim), &
         rlatv, tmp_var) 
    tsol_2d = gr_int_dyn(tmp_var)

    ! Soil moisture
    call nf95_inq_varid(ncid, 'CDSW', varid)
    call nf95_get_var(ncid, varid, var_ana)
    CALL conf_dat2d(lon_ini, lat_ini, lon_rad, lat_rad, var_ana)
    CALL inter_barxy(lon_rad, lat_rad(:jml_phys -1), var_ana, rlonu(:iim), &
            rlatv, tmp_var)
    qsol_2d = gr_int_dyn(tmp_var)

    call nf95_close(ncid)

  END SUBROUTINE start_init_phys

END MODULE start_init_phys_m
