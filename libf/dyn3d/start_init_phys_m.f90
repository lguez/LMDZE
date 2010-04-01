MODULE start_init_phys_m

  ! From startvar.F, version 1.4
  ! 2006/01/27 15:14:22 Fairhead

  IMPLICIT NONE

  REAL, ALLOCATABLE, SAVE, DIMENSION(:, :):: qsol_2d

CONTAINS

  SUBROUTINE start_init_phys(tsol_2d)

    USE flincom, only: flininfo, flinopen_nozoom, flinget, flinclo
    use conf_dat2d_m, only: conf_dat2d
    use inter_barxy_m, only: inter_barxy
    use gr_int_dyn_m, only: gr_int_dyn
    use comgeom, only: rlonu, rlatv
    use dimens_m, only: iim, jjm

    REAL, intent(out):: tsol_2d(:, :)

    !  LOCAL

    INTEGER fid_phys, iml_phys, jml_phys
    REAL, ALLOCATABLE, DIMENSION(:, :):: lon_phys, lat_phys
    REAL:: date, dt
    REAL, DIMENSION(:), ALLOCATABLE:: levphys_ini
    !ac
    INTEGER:: itau(1)
    INTEGER::  llm_tmp, ttm_tmp

    CHARACTER(len=120) physfname
    LOGICAL:: check=.TRUE.

    REAL, ALLOCATABLE:: lon_rad(:), lat_rad(:)
    REAL, ALLOCATABLE:: lon_ini(:), lat_ini(:)
    REAL, ALLOCATABLE:: var_ana(:, :)
    real tmp_var(iim, jjm + 1)

    !-----------------------------------

    print *, "Call sequence information: start_init_phys"
    if (any(shape(tsol_2d) /= (/iim + 1, jjm + 1/))) stop "start_init_phys"
    physfname = 'ECPHY.nc'
    IF ( check ) print *, 'Opening the surface analysis'
    CALL flininfo(physfname, iml_phys, jml_phys, llm_tmp, ttm_tmp, fid_phys)

    ALLOCATE(lat_phys(iml_phys, jml_phys))
    ALLOCATE(lon_phys(iml_phys, jml_phys))
    ALLOCATE(levphys_ini(llm_tmp))

    CALL flinopen_nozoom(physfname, iml_phys, jml_phys,  &
         llm_tmp, lon_phys, lat_phys, levphys_ini, ttm_tmp,  &
         itau, date, dt, fid_phys)

    DEALLOCATE(levphys_ini)

    ! Allocate the space we will need to get the data out of this file
    ALLOCATE(var_ana(iml_phys, jml_phys))

    !   In case we have a file which is in degrees we do the transformation
    ALLOCATE(lon_rad(iml_phys))
    ALLOCATE(lon_ini(iml_phys))

    IF ( MAXVAL(lon_phys(:, :)) > 2.0 * ASIN(1.0) ) THEN
       lon_ini(:) = lon_phys(:, 1) * 2.0 * ASIN(1.0) / 180.0
    ELSE
       lon_ini(:) = lon_phys(:, 1) 
    ENDIF

    ALLOCATE(lat_rad(jml_phys))
    ALLOCATE(lat_ini(jml_phys))

    IF ( MAXVAL(lat_phys(:, :)) > 2.0 * ASIN(1.0) ) THEN
       lat_ini(:) = lat_phys(1, :) * 2.0 * ASIN(1.0) / 180.0
    ELSE
       lat_ini(:) = lat_phys(1, :) 
    ENDIF

    !   We get the two standard varibales
    !   Surface temperature
    ! 'ST'            : Surface temperature
    CALL flinget(fid_phys, 'ST', iml_phys, jml_phys,  &
         llm_tmp, ttm_tmp, 1, 1, var_ana)
    CALL conf_dat2d(lon_ini, lat_ini, lon_rad, lat_rad, var_ana)
    CALL inter_barxy(lon_rad, lat_rad(:jml_phys -1), var_ana, rlonu(:iim), &
         rlatv, tmp_var) 

    tsol_2d(:, :) = gr_int_dyn(tmp_var)

    ALLOCATE(qsol_2d(iim + 1, jjm + 1))
    ! Soil moisture
    CALL flinget(fid_phys, 'CDSW', iml_phys, jml_phys, &
         llm_tmp, ttm_tmp, 1, 1, var_ana)
    CALL conf_dat2d(lon_ini, lat_ini, lon_rad, lat_rad, var_ana)
    CALL inter_barxy(lon_rad, lat_rad(:jml_phys -1), var_ana, rlonu(:iim), &
            rlatv, tmp_var)
    qsol_2d(:, :) = gr_int_dyn(tmp_var)

    CALL flinclo(fid_phys)

  END SUBROUTINE start_init_phys

END MODULE start_init_phys_m
