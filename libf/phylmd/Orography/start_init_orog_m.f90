MODULE start_init_orog_m

  ! From startvar.F, version 1.4
  ! 2006/01/27 15:14:22 Fairhead

  IMPLICIT NONE

  REAL, ALLOCATABLE, SAVE:: mask(:, :) ! fraction of land (iim + 1, jjm + 1)
  REAL, ALLOCATABLE, SAVE:: phis(:, :) ! surface geopotential, in m2 s-2

CONTAINS

  SUBROUTINE start_init_orog(relief, zstd_2d, zsig_2d, zgam_2d, zthe_2d, &
       zpic_2d, zval_2d)

    USE flincom, only: flininfo, flinopen_nozoom, flinclo
    use flinget_m, only: flinget
    use conf_dat2d_m, only: conf_dat2d
    use comgeom, only: rlatu, rlonv
    use dimens_m, only: iim, jjm
    use indicesol, only: epsfra
    use comconst, only: pi
    use grid_noro_m, only: grid_noro

    REAL, intent(out):: relief(:, :) ! orographie moyenne

    REAL, intent(out):: zstd_2d(:, :)
    ! (deviation standard de l'orographie sous-maille)

    REAL, intent(out):: zsig_2d(:, :)
    ! (pente de l'orographie sous-maille)
    
    REAL, intent(out):: zgam_2d(:, :)
    ! (anisotropie de l'orographie sous maille)

    REAL, intent(out):: zthe_2d(:, :)
    ! (orientation de l'axe oriente dans la direction de plus grande
    ! pente de l'orographie sous maille)

    REAL, intent(out):: zpic_2d(:, :) ! hauteur pics de la SSO
    REAL, intent(out):: zval_2d(:, :) ! hauteur vallees de la SSO

    ! Local:

    INTEGER, SAVE:: iml_rel
    INTEGER, SAVE:: jml_rel
    REAL lev(1), date, dt
    INTEGER itau(1), fid
    INTEGER  llm_tmp, ttm_tmp
    REAL, ALLOCATABLE:: relief_hi(:, :)
    REAL, ALLOCATABLE:: lon_rad(:), lat_rad(:)
    REAL, ALLOCATABLE:: lon_ini(:), lat_ini(:)
    REAL, ALLOCATABLE:: lon_rel(:, :), lat_rel(:, :)

    CHARACTER(len=120) orogfname

    !-----------------------------------

    print *, "Call sequence information: start_init_orog"

    if (any((/size(relief, 1), size(zstd_2d, 1), size(zsig_2d, 1), &
         size(zgam_2d, 1), size(zthe_2d, 1), size(zpic_2d, 1), &
         size(zval_2d, 1)/) /= iim + 1)) stop "start_init_orog size 1"
    if (any((/size(relief, 2), size(zstd_2d, 2), size(zsig_2d, 2), &
         size(zgam_2d, 2), size(zthe_2d, 2), size(zpic_2d, 2), &
         size(zval_2d, 2)/) /= jjm + 1)) stop "start_init_orog size 2"

    orogfname = 'Relief.nc'
    print *, 'Reading the high resolution orography'

    CALL flininfo(orogfname, iml_rel, jml_rel, llm_tmp, ttm_tmp, fid)

    ALLOCATE(lat_rel(iml_rel, jml_rel))
    ALLOCATE(lon_rel(iml_rel, jml_rel))
    ALLOCATE(relief_hi(iml_rel, jml_rel))

    CALL flinopen_nozoom(iml_rel, jml_rel, llm_tmp, &
         lon_rel, lat_rel, lev, ttm_tmp, itau, date, dt, fid)
    ! 'RELIEF': high resolution orography 
    CALL flinget(fid, 'RELIEF', iml_rel, jml_rel, llm_tmp, ttm_tmp, 1, 1, &
         relief_hi)
    CALL flinclo(fid)

    ! In case we have a file which is in degrees we do the transformation:

    ALLOCATE(lon_rad(iml_rel))
    ALLOCATE(lon_ini(iml_rel))

    IF (MAXVAL(lon_rel(:, :)) > pi) THEN
       lon_ini(:) = lon_rel(:, 1) * pi / 180.
    ELSE
       lon_ini(:) = lon_rel(:, 1) 
    ENDIF

    ALLOCATE(lat_rad(jml_rel))
    ALLOCATE(lat_ini(jml_rel))

    IF (MAXVAL(lat_rel(:, :)) > pi) THEN
       lat_ini(:) = lat_rel(1, :) * pi / 180.
    ELSE
       lat_ini(:) = lat_rel(1, :) 
    ENDIF

    CALL conf_dat2d(lon_ini, lat_ini, lon_rad, lat_rad, relief_hi , &
         interbar=.FALSE.)

    print *, 'Compute all the parameters needed for the gravity wave drag code'

    ! Allocate the data we need to put in the interpolated fields:
    ALLOCATE(phis(iim + 1, jjm + 1))
    ALLOCATE(mask(iim + 1, jjm + 1))

    CALL grid_noro(lon_rad, lat_rad, relief_hi, rlonv, rlatu, phis, relief, &
         zstd_2d, zsig_2d, zgam_2d, zthe_2d, zpic_2d, zval_2d, mask)

    phis(iim + 1, :) = phis(1, :)
    phis(:, :) = phis(:, :) * 9.81

    mask(2:, 1) = mask(1, 1) ! north pole
    mask(2:, jjm + 1) = mask(1, jjm + 1) ! south pole
    mask(iim + 1, 2:jjm) = mask(1, 2:jjm) ! Greenwich
    WHERE (mask < EPSFRA)
       mask = 0.
    elsewhere (1. - mask < EPSFRA)
       mask = 1.
    endwhere

  END SUBROUTINE start_init_orog

END MODULE start_init_orog_m
