MODULE start_init_orog_m

  ! From startvar.F, version 1.4
  ! 2006/01/27 15:14:22 Fairhead

  use dimensions, only: iim, jjm

  IMPLICIT NONE

  REAL, SAVE:: mask(iim + 1, jjm + 1) ! interpolated fraction of land

  private iim, jjm

CONTAINS

  SUBROUTINE start_init_orog(phis, zmea_2d, zstd_2d, zsig_2d, zgam_2d, &
       zthe_2d, zpic_2d, zval_2d)

    use conf_dat2d_m, only: conf_dat2d
    use dynetat0_m, only: rlatu, rlonv
    use grid_noro_m, only: grid_noro
    use indicesol, only: epsfra
    use netcdf, only: nf90_nowrite
    use netcdf95, only: nf95_open, nf95_gw_var, nf95_inq_varid, nf95_close
    use nr_util, only: pi, assert

    REAL, intent(out):: phis(:, :) ! (iim + 1, jjm + 1)
    ! surface geopotential, in m2 s-2

    REAL, intent(out):: zmea_2d(:, :) ! (iim + 1, jjm + 1) orographie moyenne

    REAL, intent(out):: zstd_2d(:, :) ! (iim + 1, jjm + 1)
    ! (deviation standard de l'orographie sous-maille)

    REAL, intent(out):: zsig_2d(:, :) ! (iim + 1, jjm + 1)
    ! (pente de l'orographie sous-maille)
    
    REAL, intent(out):: zgam_2d(:, :) ! (iim + 1, jjm + 1)
    ! (anisotropie de l'orographie sous maille)

    REAL, intent(out):: zthe_2d(:, :) ! (iim + 1, jjm + 1)
    ! (orientation de l'axe oriente dans la direction de plus grande
    ! pente de l'orographie sous maille)

    REAL, intent(out):: zpic_2d(:, :) ! (iim + 1, jjm + 1)
    ! hauteur pics de la SSO

    REAL, intent(out):: zval_2d(:, :) ! (iim + 1, jjm + 1)
    ! hauteur vallees de la SSO

    ! Local:

    INTEGER iml_rel
    INTEGER jml_rel
    INTEGER ncid, varid
    REAL, ALLOCATABLE:: relief(:, :) ! in m
    REAL, ALLOCATABLE:: lon_rad(:), lat_rad(:)
    REAL, ALLOCATABLE:: lon_ini(:), lat_ini(:)

    !-----------------------------------

    print *, "Call sequence information: start_init_orog"

    call assert((/size(phis, 1), size(zmea_2d, 1), size(zstd_2d, 1), &
         size(zsig_2d, 1), size(zgam_2d, 1), size(zthe_2d, 1), &
         size(zpic_2d, 1), size(zval_2d, 1)/) == iim + 1, "start_init_orog iim")
    call assert((/size(phis, 2), size(zmea_2d, 2), size(zstd_2d, 2), &
         size(zsig_2d, 2), size(zgam_2d, 2), size(zthe_2d, 2), &
         size(zpic_2d, 2), size(zval_2d, 2)/) == jjm + 1, "start_init_orog jjm")

    print *, 'Reading the high resolution orography...'

    call nf95_open('Relief.nc', nf90_nowrite, ncid)

    call nf95_inq_varid(ncid, "longitude", varid)
    call nf95_gw_var(ncid, varid, lon_ini)
    lon_ini = lon_ini * pi / 180. ! convert to rad
    iml_rel = size(lon_ini)

    call nf95_inq_varid(ncid, "latitude", varid)
    call nf95_gw_var(ncid, varid, lat_ini)
    lat_ini = lat_ini * pi / 180. ! convert to rad
    jml_rel = size(lat_ini)

    call nf95_inq_varid(ncid, "RELIEF", varid)
    call nf95_gw_var(ncid, varid, relief)

    call nf95_close(ncid)

    ALLOCATE(lon_rad(iml_rel))
    ALLOCATE(lat_rad(jml_rel))

    CALL conf_dat2d(lon_ini, lat_ini, lon_rad, lat_rad, relief, &
         interbar = .FALSE.)

    print *, 'Compute all the parameters needed for the gravity wave drag code'
    CALL grid_noro(lon_rad, lat_rad, relief, rlonv, rlatu, phis, zmea_2d, &
         zstd_2d, zsig_2d, zgam_2d, zthe_2d, zpic_2d, zval_2d, mask)
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
