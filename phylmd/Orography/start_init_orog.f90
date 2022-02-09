MODULE start_init_orog_m

  ! From startvar.F, version 1.4
  ! 2006/01/27

  IMPLICIT NONE

CONTAINS

  SUBROUTINE start_init_orog(zmea, zstd, zsig, zgam, zthe, zpic, zval)

    ! Libraries:
    use jumble, only: pi
    use netcdf, only: nf90_nowrite
    use netcdf95, only: nf95_open, nf95_gw_var, nf95_inq_varid, nf95_close

    use conf_dat2d_m, only: conf_dat2d
    use dimensions, only: iim, jjm
    use dynetat0_m, only: rlatu, rlonv
    use grid_change, only: dyn_phy
    use grid_noro_m, only: grid_noro

    REAL, intent(out):: zmea(:) ! (klon) orographie moyenne

    REAL, intent(out):: zstd(:) ! (klon)
    ! (deviation standard de l'orographie sous-maille)

    REAL, intent(out):: zsig(:) ! (klon) pente de l'orographie sous-maille
    REAL, intent(out):: zgam(:) ! (klon) anisotropie de l'orographie sous maille

    REAL, intent(out):: zthe(:) ! (klon)
    ! (orientation de l'axe oriente dans la direction de plus grande
    ! pente de l'orographie sous maille)

    REAL, intent(out):: zpic(:) ! (klon) hauteur pics de la SSO
    REAL, intent(out):: zval(:) ! (klon) hauteur vallees de la SSO

    ! Local:

    INTEGER iml_rel
    INTEGER jml_rel
    INTEGER ncid, varid
    REAL, ALLOCATABLE:: relief(:, :) ! in m
    REAL, ALLOCATABLE:: lon_rad(:), lat_rad(:)
    REAL, ALLOCATABLE:: lon_ini(:), lat_ini(:)
    real, dimension(iim + 1, jjm + 1):: zmea_2d, zstd_2d, zsig_2d, zgam_2d
    real, dimension(iim + 1, jjm + 1):: zthe_2d, zpic_2d, zval_2d

    !-----------------------------------

    print *, "Call sequence information: start_init_orog"

    ! Read the high resolution orography:

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

    ALLOCATE(lon_rad(iml_rel), lat_rad(jml_rel))
    CALL conf_dat2d(lon_ini, lat_ini, lon_rad, lat_rad, relief, &
         interbar = .FALSE.)
    CALL grid_noro(lon_rad, lat_rad, relief, rlonv, rlatu, zmea_2d, zstd_2d, &
         zsig_2d, zgam_2d, zthe_2d, zpic_2d, zval_2d)
    zmea = pack(zmea_2d, dyn_phy)
    zstd = pack(zstd_2d, dyn_phy)
    zsig = pack(zsig_2d, dyn_phy)
    zgam = pack(zgam_2d, dyn_phy)
    zthe = pack(zthe_2d, dyn_phy)
    zpic = pack(zpic_2d, dyn_phy)
    zval = pack(zval_2d, dyn_phy)

  END SUBROUTINE start_init_orog

END MODULE start_init_orog_m
