module grilles_gcm_netcdf_sub_m

  IMPLICIT NONE

contains

  SUBROUTINE grilles_gcm_netcdf_sub(phis)

    ! This subroutine creates the file "grilles_gcm.nc" containg
    ! longitudes and latitudes in degrees for grids u and v.

    USE comconst, ONLY: g
    USE comgeom, ONLY: aire_2d
    USE dimens_m, ONLY: iim, jjm
    use dynetat0_m, only: rlatu, rlatv, rlonu, rlonv
    USE netcdf, ONLY: nf90_clobber, nf90_float, nf90_int
    USE netcdf95, ONLY: nf95_close, nf95_create, nf95_def_dim, nf95_def_var, &
         nf95_enddef, nf95_put_att, nf95_put_var
    USE nr_util, ONLY: pi
    USE start_init_orog_m, ONLY: mask

    REAL, INTENT(IN):: phis(:, :) ! (iim + 1, jjm + 1) 
    ! surface geopotential, in m2 s-2

    ! Local:

    REAL grille_s(iim + 1, jjm + 1)
    INTEGER i, j

    ! For Netcdf:
    INTEGER ncid
    INTEGER dimid_lonu, dimid_lonv, dimid_latu, dimid_latv
    INTEGER varid_lonu, varid_lonv, varid_latu, varid_latv
    INTEGER varid
    INTEGER varid_orog, varid_area, varid_mask

    !-----------------------------------------------------------------

    print *, "Call sequence information: grilles_gcm_netcdf_sub"

    call NF95_CREATE('grilles_gcm.nc', NF90_CLOBBER, ncid)

    call NF95_DEF_DIM(ncid, 'lonu', iim + 1, dimid_lonu)
    call NF95_DEF_DIM(ncid, 'lonv', iim + 1, dimid_lonv)
    call NF95_DEF_DIM(ncid, 'latu', jjm + 1, dimid_latu)
    call NF95_DEF_DIM(ncid, 'latv', jjm, dimid_latv)

    call NF95_DEF_VAR(ncid, 'lonu', NF90_FLOAT, dimid_lonu, varid_lonu)
    call NF95_PUT_ATT(ncid, varid_lonu, 'units', 'degrees_east')
    call NF95_PUT_ATT(ncid, varid_lonu, 'long_name', 'Longitude en u')

    call NF95_DEF_VAR(ncid, 'lonv', NF90_FLOAT, dimid_lonv, varid_lonv)
    call NF95_PUT_ATT(ncid, varid_lonv, 'units', 'degrees_east')
    call NF95_PUT_ATT(ncid, varid_lonv, 'long_name', 'Longitude en v')

    call NF95_DEF_VAR(ncid, 'latu', NF90_FLOAT, dimid_latu, varid_latu)
    call NF95_PUT_ATT(ncid, varid_latu, 'units', 'degrees_north')
    call NF95_PUT_ATT(ncid, varid_latu, 'long_name', 'Latitude en u')

    call NF95_DEF_VAR(ncid, 'latv', NF90_FLOAT, dimid_latv, varid_latv)
    call NF95_PUT_ATT(ncid, varid_latv, 'units', 'degrees_north')
    call NF95_PUT_ATT(ncid, varid_latv, 'long_name', 'Latitude en v')

    call NF95_DEF_VAR(ncid, 'grille_u', NF90_FLOAT, (/dimid_lonu, &
         dimid_latu/), varid)
    call NF95_PUT_ATT(ncid, varid, 'long_name', 'Grille aux points u')

    call NF95_DEF_VAR(ncid, 'grille_v', NF90_FLOAT, (/dimid_lonv, &
         dimid_latv/), varid)
    call NF95_PUT_ATT(ncid, varid, 'long_name', 'Grille aux points v')

    call NF95_DEF_VAR(ncid, 'grille_s', NF90_FLOAT, (/dimid_lonv, &
         dimid_latu/), varid)
    call NF95_PUT_ATT(ncid, varid, 'long_name', 'Grille aux points scalaires')

    call NF95_DEF_VAR(ncid, 'orog', NF90_FLOAT, (/dimid_lonv, dimid_latu/), &
         varid_orog)
    call NF95_PUT_ATT(ncid, varid_orog, 'units', 'm')
    call NF95_PUT_ATT(ncid, varid_orog, 'standard_name', 'surface_altitude')

    call NF95_DEF_VAR(ncid, 'aire', NF90_FLOAT, (/dimid_lonv, dimid_latu/), &
         varid_area)
    call NF95_PUT_ATT(ncid, varid_area, 'units', 'm2')

    call NF95_DEF_VAR(ncid, 'mask', NF90_INT, (/dimid_lonv, dimid_latu/), &
         varid_mask)

    call NF95_ENDDEF(ncid)

    call NF95_PUT_VAR(ncid, varid_lonu, rlonu * 180. / pi)
    call NF95_PUT_VAR(ncid, varid_lonv, rlonv * 180. / pi)
    call NF95_PUT_VAR(ncid, varid_latu, rlatu * 180. / pi)
    call NF95_PUT_VAR(ncid, varid_latv, rlatv * 180. / pi)

    forall (i = 1:iim + 1, j = 1:jjm + 1) grille_s(i, j) = MOD(i, 2) + MOD(j, 2)
    call NF95_PUT_VAR(ncid, varid, grille_s)

    call NF95_PUT_VAR(ncid, varid_orog, phis / g)
    call NF95_PUT_VAR(ncid, varid_area, aire_2d)
    call NF95_PUT_VAR(ncid, varid_mask, nINT(mask))

    CALL nf95_close(ncid)

  END SUBROUTINE grilles_gcm_netcdf_sub

end module grilles_gcm_netcdf_sub_m
