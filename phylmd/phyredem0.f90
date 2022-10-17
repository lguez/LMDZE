module phyredem0_m

  IMPLICIT NONE

  INTEGER, protected:: ncid_restartphy

contains

  SUBROUTINE phyredem0(itau_phy_redem)

    ! From phylmd/phyredem.F, version 1.3, 2005/05/25 13:10:09
    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1993/08/18

    ! Objet : \'ecriture de l'\'etat de d\'emarrage ou red\'emarrage
    ! pour la physique

    USE netcdf, ONLY: nf90_clobber, nf90_global, nf90_float
    USE netcdf95, ONLY: nf95_create, nf95_put_att, nf95_def_dim, &
         nf95_def_var, nf95_enddef, nf95_put_var

    USE dimphy, ONLY: klev, klon
    USE dimsoil, ONLY: nsoilmx
    USE indicesol, ONLY: nbsrf
    use phyetat0_m, only: rlat, rlon

    integer, intent(in):: itau_phy_redem

    ! Local:

    INTEGER idim2, idim3, dimid_nbsrf, dimid_nsoilmx
    integer varid, varid_rlon, varid_rlat

    !------------------------------------------------------------

    PRINT *, 'Call sequence information: phyredem0'
    CALL nf95_create("restartphy.nc", nf90_clobber, ncid_restartphy)

    call nf95_put_att(ncid_restartphy, nf90_global, 'title', &
         'start file for the physics code')
    call nf95_put_att(ncid_restartphy, nf90_global, "itau_phy", itau_phy_redem)

    call nf95_def_dim(ncid_restartphy, 'points_physiques', klon, idim2)
    call nf95_def_dim(ncid_restartphy, 'klev', klev, idim3)
    call nf95_def_dim(ncid_restartphy, 'nbsrf', nbsrf, dimid_nbsrf)
    call nf95_def_dim(ncid_restartphy, 'nsoilmx', nsoilmx, dimid_nsoilmx)

    call nf95_def_var(ncid_restartphy, 'longitude', nf90_float, idim2, &
         varid_rlon)
    call nf95_def_var(ncid_restartphy, 'latitude', nf90_float, idim2, &
         varid_rlat)

    call nf95_def_var(ncid_restartphy, 'masque', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', 'masque terre mer')

    ! Fractions de chaque sous-surface

    call nf95_def_var(ncid_restartphy, 'pctsrf', nf90_float, &
         [idim2, dimid_nbsrf], varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', 'surface fraction')

    call nf95_def_var(ncid_restartphy, 'TS', nf90_float, &
         [idim2, dimid_nbsrf], varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', &
         'surface temperature')

    call nf95_def_var(ncid_restartphy, 'Tsoil', nf90_float, &
         [idim2, dimid_nsoilmx, dimid_nbsrf], varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', 'soil temperature')

    call nf95_def_var(ncid_restartphy, 'QS', nf90_float, &
         [idim2, dimid_nbsrf], varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', &
         'Humidite de surface')

    call nf95_def_var(ncid_restartphy, 'QSOL', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', &
         'Eau dans le sol (mm)')

    call nf95_def_var(ncid_restartphy, 'ALBE', nf90_float, &
         [idim2, dimid_nbsrf], varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', 'albedo de surface')

    call nf95_def_var(ncid_restartphy, 'SNOW', nf90_float, &
         [idim2, dimid_nbsrf], varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', 'Neige de surface')

    call nf95_def_var(ncid_restartphy, 'RADS', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', &
         'Rayonnement net a la surface')

    call nf95_def_var(ncid_restartphy, 'solsw', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', &
         'Rayonnement solaire a la surface')

    call nf95_def_var(ncid_restartphy, 'sollw', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', &
         'Rayonnement IR a la surface')

    call nf95_def_var(ncid_restartphy, 'fder', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', 'Derive de flux')

    call nf95_def_var(ncid_restartphy, 'rain_f', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', &
         'precipitation liquide')

    call nf95_def_var(ncid_restartphy, 'snow_f', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', &
         'precipitation solide')

    call nf95_def_var(ncid_restartphy, 'RUG', nf90_float, &
         [idim2, dimid_nbsrf], varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', &
         'rugosite de surface')

    call nf95_def_var(ncid_restartphy, 'AGESNO', nf90_float, &
         [idim2, dimid_nbsrf], varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', &
         'Age de la neige surface')

    call nf95_def_var(ncid_restartphy, 'ZMEA', nf90_float, idim2, varid)
    call nf95_def_var(ncid_restartphy, 'ZSTD', nf90_float, idim2, varid)
    call nf95_def_var(ncid_restartphy, 'ZSIG', nf90_float, idim2, varid)
    call nf95_def_var(ncid_restartphy, 'ZGAM', nf90_float, idim2, varid)
    call nf95_def_var(ncid_restartphy, 'ZTHE', nf90_float, idim2, varid)
    call nf95_def_var(ncid_restartphy, 'ZPIC', nf90_float, idim2, varid)
    call nf95_def_var(ncid_restartphy, 'ZVAL', nf90_float, idim2, varid)
    call nf95_def_var(ncid_restartphy, 'TANCIEN', nf90_float, &
         [idim2, idim3], varid)
    call nf95_def_var(ncid_restartphy, 'QANCIEN', nf90_float, &
         [idim2, idim3], varid)

    call nf95_def_var(ncid_restartphy, 'RUGMER', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', &
         'Longueur de rugosite sur mer')

    call nf95_def_var(ncid_restartphy, 'CLWCON', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', &
         'Eau liquide convective')

    call nf95_def_var(ncid_restartphy, 'RNEBCON', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', &
         'Nebulosite convective')

    call nf95_def_var(ncid_restartphy, 'RATQS', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', 'Ratqs')

    call nf95_def_var(ncid_restartphy, 'RUNOFFLIC0', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', 'Runofflic0')

    call nf95_def_var(ncid_restartphy, 'sig1', nf90_float, [idim2, idim3], &
         varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', &
         'section adiabatic updraft')

    call nf95_def_var(ncid_restartphy, 'w01', nf90_float, [idim2, idim3], &
         varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', &
         'vertical velocity within adiabatic updraft')

    call nf95_def_var(ncid_restartphy, 'trs', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', &
         'radon concentation in soil')

    call nf95_enddef(ncid_restartphy)

    call nf95_put_var(ncid_restartphy, varid_rlon, rlon)
    call nf95_put_var(ncid_restartphy, varid_rlat, rlat)

  END SUBROUTINE phyredem0

end module phyredem0_m
