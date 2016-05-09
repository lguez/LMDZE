module phyredem0_m

  IMPLICIT NONE

  INTEGER ncid_restartphy

contains

  SUBROUTINE phyredem0(lmt_pas)

    ! From phylmd/phyredem.F, version 1.3, 2005/05/25 13:10:09
    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1993/08/18

    ! Objet : \'ecriture de l'\'etat de d\'emarrage ou red\'emarrage
    ! pour la physique

    use conf_gcm_m, only: nday
    USE dimphy, ONLY: klev, klon
    USE dimsoil, ONLY: nsoilmx
    USE indicesol, ONLY: nbsrf
    USE netcdf, ONLY: nf90_clobber, nf90_global, nf90_float
    USE netcdf95, ONLY: nf95_create, nf95_put_att, nf95_def_dim, &
         nf95_def_var, nf95_enddef, nf95_put_var
    use phyetat0_m, only: rlat, rlon, itau_phy

    INTEGER, intent(in):: lmt_pas ! number of time steps of "physics" per day

    ! Local:

    INTEGER idim2, idim3, dimid_nbsrf, dimid_nsoilmx
    integer varid, varid_rlon, varid_rlat

    !------------------------------------------------------------

    PRINT *, 'Call sequence information: phyredem0'
    CALL nf95_create("restartphy.nc", nf90_clobber, ncid_restartphy)

    call nf95_put_att(ncid_restartphy, nf90_global, 'title', &
         'start file for the physics code')
    call nf95_put_att(ncid_restartphy, nf90_global, "itau_phy", &
         itau_phy + nday * lmt_pas)

    call nf95_def_dim(ncid_restartphy, 'points_physiques', klon, idim2)
    call nf95_def_dim(ncid_restartphy, 'klev', klev, idim3)
    call nf95_def_dim(ncid_restartphy, 'nbsrf', nbsrf, dimid_nbsrf)
    call nf95_def_dim(ncid_restartphy, 'nsoilmx', nsoilmx, dimid_nsoilmx)

    call nf95_def_var(ncid_restartphy, 'longitude', nf90_float, idim2, &
         varid_rlon)
    call nf95_def_var(ncid_restartphy, 'latitude', nf90_float, idim2, &
         varid_rlat)

    call nf95_def_var(ncid_restartphy, 'masque', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'masque terre mer')

    ! Fractions de chaque sous-surface

    call nf95_def_var(ncid_restartphy, 'FTER', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'fraction de continent')

    call nf95_def_var(ncid_restartphy, 'FLIC', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', &
         'fraction glace de terre')

    call nf95_def_var(ncid_restartphy, 'FOCE', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'fraction ocean')

    call nf95_def_var(ncid_restartphy, 'FSIC', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'fraction glace mer')

    call nf95_def_var(ncid_restartphy, 'TS', nf90_float, &
         (/idim2, dimid_nbsrf/), varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'surface temperature')

    call nf95_def_var(ncid_restartphy, 'Tsoil', nf90_float, &
         (/idim2, dimid_nsoilmx, dimid_nbsrf/), varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'soil temperature')

    call nf95_def_var(ncid_restartphy, 'QS', nf90_float, &
         (/idim2, dimid_nbsrf/), varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'Humidite de surface')

    call nf95_def_var(ncid_restartphy, 'QSOL', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'Eau dans le sol (mm)')

    call nf95_def_var(ncid_restartphy, 'ALBE', nf90_float, &
         (/idim2, dimid_nbsrf/), varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'albedo de surface')

    call nf95_def_var(ncid_restartphy, 'EVAP', nf90_float, &
         (/idim2, dimid_nbsrf/), varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'Evaporation de surface')

    call nf95_def_var(ncid_restartphy, 'SNOW', nf90_float, &
         (/idim2, dimid_nbsrf/), varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'Neige de surface')

    call nf95_def_var(ncid_restartphy, 'RADS', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', &
         'Rayonnement net a la surface')

    call nf95_def_var(ncid_restartphy, 'solsw', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', &
         'Rayonnement solaire a la surface')

    call nf95_def_var(ncid_restartphy, 'sollw', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', &
         'Rayonnement IF a la surface')

    call nf95_def_var(ncid_restartphy, 'fder', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'Derive de flux')

    call nf95_def_var(ncid_restartphy, 'rain_f', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'precipitation liquide')

    call nf95_def_var(ncid_restartphy, 'snow_f', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'precipitation solide')

    call nf95_def_var(ncid_restartphy, 'RUG', nf90_float, &
         (/idim2, dimid_nbsrf/), varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'rugosite de surface')

    call nf95_def_var(ncid_restartphy, 'AGESNO', nf90_float, &
         (/idim2, dimid_nbsrf/), varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', &
         'Age de la neige surface')

    call nf95_def_var(ncid_restartphy, 'ZMEA', nf90_float, idim2, varid)
    call nf95_def_var(ncid_restartphy, 'ZSTD', nf90_float, idim2, varid)
    call nf95_def_var(ncid_restartphy, 'ZSIG', nf90_float, idim2, varid)
    call nf95_def_var(ncid_restartphy, 'ZGAM', nf90_float, idim2, varid)
    call nf95_def_var(ncid_restartphy, 'ZTHE', nf90_float, idim2, varid)
    call nf95_def_var(ncid_restartphy, 'ZPIC', nf90_float, idim2, varid)
    call nf95_def_var(ncid_restartphy, 'ZVAL', nf90_float, idim2, varid)
    call nf95_def_var(ncid_restartphy, 'TANCIEN', nf90_float, &
         (/idim2, idim3/), varid)
    call nf95_def_var(ncid_restartphy, 'QANCIEN', nf90_float, &
         (/idim2, idim3/), varid)

    call nf95_def_var(ncid_restartphy, 'RUGMER', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', &
         'Longueur de rugosite sur mer')

    call nf95_def_var(ncid_restartphy, 'CLWCON', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'Eau liquide convective')

    call nf95_def_var(ncid_restartphy, 'RNEBCON', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'Nebulosite convective')

    call nf95_def_var(ncid_restartphy, 'RATQS', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'Ratqs')

    call nf95_def_var(ncid_restartphy, 'RUNOFFLIC0', nf90_float, idim2, varid)
    call nf95_put_att(ncid_restartphy, varid, 'title', 'Runofflic0')

    call nf95_def_var(ncid_restartphy, 'sig1', nf90_float, (/idim2, idim3/), &
         varid)
    call nf95_put_att(ncid_restartphy, varid, 'long_name', &
         'section adiabatic updraft')

    call nf95_def_var(ncid_restartphy, 'w01', nf90_float, (/idim2, idim3/), &
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
