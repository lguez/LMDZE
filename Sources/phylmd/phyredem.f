module phyredem_m

  IMPLICIT NONE

contains

  SUBROUTINE phyredem(fichnom, pctsrf, tsol, tsoil, tslab, seaice, qsurf, &
       qsol, snow, albedo, evap, rain_fall, snow_fall, solsw, sollw, &
       fder, radsol, frugs, agesno, zmea, zstd, zsig, zgam, zthe, zpic, zval, &
       t_ancien, q_ancien, rnebcon, ratqs, clwcon, run_off_lic_0, sig1, w01)

    ! From phylmd/phyredem.F, version 1.3, 2005/05/25 13:10:09
    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1993/08/18

    ! Objet : \'ecriture de l'\'etat de d\'emarrage ou red\'emarrage
    ! pour la physique

    USE dimphy, ONLY: klev, klon, zmasq
    USE dimsoil, ONLY: nsoilmx
    USE indicesol, ONLY: is_lic, is_oce, is_sic, is_ter, nbsrf
    USE netcdf, ONLY: nf90_clobber, nf90_global, nf90_float
    USE netcdf95, ONLY: nf95_create, nf95_put_att, nf95_def_dim, &
         nf95_def_var, nf95_enddef, nf95_put_var, nf95_close
    use phyetat0_m, only: rlat, rlon
    USE temps, ONLY: itau_phy

    CHARACTER(len=*), INTENT(IN):: fichnom
    REAL, INTENT(IN):: pctsrf(:, :) ! (klon, nbsrf)
    REAL, INTENT(IN):: tsol(:, :) ! (klon, nbsrf)
    REAL, INTENT(IN):: tsoil(:, :, :) ! (klon, nsoilmx, nbsrf)
    REAL, INTENT(IN):: tslab(:), seaice(:) ! (klon) slab ocean
    REAL, INTENT(IN):: qsurf(:, :) ! (klon, nbsrf)

    REAL, intent(in):: qsol(:) ! (klon)
    ! column-density of water in soil, in kg m-2

    REAL, INTENT(IN):: snow(klon, nbsrf)
    REAL, INTENT(IN):: albedo(klon, nbsrf)
    REAL, INTENT(IN):: evap(klon, nbsrf)
    REAL, INTENT(IN):: rain_fall(klon)
    REAL, INTENT(IN):: snow_fall(klon)
    REAL, INTENT(IN):: solsw(klon)
    REAL, INTENT(IN):: sollw(klon)
    REAL, INTENT(IN):: fder(klon)
    REAL, INTENT(IN):: radsol(klon)
    REAL, INTENT(IN):: frugs(klon, nbsrf)
    REAL, INTENT(IN):: agesno(klon, nbsrf)
    REAL, INTENT(IN):: zmea(klon)
    REAL, intent(in):: zstd(klon)
    REAL, intent(in):: zsig(klon)
    REAL, intent(in):: zgam(klon)
    REAL, intent(in):: zthe(klon)
    REAL, intent(in):: zpic(klon)
    REAL, intent(in):: zval(klon)
    REAL, intent(in):: t_ancien(klon, klev), q_ancien(klon, klev)
    REAL, intent(in):: rnebcon(klon, klev), ratqs(klon, klev)
    REAL, intent(in):: clwcon(klon, klev)
    REAL, intent(in):: run_off_lic_0(klon)
    real, intent(in):: sig1(klon, klev) ! section adiabatic updraft

    real, intent(in):: w01(klon, klev) 
    ! vertical velocity within adiabatic updraft

    ! Local:

    INTEGER ncid, idim2, idim3, dimid_nbsrf, dimid_nsoilmx
    integer varid, varid_run_off_lic_0, varid_sig1, varid_w01, varid_rlon
    integer varid_rlat, varid_zmasq, varid_fter, varid_flic, varid_foce
    integer varid_fsic, varid_ts, varid_tsoil, varid_tslab, varid_SEAICE
    integer varid_qs, varid_QSOL, varid_albe, varid_evap, varid_snow
    integer varid_rads, varid_solsw, varid_sollw, varid_fder, varid_rain_f
    integer varid_snow_f, varid_rug, varid_agesno, varid_zmea, varid_zstd
    integer varid_zsig, varid_zgam, varid_zthe, varid_zpic, varid_zval
    integer varid_tancien, varid_qancien, varid_rugmer, varid_clwcon
    integer varid_rnebcon, varid_ratqs

    !------------------------------------------------------------

    PRINT *, 'Call sequence information: phyredem'
    CALL nf95_create(fichnom, nf90_clobber, ncid)

    call nf95_put_att(ncid, nf90_global, 'title', &
         'start file for the physics code')
    call nf95_put_att(ncid, nf90_global, "itau_phy", itau_phy)

    call nf95_def_dim(ncid, 'points_physiques', klon, idim2)
    call nf95_def_dim(ncid, 'klev', klev, idim3)
    call nf95_def_dim(ncid, 'nbsrf', nbsrf, dimid_nbsrf)
    call nf95_def_dim(ncid, 'nsoilmx', nsoilmx, dimid_nsoilmx)

    call nf95_def_var(ncid, 'longitude', nf90_float, idim2, varid_rlon)
    call nf95_def_var(ncid, 'latitude', nf90_float, idim2, varid_rlat)

    call nf95_def_var(ncid, 'masque', nf90_float, idim2, varid_zmasq)
    call nf95_put_att(ncid, varid_zmasq, 'title', 'masque terre mer')

    ! Fractions de chaque sous-surface

    call nf95_def_var(ncid, 'FTER', nf90_float, idim2, varid_fter)
    call nf95_put_att(ncid, varid_fter, 'title', 'fraction de continent')

    call nf95_def_var(ncid, 'FLIC', nf90_float, idim2, varid_flic)
    call nf95_put_att(ncid, varid_flic, 'title', 'fraction glace de terre')

    call nf95_def_var(ncid, 'FOCE', nf90_float, idim2, varid_foce)
    call nf95_put_att(ncid, varid_foce, 'title', 'fraction ocean')

    call nf95_def_var(ncid, 'FSIC', nf90_float, idim2, varid_fsic)
    call nf95_put_att(ncid, varid_fsic, 'title', 'fraction glace mer')

    call nf95_def_var(ncid, 'TS', nf90_float, (/idim2, dimid_nbsrf/), varid_ts)
    call nf95_put_att(ncid, varid_ts, 'title', 'surface temperature')

    call nf95_def_var(ncid, 'Tsoil', nf90_float, (/idim2, dimid_nsoilmx, &
         dimid_nbsrf/), varid_tsoil)
    call nf95_put_att(ncid, varid_tsoil, 'title', 'soil temperature')

    call nf95_def_var(ncid, 'TSLAB', nf90_float, idim2, varid_tslab)
    call nf95_put_att(ncid, varid_tslab, 'title', &
         'Ecart de la SST (pour slab-ocean)')

    call nf95_def_var(ncid, 'SEAICE', nf90_float, idim2, varid_SEAICE)
    call nf95_put_att(ncid, varid_SEAICE, 'title', &
         'Glace de mer kg/m2 (pour slab-ocean)')

    call nf95_def_var(ncid, 'QS', nf90_float, (/idim2, dimid_nbsrf/), varid_qs)
    call nf95_put_att(ncid, varid, 'title', 'Humidite de surface')

    call nf95_def_var(ncid, 'QSOL', nf90_float, idim2, varid_QSOL)
    call nf95_put_att(ncid, varid_QSOL, 'title', 'Eau dans le sol (mm)')

    call nf95_def_var(ncid, 'ALBE', nf90_float, (/idim2, dimid_nbsrf/), &
         varid_albe)
    call nf95_put_att(ncid, varid_albe, 'title', 'albedo de surface')

    call nf95_def_var(ncid, 'EVAP', nf90_float, (/idim2, dimid_nbsrf/), &
         varid_evap)
    call nf95_put_att(ncid, varid_evap, 'title', 'Evaporation de surface')

    call nf95_def_var(ncid, 'SNOW', nf90_float, (/idim2, dimid_nbsrf/), &
         varid_snow)
    call nf95_put_att(ncid, varid_snow, 'title', 'Neige de surface')

    call nf95_def_var(ncid, 'RADS', nf90_float, idim2, varid_rads)
    call nf95_put_att(ncid, varid_rads, 'title', &
         'Rayonnement net a la surface')

    call nf95_def_var(ncid, 'solsw', nf90_float, idim2, varid_solsw)
    call nf95_put_att(ncid, varid_solsw, 'title', &
         'Rayonnement solaire a la surface')

    call nf95_def_var(ncid, 'sollw', nf90_float, idim2, varid_sollw)
    call nf95_put_att(ncid, varid_sollw, 'title', &
         'Rayonnement IF a la surface')

    call nf95_def_var(ncid, 'fder', nf90_float, idim2, varid_fder)
    call nf95_put_att(ncid, varid_fder, 'title', 'Derive de flux')

    call nf95_def_var(ncid, 'rain_f', nf90_float, idim2, varid_rain_f)
    call nf95_put_att(ncid, varid_rain_f, 'title', 'precipitation liquide')

    call nf95_def_var(ncid, 'snow_f', nf90_float, idim2, varid_snow_f)
    call nf95_put_att(ncid, varid_snow_f, 'title', 'precipitation solide')

    call nf95_def_var(ncid, 'RUG', nf90_float, (/idim2, dimid_nbsrf/), &
         varid_rug)
    call nf95_put_att(ncid, varid_rug, 'title', 'rugosite de surface')

    call nf95_def_var(ncid, 'AGESNO', nf90_float, (/idim2, dimid_nbsrf/), &
         varid_agesno)
    call nf95_put_att(ncid, varid_agesno, 'title', 'Age de la neige surface')

    call nf95_def_var(ncid, 'ZMEA', nf90_float, idim2, varid_zmea)
    call nf95_def_var(ncid, 'ZSTD', nf90_float, idim2, varid_zstd)
    call nf95_def_var(ncid, 'ZSIG', nf90_float, idim2, varid_zsig)
    call nf95_def_var(ncid, 'ZGAM', nf90_float, idim2, varid_zgam)
    call nf95_def_var(ncid, 'ZTHE', nf90_float, idim2, varid_zthe)
    call nf95_def_var(ncid, 'ZPIC', nf90_float, idim2, varid_zpic)
    call nf95_def_var(ncid, 'ZVAL', nf90_float, idim2, varid_zval)
    call nf95_def_var(ncid, 'TANCIEN', nf90_float, (/idim2, idim3/), &
         varid_tancien)
    call nf95_def_var(ncid, 'QANCIEN', nf90_float, (/idim2, idim3/), &
         varid_qancien)

    call nf95_def_var(ncid, 'RUGMER', nf90_float, idim2, varid_rugmer)
    call nf95_put_att(ncid, varid_rugmer, 'title', &
         'Longueur de rugosite sur mer')

    call nf95_def_var(ncid, 'CLWCON', nf90_float, idim2, varid_clwcon)
    call nf95_put_att(ncid, varid_clwcon, 'title', 'Eau liquide convective')

    call nf95_def_var(ncid, 'RNEBCON', nf90_float, idim2, varid_rnebcon)
    call nf95_put_att(ncid, varid_rnebcon, 'title', 'Nebulosite convective')

    call nf95_def_var(ncid, 'RATQS', nf90_float, idim2, varid_ratqs)
    call nf95_put_att(ncid, varid_ratqs, 'title', 'Ratqs')

    call nf95_def_var(ncid, 'RUNOFFLIC0', nf90_float, idim2, &
         varid_run_off_lic_0)
    call nf95_put_att(ncid, varid_run_off_lic_0, 'title', 'Runofflic0')

    call nf95_def_var(ncid, 'sig1', nf90_float, (/idim2, idim3/), varid_sig1)
    call nf95_put_att(ncid, varid_sig1, 'long_name', &
         'section adiabatic updraft')

    call nf95_def_var(ncid, 'w01', nf90_float, (/idim2, idim3/), varid_w01)
    call nf95_put_att(ncid, varid_w01, 'long_name', &
         'vertical velocity within adiabatic updraft')

    call nf95_enddef(ncid)

    call nf95_put_var(ncid, varid_rlon, rlon)
    call nf95_put_var(ncid, varid_rlat, rlat)
    call nf95_put_var(ncid, varid_zmasq, zmasq)
    call nf95_put_var(ncid, varid_fter, pctsrf(:, is_ter))
    call nf95_put_var(ncid, varid_flic, pctsrf(:, is_lic))
    call nf95_put_var(ncid, varid_foce, pctsrf(:, is_oce))
    call nf95_put_var(ncid, varid_fsic, pctsrf(:, is_sic))
    call nf95_put_var(ncid, varid_ts, tsol)
    call nf95_put_var(ncid, varid_tsoil, tsoil)
    call nf95_put_var(ncid, varid_tslab, tslab)
    call nf95_put_var(ncid, varid_SEAICE, seaice)
    call nf95_put_var(ncid, varid_qs, qsurf)
    call nf95_put_var(ncid, varid_QSOL, qsol)
    call nf95_put_var(ncid, varid_albe, albedo)
    call nf95_put_var(ncid, varid_evap, evap)
    call nf95_put_var(ncid, varid_snow, snow)
    call nf95_put_var(ncid, varid_rads, radsol)
    call nf95_put_var(ncid, varid_solsw, solsw)
    call nf95_put_var(ncid, varid_sollw, sollw)
    call nf95_put_var(ncid, varid_fder, fder)
    call nf95_put_var(ncid, varid_rain_f, rain_fall)
    call nf95_put_var(ncid, varid_snow_f, snow_fall)
    call nf95_put_var(ncid, varid_rug, frugs)
    call nf95_put_var(ncid, varid_agesno, agesno)
    call nf95_put_var(ncid, varid_zmea, zmea)
    call nf95_put_var(ncid, varid_zstd, zstd)
    call nf95_put_var(ncid, varid_zsig, zsig)
    call nf95_put_var(ncid, varid_zgam, zgam)
    call nf95_put_var(ncid, varid_zthe, zthe)
    call nf95_put_var(ncid, varid_zpic, zpic)
    call nf95_put_var(ncid, varid_zval, zval)
    call nf95_put_var(ncid, varid_tancien, t_ancien)
    call nf95_put_var(ncid, varid_qancien, q_ancien)
    call nf95_put_var(ncid, varid_rugmer, frugs(:, is_oce))
    call nf95_put_var(ncid, varid_clwcon, clwcon(:, 1))
    call nf95_put_var(ncid, varid_rnebcon, rnebcon(:, 1))
    call nf95_put_var(ncid, varid_ratqs, ratqs(:, 1))
    call nf95_put_var(ncid, varid_run_off_lic_0, run_off_lic_0)
    call nf95_put_var(ncid, varid_sig1, sig1)
    call nf95_put_var(ncid, varid_w01, w01)

    call nf95_close(ncid)

  END SUBROUTINE phyredem

end module phyredem_m
