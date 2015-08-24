module phyredem_m

  IMPLICIT NONE

contains

  SUBROUTINE phyredem(pctsrf, tsol, tsoil, tslab, seaice, qsurf, qsol, snow, &
       albedo, evap, rain_fall, snow_fall, solsw, sollw, fder, radsol, frugs, &
       agesno, zmea, zstd, zsig, zgam, zthe, zpic, zval, t_ancien, q_ancien, &
       rnebcon, ratqs, clwcon, run_off_lic_0, sig1, w01)

    ! From phylmd/phyredem.F, version 1.3, 2005/05/25 13:10:09
    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1993/08/18

    ! Objet : \'ecriture de l'\'etat de d\'emarrage ou red\'emarrage
    ! pour la physique

    USE dimphy, ONLY: klev, klon, zmasq
    USE indicesol, ONLY: is_lic, is_oce, is_sic, is_ter, nbsrf
    USE netcdf95, ONLY: nf95_inq_varid, nf95_put_var, nf95_close
    use phyredem0_m, only: ncid_restartphy

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
    integer varid

    !------------------------------------------------------------

    PRINT *, 'Call sequence information: phyredem'

    call nf95_inq_varid(ncid_restartphy, "masque", varid)
    call nf95_put_var(ncid_restartphy, varid, zmasq)

    call nf95_inq_varid(ncid_restartphy, "FTER", varid)
    call nf95_put_var(ncid_restartphy, varid, pctsrf(:, is_ter))

    call nf95_inq_varid(ncid_restartphy, "FLIC", varid)
    call nf95_put_var(ncid_restartphy, varid, pctsrf(:, is_lic))

    call nf95_inq_varid(ncid_restartphy, "FOCE", varid)
    call nf95_put_var(ncid_restartphy, varid, pctsrf(:, is_oce))

    call nf95_inq_varid(ncid_restartphy, "FSIC", varid)
    call nf95_put_var(ncid_restartphy, varid, pctsrf(:, is_sic))

    call nf95_inq_varid(ncid_restartphy, "TS", varid)
    call nf95_put_var(ncid_restartphy, varid, tsol)

    call nf95_inq_varid(ncid_restartphy, "Tsoil", varid)
    call nf95_put_var(ncid_restartphy, varid, tsoil)

    call nf95_inq_varid(ncid_restartphy, "TSLAB", varid)
    call nf95_put_var(ncid_restartphy, varid, tslab)

    call nf95_inq_varid(ncid_restartphy, "SEAICE", varid)
    call nf95_put_var(ncid_restartphy, varid, seaice)

    call nf95_inq_varid(ncid_restartphy, "QS", varid)
    call nf95_put_var(ncid_restartphy, varid, qsurf)

    call nf95_inq_varid(ncid_restartphy, "QSOL", varid)
    call nf95_put_var(ncid_restartphy, varid, qsol)

    call nf95_inq_varid(ncid_restartphy, "ALBE", varid)
    call nf95_put_var(ncid_restartphy, varid, albedo)

    call nf95_inq_varid(ncid_restartphy, "EVAP", varid)
    call nf95_put_var(ncid_restartphy, varid, evap)

    call nf95_inq_varid(ncid_restartphy, "SNOW", varid)
    call nf95_put_var(ncid_restartphy, varid, snow)

    call nf95_inq_varid(ncid_restartphy, "RADS", varid)
    call nf95_put_var(ncid_restartphy, varid, radsol)

    call nf95_inq_varid(ncid_restartphy, "solsw", varid)
    call nf95_put_var(ncid_restartphy, varid, solsw)

    call nf95_inq_varid(ncid_restartphy, "sollw", varid)
    call nf95_put_var(ncid_restartphy, varid, sollw)

    call nf95_inq_varid(ncid_restartphy, "fder", varid)
    call nf95_put_var(ncid_restartphy, varid, fder)

    call nf95_inq_varid(ncid_restartphy, "rain_f", varid)
    call nf95_put_var(ncid_restartphy, varid, rain_fall)

    call nf95_inq_varid(ncid_restartphy, "snow_f", varid)
    call nf95_put_var(ncid_restartphy, varid, snow_fall)

    call nf95_inq_varid(ncid_restartphy, "RUG", varid)
    call nf95_put_var(ncid_restartphy, varid, frugs)

    call nf95_inq_varid(ncid_restartphy, "AGESNO", varid)
    call nf95_put_var(ncid_restartphy, varid, agesno)

    call nf95_inq_varid(ncid_restartphy, "ZMEA", varid)
    call nf95_put_var(ncid_restartphy, varid, zmea)

    call nf95_inq_varid(ncid_restartphy, "ZSTD", varid)
    call nf95_put_var(ncid_restartphy, varid, zstd)

    call nf95_inq_varid(ncid_restartphy, "ZSIG", varid)
    call nf95_put_var(ncid_restartphy, varid, zsig)

    call nf95_inq_varid(ncid_restartphy, "ZGAM", varid)
    call nf95_put_var(ncid_restartphy, varid, zgam)

    call nf95_inq_varid(ncid_restartphy, "ZTHE", varid)
    call nf95_put_var(ncid_restartphy, varid, zthe)

    call nf95_inq_varid(ncid_restartphy, "ZPIC", varid)
    call nf95_put_var(ncid_restartphy, varid, zpic)

    call nf95_inq_varid(ncid_restartphy, "ZVAL", varid)
    call nf95_put_var(ncid_restartphy, varid, zval)

    call nf95_inq_varid(ncid_restartphy, "TANCIEN", varid)
    call nf95_put_var(ncid_restartphy, varid, t_ancien)

    call nf95_inq_varid(ncid_restartphy, "QANCIEN", varid)
    call nf95_put_var(ncid_restartphy, varid, q_ancien)

    call nf95_inq_varid(ncid_restartphy, "RUGMER", varid)
    call nf95_put_var(ncid_restartphy, varid, frugs(:, is_oce))

    call nf95_inq_varid(ncid_restartphy, "CLWCON", varid)
    call nf95_put_var(ncid_restartphy, varid, clwcon(:, 1))

    call nf95_inq_varid(ncid_restartphy, "RNEBCON", varid)
    call nf95_put_var(ncid_restartphy, varid, rnebcon(:, 1))

    call nf95_inq_varid(ncid_restartphy, "RATQS", varid)
    call nf95_put_var(ncid_restartphy, varid, ratqs(:, 1))

    call nf95_inq_varid(ncid_restartphy, "RUNOFFLIC0", varid)
    call nf95_put_var(ncid_restartphy, varid, run_off_lic_0)

    call nf95_inq_varid(ncid_restartphy, "sig1", varid)
    call nf95_put_var(ncid_restartphy, varid, sig1)

    call nf95_inq_varid(ncid_restartphy, "w01", varid)
    call nf95_put_var(ncid_restartphy, varid, w01)

    call nf95_close(ncid_restartphy)

  END SUBROUTINE phyredem

end module phyredem_m