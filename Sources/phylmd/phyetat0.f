module phyetat0_m

  use dimphy, only: klon

  IMPLICIT none

  REAL, save:: rlat(klon), rlon(klon)
  ! latitude and longitude of a point of the scalar grid identified
  ! by a simple index, in degrees

  private klon

contains

  SUBROUTINE phyetat0(pctsrf, tsol, tsoil, qsurf, qsol, snow, albe, evap, &
       rain_fall, snow_fall, solsw, sollw, fder, radsol, frugs, agesno, zmea, &
       zstd, zsig, zgam, zthe, zpic, zval, t_ancien, q_ancien, ancien_ok, &
       rnebcon, ratqs, clwcon, run_off_lic_0, sig1, w01, ncid_startphy, &
       itau_phy)

    ! From phylmd/phyetat0.F, version 1.4 2005/06/03 10:03:07
    ! Author: Z.X. Li (LMD/CNRS)
    ! Date: 1993/08/18
    ! Objet : lecture de l'état initial pour la physique

    use dimphy, only: zmasq, klev
    USE dimsoil, ONLY : nsoilmx
    USE indicesol, ONLY : epsfra, is_lic, is_oce, is_sic, is_ter, nbsrf
    use netcdf, only: nf90_global, nf90_inq_varid, NF90_NOERR, NF90_NOWRITE
    use netcdf95, only: nf95_get_att, nf95_get_var, nf95_inq_varid, &
         nf95_inquire_variable, NF95_OPEN

    REAL, intent(out):: pctsrf(klon, nbsrf)
    REAL, intent(out):: tsol(klon, nbsrf)
    REAL, intent(out):: tsoil(klon, nsoilmx, nbsrf)
    REAL, intent(out):: qsurf(klon, nbsrf)
    REAL, intent(out):: qsol(:) ! (klon)
    REAL, intent(out):: snow(klon, nbsrf)
    REAL, intent(out):: albe(klon, nbsrf)
    REAL, intent(out):: evap(klon, nbsrf)
    REAL, intent(out):: rain_fall(klon)
    REAL, intent(out):: snow_fall(klon)
    real, intent(out):: solsw(klon)
    REAL, intent(out):: sollw(klon)
    real, intent(out):: fder(klon)
    REAL, intent(out):: radsol(klon)
    REAL, intent(out):: frugs(klon, nbsrf)
    REAL, intent(out):: agesno(klon, nbsrf)
    REAL, intent(out):: zmea(klon)
    REAL, intent(out):: zstd(klon)
    REAL, intent(out):: zsig(klon)
    REAL, intent(out):: zgam(klon)
    REAL, intent(out):: zthe(klon)
    REAL, intent(out):: zpic(klon)
    REAL, intent(out):: zval(klon)
    REAL, intent(out):: t_ancien(klon, klev), q_ancien(klon, klev)
    LOGICAL, intent(out):: ancien_ok
    real, intent(out):: rnebcon(klon, klev), ratqs(klon, klev)
    REAL, intent(out):: clwcon(klon, klev), run_off_lic_0(klon)
    real, intent(out):: sig1(klon, klev) ! section adiabatic updraft

    real, intent(out):: w01(klon, klev) 
    ! vertical velocity within adiabatic updraft

    integer, intent(out):: ncid_startphy, itau_phy

    ! Local:
    REAL fractint(klon)
    INTEGER varid, ndims
    INTEGER ierr, i

    !---------------------------------------------------------------

    print *, "Call sequence information: phyetat0"

    ! Fichier contenant l'état initial :
    call NF95_OPEN("startphy.nc", NF90_NOWRITE, ncid_startphy)

    call nf95_get_att(ncid_startphy, nf90_global, "itau_phy", itau_phy)

    ! Lecture des latitudes (coordonnees):

    call NF95_INQ_VARID(ncid_startphy, "latitude", varid)
    call NF95_GET_VAR(ncid_startphy, varid, rlat)

    ! Lecture des longitudes (coordonnees):

    call NF95_INQ_VARID(ncid_startphy, "longitude", varid)
    call NF95_GET_VAR(ncid_startphy, varid, rlon)

    ! Lecture du masque terre mer

    call NF95_INQ_VARID(ncid_startphy, "masque", varid)
    call nf95_get_var(ncid_startphy, varid, zmasq)

    ! Lecture des fractions pour chaque sous-surface

    ! initialisation des sous-surfaces

    pctsrf = 0.

    ! fraction de terre

    ierr = NF90_INQ_VARID(ncid_startphy, "FTER", varid)
    IF (ierr == NF90_NOERR) THEN
       call nf95_get_var(ncid_startphy, varid, pctsrf(:, is_ter))
    else
       PRINT *, 'phyetat0: Le champ <FTER> est absent'
    ENDIF

    ! fraction de glace de terre

    ierr = NF90_INQ_VARID(ncid_startphy, "FLIC", varid)
    IF (ierr == NF90_NOERR) THEN
       call nf95_get_var(ncid_startphy, varid, pctsrf(:, is_lic))
    else
       PRINT *, 'phyetat0: Le champ <FLIC> est absent'
    ENDIF

    ! fraction d'ocean

    ierr = NF90_INQ_VARID(ncid_startphy, "FOCE", varid)
    IF (ierr == NF90_NOERR) THEN
       call nf95_get_var(ncid_startphy, varid, pctsrf(:, is_oce))
    else
       PRINT *, 'phyetat0: Le champ <FOCE> est absent'
    ENDIF

    ! fraction glace de mer

    ierr = NF90_INQ_VARID(ncid_startphy, "FSIC", varid)
    IF (ierr == NF90_NOERR) THEN
       call nf95_get_var(ncid_startphy, varid, pctsrf(:, is_sic))
    else
       PRINT *, 'phyetat0: Le champ <FSIC> est absent'
    ENDIF

    ! Verification de l'adequation entre le masque et les sous-surfaces

    fractint = pctsrf(:, is_ter) + pctsrf(:, is_lic)
    DO i = 1 , klon
       IF ( abs(fractint(i) - zmasq(i) ) > EPSFRA ) THEN
          WRITE(*, *) 'phyetat0: attention fraction terre pas ', &
               'coherente ', i, zmasq(i), pctsrf(i, is_ter) &
               , pctsrf(i, is_lic)
       ENDIF
    END DO
    fractint = pctsrf(:, is_oce) + pctsrf(:, is_sic)
    DO i = 1 , klon
       IF ( abs( fractint(i) - (1. - zmasq(i))) > EPSFRA ) THEN
          WRITE(*, *) 'phyetat0 attention fraction ocean pas ', &
               'coherente ', i, zmasq(i) , pctsrf(i, is_oce) &
               , pctsrf(i, is_sic)
       ENDIF
    END DO

    ! Lecture des temperatures du sol:
    call NF95_INQ_VARID(ncid_startphy, "TS", varid)
    call nf95_inquire_variable(ncid_startphy, varid, ndims = ndims)
    if (ndims == 2) then
       call NF95_GET_VAR(ncid_startphy, varid, tsol)
    else
       print *, "Found only one surface type for soil temperature."
       call nf95_get_var(ncid_startphy, varid, tsol(:, 1))
       tsol(:, 2:nbsrf) = spread(tsol(:, 1), dim = 2, ncopies = nbsrf - 1)
    end if

    ! Lecture des temperatures du sol profond:

    call NF95_INQ_VARID(ncid_startphy, 'Tsoil', varid)
    call NF95_GET_VAR(ncid_startphy, varid, tsoil)

    ! Lecture de l'humidite de l'air juste au dessus du sol:

    call NF95_INQ_VARID(ncid_startphy, "QS", varid)
    call nf95_get_var(ncid_startphy, varid, qsurf)

    ! Eau dans le sol (pour le modele de sol "bucket")

    ierr = NF90_INQ_VARID(ncid_startphy, "QSOL", varid)
    IF (ierr == NF90_NOERR) THEN
       call nf95_get_var(ncid_startphy, varid, qsol)
    else
       PRINT *, 'phyetat0: Le champ <QSOL> est absent'
       PRINT *, ' Valeur par defaut nulle'
       qsol = 0.
    ENDIF

    ! Lecture de neige au sol:

    call NF95_INQ_VARID(ncid_startphy, "SNOW", varid)
    call nf95_get_var(ncid_startphy, varid, snow)

    ! Lecture de albedo au sol:

    call NF95_INQ_VARID(ncid_startphy, "ALBE", varid)
    call nf95_get_var(ncid_startphy, varid, albe)

    ! Lecture de evaporation: 

    call NF95_INQ_VARID(ncid_startphy, "EVAP", varid)
    call nf95_get_var(ncid_startphy, varid, evap)

    ! Lecture precipitation liquide:

    call NF95_INQ_VARID(ncid_startphy, "rain_f", varid)
    call NF95_GET_VAR(ncid_startphy, varid, rain_fall)

    ! Lecture precipitation solide:

    call NF95_INQ_VARID(ncid_startphy, "snow_f", varid)
    call NF95_GET_VAR(ncid_startphy, varid, snow_fall)

    ! Lecture rayonnement solaire au sol:

    ierr = NF90_INQ_VARID(ncid_startphy, "solsw", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, 'phyetat0: Le champ <solsw> est absent'
       PRINT *, 'mis a zero'
       solsw = 0.
    ELSE
       call nf95_get_var(ncid_startphy, varid, solsw)
    ENDIF

    ! Lecture rayonnement IF au sol:

    ierr = NF90_INQ_VARID(ncid_startphy, "sollw", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, 'phyetat0: Le champ <sollw> est absent'
       PRINT *, 'mis a zero'
       sollw = 0.
    ELSE
       call nf95_get_var(ncid_startphy, varid, sollw)
    ENDIF

    ! Lecture derive des flux:

    ierr = NF90_INQ_VARID(ncid_startphy, "fder", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, 'phyetat0: Le champ <fder> est absent'
       PRINT *, 'mis a zero'
       fder = 0.
    ELSE
       call nf95_get_var(ncid_startphy, varid, fder)
    ENDIF

    ! Lecture du rayonnement net au sol:

    call NF95_INQ_VARID(ncid_startphy, "RADS", varid)
    call NF95_GET_VAR(ncid_startphy, varid, radsol)

    ! Lecture de la longueur de rugosite 

    call NF95_INQ_VARID(ncid_startphy, "RUG", varid)
    call nf95_get_var(ncid_startphy, varid, frugs)

    ! Lecture de l'age de la neige:

    call NF95_INQ_VARID(ncid_startphy, "AGESNO", varid)
    call nf95_get_var(ncid_startphy, varid, agesno)

    call NF95_INQ_VARID(ncid_startphy, "ZMEA", varid)
    call NF95_GET_VAR(ncid_startphy, varid, zmea)

    call NF95_INQ_VARID(ncid_startphy, "ZSTD", varid)
    call NF95_GET_VAR(ncid_startphy, varid, zstd)

    call NF95_INQ_VARID(ncid_startphy, "ZSIG", varid)
    call NF95_GET_VAR(ncid_startphy, varid, zsig)

    call NF95_INQ_VARID(ncid_startphy, "ZGAM", varid)
    call NF95_GET_VAR(ncid_startphy, varid, zgam)

    call NF95_INQ_VARID(ncid_startphy, "ZTHE", varid)
    call NF95_GET_VAR(ncid_startphy, varid, zthe)

    call NF95_INQ_VARID(ncid_startphy, "ZPIC", varid)
    call NF95_GET_VAR(ncid_startphy, varid, zpic)

    call NF95_INQ_VARID(ncid_startphy, "ZVAL", varid)
    call NF95_GET_VAR(ncid_startphy, varid, zval)

    ancien_ok = .TRUE.

    ierr = NF90_INQ_VARID(ncid_startphy, "TANCIEN", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, "phyetat0: Le champ <TANCIEN> est absent"
       PRINT *, "Depart legerement fausse. Mais je continue"
       ancien_ok = .FALSE.
    ELSE
       call nf95_get_var(ncid_startphy, varid, t_ancien)
    ENDIF

    ierr = NF90_INQ_VARID(ncid_startphy, "QANCIEN", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, "phyetat0: Le champ <QANCIEN> est absent"
       PRINT *, "Depart legerement fausse. Mais je continue"
       ancien_ok = .FALSE.
    ELSE
       call nf95_get_var(ncid_startphy, varid, q_ancien)
    ENDIF

    ierr = NF90_INQ_VARID(ncid_startphy, "CLWCON", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, "phyetat0: Le champ CLWCON est absent"
       PRINT *, "Depart legerement fausse. Mais je continue"
       clwcon = 0.
    ELSE
       call nf95_get_var(ncid_startphy, varid, clwcon(:, 1))
       clwcon(:, 2:) = 0.
    ENDIF

    ierr = NF90_INQ_VARID(ncid_startphy, "RNEBCON", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, "phyetat0: Le champ RNEBCON est absent"
       PRINT *, "Depart legerement fausse. Mais je continue"
       rnebcon = 0.
    ELSE
       call nf95_get_var(ncid_startphy, varid, rnebcon(:, 1))
       rnebcon(:, 2:) = 0.
    ENDIF

    ! Lecture ratqs

    ierr = NF90_INQ_VARID(ncid_startphy, "RATQS", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, "phyetat0: Le champ <RATQS> est absent"
       PRINT *, "Depart legerement fausse. Mais je continue"
       ratqs = 0.
    ELSE
       call nf95_get_var(ncid_startphy, varid, ratqs(:, 1))
       ratqs(:, 2:) = 0.
    ENDIF

    ! Lecture run_off_lic_0

    ierr = NF90_INQ_VARID(ncid_startphy, "RUNOFFLIC0", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, "phyetat0: Le champ <RUNOFFLIC0> est absent"
       PRINT *, "Depart legerement fausse. Mais je continue"
       run_off_lic_0 = 0.
    ELSE
       call nf95_get_var(ncid_startphy, varid, run_off_lic_0)
    ENDIF

    call nf95_inq_varid(ncid_startphy, "sig1", varid)
    call nf95_get_var(ncid_startphy, varid, sig1)

    call nf95_inq_varid(ncid_startphy, "w01", varid)
    call nf95_get_var(ncid_startphy, varid, w01)

  END SUBROUTINE phyetat0

end module phyetat0_m
