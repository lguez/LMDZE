module phyetat0_m

  use dimphy, only: klon

  IMPLICIT none

  REAL, save:: rlat(klon), rlon(klon) ! latitude and longitude, in degrees

  private klon

contains

  SUBROUTINE phyetat0(fichnom, pctsrf, tsol, tsoil, ocean, tslab, seaice, &
       qsurf, qsol, snow, albe, alblw, evap, rain_fall, snow_fall, solsw, &
       sollw, fder, radsol, frugs, agesno, zmea, zstd, zsig, zgam, zthe, &
       zpic, zval, t_ancien, q_ancien, ancien_ok, rnebcon, ratqs, clwcon, &
       run_off_lic_0)

    ! From phylmd/phyetat0.F, version 1.4 2005/06/03 10:03:07
    ! Author: Z.X. Li (LMD/CNRS)
    ! Date: 1993/08/18
    ! Objet : Lecture de l'état initial pour la physique

    USE indicesol, ONLY : epsfra, is_lic, is_oce, is_sic, is_ter, nbsrf
    USE dimsoil, ONLY : nsoilmx
    USE temps, ONLY : itau_phy
    use netcdf, only: nf90_get_att, nf90_global, nf90_inq_varid, NF90_NOERR, &
         NF90_NOWRITE
    use netcdf95, only: handle_err, nf95_get_var, nf95_close, NF95_OPEN, &
         nf95_inq_varid
    use dimphy, only: zmasq, klev

    CHARACTER(len=*), intent(in):: fichnom
    REAL tsol(klon, nbsrf)
    REAL tsoil(klon, nsoilmx, nbsrf)
    REAL tslab(klon), seaice(klon)
    REAL qsurf(klon, nbsrf)
    REAL qsol(klon)
    REAL snow(klon, nbsrf)
    REAL albe(klon, nbsrf)
    REAL alblw(klon, nbsrf)
    REAL evap(klon, nbsrf)
    REAL radsol(klon)
    REAL rain_fall(klon)
    REAL snow_fall(klon)
    REAL sollw(klon)
    real solsw(klon)
    real fder(klon)
    REAL frugs(klon, nbsrf)
    REAL agesno(klon, nbsrf)
    REAL zmea(klon)
    REAL, intent(out):: zstd(klon)
    REAL, intent(out):: zsig(klon)
    REAL zgam(klon)
    REAL zthe(klon)
    REAL zpic(klon)
    REAL zval(klon)
    REAL pctsrf(klon, nbsrf)
    REAL fractint(klon)
    REAL run_off_lic_0(klon)

    REAL t_ancien(klon, klev), q_ancien(klon, klev)
    real rnebcon(klon, klev), clwcon(klon, klev), ratqs(klon, klev)
    LOGICAL, intent(out):: ancien_ok

    CHARACTER(len=*), intent(in):: ocean

    REAL xmin, xmax

    INTEGER ncid, varid
    INTEGER ierr, i, nsrf, isoil 
    CHARACTER*7 str7
    CHARACTER*2 str2

    !---------------------------------------------------------------

    print *, "Call sequence information: phyetat0"

    ! Ouvrir le fichier contenant l'etat initial:
    print *, 'fichnom = ', fichnom
    call NF95_OPEN(fichnom, NF90_NOWRITE, ncid)

    ierr = nf90_get_att(ncid, nf90_global, "itau_phy", itau_phy)
    call handle_err("phyetat0 itau_phy", ierr, ncid, nf90_global)

    ! Lecture des latitudes (coordonnees):

    call NF95_INQ_VARID(ncid, "latitude", varid)
    call NF95_GET_VAR(ncid, varid, rlat)

    ! Lecture des longitudes (coordonnees):

    call NF95_INQ_VARID(ncid, "longitude", varid)
    call NF95_GET_VAR(ncid, varid, rlon)

    ! Lecture du masque terre mer

    ierr = NF90_INQ_VARID(ncid, "masque", varid)
    IF (ierr == NF90_NOERR) THEN
       call nf95_get_var(ncid, varid, zmasq)
    else
       PRINT *, 'phyetat0: Le champ <masque> est absent'
       PRINT *, 'fichier startphy non compatible avec phyetat0'
    ENDIF
    ! Lecture des fractions pour chaque sous-surface

    ! initialisation des sous-surfaces

    pctsrf = 0.

    ! fraction de terre

    ierr = NF90_INQ_VARID(ncid, "FTER", varid)
    IF (ierr == NF90_NOERR) THEN
       call nf95_get_var(ncid, varid, pctsrf(:, is_ter))
    else
       PRINT *, 'phyetat0: Le champ <FTER> est absent'
    ENDIF

    ! fraction de glace de terre

    ierr = NF90_INQ_VARID(ncid, "FLIC", varid)
    IF (ierr == NF90_NOERR) THEN
       call nf95_get_var(ncid, varid, pctsrf(:, is_lic))
    else
       PRINT *, 'phyetat0: Le champ <FLIC> est absent'
    ENDIF

    ! fraction d'ocean

    ierr = NF90_INQ_VARID(ncid, "FOCE", varid)
    IF (ierr == NF90_NOERR) THEN
       call nf95_get_var(ncid, varid, pctsrf(:, is_oce))
    else
       PRINT *, 'phyetat0: Le champ <FOCE> est absent'
    ENDIF

    ! fraction glace de mer

    ierr = NF90_INQ_VARID(ncid, "FSIC", varid)
    IF (ierr == NF90_NOERR) THEN
       call nf95_get_var(ncid, varid, pctsrf(:, is_sic))
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

    ierr = NF90_INQ_VARID(ncid, "TS", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, 'phyetat0 : Le champ <TS> est absent'
       PRINT *, ' Mais je vais essayer de lire TS**'
       DO nsrf = 1, nbsrf
          IF (nsrf > 99) THEN
             PRINT *, "Trop de sous-mailles"
             stop 1
          ENDIF
          WRITE(str2, '(i2.2)') nsrf
          call NF95_INQ_VARID(ncid, "TS"//str2, varid)
          call NF95_GET_VAR(ncid, varid, tsol(:, nsrf))
          xmin = 1.0E+20
          xmax = -1.0E+20
          DO i = 1, klon
             xmin = MIN(tsol(i, nsrf), xmin)
             xmax = MAX(tsol(i, nsrf), xmax)
          ENDDO
          PRINT *, 'Temperature du sol TS**:', nsrf, xmin, xmax
       ENDDO
    ELSE
       PRINT *, 'phyetat0: Le champ <TS> est present'
       PRINT *, ' J ignore donc les autres temperatures TS**'
       call nf95_get_var(ncid, varid, tsol(:, 1))
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(tsol(i, 1), xmin)
          xmax = MAX(tsol(i, 1), xmax)
       ENDDO
       PRINT *, 'Temperature du sol <TS>', xmin, xmax
       DO nsrf = 2, nbsrf
          DO i = 1, klon
             tsol(i, nsrf) = tsol(i, 1)
          ENDDO
       ENDDO
    ENDIF

    ! Lecture des temperatures du sol profond:

    DO nsrf = 1, nbsrf
       DO isoil=1, nsoilmx
          IF (isoil > 99 .AND. nsrf > 99) THEN
             PRINT *, "Trop de couches ou sous-mailles"
             stop 1
          ENDIF
          WRITE(str7, '(i2.2, "srf", i2.2)') isoil, nsrf
          ierr = NF90_INQ_VARID(ncid, 'Tsoil'//str7, varid)
          IF (ierr /= NF90_NOERR) THEN
             PRINT *, "phyetat0: Le champ <Tsoil"//str7//"> est absent"
             PRINT *, " Il prend donc la valeur de surface"
             DO i=1, klon
                tsoil(i, isoil, nsrf)=tsol(i, nsrf)
             ENDDO
          ELSE
             call NF95_GET_VAR(ncid, varid, tsoil(:, isoil, nsrf))
          ENDIF
       ENDDO
    ENDDO

    !IM "slab" ocean 

    ! Lecture de tslab (pour slab ocean seulement): 

    IF (ocean .eq. 'slab ') then
       call NF95_INQ_VARID(ncid, "TSLAB", varid)
       call nf95_get_var(ncid, varid, tslab)
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(tslab(i), xmin)
          xmax = MAX(tslab(i), xmax)
       ENDDO
       PRINT *, 'Ecart de la SST tslab:', xmin, xmax

       ! Lecture de seaice (pour slab ocean seulement):

       call NF95_INQ_VARID(ncid, "SEAICE", varid)
       call nf95_get_var(ncid, varid, seaice)
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(seaice(i), xmin)
          xmax = MAX(seaice(i), xmax)
       ENDDO
       PRINT *, 'Masse de la glace de mer seaice:', xmin, xmax
    ELSE
       tslab = 0.
       seaice = 0.
    ENDIF

    ! Lecture de l'humidite de l'air juste au dessus du sol:

    ierr = NF90_INQ_VARID(ncid, "QS", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, 'phyetat0: Le champ <QS> est absent'
       PRINT *, ' Mais je vais essayer de lire QS**'
       DO nsrf = 1, nbsrf
          IF (nsrf > 99) THEN
             PRINT *, "Trop de sous-mailles"
             stop 1
          ENDIF
          WRITE(str2, '(i2.2)') nsrf
          call NF95_INQ_VARID(ncid, "QS"//str2, varid)
          call NF95_GET_VAR(ncid, varid, qsurf(:, nsrf))
          xmin = 1.0E+20
          xmax = -1.0E+20
          DO i = 1, klon
             xmin = MIN(qsurf(i, nsrf), xmin)
             xmax = MAX(qsurf(i, nsrf), xmax)
          ENDDO
          PRINT *, 'Humidite pres du sol QS**:', nsrf, xmin, xmax
       ENDDO
    ELSE
       PRINT *, 'phyetat0: Le champ <QS> est present'
       PRINT *, ' J ignore donc les autres humidites QS**'
       call nf95_get_var(ncid, varid, qsurf(:, 1))
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(qsurf(i, 1), xmin)
          xmax = MAX(qsurf(i, 1), xmax)
       ENDDO
       PRINT *, 'Humidite pres du sol <QS>', xmin, xmax
       DO nsrf = 2, nbsrf
          DO i = 1, klon
             qsurf(i, nsrf) = qsurf(i, 1)
          ENDDO
       ENDDO
    ENDIF

    ! Eau dans le sol (pour le modele de sol "bucket")

    ierr = NF90_INQ_VARID(ncid, "QSOL", varid)
    IF (ierr == NF90_NOERR) THEN
       call nf95_get_var(ncid, varid, qsol)
    else
       PRINT *, 'phyetat0: Le champ <QSOL> est absent'
       PRINT *, ' Valeur par defaut nulle'
       qsol = 0.
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(qsol(i), xmin)
       xmax = MAX(qsol(i), xmax)
    ENDDO
    PRINT *, 'Eau dans le sol (mm) <QSOL>', xmin, xmax

    ! Lecture de neige au sol:

    ierr = NF90_INQ_VARID(ncid, "SNOW", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, 'phyetat0: Le champ <SNOW> est absent'
       PRINT *, ' Mais je vais essayer de lire SNOW**'
       DO nsrf = 1, nbsrf
          IF (nsrf > 99) THEN
             PRINT *, "Trop de sous-mailles"
             stop 1
          ENDIF
          WRITE(str2, '(i2.2)') nsrf
          call NF95_INQ_VARID(ncid, "SNOW"//str2, varid)
          call NF95_GET_VAR(ncid, varid, snow(:, nsrf))
          xmin = 1.0E+20
          xmax = -1.0E+20
          DO i = 1, klon
             xmin = MIN(snow(i, nsrf), xmin)
             xmax = MAX(snow(i, nsrf), xmax)
          ENDDO
          PRINT *, 'Neige du sol SNOW**:', nsrf, xmin, xmax
       ENDDO
    ELSE
       PRINT *, 'phyetat0: Le champ <SNOW> est present'
       PRINT *, ' J ignore donc les autres neiges SNOW**'
       call nf95_get_var(ncid, varid, snow(:, 1))
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(snow(i, 1), xmin)
          xmax = MAX(snow(i, 1), xmax)
       ENDDO
       PRINT *, 'Neige du sol <SNOW>', xmin, xmax
       DO nsrf = 2, nbsrf
          DO i = 1, klon
             snow(i, nsrf) = snow(i, 1)
          ENDDO
       ENDDO
    ENDIF

    ! Lecture de albedo au sol:

    ierr = NF90_INQ_VARID(ncid, "ALBE", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, 'phyetat0: Le champ <ALBE> est absent'
       PRINT *, ' Mais je vais essayer de lire ALBE**'
       DO nsrf = 1, nbsrf
          IF (nsrf > 99) THEN
             PRINT *, "Trop de sous-mailles"
             stop 1
          ENDIF
          WRITE(str2, '(i2.2)') nsrf
          call NF95_INQ_VARID(ncid, "ALBE"//str2, varid)
          call NF95_GET_VAR(ncid, varid, albe(:, nsrf))
          xmin = 1.0E+20
          xmax = -1.0E+20
          DO i = 1, klon
             xmin = MIN(albe(i, nsrf), xmin)
             xmax = MAX(albe(i, nsrf), xmax)
          ENDDO
          PRINT *, 'Albedo du sol ALBE**:', nsrf, xmin, xmax
       ENDDO
    ELSE
       PRINT *, 'phyetat0: Le champ <ALBE> est present'
       PRINT *, ' J ignore donc les autres ALBE**'
       call nf95_get_var(ncid, varid, albe(:, 1))
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(albe(i, 1), xmin)
          xmax = MAX(albe(i, 1), xmax)
       ENDDO
       PRINT *, 'Neige du sol <ALBE>', xmin, xmax
       DO nsrf = 2, nbsrf
          DO i = 1, klon
             albe(i, nsrf) = albe(i, 1)
          ENDDO
       ENDDO
    ENDIF

    ! Lecture de albedo au sol LW:

    ierr = NF90_INQ_VARID(ncid, "ALBLW", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, 'phyetat0: Le champ <ALBLW> est absent'
       ! PRINT *, ' Mais je vais essayer de lire ALBLW**'
       PRINT *, ' Mais je vais prendre ALBE**'
       DO nsrf = 1, nbsrf
          DO i = 1, klon
             alblw(i, nsrf) = albe(i, nsrf)
          ENDDO
       ENDDO
    ELSE
       PRINT *, 'phyetat0: Le champ <ALBLW> est present'
       PRINT *, ' J ignore donc les autres ALBLW**'
       call nf95_get_var(ncid, varid, alblw(:, 1))
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(alblw(i, 1), xmin)
          xmax = MAX(alblw(i, 1), xmax)
       ENDDO
       PRINT *, 'Neige du sol <ALBLW>', xmin, xmax
       DO nsrf = 2, nbsrf
          DO i = 1, klon
             alblw(i, nsrf) = alblw(i, 1)
          ENDDO
       ENDDO
    ENDIF

    ! Lecture de evaporation: 

    ierr = NF90_INQ_VARID(ncid, "EVAP", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, 'phyetat0: Le champ <EVAP> est absent'
       PRINT *, ' Mais je vais essayer de lire EVAP**'
       DO nsrf = 1, nbsrf
          IF (nsrf > 99) THEN
             PRINT *, "Trop de sous-mailles"
             stop 1
          ENDIF
          WRITE(str2, '(i2.2)') nsrf
          call NF95_INQ_VARID(ncid, "EVAP"//str2, varid)
          call NF95_GET_VAR(ncid, varid, evap(:, nsrf))
          xmin = 1.0E+20
          xmax = -1.0E+20
          DO i = 1, klon
             xmin = MIN(evap(i, nsrf), xmin)
             xmax = MAX(evap(i, nsrf), xmax)
          ENDDO
          PRINT *, 'evap du sol EVAP**:', nsrf, xmin, xmax
       ENDDO
    ELSE
       PRINT *, 'phyetat0: Le champ <EVAP> est present'
       PRINT *, ' J ignore donc les autres EVAP**'
       call nf95_get_var(ncid, varid, evap(:, 1))
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(evap(i, 1), xmin)
          xmax = MAX(evap(i, 1), xmax)
       ENDDO
       PRINT *, 'Evap du sol <EVAP>', xmin, xmax
       DO nsrf = 2, nbsrf
          DO i = 1, klon
             evap(i, nsrf) = evap(i, 1)
          ENDDO
       ENDDO
    ENDIF

    ! Lecture precipitation liquide:

    call NF95_INQ_VARID(ncid, "rain_f", varid)
    call NF95_GET_VAR(ncid, varid, rain_fall)
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(rain_fall(i), xmin)
       xmax = MAX(rain_fall(i), xmax)
    ENDDO
    PRINT *, 'Precipitation liquide rain_f:', xmin, xmax

    ! Lecture precipitation solide:

    call NF95_INQ_VARID(ncid, "snow_f", varid)
    call NF95_GET_VAR(ncid, varid, snow_fall)
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(snow_fall(i), xmin)
       xmax = MAX(snow_fall(i), xmax)
    ENDDO
    PRINT *, 'Precipitation solide snow_f:', xmin, xmax

    ! Lecture rayonnement solaire au sol:

    ierr = NF90_INQ_VARID(ncid, "solsw", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, 'phyetat0: Le champ <solsw> est absent'
       PRINT *, 'mis a zero'
       solsw = 0.
    ELSE
       call nf95_get_var(ncid, varid, solsw)
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(solsw(i), xmin)
       xmax = MAX(solsw(i), xmax)
    ENDDO
    PRINT *, 'Rayonnement solaire au sol solsw:', xmin, xmax

    ! Lecture rayonnement IF au sol:

    ierr = NF90_INQ_VARID(ncid, "sollw", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, 'phyetat0: Le champ <sollw> est absent'
       PRINT *, 'mis a zero'
       sollw = 0.
    ELSE
       call nf95_get_var(ncid, varid, sollw)
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(sollw(i), xmin)
       xmax = MAX(sollw(i), xmax)
    ENDDO
    PRINT *, 'Rayonnement IF au sol sollw:', xmin, xmax

    ! Lecture derive des flux:

    ierr = NF90_INQ_VARID(ncid, "fder", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, 'phyetat0: Le champ <fder> est absent'
       PRINT *, 'mis a zero'
       fder = 0.
    ELSE
       call nf95_get_var(ncid, varid, fder)
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(fder(i), xmin)
       xmax = MAX(fder(i), xmax)
    ENDDO
    PRINT *, 'Derive des flux fder:', xmin, xmax

    ! Lecture du rayonnement net au sol:

    call NF95_INQ_VARID(ncid, "RADS", varid)
    call NF95_GET_VAR(ncid, varid, radsol)
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(radsol(i), xmin)
       xmax = MAX(radsol(i), xmax)
    ENDDO
    PRINT *, 'Rayonnement net au sol radsol:', xmin, xmax

    ! Lecture de la longueur de rugosite 

    ierr = NF90_INQ_VARID(ncid, "RUG", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, 'phyetat0: Le champ <RUG> est absent'
       PRINT *, ' Mais je vais essayer de lire RUG**'
       DO nsrf = 1, nbsrf
          IF (nsrf > 99) THEN
             PRINT *, "Trop de sous-mailles"
             stop 1
          ENDIF
          WRITE(str2, '(i2.2)') nsrf
          call NF95_INQ_VARID(ncid, "RUG"//str2, varid)
          call NF95_GET_VAR(ncid, varid, frugs(:, nsrf))
          xmin = 1.0E+20
          xmax = -1.0E+20
          DO i = 1, klon
             xmin = MIN(frugs(i, nsrf), xmin)
             xmax = MAX(frugs(i, nsrf), xmax)
          ENDDO
          PRINT *, 'rugosite du sol RUG**:', nsrf, xmin, xmax
       ENDDO
    ELSE
       PRINT *, 'phyetat0: Le champ <RUG> est present'
       PRINT *, ' J ignore donc les autres RUG**'
       call nf95_get_var(ncid, varid, frugs(:, 1))
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(frugs(i, 1), xmin)
          xmax = MAX(frugs(i, 1), xmax)
       ENDDO
       PRINT *, 'rugosite <RUG>', xmin, xmax
       DO nsrf = 2, nbsrf
          DO i = 1, klon
             frugs(i, nsrf) = frugs(i, 1)
          ENDDO
       ENDDO
    ENDIF

    ! Lecture de l'age de la neige:

    ierr = NF90_INQ_VARID(ncid, "AGESNO", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, 'phyetat0: Le champ <AGESNO> est absent'
       PRINT *, ' Mais je vais essayer de lire AGESNO**'
       DO nsrf = 1, nbsrf
          IF (nsrf > 99) THEN
             PRINT *, "Trop de sous-mailles"
             stop 1
          ENDIF
          WRITE(str2, '(i2.2)') nsrf
          ierr = NF90_INQ_VARID(ncid, "AGESNO"//str2, varid)
          IF (ierr /= NF90_NOERR) THEN
             PRINT *, "phyetat0: Le champ <AGESNO"//str2//"> est absent"
             agesno = 50.0
          ENDIF
          call NF95_GET_VAR(ncid, varid, agesno(:, nsrf))
          xmin = 1.0E+20
          xmax = -1.0E+20
          DO i = 1, klon
             xmin = MIN(agesno(i, nsrf), xmin)
             xmax = MAX(agesno(i, nsrf), xmax)
          ENDDO
          PRINT *, 'Age de la neige AGESNO**:', nsrf, xmin, xmax
       ENDDO
    ELSE
       PRINT *, 'phyetat0: Le champ <AGESNO> est present'
       PRINT *, ' J ignore donc les autres AGESNO**'
       call nf95_get_var(ncid, varid, agesno(:, 1))
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(agesno(i, 1), xmin)
          xmax = MAX(agesno(i, 1), xmax)
       ENDDO
       PRINT *, 'Age de la neige <AGESNO>', xmin, xmax
       DO nsrf = 2, nbsrf
          DO i = 1, klon
             agesno(i, nsrf) = agesno(i, 1)
          ENDDO
       ENDDO
    ENDIF

    call NF95_INQ_VARID(ncid, "ZMEA", varid)
    call NF95_GET_VAR(ncid, varid, zmea)
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(zmea(i), xmin)
       xmax = MAX(zmea(i), xmax)
    ENDDO
    PRINT *, 'OROGRAPHIE SOUS-MAILLE zmea:', xmin, xmax

    call NF95_INQ_VARID(ncid, "ZSTD", varid)
    call NF95_GET_VAR(ncid, varid, zstd)
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(zstd(i), xmin)
       xmax = MAX(zstd(i), xmax)
    ENDDO
    PRINT *, 'OROGRAPHIE SOUS-MAILLE zstd:', xmin, xmax

    call NF95_INQ_VARID(ncid, "ZSIG", varid)
    call NF95_GET_VAR(ncid, varid, zsig)
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(zsig(i), xmin)
       xmax = MAX(zsig(i), xmax)
    ENDDO
    PRINT *, 'OROGRAPHIE SOUS-MAILLE zsig:', xmin, xmax

    call NF95_INQ_VARID(ncid, "ZGAM", varid)
    call NF95_GET_VAR(ncid, varid, zgam)
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(zgam(i), xmin)
       xmax = MAX(zgam(i), xmax)
    ENDDO
    PRINT *, 'OROGRAPHIE SOUS-MAILLE zgam:', xmin, xmax

    call NF95_INQ_VARID(ncid, "ZTHE", varid)
    call NF95_GET_VAR(ncid, varid, zthe)
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(zthe(i), xmin)
       xmax = MAX(zthe(i), xmax)
    ENDDO
    PRINT *, 'OROGRAPHIE SOUS-MAILLE zthe:', xmin, xmax

    call NF95_INQ_VARID(ncid, "ZPIC", varid)
    call NF95_GET_VAR(ncid, varid, zpic)
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(zpic(i), xmin)
       xmax = MAX(zpic(i), xmax)
    ENDDO
    PRINT *, 'OROGRAPHIE SOUS-MAILLE zpic:', xmin, xmax

    call NF95_INQ_VARID(ncid, "ZVAL", varid)
    call NF95_GET_VAR(ncid, varid, zval)
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(zval(i), xmin)
       xmax = MAX(zval(i), xmax)
    ENDDO
    PRINT *, 'OROGRAPHIE SOUS-MAILLE zval:', xmin, xmax

    ancien_ok = .TRUE.

    ierr = NF90_INQ_VARID(ncid, "TANCIEN", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, "phyetat0: Le champ <TANCIEN> est absent"
       PRINT *, "Depart legerement fausse. Mais je continue"
       ancien_ok = .FALSE.
    ELSE
       call nf95_get_var(ncid, varid, t_ancien)
    ENDIF

    ierr = NF90_INQ_VARID(ncid, "QANCIEN", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, "phyetat0: Le champ <QANCIEN> est absent"
       PRINT *, "Depart legerement fausse. Mais je continue"
       ancien_ok = .FALSE.
    ELSE
       call nf95_get_var(ncid, varid, q_ancien)
    ENDIF

    ierr = NF90_INQ_VARID(ncid, "CLWCON", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, "phyetat0: Le champ CLWCON est absent"
       PRINT *, "Depart legerement fausse. Mais je continue"
       clwcon = 0.
    ELSE
       call nf95_get_var(ncid, varid, clwcon)
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    xmin = MINval(clwcon)
    xmax = MAXval(clwcon)
    PRINT *, 'Eau liquide convective (ecart-type) clwcon:', xmin, xmax

    ierr = NF90_INQ_VARID(ncid, "RNEBCON", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, "phyetat0: Le champ RNEBCON est absent"
       PRINT *, "Depart legerement fausse. Mais je continue"
       rnebcon = 0.
    ELSE
       call nf95_get_var(ncid, varid, rnebcon)
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    xmin = MINval(rnebcon)
    xmax = MAXval(rnebcon)
    PRINT *, 'Nebulosite convective (ecart-type) rnebcon:', xmin, xmax

    ! Lecture ratqs

    ierr = NF90_INQ_VARID(ncid, "RATQS", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, "phyetat0: Le champ <RATQS> est absent"
       PRINT *, "Depart legerement fausse. Mais je continue"
       ratqs = 0.
    ELSE
       call nf95_get_var(ncid, varid, ratqs)
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    xmin = MINval(ratqs)
    xmax = MAXval(ratqs)
    PRINT *, '(ecart-type) ratqs:', xmin, xmax

    ! Lecture run_off_lic_0

    ierr = NF90_INQ_VARID(ncid, "RUNOFFLIC0", varid)
    IF (ierr /= NF90_NOERR) THEN
       PRINT *, "phyetat0: Le champ <RUNOFFLIC0> est absent"
       PRINT *, "Depart legerement fausse. Mais je continue"
       run_off_lic_0 = 0.
    ELSE
       call nf95_get_var(ncid, varid, run_off_lic_0)
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    xmin = MINval(run_off_lic_0)
    xmax = MAXval(run_off_lic_0)
    PRINT *, '(ecart-type) run_off_lic_0:', xmin, xmax

    call NF95_CLOSE(ncid)

  END SUBROUTINE phyetat0

end module phyetat0_m
