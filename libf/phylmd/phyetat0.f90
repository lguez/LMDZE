module phyetat0_m

  use dimphy, only: klon, klev, zmasq

  IMPLICIT none

  REAL, save:: rlat(klon), rlon(klon)
  ! latitude et longitude pour chaque point, in degrees

  private klon, klev, zmasq

contains

  SUBROUTINE phyetat0(fichnom, pctsrf, tsol,tsoil, ocean, tslab,seaice, &
       qsurf,qsol,snow, &
       albe, alblw, evap, rain_fall, snow_fall, solsw, sollw, &
       fder,radsol,frugs,agesno,clesphy0, &
       zmea,zstd,zsig,zgam,zthe,zpic,zval,rugsrel, &
       t_ancien,q_ancien,ancien_ok, rnebcon, ratqs,clwcon, &
       run_off_lic_0)

    ! From phylmd/phyetat0.F,v 1.4 2005/06/03 10:03:07
    ! Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
    ! Objet: Lecture de l'etat initial pour la physique

    USE indicesol, ONLY : epsfra, is_lic, is_oce, is_sic, is_ter, nbsrf
    USE dimsoil, ONLY : nsoilmx
    USE temps, ONLY : itau_phy
    USE clesphys2, ONLY : cycle_diurne, iflag_con, nbapp_rad, new_oliq, &
         ok_limitvrai, ok_orodr, ok_orolf, soil_model
    use netcdf, only: nf90_get_att, nf90_global
    use netcdf95, only: handle_err

    include "netcdf.inc"

    CHARACTER(len=*) fichnom
    REAL tsol(klon,nbsrf)
    REAL tsoil(klon,nsoilmx,nbsrf)
    !IM "slab" ocean
    REAL tslab(klon), seaice(klon)
    REAL qsurf(klon,nbsrf)
    REAL qsol(klon)
    REAL snow(klon,nbsrf)
    REAL albe(klon,nbsrf)
    REAL alblw(klon,nbsrf)
    REAL evap(klon,nbsrf)
    REAL radsol(klon)
    REAL rain_fall(klon)
    REAL snow_fall(klon)
    REAL sollw(klon)
    real solsw(klon)
    real fder(klon)
    REAL frugs(klon,nbsrf)
    REAL agesno(klon,nbsrf)
    REAL zmea(klon)
    REAL zstd(klon)
    REAL zsig(klon)
    REAL zgam(klon)
    REAL zthe(klon)
    REAL zpic(klon)
    REAL zval(klon)
    REAL rugsrel(klon)
    REAL pctsrf(klon, nbsrf)
    REAL fractint(klon)
    REAL run_off_lic_0(klon)

    REAL t_ancien(klon,klev), q_ancien(klon,klev)
    real rnebcon(klon,klev),clwcon(klon,klev),ratqs(klon,klev)
    LOGICAL ancien_ok

    CHARACTER(len=*), intent(in):: ocean

    INTEGER        longcles
    PARAMETER    ( longcles = 20 )
    REAL, intent(in):: clesphy0( longcles )

    REAL xmin, xmax

    INTEGER nid, nvarid
    INTEGER ierr, i, nsrf, isoil 
    INTEGER length
    PARAMETER (length=100)
    CHARACTER*7 str7
    CHARACTER*2 str2

    !---------------------------------------------------------------

    print *, "Call sequence information: phyetat0"

    ! Ouvrir le fichier contenant l'etat initial:

    print *, 'fichnom = ', fichnom
    ierr = NF_OPEN (fichnom, NF_NOWRITE,nid)
    IF (ierr.NE.NF_NOERR) THEN
       write(6,*)' Pb d''ouverture du fichier '//fichnom
       write(6,*)' ierr = ', ierr
       STOP 1
    ENDIF

    iflag_con    = clesphy0(1)
    nbapp_rad    = clesphy0(2)
    cycle_diurne  = clesphy0(3) == 1
    soil_model = clesphy0(4) == 1
    new_oliq = clesphy0(5) == 1
    ok_orodr = clesphy0(6) == 1
    ok_orolf = clesphy0(7) == 1
    ok_limitvrai = clesphy0(8) == 1

    ierr = nf90_get_att(nid, nf90_global, "itau_phy", itau_phy)
    call handle_err("phyetat0 itau_phy", ierr, nid, nf90_global)

    ! Lecture des latitudes (coordonnees):

    ierr = NF_INQ_VARID (nid, "latitude", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <latitude> est absent'
       stop 1
    ENDIF
    ierr = NF_GET_VAR_REAL(nid, nvarid, rlat)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Lecture echouee pour <latitude>'
       stop 1
    ENDIF

    ! Lecture des longitudes (coordonnees):

    ierr = NF_INQ_VARID (nid, "longitude", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <longitude> est absent'
       stop 1
    ENDIF
    ierr = NF_GET_VAR_REAL(nid, nvarid, rlon)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Lecture echouee pour <latitude>'
       stop 1
    ENDIF


    ! Lecture du masque terre mer

    ierr = NF_INQ_VARID (nid, "masque", nvarid)
    IF (ierr .EQ.  NF_NOERR) THEN
       ierr = NF_GET_VAR_REAL(nid, nvarid, zmasq)
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, 'phyetat0: Lecture echouee pour <masque>'
          stop 1
       ENDIF
    else
       PRINT*, 'phyetat0: Le champ <masque> est absent'
       PRINT*, 'fichier startphy non compatible avec phyetat0'
       !      stop 1
    ENDIF
    ! Lecture des fractions pour chaque sous-surface

    ! initialisation des sous-surfaces

    pctsrf = 0.

    ! fraction de terre

    ierr = NF_INQ_VARID (nid, "FTER", nvarid)
    IF (ierr .EQ.  NF_NOERR) THEN
       ierr = NF_GET_VAR_REAL(nid, nvarid, pctsrf(1 : klon,is_ter))
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, 'phyetat0: Lecture echouee pour <FTER>'
          stop 1
       ENDIF
    else
       PRINT*, 'phyetat0: Le champ <FTER> est absent'
       !$$$         stop 1
    ENDIF

    ! fraction de glace de terre

    ierr = NF_INQ_VARID (nid, "FLIC", nvarid)
    IF (ierr .EQ.  NF_NOERR) THEN
       ierr = NF_GET_VAR_REAL(nid, nvarid, pctsrf(1 : klon,is_lic))
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, 'phyetat0: Lecture echouee pour <FLIC>'
          stop 1
       ENDIF
    else
       PRINT*, 'phyetat0: Le champ <FLIC> est absent'
       !$$$         stop 1
    ENDIF

    ! fraction d'ocean

    ierr = NF_INQ_VARID (nid, "FOCE", nvarid)
    IF (ierr .EQ.  NF_NOERR) THEN
       ierr = NF_GET_VAR_REAL(nid, nvarid, pctsrf(1 : klon,is_oce))
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, 'phyetat0: Lecture echouee pour <FOCE>'
          stop 1
       ENDIF
    else
       PRINT*, 'phyetat0: Le champ <FOCE> est absent'
       !$$$         stop 1
    ENDIF

    ! fraction glace de mer

    ierr = NF_INQ_VARID (nid, "FSIC", nvarid)
    IF (ierr .EQ.  NF_NOERR) THEN
       ierr = NF_GET_VAR_REAL(nid, nvarid, pctsrf(1 : klon, is_sic))
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, 'phyetat0: Lecture echouee pour <FSIC>'
          stop 1
       ENDIF
    else
       PRINT*, 'phyetat0: Le champ <FSIC> est absent'
       !$$$         stop 1
    ENDIF

    !  Verification de l'adequation entre le masque et les sous-surfaces

    fractint( 1 : klon) = pctsrf(1 : klon, is_ter)  &
         + pctsrf(1 : klon, is_lic)
    DO i = 1 , klon
       IF ( abs(fractint(i) - zmasq(i) ) .GT. EPSFRA ) THEN
          WRITE(*,*) 'phyetat0: attention fraction terre pas ',  &
               'coherente ', i, zmasq(i), pctsrf(i, is_ter) &
               ,pctsrf(i, is_lic)
       ENDIF
    END DO
    fractint (1 : klon) =  pctsrf(1 : klon, is_oce)  &
         + pctsrf(1 : klon, is_sic)
    DO i = 1 , klon
       IF ( abs( fractint(i) - (1. - zmasq(i))) .GT. EPSFRA ) THEN
          WRITE(*,*) 'phyetat0 attention fraction ocean pas ',  &
               'coherente ', i, zmasq(i) , pctsrf(i, is_oce) &
               ,pctsrf(i, is_sic)
       ENDIF
    END DO

    ! Lecture des temperatures du sol:

    ierr = NF_INQ_VARID (nid, "TS", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <TS> est absent'
       PRINT*, '          Mais je vais essayer de lire TS**'
       DO nsrf = 1, nbsrf
          IF (nsrf.GT.99) THEN
             PRINT*, "Trop de sous-mailles"
             stop 1
          ENDIF
          WRITE(str2,'(i2.2)') nsrf
          ierr = NF_INQ_VARID (nid, "TS"//str2, nvarid)
          IF (ierr.NE.NF_NOERR) THEN
             PRINT*, "phyetat0: Le champ <TS"//str2//"> est absent"
             stop 1
          ENDIF
          ierr = NF_GET_VAR_REAL(nid, nvarid, tsol(1,nsrf))
          IF (ierr.NE.NF_NOERR) THEN
             PRINT*, "phyetat0: Lecture echouee pour <TS"//str2//">"
             stop 1
          ENDIF
          xmin = 1.0E+20
          xmax = -1.0E+20
          DO i = 1, klon
             xmin = MIN(tsol(i,nsrf),xmin)
             xmax = MAX(tsol(i,nsrf),xmax)
          ENDDO
          PRINT*,'Temperature du sol TS**:', nsrf, xmin, xmax
       ENDDO
    ELSE
       PRINT*, 'phyetat0: Le champ <TS> est present'
       PRINT*, '          J ignore donc les autres temperatures TS**'
       ierr = NF_GET_VAR_REAL(nid, nvarid, tsol(1,1))
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Lecture echouee pour <TS>"
          stop 1
       ENDIF
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(tsol(i,1),xmin)
          xmax = MAX(tsol(i,1),xmax)
       ENDDO
       PRINT*,'Temperature du sol <TS>', xmin, xmax
       DO nsrf = 2, nbsrf
          DO i = 1, klon
             tsol(i,nsrf) = tsol(i,1)
          ENDDO
       ENDDO
    ENDIF

    ! Lecture des temperatures du sol profond:

    DO nsrf = 1, nbsrf
       DO isoil=1, nsoilmx
          IF (isoil.GT.99 .AND. nsrf.GT.99) THEN
             PRINT*, "Trop de couches ou sous-mailles"
             stop 1
          ENDIF
          WRITE(str7,'(i2.2,"srf",i2.2)') isoil, nsrf
          ierr = NF_INQ_VARID (nid, 'Tsoil'//str7, nvarid)
          IF (ierr.NE.NF_NOERR) THEN
             PRINT*, "phyetat0: Le champ <Tsoil"//str7//"> est absent"
             PRINT*, "          Il prend donc la valeur de surface"
             DO i=1, klon
                tsoil(i,isoil,nsrf)=tsol(i,nsrf)
             ENDDO
          ELSE
             ierr = NF_GET_VAR_REAL(nid, nvarid, tsoil(1,isoil,nsrf))
             IF (ierr.NE.NF_NOERR) THEN
                PRINT*, "Lecture echouee pour <Tsoil"//str7//">"
                stop 1
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    !IM "slab" ocean 

    ! Lecture de tslab (pour slab ocean seulement):      

    IF (ocean .eq. 'slab  ') then
       ierr = NF_INQ_VARID (nid, "TSLAB", nvarid)
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Le champ <TSLAB> est absent"
          stop 1
       ENDIF
       ierr = NF_GET_VAR_REAL(nid, nvarid, tslab)
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Lecture echouee pour <TSLAB>"
          stop 1
       ENDIF
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(tslab(i),xmin)
          xmax = MAX(tslab(i),xmax)
       ENDDO
       PRINT*,'Ecart de la SST tslab:', xmin, xmax

       ! Lecture de seaice (pour slab ocean seulement):

       ierr = NF_INQ_VARID (nid, "SEAICE", nvarid)
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Le champ <SEAICE> est absent"
          stop 1
       ENDIF
       ierr = NF_GET_VAR_REAL(nid, nvarid, seaice)
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Lecture echouee pour <SEAICE>"
          stop 1
       ENDIF
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(seaice(i),xmin)
          xmax = MAX(seaice(i),xmax)
       ENDDO
       PRINT*,'Masse de la glace de mer seaice:', xmin, xmax
    ELSE
       tslab = 0.
       seaice = 0.
    ENDIF

    ! Lecture de l'humidite de l'air juste au dessus du sol:

    ierr = NF_INQ_VARID (nid, "QS", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <QS> est absent'
       PRINT*, '          Mais je vais essayer de lire QS**'
       DO nsrf = 1, nbsrf
          IF (nsrf.GT.99) THEN
             PRINT*, "Trop de sous-mailles"
             stop 1
          ENDIF
          WRITE(str2,'(i2.2)') nsrf
          ierr = NF_INQ_VARID (nid, "QS"//str2, nvarid)
          IF (ierr.NE.NF_NOERR) THEN
             PRINT*, "phyetat0: Le champ <QS"//str2//"> est absent"
             stop 1
          ENDIF
          ierr = NF_GET_VAR_REAL(nid, nvarid, qsurf(1,nsrf))
          IF (ierr.NE.NF_NOERR) THEN
             PRINT*, "phyetat0: Lecture echouee pour <QS"//str2//">"
             stop 1
          ENDIF
          xmin = 1.0E+20
          xmax = -1.0E+20
          DO i = 1, klon
             xmin = MIN(qsurf(i,nsrf),xmin)
             xmax = MAX(qsurf(i,nsrf),xmax)
          ENDDO
          PRINT*,'Humidite pres du sol QS**:', nsrf, xmin, xmax
       ENDDO
    ELSE
       PRINT*, 'phyetat0: Le champ <QS> est present'
       PRINT*, '          J ignore donc les autres humidites QS**'
       ierr = NF_GET_VAR_REAL(nid, nvarid, qsurf(1,1))
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Lecture echouee pour <QS>"
          stop 1
       ENDIF
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(qsurf(i,1),xmin)
          xmax = MAX(qsurf(i,1),xmax)
       ENDDO
       PRINT*,'Humidite pres du sol <QS>', xmin, xmax
       DO nsrf = 2, nbsrf
          DO i = 1, klon
             qsurf(i,nsrf) = qsurf(i,1)
          ENDDO
       ENDDO
    ENDIF

    ! Eau dans le sol (pour le modele de sol "bucket")

    ierr = NF_INQ_VARID (nid, "QSOL", nvarid)
    IF (ierr .EQ.  NF_NOERR) THEN
       ierr = NF_GET_VAR_REAL(nid, nvarid, qsol)
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, 'phyetat0: Lecture echouee pour <QSOL>'
          stop 1
       ENDIF
    else
       PRINT*, 'phyetat0: Le champ <QSOL> est absent'
       PRINT*, '          Valeur par defaut nulle'
       qsol(:)=0.
       !$$$         stop 1
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(qsol(i),xmin)
       xmax = MAX(qsol(i),xmax)
    ENDDO
    PRINT*,'Eau dans le sol (mm) <QSOL>', xmin, xmax

    ! Lecture de neige au sol:

    ierr = NF_INQ_VARID (nid, "SNOW", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <SNOW> est absent'
       PRINT*, '          Mais je vais essayer de lire SNOW**'
       DO nsrf = 1, nbsrf
          IF (nsrf.GT.99) THEN
             PRINT*, "Trop de sous-mailles"
             stop 1
          ENDIF
          WRITE(str2,'(i2.2)') nsrf
          ierr = NF_INQ_VARID (nid, "SNOW"//str2, nvarid)
          IF (ierr.NE.NF_NOERR) THEN
             PRINT*, "phyetat0: Le champ <SNOW"//str2//"> est absent"
             stop 1
          ENDIF
          ierr = NF_GET_VAR_REAL(nid, nvarid, snow(1,nsrf))
          IF (ierr.NE.NF_NOERR) THEN
             PRINT*, "phyetat0: Lecture echouee pour <SNOW"//str2//">"
             stop 1
          ENDIF
          xmin = 1.0E+20
          xmax = -1.0E+20
          DO i = 1, klon
             xmin = MIN(snow(i,nsrf),xmin)
             xmax = MAX(snow(i,nsrf),xmax)
          ENDDO
          PRINT*,'Neige du sol SNOW**:', nsrf, xmin, xmax
       ENDDO
    ELSE
       PRINT*, 'phyetat0: Le champ <SNOW> est present'
       PRINT*, '          J ignore donc les autres neiges SNOW**'
       ierr = NF_GET_VAR_REAL(nid, nvarid, snow(1,1))
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Lecture echouee pour <SNOW>"
          stop 1
       ENDIF
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(snow(i,1),xmin)
          xmax = MAX(snow(i,1),xmax)
       ENDDO
       PRINT*,'Neige du sol <SNOW>', xmin, xmax
       DO nsrf = 2, nbsrf
          DO i = 1, klon
             snow(i,nsrf) = snow(i,1)
          ENDDO
       ENDDO
    ENDIF

    ! Lecture de albedo au sol:

    ierr = NF_INQ_VARID (nid, "ALBE", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <ALBE> est absent'
       PRINT*, '          Mais je vais essayer de lire ALBE**'
       DO nsrf = 1, nbsrf
          IF (nsrf.GT.99) THEN
             PRINT*, "Trop de sous-mailles"
             stop 1
          ENDIF
          WRITE(str2,'(i2.2)') nsrf
          ierr = NF_INQ_VARID (nid, "ALBE"//str2, nvarid)
          IF (ierr.NE.NF_NOERR) THEN
             PRINT*, "phyetat0: Le champ <ALBE"//str2//"> est absent"
             stop 1
          ENDIF
          ierr = NF_GET_VAR_REAL(nid, nvarid, albe(1,nsrf))
          IF (ierr.NE.NF_NOERR) THEN
             PRINT*, "phyetat0: Lecture echouee pour <ALBE"//str2//">"
             stop 1
          ENDIF
          xmin = 1.0E+20
          xmax = -1.0E+20
          DO i = 1, klon
             xmin = MIN(albe(i,nsrf),xmin)
             xmax = MAX(albe(i,nsrf),xmax)
          ENDDO
          PRINT*,'Albedo du sol ALBE**:', nsrf, xmin, xmax
       ENDDO
    ELSE
       PRINT*, 'phyetat0: Le champ <ALBE> est present'
       PRINT*, '          J ignore donc les autres ALBE**'
       ierr = NF_GET_VAR_REAL(nid, nvarid, albe(1,1))
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Lecture echouee pour <ALBE>"
          stop 1
       ENDIF
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(albe(i,1),xmin)
          xmax = MAX(albe(i,1),xmax)
       ENDDO
       PRINT*,'Neige du sol <ALBE>', xmin, xmax
       DO nsrf = 2, nbsrf
          DO i = 1, klon
             albe(i,nsrf) = albe(i,1)
          ENDDO
       ENDDO
    ENDIF


    ! Lecture de albedo au sol LW:

    ierr = NF_INQ_VARID (nid, "ALBLW", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <ALBLW> est absent'
       !        PRINT*, '          Mais je vais essayer de lire ALBLW**'
       PRINT*, '          Mais je vais prendre ALBE**'
       DO nsrf = 1, nbsrf
          DO i = 1, klon
             alblw(i,nsrf) = albe(i,nsrf)
          ENDDO
       ENDDO
    ELSE
       PRINT*, 'phyetat0: Le champ <ALBLW> est present'
       PRINT*, '          J ignore donc les autres ALBLW**'
       ierr = NF_GET_VAR_REAL(nid, nvarid, alblw(1,1))
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Lecture echouee pour <ALBLW>"
          stop 1
       ENDIF
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(alblw(i,1),xmin)
          xmax = MAX(alblw(i,1),xmax)
       ENDDO
       PRINT*,'Neige du sol <ALBLW>', xmin, xmax
       DO nsrf = 2, nbsrf
          DO i = 1, klon
             alblw(i,nsrf) = alblw(i,1)
          ENDDO
       ENDDO
    ENDIF

    ! Lecture de evaporation:  

    ierr = NF_INQ_VARID (nid, "EVAP", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <EVAP> est absent'
       PRINT*, '          Mais je vais essayer de lire EVAP**'
       DO nsrf = 1, nbsrf
          IF (nsrf.GT.99) THEN
             PRINT*, "Trop de sous-mailles"
             stop 1
          ENDIF
          WRITE(str2,'(i2.2)') nsrf
          ierr = NF_INQ_VARID (nid, "EVAP"//str2, nvarid)
          IF (ierr.NE.NF_NOERR) THEN
             PRINT*, "phyetat0: Le champ <EVAP"//str2//"> est absent"
             stop 1
          ENDIF
          ierr = NF_GET_VAR_REAL(nid, nvarid, evap(1,nsrf))
          IF (ierr.NE.NF_NOERR) THEN
             PRINT*, "phyetat0: Lecture echouee pour <EVAP"//str2//">"
             stop 1
          ENDIF
          xmin = 1.0E+20
          xmax = -1.0E+20
          DO i = 1, klon
             xmin = MIN(evap(i,nsrf),xmin)
             xmax = MAX(evap(i,nsrf),xmax)
          ENDDO
          PRINT*,'evap du sol EVAP**:', nsrf, xmin, xmax
       ENDDO
    ELSE
       PRINT*, 'phyetat0: Le champ <EVAP> est present'
       PRINT*, '          J ignore donc les autres EVAP**'
       ierr = NF_GET_VAR_REAL(nid, nvarid, evap(1,1))
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Lecture echouee pour <EVAP>"
          stop 1
       ENDIF
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(evap(i,1),xmin)
          xmax = MAX(evap(i,1),xmax)
       ENDDO
       PRINT*,'Evap du sol <EVAP>', xmin, xmax
       DO nsrf = 2, nbsrf
          DO i = 1, klon
             evap(i,nsrf) = evap(i,1)
          ENDDO
       ENDDO
    ENDIF

    ! Lecture precipitation liquide:

    ierr = NF_INQ_VARID (nid, "rain_f", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <rain_f> est absent'
       stop 1
    ENDIF
    ierr = NF_GET_VAR_REAL(nid, nvarid, rain_fall)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Lecture echouee pour <rain_f>'
       stop 1
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(rain_fall(i),xmin)
       xmax = MAX(rain_fall(i),xmax)
    ENDDO
    PRINT*,'Precipitation liquide rain_f:', xmin, xmax

    ! Lecture precipitation solide:

    ierr = NF_INQ_VARID (nid, "snow_f", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <snow_f> est absent'
       stop 1
    ENDIF
    ierr = NF_GET_VAR_REAL(nid, nvarid, snow_fall)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Lecture echouee pour <snow_f>'
       stop 1
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(snow_fall(i),xmin)
       xmax = MAX(snow_fall(i),xmax)
    ENDDO
    PRINT*,'Precipitation solide snow_f:', xmin, xmax

    ! Lecture rayonnement solaire au sol:

    ierr = NF_INQ_VARID (nid, "solsw", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <solsw> est absent'
       PRINT*, 'mis a zero'
       solsw = 0.
    ELSE
       ierr = NF_GET_VAR_REAL(nid, nvarid, solsw)
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, 'phyetat0: Lecture echouee pour <solsw>'
          stop 1
       ENDIF
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(solsw(i),xmin)
       xmax = MAX(solsw(i),xmax)
    ENDDO
    PRINT*,'Rayonnement solaire au sol solsw:', xmin, xmax

    ! Lecture rayonnement IF au sol:

    ierr = NF_INQ_VARID (nid, "sollw", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <sollw> est absent'
       PRINT*, 'mis a zero'
       sollw = 0.
    ELSE
       ierr = NF_GET_VAR_REAL(nid, nvarid, sollw)
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, 'phyetat0: Lecture echouee pour <sollw>'
          stop 1
       ENDIF
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(sollw(i),xmin)
       xmax = MAX(sollw(i),xmax)
    ENDDO
    PRINT*,'Rayonnement IF au sol sollw:', xmin, xmax


    ! Lecture derive des flux:

    ierr = NF_INQ_VARID (nid, "fder", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <fder> est absent'
       PRINT*, 'mis a zero'
       fder = 0.
    ELSE
       ierr = NF_GET_VAR_REAL(nid, nvarid, fder)
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, 'phyetat0: Lecture echouee pour <fder>'
          stop 1
       ENDIF
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(fder(i),xmin)
       xmax = MAX(fder(i),xmax)
    ENDDO
    PRINT*,'Derive des flux fder:', xmin, xmax


    ! Lecture du rayonnement net au sol:

    ierr = NF_INQ_VARID (nid, "RADS", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <RADS> est absent'
       stop 1
    ENDIF
    ierr = NF_GET_VAR_REAL(nid, nvarid, radsol)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Lecture echouee pour <RADS>'
       stop 1
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(radsol(i),xmin)
       xmax = MAX(radsol(i),xmax)
    ENDDO
    PRINT*,'Rayonnement net au sol radsol:', xmin, xmax

    ! Lecture de la longueur de rugosite 


    ierr = NF_INQ_VARID (nid, "RUG", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <RUG> est absent'
       PRINT*, '          Mais je vais essayer de lire RUG**'
       DO nsrf = 1, nbsrf
          IF (nsrf.GT.99) THEN
             PRINT*, "Trop de sous-mailles"
             stop 1
          ENDIF
          WRITE(str2,'(i2.2)') nsrf
          ierr = NF_INQ_VARID (nid, "RUG"//str2, nvarid)
          IF (ierr.NE.NF_NOERR) THEN
             PRINT*, "phyetat0: Le champ <RUG"//str2//"> est absent"
             stop 1
          ENDIF
          ierr = NF_GET_VAR_REAL(nid, nvarid, frugs(1,nsrf))
          IF (ierr.NE.NF_NOERR) THEN
             PRINT*, "phyetat0: Lecture echouee pour <RUG"//str2//">"
             stop 1
          ENDIF
          xmin = 1.0E+20
          xmax = -1.0E+20
          DO i = 1, klon
             xmin = MIN(frugs(i,nsrf),xmin)
             xmax = MAX(frugs(i,nsrf),xmax)
          ENDDO
          PRINT*,'rugosite du sol RUG**:', nsrf, xmin, xmax
       ENDDO
    ELSE
       PRINT*, 'phyetat0: Le champ <RUG> est present'
       PRINT*, '          J ignore donc les autres RUG**'
       ierr = NF_GET_VAR_REAL(nid, nvarid, frugs(1,1))
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Lecture echouee pour <RUG>"
          stop 1
       ENDIF
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(frugs(i,1),xmin)
          xmax = MAX(frugs(i,1),xmax)
       ENDDO
       PRINT*,'rugosite <RUG>', xmin, xmax
       DO nsrf = 2, nbsrf
          DO i = 1, klon
             frugs(i,nsrf) = frugs(i,1)
          ENDDO
       ENDDO
    ENDIF


    ! Lecture de l'age de la neige:

    ierr = NF_INQ_VARID (nid, "AGESNO", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <AGESNO> est absent'
       PRINT*, '          Mais je vais essayer de lire AGESNO**'
       DO nsrf = 1, nbsrf
          IF (nsrf.GT.99) THEN
             PRINT*, "Trop de sous-mailles"
             stop 1
          ENDIF
          WRITE(str2,'(i2.2)') nsrf
          ierr = NF_INQ_VARID (nid, "AGESNO"//str2, nvarid)
          IF (ierr.NE.NF_NOERR) THEN
             PRINT*, "phyetat0: Le champ <AGESNO"//str2//"> est absent"
             agesno = 50.0
          ENDIF
          ierr = NF_GET_VAR_REAL(nid, nvarid, agesno(1,nsrf))
          IF (ierr.NE.NF_NOERR) THEN
             PRINT*, "phyetat0: Lecture echouee pour <AGESNO"//str2//">"
             stop 1
          ENDIF
          xmin = 1.0E+20
          xmax = -1.0E+20
          DO i = 1, klon
             xmin = MIN(agesno(i,nsrf),xmin)
             xmax = MAX(agesno(i,nsrf),xmax)
          ENDDO
          PRINT*,'Age de la neige AGESNO**:', nsrf, xmin, xmax
       ENDDO
    ELSE
       PRINT*, 'phyetat0: Le champ <AGESNO> est present'
       PRINT*, '          J ignore donc les autres AGESNO**'
       ierr = NF_GET_VAR_REAL(nid, nvarid, agesno(1,1))
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Lecture echouee pour <AGESNO>"
          stop 1
       ENDIF
       xmin = 1.0E+20
       xmax = -1.0E+20
       DO i = 1, klon
          xmin = MIN(agesno(i,1),xmin)
          xmax = MAX(agesno(i,1),xmax)
       ENDDO
       PRINT*,'Age de la neige <AGESNO>', xmin, xmax
       DO nsrf = 2, nbsrf
          DO i = 1, klon
             agesno(i,nsrf) = agesno(i,1)
          ENDDO
       ENDDO
    ENDIF


    ierr = NF_INQ_VARID (nid, "ZMEA", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <ZMEA> est absent'
       stop 1
    ENDIF
    ierr = NF_GET_VAR_REAL(nid, nvarid, zmea)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Lecture echouee pour <ZMEA>'
       stop 1
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(zmea(i),xmin)
       xmax = MAX(zmea(i),xmax)
    ENDDO
    PRINT*,'OROGRAPHIE SOUS-MAILLE zmea:', xmin, xmax


    ierr = NF_INQ_VARID (nid, "ZSTD", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <ZSTD> est absent'
       stop 1
    ENDIF
    ierr = NF_GET_VAR_REAL(nid, nvarid, zstd)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Lecture echouee pour <ZSTD>'
       stop 1
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(zstd(i),xmin)
       xmax = MAX(zstd(i),xmax)
    ENDDO
    PRINT*,'OROGRAPHIE SOUS-MAILLE zstd:', xmin, xmax


    ierr = NF_INQ_VARID (nid, "ZSIG", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <ZSIG> est absent'
       stop 1
    ENDIF
    ierr = NF_GET_VAR_REAL(nid, nvarid, zsig)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Lecture echouee pour <ZSIG>'
       stop 1
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(zsig(i),xmin)
       xmax = MAX(zsig(i),xmax)
    ENDDO
    PRINT*,'OROGRAPHIE SOUS-MAILLE zsig:', xmin, xmax


    ierr = NF_INQ_VARID (nid, "ZGAM", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <ZGAM> est absent'
       stop 1
    ENDIF
    ierr = NF_GET_VAR_REAL(nid, nvarid, zgam)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Lecture echouee pour <ZGAM>'
       stop 1
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(zgam(i),xmin)
       xmax = MAX(zgam(i),xmax)
    ENDDO
    PRINT*,'OROGRAPHIE SOUS-MAILLE zgam:', xmin, xmax


    ierr = NF_INQ_VARID (nid, "ZTHE", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <ZTHE> est absent'
       stop 1
    ENDIF
    ierr = NF_GET_VAR_REAL(nid, nvarid, zthe)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Lecture echouee pour <ZTHE>'
       stop 1
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(zthe(i),xmin)
       xmax = MAX(zthe(i),xmax)
    ENDDO
    PRINT*,'OROGRAPHIE SOUS-MAILLE zthe:', xmin, xmax


    ierr = NF_INQ_VARID (nid, "ZPIC", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <ZPIC> est absent'
       stop 1
    ENDIF
    ierr = NF_GET_VAR_REAL(nid, nvarid, zpic)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Lecture echouee pour <ZPIC>'
       stop 1
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(zpic(i),xmin)
       xmax = MAX(zpic(i),xmax)
    ENDDO
    PRINT*,'OROGRAPHIE SOUS-MAILLE zpic:', xmin, xmax

    ierr = NF_INQ_VARID (nid, "ZVAL", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <ZVAL> est absent'
       stop 1
    ENDIF
    ierr = NF_GET_VAR_REAL(nid, nvarid, zval)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Lecture echouee pour <ZVAL>'
       stop 1
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(zval(i),xmin)
       xmax = MAX(zval(i),xmax)
    ENDDO
    PRINT*,'OROGRAPHIE SOUS-MAILLE zval:', xmin, xmax


    ierr = NF_INQ_VARID (nid, "RUGSREL", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Le champ <RUGSREL> est absent'
       stop 1
    ENDIF
    ierr = NF_GET_VAR_REAL(nid, nvarid, rugsrel)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, 'phyetat0: Lecture echouee pour <RUGSREL>'
       stop 1
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    DO i = 1, klon
       xmin = MIN(rugsrel(i),xmin)
       xmax = MAX(rugsrel(i),xmax)
    ENDDO
    PRINT*,'Rugosite relief (ecart-type) rugsrel:', xmin, xmax


    ancien_ok = .TRUE.

    ierr = NF_INQ_VARID (nid, "TANCIEN", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, "phyetat0: Le champ <TANCIEN> est absent"
       PRINT*, "Depart legerement fausse. Mais je continue"
       ancien_ok = .FALSE.
    ELSE
       ierr = NF_GET_VAR_REAL(nid, nvarid, t_ancien)
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Lecture echouee pour <TANCIEN>"
          stop 1
       ENDIF
    ENDIF

    ierr = NF_INQ_VARID (nid, "QANCIEN", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, "phyetat0: Le champ <QANCIEN> est absent"
       PRINT*, "Depart legerement fausse. Mais je continue"
       ancien_ok = .FALSE.
    ELSE
       ierr = NF_GET_VAR_REAL(nid, nvarid, q_ancien)
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Lecture echouee pour <QANCIEN>"
          stop 1
       ENDIF
    ENDIF

    ierr = NF_INQ_VARID (nid, "CLWCON", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, "phyetat0: Le champ CLWCON est absent"
       PRINT*, "Depart legerement fausse. Mais je continue"
       clwcon = 0.
    ELSE
       ierr = NF_GET_VAR_REAL(nid, nvarid, clwcon)
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Lecture echouee pour <CLWCON>"
          stop 1
       ENDIF
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    xmin = MINval(clwcon)
    xmax = MAXval(clwcon)
    PRINT*,'Eau liquide convective (ecart-type) clwcon:', xmin, xmax

    ierr = NF_INQ_VARID (nid, "RNEBCON", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, "phyetat0: Le champ RNEBCON est absent"
       PRINT*, "Depart legerement fausse. Mais je continue"
       rnebcon = 0.
    ELSE
       ierr = NF_GET_VAR_REAL(nid, nvarid, rnebcon)
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Lecture echouee pour <RNEBCON>"
          stop 1
       ENDIF
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    xmin = MINval(rnebcon)
    xmax = MAXval(rnebcon)
    PRINT*,'Nebulosite convective (ecart-type) rnebcon:', xmin, xmax


    ierr = NF_INQ_VARID (nid, "QANCIEN", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, "phyetat0: Le champ <QANCIEN> est absent"
       PRINT*, "Depart legerement fausse. Mais je continue"
       ancien_ok = .FALSE.
    ELSE
       ierr = NF_GET_VAR_REAL(nid, nvarid, q_ancien)
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Lecture echouee pour <QANCIEN>"
          stop 1
       ENDIF
    ENDIF

    ! Lecture ratqs

    ierr = NF_INQ_VARID (nid, "RATQS", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, "phyetat0: Le champ <RATQS> est absent"
       PRINT*, "Depart legerement fausse. Mais je continue"
       ratqs = 0.
    ELSE
       ierr = NF_GET_VAR_REAL(nid, nvarid, ratqs)
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Lecture echouee pour <RATQS>"
          stop 1
       ENDIF
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    xmin = MINval(ratqs)
    xmax = MAXval(ratqs)
    PRINT*,'(ecart-type) ratqs:', xmin, xmax

    ! Lecture run_off_lic_0

    ierr = NF_INQ_VARID (nid, "RUNOFFLIC0", nvarid)
    IF (ierr.NE.NF_NOERR) THEN
       PRINT*, "phyetat0: Le champ <RUNOFFLIC0> est absent"
       PRINT*, "Depart legerement fausse. Mais je continue"
       run_off_lic_0 = 0.
    ELSE
       ierr = NF_GET_VAR_REAL(nid, nvarid, run_off_lic_0)
       IF (ierr.NE.NF_NOERR) THEN
          PRINT*, "phyetat0: Lecture echouee pour <RUNOFFLIC0>"
          stop 1
       ENDIF
    ENDIF
    xmin = 1.0E+20
    xmax = -1.0E+20
    xmin = MINval(run_off_lic_0)
    xmax = MAXval(run_off_lic_0)
    PRINT*,'(ecart-type) run_off_lic_0:', xmin, xmax

    ! Fermer le fichier:

    ierr = NF_CLOSE(nid)

  END SUBROUTINE phyetat0

end module phyetat0_m
