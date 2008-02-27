!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/phyredem.F,v 1.3 2005/05/25 13:10:09 fairhead Exp $
!
c
      SUBROUTINE phyredem (fichnom,dtime,radpas,
     .           rlat,rlon, pctsrf,tsol,tsoil,
cIM "slab" ocean
     .           tslab,seaice,
     .           qsurf,qsol,snow,
     .           albedo, alblw, evap, rain_fall, snow_fall,
     .           solsw, sollw,fder,
     .           radsol,frugs,agesno,
     .           zmea,zstd,zsig,zgam,zthe,zpic,zval,rugsrel,
     .           t_ancien, q_ancien, rnebcon, ratqs, clwcon,
     .           run_off_lic_0)
      use dimens_m
      use indicesol
      use dimphy
      use conf_gcm_m
      use dimsoil
      use temps
      use clesphys
      IMPLICIT none
c======================================================================
c Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
c Objet: Ecriture de l'etat de redemarrage pour la physique
c======================================================================
      include "netcdf.inc"
c======================================================================
      CHARACTER*(*) fichnom
      REAL dtime
      INTEGER radpas
      REAL, intent(in):: rlat(klon), rlon(klon)
      REAL tsol(klon,nbsrf)
      REAL tsoil(klon,nsoilmx,nbsrf)
cIM "slab" ocean
      REAL tslab(klon), seaice(klon)
      REAL qsurf(klon,nbsrf)
      REAL qsol(klon)
      REAL snow(klon,nbsrf)
      REAL albedo(klon,nbsrf)
cIM BEG
      REAL alblw(klon,nbsrf)
cIM END
      REAL evap(klon,nbsrf)
      REAL rain_fall(klon)
      REAL snow_fall(klon)
      real solsw(klon)
      real sollw(klon)
      real fder(klon)
      REAL radsol(klon)
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
      REAL t_ancien(klon,klev), q_ancien(klon,klev)
      real clwcon(klon,klev),rnebcon(klon,klev),ratqs(klon,klev)
      REAL run_off_lic_0(klon)
c
      INTEGER nid, nvarid, idim1, idim2, idim3
      INTEGER ierr
      INTEGER length
      PARAMETER (length=100)
      REAL tab_cntrl(length)
c
      INTEGER isoil, nsrf
      CHARACTER*7 str7
      CHARACTER*2 str2
c
      print *, "Call sequence information: phyredem"
      ierr = NF_CREATE(fichnom, NF_CLOBBER, nid)
      IF (ierr.NE.NF_NOERR) THEN
        write(6,*)' Pb d''ouverture du fichier '//fichnom
        write(6,*)' ierr = ', ierr
        STOP 1
      ENDIF
c
      ierr = NF_PUT_ATT_TEXT (nid, NF_GLOBAL, "title", 28,
     .                       "Fichier redemmarage physique")
c
      ierr = NF_DEF_DIM (nid, "index", length, idim1)
      ierr = NF_DEF_DIM (nid, "points_physiques", klon, idim2)
      ierr = NF_DEF_DIM (nid, "horizon_vertical", klon*klev, idim3)
c
      ierr = NF_ENDDEF(nid)
c
      DO ierr = 1, length
         tab_cntrl(ierr) = 0.0
      ENDDO
      tab_cntrl(1) = dtime
      tab_cntrl(2) = radpas
      tab_cntrl(3) = co2_ppm
      tab_cntrl(4) = solaire
      tab_cntrl(5) = iflag_con
      tab_cntrl(6) = nbapp_rad

      IF( cycle_diurne ) tab_cntrl( 7 ) = 1.
      IF(   soil_model ) tab_cntrl( 8 ) = 1.
      IF(     new_oliq ) tab_cntrl( 9 ) = 1.
      IF(     ok_orodr ) tab_cntrl(10 ) = 1.
      IF(     ok_orolf ) tab_cntrl(11 ) = 1.

      tab_cntrl(13) = day_end
      tab_cntrl(14) = annee_ref
      tab_cntrl(15) = itau_phy
c
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "controle", NF_FLOAT, 1, idim1,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 22,
     .                        "Parametres de controle")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,tab_cntrl)
c
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "longitude", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 32,
     .                        "Longitudes de la grille physique")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlon)
c
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "latitude", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 31,
     .                        "Latitudes de la grille physique")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlat)
c
C PB ajout du masque terre/mer
C
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "masque", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 16,
     .                        "masque terre mer")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zmasq)
c BP ajout des fraction de chaque sous-surface
C
C 1. fraction de terre 
C
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "FTER", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 21,
     .                        "fraction de continent")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,pctsrf(1 : klon, is_ter))
C 
C 2. Fraction de glace de terre
C 
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "FLIC", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 24,
     .                        "fraction glace de terre")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,pctsrf(1 : klon, is_lic))
C
C 3. fraction ocean
C
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "FOCE", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 14,
     .                        "fraction ocean")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,pctsrf(1 : klon, is_oce))
C
C 4. Fraction glace de mer
C
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "FSIC", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 18,
     .                        "fraction glace mer")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,pctsrf(1 : klon, is_sic))
C
C
c
      DO nsrf = 1, nbsrf
        IF (nsrf.LE.99) THEN
        WRITE(str2,'(i2.2)') nsrf
        ierr = NF_REDEF (nid)
        ierr = NF_DEF_VAR (nid, "TS"//str2, NF_FLOAT, 1, idim2,nvarid)
        ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 28,
     .                        "Temperature de surface No."//str2)
        ierr = NF_ENDDEF(nid)
        ELSE
        PRINT*, "Trop de sous-mailles"
        stop 1
        ENDIF
        ierr = NF_PUT_VAR_REAL (nid,nvarid,tsol(1,nsrf))
      ENDDO
c
      DO nsrf = 1, nbsrf
      DO isoil=1, nsoilmx
        IF (isoil.LE.99 .AND. nsrf.LE.99) THEN
        WRITE(str7,'(i2.2,"srf",i2.2)') isoil,nsrf
        ierr = NF_REDEF (nid)
        ierr = NF_DEF_VAR (nid, "Tsoil"//str7,NF_FLOAT,1,idim2,nvarid)
        ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 29,
     .                        "Temperature du sol No."//str7)
        ierr = NF_ENDDEF(nid)
        ELSE
        PRINT*, "Trop de couches"
        stop 1
        ENDIF
        ierr = NF_PUT_VAR_REAL (nid,nvarid,tsoil(1,isoil,nsrf))
      ENDDO
      ENDDO
c
cIM "slab" ocean
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "TSLAB", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 33,
     .                        "Ecart de la SST (pour slab-ocean)")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,tslab)
c
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "SEAICE", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 33,
     .                        "Glace de mer kg/m2 (pour slab-ocean)")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,seaice)
c
      DO nsrf = 1, nbsrf
        IF (nsrf.LE.99) THEN
        WRITE(str2,'(i2.2)') nsrf
        ierr = NF_REDEF (nid)
        ierr = NF_DEF_VAR (nid,"QS"//str2,NF_FLOAT,1,idim2,nvarid)
        ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 25,
     .                        "Humidite de surface No."//str2)
        ierr = NF_ENDDEF(nid)
        ELSE
        PRINT*, "Trop de sous-mailles"
        stop 1
        ENDIF
      ierr = NF_PUT_VAR_REAL (nid,nvarid,qsurf(1,nsrf))
      END DO
C
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid,"QSOL",NF_FLOAT,1,idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 20,
     .    "Eau dans le sol (mm)")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,qsol)
c
      DO nsrf = 1, nbsrf
        IF (nsrf.LE.99) THEN
        WRITE(str2,'(i2.2)') nsrf
        ierr = NF_REDEF (nid)
        ierr = NF_DEF_VAR (nid,"ALBE"//str2,NF_FLOAT,1,idim2,nvarid)
        ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 23,
     .                        "albedo de surface No."//str2)
        ierr = NF_ENDDEF(nid)
        ELSE
        PRINT*, "Trop de sous-mailles"
        stop 1
        ENDIF
      ierr = NF_PUT_VAR_REAL (nid,nvarid,albedo(1,nsrf))
      ENDDO

cIM BEG albedo LW
        DO nsrf = 1, nbsrf
        IF (nsrf.LE.99) THEN
        WRITE(str2,'(i2.2)') nsrf
        ierr = NF_REDEF (nid)
        ierr = NF_DEF_VAR (nid,"ALBLW"//str2,NF_FLOAT,1,idim2,nvarid)
        ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 23,
     .                        "albedo LW de surface No."//str2)
        ierr = NF_ENDDEF(nid)
        ELSE
        PRINT*, "Trop de sous-mailles"
        stop 1
        ENDIF
      ierr = NF_PUT_VAR_REAL (nid,nvarid,alblw(1,nsrf))
      ENDDO
cIM END albedo LW
c
      DO nsrf = 1, nbsrf
        IF (nsrf.LE.99) THEN
        WRITE(str2,'(i2.2)') nsrf
        ierr = NF_REDEF (nid)
        ierr = NF_DEF_VAR (nid,"EVAP"//str2,NF_FLOAT,1,idim2,nvarid)
        ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 28,
     .                        "Evaporation de surface No."//str2)
        ierr = NF_ENDDEF(nid)
        ELSE
        PRINT*, "Trop de sous-mailles"
        stop 1
        ENDIF
      ierr = NF_PUT_VAR_REAL (nid,nvarid,evap(1,nsrf))
      ENDDO

c
      DO nsrf = 1, nbsrf
        IF (nsrf.LE.99) THEN
        WRITE(str2,'(i2.2)') nsrf
        ierr = NF_REDEF (nid)
        ierr = NF_DEF_VAR (nid,"SNOW"//str2,NF_FLOAT,1,idim2,nvarid)
        ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 22,
     .                        "Neige de surface No."//str2)
        ierr = NF_ENDDEF(nid)
        ELSE
        PRINT*, "Trop de sous-mailles"
        stop 1
        ENDIF
      ierr = NF_PUT_VAR_REAL (nid,nvarid,snow(1,nsrf))
      ENDDO

c
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "RADS", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 28,
     .                        "Rayonnement net a la surface")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,radsol)
c
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "solsw", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 32,
     .                        "Rayonnement solaire a la surface")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,solsw)
c
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "sollw", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 27,
     .                        "Rayonnement IF a la surface")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,sollw)
c
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "fder", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 14,
     .                        "Derive de flux")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,fder)
c
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "rain_f", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 21,
     .                        "precipitation liquide")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rain_fall)
c
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "snow_f", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 20,
     .                        "precipitation solide")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,snow_fall)
c
      DO nsrf = 1, nbsrf
        IF (nsrf.LE.99) THEN
        WRITE(str2,'(i2.2)') nsrf
        ierr = NF_REDEF (nid)
        ierr = NF_DEF_VAR (nid,"RUG"//str2,NF_FLOAT,1,idim2,nvarid)
        ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 23,
     .                        "rugosite de surface No."//str2)
        ierr = NF_ENDDEF(nid)
        ELSE
        PRINT*, "Trop de sous-mailles"
        stop 1
        ENDIF
      ierr = NF_PUT_VAR_REAL (nid,nvarid,frugs(1,nsrf))
      ENDDO
c
      DO nsrf = 1, nbsrf
        IF (nsrf.LE.99) THEN
            WRITE(str2,'(i2.2)') nsrf
            ierr = NF_REDEF (nid)
            ierr = NF_DEF_VAR (nid,"AGESNO"//str2,NF_FLOAT,1,idim2
     $          ,nvarid)
            ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 15,
     .                        "Age de la neige surface No."//str2)
            ierr = NF_ENDDEF(nid)
        ELSE
            PRINT*, "Trop de sous-mailles"
            stop 1
        ENDIF
      ierr = NF_PUT_VAR_REAL (nid,nvarid,agesno(1,nsrf))
      ENDDO
c
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "ZMEA", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zmea)
c
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "ZSTD", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zstd)
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "ZSIG", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zsig)
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "ZGAM", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zgam)
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "ZTHE", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zthe)
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "ZPIC", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zpic)
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "ZVAL", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zval)
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "RUGSREL", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rugsrel)
c
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "TANCIEN", NF_FLOAT, 1, idim3,nvarid)
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,t_ancien)
c
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "QANCIEN", NF_FLOAT, 1, idim3,nvarid)
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,q_ancien)
c
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "RUGMER", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 28,
     .                        "Longueur de rugosite sur mer")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,frugs(1,is_oce))
c
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "CLWCON", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 28,
     .                        "Eau liquide convective")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,clwcon)
c
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "RNEBCON", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 28,
     .                        "Nebulosite convective")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rnebcon)
c
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid, "RATQS", NF_FLOAT, 1, idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT(nid,nvarid,"title", 5,
     .                        "Ratqs")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,ratqs)
c
c run_off_lic_0
c
      ierr = NF_REDEF (nid)
      ierr=NF_DEF_VAR(nid,"RUNOFFLIC0",NF_FLOAT, 1,idim2,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 10,
     .                        "Runofflic0")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_REAL (nid,nvarid,run_off_lic_0)
c
c
      ierr = NF_CLOSE(nid)
c
      RETURN
      END
