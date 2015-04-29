SUBROUTINE nuage (paprs, pplay, &
     t, pqlwp, pclc, pcltau, pclemi, &
     pch, pcl, pcm, pct, pctlwp, &
     ok_aie, &
     sulfate, sulfate_pi,  &
     bl95_b0, bl95_b1, &
     cldtaupi, re, fl)
  !
  ! From LMDZ4/libf/phylmd/nuage.F, version 1.1.1.1 2004/05/19 12:53:07
  !
  use dimens_m
  use dimphy
  use SUPHEC_M
  IMPLICIT none
  !======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) date: 19930910
  ! Objet: Calculer epaisseur optique et emmissivite des nuages
  !======================================================================
  ! Arguments:
  ! t-------input-R-temperature
  ! pqlwp---input-R-eau liquide nuageuse dans l'atmosphere (kg/kg)
  ! pclc----input-R-couverture nuageuse pour le rayonnement (0 a 1)
  ! ok_aie--input-L-apply aerosol indirect effect or not
  ! sulfate-input-R-sulfate aerosol mass concentration [um/m^3]
  ! sulfate_pi-input-R-dito, pre-industrial value
  ! bl95_b0-input-R-a parameter, may be varied for tests (s-sea, l-land)
  ! bl95_b1-input-R-a parameter, may be varied for tests (    -"-      )
  !
  ! cldtaupi-output-R-pre-industrial value of cloud optical thickness,
  !                   needed for the diagnostics of the aerosol indirect
  !                   radiative forcing (see radlwsw)
  ! re------output-R-Cloud droplet effective radius multiplied by fl [um]
  ! fl------output-R-Denominator to re, introduced to avoid problems in
  !                  the averaging of the output. fl is the fraction of liquid
  !                  water clouds within a grid cell
  !
  ! pcltau--output-R-epaisseur optique des nuages
  ! pclemi--output-R-emissivite des nuages (0 a 1)
  !======================================================================
  !
  !
  REAL, intent(in):: paprs(klon,klev+1)
  real, intent(in):: pplay(klon,klev)
  REAL, intent(in):: t(klon,klev)
  !
  REAL pclc(klon,klev)
  REAL pqlwp(klon,klev)
  REAL pcltau(klon,klev), pclemi(klon,klev)
  !
  REAL pct(klon), pctlwp(klon), pch(klon), pcl(klon), pcm(klon)
  !
  LOGICAL lo
  !
  REAL cetahb, cetamb
  PARAMETER (cetahb = 0.45, cetamb = 0.80)
  !
  INTEGER i, k
  REAL zflwp, zfice
  !
  REAL radius, rad_froid, rad_chaud, rad_chau1, rad_chau2
  PARAMETER (rad_chau1=13.0, rad_chau2=9.0, rad_froid=35.0)
  !cc      PARAMETER (rad_chaud=15.0, rad_froid=35.0)
  ! sintex initial      PARAMETER (rad_chaud=10.0, rad_froid=30.0)
  REAL coef, coef_froi, coef_chau
  PARAMETER (coef_chau=0.13, coef_froi=0.09)
  REAL seuil_neb, t_glace
  PARAMETER (seuil_neb=0.001, t_glace=273.0-15.0)
  INTEGER nexpo ! exponentiel pour glace/eau
  PARAMETER (nexpo=6)

  !jq for the aerosol indirect effect
  !jq introduced by Johannes Quaas (quaas@lmd.jussieu.fr), 27/11/2003
  !jq
  LOGICAL ok_aie            ! Apply AIE or not?

  REAL sulfate(klon, klev)  ! sulfate aerosol mass concentration [ug m-3]
  REAL cdnc(klon, klev)     ! cloud droplet number concentration [m-3]
  REAL re(klon, klev)       ! cloud droplet effective radius [um]
  REAL sulfate_pi(klon, klev)  ! sulfate aerosol mass concentration [ug m-3] (pre-industrial value)
  REAL cdnc_pi(klon, klev)     ! cloud droplet number concentration [m-3] (pi value)

  REAL fl(klon, klev)  ! xliq * rneb (denominator to re ; fraction of liquid water clouds within the grid cell)

  REAL bl95_b0, bl95_b1     ! Parameter in B&L 95-Formula

  REAL cldtaupi(klon, klev) ! pre-industrial cloud opt thickness for diag

  !cc      PARAMETER (nexpo=1)
  !
  ! Calculer l'epaisseur optique et l'emmissivite des nuages
  !
  DO k = 1, klev
     DO i = 1, klon
        rad_chaud = rad_chau1
        IF (k.LE.3) rad_chaud = rad_chau2

        pclc(i,k) = MAX(pclc(i,k), seuil_neb)
        zflwp = 1000.*pqlwp(i,k)/RG/pclc(i,k) &
             *(paprs(i,k)-paprs(i,k+1))
        zfice = 1.0 - (t(i,k)-t_glace) / (273.13-t_glace)
        zfice = MIN(MAX(zfice,0.0),1.0)
        zfice = zfice**nexpo

        IF (ok_aie) THEN
           ! Formula "D" of Boucher and Lohmann, Tellus, 1995
           !
           cdnc(i,k) = 10.**(bl95_b0+bl95_b1* &
                log(MAX(sulfate(i,k),1.e-4))/log(10.))*1.e6 !-m-3
           ! Cloud droplet number concentration (CDNC) is restricted
           ! to be within [20, 1000 cm^3]
           !
           cdnc(i,k)=MIN(1000.e6,MAX(20.e6,cdnc(i,k)))
           cdnc_pi(i,k) = 10.**(bl95_b0+bl95_b1* &
                log(MAX(sulfate_pi(i,k),1.e-4))/log(10.))*1.e6 !-m-3
           cdnc_pi(i,k)=MIN(1000.e6,MAX(20.e6,cdnc_pi(i,k)))
           !
           !
           ! air density: pplay(i,k) / (RD * zT(i,k))
           ! factor 1.1: derive effective radius from volume-mean radius
           ! factor 1000 is the water density
           ! _chaud means that this is the CDR for liquid water clouds
           !
           rad_chaud =  &
                1.1 * ( (pqlwp(i,k) * pplay(i,k) / (RD * T(i,k)) )   &
                / (4./3. * RPI * 1000. * cdnc(i,k)) )**(1./3.)
           !
           ! Convert to um. CDR shall be at least 3 um.
           !
           rad_chaud = MAX(rad_chaud*1.e6, 3.)

           ! For output diagnostics
           !
           ! Cloud droplet effective radius [um]
           !
           ! we multiply here with f * xl (fraction of liquid water
           ! clouds in the grid cell) to avoid problems in the
           ! averaging of the output.
           ! In the output of IOIPSL, derive the real cloud droplet
           ! effective radius as re/fl
           !
           fl(i,k) = pclc(i,k)*(1.-zfice)
           re(i,k) = rad_chaud*fl(i,k)

           ! Pre-industrial cloud opt thickness
           !
           ! "radius" is calculated as rad_chaud above (plus the
           ! ice cloud contribution) but using cdnc_pi instead of
           ! cdnc.
           radius = MAX(1.1e6 * ( (pqlwp(i,k)*pplay(i,k)/(RD*T(i,k)))   &
                / (4./3.*RPI*1000.*cdnc_pi(i,k)) )**(1./3.),  &
                3.) * (1.-zfice) + rad_froid * zfice
           cldtaupi(i,k) = 3.0/2.0 * zflwp / radius

        END IF                  ! ok_aie

        radius = rad_chaud * (1.-zfice) + rad_froid * zfice
        coef = coef_chau * (1.-zfice) + coef_froi * zfice
        pcltau(i,k) = 3.0/2.0 * zflwp / radius
        pclemi(i,k) = 1.0 - EXP( - coef * zflwp)
        lo = (pclc(i,k) .LE. seuil_neb)
        IF (lo) pclc(i,k) = 0.0
        IF (lo) pcltau(i,k) = 0.0
        IF (lo) pclemi(i,k) = 0.0

        IF (.NOT.ok_aie) cldtaupi(i,k)=pcltau(i,k)
     END DO
  END DO
  !
  ! COMPUTE CLOUD LIQUID PATH AND TOTAL CLOUDINESS
  !
  DO i = 1, klon
     pct(i)=1.0
     pch(i)=1.0
     pcm(i) = 1.0
     pcl(i) = 1.0
     pctlwp(i) = 0.0
  END DO
  !
  DO k = klev, 1, -1
     DO i = 1, klon
        pctlwp(i) = pctlwp(i)  &
             + pqlwp(i,k)*(paprs(i,k)-paprs(i,k+1))/RG
        pct(i) = pct(i)*(1.0-pclc(i,k))
        if (pplay(i,k).LE.cetahb*paprs(i,1)) &
             pch(i) = pch(i)*(1.0-pclc(i,k))
        if (pplay(i,k).GT.cetahb*paprs(i,1) .AND. &
             pplay(i,k).LE.cetamb*paprs(i,1))  &
             pcm(i) = pcm(i)*(1.0-pclc(i,k))
        if (pplay(i,k).GT.cetamb*paprs(i,1)) &
             pcl(i) = pcl(i)*(1.0-pclc(i,k))
     END DO
  END DO
  !
  DO i = 1, klon
     pct(i)=1.-pct(i)
     pch(i)=1.-pch(i)
     pcm(i)=1.-pcm(i)
     pcl(i)=1.-pcl(i)
  END DO
  !
END SUBROUTINE nuage
