!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/radlwsw.F,v 1.4 2005/06/06 13:16:33 fairhead Exp $
!
      SUBROUTINE radlwsw(dist, rmu0, fract, 
     .                  paprs, pplay,tsol,albedo, alblw, t,q,wo,
     .                  cldfra, cldemi, cldtaupd,
     .                  heat,heat0,cool,cool0,radsol,albpla,
     .                  topsw,toplw,solsw,sollw,
     .                  sollwdown,
     .                  topsw0,toplw0,solsw0,sollw0,
     .                  lwdn0, lwdn, lwup0, lwup,
     .                  swdn0, swdn, swup0, swup,
     .                  ok_ade, ok_aie,
     .                  tau_ae, piz_ae, cg_ae,
     .                  topswad, solswad,
     .                  cldtaupi, topswai, solswai)
c      
      use dimphy
      use clesphys 
      use YOMCST
      use raddim, only: kflev, kdlon
      use yoethf
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19960719
c Objet: interface entre le modele et les rayonnements
c Arguments:
c dist-----input-R- distance astronomique terre-soleil
c rmu0-----input-R- cosinus de l'angle zenithal
c fract----input-R- duree d'ensoleillement normalisee
c co2_ppm--input-R- concentration du gaz carbonique (en ppm)
c solaire--input-R- constante solaire (W/m**2)
c paprs----input-R- pression a inter-couche (Pa)
c pplay----input-R- pression au milieu de couche (Pa)
c tsol-----input-R- temperature du sol (en K)
c albedo---input-R- albedo du sol (entre 0 et 1)
c t--------input-R- temperature (K)
c q--------input-R- vapeur d'eau (en kg/kg)
c wo-------input-R- contenu en ozone (en kg/kg) correction MPL 100505
c cldfra---input-R- fraction nuageuse (entre 0 et 1)
c cldtaupd---input-R- epaisseur optique des nuages dans le visible (present-day value)
c cldemi---input-R- emissivite des nuages dans l'IR (entre 0 et 1)
c ok_ade---input-L- apply the Aerosol Direct Effect or not?
c ok_aie---input-L- apply the Aerosol Indirect Effect or not?
c tau_ae, piz_ae, cg_ae-input-R- aerosol optical properties (calculated in aeropt.F)
c cldtaupi-input-R- epaisseur optique des nuages dans le visible
c                   calculated for pre-industrial (pi) aerosol concentrations, i.e. with smaller
c                   droplet concentration, thus larger droplets, thus generally cdltaupi cldtaupd
c                   it is needed for the diagnostics of the aerosol indirect radiative forcing      
c
c heat-----output-R- echauffement atmospherique (visible) (K/jour)
c cool-----output-R- refroidissement dans l'IR (K/jour)
c radsol---output-R- bilan radiatif net au sol (W/m**2) (+ vers le bas)
c albpla---output-R- albedo planetaire (entre 0 et 1)
c topsw----output-R- flux solaire net au sommet de l'atm.
c toplw----output-R- ray. IR montant au sommet de l'atmosphere
c solsw----output-R- flux solaire net a la surface
c sollw----output-R- ray. IR montant a la surface
c solswad---output-R- ray. solaire net absorbe a la surface (aerosol dir)
c topswad---output-R- ray. solaire absorbe au sommet de l'atm. (aerosol dir)
c solswai---output-R- ray. solaire net absorbe a la surface (aerosol ind)
c topswai---output-R- ray. solaire absorbe au sommet de l'atm. (aerosol ind)
c
c ATTENTION: swai and swad have to be interpreted in the following manner:
c ---------
c ok_ade=F & ok_aie=F -both are zero
c ok_ade=T & ok_aie=F -aerosol direct forcing is F_{AD} = topsw-topswad
c                        indirect is zero
c ok_ade=F & ok_aie=T -aerosol indirect forcing is F_{AI} = topsw-topswai
c                        direct is zero
c ok_ade=T & ok_aie=T -aerosol indirect forcing is F_{AI} = topsw-topswai
c                        aerosol direct forcing is F_{AD} = topswai-topswad
c
      
c======================================================================
c
      real rmu0(klon), fract(klon), dist
cIM   real co2_ppm
cIM   real solaire
c
      real, intent(in):: paprs(klon,klev+1)
      real pplay(klon,klev)
      real albedo(klon), alblw(klon), tsol(klon)
      real t(klon,klev), q(klon,klev)
      real, intent(in):: wo(klon,klev)
      real cldfra(klon,klev), cldemi(klon,klev), cldtaupd(klon,klev)
      real heat(klon,klev), cool(klon,klev)
      real heat0(klon,klev), cool0(klon,klev)
      real radsol(klon), topsw(klon), toplw(klon)
      real solsw(klon), sollw(klon), albpla(klon)
      real topsw0(klon), toplw0(klon), solsw0(klon), sollw0(klon)
      real sollwdown(klon)
cIM output 3D 
      REAL*8 ZFSUP(KDLON,KFLEV+1)
      REAL*8 ZFSDN(KDLON,KFLEV+1)
      REAL*8 ZFSUP0(KDLON,KFLEV+1)
      REAL*8 ZFSDN0(KDLON,KFLEV+1)
c
      REAL*8 ZFLUP(KDLON,KFLEV+1)
      REAL*8 ZFLDN(KDLON,KFLEV+1)
      REAL*8 ZFLUP0(KDLON,KFLEV+1)
      REAL*8 ZFLDN0(KDLON,KFLEV+1)
c
      REAL*8 zx_alpha1, zx_alpha2
c
c
      INTEGER k, kk, i, j, iof, nb_gr
      EXTERNAL lw, sw
c
cIM ctes ds clesphys.h  REAL*8 RCO2, RCH4, RN2O, RCFC11, RCFC12
      REAL*8 PSCT
c
      REAL*8 PALBD(kdlon,2), PALBP(kdlon,2)
      REAL*8 PEMIS(kdlon), PDT0(kdlon), PVIEW(kdlon)
      REAL*8 PPSOL(kdlon), PDP(kdlon,klev)
      REAL*8 PTL(kdlon,kflev+1), PPMB(kdlon,kflev+1)
      REAL*8 PTAVE(kdlon,kflev)
      REAL*8 PWV(kdlon,kflev), PQS(kdlon,kflev), POZON(kdlon,kflev)
      REAL*8 PAER(kdlon,kflev,5)
      REAL*8 PCLDLD(kdlon,kflev)
      REAL*8 PCLDLU(kdlon,kflev)
      REAL*8 PCLDSW(kdlon,kflev)
      REAL*8 PTAU(kdlon,2,kflev)
      REAL*8 POMEGA(kdlon,2,kflev)
      REAL*8 PCG(kdlon,2,kflev)
c
      REAL*8 zfract(kdlon), zrmu0(kdlon), zdist
c
      REAL*8 zheat(kdlon,kflev), zcool(kdlon,kflev)
      REAL*8 zheat0(kdlon,kflev), zcool0(kdlon,kflev)
      REAL*8 ztopsw(kdlon), ztoplw(kdlon)
      REAL*8 zsolsw(kdlon), zsollw(kdlon), zalbpla(kdlon)
cIM
      REAL*8 zsollwdown(kdlon)
c
      REAL*8 ztopsw0(kdlon), ztoplw0(kdlon)
      REAL*8 zsolsw0(kdlon), zsollw0(kdlon)
      REAL*8 zznormcp
cIM output 3D : SWup, SWdn, LWup, LWdn
      REAL swdn(klon,kflev+1),swdn0(klon,kflev+1)
      REAL swup(klon,kflev+1),swup0(klon,kflev+1)
      REAL lwdn(klon,kflev+1),lwdn0(klon,kflev+1)
      REAL lwup(klon,kflev+1),lwup0(klon,kflev+1)
c-OB
cjq the following quantities are needed for the aerosol radiative forcings

      real topswad(klon), solswad(klon) ! output: aerosol direct forcing at TOA and surface
      real topswai(klon), solswai(klon) ! output: aerosol indirect forcing atTOA and surface
      real tau_ae(klon,klev,2), piz_ae(klon,klev,2), cg_ae(klon,klev,2) ! aerosol optical properties (see aeropt.F)
      real cldtaupi(klon,klev)  ! cloud optical thickness for pre-industrial aerosol concentrations
                                ! (i.e., with a smaller droplet concentrationand thus larger droplet radii)
      logical ok_ade, ok_aie    ! switches whether to use aerosol direct (indirect) effects or not
      real*8 tauae(kdlon,kflev,2) ! aer opt properties
      real*8 pizae(kdlon,kflev,2)
      real*8 cgae(kdlon,kflev,2)
      REAL*8 PTAUA(kdlon,2,kflev) ! present-day value of cloud opt thickness (PTAU is pre-industrial value), local use
      REAL*8 POMEGAA(kdlon,2,kflev) ! dito for single scatt albedo
      REAL*8 ztopswad(kdlon), zsolswad(kdlon) ! Aerosol direct forcing at TOAand surface
      REAL*8 ztopswai(kdlon), zsolswai(kdlon) ! dito, indirect
cjq-end
!rv
      tauae(:,:,:)=0.
      pizae(:,:,:)=0.
      cgae(:,:,:)=0.
!rv
      
c
c-------------------------------------------
      nb_gr = klon / kdlon
      IF (nb_gr*kdlon .NE. klon) THEN
         PRINT*, "kdlon mauvais:", klon, kdlon, nb_gr
         stop 1
      ENDIF
      IF (kflev .NE. klev) THEN
          PRINT*, "kflev differe de klev, kflev, klev"
          stop 1
      ENDIF
c-------------------------------------------
      DO k = 1, klev
      DO i = 1, klon
         heat(i,k)=0.
         cool(i,k)=0.
         heat0(i,k)=0.
         cool0(i,k)=0.
      ENDDO
      ENDDO
c
      zdist = dist
c
cIM anciennes valeurs
c     RCO2 = co2_ppm * 1.0e-06  * 44.011/28.97
c
cIM : on met RCO2, RCH4, RN2O, RCFC11 et RCFC12 dans clesphys.h /lecture ds conf_phys.F90
c     RCH4 = 1.65E-06* 16.043/28.97
c     RN2O = 306.E-09* 44.013/28.97
c     RCFC11 = 280.E-12* 137.3686/28.97
c     RCFC12 = 484.E-12* 120.9140/28.97
cIM anciennes valeurs
c     RCH4 = 1.72E-06* 16.043/28.97
c     RN2O = 310.E-09* 44.013/28.97
c
c     PRINT*,'IMradlwsw : solaire, co2= ', solaire, co2_ppm
      PSCT = solaire/zdist/zdist
c
      DO 99999 j = 1, nb_gr
      iof = kdlon*(j-1)
c
      DO i = 1, kdlon
         zfract(i) = fract(iof+i)
         zrmu0(i) = rmu0(iof+i)
         PALBD(i,1) = albedo(iof+i)
!         PALBD(i,2) = albedo(iof+i)
         PALBD(i,2) = alblw(iof+i)
         PALBP(i,1) = albedo(iof+i)
!         PALBP(i,2) = albedo(iof+i)
         PALBP(i,2) = alblw(iof+i)
cIM cf. JLD pour etre en accord avec ORCHIDEE il faut mettre PEMIS(i) = 0.96
         PEMIS(i) = 1.0 
         PVIEW(i) = 1.66
         PPSOL(i) = paprs(iof+i,1)
         zx_alpha1 = (paprs(iof+i,1)-pplay(iof+i,2)) 
     .             / (pplay(iof+i,1)-pplay(iof+i,2))
         zx_alpha2 = 1.0 - zx_alpha1
         PTL(i,1) = t(iof+i,1) * zx_alpha1 + t(iof+i,2) * zx_alpha2
         PTL(i,klev+1) = t(iof+i,klev)
         PDT0(i) = tsol(iof+i) - PTL(i,1)
      ENDDO
      DO k = 2, kflev
      DO i = 1, kdlon
         PTL(i,k) = (t(iof+i,k)+t(iof+i,k-1))*0.5
      ENDDO
      ENDDO
      DO k = 1, kflev
      DO i = 1, kdlon
         PDP(i,k) = paprs(iof+i,k)-paprs(iof+i,k+1)
         PTAVE(i,k) = t(iof+i,k)
         PWV(i,k) = MAX (q(iof+i,k), 1.0e-12)
         PQS(i,k) = PWV(i,k)
c wo:    cm.atm (epaisseur en cm dans la situation standard)
c POZON: kg/kg
         IF (bug_ozone) then
           POZON(i,k) = MAX(wo(iof+i,k),1.0e-12)*RG/46.6968
     .               /(paprs(iof+i,k)-paprs(iof+i,k+1))
     .               *(paprs(iof+i,1)/101325.0)
         ELSE
c le calcul qui suit est maintenant fait dans ozonecm (MPL)
           POZON(i,k) = wo(i,k)
         ENDIF
         PCLDLD(i,k) = cldfra(iof+i,k)*cldemi(iof+i,k)
         PCLDLU(i,k) = cldfra(iof+i,k)*cldemi(iof+i,k)
         PCLDSW(i,k) = cldfra(iof+i,k)
         PTAU(i,1,k) = MAX(cldtaupi(iof+i,k), 1.0e-05)! 1e-12 serait instable
         PTAU(i,2,k) = MAX(cldtaupi(iof+i,k), 1.0e-05)! pour 32-bit machines
         POMEGA(i,1,k) = 0.9999 - 5.0e-04 * EXP(-0.5 * PTAU(i,1,k))
         POMEGA(i,2,k) = 0.9988 - 2.5e-03 * EXP(-0.05 * PTAU(i,2,k))
         PCG(i,1,k) = 0.865
         PCG(i,2,k) = 0.910
c-OB
cjq Introduced for aerosol indirect forcings.
cjq The following values use the cloud optical thickness calculated from
cjq present-day aerosol concentrations whereas the quantities without the
cjq "A" at the end are for pre-industial (natural-only) aerosol concentrations
cjq
         PTAUA(i,1,k) = MAX(cldtaupd(iof+i,k), 1.0e-05)! 1e-12 serait instable
         PTAUA(i,2,k) = MAX(cldtaupd(iof+i,k), 1.0e-05)! pour 32-bit machines
         POMEGAA(i,1,k) = 0.9999 - 5.0e-04 * EXP(-0.5 * PTAUA(i,1,k))
         POMEGAA(i,2,k) = 0.9988 - 2.5e-03 * EXP(-0.05 * PTAUA(i,2,k))
cjq-end
      ENDDO
      ENDDO
c
      DO k = 1, kflev+1
      DO i = 1, kdlon
         PPMB(i,k) = paprs(iof+i,k)/100.0
      ENDDO
      ENDDO
c
      DO kk = 1, 5
      DO k = 1, kflev
      DO i = 1, kdlon
         PAER(i,k,kk) = 1.0E-15
      ENDDO
      ENDDO
      ENDDO
c-OB
      DO k = 1, kflev
      DO i = 1, kdlon
        tauae(i,k,1)=tau_ae(iof+i,k,1)
        pizae(i,k,1)=piz_ae(iof+i,k,1)
        cgae(i,k,1) =cg_ae(iof+i,k,1)
        tauae(i,k,2)=tau_ae(iof+i,k,2)
        pizae(i,k,2)=piz_ae(iof+i,k,2)
        cgae(i,k,2) =cg_ae(iof+i,k,2)
      ENDDO
      ENDDO
c
c======================================================================
cIM ctes ds clesphys.h   CALL LW(RCO2,RCH4,RN2O,RCFC11,RCFC12,
      CALL LW(
     .        PPMB, PDP,
     .        PPSOL,PDT0,PEMIS,
     .        PTL, PTAVE, PWV, POZON, PAER,
     .        PCLDLD,PCLDLU,
     .        PVIEW,
     .        zcool, zcool0,
     .        ztoplw,zsollw,ztoplw0,zsollw0,
     .        zsollwdown,
     .        ZFLUP, ZFLDN, ZFLUP0,ZFLDN0)
cIM ctes ds clesphys.h   CALL SW(PSCT, RCO2, zrmu0, zfract,
      CALL SW(PSCT, zrmu0, zfract,
     S        PPMB, PDP,
     S        PPSOL, PALBD, PALBP,
     S        PTAVE, PWV, PQS, POZON, PAER,
     S        PCLDSW, PTAU, POMEGA, PCG,
     S        zheat, zheat0,
     S        zalbpla,ztopsw,zsolsw,ztopsw0,zsolsw0,
     S        ZFSUP,ZFSDN,ZFSUP0,ZFSDN0,
     S        tauae, pizae, cgae, ! aerosol optical properties
     s        PTAUA, POMEGAA,
     s        ztopswad,zsolswad,ztopswai,zsolswai, ! diagnosed aerosol forcing
     J        ok_ade, ok_aie) ! apply aerosol effects or not?

c======================================================================
      DO i = 1, kdlon
         radsol(iof+i) = zsolsw(i) + zsollw(i)
         topsw(iof+i) = ztopsw(i)
         toplw(iof+i) = ztoplw(i)
         solsw(iof+i) = zsolsw(i)
         sollw(iof+i) = zsollw(i)
         sollwdown(iof+i) = zsollwdown(i)
cIM
         DO k = 1, kflev+1
         lwdn0 ( iof+i,k)   = ZFLDN0 ( i,k)
         lwdn  ( iof+i,k)   = ZFLDN  ( i,k)
         lwup0 ( iof+i,k)   = ZFLUP0 ( i,k)
         lwup  ( iof+i,k)   = ZFLUP  ( i,k)
         ENDDO
c
         topsw0(iof+i) = ztopsw0(i)
         toplw0(iof+i) = ztoplw0(i)
         solsw0(iof+i) = zsolsw0(i)
         sollw0(iof+i) = zsollw0(i)
         albpla(iof+i) = zalbpla(i)
cIM
         DO k = 1, kflev+1
         swdn0 ( iof+i,k)   = ZFSDN0 ( i,k)
         swdn  ( iof+i,k)   = ZFSDN  ( i,k)
         swup0 ( iof+i,k)   = ZFSUP0 ( i,k)
         swup  ( iof+i,k)   = ZFSUP  ( i,k)
         ENDDO !k=1, kflev+1
      ENDDO
cjq-transform the aerosol forcings, if they have
cjq to be calculated
      IF (ok_ade) THEN
      DO i = 1, kdlon
         topswad(iof+i) = ztopswad(i)
         solswad(iof+i) = zsolswad(i)
      ENDDO
      ELSE
      DO i = 1, kdlon
         topswad(iof+i) = 0.0
         solswad(iof+i) = 0.0
      ENDDO
      ENDIF
      IF (ok_aie) THEN
      DO i = 1, kdlon
         topswai(iof+i) = ztopswai(i)
         solswai(iof+i) = zsolswai(i)
      ENDDO
      ELSE
      DO i = 1, kdlon
         topswai(iof+i) = 0.0
         solswai(iof+i) = 0.0
      ENDDO
      ENDIF
cjq-end
      DO k = 1, kflev
c      DO i = 1, kdlon
c         heat(iof+i,k) = zheat(i,k)
c         cool(iof+i,k) = zcool(i,k)
c         heat0(iof+i,k) = zheat0(i,k)
c         cool0(iof+i,k) = zcool0(i,k)
c      ENDDO
      DO i = 1, kdlon
C        scale factor to take into account the difference between
C        dry air and watter vapour scpecific heat capacity
         zznormcp=1.0+RVTMP2*PWV(i,k)
         heat(iof+i,k) = zheat(i,k)/zznormcp
         cool(iof+i,k) = zcool(i,k)/zznormcp
         heat0(iof+i,k) = zheat0(i,k)/zznormcp
         cool0(iof+i,k) = zcool0(i,k)/zznormcp
      ENDDO
      ENDDO
c
99999 CONTINUE
      RETURN
      END
cIM ctes ds clesphys.h   SUBROUTINE SW(PSCT, RCO2, PRMU0, PFRAC, 
      SUBROUTINE SW(PSCT, PRMU0, PFRAC, 
     S              PPMB, PDP, 
     S              PPSOL, PALBD, PALBP,
     S              PTAVE, PWV, PQS, POZON, PAER,
     S              PCLDSW, PTAU, POMEGA, PCG,
     S              PHEAT, PHEAT0,
     S              PALBPLA,PTOPSW,PSOLSW,PTOPSW0,PSOLSW0,
     S              ZFSUP,ZFSDN,ZFSUP0,ZFSDN0,
     S              tauae, pizae, cgae,
     s              PTAUA, POMEGAA,
     S              PTOPSWAD,PSOLSWAD,PTOPSWAI,PSOLSWAI,
     J              ok_ade, ok_aie )
      
      use dimens_m
      use dimphy
      use clesphys
      use YOMCST
      use raddim
      IMPLICIT none

C
C     ------------------------------------------------------------------
C
C     PURPOSE.
C     --------
C
C          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
C     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).
C
C     METHOD.
C     -------
C
C          1. COMPUTES ABSORBER AMOUNTS                 (SWU)
C          2. COMPUTES FLUXES IN 1ST SPECTRAL INTERVAL  (SW1S)
C          3. COMPUTES FLUXES IN 2ND SPECTRAL INTERVAL  (SW2S)
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
C        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C        95-01-01   J.-J. MORCRETTE  Direct/Diffuse Albedo
c        03-11-27   J. QUAAS Introduce aerosol forcings (based on BOUCHER)
C     ------------------------------------------------------------------
C
C* ARGUMENTS:
C
      REAL*8 PSCT  ! constante solaire (valeur conseillee: 1370)
cIM ctes ds clesphys.h   REAL*8 RCO2  ! concentration CO2 (IPCC: 353.E-06*44.011/28.97)
C
      REAL*8 PPSOL(KDLON)        ! SURFACE PRESSURE (PA)
      REAL*8 PDP(KDLON,KFLEV)    ! LAYER THICKNESS (PA)
      REAL*8 PPMB(KDLON,KFLEV+1) ! HALF-LEVEL PRESSURE (MB)
C
      REAL*8 PRMU0(KDLON)  ! COSINE OF ZENITHAL ANGLE
      REAL*8 PFRAC(KDLON)  ! fraction de la journee
C
      REAL*8 PTAVE(KDLON,KFLEV)  ! LAYER TEMPERATURE (K)
      REAL*8 PWV(KDLON,KFLEV)    ! SPECIFIC HUMIDITY (KG/KG)
      REAL*8 PQS(KDLON,KFLEV)    ! SATURATED WATER VAPOUR (KG/KG)
      REAL*8 POZON(KDLON,KFLEV)  ! OZONE CONCENTRATION (KG/KG)
      REAL*8 PAER(KDLON,KFLEV,5) ! AEROSOLS' OPTICAL THICKNESS
C
      REAL*8 PALBD(KDLON,2)  ! albedo du sol (lumiere diffuse)
      REAL*8 PALBP(KDLON,2)  ! albedo du sol (lumiere parallele)
C
      REAL*8 PCLDSW(KDLON,KFLEV)    ! CLOUD FRACTION
      REAL*8 PTAU(KDLON,2,KFLEV)    ! CLOUD OPTICAL THICKNESS
      REAL*8 PCG(KDLON,2,KFLEV)     ! ASYMETRY FACTOR
      REAL*8 POMEGA(KDLON,2,KFLEV)  ! SINGLE SCATTERING ALBEDO
C
      REAL*8 PHEAT(KDLON,KFLEV) ! SHORTWAVE HEATING (K/DAY)
      REAL*8 PHEAT0(KDLON,KFLEV)! SHORTWAVE HEATING (K/DAY) clear-sky
      REAL*8 PALBPLA(KDLON)     ! PLANETARY ALBEDO
      REAL*8 PTOPSW(KDLON)      ! SHORTWAVE FLUX AT T.O.A.
      REAL*8 PSOLSW(KDLON)      ! SHORTWAVE FLUX AT SURFACE
      REAL*8 PTOPSW0(KDLON)     ! SHORTWAVE FLUX AT T.O.A. (CLEAR-SKY)
      REAL*8 PSOLSW0(KDLON)     ! SHORTWAVE FLUX AT SURFACE (CLEAR-SKY)
C
C* LOCAL VARIABLES:
C
      REAL*8 ZOZ(KDLON,KFLEV)
      REAL*8 ZAKI(KDLON,2)     
      REAL*8 ZCLD(KDLON,KFLEV)
      REAL*8 ZCLEAR(KDLON) 
      REAL*8 ZDSIG(KDLON,KFLEV)
      REAL*8 ZFACT(KDLON)
      REAL*8 ZFD(KDLON,KFLEV+1)
      REAL*8 ZFDOWN(KDLON,KFLEV+1)
      REAL*8 ZFU(KDLON,KFLEV+1)
      REAL*8 ZFUP(KDLON,KFLEV+1)
      REAL*8 ZRMU(KDLON)
      REAL*8 ZSEC(KDLON)
      REAL*8 ZUD(KDLON,5,KFLEV+1)
      REAL*8 ZCLDSW0(KDLON,KFLEV)
c
      REAL*8 ZFSUP(KDLON,KFLEV+1)
      REAL*8 ZFSDN(KDLON,KFLEV+1)
      REAL*8 ZFSUP0(KDLON,KFLEV+1)
      REAL*8 ZFSDN0(KDLON,KFLEV+1)
C
      INTEGER inu, jl, jk, i, k, kpl1
c
      INTEGER swpas  ! Every swpas steps, sw is calculated
      PARAMETER(swpas=1)
c
      INTEGER itapsw
      LOGICAL appel1er
      DATA itapsw /0/
      DATA appel1er /.TRUE./
cjq-Introduced for aerosol forcings
      real*8 flag_aer
      logical ok_ade, ok_aie    ! use aerosol forcings or not?
      real*8 tauae(kdlon,kflev,2)  ! aerosol optical properties
      real*8 pizae(kdlon,kflev,2)  ! (see aeropt.F)
      real*8 cgae(kdlon,kflev,2)   ! -"-
      REAL*8 PTAUA(KDLON,2,KFLEV)    ! CLOUD OPTICAL THICKNESS (pre-industrial value)
      REAL*8 POMEGAA(KDLON,2,KFLEV)  ! SINGLE SCATTERING ALBEDO
      REAL*8 PTOPSWAD(KDLON)     ! SHORTWAVE FLUX AT T.O.A.(+AEROSOL DIR)
      REAL*8 PSOLSWAD(KDLON)     ! SHORTWAVE FLUX AT SURFACE(+AEROSOL DIR)
      REAL*8 PTOPSWAI(KDLON)     ! SHORTWAVE FLUX AT T.O.A.(+AEROSOL IND)
      REAL*8 PSOLSWAI(KDLON)     ! SHORTWAVE FLUX AT SURFACE(+AEROSOL IND)
cjq - Fluxes including aerosol effects
      REAL*8 ZFSUPAD(KDLON,KFLEV+1)
      REAL*8 ZFSDNAD(KDLON,KFLEV+1)
      REAL*8 ZFSUPAI(KDLON,KFLEV+1)
      REAL*8 ZFSDNAI(KDLON,KFLEV+1)
      logical initialized
      SAVE ZFSUPAD, ZFSDNAD, ZFSUPAI, ZFSDNAI ! aerosol fluxes
!rv
      save flag_aer
      data initialized/.false./
cjq-end
      if(.not.initialized) then
        flag_aer=0.
        initialized=.TRUE.
      endif
!rv
      
c
      IF (appel1er) THEN
         PRINT*, 'SW calling frequency : ', swpas
         PRINT*, "   In general, it should be 1"
         appel1er = .FALSE.
      ENDIF
C     ------------------------------------------------------------------
      IF (MOD(itapsw,swpas).EQ.0) THEN
c
      DO JK = 1 , KFLEV
      DO JL = 1, KDLON
         ZCLDSW0(JL,JK) = 0.0
         IF (bug_ozone) then
           ZOZ(JL,JK) = POZON(JL,JK)*46.6968/RG
     .               *PDP(JL,JK)*(101325.0/PPSOL(JL))
         ELSE
c        Correction MPL 100505
           ZOZ(JL,JK) = POZON(JL,JK)*RMD/RMO3*46.6968/RG*PDP(JL,JK)
         ENDIF           
      ENDDO
      ENDDO
C
C
c clear-sky:
cIM ctes ds clesphys.h  CALL SWU(PSCT,RCO2,ZCLDSW0,PPMB,PPSOL,
      CALL SWU(PSCT,ZCLDSW0,PPMB,PPSOL,
     S         PRMU0,PFRAC,PTAVE,PWV,
     S         ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
      INU = 1
      CALL SW1S(INU,
     S     PAER, flag_aer, tauae, pizae, cgae,
     S     PALBD, PALBP, PCG, ZCLD, ZCLEAR, ZCLDSW0,
     S     ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,
     S     ZFD, ZFU)
      INU = 2
      CALL SW2S(INU,
     S     PAER, flag_aer, tauae, pizae, cgae,
     S     ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, ZCLDSW0,
     S     ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,
     S     PWV, PQS,
     S     ZFDOWN, ZFUP)
      DO JK = 1 , KFLEV+1
      DO JL = 1, KDLON
         ZFSUP0(JL,JK) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
         ZFSDN0(JL,JK) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
      ENDDO
      ENDDO
      
      flag_aer=0.0
      CALL SWU(PSCT,PCLDSW,PPMB,PPSOL,
     S         PRMU0,PFRAC,PTAVE,PWV,
     S         ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
      INU = 1
      CALL SW1S(INU,
     S     PAER, flag_aer, tauae, pizae, cgae,
     S     PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,
     S     ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,
     S     ZFD, ZFU)
      INU = 2
      CALL SW2S(INU,
     S     PAER, flag_aer, tauae, pizae, cgae,
     S     ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,
     S     ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,
     S     PWV, PQS,
     S    ZFDOWN, ZFUP)

c cloudy-sky:
      
      DO JK = 1 , KFLEV+1
      DO JL = 1, KDLON
         ZFSUP(JL,JK) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
         ZFSDN(JL,JK) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
      ENDDO
      ENDDO
      
c      
      IF (ok_ade) THEN
c
c cloudy-sky + aerosol dir OB
      flag_aer=1.0
      CALL SWU(PSCT,PCLDSW,PPMB,PPSOL,
     S         PRMU0,PFRAC,PTAVE,PWV,
     S         ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
      INU = 1
      CALL SW1S(INU,
     S     PAER, flag_aer, tauae, pizae, cgae,
     S     PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,
     S     ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,
     S     ZFD, ZFU)
      INU = 2
      CALL SW2S(INU,
     S     PAER, flag_aer, tauae, pizae, cgae,
     S     ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,
     S     ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,
     S     PWV, PQS,
     S    ZFDOWN, ZFUP)
      DO JK = 1 , KFLEV+1
      DO JL = 1, KDLON
         ZFSUPAD(JL,JK) = ZFSUP(JL,JK) 
         ZFSDNAD(JL,JK) = ZFSDN(JL,JK) 
         ZFSUP(JL,JK) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
         ZFSDN(JL,JK) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
      ENDDO
      ENDDO 
      
      ENDIF ! ok_ade
      
      IF (ok_aie) THEN
         
cjq   cloudy-sky + aerosol direct + aerosol indirect
      flag_aer=1.0
      CALL SWU(PSCT,PCLDSW,PPMB,PPSOL,
     S         PRMU0,PFRAC,PTAVE,PWV,
     S         ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
      INU = 1
      CALL SW1S(INU,
     S     PAER, flag_aer, tauae, pizae, cgae,
     S     PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,
     S     ZDSIG, POMEGAA, ZOZ, ZRMU, ZSEC, PTAUA, ZUD,
     S     ZFD, ZFU)
      INU = 2
      CALL SW2S(INU,
     S     PAER, flag_aer, tauae, pizae, cgae,
     S     ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,
     S     ZDSIG, POMEGAA, ZOZ, ZRMU, ZSEC, PTAUA, ZUD,
     S     PWV, PQS,
     S    ZFDOWN, ZFUP)
      DO JK = 1 , KFLEV+1
      DO JL = 1, KDLON
         ZFSUPAI(JL,JK) = ZFSUP(JL,JK) 
         ZFSDNAI(JL,JK) = ZFSDN(JL,JK)          
         ZFSUP(JL,JK) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
         ZFSDN(JL,JK) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
      ENDDO
      ENDDO
      ENDIF ! ok_aie      
cjq -end
      
      itapsw = 0
      ENDIF
      itapsw = itapsw + 1
C
      DO k = 1, KFLEV
         kpl1 = k+1
         DO i = 1, KDLON
            PHEAT(i,k) = -(ZFSUP(i,kpl1)-ZFSUP(i,k))
     .                     -(ZFSDN(i,k)-ZFSDN(i,kpl1))
            PHEAT(i,k) = PHEAT(i,k) * RDAY*RG/RCPD / PDP(i,k)
            PHEAT0(i,k) = -(ZFSUP0(i,kpl1)-ZFSUP0(i,k))
     .                     -(ZFSDN0(i,k)-ZFSDN0(i,kpl1))
            PHEAT0(i,k) = PHEAT0(i,k) * RDAY*RG/RCPD / PDP(i,k)
         ENDDO
      ENDDO
      DO i = 1, KDLON
         PALBPLA(i) = ZFSUP(i,KFLEV+1)/(ZFSDN(i,KFLEV+1)+1.0e-20)
c
         PSOLSW(i) = ZFSDN(i,1) - ZFSUP(i,1)
         PTOPSW(i) = ZFSDN(i,KFLEV+1) - ZFSUP(i,KFLEV+1)
c
         PSOLSW0(i) = ZFSDN0(i,1) - ZFSUP0(i,1)
         PTOPSW0(i) = ZFSDN0(i,KFLEV+1) - ZFSUP0(i,KFLEV+1)
c-OB
         PSOLSWAD(i) = ZFSDNAD(i,1) - ZFSUPAD(i,1)
         PTOPSWAD(i) = ZFSDNAD(i,KFLEV+1) - ZFSUPAD(i,KFLEV+1)
c
         PSOLSWAI(i) = ZFSDNAI(i,1) - ZFSUPAI(i,1)
         PTOPSWAI(i) = ZFSDNAI(i,KFLEV+1) - ZFSUPAI(i,KFLEV+1)
c-fin 
      ENDDO
C
      RETURN
      END
c
cIM ctes ds clesphys.h   SUBROUTINE SWU (PSCT,RCO2,PCLDSW,PPMB,PPSOL,PRMU0,PFRAC,
      SUBROUTINE SWU (PSCT,PCLDSW,PPMB,PPSOL,PRMU0,PFRAC,
     S                PTAVE,PWV,PAKI,PCLD,PCLEAR,PDSIG,PFACT,
     S                PRMU,PSEC,PUD)
      use dimens_m
      use dimphy
      use clesphys
      use YOMCST
      use raddim
      use radepsi
      use radopt
      IMPLICIT none
C
C* ARGUMENTS:
C
      REAL*8 PSCT
cIM ctes ds clesphys.h   REAL*8 RCO2
      REAL*8 PCLDSW(KDLON,KFLEV)
      REAL*8 PPMB(KDLON,KFLEV+1)
      REAL*8 PPSOL(KDLON)
      REAL*8 PRMU0(KDLON)
      REAL*8 PFRAC(KDLON)
      REAL*8 PTAVE(KDLON,KFLEV)
      REAL*8 PWV(KDLON,KFLEV)
C
      REAL*8 PAKI(KDLON,2)
      REAL*8 PCLD(KDLON,KFLEV)
      REAL*8 PCLEAR(KDLON)
      REAL*8 PDSIG(KDLON,KFLEV)
      REAL*8 PFACT(KDLON)
      REAL*8 PRMU(KDLON)
      REAL*8 PSEC(KDLON)
      REAL*8 PUD(KDLON,5,KFLEV+1)
C
C* LOCAL VARIABLES:
C
      INTEGER IIND(2)
      REAL*8 ZC1J(KDLON,KFLEV+1)
      REAL*8 ZCLEAR(KDLON)
      REAL*8 ZCLOUD(KDLON)
      REAL*8 ZN175(KDLON)
      REAL*8 ZN190(KDLON)
      REAL*8 ZO175(KDLON)
      REAL*8 ZO190(KDLON)
      REAL*8 ZSIGN(KDLON)
      REAL*8 ZR(KDLON,2) 
      REAL*8 ZSIGO(KDLON)
      REAL*8 ZUD(KDLON,2)
      REAL*8 ZRTH, ZRTU, ZWH2O, ZDSCO2, ZDSH2O, ZFPPW
      INTEGER jl, jk, jkp1, jkl, jklp1, ja
C
C* Prescribed Data:
c
      REAL*8 ZPDH2O,ZPDUMG
      SAVE ZPDH2O,ZPDUMG
      REAL*8 ZPRH2O,ZPRUMG
      SAVE ZPRH2O,ZPRUMG
      REAL*8 RTDH2O,RTDUMG
      SAVE RTDH2O,RTDUMG
      REAL*8 RTH2O ,RTUMG
      SAVE RTH2O ,RTUMG
      DATA ZPDH2O,ZPDUMG / 0.8   , 0.75 /
      DATA ZPRH2O,ZPRUMG / 30000., 30000. /
      DATA RTDH2O,RTDUMG /  0.40  , 0.375 /
      DATA RTH2O ,RTUMG  /  240.  , 240.  /
C     ------------------------------------------------------------------
C
C*         1.     COMPUTES AMOUNTS OF ABSORBERS
C                 -----------------------------
C
 100  CONTINUE
C
      IIND(1)=1
      IIND(2)=2
C      
C
C*         1.1    INITIALIZES QUANTITIES
C                 ----------------------
C
 110  CONTINUE
C
      DO 111 JL = 1, KDLON
      PUD(JL,1,KFLEV+1)=0.
      PUD(JL,2,KFLEV+1)=0.
      PUD(JL,3,KFLEV+1)=0.
      PUD(JL,4,KFLEV+1)=0.
      PUD(JL,5,KFLEV+1)=0.
      PFACT(JL)= PRMU0(JL) * PFRAC(JL) * PSCT
      PRMU(JL)=SQRT(1224.* PRMU0(JL) * PRMU0(JL) + 1.) / 35.
      PSEC(JL)=1./PRMU(JL)
      ZC1J(JL,KFLEV+1)=0.
 111  CONTINUE
C
C*          1.3    AMOUNTS OF ABSORBERS
C                  --------------------
C
 130  CONTINUE
C
      DO 131 JL= 1, KDLON
      ZUD(JL,1) = 0.
      ZUD(JL,2) = 0.
      ZO175(JL) = PPSOL(JL)** (ZPDUMG+1.)
      ZO190(JL) = PPSOL(JL)** (ZPDH2O+1.)
      ZSIGO(JL) = PPSOL(JL)
      ZCLEAR(JL)=1.
      ZCLOUD(JL)=0.
 131  CONTINUE
C
      DO 133 JK = 1 , KFLEV
      JKP1 = JK + 1
      JKL = KFLEV+1 - JK
      JKLP1 = JKL+1
      DO 132 JL = 1, KDLON
      ZRTH=(RTH2O/PTAVE(JL,JK))**RTDH2O
      ZRTU=(RTUMG/PTAVE(JL,JK))**RTDUMG
      ZWH2O = MAX (PWV(JL,JK) , ZEPSCQ )
      ZSIGN(JL) = 100. * PPMB(JL,JKP1)
      PDSIG(JL,JK) = (ZSIGO(JL) - ZSIGN(JL))/PPSOL(JL)
      ZN175(JL) = ZSIGN(JL) ** (ZPDUMG+1.)
      ZN190(JL) = ZSIGN(JL) ** (ZPDH2O+1.)
      ZDSCO2 = ZO175(JL) - ZN175(JL)
      ZDSH2O = ZO190(JL) - ZN190(JL)
      PUD(JL,1,JK) = 1./( 10.* RG * (ZPDH2O+1.) )/(ZPRH2O**ZPDH2O)
     .             * ZDSH2O * ZWH2O  * ZRTH
      PUD(JL,2,JK) = 1./( 10.* RG * (ZPDUMG+1.) )/(ZPRUMG**ZPDUMG)
     .             * ZDSCO2 * RCO2 * ZRTU
      ZFPPW=1.6078*ZWH2O/(1.+0.608*ZWH2O)
      PUD(JL,4,JK)=PUD(JL,1,JK)*ZFPPW
      PUD(JL,5,JK)=PUD(JL,1,JK)*(1.-ZFPPW)
      ZUD(JL,1) = ZUD(JL,1) + PUD(JL,1,JK)
      ZUD(JL,2) = ZUD(JL,2) + PUD(JL,2,JK)
      ZSIGO(JL) = ZSIGN(JL)
      ZO175(JL) = ZN175(JL)
      ZO190(JL) = ZN190(JL)
C      
      IF (NOVLP.EQ.1) THEN
         ZCLEAR(JL)=ZCLEAR(JL)
     S               *(1.-MAX(PCLDSW(JL,JKL),ZCLOUD(JL)))
     S               /(1.-MIN(ZCLOUD(JL),1.-ZEPSEC))
         ZC1J(JL,JKL)= 1.0 - ZCLEAR(JL)
         ZCLOUD(JL) = PCLDSW(JL,JKL)
      ELSE IF (NOVLP.EQ.2) THEN
         ZCLOUD(JL) = MAX(PCLDSW(JL,JKL),ZCLOUD(JL))
         ZC1J(JL,JKL) = ZCLOUD(JL)
      ELSE IF (NOVLP.EQ.3) THEN
         ZCLEAR(JL) = ZCLEAR(JL)*(1.-PCLDSW(JL,JKL))
         ZCLOUD(JL) = 1.0 - ZCLEAR(JL)
         ZC1J(JL,JKL) = ZCLOUD(JL)
      END IF
 132  CONTINUE
 133  CONTINUE
      DO 134 JL=1, KDLON
      PCLEAR(JL)=1.-ZC1J(JL,1)
 134  CONTINUE
      DO 136 JK=1,KFLEV
      DO 135 JL=1, KDLON
      IF (PCLEAR(JL).LT.1.) THEN
         PCLD(JL,JK)=PCLDSW(JL,JK)/(1.-PCLEAR(JL))
      ELSE
         PCLD(JL,JK)=0.
      END IF
 135  CONTINUE
 136  CONTINUE           
C      
C
C*         1.4    COMPUTES CLEAR-SKY GREY ABSORPTION COEFFICIENTS
C                 -----------------------------------------------
C
 140  CONTINUE
C
      DO 142 JA = 1,2
      DO 141 JL = 1, KDLON
      ZUD(JL,JA) = ZUD(JL,JA) * PSEC(JL)
 141  CONTINUE
 142  CONTINUE
C
      CALL SWTT1(2, 2, IIND, ZUD, ZR)
C
      DO 144 JA = 1,2
      DO 143 JL = 1, KDLON
      PAKI(JL,JA) = -LOG( ZR(JL,JA) ) / ZUD(JL,JA)
 143  CONTINUE
 144  CONTINUE
C
C
C     ------------------------------------------------------------------
C
      RETURN
      END
      SUBROUTINE SW1S ( KNU
     S  ,  PAER  , flag_aer, tauae, pizae, cgae
     S  ,  PALBD , PALBP, PCG  , PCLD , PCLEAR, PCLDSW
     S  ,  PDSIG , POMEGA, POZ  , PRMU , PSEC , PTAU  , PUD  
     S  ,  PFD   , PFU)
      use dimens_m
      use dimphy
      use raddim
      IMPLICIT none
C
C     ------------------------------------------------------------------
C     PURPOSE.
C     --------
C
C          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
C     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).
C
C     METHOD.
C     -------
C
C          1. COMPUTES UPWARD AND DOWNWARD FLUXES CORRESPONDING TO
C     CONTINUUM SCATTERING
C          2. MULTIPLY BY OZONE TRANSMISSION FUNCTION
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
C        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C        94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO
C     ------------------------------------------------------------------
C
C* ARGUMENTS:
C
      INTEGER KNU
c-OB
      real*8 flag_aer
      real*8 tauae(kdlon,kflev,2)
      real*8 pizae(kdlon,kflev,2)
      real*8 cgae(kdlon,kflev,2)
      REAL*8 PAER(KDLON,KFLEV,5)
      REAL*8 PALBD(KDLON,2)
      REAL*8 PALBP(KDLON,2)
      REAL*8 PCG(KDLON,2,KFLEV)  
      REAL*8 PCLD(KDLON,KFLEV)
      REAL*8 PCLDSW(KDLON,KFLEV)
      REAL*8 PCLEAR(KDLON)
      REAL*8 PDSIG(KDLON,KFLEV)
      REAL*8 POMEGA(KDLON,2,KFLEV)
      REAL*8 POZ(KDLON,KFLEV)
      REAL*8 PRMU(KDLON)
      REAL*8 PSEC(KDLON)
      REAL*8 PTAU(KDLON,2,KFLEV)
      REAL*8 PUD(KDLON,5,KFLEV+1)
C
      REAL*8 PFD(KDLON,KFLEV+1)
      REAL*8 PFU(KDLON,KFLEV+1)
C
C* LOCAL VARIABLES:
C
      INTEGER IIND(4)
C      
      REAL*8 ZCGAZ(KDLON,KFLEV) 
      REAL*8 ZDIFF(KDLON)
      REAL*8 ZDIRF(KDLON)        
      REAL*8 ZPIZAZ(KDLON,KFLEV)
      REAL*8 ZRAYL(KDLON)
      REAL*8 ZRAY1(KDLON,KFLEV+1)
      REAL*8 ZRAY2(KDLON,KFLEV+1)
      REAL*8 ZREFZ(KDLON,2,KFLEV+1)
      REAL*8 ZRJ(KDLON,6,KFLEV+1)
      REAL*8 ZRJ0(KDLON,6,KFLEV+1)
      REAL*8 ZRK(KDLON,6,KFLEV+1)
      REAL*8 ZRK0(KDLON,6,KFLEV+1)
      REAL*8 ZRMUE(KDLON,KFLEV+1)
      REAL*8 ZRMU0(KDLON,KFLEV+1)
      REAL*8 ZR(KDLON,4)
      REAL*8 ZTAUAZ(KDLON,KFLEV)
      REAL*8 ZTRA1(KDLON,KFLEV+1)
      REAL*8 ZTRA2(KDLON,KFLEV+1)
      REAL*8 ZW(KDLON,4)
C
      INTEGER jl, jk, k, jaj, ikm1, ikl
c
c Prescribed Data:
c
      REAL*8 RSUN(2)
      SAVE RSUN
      REAL*8 RRAY(2,6)
      SAVE RRAY
      DATA RSUN(1) / 0.441676 /
      DATA RSUN(2) / 0.558324 /
      DATA (RRAY(1,K),K=1,6) /
     S .428937E-01, .890743E+00,-.288555E+01,
     S .522744E+01,-.469173E+01, .161645E+01/
      DATA (RRAY(2,K),K=1,6) /
     S .697200E-02, .173297E-01,-.850903E-01,
     S .248261E+00,-.302031E+00, .129662E+00/
C     ------------------------------------------------------------------
C
C*         1.     FIRST SPECTRAL INTERVAL (0.25-0.68 MICRON)
C                 ----------------------- ------------------
C
 100  CONTINUE
C
C
C*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
C                 -----------------------------------------
C
 110  CONTINUE
C
      DO 111 JL = 1, KDLON
      ZRAYL(JL) =  RRAY(KNU,1) + PRMU(JL) * (RRAY(KNU,2) + PRMU(JL)
     S          * (RRAY(KNU,3) + PRMU(JL) * (RRAY(KNU,4) + PRMU(JL)
     S          * (RRAY(KNU,5) + PRMU(JL) *  RRAY(KNU,6)       ))))
 111  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         2.    CONTINUUM SCATTERING CALCULATIONS
C                ---------------------------------
C
 200  CONTINUE
C
C*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
C                --------------------------------
C  
 210  CONTINUE
C
      CALL SWCLR ( KNU
     S  , PAER   , flag_aer, tauae, pizae, cgae
     S  , PALBP  , PDSIG , ZRAYL, PSEC
     S  , ZCGAZ  , ZPIZAZ, ZRAY1 , ZRAY2, ZREFZ, ZRJ0
     S  , ZRK0   , ZRMU0 , ZTAUAZ, ZTRA1, ZTRA2)
C
C
C*         2.2   CLOUDY FRACTION OF THE COLUMN
C                -----------------------------
C
 220  CONTINUE
C
      CALL SWR ( KNU
     S  , PALBD ,PCG   ,PCLD  ,PDSIG ,POMEGA,ZRAYL
     S  , PSEC  ,PTAU
     S  , ZCGAZ ,ZPIZAZ,ZRAY1 ,ZRAY2 ,ZREFZ ,ZRJ  ,ZRK,ZRMUE
     S  , ZTAUAZ,ZTRA1 ,ZTRA2)
C
C
C     ------------------------------------------------------------------
C
C*         3.    OZONE ABSORPTION
C                ----------------
C
 300  CONTINUE
C
      IIND(1)=1
      IIND(2)=3
      IIND(3)=1
      IIND(4)=3
C      
C
C*         3.1   DOWNWARD FLUXES
C                ---------------
C
 310  CONTINUE
C
      JAJ = 2
C
      DO 311 JL = 1, KDLON
      ZW(JL,1)=0.
      ZW(JL,2)=0.
      ZW(JL,3)=0.
      ZW(JL,4)=0.
      PFD(JL,KFLEV+1)=((1.-PCLEAR(JL))*ZRJ(JL,JAJ,KFLEV+1)
     S     + PCLEAR(JL) *ZRJ0(JL,JAJ,KFLEV+1)) * RSUN(KNU)
 311  CONTINUE
      DO 314 JK = 1 , KFLEV
      IKL = KFLEV+1-JK
      DO 312 JL = 1, KDLON
      ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKL)/ZRMUE(JL,IKL)
      ZW(JL,2)=ZW(JL,2)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
      ZW(JL,3)=ZW(JL,3)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
      ZW(JL,4)=ZW(JL,4)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
 312  CONTINUE
C
      CALL SWTT1(KNU, 4, IIND, ZW, ZR)
C
      DO 313 JL = 1, KDLON
      ZDIFF(JL) = ZR(JL,1)*ZR(JL,2)*ZRJ(JL,JAJ,IKL)
      ZDIRF(JL) = ZR(JL,3)*ZR(JL,4)*ZRJ0(JL,JAJ,IKL)
      PFD(JL,IKL) = ((1.-PCLEAR(JL)) * ZDIFF(JL)
     S                  +PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
 313  CONTINUE
 314  CONTINUE
C
C
C*         3.2   UPWARD FLUXES
C                -------------
C
 320  CONTINUE
C
      DO 325 JL = 1, KDLON
      PFU(JL,1) = ((1.-PCLEAR(JL))*ZDIFF(JL)*PALBD(JL,KNU)
     S               + PCLEAR(JL) *ZDIRF(JL)*PALBP(JL,KNU))
     S          * RSUN(KNU)
 325  CONTINUE
C
      DO 328 JK = 2 , KFLEV+1
      IKM1=JK-1
      DO 326 JL = 1, KDLON
      ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKM1)*1.66
      ZW(JL,2)=ZW(JL,2)+POZ(JL,  IKM1)*1.66
      ZW(JL,3)=ZW(JL,3)+PUD(JL,1,IKM1)*1.66
      ZW(JL,4)=ZW(JL,4)+POZ(JL,  IKM1)*1.66
 326  CONTINUE
C
      CALL SWTT1(KNU, 4, IIND, ZW, ZR)
C
      DO 327 JL = 1, KDLON
      ZDIFF(JL) = ZR(JL,1)*ZR(JL,2)*ZRK(JL,JAJ,JK)
      ZDIRF(JL) = ZR(JL,3)*ZR(JL,4)*ZRK0(JL,JAJ,JK)
      PFU(JL,JK) = ((1.-PCLEAR(JL)) * ZDIFF(JL)
     S                 +PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
 327  CONTINUE
 328  CONTINUE
C
C     ------------------------------------------------------------------
C
      RETURN
      END
      SUBROUTINE SW2S ( KNU
     S  ,  PAER  , flag_aer, tauae, pizae, cgae
     S  ,  PAKI, PALBD, PALBP, PCG   , PCLD, PCLEAR, PCLDSW
     S  ,  PDSIG ,POMEGA,POZ , PRMU , PSEC  , PTAU
     S  ,  PUD   ,PWV , PQS
     S  ,  PFDOWN,PFUP                                            )
      use dimens_m
      use dimphy
      use raddim
      use radepsi
      IMPLICIT none
C
C     ------------------------------------------------------------------
C     PURPOSE.
C     --------
C
C          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN THE
C     SECOND SPECTRAL INTERVAL FOLLOWING FOUQUART AND BONNEL (1980).
C
C     METHOD.
C     -------
C
C          1. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING TO
C     CONTINUUM SCATTERING
C          2. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING FOR
C     A GREY MOLECULAR ABSORPTION
C          3. LAPLACE TRANSFORM ON THE PREVIOUS TO GET EFFECTIVE AMOUNTS
C     OF ABSORBERS
C          4. APPLY H2O AND U.M.G. TRANSMISSION FUNCTIONS
C          5. MULTIPLY BY OZONE TRANSMISSION FUNCTION
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
C        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C        94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO
C     ------------------------------------------------------------------
C* ARGUMENTS:
C
      INTEGER KNU
c-OB
      real*8 flag_aer
      real*8 tauae(kdlon,kflev,2)
      real*8 pizae(kdlon,kflev,2)
      real*8 cgae(kdlon,kflev,2)
      REAL*8 PAER(KDLON,KFLEV,5)
      REAL*8 PAKI(KDLON,2)
      REAL*8 PALBD(KDLON,2)
      REAL*8 PALBP(KDLON,2)
      REAL*8 PCG(KDLON,2,KFLEV)
      REAL*8 PCLD(KDLON,KFLEV)
      REAL*8 PCLDSW(KDLON,KFLEV)
      REAL*8 PCLEAR(KDLON)
      REAL*8 PDSIG(KDLON,KFLEV)
      REAL*8 POMEGA(KDLON,2,KFLEV)
      REAL*8 POZ(KDLON,KFLEV)
      REAL*8 PQS(KDLON,KFLEV)
      REAL*8 PRMU(KDLON)
      REAL*8 PSEC(KDLON)
      REAL*8 PTAU(KDLON,2,KFLEV)
      REAL*8 PUD(KDLON,5,KFLEV+1)
      REAL*8 PWV(KDLON,KFLEV)
C
      REAL*8 PFDOWN(KDLON,KFLEV+1)
      REAL*8 PFUP(KDLON,KFLEV+1)
C
C* LOCAL VARIABLES:
C
      INTEGER IIND2(2), IIND3(3)
      REAL*8 ZCGAZ(KDLON,KFLEV)
      REAL*8 ZFD(KDLON,KFLEV+1)
      REAL*8 ZFU(KDLON,KFLEV+1) 
      REAL*8 ZG(KDLON)
      REAL*8 ZGG(KDLON)
      REAL*8 ZPIZAZ(KDLON,KFLEV)
      REAL*8 ZRAYL(KDLON)
      REAL*8 ZRAY1(KDLON,KFLEV+1)
      REAL*8 ZRAY2(KDLON,KFLEV+1)
      REAL*8 ZREF(KDLON)
      REAL*8 ZREFZ(KDLON,2,KFLEV+1)
      REAL*8 ZRE1(KDLON)
      REAL*8 ZRE2(KDLON)
      REAL*8 ZRJ(KDLON,6,KFLEV+1)
      REAL*8 ZRJ0(KDLON,6,KFLEV+1)
      REAL*8 ZRK(KDLON,6,KFLEV+1)
      REAL*8 ZRK0(KDLON,6,KFLEV+1)
      REAL*8 ZRL(KDLON,8)
      REAL*8 ZRMUE(KDLON,KFLEV+1)
      REAL*8 ZRMU0(KDLON,KFLEV+1)
      REAL*8 ZRMUZ(KDLON)
      REAL*8 ZRNEB(KDLON)
      REAL*8 ZRUEF(KDLON,8)
      REAL*8 ZR1(KDLON) 
      REAL*8 ZR2(KDLON,2)
      REAL*8 ZR3(KDLON,3)
      REAL*8 ZR4(KDLON)
      REAL*8 ZR21(KDLON)
      REAL*8 ZR22(KDLON)
      REAL*8 ZS(KDLON)
      REAL*8 ZTAUAZ(KDLON,KFLEV)
      REAL*8 ZTO1(KDLON)
      REAL*8 ZTR(KDLON,2,KFLEV+1)
      REAL*8 ZTRA1(KDLON,KFLEV+1)
      REAL*8 ZTRA2(KDLON,KFLEV+1)
      REAL*8 ZTR1(KDLON)
      REAL*8 ZTR2(KDLON)
      REAL*8 ZW(KDLON)   
      REAL*8 ZW1(KDLON)
      REAL*8 ZW2(KDLON,2)
      REAL*8 ZW3(KDLON,3)
      REAL*8 ZW4(KDLON)
      REAL*8 ZW5(KDLON)
C
      INTEGER jl, jk, k, jaj, ikm1, ikl, jn, jabs, jkm1
      INTEGER jref, jkl, jklp1, jajp, jkki, jkkp4, jn2j, iabs
      REAL*8 ZRMUM1, ZWH2O, ZCNEB, ZAA, ZBB, ZRKI, ZRE11
C
C* Prescribed Data:
C
      REAL*8 RSUN(2)
      SAVE RSUN
      REAL*8 RRAY(2,6)
      SAVE RRAY
      DATA RSUN(1) / 0.441676 /
      DATA RSUN(2) / 0.558324 /
      DATA (RRAY(1,K),K=1,6) /
     S .428937E-01, .890743E+00,-.288555E+01,
     S .522744E+01,-.469173E+01, .161645E+01/
      DATA (RRAY(2,K),K=1,6) /
     S .697200E-02, .173297E-01,-.850903E-01,
     S .248261E+00,-.302031E+00, .129662E+00/
C
C     ------------------------------------------------------------------
C
C*         1.     SECOND SPECTRAL INTERVAL (0.68-4.00 MICRON)
C                 -------------------------------------------
C
 100  CONTINUE
C
C
C*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
C                 -----------------------------------------
C
 110  CONTINUE
C
      DO 111 JL = 1, KDLON
      ZRMUM1 = 1. - PRMU(JL)
      ZRAYL(JL) =  RRAY(KNU,1) + ZRMUM1   * (RRAY(KNU,2) + ZRMUM1
     S          * (RRAY(KNU,3) + ZRMUM1   * (RRAY(KNU,4) + ZRMUM1
     S          * (RRAY(KNU,5) + ZRMUM1   *  RRAY(KNU,6)     ))))
 111  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         2.    CONTINUUM SCATTERING CALCULATIONS
C                ---------------------------------
C
 200  CONTINUE
C
C*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
C                --------------------------------
C  
 210  CONTINUE
C
      CALL SWCLR ( KNU
     S  , PAER   , flag_aer, tauae, pizae, cgae
     S  , PALBP  , PDSIG , ZRAYL, PSEC 
     S  , ZCGAZ  , ZPIZAZ, ZRAY1 , ZRAY2, ZREFZ, ZRJ0
     S  , ZRK0   , ZRMU0 , ZTAUAZ, ZTRA1, ZTRA2)
C
C
C*         2.2   CLOUDY FRACTION OF THE COLUMN
C                -----------------------------
C
 220  CONTINUE
C
      CALL SWR ( KNU
     S  , PALBD , PCG   , PCLD , PDSIG, POMEGA, ZRAYL
     S  , PSEC  , PTAU
     S  , ZCGAZ , ZPIZAZ, ZRAY1, ZRAY2, ZREFZ , ZRJ  , ZRK, ZRMUE
     S  , ZTAUAZ, ZTRA1 , ZTRA2)
C
C
C     ------------------------------------------------------------------
C
C*         3.    SCATTERING CALCULATIONS WITH GREY MOLECULAR ABSORPTION
C                ------------------------------------------------------
C
 300  CONTINUE
C
      JN = 2
C
      DO 361 JABS=1,2
C
C
C*         3.1  SURFACE CONDITIONS
C               ------------------
C
 310  CONTINUE
C
      DO 311 JL = 1, KDLON
      ZREFZ(JL,2,1) = PALBD(JL,KNU)
      ZREFZ(JL,1,1) = PALBD(JL,KNU)
 311  CONTINUE
C
C
C*         3.2  INTRODUCING CLOUD EFFECTS
C               -------------------------
C
 320  CONTINUE
C
      DO 324 JK = 2 , KFLEV+1
      JKM1 = JK - 1
      IKL=KFLEV+1-JKM1
      DO 322 JL = 1, KDLON
      ZRNEB(JL) = PCLD(JL,JKM1)
      IF (JABS.EQ.1 .AND. ZRNEB(JL).GT.2.*ZEELOG) THEN
         ZWH2O=MAX(PWV(JL,JKM1),ZEELOG)
         ZCNEB=MAX(ZEELOG,MIN(ZRNEB(JL),1.-ZEELOG))
         ZBB=PUD(JL,JABS,JKM1)*PQS(JL,JKM1)/ZWH2O
         ZAA=MAX((PUD(JL,JABS,JKM1)-ZCNEB*ZBB)/(1.-ZCNEB),ZEELOG)
      ELSE
         ZAA=PUD(JL,JABS,JKM1)
         ZBB=ZAA
      END IF
      ZRKI = PAKI(JL,JABS)
      ZS(JL) = EXP(-ZRKI * ZAA * 1.66)
      ZG(JL) = EXP(-ZRKI * ZAA / ZRMUE(JL,JK))
      ZTR1(JL) = 0.
      ZRE1(JL) = 0.
      ZTR2(JL) = 0.
      ZRE2(JL) = 0.
C
      ZW(JL)= POMEGA(JL,KNU,JKM1)
      ZTO1(JL) = PTAU(JL,KNU,JKM1) / ZW(JL)
     S               + ZTAUAZ(JL,JKM1) / ZPIZAZ(JL,JKM1)
     S               + ZBB * ZRKI

      ZR21(JL) = PTAU(JL,KNU,JKM1) + ZTAUAZ(JL,JKM1)
      ZR22(JL) = PTAU(JL,KNU,JKM1) / ZR21(JL)
      ZGG(JL) = ZR22(JL) * PCG(JL,KNU,JKM1)
     S              + (1. - ZR22(JL)) * ZCGAZ(JL,JKM1)
      ZW(JL) = ZR21(JL) / ZTO1(JL)
      ZREF(JL) = ZREFZ(JL,1,JKM1)
      ZRMUZ(JL) = ZRMUE(JL,JK)
 322  CONTINUE
C
      CALL SWDE(ZGG, ZREF, ZRMUZ, ZTO1, ZW,
     S          ZRE1, ZRE2, ZTR1, ZTR2)
C
      DO 323 JL = 1, KDLON
C
      ZREFZ(JL,2,JK) = (1.-ZRNEB(JL)) * (ZRAY1(JL,JKM1)
     S               + ZREFZ(JL,2,JKM1) * ZTRA1(JL,JKM1)
     S               * ZTRA2(JL,JKM1) ) * ZG(JL) * ZS(JL)
     S               + ZRNEB(JL) * ZRE1(JL)
C
      ZTR(JL,2,JKM1)=ZRNEB(JL)*ZTR1(JL)
     S              + (ZTRA1(JL,JKM1)) * ZG(JL) * (1.-ZRNEB(JL))
C
      ZREFZ(JL,1,JK)=(1.-ZRNEB(JL))*(ZRAY1(JL,JKM1)
     S                  +ZREFZ(JL,1,JKM1)*ZTRA1(JL,JKM1)*ZTRA2(JL,JKM1)
     S             /(1.-ZRAY2(JL,JKM1)*ZREFZ(JL,1,JKM1)))*ZG(JL)*ZS(JL)
     S             + ZRNEB(JL) * ZRE2(JL)
C
      ZTR(JL,1,JKM1)= ZRNEB(JL) * ZTR2(JL)
     S              + (ZTRA1(JL,JKM1)/(1.-ZRAY2(JL,JKM1)
     S              * ZREFZ(JL,1,JKM1)))
     S              * ZG(JL) * (1. -ZRNEB(JL))
C
 323  CONTINUE
 324  CONTINUE
C
C*         3.3  REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
C               -------------------------------------------------
C
 330  CONTINUE
C
      DO 351 JREF=1,2
C
      JN = JN + 1
C
      DO 331 JL = 1, KDLON
      ZRJ(JL,JN,KFLEV+1) = 1.
      ZRK(JL,JN,KFLEV+1) = ZREFZ(JL,JREF,KFLEV+1)
 331  CONTINUE
C
      DO 333 JK = 1 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 332 JL = 1, KDLON
      ZRE11 = ZRJ(JL,JN,JKLP1) * ZTR(JL,JREF,JKL)
      ZRJ(JL,JN,JKL) = ZRE11
      ZRK(JL,JN,JKL) = ZRE11 * ZREFZ(JL,JREF,JKL)
 332  CONTINUE
 333  CONTINUE
 351  CONTINUE
 361  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         4.    INVERT GREY AND CONTINUUM FLUXES
C                --------------------------------
C
 400  CONTINUE
C
C
C*         4.1   UPWARD (ZRK) AND DOWNWARD (ZRJ) PSEUDO-FLUXES
C                ---------------------------------------------
C
 410  CONTINUE
C
      DO 414 JK = 1 , KFLEV+1
      DO 413 JAJ = 1 , 5 , 2
      JAJP = JAJ + 1
      DO 412 JL = 1, KDLON
      ZRJ(JL,JAJ,JK)=        ZRJ(JL,JAJ,JK) - ZRJ(JL,JAJP,JK)
      ZRK(JL,JAJ,JK)=        ZRK(JL,JAJ,JK) - ZRK(JL,JAJP,JK)
      ZRJ(JL,JAJ,JK)= MAX( ZRJ(JL,JAJ,JK) , ZEELOG )
      ZRK(JL,JAJ,JK)= MAX( ZRK(JL,JAJ,JK) , ZEELOG )
 412  CONTINUE
 413  CONTINUE
 414  CONTINUE
C
      DO 417 JK = 1 , KFLEV+1
      DO 416 JAJ = 2 , 6 , 2
      DO 415 JL = 1, KDLON
      ZRJ(JL,JAJ,JK)= MAX( ZRJ(JL,JAJ,JK) , ZEELOG )
      ZRK(JL,JAJ,JK)= MAX( ZRK(JL,JAJ,JK) , ZEELOG )
 415  CONTINUE
 416  CONTINUE
 417  CONTINUE
C
C*         4.2    EFFECTIVE ABSORBER AMOUNTS BY INVERSE LAPLACE
C                 ---------------------------------------------
C
 420  CONTINUE
C
      DO 437 JK = 1 , KFLEV+1
      JKKI = 1
      DO 425 JAJ = 1 , 2
      IIND2(1)=JAJ
      IIND2(2)=JAJ
      DO 424 JN = 1 , 2
      JN2J = JN + 2 * JAJ
      JKKP4 = JKKI + 4
C
C*         4.2.1  EFFECTIVE ABSORBER AMOUNTS
C                 --------------------------
C
 4210 CONTINUE
C
      DO 4211 JL = 1, KDLON
      ZW2(JL,1) = LOG( ZRJ(JL,JN,JK) / ZRJ(JL,JN2J,JK))
     S                               / PAKI(JL,JAJ)
      ZW2(JL,2) = LOG( ZRK(JL,JN,JK) / ZRK(JL,JN2J,JK))
     S                               / PAKI(JL,JAJ)
 4211 CONTINUE
C
C*         4.2.2  TRANSMISSION FUNCTION
C                 ---------------------
C
 4220 CONTINUE
C
      CALL SWTT1(KNU, 2, IIND2, ZW2, ZR2)
C
      DO 4221 JL = 1, KDLON
      ZRL(JL,JKKI) = ZR2(JL,1)
      ZRUEF(JL,JKKI) = ZW2(JL,1)
      ZRL(JL,JKKP4) = ZR2(JL,2)
      ZRUEF(JL,JKKP4) = ZW2(JL,2)
 4221 CONTINUE
C
      JKKI=JKKI+1
 424  CONTINUE
 425  CONTINUE
C
C*         4.3    UPWARD AND DOWNWARD FLUXES WITH H2O AND UMG ABSORPTION
C                 ------------------------------------------------------
C
 430  CONTINUE
C
      DO 431 JL = 1, KDLON
      PFDOWN(JL,JK) = ZRJ(JL,1,JK) * ZRL(JL,1) * ZRL(JL,3)
     S              + ZRJ(JL,2,JK) * ZRL(JL,2) * ZRL(JL,4)
      PFUP(JL,JK)   = ZRK(JL,1,JK) * ZRL(JL,5) * ZRL(JL,7)
     S              + ZRK(JL,2,JK) * ZRL(JL,6) * ZRL(JL,8)
 431  CONTINUE
 437  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         5.    MOLECULAR ABSORPTION ON CLEAR-SKY FLUXES
C                ----------------------------------------
C
 500  CONTINUE
C
C
C*         5.1   DOWNWARD FLUXES
C                ---------------
C
 510  CONTINUE
C
      JAJ = 2
      IIND3(1)=1
      IIND3(2)=2
      IIND3(3)=3
C      
      DO 511 JL = 1, KDLON
      ZW3(JL,1)=0.
      ZW3(JL,2)=0.
      ZW3(JL,3)=0.
      ZW4(JL)  =0.
      ZW5(JL)  =0.
      ZR4(JL)  =1.
      ZFD(JL,KFLEV+1)= ZRJ0(JL,JAJ,KFLEV+1)
 511  CONTINUE
      DO 514 JK = 1 , KFLEV
      IKL = KFLEV+1-JK
      DO 512 JL = 1, KDLON
      ZW3(JL,1)=ZW3(JL,1)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
      ZW3(JL,2)=ZW3(JL,2)+PUD(JL,2,IKL)/ZRMU0(JL,IKL)
      ZW3(JL,3)=ZW3(JL,3)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
      ZW4(JL)  =ZW4(JL)  +PUD(JL,4,IKL)/ZRMU0(JL,IKL)
      ZW5(JL)  =ZW5(JL)  +PUD(JL,5,IKL)/ZRMU0(JL,IKL)
 512  CONTINUE
C
      CALL SWTT1(KNU, 3, IIND3, ZW3, ZR3)
C
      DO 513 JL = 1, KDLON
C     ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
      ZFD(JL,IKL) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)
     S            * ZRJ0(JL,JAJ,IKL)
 513  CONTINUE
 514  CONTINUE
C
C
C*         5.2   UPWARD FLUXES
C                -------------
C
 520  CONTINUE
C
      DO 525 JL = 1, KDLON
      ZFU(JL,1) = ZFD(JL,1)*PALBP(JL,KNU)
 525  CONTINUE
C
      DO 528 JK = 2 , KFLEV+1
      IKM1=JK-1
      DO 526 JL = 1, KDLON
      ZW3(JL,1)=ZW3(JL,1)+PUD(JL,1,IKM1)*1.66
      ZW3(JL,2)=ZW3(JL,2)+PUD(JL,2,IKM1)*1.66
      ZW3(JL,3)=ZW3(JL,3)+POZ(JL,  IKM1)*1.66
      ZW4(JL)  =ZW4(JL)  +PUD(JL,4,IKM1)*1.66
      ZW5(JL)  =ZW5(JL)  +PUD(JL,5,IKM1)*1.66
 526  CONTINUE
C
      CALL SWTT1(KNU, 3, IIND3, ZW3, ZR3)
C
      DO 527 JL = 1, KDLON
C     ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
      ZFU(JL,JK) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)
     S           * ZRK0(JL,JAJ,JK)
 527  CONTINUE
 528  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         6.     INTRODUCTION OF OZONE AND H2O CONTINUUM ABSORPTION
C                 --------------------------------------------------
C
 600  CONTINUE
      IABS=3
C
C*         6.1    DOWNWARD FLUXES
C                 ---------------
C
 610  CONTINUE
      DO 611 JL = 1, KDLON
      ZW1(JL)=0.
      ZW4(JL)=0.
      ZW5(JL)=0.
      ZR1(JL)=0.
      PFDOWN(JL,KFLEV+1) = ((1.-PCLEAR(JL))*PFDOWN(JL,KFLEV+1)
     S                   + PCLEAR(JL) * ZFD(JL,KFLEV+1)) * RSUN(KNU)
 611  CONTINUE
C
      DO 614 JK = 1 , KFLEV
      IKL=KFLEV+1-JK
      DO 612 JL = 1, KDLON
      ZW1(JL) = ZW1(JL)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
      ZW4(JL) = ZW4(JL)+PUD(JL,4,IKL)/ZRMUE(JL,IKL)
      ZW5(JL) = ZW5(JL)+PUD(JL,5,IKL)/ZRMUE(JL,IKL)
C     ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
 612  CONTINUE
C
      CALL SWTT(KNU, IABS, ZW1, ZR1)
C
      DO 613 JL = 1, KDLON
      PFDOWN(JL,IKL) = ((1.-PCLEAR(JL))*ZR1(JL)*ZR4(JL)*PFDOWN(JL,IKL)
     S                     +PCLEAR(JL)*ZFD(JL,IKL)) * RSUN(KNU)
 613  CONTINUE
 614  CONTINUE
C
C
C*         6.2    UPWARD FLUXES
C                 -------------
C
 620  CONTINUE
      DO 621 JL = 1, KDLON
      PFUP(JL,1) = ((1.-PCLEAR(JL))*ZR1(JL)*ZR4(JL) * PFUP(JL,1)
     S                 +PCLEAR(JL)*ZFU(JL,1)) * RSUN(KNU)
 621  CONTINUE
C
      DO 624 JK = 2 , KFLEV+1
      IKM1=JK-1
      DO 622 JL = 1, KDLON
      ZW1(JL) = ZW1(JL)+POZ(JL  ,IKM1)*1.66
      ZW4(JL) = ZW4(JL)+PUD(JL,4,IKM1)*1.66
      ZW5(JL) = ZW5(JL)+PUD(JL,5,IKM1)*1.66
C     ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
 622  CONTINUE
C
      CALL SWTT(KNU, IABS, ZW1, ZR1)
C
      DO 623 JL = 1, KDLON
      PFUP(JL,JK) = ((1.-PCLEAR(JL))*ZR1(JL)*ZR4(JL) * PFUP(JL,JK)
     S                 +PCLEAR(JL)*ZFU(JL,JK)) * RSUN(KNU)
 623  CONTINUE
 624  CONTINUE
C
C     ------------------------------------------------------------------
C
      RETURN
      END
      SUBROUTINE SWCLR  ( KNU
     S  , PAER  , flag_aer, tauae, pizae, cgae
     S  , PALBP , PDSIG , PRAYL , PSEC
     S  , PCGAZ , PPIZAZ, PRAY1 , PRAY2 , PREFZ , PRJ  
     S  , PRK   , PRMU0 , PTAUAZ, PTRA1 , PTRA2                   )
      use dimens_m
      use dimphy
      use raddim
      use radepsi
      use radopt
      IMPLICIT none
C
C     ------------------------------------------------------------------
C     PURPOSE.
C     --------
C           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
C     CLEAR-SKY COLUMN
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
C        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 94-11-15
C     ------------------------------------------------------------------
C* ARGUMENTS:
C
      INTEGER KNU
c-OB
      real*8 flag_aer
      real*8 tauae(kdlon,kflev,2)
      real*8 pizae(kdlon,kflev,2)
      real*8 cgae(kdlon,kflev,2)
      REAL*8 PAER(KDLON,KFLEV,5)
      REAL*8 PALBP(KDLON,2)
      REAL*8 PDSIG(KDLON,KFLEV)
      REAL*8 PRAYL(KDLON)
      REAL*8 PSEC(KDLON)
C
      REAL*8 PCGAZ(KDLON,KFLEV)     
      REAL*8 PPIZAZ(KDLON,KFLEV)
      REAL*8 PRAY1(KDLON,KFLEV+1)
      REAL*8 PRAY2(KDLON,KFLEV+1)
      REAL*8 PREFZ(KDLON,2,KFLEV+1)
      REAL*8 PRJ(KDLON,6,KFLEV+1)
      REAL*8 PRK(KDLON,6,KFLEV+1)
      REAL*8 PRMU0(KDLON,KFLEV+1)
      REAL*8 PTAUAZ(KDLON,KFLEV)
      REAL*8 PTRA1(KDLON,KFLEV+1)
      REAL*8 PTRA2(KDLON,KFLEV+1)
C
C* LOCAL VARIABLES:
C
      REAL*8 ZC0I(KDLON,KFLEV+1)       
      REAL*8 ZCLE0(KDLON,KFLEV)
      REAL*8 ZCLEAR(KDLON)
      REAL*8 ZR21(KDLON)
      REAL*8 ZR23(KDLON)
      REAL*8 ZSS0(KDLON)
      REAL*8 ZSCAT(KDLON)
      REAL*8 ZTR(KDLON,2,KFLEV+1)
C
      INTEGER jl, jk, ja, jae, jkl, jklp1, jaj, jkm1, in
      REAL*8 ZTRAY, ZGAR, ZRATIO, ZFF, ZFACOA, ZCORAE
      REAL*8 ZMUE, ZGAP, ZWW, ZTO, ZDEN, ZMU1, ZDEN1
      REAL*8 ZBMU0, ZBMU1, ZRE11
C
C* Prescribed Data for Aerosols:
C
      REAL*8 TAUA(2,5), RPIZA(2,5), RCGA(2,5)
      SAVE TAUA, RPIZA, RCGA
      DATA ((TAUA(IN,JA),JA=1,5),IN=1,2) /
     S .730719, .912819, .725059, .745405, .682188 ,
     S .730719, .912819, .725059, .745405, .682188 /
      DATA ((RPIZA(IN,JA),JA=1,5),IN=1,2) /
     S .872212, .982545, .623143, .944887, .997975 ,
     S .872212, .982545, .623143, .944887, .997975 /
      DATA ((RCGA (IN,JA),JA=1,5),IN=1,2) /
     S .647596, .739002, .580845, .662657, .624246 ,
     S .647596, .739002, .580845, .662657, .624246 /
C     ------------------------------------------------------------------
C
C*         1.    OPTICAL PARAMETERS FOR AEROSOLS AND RAYLEIGH
C                --------------------------------------------
C
 100  CONTINUE
C
      DO 103 JK = 1 , KFLEV+1
      DO 102 JA = 1 , 6
      DO 101 JL = 1, KDLON
      PRJ(JL,JA,JK) = 0.
      PRK(JL,JA,JK) = 0.
 101  CONTINUE
 102  CONTINUE
 103  CONTINUE
C
      DO 108 JK = 1 , KFLEV
c-OB
c      DO 104 JL = 1, KDLON
c      PCGAZ(JL,JK) = 0.
c      PPIZAZ(JL,JK) =  0.
c      PTAUAZ(JL,JK) = 0.
c 104  CONTINUE
c-OB
c      DO 106 JAE=1,5
c      DO 105 JL = 1, KDLON
c      PTAUAZ(JL,JK)=PTAUAZ(JL,JK)
c     S        +PAER(JL,JK,JAE)*TAUA(KNU,JAE)
c      PPIZAZ(JL,JK)=PPIZAZ(JL,JK)+PAER(JL,JK,JAE)
c     S        * TAUA(KNU,JAE)*RPIZA(KNU,JAE)
c      PCGAZ(JL,JK) =  PCGAZ(JL,JK) +PAER(JL,JK,JAE)
c     S        * TAUA(KNU,JAE)*RPIZA(KNU,JAE)*RCGA(KNU,JAE)
c 105  CONTINUE
c 106  CONTINUE
c-OB
      DO 105 JL = 1, KDLON
      PTAUAZ(JL,JK)=flag_aer * tauae(JL,JK,KNU)
      PPIZAZ(JL,JK)=flag_aer * pizae(JL,JK,KNU)
      PCGAZ (JL,JK)=flag_aer * cgae(JL,JK,KNU)
 105  CONTINUE
C
      IF (flag_aer.GT.0) THEN
c-OB
      DO 107 JL = 1, KDLON
c         PCGAZ(JL,JK)=PCGAZ(JL,JK)/PPIZAZ(JL,JK)
c         PPIZAZ(JL,JK)=PPIZAZ(JL,JK)/PTAUAZ(JL,JK)
         ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
         ZRATIO = ZTRAY / (ZTRAY + PTAUAZ(JL,JK))
         ZGAR = PCGAZ(JL,JK)
         ZFF = ZGAR * ZGAR
         PTAUAZ(JL,JK)=ZTRAY+PTAUAZ(JL,JK)*(1.-PPIZAZ(JL,JK)*ZFF)
         PCGAZ(JL,JK) = ZGAR * (1. - ZRATIO) / (1. + ZGAR)
         PPIZAZ(JL,JK) =ZRATIO+(1.-ZRATIO)*PPIZAZ(JL,JK)*(1.-ZFF)
     S                       / (1. - PPIZAZ(JL,JK) * ZFF)
 107  CONTINUE
      ELSE
      DO JL = 1, KDLON
         ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
         PTAUAZ(JL,JK) = ZTRAY
         PCGAZ(JL,JK) = 0.
         PPIZAZ(JL,JK) = 1.-REPSCT
      END DO
      END IF   ! check flag_aer
c     107  CONTINUE
c      PRINT 9107,JK,((PAER(JL,JK,JAE),JAE=1,5)
c     $ ,PTAUAZ(JL,JK),PPIZAZ(JL,JK),PCGAZ(JL,JK),JL=1,KDLON)
c 9107 FORMAT(1X,'SWCLR_107',I3,8E12.5)
C
 108  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
C                ----------------------------------------------
C
 200  CONTINUE
C
      DO 201 JL = 1, KDLON
      ZR23(JL) = 0.
      ZC0I(JL,KFLEV+1) = 0.
      ZCLEAR(JL) = 1.
      ZSCAT(JL) = 0.
 201  CONTINUE
C
      JK = 1
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 202 JL = 1, KDLON
      ZFACOA = 1. - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
      ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
      ZR21(JL) = EXP(-ZCORAE   )
      ZSS0(JL) = 1.-ZR21(JL)
      ZCLE0(JL,JKL) = ZSS0(JL)
C
      IF (NOVLP.EQ.1) THEN
c* maximum-random
         ZCLEAR(JL) = ZCLEAR(JL)
     S                  *(1.0-MAX(ZSS0(JL),ZSCAT(JL)))
     S                  /(1.0-MIN(ZSCAT(JL),1.-ZEPSEC))
         ZC0I(JL,JKL) = 1.0 - ZCLEAR(JL)
         ZSCAT(JL) = ZSS0(JL)
      ELSE IF (NOVLP.EQ.2) THEN
C* maximum
         ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
         ZC0I(JL,JKL) = ZSCAT(JL)
      ELSE IF (NOVLP.EQ.3) THEN
c* random
         ZCLEAR(JL)=ZCLEAR(JL)*(1.0-ZSS0(JL))
         ZSCAT(JL) = 1.0 - ZCLEAR(JL)
         ZC0I(JL,JKL) = ZSCAT(JL)
      END IF
 202  CONTINUE
C
      DO 205 JK = 2 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 204 JL = 1, KDLON
      ZFACOA = 1. - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
      ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
      ZR21(JL) = EXP(-ZCORAE   )
      ZSS0(JL) = 1.-ZR21(JL)
      ZCLE0(JL,JKL) = ZSS0(JL)
c     
      IF (NOVLP.EQ.1) THEN
c* maximum-random
         ZCLEAR(JL) = ZCLEAR(JL)
     S                  *(1.0-MAX(ZSS0(JL),ZSCAT(JL)))
     S                  /(1.0-MIN(ZSCAT(JL),1.-ZEPSEC))
         ZC0I(JL,JKL) = 1.0 - ZCLEAR(JL)
         ZSCAT(JL) = ZSS0(JL)
      ELSE IF (NOVLP.EQ.2) THEN
C* maximum
         ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
         ZC0I(JL,JKL) = ZSCAT(JL)
      ELSE IF (NOVLP.EQ.3) THEN
c* random
         ZCLEAR(JL)=ZCLEAR(JL)*(1.0-ZSS0(JL))
         ZSCAT(JL) = 1.0 - ZCLEAR(JL)
         ZC0I(JL,JKL) = ZSCAT(JL)
      END IF                  
 204  CONTINUE
 205  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
C                -----------------------------------------------
C
 300  CONTINUE
C
      DO 301 JL = 1, KDLON
      PRAY1(JL,KFLEV+1) = 0.
      PRAY2(JL,KFLEV+1) = 0.
      PREFZ(JL,2,1) = PALBP(JL,KNU)
      PREFZ(JL,1,1) = PALBP(JL,KNU)
      PTRA1(JL,KFLEV+1) = 1.
      PTRA2(JL,KFLEV+1) = 1.
 301  CONTINUE
C
      DO 346 JK = 2 , KFLEV+1
      JKM1 = JK-1
      DO 342 JL = 1, KDLON
C
C
C     ------------------------------------------------------------------
C
C*         3.1  EQUIVALENT ZENITH ANGLE
C               -----------------------
C
 310  CONTINUE
C
      ZMUE = (1.-ZC0I(JL,JK)) * PSEC(JL)
     S            + ZC0I(JL,JK) * 1.66
      PRMU0(JL,JK) = 1./ZMUE
C
C
C     ------------------------------------------------------------------
C
C*         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
C               ----------------------------------------------------
C
 320  CONTINUE
C
      ZGAP = PCGAZ(JL,JKM1)
      ZBMU0 = 0.5 - 0.75 * ZGAP / ZMUE
      ZWW = PPIZAZ(JL,JKM1)
      ZTO = PTAUAZ(JL,JKM1)
      ZDEN = 1. + (1. - ZWW + ZBMU0 * ZWW) * ZTO * ZMUE
     S       + (1-ZWW) * (1. - ZWW +2.*ZBMU0*ZWW)*ZTO*ZTO*ZMUE*ZMUE
      PRAY1(JL,JKM1) = ZBMU0 * ZWW * ZTO * ZMUE / ZDEN
      PTRA1(JL,JKM1) = 1. / ZDEN
C
      ZMU1 = 0.5
      ZBMU1 = 0.5 - 0.75 * ZGAP * ZMU1
      ZDEN1= 1. + (1. - ZWW + ZBMU1 * ZWW) * ZTO / ZMU1
     S       + (1-ZWW) * (1. - ZWW +2.*ZBMU1*ZWW)*ZTO*ZTO/ZMU1/ZMU1
      PRAY2(JL,JKM1) = ZBMU1 * ZWW * ZTO / ZMU1 / ZDEN1
      PTRA2(JL,JKM1) = 1. / ZDEN1
C
C
C
      PREFZ(JL,1,JK) = (PRAY1(JL,JKM1)
     S               + PREFZ(JL,1,JKM1) * PTRA1(JL,JKM1)
     S               * PTRA2(JL,JKM1)
     S               / (1.-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))
C
      ZTR(JL,1,JKM1) = (PTRA1(JL,JKM1)
     S               / (1.-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))
C
      PREFZ(JL,2,JK) = (PRAY1(JL,JKM1)
     S               + PREFZ(JL,2,JKM1) * PTRA1(JL,JKM1)
     S               * PTRA2(JL,JKM1) )
C
      ZTR(JL,2,JKM1) = PTRA1(JL,JKM1) 
C
 342  CONTINUE
 346  CONTINUE
      DO 347 JL = 1, KDLON
      ZMUE = (1.-ZC0I(JL,1))*PSEC(JL)+ZC0I(JL,1)*1.66
      PRMU0(JL,1)=1./ZMUE
 347  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
C                 -------------------------------------------------
C
 350  CONTINUE
C
      IF (KNU.EQ.1) THEN
      JAJ = 2
      DO 351 JL = 1, KDLON
      PRJ(JL,JAJ,KFLEV+1) = 1.
      PRK(JL,JAJ,KFLEV+1) = PREFZ(JL, 1,KFLEV+1)
 351  CONTINUE
C
      DO 353 JK = 1 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 352 JL = 1, KDLON
      ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,  1,JKL)
      PRJ(JL,JAJ,JKL) = ZRE11
      PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,  1,JKL)
 352  CONTINUE
 353  CONTINUE
 354  CONTINUE
C
      ELSE
C
      DO 358 JAJ = 1 , 2
      DO 355 JL = 1, KDLON
      PRJ(JL,JAJ,KFLEV+1) = 1.
      PRK(JL,JAJ,KFLEV+1) = PREFZ(JL,JAJ,KFLEV+1)
 355  CONTINUE
C
      DO 357 JK = 1 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 356 JL = 1, KDLON
      ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,JAJ,JKL)
      PRJ(JL,JAJ,JKL) = ZRE11
      PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,JAJ,JKL)
 356  CONTINUE
 357  CONTINUE
 358  CONTINUE
C
      END IF
C
C     ------------------------------------------------------------------
C
      RETURN
      END
      SUBROUTINE SWR ( KNU
     S  , PALBD , PCG   , PCLD , PDSIG, POMEGA, PRAYL
     S  , PSEC  , PTAU
     S  , PCGAZ , PPIZAZ, PRAY1, PRAY2, PREFZ , PRJ  , PRK , PRMUE
     S  , PTAUAZ, PTRA1 , PTRA2 )
      use dimens_m
      use dimphy
      use raddim
      use radepsi
      use radopt
      IMPLICIT none
C
C     ------------------------------------------------------------------
C     PURPOSE.
C     --------
C           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
C     CONTINUUM SCATTERING
C
C     METHOD.
C     -------
C
C          1. COMPUTES CONTINUUM FLUXES CORRESPONDING TO AEROSOL
C     OR/AND RAYLEIGH SCATTERING (NO MOLECULAR GAS ABSORPTION)
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
C        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C     ------------------------------------------------------------------
C* ARGUMENTS:
C
      INTEGER KNU
      REAL*8 PALBD(KDLON,2)
      REAL*8 PCG(KDLON,2,KFLEV)
      REAL*8 PCLD(KDLON,KFLEV)
      REAL*8 PDSIG(KDLON,KFLEV)
      REAL*8 POMEGA(KDLON,2,KFLEV)
      REAL*8 PRAYL(KDLON)
      REAL*8 PSEC(KDLON)
      REAL*8 PTAU(KDLON,2,KFLEV)
C
      REAL*8 PRAY1(KDLON,KFLEV+1)
      REAL*8 PRAY2(KDLON,KFLEV+1)
      REAL*8 PREFZ(KDLON,2,KFLEV+1)
      REAL*8 PRJ(KDLON,6,KFLEV+1)
      REAL*8 PRK(KDLON,6,KFLEV+1)
      REAL*8 PRMUE(KDLON,KFLEV+1)
      REAL*8 PCGAZ(KDLON,KFLEV)
      REAL*8 PPIZAZ(KDLON,KFLEV)
      REAL*8 PTAUAZ(KDLON,KFLEV)
      REAL*8 PTRA1(KDLON,KFLEV+1)
      REAL*8 PTRA2(KDLON,KFLEV+1)
C
C* LOCAL VARIABLES:
C
      REAL*8 ZC1I(KDLON,KFLEV+1)
      REAL*8 ZCLEQ(KDLON,KFLEV)
      REAL*8 ZCLEAR(KDLON)
      REAL*8 ZCLOUD(KDLON)
      REAL*8 ZGG(KDLON)
      REAL*8 ZREF(KDLON)
      REAL*8 ZRE1(KDLON)
      REAL*8 ZRE2(KDLON)
      REAL*8 ZRMUZ(KDLON)
      REAL*8 ZRNEB(KDLON)
      REAL*8 ZR21(KDLON)
      REAL*8 ZR22(KDLON)
      REAL*8 ZR23(KDLON)
      REAL*8 ZSS1(KDLON)
      REAL*8 ZTO1(KDLON)
      REAL*8 ZTR(KDLON,2,KFLEV+1)
      REAL*8 ZTR1(KDLON)
      REAL*8 ZTR2(KDLON)
      REAL*8 ZW(KDLON)
C
      INTEGER jk, jl, ja, jkl, jklp1, jkm1, jaj
      REAL*8 ZFACOA, ZFACOC, ZCORAE, ZCORCD
      REAL*8 ZMUE, ZGAP, ZWW, ZTO, ZDEN, ZDEN1
      REAL*8 ZMU1, ZRE11, ZBMU0, ZBMU1
C
C     ------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
 100  CONTINUE
C
      DO 103 JK = 1 , KFLEV+1
      DO 102 JA = 1 , 6
      DO 101 JL = 1, KDLON
      PRJ(JL,JA,JK) = 0.
      PRK(JL,JA,JK) = 0.
 101  CONTINUE
 102  CONTINUE
 103  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
C                ----------------------------------------------
C
 200  CONTINUE
C
      DO 201 JL = 1, KDLON
      ZR23(JL) = 0.
      ZC1I(JL,KFLEV+1) = 0.
      ZCLEAR(JL) = 1.
      ZCLOUD(JL) = 0.
 201  CONTINUE
C
      JK = 1
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 202 JL = 1, KDLON
      ZFACOA = 1. - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
      ZFACOC = 1. - POMEGA(JL,KNU,JKL) * PCG(JL,KNU,JKL)
     S                                 * PCG(JL,KNU,JKL)
      ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
      ZCORCD = ZFACOC * PTAU(JL,KNU,JKL) * PSEC(JL)
      ZR21(JL) = EXP(-ZCORAE   )
      ZR22(JL) = EXP(-ZCORCD   )
      ZSS1(JL) = PCLD(JL,JKL)*(1.0-ZR21(JL)*ZR22(JL))
     S               + (1.0-PCLD(JL,JKL))*(1.0-ZR21(JL))
      ZCLEQ(JL,JKL) = ZSS1(JL)
C
      IF (NOVLP.EQ.1) THEN
c* maximum-random
         ZCLEAR(JL) = ZCLEAR(JL)
     S                  *(1.0-MAX(ZSS1(JL),ZCLOUD(JL)))
     S                  /(1.0-MIN(ZCLOUD(JL),1.-ZEPSEC))
         ZC1I(JL,JKL) = 1.0 - ZCLEAR(JL)
         ZCLOUD(JL) = ZSS1(JL)
      ELSE IF (NOVLP.EQ.2) THEN
C* maximum
         ZCLOUD(JL) = MAX( ZSS1(JL) , ZCLOUD(JL) )
         ZC1I(JL,JKL) = ZCLOUD(JL)
      ELSE IF (NOVLP.EQ.3) THEN
c* random
         ZCLEAR(JL) = ZCLEAR(JL)*(1.0 - ZSS1(JL))
         ZCLOUD(JL) = 1.0 - ZCLEAR(JL)
         ZC1I(JL,JKL) = ZCLOUD(JL)
      END IF
 202  CONTINUE
C
      DO 205 JK = 2 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 204 JL = 1, KDLON
      ZFACOA = 1. - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
      ZFACOC = 1. - POMEGA(JL,KNU,JKL) * PCG(JL,KNU,JKL)
     S                                 * PCG(JL,KNU,JKL)
      ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
      ZCORCD = ZFACOC * PTAU(JL,KNU,JKL) * PSEC(JL)
      ZR21(JL) = EXP(-ZCORAE   )
      ZR22(JL) = EXP(-ZCORCD   )
      ZSS1(JL) = PCLD(JL,JKL)*(1.0-ZR21(JL)*ZR22(JL))
     S               + (1.0-PCLD(JL,JKL))*(1.0-ZR21(JL))
      ZCLEQ(JL,JKL) = ZSS1(JL)
c     
      IF (NOVLP.EQ.1) THEN
c* maximum-random
         ZCLEAR(JL) = ZCLEAR(JL)
     S                  *(1.0-MAX(ZSS1(JL),ZCLOUD(JL)))
     S                  /(1.0-MIN(ZCLOUD(JL),1.-ZEPSEC))
         ZC1I(JL,JKL) = 1.0 - ZCLEAR(JL)
         ZCLOUD(JL) = ZSS1(JL)
      ELSE IF (NOVLP.EQ.2) THEN
C* maximum
         ZCLOUD(JL) = MAX( ZSS1(JL) , ZCLOUD(JL) )
         ZC1I(JL,JKL) = ZCLOUD(JL)
      ELSE IF (NOVLP.EQ.3) THEN
c* random
         ZCLEAR(JL) = ZCLEAR(JL)*(1.0 - ZSS1(JL))
         ZCLOUD(JL) = 1.0 - ZCLEAR(JL)
         ZC1I(JL,JKL) = ZCLOUD(JL)
      END IF
 204  CONTINUE
 205  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
C                -----------------------------------------------
C
 300  CONTINUE
C
      DO 301 JL = 1, KDLON
      PRAY1(JL,KFLEV+1) = 0.
      PRAY2(JL,KFLEV+1) = 0.
      PREFZ(JL,2,1) = PALBD(JL,KNU)
      PREFZ(JL,1,1) = PALBD(JL,KNU)
      PTRA1(JL,KFLEV+1) = 1.
      PTRA2(JL,KFLEV+1) = 1.
 301  CONTINUE
C
      DO 346 JK = 2 , KFLEV+1
      JKM1 = JK-1
      DO 342 JL = 1, KDLON
      ZRNEB(JL)= PCLD(JL,JKM1)
      ZRE1(JL)=0.
      ZTR1(JL)=0.
      ZRE2(JL)=0.
      ZTR2(JL)=0.
C
C
C     ------------------------------------------------------------------
C
C*         3.1  EQUIVALENT ZENITH ANGLE
C               -----------------------
C
 310  CONTINUE
C
      ZMUE = (1.-ZC1I(JL,JK)) * PSEC(JL)
     S            + ZC1I(JL,JK) * 1.66
      PRMUE(JL,JK) = 1./ZMUE
C
C
C     ------------------------------------------------------------------
C
C*         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
C               ----------------------------------------------------
C
 320  CONTINUE
C
      ZGAP = PCGAZ(JL,JKM1)
      ZBMU0 = 0.5 - 0.75 * ZGAP / ZMUE
      ZWW = PPIZAZ(JL,JKM1)
      ZTO = PTAUAZ(JL,JKM1)
      ZDEN = 1. + (1. - ZWW + ZBMU0 * ZWW) * ZTO * ZMUE
     S       + (1-ZWW) * (1. - ZWW +2.*ZBMU0*ZWW)*ZTO*ZTO*ZMUE*ZMUE
      PRAY1(JL,JKM1) = ZBMU0 * ZWW * ZTO * ZMUE / ZDEN
      PTRA1(JL,JKM1) = 1. / ZDEN
c      PRINT *,' LOOP 342 ** 3 ** JL=',JL,PRAY1(JL,JKM1),PTRA1(JL,JKM1)
C
      ZMU1 = 0.5
      ZBMU1 = 0.5 - 0.75 * ZGAP * ZMU1
      ZDEN1= 1. + (1. - ZWW + ZBMU1 * ZWW) * ZTO / ZMU1
     S       + (1-ZWW) * (1. - ZWW +2.*ZBMU1*ZWW)*ZTO*ZTO/ZMU1/ZMU1
      PRAY2(JL,JKM1) = ZBMU1 * ZWW * ZTO / ZMU1 / ZDEN1
      PTRA2(JL,JKM1) = 1. / ZDEN1
C
C
C     ------------------------------------------------------------------
C
C*         3.3  EFFECT OF CLOUD LAYER
C               ---------------------
C
 330  CONTINUE
C
      ZW(JL) = POMEGA(JL,KNU,JKM1)
      ZTO1(JL) = PTAU(JL,KNU,JKM1)/ZW(JL)
     S         + PTAUAZ(JL,JKM1)/PPIZAZ(JL,JKM1)
      ZR21(JL) = PTAU(JL,KNU,JKM1) + PTAUAZ(JL,JKM1)
      ZR22(JL) = PTAU(JL,KNU,JKM1) / ZR21(JL)
      ZGG(JL) = ZR22(JL) * PCG(JL,KNU,JKM1)
     S              + (1. - ZR22(JL)) * PCGAZ(JL,JKM1)
C Modif PhD - JJM 19/03/96 pour erreurs arrondis
C machine
C PHD PROTECTION ZW(JL) = ZR21(JL) / ZTO1(JL)
      IF (ZW(JL).EQ.1. .AND. PPIZAZ(JL,JKM1).EQ.1.) THEN
         ZW(JL)=1.
      ELSE
         ZW(JL) = ZR21(JL) / ZTO1(JL)
      END IF
      ZREF(JL) = PREFZ(JL,1,JKM1)
      ZRMUZ(JL) = PRMUE(JL,JK)
 342  CONTINUE
C
      CALL SWDE(ZGG  , ZREF  , ZRMUZ , ZTO1 , ZW,
     S          ZRE1 , ZRE2  , ZTR1  , ZTR2)
C
      DO 345 JL = 1, KDLON
C
      PREFZ(JL,1,JK) = (1.-ZRNEB(JL)) * (PRAY1(JL,JKM1)
     S               + PREFZ(JL,1,JKM1) * PTRA1(JL,JKM1)
     S               * PTRA2(JL,JKM1)
     S               / (1.-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))
     S               + ZRNEB(JL) * ZRE2(JL)
C
      ZTR(JL,1,JKM1) = ZRNEB(JL) * ZTR2(JL) + (PTRA1(JL,JKM1)
     S               / (1.-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))
     S               * (1.-ZRNEB(JL))
C
      PREFZ(JL,2,JK) = (1.-ZRNEB(JL)) * (PRAY1(JL,JKM1)
     S               + PREFZ(JL,2,JKM1) * PTRA1(JL,JKM1)
     S               * PTRA2(JL,JKM1) )
     S               + ZRNEB(JL) * ZRE1(JL)
C
      ZTR(JL,2,JKM1) = ZRNEB(JL) * ZTR1(JL)
     S               + PTRA1(JL,JKM1) * (1.-ZRNEB(JL))
C
 345  CONTINUE
 346  CONTINUE
      DO 347 JL = 1, KDLON
      ZMUE = (1.-ZC1I(JL,1))*PSEC(JL)+ZC1I(JL,1)*1.66
      PRMUE(JL,1)=1./ZMUE
 347  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
C                 -------------------------------------------------
C
 350  CONTINUE
C
      IF (KNU.EQ.1) THEN
      JAJ = 2
      DO 351 JL = 1, KDLON
      PRJ(JL,JAJ,KFLEV+1) = 1.
      PRK(JL,JAJ,KFLEV+1) = PREFZ(JL, 1,KFLEV+1)
 351  CONTINUE
C
      DO 353 JK = 1 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 352 JL = 1, KDLON
      ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,  1,JKL)
      PRJ(JL,JAJ,JKL) = ZRE11
      PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,  1,JKL)
 352  CONTINUE
 353  CONTINUE
 354  CONTINUE
C
      ELSE
C
      DO 358 JAJ = 1 , 2
      DO 355 JL = 1, KDLON
      PRJ(JL,JAJ,KFLEV+1) = 1.
      PRK(JL,JAJ,KFLEV+1) = PREFZ(JL,JAJ,KFLEV+1)
 355  CONTINUE
C
      DO 357 JK = 1 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 356 JL = 1, KDLON
      ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,JAJ,JKL)
      PRJ(JL,JAJ,JKL) = ZRE11
      PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,JAJ,JKL)
 356  CONTINUE
 357  CONTINUE
 358  CONTINUE
C
      END IF
C
C     ------------------------------------------------------------------
C
      RETURN
      END
      SUBROUTINE SWDE (PGG,PREF,PRMUZ,PTO1,PW,
     S                 PRE1,PRE2,PTR1,PTR2)
      use dimens_m
      use dimphy
      use raddim
      IMPLICIT none
C
C     ------------------------------------------------------------------
C     PURPOSE.
C     --------
C           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLOUDY
C     LAYER USING THE DELTA-EDDINGTON'S APPROXIMATION.
C
C     METHOD.
C     -------
C
C          STANDARD DELTA-EDDINGTON LAYER CALCULATIONS.
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 88-12-15
C     ------------------------------------------------------------------
C* ARGUMENTS:
C
      REAL*8 PGG(KDLON)   ! ASSYMETRY FACTOR
      REAL*8 PREF(KDLON)  ! REFLECTIVITY OF THE UNDERLYING LAYER
      REAL*8 PRMUZ(KDLON) ! COSINE OF SOLAR ZENITH ANGLE
      REAL*8 PTO1(KDLON)  ! OPTICAL THICKNESS
      REAL*8 PW(KDLON)    ! SINGLE SCATTERING ALBEDO
      REAL*8 PRE1(KDLON)  ! LAYER REFLECTIVITY (NO UNDERLYING-LAYER REFLECTION)
      REAL*8 PRE2(KDLON)  ! LAYER REFLECTIVITY
      REAL*8 PTR1(KDLON)  ! LAYER TRANSMISSIVITY (NO UNDERLYING-LAYER REFLECTION)
      REAL*8 PTR2(KDLON)  ! LAYER TRANSMISSIVITY
C
C* LOCAL VARIABLES:
C
      INTEGER jl
      REAL*8 ZFF, ZGP, ZTOP, ZWCP, ZDT, ZX1, ZWM
      REAL*8 ZRM2, ZRK, ZX2, ZRP, ZALPHA, ZBETA, ZARG
      REAL*8 ZEXMU0, ZARG2, ZEXKP, ZEXKM, ZXP2P, ZXM2P, ZAP2B, ZAM2B
      REAL*8 ZA11, ZA12, ZA13, ZA21, ZA22, ZA23
      REAL*8 ZDENA, ZC1A, ZC2A, ZRI0A, ZRI1A
      REAL*8 ZRI0B, ZRI1B
      REAL*8 ZB21, ZB22, ZB23, ZDENB, ZC1B, ZC2B
      REAL*8 ZRI0C, ZRI1C, ZRI0D, ZRI1D
C     ------------------------------------------------------------------
C
C*         1.      DELTA-EDDINGTON CALCULATIONS
C
 100  CONTINUE
C
      DO 131 JL   =   1, KDLON
C
C*         1.1     SET UP THE DELTA-MODIFIED PARAMETERS
C
 110  CONTINUE
C
      ZFF = PGG(JL)*PGG(JL)
      ZGP = PGG(JL)/(1.+PGG(JL))
      ZTOP = (1.- PW(JL) * ZFF) * PTO1(JL)
      ZWCP = (1-ZFF)* PW(JL) /(1.- PW(JL) * ZFF)
      ZDT = 2./3.
      ZX1 = 1.-ZWCP*ZGP
      ZWM = 1.-ZWCP
      ZRM2 =  PRMUZ(JL) * PRMUZ(JL)
      ZRK = SQRT(3.*ZWM*ZX1)
      ZX2 = 4.*(1.-ZRK*ZRK*ZRM2)
      ZRP=ZRK/ZX1
      ZALPHA = 3.*ZWCP*ZRM2*(1.+ZGP*ZWM)/ZX2
      ZBETA = 3.*ZWCP* PRMUZ(JL) *(1.+3.*ZGP*ZRM2*ZWM)/ZX2
CMAF      ZARG=MIN(ZTOP/PRMUZ(JL),200.)
      ZARG=MIN(ZTOP/PRMUZ(JL),2.0d+2)
      ZEXMU0=EXP(-ZARG)
CMAF      ZARG2=MIN(ZRK*ZTOP,200.)
      ZARG2=MIN(ZRK*ZTOP,2.0d+2)
      ZEXKP=EXP(ZARG2)
      ZEXKM = 1./ZEXKP
      ZXP2P = 1.+ZDT*ZRP
      ZXM2P = 1.-ZDT*ZRP
      ZAP2B = ZALPHA+ZDT*ZBETA
      ZAM2B = ZALPHA-ZDT*ZBETA
C
C*         1.2     WITHOUT REFLECTION FROM THE UNDERLYING LAYER
C
 120  CONTINUE
C
      ZA11 = ZXP2P
      ZA12 = ZXM2P
      ZA13 = ZAP2B
      ZA22 = ZXP2P*ZEXKP
      ZA21 = ZXM2P*ZEXKM
      ZA23 = ZAM2B*ZEXMU0
      ZDENA = ZA11 * ZA22 - ZA21 * ZA12
      ZC1A = (ZA22*ZA13-ZA12*ZA23)/ZDENA
      ZC2A = (ZA11*ZA23-ZA21*ZA13)/ZDENA
      ZRI0A = ZC1A+ZC2A-ZALPHA
      ZRI1A = ZRP*(ZC1A-ZC2A)-ZBETA
      PRE1(JL) = (ZRI0A-ZDT*ZRI1A)/ PRMUZ(JL)
      ZRI0B = ZC1A*ZEXKM+ZC2A*ZEXKP-ZALPHA*ZEXMU0
      ZRI1B = ZRP*(ZC1A*ZEXKM-ZC2A*ZEXKP)-ZBETA*ZEXMU0
      PTR1(JL) = ZEXMU0+(ZRI0B+ZDT*ZRI1B)/ PRMUZ(JL)
C
C*         1.3     WITH REFLECTION FROM THE UNDERLYING LAYER
C
 130  CONTINUE
C
      ZB21 = ZA21- PREF(JL) *ZXP2P*ZEXKM
      ZB22 = ZA22- PREF(JL) *ZXM2P*ZEXKP
      ZB23 = ZA23- PREF(JL) *ZEXMU0*(ZAP2B - PRMUZ(JL) )
      ZDENB = ZA11 * ZB22 - ZB21 * ZA12
      ZC1B = (ZB22*ZA13-ZA12*ZB23)/ZDENB
      ZC2B = (ZA11*ZB23-ZB21*ZA13)/ZDENB
      ZRI0C = ZC1B+ZC2B-ZALPHA
      ZRI1C = ZRP*(ZC1B-ZC2B)-ZBETA
      PRE2(JL) = (ZRI0C-ZDT*ZRI1C) / PRMUZ(JL)
      ZRI0D = ZC1B*ZEXKM + ZC2B*ZEXKP - ZALPHA*ZEXMU0
      ZRI1D = ZRP * (ZC1B*ZEXKM - ZC2B*ZEXKP) - ZBETA*ZEXMU0
      PTR2(JL) = ZEXMU0 + (ZRI0D + ZDT*ZRI1D) / PRMUZ(JL)
C
 131  CONTINUE
      RETURN
      END
      SUBROUTINE SWTT (KNU,KA,PU,PTR)
      use dimens_m
      use dimphy
      use raddim
      IMPLICIT none
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
C     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN THE TWO SPECTRAL
C     INTERVALS.
C
C     METHOD.
C     -------
C
C          TRANSMISSION FUNCTION ARE COMPUTED USING PADE APPROXIMANTS
C     AND HORNER'S ALGORITHM.
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 88-12-15
C-----------------------------------------------------------------------
C
C* ARGUMENTS
C
      INTEGER KNU     ! INDEX OF THE SPECTRAL INTERVAL
      INTEGER KA      ! INDEX OF THE ABSORBER
      REAL*8 PU(KDLON)  ! ABSORBER AMOUNT
C
      REAL*8 PTR(KDLON) ! TRANSMISSION FUNCTION
C
C* LOCAL VARIABLES:
C
      REAL*8 ZR1(KDLON), ZR2(KDLON)
      INTEGER jl, i,j
C
C* Prescribed Data:
C
      REAL*8 APAD(2,3,7), BPAD(2,3,7), D(2,3)
      SAVE APAD, BPAD, D
      DATA ((APAD(1,I,J),I=1,3),J=1,7) /
     S 0.912418292E+05, 0.000000000E-00, 0.925887084E-04,
     S 0.723613782E+05, 0.000000000E-00, 0.129353723E-01,
     S 0.596037057E+04, 0.000000000E-00, 0.800821928E+00,
     S 0.000000000E-00, 0.000000000E-00, 0.242715973E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.878331486E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.191559725E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.000000000E+00 /
      DATA ((APAD(2,I,J),I=1,3),J=1,7) /
     S 0.376655383E-08, 0.739646016E-08, 0.410177786E+03,
     S 0.978576773E-04, 0.131849595E-03, 0.672595424E+02,
     S 0.387714006E+00, 0.437772681E+00, 0.000000000E-00,
     S 0.118461660E+03, 0.151345118E+03, 0.000000000E-00,
     S 0.119079797E+04, 0.233628890E+04, 0.000000000E-00,
     S 0.293353397E+03, 0.797219934E+03, 0.000000000E-00,
     S 0.000000000E+00, 0.000000000E+00, 0.000000000E+00 /
C
      DATA ((BPAD(1,I,J),I=1,3),J=1,7) /
     S 0.912418292E+05, 0.000000000E-00, 0.925887084E-04,
     S 0.724555318E+05, 0.000000000E-00, 0.131812683E-01,
     S 0.602593328E+04, 0.000000000E-00, 0.812706117E+00,
     S 0.100000000E+01, 0.000000000E-00, 0.249863591E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.931071925E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.252233437E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.100000000E+01 /
      DATA ((BPAD(2,I,J),I=1,3),J=1,7) /
     S 0.376655383E-08, 0.739646016E-08, 0.410177786E+03,
     S 0.979023421E-04, 0.131861712E-03, 0.731185438E+02,
     S 0.388611139E+00, 0.437949001E+00, 0.100000000E+01,
     S 0.120291383E+03, 0.151692730E+03, 0.000000000E+00,
     S 0.130531005E+04, 0.237071130E+04, 0.000000000E+00,
     S 0.415049409E+03, 0.867914360E+03, 0.000000000E+00,
     S 0.100000000E+01, 0.100000000E+01, 0.000000000E+00 /
c
      DATA (D(1,I),I=1,3) / 0.00, 0.00, 0.00 /
      DATA (D(2,I),I=1,3) / 0.000000000, 0.000000000, 0.800000000 /
C
C-----------------------------------------------------------------------
C
C*         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION
C
 100  CONTINUE
C
      DO 201 JL = 1, KDLON
      ZR1(JL) = APAD(KNU,KA,1) + PU(JL) * (APAD(KNU,KA,2) + PU(JL)
     S      * ( APAD(KNU,KA,3) + PU(JL) * (APAD(KNU,KA,4) + PU(JL)
     S      * ( APAD(KNU,KA,5) + PU(JL) * (APAD(KNU,KA,6) + PU(JL)
     S      * ( APAD(KNU,KA,7) ))))))
C
      ZR2(JL) = BPAD(KNU,KA,1) + PU(JL) * (BPAD(KNU,KA,2) + PU(JL)
     S      * ( BPAD(KNU,KA,3) + PU(JL) * (BPAD(KNU,KA,4) + PU(JL)
     S      * ( BPAD(KNU,KA,5) + PU(JL) * (BPAD(KNU,KA,6) + PU(JL)
     S      * ( BPAD(KNU,KA,7) ))))))
C     
C
C*         2.      ADD THE BACKGROUND TRANSMISSION
C
 200  CONTINUE
C
C
      PTR(JL) = (ZR1(JL) / ZR2(JL)) * (1. - D(KNU,KA)) + D(KNU,KA)
 201  CONTINUE
C
      RETURN
      END
      SUBROUTINE SWTT1(KNU,KABS,KIND, PU, PTR)
      use dimens_m
      use dimphy
      use raddim
      IMPLICIT none
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
C     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN THE TWO SPECTRAL
C     INTERVALS.
C
C     METHOD.
C     -------
C
C          TRANSMISSION FUNCTION ARE COMPUTED USING PADE APPROXIMANTS
C     AND HORNER'S ALGORITHM.
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 95-01-20
C-----------------------------------------------------------------------
C* ARGUMENTS:
C
      INTEGER KNU          ! INDEX OF THE SPECTRAL INTERVAL
      INTEGER KABS         ! NUMBER OF ABSORBERS
      INTEGER KIND(KABS)   ! INDICES OF THE ABSORBERS
      REAL*8 PU(KDLON,KABS)  ! ABSORBER AMOUNT
C
      REAL*8 PTR(KDLON,KABS) ! TRANSMISSION FUNCTION
C
C* LOCAL VARIABLES:
C
      REAL*8 ZR1(KDLON)
      REAL*8 ZR2(KDLON)
      REAL*8 ZU(KDLON)
      INTEGER jl, ja, i, j, ia
C
C* Prescribed Data:
C
      REAL*8 APAD(2,3,7), BPAD(2,3,7), D(2,3)
      SAVE APAD, BPAD, D
      DATA ((APAD(1,I,J),I=1,3),J=1,7) /
     S 0.912418292E+05, 0.000000000E-00, 0.925887084E-04,
     S 0.723613782E+05, 0.000000000E-00, 0.129353723E-01,
     S 0.596037057E+04, 0.000000000E-00, 0.800821928E+00,
     S 0.000000000E-00, 0.000000000E-00, 0.242715973E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.878331486E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.191559725E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.000000000E+00 /
      DATA ((APAD(2,I,J),I=1,3),J=1,7) /
     S 0.376655383E-08, 0.739646016E-08, 0.410177786E+03,
     S 0.978576773E-04, 0.131849595E-03, 0.672595424E+02,
     S 0.387714006E+00, 0.437772681E+00, 0.000000000E-00,
     S 0.118461660E+03, 0.151345118E+03, 0.000000000E-00,
     S 0.119079797E+04, 0.233628890E+04, 0.000000000E-00,
     S 0.293353397E+03, 0.797219934E+03, 0.000000000E-00,
     S 0.000000000E+00, 0.000000000E+00, 0.000000000E+00 /
C
      DATA ((BPAD(1,I,J),I=1,3),J=1,7) /
     S 0.912418292E+05, 0.000000000E-00, 0.925887084E-04,
     S 0.724555318E+05, 0.000000000E-00, 0.131812683E-01,
     S 0.602593328E+04, 0.000000000E-00, 0.812706117E+00,
     S 0.100000000E+01, 0.000000000E-00, 0.249863591E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.931071925E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.252233437E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.100000000E+01 /
      DATA ((BPAD(2,I,J),I=1,3),J=1,7) /
     S 0.376655383E-08, 0.739646016E-08, 0.410177786E+03,
     S 0.979023421E-04, 0.131861712E-03, 0.731185438E+02,
     S 0.388611139E+00, 0.437949001E+00, 0.100000000E+01,
     S 0.120291383E+03, 0.151692730E+03, 0.000000000E+00,
     S 0.130531005E+04, 0.237071130E+04, 0.000000000E+00,
     S 0.415049409E+03, 0.867914360E+03, 0.000000000E+00,
     S 0.100000000E+01, 0.100000000E+01, 0.000000000E+00 /
c
      DATA (D(1,I),I=1,3) / 0.00, 0.00, 0.00 /
      DATA (D(2,I),I=1,3) / 0.000000000, 0.000000000, 0.800000000 /
C-----------------------------------------------------------------------
C
C*         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION
C
 100  CONTINUE
C
      DO 202 JA = 1,KABS
      IA=KIND(JA)
      DO 201 JL = 1, KDLON
      ZU(JL) = PU(JL,JA)
      ZR1(JL) = APAD(KNU,IA,1) + ZU(JL) * (APAD(KNU,IA,2) + ZU(JL)
     S      * ( APAD(KNU,IA,3) + ZU(JL) * (APAD(KNU,IA,4) + ZU(JL)
     S      * ( APAD(KNU,IA,5) + ZU(JL) * (APAD(KNU,IA,6) + ZU(JL)
     S      * ( APAD(KNU,IA,7) ))))))
C
      ZR2(JL) = BPAD(KNU,IA,1) + ZU(JL) * (BPAD(KNU,IA,2) + ZU(JL)
     S      * ( BPAD(KNU,IA,3) + ZU(JL) * (BPAD(KNU,IA,4) + ZU(JL)
     S      * ( BPAD(KNU,IA,5) + ZU(JL) * (BPAD(KNU,IA,6) + ZU(JL)
     S      * ( BPAD(KNU,IA,7) ))))))
C     
C
C*         2.      ADD THE BACKGROUND TRANSMISSION
C
 200  CONTINUE
C
      PTR(JL,JA) = (ZR1(JL)/ZR2(JL)) * (1.-D(KNU,IA)) + D(KNU,IA) 
 201  CONTINUE
 202  CONTINUE
C
      RETURN
      END
cIM ctes ds clesphys.h   SUBROUTINE LW(RCO2,RCH4,RN2O,RCFC11,RCFC12,
      SUBROUTINE LW(
     .              PPMB, PDP,
     .              PPSOL,PDT0,PEMIS,
     .              PTL, PTAVE, PWV, POZON, PAER,
     .              PCLDLD,PCLDLU,
     .              PVIEW,
     .              PCOLR, PCOLR0,
     .              PTOPLW,PSOLLW,PTOPLW0,PSOLLW0,
     .              psollwdown,
     .              plwup, plwdn, plwup0, plwdn0)
      use dimens_m
      use dimphy
      use clesphys
      use YOMCST
      use raddim
      IMPLICIT none
      include "raddimlw.h"
C
C-----------------------------------------------------------------------
C     METHOD.
C     -------
C
C          1. COMPUTES THE PRESSURE AND TEMPERATURE WEIGHTED AMOUNTS OF
C     ABSORBERS.
C          2. COMPUTES THE PLANCK FUNCTIONS ON THE INTERFACES AND THE
C     GRADIENT OF PLANCK FUNCTIONS IN THE LAYERS.
C          3. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING THE CON-
C     TRIBUTIONS OF THE ADJACENT AND DISTANT LAYERS AND THOSE FROM THE
C     BOUNDARIES.
C          4. COMPUTES THE CLEAR-SKY DOWNWARD AND UPWARD EMISSIVITIES.
C          5. INTRODUCES THE EFFECTS OF THE CLOUDS ON THE FLUXES.
C
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C-----------------------------------------------------------------------
cIM ctes ds clesphys.h
c     REAL*8 RCO2   ! CO2 CONCENTRATION (IPCC:353.E-06* 44.011/28.97)
c     REAL*8 RCH4   ! CH4 CONCENTRATION (IPCC: 1.72E-06* 16.043/28.97)
c     REAL*8 RN2O   ! N2O CONCENTRATION (IPCC: 310.E-09* 44.013/28.97)
c     REAL*8 RCFC11 ! CFC11 CONCENTRATION (IPCC: 280.E-12* 137.3686/28.97)
c     REAL*8 RCFC12 ! CFC12 CONCENTRATION (IPCC: 484.E-12* 120.9140/28.97)
      REAL*8 PCLDLD(KDLON,KFLEV)  ! DOWNWARD EFFECTIVE CLOUD COVER
      REAL*8 PCLDLU(KDLON,KFLEV)  ! UPWARD EFFECTIVE CLOUD COVER
      REAL*8 PDP(KDLON,KFLEV)     ! LAYER PRESSURE THICKNESS (Pa)
      REAL*8 PDT0(KDLON)          ! SURFACE TEMPERATURE DISCONTINUITY (K)
      REAL*8 PEMIS(KDLON)         ! SURFACE EMISSIVITY
      REAL*8 PPMB(KDLON,KFLEV+1)  ! HALF LEVEL PRESSURE (mb)
      REAL*8 PPSOL(KDLON)         ! SURFACE PRESSURE (Pa)
      REAL*8 POZON(KDLON,KFLEV)   ! O3 CONCENTRATION (kg/kg)
      REAL*8 PTL(KDLON,KFLEV+1)   ! HALF LEVEL TEMPERATURE (K)
      REAL*8 PAER(KDLON,KFLEV,5)  ! OPTICAL THICKNESS OF THE AEROSOLS
      REAL*8 PTAVE(KDLON,KFLEV)   ! LAYER TEMPERATURE (K)
      REAL*8 PVIEW(KDLON)         ! COSECANT OF VIEWING ANGLE
      REAL*8 PWV(KDLON,KFLEV)     ! SPECIFIC HUMIDITY (kg/kg)
C
      REAL*8 PCOLR(KDLON,KFLEV)   ! LONG-WAVE TENDENCY (K/day)
      REAL*8 PCOLR0(KDLON,KFLEV)  ! LONG-WAVE TENDENCY (K/day) clear-sky
      REAL*8 PTOPLW(KDLON)        ! LONGWAVE FLUX AT T.O.A.
      REAL*8 PSOLLW(KDLON)        ! LONGWAVE FLUX AT SURFACE
      REAL*8 PTOPLW0(KDLON)       ! LONGWAVE FLUX AT T.O.A. (CLEAR-SKY)
      REAL*8 PSOLLW0(KDLON)       ! LONGWAVE FLUX AT SURFACE (CLEAR-SKY)
c Rajout LF
      real*8 psollwdown(kdlon)    ! LONGWAVE downwards flux at surface
cIM
      REAL*8 plwup(KDLON,KFLEV+1)  ! LW up total sky
      REAL*8 plwup0(KDLON,KFLEV+1) ! LW up clear sky
      REAL*8 plwdn(KDLON,KFLEV+1)  ! LW down total sky
      REAL*8 plwdn0(KDLON,KFLEV+1) ! LW down clear sky
C-------------------------------------------------------------------------
      REAL*8 ZABCU(KDLON,NUA,3*KFLEV+1)
      REAL*8 ZOZ(KDLON,KFLEV)
c
      REAL*8 ZFLUX(KDLON,2,KFLEV+1) ! RADIATIVE FLUXES (1:up; 2:down)
      REAL*8 ZFLUC(KDLON,2,KFLEV+1) ! CLEAR-SKY RADIATIVE FLUXES
      REAL*8 ZBINT(KDLON,KFLEV+1)            ! Intermediate variable
      REAL*8 ZBSUI(KDLON)                    ! Intermediate variable
      REAL*8 ZCTS(KDLON,KFLEV)               ! Intermediate variable
      REAL*8 ZCNTRB(KDLON,KFLEV+1,KFLEV+1)   ! Intermediate variable
      SAVE ZFLUX, ZFLUC, ZBINT, ZBSUI, ZCTS, ZCNTRB
c
      INTEGER ilim, i, k, kpl1
C
      INTEGER lw0pas ! Every lw0pas steps, clear-sky is done
      PARAMETER (lw0pas=1)
      INTEGER lwpas  ! Every lwpas steps, cloudy-sky is done
      PARAMETER (lwpas=1)
c
      INTEGER itaplw0, itaplw
      LOGICAL appel1er
      SAVE appel1er, itaplw0, itaplw
      DATA appel1er /.TRUE./
      DATA itaplw0,itaplw /0,0/
C     ------------------------------------------------------------------
      IF (appel1er) THEN
         PRINT*, "LW clear-sky calling frequency: ", lw0pas
         PRINT*, "LW cloudy-sky calling frequency: ", lwpas
         PRINT*, "   In general, they should be 1"
         appel1er=.FALSE.
      ENDIF
C
      IF (MOD(itaplw0,lw0pas).EQ.0) THEN
      DO k = 1, KFLEV  ! convertir ozone de kg/kg en pa/pa
      DO i = 1, KDLON
c convertir ozone de kg/kg en pa (modif MPL 100505)
         ZOZ(i,k) = POZON(i,k)*PDP(i,k) * RMD/RMO3
c        print *,'LW: ZOZ*10**6=',ZOZ(i,k)*1000000.
      ENDDO
      ENDDO
cIM ctes ds clesphys.h   CALL LWU(RCO2,RCH4, RN2O, RCFC11, RCFC12,
      CALL LWU(
     S         PAER,PDP,PPMB,PPSOL,ZOZ,PTAVE,PVIEW,PWV,ZABCU)
      CALL LWBV(ILIM,PDP,PDT0,PEMIS,PPMB,PTL,PTAVE,ZABCU,
     S          ZFLUC,ZBINT,ZBSUI,ZCTS,ZCNTRB)
      itaplw0 = 0
      ENDIF
      itaplw0 = itaplw0 + 1
C
      IF (MOD(itaplw,lwpas).EQ.0) THEN
      CALL LWC(ILIM,PCLDLD,PCLDLU,PEMIS,
     S         ZFLUC,ZBINT,ZBSUI,ZCTS,ZCNTRB,
     S         ZFLUX)
      itaplw = 0
      ENDIF
      itaplw = itaplw + 1
C
      DO k = 1, KFLEV
         kpl1 = k+1
         DO i = 1, KDLON
            PCOLR(i,k) = ZFLUX(i,1,kpl1)+ZFLUX(i,2,kpl1)
     .                 - ZFLUX(i,1,k)-   ZFLUX(i,2,k)
            PCOLR(i,k) = PCOLR(i,k) * RDAY*RG/RCPD / PDP(i,k)
            PCOLR0(i,k) = ZFLUC(i,1,kpl1)+ZFLUC(i,2,kpl1)
     .                 - ZFLUC(i,1,k)-   ZFLUC(i,2,k)
            PCOLR0(i,k) = PCOLR0(i,k) * RDAY*RG/RCPD / PDP(i,k)
         ENDDO
      ENDDO
      DO i = 1, KDLON
         PSOLLW(i) = -ZFLUX(i,1,1)-ZFLUX(i,2,1)
         PTOPLW(i) = ZFLUX(i,1,KFLEV+1) + ZFLUX(i,2,KFLEV+1)
c
         PSOLLW0(i) = -ZFLUC(i,1,1)-ZFLUC(i,2,1)
         PTOPLW0(i) = ZFLUC(i,1,KFLEV+1) + ZFLUC(i,2,KFLEV+1)
         psollwdown(i) = -ZFLUX(i,2,1)
c
cIM attention aux signes !; LWtop >0, LWdn < 0
         DO k = 1, KFLEV+1
           plwup(i,k) = ZFLUX(i,1,k)
           plwup0(i,k) = ZFLUC(i,1,k)
           plwdn(i,k) = ZFLUX(i,2,k)
           plwdn0(i,k) = ZFLUC(i,2,k)
         ENDDO
      ENDDO
C     ------------------------------------------------------------------
      RETURN
      END
cIM ctes ds clesphys.h   SUBROUTINE LWU(RCO2, RCH4, RN2O, RCFC11, RCFC12,
      SUBROUTINE LWU(
     S               PAER,PDP,PPMB,PPSOL,POZ,PTAVE,PVIEW,PWV,
     S               PABCU)
      use dimens_m
      use dimphy
      use clesphys
      use YOMCST
      use raddim
      use radepsi
      use radopt
      IMPLICIT none
      include "raddimlw.h"
C
C     PURPOSE.
C     --------
C           COMPUTES ABSORBER AMOUNTS INCLUDING PRESSURE AND
C           TEMPERATURE EFFECTS
C
C     METHOD.
C     -------
C
C          1. COMPUTES THE PRESSURE AND TEMPERATURE WEIGHTED AMOUNTS OF
C     ABSORBERS.
C
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C        Voigt lines (loop 404 modified) - JJM & PhD - 01/96
C-----------------------------------------------------------------------
C* ARGUMENTS:
cIM ctes ds clesphys.h
c     REAL*8 RCO2
c     REAL*8 RCH4, RN2O, RCFC11, RCFC12
      REAL*8 PAER(KDLON,KFLEV,5)
      REAL*8 PDP(KDLON,KFLEV)
      REAL*8 PPMB(KDLON,KFLEV+1)
      REAL*8 PPSOL(KDLON)
      REAL*8 POZ(KDLON,KFLEV)
      REAL*8 PTAVE(KDLON,KFLEV)
      REAL*8 PVIEW(KDLON)
      REAL*8 PWV(KDLON,KFLEV)
C
      REAL*8 PABCU(KDLON,NUA,3*KFLEV+1) ! EFFECTIVE ABSORBER AMOUNTS
C
C-----------------------------------------------------------------------
C* LOCAL VARIABLES:
      REAL*8 ZABLY(KDLON,NUA,3*KFLEV+1)
      REAL*8 ZDUC(KDLON,3*KFLEV+1)
      REAL*8 ZPHIO(KDLON)
      REAL*8 ZPSC2(KDLON)
      REAL*8 ZPSC3(KDLON)
      REAL*8 ZPSH1(KDLON)
      REAL*8 ZPSH2(KDLON)
      REAL*8 ZPSH3(KDLON)
      REAL*8 ZPSH4(KDLON)
      REAL*8 ZPSH5(KDLON)
      REAL*8 ZPSH6(KDLON)
      REAL*8 ZPSIO(KDLON)
      REAL*8 ZTCON(KDLON)
      REAL*8 ZPHM6(KDLON)
      REAL*8 ZPSM6(KDLON)
      REAL*8 ZPHN6(KDLON)
      REAL*8 ZPSN6(KDLON)
      REAL*8 ZSSIG(KDLON,3*KFLEV+1)
      REAL*8 ZTAVI(KDLON)
      REAL*8 ZUAER(KDLON,Ninter)
      REAL*8 ZXOZ(KDLON)
      REAL*8 ZXWV(KDLON)
C
      INTEGER jl, jk, jkj, jkjr, jkjp, ig1
      INTEGER jki, jkip1, ja, jj
      INTEGER jkl, jkp1, jkk, jkjpn
      INTEGER jae1, jae2, jae3, jae, jjpn
      INTEGER ir, jc, jcp1
      REAL*8 zdpm, zupm, zupmh2o, zupmco2, zupmo3, zu6, zup
      REAL*8 zfppw, ztx, ztx2, zzably
      REAL*8 zcah1, zcbh1, zcah2, zcbh2, zcah3, zcbh3
      REAL*8 zcah4, zcbh4, zcah5, zcbh5, zcah6, zcbh6
      REAL*8 zcac8, zcbc8
      REAL*8 zalup, zdiff
c
      REAL*8 PVGCO2, PVGH2O, PVGO3
C
      REAL*8 R10E  ! DECIMAL/NATURAL LOG.FACTOR
      PARAMETER (R10E=0.4342945)
c
c Used Data Block:
c
      REAL*8 TREF
      SAVE TREF
      REAL*8 RT1(2)
      SAVE RT1
      REAL*8 RAER(5,5)
      SAVE RAER
      REAL*8 AT(8,3), BT(8,3)
      SAVE AT, BT
      REAL*8 OCT(4)
      SAVE OCT
      DATA TREF /250.0/
      DATA (RT1(IG1),IG1=1,2) / -0.577350269, +0.577350269 /
      DATA RAER / .038520, .037196, .040532, .054934, .038520
     1          , .12613 , .18313 , .10357 , .064106, .126130
     2          , .012579, .013649, .018652, .025181, .012579
     3          , .011890, .016142, .021105, .028908, .011890
     4          , .013792, .026810, .052203, .066338, .013792 /
      DATA (AT(1,IR),IR=1,3) /
     S 0.298199E-02,-.394023E-03,0.319566E-04 /
      DATA (BT(1,IR),IR=1,3) /
     S-0.106432E-04,0.660324E-06,0.174356E-06 /
      DATA (AT(2,IR),IR=1,3) /
     S 0.143676E-01,0.366501E-02,-.160822E-02 /
      DATA (BT(2,IR),IR=1,3) /
     S-0.553979E-04,-.101701E-04,0.920868E-05 /
      DATA (AT(3,IR),IR=1,3) /
     S 0.197861E-01,0.315541E-02,-.174547E-02 /
      DATA (BT(3,IR),IR=1,3) /
     S-0.877012E-04,0.513302E-04,0.523138E-06 /
      DATA (AT(4,IR),IR=1,3) /
     S 0.289560E-01,-.208807E-02,-.121943E-02 /
      DATA (BT(4,IR),IR=1,3) /
     S-0.165960E-03,0.157704E-03,-.146427E-04 /
      DATA (AT(5,IR),IR=1,3) /
     S 0.103800E-01,0.436296E-02,-.161431E-02 /
      DATA (BT(5,IR),IR=1,3) /
     S -.276744E-04,-.327381E-04,0.127646E-04 /
      DATA (AT(6,IR),IR=1,3) /
     S 0.868859E-02,-.972752E-03,0.000000E-00 /
      DATA (BT(6,IR),IR=1,3) /
     S -.278412E-04,-.713940E-06,0.117469E-05 /
      DATA (AT(7,IR),IR=1,3) /
     S 0.250073E-03,0.455875E-03,0.109242E-03 /
      DATA (BT(7,IR),IR=1,3) /
     S 0.199846E-05,-.216313E-05,0.175991E-06 /
      DATA (AT(8,IR),IR=1,3) /
     S 0.307423E-01,0.110879E-02,-.322172E-03 /
      DATA (BT(8,IR),IR=1,3) /
     S-0.108482E-03,0.258096E-05,-.814575E-06 /
c
      DATA OCT /-.326E-03, -.102E-05, .137E-02, -.535E-05/
C-----------------------------------------------------------------------
c
      IF (LEVOIGT) THEN
         PVGCO2= 60.
         PVGH2O= 30.
         PVGO3 =400.
      ELSE
         PVGCO2= 0.
         PVGH2O= 0.
         PVGO3 = 0.
      ENDIF
C
C
C*         2.    PRESSURE OVER GAUSS SUB-LEVELS
C                ------------------------------
C
 200  CONTINUE
C
      DO 201 JL = 1, KDLON
      ZSSIG(JL, 1 ) = PPMB(JL,1) * 100.
 201  CONTINUE
C
      DO 206 JK = 1 , KFLEV
      JKJ=(JK-1)*NG1P1+1
      JKJR = JKJ
      JKJP = JKJ + NG1P1
      DO 203 JL = 1, KDLON
      ZSSIG(JL,JKJP)=PPMB(JL,JK+1)* 100.
 203  CONTINUE
      DO 205 IG1=1,NG1
      JKJ=JKJ+1
      DO 204 JL = 1, KDLON
      ZSSIG(JL,JKJ)= (ZSSIG(JL,JKJR)+ZSSIG(JL,JKJP))*0.5
     S  + RT1(IG1) * (ZSSIG(JL,JKJP) - ZSSIG(JL,JKJR)) * 0.5
 204  CONTINUE
 205  CONTINUE
 206  CONTINUE
C
C-----------------------------------------------------------------------
C
C
C*         4.    PRESSURE THICKNESS AND MEAN PRESSURE OF SUB-LAYERS
C                --------------------------------------------------
C
 400  CONTINUE
C
      DO 402 JKI=1,3*KFLEV
      JKIP1=JKI+1
      DO 401 JL = 1, KDLON
      ZABLY(JL,5,JKI)=(ZSSIG(JL,JKI)+ZSSIG(JL,JKIP1))*0.5
      ZABLY(JL,3,JKI)=(ZSSIG(JL,JKI)-ZSSIG(JL,JKIP1))
     S                                 /(10.*RG)
 401  CONTINUE
 402  CONTINUE
C
      DO 406 JK = 1 , KFLEV
      JKP1=JK+1
      JKL = KFLEV+1 - JK
      DO 403 JL = 1, KDLON
      ZXWV(JL) = MAX (PWV(JL,JK) , ZEPSCQ )
      ZXOZ(JL) = MAX (POZ(JL,JK) / PDP(JL,JK) , ZEPSCO )
 403  CONTINUE
      JKJ=(JK-1)*NG1P1+1
      JKJPN=JKJ+NG1
      DO 405 JKK=JKJ,JKJPN
      DO 404 JL = 1, KDLON
      ZDPM = ZABLY(JL,3,JKK)
      ZUPM = ZABLY(JL,5,JKK)             * ZDPM / 101325.
      ZUPMCO2 = ( ZABLY(JL,5,JKK) + PVGCO2 ) * ZDPM / 101325.
      ZUPMH2O = ( ZABLY(JL,5,JKK) + PVGH2O ) * ZDPM / 101325.
      ZUPMO3  = ( ZABLY(JL,5,JKK) + PVGO3  ) * ZDPM / 101325.
      ZDUC(JL,JKK) = ZDPM
      ZABLY(JL,12,JKK) = ZXOZ(JL) * ZDPM
      ZABLY(JL,13,JKK) = ZXOZ(JL) * ZUPMO3
      ZU6 = ZXWV(JL) * ZUPM
      ZFPPW = 1.6078 * ZXWV(JL) / (1.+0.608*ZXWV(JL))
      ZABLY(JL,6,JKK) = ZXWV(JL) * ZUPMH2O
      ZABLY(JL,11,JKK) = ZU6 * ZFPPW
      ZABLY(JL,10,JKK) = ZU6 * (1.-ZFPPW)
      ZABLY(JL,9,JKK) = RCO2 * ZUPMCO2
      ZABLY(JL,8,JKK) = RCO2 * ZDPM
 404  CONTINUE
 405  CONTINUE
 406  CONTINUE
C
C-----------------------------------------------------------------------
C
C
C*         5.    CUMULATIVE ABSORBER AMOUNTS FROM TOP OF ATMOSPHERE
C                --------------------------------------------------
C
 500  CONTINUE
C
      DO 502 JA = 1, NUA
      DO 501 JL = 1, KDLON
      PABCU(JL,JA,3*KFLEV+1) = 0.
  501 CONTINUE
  502 CONTINUE
C
      DO 529 JK = 1 , KFLEV
      JJ=(JK-1)*NG1P1+1
      JJPN=JJ+NG1
      JKL=KFLEV+1-JK
C
C
C*         5.1  CUMULATIVE AEROSOL AMOUNTS FROM TOP OF ATMOSPHERE
C               --------------------------------------------------
C
 510  CONTINUE
C
      JAE1=3*KFLEV+1-JJ
      JAE2=3*KFLEV+1-(JJ+1)
      JAE3=3*KFLEV+1-JJPN
      DO 512 JAE=1,5
      DO 511 JL = 1, KDLON
      ZUAER(JL,JAE) = (RAER(JAE,1)*PAER(JL,JKL,1)
     S      +RAER(JAE,2)*PAER(JL,JKL,2)+RAER(JAE,3)*PAER(JL,JKL,3)
     S      +RAER(JAE,4)*PAER(JL,JKL,4)+RAER(JAE,5)*PAER(JL,JKL,5))
     S      /(ZDUC(JL,JAE1)+ZDUC(JL,JAE2)+ZDUC(JL,JAE3))
 511  CONTINUE
 512  CONTINUE
C
C
C
C*         5.2  INTRODUCES TEMPERATURE EFFECTS ON ABSORBER AMOUNTS
C               --------------------------------------------------
C
 520  CONTINUE
C
      DO 521 JL = 1, KDLON
      ZTAVI(JL)=PTAVE(JL,JKL)
      ZTCON(JL)=EXP(6.08*(296./ZTAVI(JL)-1.))
      ZTX=ZTAVI(JL)-TREF
      ZTX2=ZTX*ZTX
      ZZABLY = ZABLY(JL,6,JAE1)+ZABLY(JL,6,JAE2)+ZABLY(JL,6,JAE3)
CMAF      ZUP=MIN( MAX( 0.5*R10E*LOG( ZZABLY ) + 5., 0.), 6.0)
      ZUP=MIN( MAX( 0.5*R10E*LOG( ZZABLY ) + 5., 0.d+0), 6.d+0)
      ZCAH1=AT(1,1)+ZUP*(AT(1,2)+ZUP*(AT(1,3)))
      ZCBH1=BT(1,1)+ZUP*(BT(1,2)+ZUP*(BT(1,3)))
      ZPSH1(JL)=EXP( ZCAH1 * ZTX + ZCBH1 * ZTX2 )
      ZCAH2=AT(2,1)+ZUP*(AT(2,2)+ZUP*(AT(2,3)))
      ZCBH2=BT(2,1)+ZUP*(BT(2,2)+ZUP*(BT(2,3)))
      ZPSH2(JL)=EXP( ZCAH2 * ZTX + ZCBH2 * ZTX2 )
      ZCAH3=AT(3,1)+ZUP*(AT(3,2)+ZUP*(AT(3,3)))
      ZCBH3=BT(3,1)+ZUP*(BT(3,2)+ZUP*(BT(3,3)))
      ZPSH3(JL)=EXP( ZCAH3 * ZTX + ZCBH3 * ZTX2 )
      ZCAH4=AT(4,1)+ZUP*(AT(4,2)+ZUP*(AT(4,3)))
      ZCBH4=BT(4,1)+ZUP*(BT(4,2)+ZUP*(BT(4,3)))
      ZPSH4(JL)=EXP( ZCAH4 * ZTX + ZCBH4 * ZTX2 )
      ZCAH5=AT(5,1)+ZUP*(AT(5,2)+ZUP*(AT(5,3)))
      ZCBH5=BT(5,1)+ZUP*(BT(5,2)+ZUP*(BT(5,3)))
      ZPSH5(JL)=EXP( ZCAH5 * ZTX + ZCBH5 * ZTX2 )
      ZCAH6=AT(6,1)+ZUP*(AT(6,2)+ZUP*(AT(6,3)))
      ZCBH6=BT(6,1)+ZUP*(BT(6,2)+ZUP*(BT(6,3)))
      ZPSH6(JL)=EXP( ZCAH6 * ZTX + ZCBH6 * ZTX2 )
      ZPHM6(JL)=EXP(-5.81E-4 * ZTX - 1.13E-6 * ZTX2 )
      ZPSM6(JL)=EXP(-5.57E-4 * ZTX - 3.30E-6 * ZTX2 )
      ZPHN6(JL)=EXP(-3.46E-5 * ZTX + 2.05E-7 * ZTX2 )
      ZPSN6(JL)=EXP( 3.70E-3 * ZTX - 2.30E-6 * ZTX2 )
 521  CONTINUE
C
      DO 522 JL = 1, KDLON
      ZTAVI(JL)=PTAVE(JL,JKL)
      ZTX=ZTAVI(JL)-TREF
      ZTX2=ZTX*ZTX
      ZZABLY = ZABLY(JL,9,JAE1)+ZABLY(JL,9,JAE2)+ZABLY(JL,9,JAE3)
      ZALUP = R10E * LOG ( ZZABLY )
CMAF      ZUP   = MAX( 0.0 , 5.0 + 0.5 * ZALUP )
      ZUP   = MAX( 0.d+0 , 5.0 + 0.5 * ZALUP )
      ZPSC2(JL) = (ZTAVI(JL)/TREF) ** ZUP
      ZCAC8=AT(8,1)+ZUP*(AT(8,2)+ZUP*(AT(8,3)))
      ZCBC8=BT(8,1)+ZUP*(BT(8,2)+ZUP*(BT(8,3)))
      ZPSC3(JL)=EXP( ZCAC8 * ZTX + ZCBC8 * ZTX2 )
      ZPHIO(JL) = EXP( OCT(1) * ZTX + OCT(2) * ZTX2)
      ZPSIO(JL) = EXP( 2.* (OCT(3)*ZTX+OCT(4)*ZTX2))
 522  CONTINUE
C
      DO 524 JKK=JJ,JJPN
      JC=3*KFLEV+1-JKK
      JCP1=JC+1
      DO 523 JL = 1, KDLON
      ZDIFF = PVIEW(JL)
      PABCU(JL,10,JC)=PABCU(JL,10,JCP1)
     S                +ZABLY(JL,10,JC)           *ZDIFF
      PABCU(JL,11,JC)=PABCU(JL,11,JCP1)
     S                +ZABLY(JL,11,JC)*ZTCON(JL)*ZDIFF
C
      PABCU(JL,12,JC)=PABCU(JL,12,JCP1)
     S                +ZABLY(JL,12,JC)*ZPHIO(JL)*ZDIFF
      PABCU(JL,13,JC)=PABCU(JL,13,JCP1)
     S                +ZABLY(JL,13,JC)*ZPSIO(JL)*ZDIFF
C
      PABCU(JL,7,JC)=PABCU(JL,7,JCP1)
     S               +ZABLY(JL,9,JC)*ZPSC2(JL)*ZDIFF
      PABCU(JL,8,JC)=PABCU(JL,8,JCP1)
     S               +ZABLY(JL,9,JC)*ZPSC3(JL)*ZDIFF
      PABCU(JL,9,JC)=PABCU(JL,9,JCP1)
     S               +ZABLY(JL,9,JC)*ZPSC3(JL)*ZDIFF
C
      PABCU(JL,1,JC)=PABCU(JL,1,JCP1)
     S               +ZABLY(JL,6,JC)*ZPSH1(JL)*ZDIFF
      PABCU(JL,2,JC)=PABCU(JL,2,JCP1)
     S               +ZABLY(JL,6,JC)*ZPSH2(JL)*ZDIFF
      PABCU(JL,3,JC)=PABCU(JL,3,JCP1)
     S               +ZABLY(JL,6,JC)*ZPSH5(JL)*ZDIFF
      PABCU(JL,4,JC)=PABCU(JL,4,JCP1)
     S               +ZABLY(JL,6,JC)*ZPSH3(JL)*ZDIFF
      PABCU(JL,5,JC)=PABCU(JL,5,JCP1)
     S               +ZABLY(JL,6,JC)*ZPSH4(JL)*ZDIFF
      PABCU(JL,6,JC)=PABCU(JL,6,JCP1)
     S               +ZABLY(JL,6,JC)*ZPSH6(JL)*ZDIFF
C
      PABCU(JL,14,JC)=PABCU(JL,14,JCP1)
     S                +ZUAER(JL,1)    *ZDUC(JL,JC)*ZDIFF
      PABCU(JL,15,JC)=PABCU(JL,15,JCP1)
     S                +ZUAER(JL,2)    *ZDUC(JL,JC)*ZDIFF
      PABCU(JL,16,JC)=PABCU(JL,16,JCP1)
     S                +ZUAER(JL,3)    *ZDUC(JL,JC)*ZDIFF
      PABCU(JL,17,JC)=PABCU(JL,17,JCP1)
     S                +ZUAER(JL,4)    *ZDUC(JL,JC)*ZDIFF
      PABCU(JL,18,JC)=PABCU(JL,18,JCP1)
     S                +ZUAER(JL,5)    *ZDUC(JL,JC)*ZDIFF
C
      PABCU(JL,19,JC)=PABCU(JL,19,JCP1)
     S               +ZABLY(JL,8,JC)*RCH4/RCO2*ZPHM6(JL)*ZDIFF
      PABCU(JL,20,JC)=PABCU(JL,20,JCP1)
     S               +ZABLY(JL,9,JC)*RCH4/RCO2*ZPSM6(JL)*ZDIFF
      PABCU(JL,21,JC)=PABCU(JL,21,JCP1)
     S               +ZABLY(JL,8,JC)*RN2O/RCO2*ZPHN6(JL)*ZDIFF
      PABCU(JL,22,JC)=PABCU(JL,22,JCP1)
     S               +ZABLY(JL,9,JC)*RN2O/RCO2*ZPSN6(JL)*ZDIFF
C
      PABCU(JL,23,JC)=PABCU(JL,23,JCP1)
     S               +ZABLY(JL,8,JC)*RCFC11/RCO2         *ZDIFF
      PABCU(JL,24,JC)=PABCU(JL,24,JCP1)
     S               +ZABLY(JL,8,JC)*RCFC12/RCO2         *ZDIFF
 523  CONTINUE
 524  CONTINUE
C
 529  CONTINUE
C
C
      RETURN
      END
      SUBROUTINE LWBV(KLIM,PDP,PDT0,PEMIS,PPMB,PTL,PTAVE,PABCU,
     S                PFLUC,PBINT,PBSUI,PCTS,PCNTRB)
      use dimens_m
      use dimphy
      use YOMCST
      use raddim
      IMPLICIT none
      include "raddimlw.h"
C
C     PURPOSE.
C     --------
C           TO COMPUTE THE PLANCK FUNCTION AND PERFORM THE
C           VERTICAL INTEGRATION. SPLIT OUT FROM LW FOR MEMORY
C           SAVING
C
C     METHOD.
C     -------
C
C          1. COMPUTES THE PLANCK FUNCTIONS ON THE INTERFACES AND THE
C     GRADIENT OF PLANCK FUNCTIONS IN THE LAYERS.
C          2. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING THE CON-
C     TRIBUTIONS OF THE ADJACENT AND DISTANT LAYERS AND THOSE FROM THE
C     BOUNDARIES.
C          3. COMPUTES THE CLEAR-SKY COOLING RATES.
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C        MODIFICATION : 93-10-15 M.HAMRUD (SPLIT OUT FROM LW TO SAVE
C                                          MEMORY)
C-----------------------------------------------------------------------
C* ARGUMENTS:
      INTEGER KLIM
C
      REAL*8 PDP(KDLON,KFLEV)
      REAL*8 PDT0(KDLON)
      REAL*8 PEMIS(KDLON)
      REAL*8 PPMB(KDLON,KFLEV+1)
      REAL*8 PTL(KDLON,KFLEV+1)
      REAL*8 PTAVE(KDLON,KFLEV)
C
      REAL*8 PFLUC(KDLON,2,KFLEV+1)
C     
      REAL*8 PABCU(KDLON,NUA,3*KFLEV+1)
      REAL*8 PBINT(KDLON,KFLEV+1)
      REAL*8 PBSUI(KDLON)
      REAL*8 PCTS(KDLON,KFLEV)
      REAL*8 PCNTRB(KDLON,KFLEV+1,KFLEV+1)
C
C-------------------------------------------------------------------------
C
C* LOCAL VARIABLES:
      REAL*8 ZB(KDLON,Ninter,KFLEV+1)
      REAL*8 ZBSUR(KDLON,Ninter)
      REAL*8 ZBTOP(KDLON,Ninter)
      REAL*8 ZDBSL(KDLON,Ninter,KFLEV*2)
      REAL*8 ZGA(KDLON,8,2,KFLEV)
      REAL*8 ZGB(KDLON,8,2,KFLEV)
      REAL*8 ZGASUR(KDLON,8,2)
      REAL*8 ZGBSUR(KDLON,8,2)
      REAL*8 ZGATOP(KDLON,8,2)
      REAL*8 ZGBTOP(KDLON,8,2)
C
      INTEGER nuaer, ntraer
C     ------------------------------------------------------------------
C* COMPUTES PLANCK FUNCTIONS:
       CALL LWB(PDT0,PTAVE,PTL,
     S          ZB,PBINT,PBSUI,ZBSUR,ZBTOP,ZDBSL,
     S          ZGA,ZGB,ZGASUR,ZGBSUR,ZGATOP,ZGBTOP)
C     ------------------------------------------------------------------
C* PERFORMS THE VERTICAL INTEGRATION:
      NUAER = NUA
      NTRAER = NTRA
      CALL LWV(NUAER,NTRAER, KLIM
     R  , PABCU,ZB,PBINT,PBSUI,ZBSUR,ZBTOP,ZDBSL,PEMIS,PPMB,PTAVE
     R  , ZGA,ZGB,ZGASUR,ZGBSUR,ZGATOP,ZGBTOP
     S  , PCNTRB,PCTS,PFLUC)
C     ------------------------------------------------------------------
      RETURN
      END
      SUBROUTINE LWC(KLIM,PCLDLD,PCLDLU,PEMIS,PFLUC,
     R               PBINT,PBSUIN,PCTS,PCNTRB,
     S               PFLUX)
      use dimens_m
      use dimphy
      use raddim
      use radepsi
      use radopt
      IMPLICIT none
C
C     PURPOSE.
C     --------
C           INTRODUCES CLOUD EFFECTS ON LONGWAVE FLUXES OR
C           RADIANCES
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C     ==== INPUTS ===
C PBINT  : (KDLON,0:KFLEV)     ; HALF LEVEL PLANCK FUNCTION
C PBSUIN : (KDLON)             ; SURFACE PLANCK FUNCTION
C PCLDLD : (KDLON,KFLEV)       ; DOWNWARD EFFECTIVE CLOUD FRACTION
C PCLDLU : (KDLON,KFLEV)       ; UPWARD EFFECTIVE CLOUD FRACTION
C PCNTRB : (KDLON,KFLEV+1,KFLEV+1); CLEAR-SKY ENERGY EXCHANGE
C PCTS   : (KDLON,KFLEV)       ; CLEAR-SKY LAYER COOLING-TO-SPACE
C PEMIS  : (KDLON)             ; SURFACE EMISSIVITY
C PFLUC
C     ==== OUTPUTS ===
C PFLUX(KDLON,2,KFLEV)         ; RADIATIVE FLUXES :
C                     1  ==>  UPWARD   FLUX TOTAL
C                     2  ==>  DOWNWARD FLUX TOTAL
C
C     METHOD.
C     -------
C
C          1. INITIALIZES ALL FLUXES TO CLEAR-SKY VALUES
C          2. EFFECT OF ONE OVERCAST UNITY EMISSIVITY CLOUD LAYER
C          3. EFFECT OF SEMI-TRANSPARENT, PARTIAL OR MULTI-LAYERED
C     CLOUDS
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C        Voigt lines (loop 231 to 233)  - JJM & PhD - 01/96
C-----------------------------------------------------------------------
C* ARGUMENTS:
      INTEGER klim
      REAL*8 PFLUC(KDLON,2,KFLEV+1) ! CLEAR-SKY RADIATIVE FLUXES
      REAL*8 PBINT(KDLON,KFLEV+1)   ! HALF LEVEL PLANCK FUNCTION
      REAL*8 PBSUIN(KDLON)          ! SURFACE PLANCK FUNCTION
      REAL*8 PCNTRB(KDLON,KFLEV+1,KFLEV+1) !CLEAR-SKY ENERGY EXCHANGE
      REAL*8 PCTS(KDLON,KFLEV)      ! CLEAR-SKY LAYER COOLING-TO-SPACE
c
      REAL*8 PCLDLD(KDLON,KFLEV)
      REAL*8 PCLDLU(KDLON,KFLEV)
      REAL*8 PEMIS(KDLON)
C
      REAL*8 PFLUX(KDLON,2,KFLEV+1)
C-----------------------------------------------------------------------
C* LOCAL VARIABLES:
      INTEGER IMX(KDLON), IMXP(KDLON)
C
      REAL*8 ZCLEAR(KDLON),ZCLOUD(KDLON),ZDNF(KDLON,KFLEV+1,KFLEV+1)
     S  , ZFD(KDLON), ZFN10(KDLON), ZFU(KDLON)
     S  , ZUPF(KDLON,KFLEV+1,KFLEV+1)
      REAL*8 ZCLM(KDLON,KFLEV+1,KFLEV+1)
C
      INTEGER jk, jl, imaxc, imx1, imx2, jkj, jkp1, jkm1
      INTEGER jk1, jk2, jkc, jkcp1, jcloud
      INTEGER imxm1, imxp1
      REAL*8 zcfrac
C     ------------------------------------------------------------------
C
C*         1.     INITIALIZATION
C                 --------------
C
 100  CONTINUE
C
      IMAXC = 0
C
      DO 101 JL = 1, KDLON
      IMX(JL)=0
      IMXP(JL)=0
      ZCLOUD(JL) = 0.
 101  CONTINUE
C
C*         1.1    SEARCH THE LAYER INDEX OF THE HIGHEST CLOUD
C                 -------------------------------------------
C
 110  CONTINUE
C
      DO 112 JK = 1 , KFLEV
      DO 111 JL = 1, KDLON
      IMX1=IMX(JL)
      IMX2=JK
      IF (PCLDLU(JL,JK).GT.ZEPSC) THEN
         IMXP(JL)=IMX2
      ELSE
         IMXP(JL)=IMX1
      END IF
      IMAXC=MAX(IMXP(JL),IMAXC)
      IMX(JL)=IMXP(JL)
 111  CONTINUE
 112  CONTINUE
CGM*******
      IMAXC=KFLEV
CGM*******
C
      DO 114 JK = 1 , KFLEV+1
      DO 113 JL = 1, KDLON
      PFLUX(JL,1,JK) = PFLUC(JL,1,JK)
      PFLUX(JL,2,JK) = PFLUC(JL,2,JK)
 113  CONTINUE
 114  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.      EFFECT OF CLOUDINESS ON LONGWAVE FLUXES
C                  ---------------------------------------
C
      IF (IMAXC.GT.0) THEN
C
         IMXP1 = IMAXC + 1
         IMXM1 = IMAXC - 1
C
C*         2.0     INITIALIZE TO CLEAR-SKY FLUXES
C                  ------------------------------
C
 200  CONTINUE
C
         DO 203 JK1=1,KFLEV+1
         DO 202 JK2=1,KFLEV+1
         DO 201 JL = 1, KDLON
         ZUPF(JL,JK2,JK1)=PFLUC(JL,1,JK1)
         ZDNF(JL,JK2,JK1)=PFLUC(JL,2,JK1)
 201     CONTINUE
 202     CONTINUE
 203     CONTINUE
C
C*         2.1     FLUXES FOR ONE OVERCAST UNITY EMISSIVITY CLOUD
C                  ----------------------------------------------
C
 210  CONTINUE
C
         DO 213 JKC = 1 , IMAXC
         JCLOUD=JKC
         JKCP1=JCLOUD+1
C
C*         2.1.1   ABOVE THE CLOUD
C                  ---------------
C
 2110 CONTINUE
C
         DO 2115 JK=JKCP1,KFLEV+1
         JKM1=JK-1
         DO 2111 JL = 1, KDLON
         ZFU(JL)=0.
 2111    CONTINUE
         IF (JK .GT. JKCP1) THEN
            DO 2113 JKJ=JKCP1,JKM1
            DO 2112 JL = 1, KDLON
            ZFU(JL) = ZFU(JL) + PCNTRB(JL,JK,JKJ)
 2112       CONTINUE
 2113       CONTINUE
         END IF
C
         DO 2114 JL = 1, KDLON
         ZUPF(JL,JKCP1,JK)=PBINT(JL,JK)-ZFU(JL)
 2114    CONTINUE
 2115    CONTINUE
C
C*         2.1.2   BELOW THE CLOUD
C                  ---------------
C
 2120 CONTINUE
C
         DO 2125 JK=1,JCLOUD
         JKP1=JK+1
         DO 2121 JL = 1, KDLON
         ZFD(JL)=0.
 2121    CONTINUE
C
         IF (JK .LT. JCLOUD) THEN
            DO 2123 JKJ=JKP1,JCLOUD
            DO 2122 JL = 1, KDLON
            ZFD(JL) = ZFD(JL) + PCNTRB(JL,JK,JKJ)
 2122       CONTINUE
 2123       CONTINUE
         END IF
         DO 2124 JL = 1, KDLON
         ZDNF(JL,JKCP1,JK)=-PBINT(JL,JK)-ZFD(JL)
 2124    CONTINUE
 2125    CONTINUE
C
 213     CONTINUE
C
C
C*         2.2     CLOUD COVER MATRIX
C                  ------------------
C
C*    ZCLM(JK1,JK2) IS THE OBSCURATION FACTOR BY CLOUD LAYERS BETWEEN
C     HALF-LEVELS JK1 AND JK2 AS SEEN FROM JK1
C
 220  CONTINUE
C
      DO 223 JK1 = 1 , KFLEV+1
      DO 222 JK2 = 1 , KFLEV+1
      DO 221 JL = 1, KDLON
      ZCLM(JL,JK1,JK2) = 0.
 221  CONTINUE
 222  CONTINUE
 223  CONTINUE
C
C
C
C*         2.4     CLOUD COVER BELOW THE LEVEL OF CALCULATION
C                  ------------------------------------------
C
 240  CONTINUE
C
      DO 244 JK1 = 2 , KFLEV+1
      DO 241 JL = 1, KDLON
      ZCLEAR(JL)=1.
      ZCLOUD(JL)=0.
 241  CONTINUE
      DO 243 JK = JK1 - 1 , 1 , -1
      DO 242 JL = 1, KDLON
      IF (NOVLP.EQ.1) THEN
c* maximum-random       
         ZCLEAR(JL)=ZCLEAR(JL)*(1.0-MAX(PCLDLU(JL,JK),ZCLOUD(JL)))
     *                        /(1.0-MIN(ZCLOUD(JL),1.-ZEPSEC))
         ZCLM(JL,JK1,JK) = 1.0 - ZCLEAR(JL)
         ZCLOUD(JL) = PCLDLU(JL,JK)
      ELSE IF (NOVLP.EQ.2) THEN 
c* maximum      
         ZCLOUD(JL) = MAX(ZCLOUD(JL) , PCLDLU(JL,JK))
         ZCLM(JL,JK1,JK) = ZCLOUD(JL)
      ELSE IF (NOVLP.EQ.3) THEN
c* random      
         ZCLEAR(JL) = ZCLEAR(JL)*(1.0 - PCLDLU(JL,JK))
         ZCLOUD(JL) = 1.0 - ZCLEAR(JL)
         ZCLM(JL,JK1,JK) = ZCLOUD(JL)
      END IF
 242  CONTINUE
 243  CONTINUE
 244  CONTINUE
C
C
C*         2.5     CLOUD COVER ABOVE THE LEVEL OF CALCULATION
C                  ------------------------------------------
C
 250  CONTINUE
C
      DO 254 JK1 = 1 , KFLEV
      DO 251 JL = 1, KDLON
      ZCLEAR(JL)=1.
      ZCLOUD(JL)=0.
 251  CONTINUE
      DO 253 JK = JK1 , KFLEV
      DO 252 JL = 1, KDLON
      IF (NOVLP.EQ.1) THEN
c* maximum-random       
         ZCLEAR(JL)=ZCLEAR(JL)*(1.0-MAX(PCLDLD(JL,JK),ZCLOUD(JL)))
     *                        /(1.0-MIN(ZCLOUD(JL),1.-ZEPSEC))
         ZCLM(JL,JK1,JK) = 1.0 - ZCLEAR(JL)
         ZCLOUD(JL) = PCLDLD(JL,JK)
      ELSE IF (NOVLP.EQ.2) THEN 
c* maximum      
         ZCLOUD(JL) = MAX(ZCLOUD(JL) , PCLDLD(JL,JK))
         ZCLM(JL,JK1,JK) = ZCLOUD(JL)
      ELSE IF (NOVLP.EQ.3) THEN
c* random      
         ZCLEAR(JL) = ZCLEAR(JL)*(1.0 - PCLDLD(JL,JK))
         ZCLOUD(JL) = 1.0 - ZCLEAR(JL)
         ZCLM(JL,JK1,JK) = ZCLOUD(JL)
      END IF
 252  CONTINUE
 253  CONTINUE
 254  CONTINUE
C
C
C
C*         3.      FLUXES FOR PARTIAL/MULTIPLE LAYERED CLOUDINESS
C                  ----------------------------------------------
C
 300  CONTINUE
C
C*         3.1     DOWNWARD FLUXES
C                  ---------------
C
 310  CONTINUE
C
      DO 311 JL = 1, KDLON
      PFLUX(JL,2,KFLEV+1) = 0.
 311  CONTINUE
C
      DO 317 JK1 = KFLEV , 1 , -1
C
C*                 CONTRIBUTION FROM CLEAR-SKY FRACTION
C
      DO 312 JL = 1, KDLON
      ZFD (JL) = (1. - ZCLM(JL,JK1,KFLEV)) * ZDNF(JL,1,JK1)
 312  CONTINUE
C
C*                 CONTRIBUTION FROM ADJACENT CLOUD
C
      DO 313 JL = 1, KDLON
      ZFD(JL) = ZFD(JL) + ZCLM(JL,JK1,JK1) * ZDNF(JL,JK1+1,JK1)
 313  CONTINUE
C
C*                 CONTRIBUTION FROM OTHER CLOUDY FRACTIONS
C
      DO 315 JK = KFLEV-1 , JK1 , -1
      DO 314 JL = 1, KDLON
      ZCFRAC = ZCLM(JL,JK1,JK+1) - ZCLM(JL,JK1,JK)
      ZFD(JL) =  ZFD(JL) + ZCFRAC * ZDNF(JL,JK+2,JK1)
 314  CONTINUE
 315  CONTINUE
C
      DO 316 JL = 1, KDLON
      PFLUX(JL,2,JK1) = ZFD (JL)
 316  CONTINUE
C
 317  CONTINUE
C
C
C
C
C*         3.2     UPWARD FLUX AT THE SURFACE
C                  --------------------------
C
 320  CONTINUE
C
      DO 321 JL = 1, KDLON
      PFLUX(JL,1,1) = PEMIS(JL)*PBSUIN(JL)-(1.-PEMIS(JL))*PFLUX(JL,2,1)
 321  CONTINUE
C
C
C
C*         3.3     UPWARD FLUXES
C                  -------------
C
 330  CONTINUE
C
      DO 337 JK1 = 2 , KFLEV+1
C
C*                 CONTRIBUTION FROM CLEAR-SKY FRACTION
C
      DO 332 JL = 1, KDLON
      ZFU (JL) = (1. - ZCLM(JL,JK1,1)) * ZUPF(JL,1,JK1)
 332  CONTINUE
C
C*                 CONTRIBUTION FROM ADJACENT CLOUD
C
      DO 333 JL = 1, KDLON
      ZFU(JL) =  ZFU(JL) + ZCLM(JL,JK1,JK1-1) * ZUPF(JL,JK1,JK1)
 333  CONTINUE
C
C*                 CONTRIBUTION FROM OTHER CLOUDY FRACTIONS
C
      DO 335 JK = 2 , JK1-1
      DO 334 JL = 1, KDLON
      ZCFRAC = ZCLM(JL,JK1,JK-1) - ZCLM(JL,JK1,JK)
      ZFU(JL) =  ZFU(JL) + ZCFRAC * ZUPF(JL,JK  ,JK1)
 334  CONTINUE
 335  CONTINUE
C
      DO 336 JL = 1, KDLON
      PFLUX(JL,1,JK1) = ZFU (JL)
 336  CONTINUE
C
 337  CONTINUE
C
C
      END IF
C
C
C*         2.3     END OF CLOUD EFFECT COMPUTATIONS
C
 230  CONTINUE
C
      IF (.NOT.LEVOIGT) THEN
        DO 231 JL = 1, KDLON
        ZFN10(JL) = PFLUX(JL,1,KLIM) + PFLUX(JL,2,KLIM)
 231    CONTINUE
        DO 233 JK = KLIM+1 , KFLEV+1
        DO 232 JL = 1, KDLON
        ZFN10(JL) = ZFN10(JL) + PCTS(JL,JK-1)
        PFLUX(JL,1,JK) = ZFN10(JL)
        PFLUX(JL,2,JK) = 0.0
 232    CONTINUE
 233    CONTINUE
      ENDIF
C
      RETURN
      END
      SUBROUTINE LWB(PDT0,PTAVE,PTL
     S  , PB,PBINT,PBSUIN,PBSUR,PBTOP,PDBSL
     S  , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP)
      use dimens_m
      use dimphy
      use raddim
      IMPLICIT none
      include "raddimlw.h"
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           COMPUTES PLANCK FUNCTIONS
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C     ==== INPUTS ===
C PDT0   : (KDLON)             ; SURFACE TEMPERATURE DISCONTINUITY
C PTAVE  : (KDLON,KFLEV)       ; TEMPERATURE
C PTL    : (KDLON,0:KFLEV)     ; HALF LEVEL TEMPERATURE
C     ==== OUTPUTS ===
C PB     : (KDLON,Ninter,KFLEV+1); SPECTRAL HALF LEVEL PLANCK FUNCTION
C PBINT  : (KDLON,KFLEV+1)     ; HALF LEVEL PLANCK FUNCTION
C PBSUIN : (KDLON)             ; SURFACE PLANCK FUNCTION
C PBSUR  : (KDLON,Ninter)        ; SURFACE SPECTRAL PLANCK FUNCTION
C PBTOP  : (KDLON,Ninter)        ; TOP SPECTRAL PLANCK FUNCTION
C PDBSL  : (KDLON,Ninter,KFLEV*2); SUB-LAYER PLANCK FUNCTION GRADIENT
C PGA    : (KDLON,8,2,KFLEV); dB/dT-weighted LAYER PADE APPROXIMANTS
C PGB    : (KDLON,8,2,KFLEV); dB/dT-weighted LAYER PADE APPROXIMANTS
C PGASUR, PGBSUR (KDLON,8,2)   ; SURFACE PADE APPROXIMANTS
C PGATOP, PGBTOP (KDLON,8,2)   ; T.O.A. PADE APPROXIMANTS
C
C        IMPLICIT ARGUMENTS :   NONE
C        --------------------
C
C     METHOD.
C     -------
C
C          1. COMPUTES THE PLANCK FUNCTION ON ALL LEVELS AND HALF LEVELS
C     FROM A POLYNOMIAL DEVELOPMENT OF PLANCK FUNCTION
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS           "
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C
C-----------------------------------------------------------------------
C
C ARGUMENTS:
C
      REAL*8 PDT0(KDLON)
      REAL*8 PTAVE(KDLON,KFLEV)
      REAL*8 PTL(KDLON,KFLEV+1)
C
      REAL*8 PB(KDLON,Ninter,KFLEV+1) ! SPECTRAL HALF LEVEL PLANCK FUNCTION
      REAL*8 PBINT(KDLON,KFLEV+1) ! HALF LEVEL PLANCK FUNCTION
      REAL*8 PBSUIN(KDLON) ! SURFACE PLANCK FUNCTION
      REAL*8 PBSUR(KDLON,Ninter) ! SURFACE SPECTRAL PLANCK FUNCTION
      REAL*8 PBTOP(KDLON,Ninter) ! TOP SPECTRAL PLANCK FUNCTION
      REAL*8 PDBSL(KDLON,Ninter,KFLEV*2) ! SUB-LAYER PLANCK FUNCTION GRADIENT
      REAL*8 PGA(KDLON,8,2,KFLEV) ! dB/dT-weighted LAYER PADE APPROXIMANTS
      REAL*8 PGB(KDLON,8,2,KFLEV) ! dB/dT-weighted LAYER PADE APPROXIMANTS
      REAL*8 PGASUR(KDLON,8,2) ! SURFACE PADE APPROXIMANTS
      REAL*8 PGBSUR(KDLON,8,2) ! SURFACE PADE APPROXIMANTS
      REAL*8 PGATOP(KDLON,8,2) ! T.O.A. PADE APPROXIMANTS
      REAL*8 PGBTOP(KDLON,8,2) ! T.O.A. PADE APPROXIMANTS
C
C-------------------------------------------------------------------------
C*  LOCAL VARIABLES:
      INTEGER INDB(KDLON),INDS(KDLON)
      REAL*8 ZBLAY(KDLON,KFLEV),ZBLEV(KDLON,KFLEV+1)
      REAL*8 ZRES(KDLON),ZRES2(KDLON),ZTI(KDLON),ZTI2(KDLON)
c
      INTEGER jk, jl, ic, jnu, jf, jg
      INTEGER jk1, jk2
      INTEGER k, j, ixtox, indto, ixtx, indt
      INTEGER indsu, indtp
      REAL*8 zdsto1, zdstox, zdst1, zdstx
c
C* Quelques parametres:
      REAL*8 TSTAND
      PARAMETER (TSTAND=250.0)
      REAL*8 TSTP
      PARAMETER (TSTP=12.5)
      INTEGER MXIXT
      PARAMETER (MXIXT=10)
C
C* Used Data Block:
      REAL*8 TINTP(11)
      SAVE TINTP
      REAL*8 GA(11,16,3), GB(11,16,3)
      SAVE GA, GB
      REAL*8 XP(6,6)
      SAVE XP
c
      DATA TINTP / 187.5, 200., 212.5, 225., 237.5, 250.,
     S             262.5, 275., 287.5, 300., 312.5 /
C-----------------------------------------------------------------------
C-- WATER VAPOR -- INT.1 -- 0- 500 CM-1 -- FROM ABS225 ----------------
C
C
C
C
C-- R.D. -- G = - 0.2 SLA
C
C
C----- INTERVAL = 1 ----- T =  187.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 1, 1,IC),IC=1,3) /
     S 0.63499072E-02,-0.99506586E-03, 0.00000000E+00/
      DATA (GB( 1, 1,IC),IC=1,3) /
     S 0.63499072E-02, 0.97222852E-01, 0.10000000E+01/
      DATA (GA( 1, 2,IC),IC=1,3) /
     S 0.77266491E-02,-0.11661515E-02, 0.00000000E+00/
      DATA (GB( 1, 2,IC),IC=1,3) /
     S 0.77266491E-02, 0.10681591E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  200.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 2, 1,IC),IC=1,3) /
     S 0.65566348E-02,-0.10184169E-02, 0.00000000E+00/
      DATA (GB( 2, 1,IC),IC=1,3) /
     S 0.65566348E-02, 0.98862238E-01, 0.10000000E+01/
      DATA (GA( 2, 2,IC),IC=1,3) /
     S 0.81323287E-02,-0.11886130E-02, 0.00000000E+00/
      DATA (GB( 2, 2,IC),IC=1,3) /
     S 0.81323287E-02, 0.10921298E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  212.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 3, 1,IC),IC=1,3) /
     S 0.67849730E-02,-0.10404730E-02, 0.00000000E+00/
      DATA (GB( 3, 1,IC),IC=1,3) /
     S 0.67849730E-02, 0.10061504E+00, 0.10000000E+01/
      DATA (GA( 3, 2,IC),IC=1,3) /
     S 0.86507620E-02,-0.12139929E-02, 0.00000000E+00/
      DATA (GB( 3, 2,IC),IC=1,3) /
     S 0.86507620E-02, 0.11198225E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  225.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 4, 1,IC),IC=1,3) /
     S 0.70481947E-02,-0.10621792E-02, 0.00000000E+00/
      DATA (GB( 4, 1,IC),IC=1,3) /
     S 0.70481947E-02, 0.10256222E+00, 0.10000000E+01/
      DATA (GA( 4, 2,IC),IC=1,3) /
     S 0.92776391E-02,-0.12445811E-02, 0.00000000E+00/
      DATA (GB( 4, 2,IC),IC=1,3) /
     S 0.92776391E-02, 0.11487826E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  237.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 5, 1,IC),IC=1,3) /
     S 0.73585943E-02,-0.10847662E-02, 0.00000000E+00/
      DATA (GB( 5, 1,IC),IC=1,3) /
     S 0.73585943E-02, 0.10475952E+00, 0.10000000E+01/
      DATA (GA( 5, 2,IC),IC=1,3) /
     S 0.99806312E-02,-0.12807672E-02, 0.00000000E+00/
      DATA (GB( 5, 2,IC),IC=1,3) /
     S 0.99806312E-02, 0.11751113E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  250.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 6, 1,IC),IC=1,3) /
     S 0.77242818E-02,-0.11094726E-02, 0.00000000E+00/
      DATA (GB( 6, 1,IC),IC=1,3) /
     S 0.77242818E-02, 0.10720986E+00, 0.10000000E+01/
      DATA (GA( 6, 2,IC),IC=1,3) /
     S 0.10709803E-01,-0.13208251E-02, 0.00000000E+00/
      DATA (GB( 6, 2,IC),IC=1,3) /
     S 0.10709803E-01, 0.11951535E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  262.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 7, 1,IC),IC=1,3) /
     S 0.81472693E-02,-0.11372949E-02, 0.00000000E+00/
      DATA (GB( 7, 1,IC),IC=1,3) /
     S 0.81472693E-02, 0.10985370E+00, 0.10000000E+01/
      DATA (GA( 7, 2,IC),IC=1,3) /
     S 0.11414739E-01,-0.13619034E-02, 0.00000000E+00/
      DATA (GB( 7, 2,IC),IC=1,3) /
     S 0.11414739E-01, 0.12069945E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  275.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 8, 1,IC),IC=1,3) /
     S 0.86227527E-02,-0.11687683E-02, 0.00000000E+00/
      DATA (GB( 8, 1,IC),IC=1,3) /
     S 0.86227527E-02, 0.11257633E+00, 0.10000000E+01/
      DATA (GA( 8, 2,IC),IC=1,3) /
     S 0.12058772E-01,-0.14014165E-02, 0.00000000E+00/
      DATA (GB( 8, 2,IC),IC=1,3) /
     S 0.12058772E-01, 0.12108524E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  287.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 9, 1,IC),IC=1,3) /
     S 0.91396814E-02,-0.12038314E-02, 0.00000000E+00/
      DATA (GB( 9, 1,IC),IC=1,3) /
     S 0.91396814E-02, 0.11522980E+00, 0.10000000E+01/
      DATA (GA( 9, 2,IC),IC=1,3) /
     S 0.12623992E-01,-0.14378639E-02, 0.00000000E+00/
      DATA (GB( 9, 2,IC),IC=1,3) /
     S 0.12623992E-01, 0.12084229E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  300.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA(10, 1,IC),IC=1,3) /
     S 0.96825438E-02,-0.12418367E-02, 0.00000000E+00/
      DATA (GB(10, 1,IC),IC=1,3) /
     S 0.96825438E-02, 0.11766343E+00, 0.10000000E+01/
      DATA (GA(10, 2,IC),IC=1,3) /
     S 0.13108146E-01,-0.14708488E-02, 0.00000000E+00/
      DATA (GB(10, 2,IC),IC=1,3) /
     S 0.13108146E-01, 0.12019005E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  312.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA(11, 1,IC),IC=1,3) /
     S 0.10233955E-01,-0.12817135E-02, 0.00000000E+00/
      DATA (GB(11, 1,IC),IC=1,3) /
     S 0.10233955E-01, 0.11975320E+00, 0.10000000E+01/
      DATA (GA(11, 2,IC),IC=1,3) /
     S 0.13518390E-01,-0.15006791E-02, 0.00000000E+00/
      DATA (GB(11, 2,IC),IC=1,3) /
     S 0.13518390E-01, 0.11932684E+00, 0.10000000E+01/
C
C
C
C--- WATER VAPOR --- INTERVAL 2 -- 500-800 CM-1--- FROM ABS225 ---------
C
C
C
C
C--- R.D.  ---  G = 0.02 + 0.50 / ( 1 + 4.5 U )
C
C
C----- INTERVAL = 2 ----- T =  187.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 1, 3,IC),IC=1,3) /
     S 0.11644593E+01, 0.41243390E+00, 0.00000000E+00/
      DATA (GB( 1, 3,IC),IC=1,3) /
     S 0.11644593E+01, 0.10346097E+01, 0.10000000E+01/
      DATA (GA( 1, 4,IC),IC=1,3) /
     S 0.12006968E+01, 0.48318936E+00, 0.00000000E+00/
      DATA (GB( 1, 4,IC),IC=1,3) /
     S 0.12006968E+01, 0.10626130E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  200.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 2, 3,IC),IC=1,3) /
     S 0.11747203E+01, 0.43407282E+00, 0.00000000E+00/
      DATA (GB( 2, 3,IC),IC=1,3) /
     S 0.11747203E+01, 0.10433655E+01, 0.10000000E+01/
      DATA (GA( 2, 4,IC),IC=1,3) /
     S 0.12108196E+01, 0.50501827E+00, 0.00000000E+00/
      DATA (GB( 2, 4,IC),IC=1,3) /
     S 0.12108196E+01, 0.10716026E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  212.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 3, 3,IC),IC=1,3) /
     S 0.11837872E+01, 0.45331413E+00, 0.00000000E+00/
      DATA (GB( 3, 3,IC),IC=1,3) /
     S 0.11837872E+01, 0.10511933E+01, 0.10000000E+01/
      DATA (GA( 3, 4,IC),IC=1,3) /
     S 0.12196717E+01, 0.52409502E+00, 0.00000000E+00/
      DATA (GB( 3, 4,IC),IC=1,3) /
     S 0.12196717E+01, 0.10795108E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  225.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 4, 3,IC),IC=1,3) /
     S 0.11918561E+01, 0.47048604E+00, 0.00000000E+00/
      DATA (GB( 4, 3,IC),IC=1,3) /
     S 0.11918561E+01, 0.10582150E+01, 0.10000000E+01/
      DATA (GA( 4, 4,IC),IC=1,3) /
     S 0.12274493E+01, 0.54085277E+00, 0.00000000E+00/
      DATA (GB( 4, 4,IC),IC=1,3) /
     S 0.12274493E+01, 0.10865006E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  237.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 5, 3,IC),IC=1,3) /
     S 0.11990757E+01, 0.48586286E+00, 0.00000000E+00/
      DATA (GB( 5, 3,IC),IC=1,3) /
     S 0.11990757E+01, 0.10645317E+01, 0.10000000E+01/
      DATA (GA( 5, 4,IC),IC=1,3) /
     S 0.12343189E+01, 0.55565422E+00, 0.00000000E+00/
      DATA (GB( 5, 4,IC),IC=1,3) /
     S 0.12343189E+01, 0.10927103E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  250.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 6, 3,IC),IC=1,3) /
     S 0.12055643E+01, 0.49968044E+00, 0.00000000E+00/
      DATA (GB( 6, 3,IC),IC=1,3) /
     S 0.12055643E+01, 0.10702313E+01, 0.10000000E+01/
      DATA (GA( 6, 4,IC),IC=1,3) /
     S 0.12404147E+01, 0.56878618E+00, 0.00000000E+00/
      DATA (GB( 6, 4,IC),IC=1,3) /
     S 0.12404147E+01, 0.10982489E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  262.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 7, 3,IC),IC=1,3) /
     S 0.12114186E+01, 0.51214132E+00, 0.00000000E+00/
      DATA (GB( 7, 3,IC),IC=1,3) /
     S 0.12114186E+01, 0.10753907E+01, 0.10000000E+01/
      DATA (GA( 7, 4,IC),IC=1,3) /
     S 0.12458431E+01, 0.58047395E+00, 0.00000000E+00/
      DATA (GB( 7, 4,IC),IC=1,3) /
     S 0.12458431E+01, 0.11032019E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  275.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 8, 3,IC),IC=1,3) /
     S 0.12167192E+01, 0.52341830E+00, 0.00000000E+00/
      DATA (GB( 8, 3,IC),IC=1,3) /
     S 0.12167192E+01, 0.10800762E+01, 0.10000000E+01/
      DATA (GA( 8, 4,IC),IC=1,3) /
     S 0.12506907E+01, 0.59089894E+00, 0.00000000E+00/
      DATA (GB( 8, 4,IC),IC=1,3) /
     S 0.12506907E+01, 0.11076379E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  287.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 9, 3,IC),IC=1,3) /
     S 0.12215344E+01, 0.53365803E+00, 0.00000000E+00/
      DATA (GB( 9, 3,IC),IC=1,3) /
     S 0.12215344E+01, 0.10843446E+01, 0.10000000E+01/
      DATA (GA( 9, 4,IC),IC=1,3) /
     S 0.12550299E+01, 0.60021475E+00, 0.00000000E+00/
      DATA (GB( 9, 4,IC),IC=1,3) /
     S 0.12550299E+01, 0.11116160E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  300.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA(10, 3,IC),IC=1,3) /
     S 0.12259226E+01, 0.54298448E+00, 0.00000000E+00/
      DATA (GB(10, 3,IC),IC=1,3) /
     S 0.12259226E+01, 0.10882439E+01, 0.10000000E+01/
      DATA (GA(10, 4,IC),IC=1,3) /
     S 0.12589256E+01, 0.60856112E+00, 0.00000000E+00/
      DATA (GB(10, 4,IC),IC=1,3) /
     S 0.12589256E+01, 0.11151910E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  312.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA(11, 3,IC),IC=1,3) /
     S 0.12299344E+01, 0.55150227E+00, 0.00000000E+00/
      DATA (GB(11, 3,IC),IC=1,3) /
     S 0.12299344E+01, 0.10918144E+01, 0.10000000E+01/
      DATA (GA(11, 4,IC),IC=1,3) /
     S 0.12624402E+01, 0.61607594E+00, 0.00000000E+00/
      DATA (GB(11, 4,IC),IC=1,3) /
     S 0.12624402E+01, 0.11184188E+01, 0.10000000E+01/
C
C
C
C
C
C
C- WATER VAPOR - INT. 3 -- 800-970 + 1110-1250 CM-1 -- FIT FROM 215 IS -
C
C
C-- WATER VAPOR LINES IN THE WINDOW REGION (800-1250 CM-1)
C
C
C
C--- G = 3.875E-03 ---------------
C
C----- INTERVAL = 3 ----- T =  187.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 1, 7,IC),IC=1,3) /
     S 0.10192131E+02, 0.80737799E+01, 0.00000000E+00/
      DATA (GB( 1, 7,IC),IC=1,3) /
     S 0.10192131E+02, 0.82623280E+01, 0.10000000E+01/
      DATA (GA( 1, 8,IC),IC=1,3) /
     S 0.92439050E+01, 0.77425778E+01, 0.00000000E+00/
      DATA (GB( 1, 8,IC),IC=1,3) /
     S 0.92439050E+01, 0.79342219E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  200.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 2, 7,IC),IC=1,3) /
     S 0.97258602E+01, 0.79171158E+01, 0.00000000E+00/
      DATA (GB( 2, 7,IC),IC=1,3) /
     S 0.97258602E+01, 0.81072291E+01, 0.10000000E+01/
      DATA (GA( 2, 8,IC),IC=1,3) /
     S 0.87567422E+01, 0.75443460E+01, 0.00000000E+00/
      DATA (GB( 2, 8,IC),IC=1,3) /
     S 0.87567422E+01, 0.77373458E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  212.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 3, 7,IC),IC=1,3) /
     S 0.92992890E+01, 0.77609605E+01, 0.00000000E+00/
      DATA (GB( 3, 7,IC),IC=1,3) /
     S 0.92992890E+01, 0.79523834E+01, 0.10000000E+01/
      DATA (GA( 3, 8,IC),IC=1,3) /
     S 0.83270144E+01, 0.73526151E+01, 0.00000000E+00/
      DATA (GB( 3, 8,IC),IC=1,3) /
     S 0.83270144E+01, 0.75467334E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  225.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 4, 7,IC),IC=1,3) /
     S 0.89154021E+01, 0.76087371E+01, 0.00000000E+00/
      DATA (GB( 4, 7,IC),IC=1,3) /
     S 0.89154021E+01, 0.78012527E+01, 0.10000000E+01/
      DATA (GA( 4, 8,IC),IC=1,3) /
     S 0.79528337E+01, 0.71711188E+01, 0.00000000E+00/
      DATA (GB( 4, 8,IC),IC=1,3) /
     S 0.79528337E+01, 0.73661786E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  237.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 5, 7,IC),IC=1,3) /
     S 0.85730084E+01, 0.74627112E+01, 0.00000000E+00/
      DATA (GB( 5, 7,IC),IC=1,3) /
     S 0.85730084E+01, 0.76561458E+01, 0.10000000E+01/
      DATA (GA( 5, 8,IC),IC=1,3) /
     S 0.76286839E+01, 0.70015571E+01, 0.00000000E+00/
      DATA (GB( 5, 8,IC),IC=1,3) /
     S 0.76286839E+01, 0.71974319E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  250.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 6, 7,IC),IC=1,3) /
     S 0.82685838E+01, 0.73239981E+01, 0.00000000E+00/
      DATA (GB( 6, 7,IC),IC=1,3) /
     S 0.82685838E+01, 0.75182174E+01, 0.10000000E+01/
      DATA (GA( 6, 8,IC),IC=1,3) /
     S 0.73477879E+01, 0.68442532E+01, 0.00000000E+00/
      DATA (GB( 6, 8,IC),IC=1,3) /
     S 0.73477879E+01, 0.70408543E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  262.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 7, 7,IC),IC=1,3) /
     S 0.79978921E+01, 0.71929934E+01, 0.00000000E+00/
      DATA (GB( 7, 7,IC),IC=1,3) /
     S 0.79978921E+01, 0.73878952E+01, 0.10000000E+01/
      DATA (GA( 7, 8,IC),IC=1,3) /
     S 0.71035818E+01, 0.66987996E+01, 0.00000000E+00/
      DATA (GB( 7, 8,IC),IC=1,3) /
     S 0.71035818E+01, 0.68960649E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  275.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 8, 7,IC),IC=1,3) /
     S 0.77568055E+01, 0.70697065E+01, 0.00000000E+00/
      DATA (GB( 8, 7,IC),IC=1,3) /
     S 0.77568055E+01, 0.72652133E+01, 0.10000000E+01/
      DATA (GA( 8, 8,IC),IC=1,3) /
     S 0.68903312E+01, 0.65644820E+01, 0.00000000E+00/
      DATA (GB( 8, 8,IC),IC=1,3) /
     S 0.68903312E+01, 0.67623672E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  287.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 9, 7,IC),IC=1,3) /
     S 0.75416266E+01, 0.69539626E+01, 0.00000000E+00/
      DATA (GB( 9, 7,IC),IC=1,3) /
     S 0.75416266E+01, 0.71500151E+01, 0.10000000E+01/
      DATA (GA( 9, 8,IC),IC=1,3) /
     S 0.67032875E+01, 0.64405267E+01, 0.00000000E+00/
      DATA (GB( 9, 8,IC),IC=1,3) /
     S 0.67032875E+01, 0.66389989E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  300.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA(10, 7,IC),IC=1,3) /
     S 0.73491694E+01, 0.68455144E+01, 0.00000000E+00/
      DATA (GB(10, 7,IC),IC=1,3) /
     S 0.73491694E+01, 0.70420667E+01, 0.10000000E+01/
      DATA (GA(10, 8,IC),IC=1,3) /
     S 0.65386461E+01, 0.63262376E+01, 0.00000000E+00/
      DATA (GB(10, 8,IC),IC=1,3) /
     S 0.65386461E+01, 0.65252707E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  312.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA(11, 7,IC),IC=1,3) /
     S 0.71767400E+01, 0.67441020E+01, 0.00000000E+00/
      DATA (GB(11, 7,IC),IC=1,3) /
     S 0.71767400E+01, 0.69411177E+01, 0.10000000E+01/
      DATA (GA(11, 8,IC),IC=1,3) /
     S 0.63934377E+01, 0.62210701E+01, 0.00000000E+00/
      DATA (GB(11, 8,IC),IC=1,3) /
     S 0.63934377E+01, 0.64206412E+01, 0.10000000E+01/
C
C
C-- WATER VAPOR -- 970-1110 CM-1 ----------------------------------------
C
C-- G = 3.6E-03
C
C----- INTERVAL = 4 ----- T =  187.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 1, 9,IC),IC=1,3) /
     S 0.24870635E+02, 0.10542131E+02, 0.00000000E+00/
      DATA (GB( 1, 9,IC),IC=1,3) /
     S 0.24870635E+02, 0.10656640E+02, 0.10000000E+01/
      DATA (GA( 1,10,IC),IC=1,3) /
     S 0.24586283E+02, 0.10490353E+02, 0.00000000E+00/
      DATA (GB( 1,10,IC),IC=1,3) /
     S 0.24586283E+02, 0.10605856E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  200.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 2, 9,IC),IC=1,3) /
     S 0.24725591E+02, 0.10515895E+02, 0.00000000E+00/
      DATA (GB( 2, 9,IC),IC=1,3) /
     S 0.24725591E+02, 0.10630910E+02, 0.10000000E+01/
      DATA (GA( 2,10,IC),IC=1,3) /
     S 0.24441465E+02, 0.10463512E+02, 0.00000000E+00/
      DATA (GB( 2,10,IC),IC=1,3) /
     S 0.24441465E+02, 0.10579514E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  212.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 3, 9,IC),IC=1,3) /
     S 0.24600320E+02, 0.10492949E+02, 0.00000000E+00/
      DATA (GB( 3, 9,IC),IC=1,3) /
     S 0.24600320E+02, 0.10608399E+02, 0.10000000E+01/
      DATA (GA( 3,10,IC),IC=1,3) /
     S 0.24311657E+02, 0.10439183E+02, 0.00000000E+00/
      DATA (GB( 3,10,IC),IC=1,3) /
     S 0.24311657E+02, 0.10555632E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  225.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 4, 9,IC),IC=1,3) /
     S 0.24487300E+02, 0.10472049E+02, 0.00000000E+00/
      DATA (GB( 4, 9,IC),IC=1,3) /
     S 0.24487300E+02, 0.10587891E+02, 0.10000000E+01/
      DATA (GA( 4,10,IC),IC=1,3) /
     S 0.24196167E+02, 0.10417324E+02, 0.00000000E+00/
      DATA (GB( 4,10,IC),IC=1,3) /
     S 0.24196167E+02, 0.10534169E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  237.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 5, 9,IC),IC=1,3) /
     S 0.24384935E+02, 0.10452961E+02, 0.00000000E+00/
      DATA (GB( 5, 9,IC),IC=1,3) /
     S 0.24384935E+02, 0.10569156E+02, 0.10000000E+01/
      DATA (GA( 5,10,IC),IC=1,3) /
     S 0.24093406E+02, 0.10397704E+02, 0.00000000E+00/
      DATA (GB( 5,10,IC),IC=1,3) /
     S 0.24093406E+02, 0.10514900E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  250.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 6, 9,IC),IC=1,3) /
     S 0.24292341E+02, 0.10435562E+02, 0.00000000E+00/
      DATA (GB( 6, 9,IC),IC=1,3) /
     S 0.24292341E+02, 0.10552075E+02, 0.10000000E+01/
      DATA (GA( 6,10,IC),IC=1,3) /
     S 0.24001597E+02, 0.10380038E+02, 0.00000000E+00/
      DATA (GB( 6,10,IC),IC=1,3) /
     S 0.24001597E+02, 0.10497547E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  262.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 7, 9,IC),IC=1,3) /
     S 0.24208572E+02, 0.10419710E+02, 0.00000000E+00/
      DATA (GB( 7, 9,IC),IC=1,3) /
     S 0.24208572E+02, 0.10536510E+02, 0.10000000E+01/
      DATA (GA( 7,10,IC),IC=1,3) /
     S 0.23919098E+02, 0.10364052E+02, 0.00000000E+00/
      DATA (GB( 7,10,IC),IC=1,3) /
     S 0.23919098E+02, 0.10481842E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  275.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 8, 9,IC),IC=1,3) /
     S 0.24132642E+02, 0.10405247E+02, 0.00000000E+00/
      DATA (GB( 8, 9,IC),IC=1,3) /
     S 0.24132642E+02, 0.10522307E+02, 0.10000000E+01/
      DATA (GA( 8,10,IC),IC=1,3) /
     S 0.23844511E+02, 0.10349509E+02, 0.00000000E+00/
      DATA (GB( 8,10,IC),IC=1,3) /
     S 0.23844511E+02, 0.10467553E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  287.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 9, 9,IC),IC=1,3) /
     S 0.24063614E+02, 0.10392022E+02, 0.00000000E+00/
      DATA (GB( 9, 9,IC),IC=1,3) /
     S 0.24063614E+02, 0.10509317E+02, 0.10000000E+01/
      DATA (GA( 9,10,IC),IC=1,3) /
     S 0.23776708E+02, 0.10336215E+02, 0.00000000E+00/
      DATA (GB( 9,10,IC),IC=1,3) /
     S 0.23776708E+02, 0.10454488E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  300.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA(10, 9,IC),IC=1,3) /
     S 0.24000649E+02, 0.10379892E+02, 0.00000000E+00/
      DATA (GB(10, 9,IC),IC=1,3) /
     S 0.24000649E+02, 0.10497402E+02, 0.10000000E+01/
      DATA (GA(10,10,IC),IC=1,3) /
     S 0.23714816E+02, 0.10324018E+02, 0.00000000E+00/
      DATA (GB(10,10,IC),IC=1,3) /
     S 0.23714816E+02, 0.10442501E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  312.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA(11, 9,IC),IC=1,3) /
     S 0.23943021E+02, 0.10368736E+02, 0.00000000E+00/
      DATA (GB(11, 9,IC),IC=1,3) /
     S 0.23943021E+02, 0.10486443E+02, 0.10000000E+01/
      DATA (GA(11,10,IC),IC=1,3) /
     S 0.23658197E+02, 0.10312808E+02, 0.00000000E+00/
      DATA (GB(11,10,IC),IC=1,3) /
     S 0.23658197E+02, 0.10431483E+02, 0.10000000E+01/
C
C
C
C-- H2O -- WEAKER PARTS OF THE STRONG BANDS  -- FROM ABS225 ----
C
C-- WATER VAPOR --- 350 - 500 CM-1
C
C-- G = - 0.2*SLA, 0.0 +0.5/(1+0.5U)
C
C----- INTERVAL = 5 ----- T =  187.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 1, 5,IC),IC=1,3) /
     S 0.15750172E+00,-0.22159303E-01, 0.00000000E+00/
      DATA (GB( 1, 5,IC),IC=1,3) /
     S 0.15750172E+00, 0.38103212E+00, 0.10000000E+01/
      DATA (GA( 1, 6,IC),IC=1,3) /
     S 0.17770551E+00,-0.24972399E-01, 0.00000000E+00/
      DATA (GB( 1, 6,IC),IC=1,3) /
     S 0.17770551E+00, 0.41646579E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  200.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 2, 5,IC),IC=1,3) /
     S 0.16174076E+00,-0.22748917E-01, 0.00000000E+00/
      DATA (GB( 2, 5,IC),IC=1,3) /
     S 0.16174076E+00, 0.38913800E+00, 0.10000000E+01/
      DATA (GA( 2, 6,IC),IC=1,3) /
     S 0.18176757E+00,-0.25537247E-01, 0.00000000E+00/
      DATA (GB( 2, 6,IC),IC=1,3) /
     S 0.18176757E+00, 0.42345095E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  212.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 3, 5,IC),IC=1,3) /
     S 0.16548628E+00,-0.23269898E-01, 0.00000000E+00/
      DATA (GB( 3, 5,IC),IC=1,3) /
     S 0.16548628E+00, 0.39613651E+00, 0.10000000E+01/
      DATA (GA( 3, 6,IC),IC=1,3) /
     S 0.18527967E+00,-0.26025624E-01, 0.00000000E+00/
      DATA (GB( 3, 6,IC),IC=1,3) /
     S 0.18527967E+00, 0.42937476E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  225.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 4, 5,IC),IC=1,3) /
     S 0.16881124E+00,-0.23732392E-01, 0.00000000E+00/
      DATA (GB( 4, 5,IC),IC=1,3) /
     S 0.16881124E+00, 0.40222421E+00, 0.10000000E+01/
      DATA (GA( 4, 6,IC),IC=1,3) /
     S 0.18833348E+00,-0.26450280E-01, 0.00000000E+00/
      DATA (GB( 4, 6,IC),IC=1,3) /
     S 0.18833348E+00, 0.43444062E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  237.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 5, 5,IC),IC=1,3) /
     S 0.17177839E+00,-0.24145123E-01, 0.00000000E+00/
      DATA (GB( 5, 5,IC),IC=1,3) /
     S 0.17177839E+00, 0.40756010E+00, 0.10000000E+01/
      DATA (GA( 5, 6,IC),IC=1,3) /
     S 0.19100108E+00,-0.26821236E-01, 0.00000000E+00/
      DATA (GB( 5, 6,IC),IC=1,3) /
     S 0.19100108E+00, 0.43880316E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  250.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 6, 5,IC),IC=1,3) /
     S 0.17443933E+00,-0.24515269E-01, 0.00000000E+00/
      DATA (GB( 6, 5,IC),IC=1,3) /
     S 0.17443933E+00, 0.41226954E+00, 0.10000000E+01/
      DATA (GA( 6, 6,IC),IC=1,3) /
     S 0.19334122E+00,-0.27146657E-01, 0.00000000E+00/
      DATA (GB( 6, 6,IC),IC=1,3) /
     S 0.19334122E+00, 0.44258354E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  262.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 7, 5,IC),IC=1,3) /
     S 0.17683622E+00,-0.24848690E-01, 0.00000000E+00/
      DATA (GB( 7, 5,IC),IC=1,3) /
     S 0.17683622E+00, 0.41645142E+00, 0.10000000E+01/
      DATA (GA( 7, 6,IC),IC=1,3) /
     S 0.19540288E+00,-0.27433354E-01, 0.00000000E+00/
      DATA (GB( 7, 6,IC),IC=1,3) /
     S 0.19540288E+00, 0.44587882E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  275.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 8, 5,IC),IC=1,3) /
     S 0.17900375E+00,-0.25150210E-01, 0.00000000E+00/
      DATA (GB( 8, 5,IC),IC=1,3) /
     S 0.17900375E+00, 0.42018474E+00, 0.10000000E+01/
      DATA (GA( 8, 6,IC),IC=1,3) /
     S 0.19722732E+00,-0.27687065E-01, 0.00000000E+00/
      DATA (GB( 8, 6,IC),IC=1,3) /
     S 0.19722732E+00, 0.44876776E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  287.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 9, 5,IC),IC=1,3) /
     S 0.18097099E+00,-0.25423873E-01, 0.00000000E+00/
      DATA (GB( 9, 5,IC),IC=1,3) /
     S 0.18097099E+00, 0.42353379E+00, 0.10000000E+01/
      DATA (GA( 9, 6,IC),IC=1,3) /
     S 0.19884918E+00,-0.27912608E-01, 0.00000000E+00/
      DATA (GB( 9, 6,IC),IC=1,3) /
     S 0.19884918E+00, 0.45131451E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  300.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA(10, 5,IC),IC=1,3) /
     S 0.18276283E+00,-0.25673139E-01, 0.00000000E+00/
      DATA (GB(10, 5,IC),IC=1,3) /
     S 0.18276283E+00, 0.42655211E+00, 0.10000000E+01/
      DATA (GA(10, 6,IC),IC=1,3) /
     S 0.20029696E+00,-0.28113944E-01, 0.00000000E+00/
      DATA (GB(10, 6,IC),IC=1,3) /
     S 0.20029696E+00, 0.45357095E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  312.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA(11, 5,IC),IC=1,3) /
     S 0.18440117E+00,-0.25901055E-01, 0.00000000E+00/
      DATA (GB(11, 5,IC),IC=1,3) /
     S 0.18440117E+00, 0.42928533E+00, 0.10000000E+01/
      DATA (GA(11, 6,IC),IC=1,3) /
     S 0.20159300E+00,-0.28294180E-01, 0.00000000E+00/
      DATA (GB(11, 6,IC),IC=1,3) /
     S 0.20159300E+00, 0.45557797E+00, 0.10000000E+01/
C
C
C
C
C- WATER VAPOR - WINGS OF VIBRATION-ROTATION BAND - 1250-1450+1880-2820 -
C--- G = 0.0
C
C
C----- INTERVAL = 6 ----- T =  187.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 1,11,IC),IC=1,3) /
     S 0.11990218E+02,-0.12823142E+01, 0.00000000E+00/
      DATA (GB( 1,11,IC),IC=1,3) /
     S 0.11990218E+02, 0.26681588E+02, 0.10000000E+01/
      DATA (GA( 1,12,IC),IC=1,3) /
     S 0.79709806E+01,-0.74805226E+00, 0.00000000E+00/
      DATA (GB( 1,12,IC),IC=1,3) /
     S 0.79709806E+01, 0.18377807E+02, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  200.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 2,11,IC),IC=1,3) /
     S 0.10904073E+02,-0.10571588E+01, 0.00000000E+00/
      DATA (GB( 2,11,IC),IC=1,3) /
     S 0.10904073E+02, 0.24728346E+02, 0.10000000E+01/
      DATA (GA( 2,12,IC),IC=1,3) /
     S 0.75400737E+01,-0.56252739E+00, 0.00000000E+00/
      DATA (GB( 2,12,IC),IC=1,3) /
     S 0.75400737E+01, 0.17643148E+02, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  212.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 3,11,IC),IC=1,3) /
     S 0.89126838E+01,-0.74864953E+00, 0.00000000E+00/
      DATA (GB( 3,11,IC),IC=1,3) /
     S 0.89126838E+01, 0.20551342E+02, 0.10000000E+01/
      DATA (GA( 3,12,IC),IC=1,3) /
     S 0.81804377E+01,-0.46188072E+00, 0.00000000E+00/
      DATA (GB( 3,12,IC),IC=1,3) /
     S 0.81804377E+01, 0.19296161E+02, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  225.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 4,11,IC),IC=1,3) /
     S 0.85622405E+01,-0.58705980E+00, 0.00000000E+00/
      DATA (GB( 4,11,IC),IC=1,3) /
     S 0.85622405E+01, 0.19955244E+02, 0.10000000E+01/
      DATA (GA( 4,12,IC),IC=1,3) /
     S 0.10564339E+02,-0.40712065E+00, 0.00000000E+00/
      DATA (GB( 4,12,IC),IC=1,3) /
     S 0.10564339E+02, 0.24951120E+02, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  237.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 5,11,IC),IC=1,3) /
     S 0.94892164E+01,-0.49305772E+00, 0.00000000E+00/
      DATA (GB( 5,11,IC),IC=1,3) /
     S 0.94892164E+01, 0.22227100E+02, 0.10000000E+01/
      DATA (GA( 5,12,IC),IC=1,3) /
     S 0.46896789E+02,-0.15295996E+01, 0.00000000E+00/
      DATA (GB( 5,12,IC),IC=1,3) /
     S 0.46896789E+02, 0.10957372E+03, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  250.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 6,11,IC),IC=1,3) /
     S 0.13580937E+02,-0.51461431E+00, 0.00000000E+00/
      DATA (GB( 6,11,IC),IC=1,3) /
     S 0.13580937E+02, 0.31770288E+02, 0.10000000E+01/
      DATA (GA( 6,12,IC),IC=1,3) /
     S-0.30926524E+01, 0.43555255E+00, 0.00000000E+00/
      DATA (GB( 6,12,IC),IC=1,3) /
     S-0.30926524E+01,-0.67432659E+01, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  262.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 7,11,IC),IC=1,3) /
     S-0.32050918E+03, 0.12373350E+02, 0.00000000E+00/
      DATA (GB( 7,11,IC),IC=1,3) /
     S-0.32050918E+03,-0.74061287E+03, 0.10000000E+01/
      DATA (GA( 7,12,IC),IC=1,3) /
     S 0.85742941E+00, 0.50380874E+00, 0.00000000E+00/
      DATA (GB( 7,12,IC),IC=1,3) /
     S 0.85742941E+00, 0.24550746E+01, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  275.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 8,11,IC),IC=1,3) /
     S-0.37133165E+01, 0.44809588E+00, 0.00000000E+00/
      DATA (GB( 8,11,IC),IC=1,3) /
     S-0.37133165E+01,-0.81329826E+01, 0.10000000E+01/
      DATA (GA( 8,12,IC),IC=1,3) /
     S 0.19164038E+01, 0.68537352E+00, 0.00000000E+00/
      DATA (GB( 8,12,IC),IC=1,3) /
     S 0.19164038E+01, 0.49089917E+01, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  287.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 9,11,IC),IC=1,3) /
     S 0.18890836E+00, 0.46548918E+00, 0.00000000E+00/
      DATA (GB( 9,11,IC),IC=1,3) /
     S 0.18890836E+00, 0.90279822E+00, 0.10000000E+01/
      DATA (GA( 9,12,IC),IC=1,3) /
     S 0.23513199E+01, 0.89437630E+00, 0.00000000E+00/
      DATA (GB( 9,12,IC),IC=1,3) /
     S 0.23513199E+01, 0.59008712E+01, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  300.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA(10,11,IC),IC=1,3) /
     S 0.14209226E+01, 0.59121475E+00, 0.00000000E+00/
      DATA (GB(10,11,IC),IC=1,3) /
     S 0.14209226E+01, 0.37532746E+01, 0.10000000E+01/
      DATA (GA(10,12,IC),IC=1,3) /
     S 0.25566644E+01, 0.11127003E+01, 0.00000000E+00/
      DATA (GB(10,12,IC),IC=1,3) /
     S 0.25566644E+01, 0.63532616E+01, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  312.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA(11,11,IC),IC=1,3) /
     S 0.19817679E+01, 0.74676119E+00, 0.00000000E+00/
      DATA (GB(11,11,IC),IC=1,3) /
     S 0.19817679E+01, 0.50437916E+01, 0.10000000E+01/
      DATA (GA(11,12,IC),IC=1,3) /
     S 0.26555181E+01, 0.13329782E+01, 0.00000000E+00/
      DATA (GB(11,12,IC),IC=1,3) /
     S 0.26555181E+01, 0.65558627E+01, 0.10000000E+01/
C
C
C
C
C
C-- END WATER VAPOR
C
C
C-- CO2 -- INT.2 -- 500-800 CM-1 --- FROM ABS225 ----------------------
C
C
C
C-- FIU = 0.8 + MAX(0.35,(7-IU)*0.9)  , X/T,  9
C
C----- INTERVAL = 2 ----- T =  187.5
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA( 1,13,IC),IC=1,3) /
     S 0.87668459E-01, 0.13845511E+01, 0.00000000E+00/
      DATA (GB( 1,13,IC),IC=1,3) /
     S 0.87668459E-01, 0.23203798E+01, 0.10000000E+01/
      DATA (GA( 1,14,IC),IC=1,3) /
     S 0.74878820E-01, 0.11718758E+01, 0.00000000E+00/
      DATA (GB( 1,14,IC),IC=1,3) /
     S 0.74878820E-01, 0.20206726E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  200.0
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA( 2,13,IC),IC=1,3) /
     S 0.83754276E-01, 0.13187042E+01, 0.00000000E+00/
      DATA (GB( 2,13,IC),IC=1,3) /
     S 0.83754276E-01, 0.22288925E+01, 0.10000000E+01/
      DATA (GA( 2,14,IC),IC=1,3) /
     S 0.71650966E-01, 0.11216131E+01, 0.00000000E+00/
      DATA (GB( 2,14,IC),IC=1,3) /
     S 0.71650966E-01, 0.19441824E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  212.5
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA( 3,13,IC),IC=1,3) /
     S 0.80460283E-01, 0.12644396E+01, 0.00000000E+00/
      DATA (GB( 3,13,IC),IC=1,3) /
     S 0.80460283E-01, 0.21515593E+01, 0.10000000E+01/
      DATA (GA( 3,14,IC),IC=1,3) /
     S 0.68979615E-01, 0.10809473E+01, 0.00000000E+00/
      DATA (GB( 3,14,IC),IC=1,3) /
     S 0.68979615E-01, 0.18807257E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  225.0
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA( 4,13,IC),IC=1,3) /
     S 0.77659686E-01, 0.12191543E+01, 0.00000000E+00/
      DATA (GB( 4,13,IC),IC=1,3) /
     S 0.77659686E-01, 0.20855896E+01, 0.10000000E+01/
      DATA (GA( 4,14,IC),IC=1,3) /
     S 0.66745345E-01, 0.10476396E+01, 0.00000000E+00/
      DATA (GB( 4,14,IC),IC=1,3) /
     S 0.66745345E-01, 0.18275618E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  237.5
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA( 5,13,IC),IC=1,3) /
     S 0.75257056E-01, 0.11809511E+01, 0.00000000E+00/
      DATA (GB( 5,13,IC),IC=1,3) /
     S 0.75257056E-01, 0.20288489E+01, 0.10000000E+01/
      DATA (GA( 5,14,IC),IC=1,3) /
     S 0.64857571E-01, 0.10200373E+01, 0.00000000E+00/
      DATA (GB( 5,14,IC),IC=1,3) /
     S 0.64857571E-01, 0.17825910E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  250.0
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA( 6,13,IC),IC=1,3) /
     S 0.73179175E-01, 0.11484154E+01, 0.00000000E+00/
      DATA (GB( 6,13,IC),IC=1,3) /
     S 0.73179175E-01, 0.19796791E+01, 0.10000000E+01/
      DATA (GA( 6,14,IC),IC=1,3) /
     S 0.63248495E-01, 0.99692726E+00, 0.00000000E+00/
      DATA (GB( 6,14,IC),IC=1,3) /
     S 0.63248495E-01, 0.17442308E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  262.5
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA( 7,13,IC),IC=1,3) /
     S 0.71369063E-01, 0.11204723E+01, 0.00000000E+00/
      DATA (GB( 7,13,IC),IC=1,3) /
     S 0.71369063E-01, 0.19367778E+01, 0.10000000E+01/
      DATA (GA( 7,14,IC),IC=1,3) /
     S 0.61866970E-01, 0.97740923E+00, 0.00000000E+00/
      DATA (GB( 7,14,IC),IC=1,3) /
     S 0.61866970E-01, 0.17112809E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  275.0
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA( 8,13,IC),IC=1,3) /
     S 0.69781812E-01, 0.10962918E+01, 0.00000000E+00/
      DATA (GB( 8,13,IC),IC=1,3) /
     S 0.69781812E-01, 0.18991112E+01, 0.10000000E+01/
      DATA (GA( 8,14,IC),IC=1,3) /
     S 0.60673632E-01, 0.96080188E+00, 0.00000000E+00/
      DATA (GB( 8,14,IC),IC=1,3) /
     S 0.60673632E-01, 0.16828137E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  287.5
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA( 9,13,IC),IC=1,3) /
     S 0.68381606E-01, 0.10752229E+01, 0.00000000E+00/
      DATA (GB( 9,13,IC),IC=1,3) /
     S 0.68381606E-01, 0.18658501E+01, 0.10000000E+01/
      DATA (GA( 9,14,IC),IC=1,3) /
     S 0.59637277E-01, 0.94657562E+00, 0.00000000E+00/
      DATA (GB( 9,14,IC),IC=1,3) /
     S 0.59637277E-01, 0.16580908E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  300.0
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA(10,13,IC),IC=1,3) /
     S 0.67139539E-01, 0.10567474E+01, 0.00000000E+00/
      DATA (GB(10,13,IC),IC=1,3) /
     S 0.67139539E-01, 0.18363226E+01, 0.10000000E+01/
      DATA (GA(10,14,IC),IC=1,3) /
     S 0.58732178E-01, 0.93430511E+00, 0.00000000E+00/
      DATA (GB(10,14,IC),IC=1,3) /
     S 0.58732178E-01, 0.16365014E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  312.5
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA(11,13,IC),IC=1,3) /
     S 0.66032012E-01, 0.10404465E+01, 0.00000000E+00/
      DATA (GB(11,13,IC),IC=1,3) /
     S 0.66032012E-01, 0.18099779E+01, 0.10000000E+01/
      DATA (GA(11,14,IC),IC=1,3) /
     S 0.57936092E-01, 0.92363528E+00, 0.00000000E+00/
      DATA (GB(11,14,IC),IC=1,3) /
     S 0.57936092E-01, 0.16175164E+01, 0.10000000E+01/
C
C
C
C
C
C
C
C
C
C
C-- CARBON DIOXIDE LINES IN THE WINDOW REGION (800-1250 CM-1)
C
C
C-- G = 0.0
C
C
C----- INTERVAL = 4 ----- T =  187.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 1,15,IC),IC=1,3) /
     S 0.13230067E+02, 0.22042132E+02, 0.00000000E+00/
      DATA (GB( 1,15,IC),IC=1,3) /
     S 0.13230067E+02, 0.22051750E+02, 0.10000000E+01/
      DATA (GA( 1,16,IC),IC=1,3) /
     S 0.13183816E+02, 0.22169501E+02, 0.00000000E+00/
      DATA (GB( 1,16,IC),IC=1,3) /
     S 0.13183816E+02, 0.22178972E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  200.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 2,15,IC),IC=1,3) /
     S 0.13213564E+02, 0.22107298E+02, 0.00000000E+00/
      DATA (GB( 2,15,IC),IC=1,3) /
     S 0.13213564E+02, 0.22116850E+02, 0.10000000E+01/
      DATA (GA( 2,16,IC),IC=1,3) /
     S 0.13189991E+02, 0.22270075E+02, 0.00000000E+00/
      DATA (GB( 2,16,IC),IC=1,3) /
     S 0.13189991E+02, 0.22279484E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  212.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 3,15,IC),IC=1,3) /
     S 0.13209140E+02, 0.22180915E+02, 0.00000000E+00/
      DATA (GB( 3,15,IC),IC=1,3) /
     S 0.13209140E+02, 0.22190410E+02, 0.10000000E+01/
      DATA (GA( 3,16,IC),IC=1,3) /
     S 0.13209485E+02, 0.22379193E+02, 0.00000000E+00/
      DATA (GB( 3,16,IC),IC=1,3) /
     S 0.13209485E+02, 0.22388551E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  225.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 4,15,IC),IC=1,3) /
     S 0.13213894E+02, 0.22259478E+02, 0.00000000E+00/
      DATA (GB( 4,15,IC),IC=1,3) /
     S 0.13213894E+02, 0.22268925E+02, 0.10000000E+01/
      DATA (GA( 4,16,IC),IC=1,3) /
     S 0.13238789E+02, 0.22492992E+02, 0.00000000E+00/
      DATA (GB( 4,16,IC),IC=1,3) /
     S 0.13238789E+02, 0.22502309E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  237.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 5,15,IC),IC=1,3) /
     S 0.13225963E+02, 0.22341039E+02, 0.00000000E+00/
      DATA (GB( 5,15,IC),IC=1,3) /
     S 0.13225963E+02, 0.22350445E+02, 0.10000000E+01/
      DATA (GA( 5,16,IC),IC=1,3) /
     S 0.13275017E+02, 0.22608508E+02, 0.00000000E+00/
      DATA (GB( 5,16,IC),IC=1,3) /
     S 0.13275017E+02, 0.22617792E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  250.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 6,15,IC),IC=1,3) /
     S 0.13243806E+02, 0.22424247E+02, 0.00000000E+00/
      DATA (GB( 6,15,IC),IC=1,3) /
     S 0.13243806E+02, 0.22433617E+02, 0.10000000E+01/
      DATA (GA( 6,16,IC),IC=1,3) /
     S 0.13316096E+02, 0.22723843E+02, 0.00000000E+00/
      DATA (GB( 6,16,IC),IC=1,3) /
     S 0.13316096E+02, 0.22733099E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  262.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 7,15,IC),IC=1,3) /
     S 0.13266104E+02, 0.22508089E+02, 0.00000000E+00/
      DATA (GB( 7,15,IC),IC=1,3) /
     S 0.13266104E+02, 0.22517429E+02, 0.10000000E+01/
      DATA (GA( 7,16,IC),IC=1,3) /
     S 0.13360555E+02, 0.22837837E+02, 0.00000000E+00/
      DATA (GB( 7,16,IC),IC=1,3) /
     S 0.13360555E+02, 0.22847071E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  275.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 8,15,IC),IC=1,3) /
     S 0.13291782E+02, 0.22591771E+02, 0.00000000E+00/
      DATA (GB( 8,15,IC),IC=1,3) /
     S 0.13291782E+02, 0.22601086E+02, 0.10000000E+01/
      DATA (GA( 8,16,IC),IC=1,3) /
     S 0.13407324E+02, 0.22949751E+02, 0.00000000E+00/
      DATA (GB( 8,16,IC),IC=1,3) /
     S 0.13407324E+02, 0.22958967E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  287.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 9,15,IC),IC=1,3) /
     S 0.13319961E+02, 0.22674661E+02, 0.00000000E+00/
      DATA (GB( 9,15,IC),IC=1,3) /
     S 0.13319961E+02, 0.22683956E+02, 0.10000000E+01/
      DATA (GA( 9,16,IC),IC=1,3) /
     S 0.13455544E+02, 0.23059032E+02, 0.00000000E+00/
      DATA (GB( 9,16,IC),IC=1,3) /
     S 0.13455544E+02, 0.23068234E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  300.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA(10,15,IC),IC=1,3) /
     S 0.13349927E+02, 0.22756246E+02, 0.00000000E+00/
      DATA (GB(10,15,IC),IC=1,3) /
     S 0.13349927E+02, 0.22765522E+02, 0.10000000E+01/
      DATA (GA(10,16,IC),IC=1,3) /
     S 0.13504450E+02, 0.23165146E+02, 0.00000000E+00/
      DATA (GB(10,16,IC),IC=1,3) /
     S 0.13504450E+02, 0.23174336E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  312.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA(11,15,IC),IC=1,3) /
     S 0.13381108E+02, 0.22836093E+02, 0.00000000E+00/
      DATA (GB(11,15,IC),IC=1,3) /
     S 0.13381108E+02, 0.22845354E+02, 0.10000000E+01/
      DATA (GA(11,16,IC),IC=1,3) /
     S 0.13553282E+02, 0.23267456E+02, 0.00000000E+00/
      DATA (GB(11,16,IC),IC=1,3) /
     S 0.13553282E+02, 0.23276638E+02, 0.10000000E+01/

C     ------------------------------------------------------------------
      DATA (( XP(  J,K),J=1,6),       K=1,6) /
     S 0.46430621E+02, 0.12928299E+03, 0.20732648E+03,
     S 0.31398411E+03, 0.18373177E+03,-0.11412303E+03,
     S 0.73604774E+02, 0.27887914E+03, 0.27076947E+03,
     S-0.57322111E+02,-0.64742459E+02, 0.87238280E+02,
     S 0.37050866E+02, 0.20498759E+03, 0.37558029E+03,
     S 0.17401171E+03,-0.13350302E+03,-0.37651795E+02,
     S 0.14930141E+02, 0.89161160E+02, 0.17793062E+03,
     S 0.93433860E+02,-0.70646020E+02,-0.26373150E+02,
     S 0.40386780E+02, 0.10855270E+03, 0.50755010E+02,
     S-0.31496190E+02, 0.12791300E+00, 0.18017770E+01,
     S 0.90811926E+01, 0.75073923E+02, 0.24654438E+03,
     S 0.39332612E+03, 0.29385281E+03, 0.89107921E+02 /
C
C
C*         1.0     PLANCK FUNCTIONS AND GRADIENTS
C                  ------------------------------
C
 100  CONTINUE
C
      DO 102 JK = 1 , KFLEV+1
      DO 101 JL = 1, KDLON
      PBINT(JL,JK) = 0.
 101  CONTINUE
 102  CONTINUE
      DO 103 JL = 1, KDLON
      PBSUIN(JL) = 0.
 103  CONTINUE
C
      DO 141 JNU=1,Ninter
C
C
C*         1.1   LEVELS FROM SURFACE TO KFLEV
C                ----------------------------
C
 110  CONTINUE
C
      DO 112 JK = 1 , KFLEV
      DO 111 JL = 1, KDLON
      ZTI(JL)=(PTL(JL,JK)-TSTAND)/TSTAND
      ZRES(JL) = XP(1,JNU)+ZTI(JL)*(XP(2,JNU)+ZTI(JL)*(XP(3,JNU)
     S       +ZTI(JL)*(XP(4,JNU)+ZTI(JL)*(XP(5,JNU)+ZTI(JL)*(XP(6,JNU)
     S       )))))
      PBINT(JL,JK)=PBINT(JL,JK)+ZRES(JL)
      PB(JL,JNU,JK)= ZRES(JL)
      ZBLEV(JL,JK) = ZRES(JL)
      ZTI2(JL)=(PTAVE(JL,JK)-TSTAND)/TSTAND
      ZRES2(JL)=XP(1,JNU)+ZTI2(JL)*(XP(2,JNU)+ZTI2(JL)*(XP(3,JNU)
     S     +ZTI2(JL)*(XP(4,JNU)+ZTI2(JL)*(XP(5,JNU)+ZTI2(JL)*(XP(6,JNU)
     S       )))))
      ZBLAY(JL,JK) = ZRES2(JL)
 111  CONTINUE
 112  CONTINUE
C
C
C*         1.2   TOP OF THE ATMOSPHERE AND SURFACE
C                ---------------------------------
C
 120  CONTINUE
C
      DO 121 JL = 1, KDLON
      ZTI(JL)=(PTL(JL,KFLEV+1)-TSTAND)/TSTAND
      ZTI2(JL) = (PTL(JL,1) + PDT0(JL) - TSTAND) / TSTAND
      ZRES(JL) = XP(1,JNU)+ZTI(JL)*(XP(2,JNU)+ZTI(JL)*(XP(3,JNU)
     S    +ZTI(JL)*(XP(4,JNU)+ZTI(JL)*(XP(5,JNU)+ZTI(JL)*(XP(6,JNU)
     S       )))))
      ZRES2(JL) = XP(1,JNU)+ZTI2(JL)*(XP(2,JNU)+ZTI2(JL)*(XP(3,JNU)
     S    +ZTI2(JL)*(XP(4,JNU)+ZTI2(JL)*(XP(5,JNU)+ZTI2(JL)*(XP(6,JNU)
     S       )))))
      PBINT(JL,KFLEV+1) = PBINT(JL,KFLEV+1)+ZRES(JL)
      PB(JL,JNU,KFLEV+1)= ZRES(JL)
      ZBLEV(JL,KFLEV+1) = ZRES(JL)
      PBTOP(JL,JNU) = ZRES(JL)
      PBSUR(JL,JNU) = ZRES2(JL)
      PBSUIN(JL) = PBSUIN(JL) + ZRES2(JL)
 121  CONTINUE
C
C
C*         1.3   GRADIENTS IN SUB-LAYERS
C                -----------------------
C
 130  CONTINUE
C
      DO 132 JK = 1 , KFLEV
      JK2 = 2 * JK
      JK1 = JK2 - 1
      DO 131 JL = 1, KDLON
      PDBSL(JL,JNU,JK1) = ZBLAY(JL,JK  ) - ZBLEV(JL,JK)
      PDBSL(JL,JNU,JK2) = ZBLEV(JL,JK+1) - ZBLAY(JL,JK)
 131  CONTINUE
 132  CONTINUE
C
 141  CONTINUE
C
C*         2.0   CHOOSE THE RELEVANT SETS OF PADE APPROXIMANTS
C                ---------------------------------------------
C
 200  CONTINUE
C
C
 210  CONTINUE
C
      DO 211 JL=1, KDLON
      ZDSTO1 = (PTL(JL,KFLEV+1)-TINTP(1)) / TSTP
      IXTOX = MAX( 1, MIN( MXIXT, INT( ZDSTO1 + 1. ) ) )
      ZDSTOX = (PTL(JL,KFLEV+1)-TINTP(IXTOX))/TSTP
      IF (ZDSTOX.LT.0.5) THEN
         INDTO=IXTOX
      ELSE
         INDTO=IXTOX+1
      END IF
      INDB(JL)=INDTO
      ZDST1 = (PTL(JL,1)-TINTP(1)) / TSTP
      IXTX = MAX( 1, MIN( MXIXT, INT( ZDST1 + 1. ) ) )
      ZDSTX = (PTL(JL,1)-TINTP(IXTX))/TSTP
      IF (ZDSTX.LT.0.5) THEN
         INDT=IXTX
      ELSE
         INDT=IXTX+1
      END IF
      INDS(JL)=INDT
 211  CONTINUE
C
      DO 214 JF=1,2
      DO 213 JG=1, 8
      DO 212 JL=1, KDLON
      INDSU=INDS(JL)
      PGASUR(JL,JG,JF)=GA(INDSU,2*JG-1,JF)
      PGBSUR(JL,JG,JF)=GB(INDSU,2*JG-1,JF)
      INDTP=INDB(JL)
      PGATOP(JL,JG,JF)=GA(INDTP,2*JG-1,JF)
      PGBTOP(JL,JG,JF)=GB(INDTP,2*JG-1,JF)
 212  CONTINUE
 213  CONTINUE
 214  CONTINUE
C
 220  CONTINUE
C
      DO 225 JK=1,KFLEV
      DO 221 JL=1, KDLON
      ZDST1 = (PTAVE(JL,JK)-TINTP(1)) / TSTP
      IXTX = MAX( 1, MIN( MXIXT, INT( ZDST1 + 1. ) ) )
      ZDSTX = (PTAVE(JL,JK)-TINTP(IXTX))/TSTP
      IF (ZDSTX.LT.0.5) THEN
         INDT=IXTX
      ELSE
         INDT=IXTX+1
      END IF
      INDB(JL)=INDT
 221  CONTINUE
C
      DO 224 JF=1,2
      DO 223 JG=1, 8
      DO 222 JL=1, KDLON
      INDT=INDB(JL)
      PGA(JL,JG,JF,JK)=GA(INDT,2*JG,JF)
      PGB(JL,JG,JF,JK)=GB(INDT,2*JG,JF)
 222  CONTINUE
 223  CONTINUE
 224  CONTINUE
 225  CONTINUE
C
C     ------------------------------------------------------------------
C
      RETURN
      END
      SUBROUTINE LWV(KUAER,KTRAER, KLIM
     R  , PABCU,PB,PBINT,PBSUIN,PBSUR,PBTOP,PDBSL,PEMIS,PPMB,PTAVE
     R  , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP
     S  , PCNTRB,PCTS,PFLUC)
      use dimens_m
      use dimphy
      use YOMCST
      use raddim
      IMPLICIT none
      include "raddimlw.h"
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           CARRIES OUT THE VERTICAL INTEGRATION TO GIVE LONGWAVE
C           FLUXES OR RADIANCES
C
C     METHOD.
C     -------
C
C          1. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING BETWEEN
C     CONTRIBUTIONS BY -  THE NEARBY LAYERS
C                      -  THE DISTANT LAYERS
C                      -  THE BOUNDARY TERMS
C          2. COMPUTES THE CLEAR-SKY DOWNWARD AND UPWARD EMISSIVITIES.
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C-----------------------------------------------------------------------
C
C* ARGUMENTS:
      INTEGER KUAER,KTRAER, KLIM
C
      REAL*8 PABCU(KDLON,NUA,3*KFLEV+1) ! EFFECTIVE ABSORBER AMOUNTS
      REAL*8 PB(KDLON,Ninter,KFLEV+1) ! SPECTRAL HALF-LEVEL PLANCK FUNCTIONS
      REAL*8 PBINT(KDLON,KFLEV+1) ! HALF-LEVEL PLANCK FUNCTIONS
      REAL*8 PBSUR(KDLON,Ninter) ! SURFACE SPECTRAL PLANCK FUNCTION
      REAL*8 PBSUIN(KDLON) ! SURFACE PLANCK FUNCTION
      REAL*8 PBTOP(KDLON,Ninter) ! T.O.A. SPECTRAL PLANCK FUNCTION
      REAL*8 PDBSL(KDLON,Ninter,KFLEV*2) ! SUB-LAYER PLANCK FUNCTION GRADIENT
      REAL*8 PEMIS(KDLON) ! SURFACE EMISSIVITY
      REAL*8 PPMB(KDLON,KFLEV+1) ! HALF-LEVEL PRESSURE (MB)
      REAL*8 PTAVE(KDLON,KFLEV) ! TEMPERATURE
      REAL*8 PGA(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
      REAL*8 PGB(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
      REAL*8 PGASUR(KDLON,8,2) ! PADE APPROXIMANTS
      REAL*8 PGBSUR(KDLON,8,2) ! PADE APPROXIMANTS
      REAL*8 PGATOP(KDLON,8,2) ! PADE APPROXIMANTS
      REAL*8 PGBTOP(KDLON,8,2) ! PADE APPROXIMANTS
C
      REAL*8 PCNTRB(KDLON,KFLEV+1,KFLEV+1) ! CLEAR-SKY ENERGY EXCHANGE MATRIX
      REAL*8 PCTS(KDLON,KFLEV) ! COOLING-TO-SPACE TERM
      REAL*8 PFLUC(KDLON,2,KFLEV+1) ! CLEAR-SKY RADIATIVE FLUXES
C-----------------------------------------------------------------------
C LOCAL VARIABLES:
      REAL*8 ZADJD(KDLON,KFLEV+1)
      REAL*8 ZADJU(KDLON,KFLEV+1)
      REAL*8 ZDBDT(KDLON,Ninter,KFLEV)
      REAL*8 ZDISD(KDLON,KFLEV+1)
      REAL*8 ZDISU(KDLON,KFLEV+1)
C
      INTEGER jk, jl
C-----------------------------------------------------------------------
C
      DO 112 JK=1,KFLEV+1
      DO 111 JL=1, KDLON
      ZADJD(JL,JK)=0.
      ZADJU(JL,JK)=0.
      ZDISD(JL,JK)=0.
      ZDISU(JL,JK)=0.
 111  CONTINUE
 112  CONTINUE
C
      DO 114 JK=1,KFLEV
      DO 113 JL=1, KDLON
      PCTS(JL,JK)=0.
 113  CONTINUE
 114  CONTINUE
C
C* CONTRIBUTION FROM ADJACENT LAYERS
C
      CALL LWVN(KUAER,KTRAER
     R  , PABCU,PDBSL,PGA,PGB
     S  , ZADJD,ZADJU,PCNTRB,ZDBDT)
C* CONTRIBUTION FROM DISTANT LAYERS
C
      CALL LWVD(KUAER,KTRAER
     R  , PABCU,ZDBDT,PGA,PGB
     S  , PCNTRB,ZDISD,ZDISU)
C
C* EXCHANGE WITH THE BOUNDARIES
C
      CALL LWVB(KUAER,KTRAER, KLIM
     R  , PABCU,ZADJD,ZADJU,PB,PBINT,PBSUIN,PBSUR,PBTOP
     R  , ZDISD,ZDISU,PEMIS,PPMB
     R  , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP
     S  , PCTS,PFLUC)
C
C
      RETURN
      END
      SUBROUTINE LWVB(KUAER,KTRAER, KLIM
     R  , PABCU,PADJD,PADJU,PB,PBINT,PBSUI,PBSUR,PBTOP
     R  , PDISD,PDISU,PEMIS,PPMB
     R  , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP
     S  , PCTS,PFLUC)
      use dimens_m
      use dimphy
      use raddim
      use radopt
      IMPLICIT none
      include "raddimlw.h"
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           INTRODUCES THE EFFECTS OF THE BOUNDARIES IN THE VERTICAL
C           INTEGRATION
C
C     METHOD.
C     -------
C
C          1. COMPUTES THE ENERGY EXCHANGE WITH TOP AND SURFACE OF THE
C     ATMOSPHERE
C          2. COMPUTES THE COOLING-TO-SPACE AND HEATING-FROM-GROUND
C     TERMS FOR THE APPROXIMATE COOLING RATE ABOVE 10 HPA
C          3. ADDS UP ALL CONTRIBUTIONS TO GET THE CLEAR-SKY FLUXES
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C        Voigt lines (loop 2413 to 2427)  - JJM & PhD - 01/96
C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
      INTEGER KUAER,KTRAER, KLIM
C
      REAL*8 PABCU(KDLON,NUA,3*KFLEV+1) ! ABSORBER AMOUNTS
      REAL*8 PADJD(KDLON,KFLEV+1) ! CONTRIBUTION BY ADJACENT LAYERS
      REAL*8 PADJU(KDLON,KFLEV+1) ! CONTRIBUTION BY ADJACENT LAYERS
      REAL*8 PB(KDLON,Ninter,KFLEV+1) ! SPECTRAL HALF-LEVEL PLANCK FUNCTIONS
      REAL*8 PBINT(KDLON,KFLEV+1) ! HALF-LEVEL PLANCK FUNCTIONS
      REAL*8 PBSUR(KDLON,Ninter) ! SPECTRAL SURFACE PLANCK FUNCTION
      REAL*8 PBSUI(KDLON) ! SURFACE PLANCK FUNCTION
      REAL*8 PBTOP(KDLON,Ninter) ! SPECTRAL T.O.A. PLANCK FUNCTION
      REAL*8 PDISD(KDLON,KFLEV+1) ! CONTRIBUTION BY DISTANT LAYERS
      REAL*8 PDISU(KDLON,KFLEV+1) ! CONTRIBUTION BY DISTANT LAYERS
      REAL*8 PEMIS(KDLON) ! SURFACE EMISSIVITY
      REAL*8 PPMB(KDLON,KFLEV+1) ! PRESSURE MB
      REAL*8 PGA(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
      REAL*8 PGB(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
      REAL*8 PGASUR(KDLON,8,2) ! SURFACE PADE APPROXIMANTS
      REAL*8 PGBSUR(KDLON,8,2) ! SURFACE PADE APPROXIMANTS
      REAL*8 PGATOP(KDLON,8,2) ! T.O.A. PADE APPROXIMANTS
      REAL*8 PGBTOP(KDLON,8,2) ! T.O.A. PADE APPROXIMANTS
C
      REAL*8 PFLUC(KDLON,2,KFLEV+1) ! CLEAR-SKY RADIATIVE FLUXES
      REAL*8 PCTS(KDLON,KFLEV) ! COOLING-TO-SPACE TERM
C
C* LOCAL VARIABLES:
C
      REAL*8 ZBGND(KDLON)
      REAL*8 ZFD(KDLON)
      REAL*8  ZFN10(KDLON)
      REAL*8 ZFU(KDLON)
      REAL*8  ZTT(KDLON,NTRA)
      REAL*8 ZTT1(KDLON,NTRA)
      REAL*8 ZTT2(KDLON,NTRA)
      REAL*8  ZUU(KDLON,NUA) 
      REAL*8 ZCNSOL(KDLON)
      REAL*8 ZCNTOP(KDLON)
C
      INTEGER jk, jl, ja
      INTEGER jstra, jstru
      INTEGER ind1, ind2, ind3, ind4, in, jlim
      REAL*8 zctstr
C-----------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
 100  CONTINUE
C
C
C*         1.2     INITIALIZE TRANSMISSION FUNCTIONS
C                  ---------------------------------
C
 120  CONTINUE
C
      DO 122 JA=1,NTRA
      DO 121 JL=1, KDLON
      ZTT (JL,JA)=1.0
      ZTT1(JL,JA)=1.0
      ZTT2(JL,JA)=1.0
 121  CONTINUE
 122  CONTINUE
C
      DO 124 JA=1,NUA
      DO 123 JL=1, KDLON
      ZUU(JL,JA)=1.0
 123  CONTINUE
 124  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.      VERTICAL INTEGRATION
C                  --------------------
C
 200  CONTINUE
C
      IND1=0
      IND3=0
      IND4=1
      IND2=1
C
C
C*         2.3     EXCHANGE WITH TOP OF THE ATMOSPHERE
C                  -----------------------------------
C
 230  CONTINUE
C
      DO 235 JK = 1 , KFLEV
      IN=(JK-1)*NG1P1+1
C
      DO 232 JA=1,KUAER
      DO 231 JL=1, KDLON
      ZUU(JL,JA)=PABCU(JL,JA,IN)
 231  CONTINUE
 232  CONTINUE
C
C
      CALL LWTT(PGATOP(1,1,1), PGBTOP(1,1,1), ZUU, ZTT)
C
      DO 234 JL = 1, KDLON
      ZCNTOP(JL)=PBTOP(JL,1)*ZTT(JL,1)          *ZTT(JL,10)
     2      +PBTOP(JL,2)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     3      +PBTOP(JL,3)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     4      +PBTOP(JL,4)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     5      +PBTOP(JL,5)*ZTT(JL,3)          *ZTT(JL,14)
     6      +PBTOP(JL,6)*ZTT(JL,6)          *ZTT(JL,15)
      ZFD(JL)=ZCNTOP(JL)-PBINT(JL,JK)-PDISD(JL,JK)-PADJD(JL,JK)
      PFLUC(JL,2,JK)=ZFD(JL)
 234  CONTINUE
C
 235  CONTINUE
C
      JK = KFLEV+1
      IN=(JK-1)*NG1P1+1
C
      DO 236 JL = 1, KDLON
      ZCNTOP(JL)= PBTOP(JL,1)
     1   + PBTOP(JL,2)
     2   + PBTOP(JL,3)
     3   + PBTOP(JL,4)
     4   + PBTOP(JL,5)
     5   + PBTOP(JL,6)
      ZFD(JL)=ZCNTOP(JL)-PBINT(JL,JK)-PDISD(JL,JK)-PADJD(JL,JK)
      PFLUC(JL,2,JK)=ZFD(JL)
 236  CONTINUE
C
C*         2.4     COOLING-TO-SPACE OF LAYERS ABOVE 10 HPA
C                  ---------------------------------------
C
 240  CONTINUE
C
C
C*         2.4.1   INITIALIZATION
C                  --------------
C
 2410 CONTINUE
C
      JLIM = KFLEV
C
      IF (.NOT.LEVOIGT) THEN
      DO 2412 JK = KFLEV,1,-1
      IF(PPMB(1,JK).LT.10.0) THEN
         JLIM=JK
      ENDIF   
 2412 CONTINUE
      ENDIF
      KLIM=JLIM
C
      IF (.NOT.LEVOIGT) THEN
        DO 2414 JA=1,KTRAER
        DO 2413 JL=1, KDLON
        ZTT1(JL,JA)=1.0
 2413   CONTINUE
 2414   CONTINUE
C
C*         2.4.2   LOOP OVER LAYERS ABOVE 10 HPA
C                  -----------------------------
C
 2420   CONTINUE
C
        DO 2427 JSTRA = KFLEV,JLIM,-1
        JSTRU=(JSTRA-1)*NG1P1+1
C
        DO 2423 JA=1,KUAER
        DO 2422 JL=1, KDLON
        ZUU(JL,JA)=PABCU(JL,JA,JSTRU)
 2422   CONTINUE
 2423   CONTINUE
C
C
        CALL LWTT(PGA(1,1,1,JSTRA), PGB(1,1,1,JSTRA), ZUU, ZTT)
C
        DO 2424 JL = 1, KDLON
        ZCTSTR =
     1   (PB(JL,1,JSTRA)+PB(JL,1,JSTRA+1))
     1       *(ZTT1(JL,1)           *ZTT1(JL,10)
     1       - ZTT (JL,1)           *ZTT (JL,10))
     2  +(PB(JL,2,JSTRA)+PB(JL,2,JSTRA+1))
     2       *(ZTT1(JL,2)*ZTT1(JL,7)*ZTT1(JL,11)
     2       - ZTT (JL,2)*ZTT (JL,7)*ZTT (JL,11))
     3  +(PB(JL,3,JSTRA)+PB(JL,3,JSTRA+1))
     3       *(ZTT1(JL,4)*ZTT1(JL,8)*ZTT1(JL,12)
     3       - ZTT (JL,4)*ZTT (JL,8)*ZTT (JL,12))
     4  +(PB(JL,4,JSTRA)+PB(JL,4,JSTRA+1))
     4       *(ZTT1(JL,5)*ZTT1(JL,9)*ZTT1(JL,13)
     4       - ZTT (JL,5)*ZTT (JL,9)*ZTT (JL,13))
     5  +(PB(JL,5,JSTRA)+PB(JL,5,JSTRA+1))
     5       *(ZTT1(JL,3)           *ZTT1(JL,14)
     5       - ZTT (JL,3)           *ZTT (JL,14))
     6  +(PB(JL,6,JSTRA)+PB(JL,6,JSTRA+1))
     6       *(ZTT1(JL,6)           *ZTT1(JL,15)
     6       - ZTT (JL,6)           *ZTT (JL,15))
        PCTS(JL,JSTRA)=ZCTSTR*0.5
 2424   CONTINUE
        DO 2426 JA=1,KTRAER
        DO 2425 JL=1, KDLON
        ZTT1(JL,JA)=ZTT(JL,JA)
 2425   CONTINUE
 2426   CONTINUE
 2427   CONTINUE
      ENDIF
C Mise a zero de securite pour PCTS en cas de LEVOIGT
      IF(LEVOIGT)THEN
        DO 2429 JSTRA = 1,KFLEV
        DO 2428 JL = 1, KDLON
          PCTS(JL,JSTRA)=0.
 2428   CONTINUE
 2429   CONTINUE
      ENDIF
C
C
C*         2.5     EXCHANGE WITH LOWER LIMIT
C                  -------------------------
C
 250  CONTINUE
C
      DO 251 JL = 1, KDLON
      ZBGND(JL)=PBSUI(JL)*PEMIS(JL)-(1.-PEMIS(JL))
     S               *PFLUC(JL,2,1)-PBINT(JL,1)
 251  CONTINUE
C
      JK = 1
      IN=(JK-1)*NG1P1+1
C
      DO 252 JL = 1, KDLON
      ZCNSOL(JL)=PBSUR(JL,1)
     1 +PBSUR(JL,2)
     2 +PBSUR(JL,3)
     3 +PBSUR(JL,4)
     4 +PBSUR(JL,5)
     5 +PBSUR(JL,6)
      ZCNSOL(JL)=ZCNSOL(JL)*ZBGND(JL)/PBSUI(JL)
      ZFU(JL)=ZCNSOL(JL)+PBINT(JL,JK)-PDISU(JL,JK)-PADJU(JL,JK)
      PFLUC(JL,1,JK)=ZFU(JL)
 252  CONTINUE
C
      DO 257 JK = 2 , KFLEV+1
      IN=(JK-1)*NG1P1+1
C
C
      DO 255 JA=1,KUAER
      DO 254 JL=1, KDLON
      ZUU(JL,JA)=PABCU(JL,JA,1)-PABCU(JL,JA,IN)
 254  CONTINUE
 255  CONTINUE
C
C
      CALL LWTT(PGASUR(1,1,1), PGBSUR(1,1,1), ZUU, ZTT)
C
      DO 256 JL = 1, KDLON
      ZCNSOL(JL)=PBSUR(JL,1)*ZTT(JL,1)          *ZTT(JL,10)
     2      +PBSUR(JL,2)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     3      +PBSUR(JL,3)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     4      +PBSUR(JL,4)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     5      +PBSUR(JL,5)*ZTT(JL,3)          *ZTT(JL,14)
     6      +PBSUR(JL,6)*ZTT(JL,6)          *ZTT(JL,15)
      ZCNSOL(JL)=ZCNSOL(JL)*ZBGND(JL)/PBSUI(JL)
      ZFU(JL)=ZCNSOL(JL)+PBINT(JL,JK)-PDISU(JL,JK)-PADJU(JL,JK)
      PFLUC(JL,1,JK)=ZFU(JL)
 256  CONTINUE
C
C
 257  CONTINUE
C
C
C
C*         2.7     CLEAR-SKY FLUXES
C                  ----------------
C
 270  CONTINUE
C
      IF (.NOT.LEVOIGT) THEN
      DO 271 JL = 1, KDLON
      ZFN10(JL) = PFLUC(JL,1,JLIM) + PFLUC(JL,2,JLIM)
 271  CONTINUE
      DO 273 JK = JLIM+1,KFLEV+1
      DO 272 JL = 1, KDLON
      ZFN10(JL) = ZFN10(JL) + PCTS(JL,JK-1)
      PFLUC(JL,1,JK) = ZFN10(JL)
      PFLUC(JL,2,JK) = 0.
 272  CONTINUE
 273  CONTINUE
      ENDIF
C
C     ------------------------------------------------------------------
C
      RETURN
      END
      SUBROUTINE LWVD(KUAER,KTRAER
     S  , PABCU,PDBDT
     R  , PGA,PGB
     S  , PCNTRB,PDISD,PDISU)
      use dimens_m
      use dimphy
      use raddim
      IMPLICIT none
      include "raddimlw.h"
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           CARRIES OUT THE VERTICAL INTEGRATION ON THE DISTANT LAYERS
C
C     METHOD.
C     -------
C
C          1. PERFORMS THE VERTICAL INTEGRATION CORRESPONDING TO THE
C     CONTRIBUTIONS OF THE DISTANT LAYERS USING TRAPEZOIDAL RULE
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C-----------------------------------------------------------------------
C* ARGUMENTS:
C
      INTEGER KUAER,KTRAER
C
      REAL*8 PABCU(KDLON,NUA,3*KFLEV+1) ! ABSORBER AMOUNTS
      REAL*8 PDBDT(KDLON,Ninter,KFLEV) ! LAYER PLANCK FUNCTION GRADIENT
      REAL*8 PGA(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
      REAL*8 PGB(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
C
      REAL*8 PCNTRB(KDLON,KFLEV+1,KFLEV+1) ! ENERGY EXCHANGE MATRIX
      REAL*8 PDISD(KDLON,KFLEV+1) !  CONTRIBUTION BY DISTANT LAYERS
      REAL*8 PDISU(KDLON,KFLEV+1) !  CONTRIBUTION BY DISTANT LAYERS
C
C* LOCAL VARIABLES:
C
      REAL*8 ZGLAYD(KDLON)
      REAL*8 ZGLAYU(KDLON)
      REAL*8 ZTT(KDLON,NTRA)
      REAL*8 ZTT1(KDLON,NTRA)
      REAL*8 ZTT2(KDLON,NTRA)
C
      INTEGER jl, jk, ja, ikp1, ikn, ikd1, jkj, ikd2
      INTEGER ikjp1, ikm1, ikj, jlk, iku1, ijkl, iku2
      INTEGER ind1, ind2, ind3, ind4, itt
      REAL*8 zww, zdzxdg, zdzxmg
C
C*         1.    INITIALIZATION
C                --------------
C
 100  CONTINUE
C
C*         1.1     INITIALIZE LAYER CONTRIBUTIONS
C                  ------------------------------
C
 110  CONTINUE
C
      DO 112 JK = 1, KFLEV+1
      DO 111 JL = 1, KDLON
      PDISD(JL,JK) = 0.
      PDISU(JL,JK) = 0.
  111 CONTINUE
  112 CONTINUE
C
C*         1.2     INITIALIZE TRANSMISSION FUNCTIONS
C                  ---------------------------------
C
 120  CONTINUE
C
C
      DO 122 JA = 1, NTRA
      DO 121 JL = 1, KDLON
      ZTT (JL,JA) = 1.0
      ZTT1(JL,JA) = 1.0
      ZTT2(JL,JA) = 1.0
  121 CONTINUE
  122 CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.      VERTICAL INTEGRATION
C                  --------------------
C
 200  CONTINUE
C
      IND1=0
      IND3=0
      IND4=1
      IND2=1
C
C
C*         2.2     CONTRIBUTION FROM DISTANT LAYERS
C                  ---------------------------------
C
 220  CONTINUE
C
C
C*         2.2.1   DISTANT AND ABOVE LAYERS
C                  ------------------------
C
 2210 CONTINUE
C
C
C
C*         2.2.2   FIRST UPPER LEVEL
C                  -----------------
C
 2220 CONTINUE
C
      DO 225 JK = 1 , KFLEV-1
      IKP1=JK+1
      IKN=(JK-1)*NG1P1+1
      IKD1= JK  *NG1P1+1
C
      CALL LWTTM(PGA(1,1,1,JK), PGB(1,1,1,JK)
     2          , PABCU(1,1,IKN),PABCU(1,1,IKD1),ZTT1)
C
C
C
C*         2.2.3   HIGHER UP
C                  ---------
C
 2230 CONTINUE
C
      ITT=1
      DO 224 JKJ=IKP1,KFLEV
      IF(ITT.EQ.1) THEN
         ITT=2
      ELSE
         ITT=1
      ENDIF
      IKJP1=JKJ+1
      IKD2= JKJ  *NG1P1+1
C
      IF(ITT.EQ.1) THEN
         CALL LWTTM(PGA(1,1,1,JKJ),PGB(1,1,1,JKJ)
     2             , PABCU(1,1,IKN),PABCU(1,1,IKD2),ZTT1)
      ELSE
         CALL LWTTM(PGA(1,1,1,JKJ),PGB(1,1,1,JKJ)
     2             , PABCU(1,1,IKN),PABCU(1,1,IKD2),ZTT2)
      ENDIF
C
      DO 2235 JA = 1, KTRAER
      DO 2234 JL = 1, KDLON
      ZTT(JL,JA) = (ZTT1(JL,JA)+ZTT2(JL,JA))*0.5
 2234 CONTINUE
 2235 CONTINUE
C
      DO 2236 JL = 1, KDLON
      ZWW=PDBDT(JL,1,JKJ)*ZTT(JL,1)          *ZTT(JL,10)
     S   +PDBDT(JL,2,JKJ)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     S   +PDBDT(JL,3,JKJ)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     S   +PDBDT(JL,4,JKJ)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     S   +PDBDT(JL,5,JKJ)*ZTT(JL,3)          *ZTT(JL,14)
     S   +PDBDT(JL,6,JKJ)*ZTT(JL,6)          *ZTT(JL,15)
      ZGLAYD(JL)=ZWW
      ZDZXDG=ZGLAYD(JL)
      PDISD(JL,JK)=PDISD(JL,JK)+ZDZXDG
      PCNTRB(JL,JK,IKJP1)=ZDZXDG
 2236 CONTINUE
C
C
 224  CONTINUE
 225  CONTINUE
C
C
C*         2.2.4   DISTANT AND BELOW LAYERS
C                  ------------------------
C
 2240 CONTINUE
C
C
C
C*         2.2.5   FIRST LOWER LEVEL
C                  -----------------
C
 2250 CONTINUE
C
      DO 228 JK=3,KFLEV+1
      IKN=(JK-1)*NG1P1+1
      IKM1=JK-1
      IKJ=JK-2
      IKU1= IKJ  *NG1P1+1
C
C
      CALL LWTTM(PGA(1,1,1,IKJ),PGB(1,1,1,IKJ)
     2          , PABCU(1,1,IKU1),PABCU(1,1,IKN),ZTT1)
C
C
C
C*         2.2.6   DOWN BELOW
C                  ----------
C
 2260 CONTINUE
C
      ITT=1
      DO 227 JLK=1,IKJ
      IF(ITT.EQ.1) THEN
         ITT=2
      ELSE
         ITT=1
      ENDIF
      IJKL=IKM1-JLK
      IKU2=(IJKL-1)*NG1P1+1
C
C
      IF(ITT.EQ.1) THEN
         CALL LWTTM(PGA(1,1,1,IJKL),PGB(1,1,1,IJKL)
     2             , PABCU(1,1,IKU2),PABCU(1,1,IKN),ZTT1)
      ELSE
         CALL LWTTM(PGA(1,1,1,IJKL),PGB(1,1,1,IJKL)
     2             , PABCU(1,1,IKU2),PABCU(1,1,IKN),ZTT2)
      ENDIF
C
      DO 2265 JA = 1, KTRAER
      DO 2264 JL = 1, KDLON
      ZTT(JL,JA) = (ZTT1(JL,JA)+ZTT2(JL,JA))*0.5
 2264 CONTINUE
 2265 CONTINUE
C
      DO 2266 JL = 1, KDLON
      ZWW=PDBDT(JL,1,IJKL)*ZTT(JL,1)          *ZTT(JL,10)
     S   +PDBDT(JL,2,IJKL)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     S   +PDBDT(JL,3,IJKL)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     S   +PDBDT(JL,4,IJKL)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     S   +PDBDT(JL,5,IJKL)*ZTT(JL,3)          *ZTT(JL,14)
     S   +PDBDT(JL,6,IJKL)*ZTT(JL,6)          *ZTT(JL,15)
      ZGLAYU(JL)=ZWW
      ZDZXMG=ZGLAYU(JL)
      PDISU(JL,JK)=PDISU(JL,JK)+ZDZXMG
      PCNTRB(JL,JK,IJKL)=ZDZXMG
 2266 CONTINUE
C
C
 227  CONTINUE
 228  CONTINUE
C
      RETURN
      END
      SUBROUTINE LWVN(KUAER,KTRAER
     R  , PABCU,PDBSL,PGA,PGB
     S  , PADJD,PADJU,PCNTRB,PDBDT)
      use dimens_m
      use dimphy
      use raddim
      IMPLICIT none
      include "raddimlw.h"
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           CARRIES OUT THE VERTICAL INTEGRATION ON NEARBY LAYERS
C           TO GIVE LONGWAVE FLUXES OR RADIANCES
C
C     METHOD.
C     -------
C
C          1. PERFORMS THE VERTICAL INTEGRATION CORRESPONDING TO THE
C     CONTRIBUTIONS OF THE ADJACENT LAYERS USING A GAUSSIAN QUADRATURE
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C-----------------------------------------------------------------------
C
C* ARGUMENTS:
C
      INTEGER KUAER,KTRAER
C
      REAL*8 PABCU(KDLON,NUA,3*KFLEV+1) ! ABSORBER AMOUNTS
      REAL*8 PDBSL(KDLON,Ninter,KFLEV*2) ! SUB-LAYER PLANCK FUNCTION GRADIENT
      REAL*8 PGA(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
      REAL*8 PGB(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
C
      REAL*8 PADJD(KDLON,KFLEV+1) ! CONTRIBUTION OF ADJACENT LAYERS
      REAL*8 PADJU(KDLON,KFLEV+1) ! CONTRIBUTION OF ADJACENT LAYERS
      REAL*8 PCNTRB(KDLON,KFLEV+1,KFLEV+1) ! CLEAR-SKY ENERGY EXCHANGE MATRIX
      REAL*8 PDBDT(KDLON,Ninter,KFLEV) !  LAYER PLANCK FUNCTION GRADIENT
C
C* LOCAL ARRAYS:
C
      REAL*8 ZGLAYD(KDLON)
      REAL*8 ZGLAYU(KDLON)
      REAL*8 ZTT(KDLON,NTRA)
      REAL*8 ZTT1(KDLON,NTRA)
      REAL*8 ZTT2(KDLON,NTRA)
      REAL*8 ZUU(KDLON,NUA)
C
      INTEGER jk, jl, ja, im12, ind, inu, ixu, jg
      INTEGER ixd, ibs, idd, imu, jk1, jk2, jnu
      REAL*8 zwtr
c
C* Data Block:
c
      REAL*8 WG1(2)
      SAVE WG1
      DATA (WG1(jk),jk=1,2) /1.0, 1.0/
C-----------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
 100  CONTINUE
C
C*         1.1     INITIALIZE LAYER CONTRIBUTIONS
C                  ------------------------------
C
 110  CONTINUE
C
      DO 112 JK = 1 , KFLEV+1
      DO 111 JL = 1, KDLON
      PADJD(JL,JK) = 0.
      PADJU(JL,JK) = 0.
 111  CONTINUE
 112  CONTINUE
C
C*         1.2     INITIALIZE TRANSMISSION FUNCTIONS
C                  ---------------------------------
C
 120  CONTINUE
C
      DO 122 JA = 1 , NTRA
      DO 121 JL = 1, KDLON
      ZTT (JL,JA) = 1.0
      ZTT1(JL,JA) = 1.0
      ZTT2(JL,JA) = 1.0
 121  CONTINUE
 122  CONTINUE
C
      DO 124 JA = 1 , NUA
      DO 123 JL = 1, KDLON
      ZUU(JL,JA) = 0.
 123  CONTINUE
 124  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.      VERTICAL INTEGRATION
C                  --------------------
C
 200  CONTINUE
C
C
C*         2.1     CONTRIBUTION FROM ADJACENT LAYERS
C                  ---------------------------------
C
 210  CONTINUE
C
      DO 215 JK = 1 , KFLEV
C
C*         2.1.1   DOWNWARD LAYERS
C                  ---------------
C
 2110 CONTINUE
C
      IM12 = 2 * (JK - 1)
      IND = (JK - 1) * NG1P1 + 1
      IXD = IND
      INU = JK * NG1P1 + 1
      IXU = IND
C
      DO 2111 JL = 1, KDLON
      ZGLAYD(JL) = 0.
      ZGLAYU(JL) = 0.
 2111 CONTINUE
C
      DO 213 JG = 1 , NG1
      IBS = IM12 + JG
      IDD = IXD + JG
      DO 2113 JA = 1 , KUAER
      DO 2112 JL = 1, KDLON
      ZUU(JL,JA) = PABCU(JL,JA,IND) - PABCU(JL,JA,IDD)
 2112 CONTINUE
 2113 CONTINUE
C
C
      CALL LWTT(PGA(1,1,1,JK), PGB(1,1,1,JK), ZUU, ZTT)
C
      DO 2114 JL = 1, KDLON
      ZWTR=PDBSL(JL,1,IBS)*ZTT(JL,1)          *ZTT(JL,10)
     S    +PDBSL(JL,2,IBS)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     S    +PDBSL(JL,3,IBS)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     S    +PDBSL(JL,4,IBS)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     S    +PDBSL(JL,5,IBS)*ZTT(JL,3)          *ZTT(JL,14)
     S    +PDBSL(JL,6,IBS)*ZTT(JL,6)          *ZTT(JL,15)
      ZGLAYD(JL)=ZGLAYD(JL)+ZWTR*WG1(JG)
 2114 CONTINUE
C
C*         2.1.2   DOWNWARD LAYERS
C                  ---------------
C
 2120 CONTINUE
C
      IMU = IXU + JG
      DO 2122 JA = 1 , KUAER
      DO 2121 JL = 1, KDLON
      ZUU(JL,JA) = PABCU(JL,JA,IMU) - PABCU(JL,JA,INU)
 2121 CONTINUE
 2122 CONTINUE
C
C
      CALL LWTT(PGA(1,1,1,JK), PGB(1,1,1,JK), ZUU, ZTT)
C
      DO 2123 JL = 1, KDLON
      ZWTR=PDBSL(JL,1,IBS)*ZTT(JL,1)          *ZTT(JL,10)
     S    +PDBSL(JL,2,IBS)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     S    +PDBSL(JL,3,IBS)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     S    +PDBSL(JL,4,IBS)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     S    +PDBSL(JL,5,IBS)*ZTT(JL,3)          *ZTT(JL,14)
     S    +PDBSL(JL,6,IBS)*ZTT(JL,6)          *ZTT(JL,15)
      ZGLAYU(JL)=ZGLAYU(JL)+ZWTR*WG1(JG)
 2123 CONTINUE
C
 213  CONTINUE
C
      DO 214 JL = 1, KDLON
      PADJD(JL,JK) = ZGLAYD(JL)
      PCNTRB(JL,JK,JK+1) = ZGLAYD(JL)
      PADJU(JL,JK+1) = ZGLAYU(JL)
      PCNTRB(JL,JK+1,JK) = ZGLAYU(JL)
      PCNTRB(JL,JK  ,JK) = 0.0
 214  CONTINUE
C
 215  CONTINUE
C
      DO 218 JK = 1 , KFLEV
      JK2 = 2 * JK
      JK1 = JK2 - 1
      DO 217 JNU = 1 , Ninter
      DO 216 JL = 1, KDLON
      PDBDT(JL,JNU,JK) = PDBSL(JL,JNU,JK1) + PDBSL(JL,JNU,JK2)
 216  CONTINUE
 217  CONTINUE
 218  CONTINUE
C
      RETURN
C
      END
      SUBROUTINE LWTT(PGA,PGB,PUU, PTT)
      use dimens_m
      use dimphy
      use raddim
      IMPLICIT none
      include "raddimlw.h"
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
C     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN ALL SIX SPECTRAL
C     INTERVALS.
C
C     METHOD.
C     -------
C
C          1. TRANSMISSION FUNCTION BY H2O AND UNIFORMLY MIXED GASES ARE
C     COMPUTED USING PADE APPROXIMANTS AND HORNER'S ALGORITHM.
C          2. TRANSMISSION BY O3 IS EVALUATED WITH MALKMUS'S BAND MODEL.
C          3. TRANSMISSION BY H2O CONTINUUM AND AEROSOLS FOLLOW AN
C     A SIMPLE EXPONENTIAL DECREASE WITH ABSORBER AMOUNT.
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 88-12-15
C
C-----------------------------------------------------------------------
      REAL*8 O1H, O2H
      PARAMETER (O1H=2230.)
      PARAMETER (O2H=100.)
      REAL*8 RPIALF0
      PARAMETER (RPIALF0=2.0)
C
C* ARGUMENTS:
C
      REAL*8 PUU(KDLON,NUA)
      REAL*8 PTT(KDLON,NTRA)
      REAL*8 PGA(KDLON,8,2)
      REAL*8 PGB(KDLON,8,2)
C
C* LOCAL VARIABLES:
C
      REAL*8 zz, zxd, zxn
      REAL*8 zpu, zpu10, zpu11, zpu12, zpu13
      REAL*8 zeu, zeu10, zeu11, zeu12, zeu13
      REAL*8 zx, zy, zsq1, zsq2, zvxy, zuxy
      REAL*8 zaercn, zto1, zto2, zxch4, zych4, zxn2o, zyn2o
      REAL*8 zsqn21, zodn21, zsqh42, zodh42
      REAL*8 zsqh41, zodh41, zsqn22, zodn22, zttf11, zttf12
      REAL*8 zuu11, zuu12, za11, za12
      INTEGER jl, ja
C     ------------------------------------------------------------------
C
C*         1.     HORNER'S ALGORITHM FOR H2O AND CO2 TRANSMISSION
C                 -----------------------------------------------
C
 100  CONTINUE
C
C
      DO 130 JA = 1 , 8
      DO 120 JL = 1, KDLON
      ZZ      =SQRT(PUU(JL,JA))
c     ZXD(JL,1)=PGB( JL, 1,1) + ZZ(JL, 1)*(PGB( JL, 1,2) + ZZ(JL, 1))
c     ZXN(JL,1)=PGA( JL, 1,1) + ZZ(JL, 1)*(PGA( JL, 1,2) )
c     PTT(JL,1)=ZXN(JL,1)/ZXD(JL,1)
      ZXD      =PGB( JL,JA,1) + ZZ       *(PGB( JL,JA,2) + ZZ       )
      ZXN      =PGA( JL,JA,1) + ZZ       *(PGA( JL,JA,2) )
      PTT(JL,JA)=ZXN      /ZXD
  120 CONTINUE
  130 CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.     CONTINUUM, OZONE AND AEROSOL TRANSMISSION FUNCTIONS
C                 ---------------------------------------------------
C
 200  CONTINUE
C
      DO 201 JL = 1, KDLON
      PTT(JL, 9) = PTT(JL, 8)
C
C-  CONTINUUM ABSORPTION: E- AND P-TYPE
C
      ZPU   = 0.002 * PUU(JL,10)
      ZPU10 = 112. * ZPU
      ZPU11 = 6.25 * ZPU
      ZPU12 = 5.00 * ZPU
      ZPU13 = 80.0 * ZPU
      ZEU   =  PUU(JL,11)
      ZEU10 =  12. * ZEU
      ZEU11 = 6.25 * ZEU
      ZEU12 = 5.00 * ZEU
      ZEU13 = 80.0 * ZEU
C
C-  OZONE ABSORPTION
C
      ZX = PUU(JL,12)
      ZY = PUU(JL,13)
      ZUXY = 4. * ZX * ZX / (RPIALF0 * ZY)
      ZSQ1 = SQRT(1. + O1H * ZUXY ) - 1.
      ZSQ2 = SQRT(1. + O2H * ZUXY ) - 1.
      ZVXY = RPIALF0 * ZY / (2. * ZX)
      ZAERCN = PUU(JL,17) + ZEU12 + ZPU12
      ZTO1 = EXP( - ZVXY * ZSQ1 - ZAERCN )
      ZTO2 = EXP( - ZVXY * ZSQ2 - ZAERCN )
C
C-- TRACE GASES (CH4, N2O, CFC-11, CFC-12)
C
C* CH4 IN INTERVAL 800-970 + 1110-1250 CM-1
C
c     NEXOTIC=1
c     IF (NEXOTIC.EQ.1) THEN
      ZXCH4 = PUU(JL,19)
      ZYCH4 = PUU(JL,20)
      ZUXY = 4. * ZXCH4*ZXCH4/(0.103*ZYCH4)
      ZSQH41 = SQRT(1. + 33.7 * ZUXY) - 1.
      ZVXY = 0.103 * ZYCH4 / (2. * ZXCH4)
      ZODH41 = ZVXY * ZSQH41
C
C* N2O IN INTERVAL 800-970 + 1110-1250 CM-1
C
      ZXN2O = PUU(JL,21)
      ZYN2O = PUU(JL,22)
      ZUXY = 4. * ZXN2O*ZXN2O/(0.416*ZYN2O)
      ZSQN21 = SQRT(1. + 21.3 * ZUXY) - 1.
      ZVXY = 0.416 * ZYN2O / (2. * ZXN2O)
      ZODN21 = ZVXY * ZSQN21
C
C* CH4 IN INTERVAL 1250-1450 + 1880-2820 CM-1
C
      ZUXY = 4. * ZXCH4*ZXCH4/(0.113*ZYCH4)
      ZSQH42 = SQRT(1. + 400. * ZUXY) - 1.
      ZVXY = 0.113 * ZYCH4 / (2. * ZXCH4)
      ZODH42 = ZVXY * ZSQH42
C
C* N2O IN INTERVAL 1250-1450 + 1880-2820 CM-1
C
      ZUXY = 4. * ZXN2O*ZXN2O/(0.197*ZYN2O)
      ZSQN22 = SQRT(1. + 2000. * ZUXY) - 1.
      ZVXY = 0.197 * ZYN2O / (2. * ZXN2O)
      ZODN22 = ZVXY * ZSQN22
C
C* CFC-11 IN INTERVAL 800-970 + 1110-1250 CM-1
C
      ZA11 = 2. * PUU(JL,23) * 4.404E+05
      ZTTF11 = 1. - ZA11 * 0.003225
C
C* CFC-12 IN INTERVAL 800-970 + 1110-1250 CM-1
C
      ZA12 = 2. * PUU(JL,24) * 6.7435E+05
      ZTTF12 = 1. - ZA12 * 0.003225
C
      ZUU11 = - PUU(JL,15) - ZEU10 - ZPU10
      ZUU12 = - PUU(JL,16) - ZEU11 - ZPU11 - ZODH41 - ZODN21
      PTT(JL,10) = EXP( - PUU(JL,14) )
      PTT(JL,11) = EXP( ZUU11 )
      PTT(JL,12) = EXP( ZUU12 ) * ZTTF11 * ZTTF12
      PTT(JL,13) = 0.7554 * ZTO1 + 0.2446 * ZTO2
      PTT(JL,14) = PTT(JL,10) * EXP( - ZEU13 - ZPU13 )
      PTT(JL,15) = EXP ( - PUU(JL,14) - ZODH42 - ZODN22 )
 201  CONTINUE
C
      RETURN
      END
      SUBROUTINE LWTTM(PGA,PGB,PUU1,PUU2, PTT)
      use dimens_m
      use dimphy
      use raddim
      IMPLICIT none
      include "raddimlw.h"
C
C     ------------------------------------------------------------------
C     PURPOSE.
C     --------
C           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
C     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN ALL SIX SPECTRAL
C     INTERVALS.
C
C     METHOD.
C     -------
C
C          1. TRANSMISSION FUNCTION BY H2O AND UNIFORMLY MIXED GASES ARE
C     COMPUTED USING PADE APPROXIMANTS AND HORNER'S ALGORITHM.
C          2. TRANSMISSION BY O3 IS EVALUATED WITH MALKMUS'S BAND MODEL.
C          3. TRANSMISSION BY H2O CONTINUUM AND AEROSOLS FOLLOW AN
C     A SIMPLE EXPONENTIAL DECREASE WITH ABSORBER AMOUNT.
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 88-12-15
C
C-----------------------------------------------------------------------
      REAL*8 O1H, O2H
      PARAMETER (O1H=2230.)
      PARAMETER (O2H=100.)
      REAL*8 RPIALF0
      PARAMETER (RPIALF0=2.0)
C
C* ARGUMENTS:
C
      REAL*8 PGA(KDLON,8,2) ! PADE APPROXIMANTS
      REAL*8 PGB(KDLON,8,2) ! PADE APPROXIMANTS
      REAL*8 PUU1(KDLON,NUA) ! ABSORBER AMOUNTS FROM TOP TO LEVEL 1
      REAL*8 PUU2(KDLON,NUA) ! ABSORBER AMOUNTS FROM TOP TO LEVEL 2
      REAL*8 PTT(KDLON,NTRA) ! TRANSMISSION FUNCTIONS
C
C* LOCAL VARIABLES:
C
      INTEGER ja, jl
      REAL*8 zz, zxd, zxn
      REAL*8 zpu, zpu10, zpu11, zpu12, zpu13
      REAL*8 zeu, zeu10, zeu11, zeu12, zeu13
      REAL*8 zx, zy, zuxy, zsq1, zsq2, zvxy, zaercn, zto1, zto2
      REAL*8 zxch4, zych4, zsqh41, zodh41
      REAL*8 zxn2o, zyn2o, zsqn21, zodn21, zsqh42, zodh42
      REAL*8 zsqn22, zodn22, za11, zttf11, za12, zttf12
      REAL*8 zuu11, zuu12
C     ------------------------------------------------------------------
C
C*         1.     HORNER'S ALGORITHM FOR H2O AND CO2 TRANSMISSION
C                 -----------------------------------------------
C
 100  CONTINUE
C
C
      DO 130 JA = 1 , 8
      DO 120 JL = 1, KDLON
      ZZ      =SQRT(PUU1(JL,JA) - PUU2(JL,JA))
      ZXD      =PGB( JL,JA,1) + ZZ       *(PGB( JL,JA,2) + ZZ       )
      ZXN      =PGA( JL,JA,1) + ZZ       *(PGA( JL,JA,2) )
      PTT(JL,JA)=ZXN      /ZXD
  120 CONTINUE
  130 CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.     CONTINUUM, OZONE AND AEROSOL TRANSMISSION FUNCTIONS
C                 ---------------------------------------------------
C
 200  CONTINUE
C
      DO 201 JL = 1, KDLON
      PTT(JL, 9) = PTT(JL, 8)
C
C-  CONTINUUM ABSORPTION: E- AND P-TYPE
C
      ZPU   = 0.002 * (PUU1(JL,10) - PUU2(JL,10))
      ZPU10 = 112. * ZPU
      ZPU11 = 6.25 * ZPU
      ZPU12 = 5.00 * ZPU
      ZPU13 = 80.0 * ZPU
      ZEU   = (PUU1(JL,11) - PUU2(JL,11))
      ZEU10 =  12. * ZEU
      ZEU11 = 6.25 * ZEU
      ZEU12 = 5.00 * ZEU
      ZEU13 = 80.0 * ZEU
C
C-  OZONE ABSORPTION
C
      ZX = (PUU1(JL,12) - PUU2(JL,12))
      ZY = (PUU1(JL,13) - PUU2(JL,13))
      ZUXY = 4. * ZX * ZX / (RPIALF0 * ZY)
      ZSQ1 = SQRT(1. + O1H * ZUXY ) - 1.
      ZSQ2 = SQRT(1. + O2H * ZUXY ) - 1.
      ZVXY = RPIALF0 * ZY / (2. * ZX)
      ZAERCN = (PUU1(JL,17) -PUU2(JL,17)) + ZEU12 + ZPU12
      ZTO1 = EXP( - ZVXY * ZSQ1 - ZAERCN )
      ZTO2 = EXP( - ZVXY * ZSQ2 - ZAERCN )
C
C-- TRACE GASES (CH4, N2O, CFC-11, CFC-12)
C
C* CH4 IN INTERVAL 800-970 + 1110-1250 CM-1
C
      ZXCH4 = (PUU1(JL,19) - PUU2(JL,19))
      ZYCH4 = (PUU1(JL,20) - PUU2(JL,20))
      ZUXY = 4. * ZXCH4*ZXCH4/(0.103*ZYCH4)
      ZSQH41 = SQRT(1. + 33.7 * ZUXY) - 1.
      ZVXY = 0.103 * ZYCH4 / (2. * ZXCH4)
      ZODH41 = ZVXY * ZSQH41
C
C* N2O IN INTERVAL 800-970 + 1110-1250 CM-1
C
      ZXN2O = (PUU1(JL,21) - PUU2(JL,21))
      ZYN2O = (PUU1(JL,22) - PUU2(JL,22))
      ZUXY = 4. * ZXN2O*ZXN2O/(0.416*ZYN2O)
      ZSQN21 = SQRT(1. + 21.3 * ZUXY) - 1.
      ZVXY = 0.416 * ZYN2O / (2. * ZXN2O)
      ZODN21 = ZVXY * ZSQN21
C
C* CH4 IN INTERVAL 1250-1450 + 1880-2820 CM-1
C
      ZUXY = 4. * ZXCH4*ZXCH4/(0.113*ZYCH4)
      ZSQH42 = SQRT(1. + 400. * ZUXY) - 1.
      ZVXY = 0.113 * ZYCH4 / (2. * ZXCH4)
      ZODH42 = ZVXY * ZSQH42
C
C* N2O IN INTERVAL 1250-1450 + 1880-2820 CM-1
C
      ZUXY = 4. * ZXN2O*ZXN2O/(0.197*ZYN2O)
      ZSQN22 = SQRT(1. + 2000. * ZUXY) - 1.
      ZVXY = 0.197 * ZYN2O / (2. * ZXN2O)
      ZODN22 = ZVXY * ZSQN22
C
C* CFC-11 IN INTERVAL 800-970 + 1110-1250 CM-1
C
      ZA11 = (PUU1(JL,23) - PUU2(JL,23)) * 4.404E+05
      ZTTF11 = 1. - ZA11 * 0.003225
C
C* CFC-12 IN INTERVAL 800-970 + 1110-1250 CM-1
C
      ZA12 = (PUU1(JL,24) - PUU2(JL,24)) * 6.7435E+05
      ZTTF12 = 1. - ZA12 * 0.003225
C
      ZUU11 = - (PUU1(JL,15) - PUU2(JL,15)) - ZEU10 - ZPU10
      ZUU12 = - (PUU1(JL,16) - PUU2(JL,16)) - ZEU11 - ZPU11 -
     S         ZODH41 - ZODN21
      PTT(JL,10) = EXP( - (PUU1(JL,14)- PUU2(JL,14)) )
      PTT(JL,11) = EXP( ZUU11 )
      PTT(JL,12) = EXP( ZUU12 ) * ZTTF11 * ZTTF12
      PTT(JL,13) = 0.7554 * ZTO1 + 0.2446 * ZTO2
      PTT(JL,14) = PTT(JL,10) * EXP( - ZEU13 - ZPU13 )
      PTT(JL,15) = EXP ( - (PUU1(JL,14) - PUU2(JL,14)) - ZODH42-ZODN22 )
 201  CONTINUE
C
      RETURN
      END
