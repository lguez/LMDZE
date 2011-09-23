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
      use SUPHEC_M
      use raddim, only: kflev, kdlon
      use yoethf_m
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
      real, intent(in):: pplay(klon,klev)
      real albedo(klon), alblw(klon), tsol(klon)
      real, intent(in):: t(klon,klev)
      real q(klon,klev)
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
