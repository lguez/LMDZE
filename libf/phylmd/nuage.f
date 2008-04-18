!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/nuage.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
      SUBROUTINE nuage (paprs, pplay,
     .                  t, pqlwp, pclc, pcltau, pclemi,
     .                  pch, pcl, pcm, pct, pctlwp,
     e                  ok_aie,
     e                  sulfate, sulfate_pi, 
     e                  bl95_b0, bl95_b1,
     s                  cldtaupi, re, fl)
      use dimens_m
      use dimphy
      use YOMCST
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19930910
c Objet: Calculer epaisseur optique et emmissivite des nuages
c======================================================================
c Arguments:
c t-------input-R-temperature
c pqlwp---input-R-eau liquide nuageuse dans l'atmosphere (kg/kg)
c pclc----input-R-couverture nuageuse pour le rayonnement (0 a 1)
c ok_aie--input-L-apply aerosol indirect effect or not
c sulfate-input-R-sulfate aerosol mass concentration [um/m^3]
c sulfate_pi-input-R-dito, pre-industrial value
c bl95_b0-input-R-a parameter, may be varied for tests (s-sea, l-land)
c bl95_b1-input-R-a parameter, may be varied for tests (    -"-      )
c      
c cldtaupi-output-R-pre-industrial value of cloud optical thickness, 
c                   needed for the diagnostics of the aerosol indirect 
c                   radiative forcing (see radlwsw)
c re------output-R-Cloud droplet effective radius multiplied by fl [um]
c fl------output-R-Denominator to re, introduced to avoid problems in
c                  the averaging of the output. fl is the fraction of liquid
c                  water clouds within a grid cell      
c 
c pcltau--output-R-epaisseur optique des nuages
c pclemi--output-R-emissivite des nuages (0 a 1)
c======================================================================
C
c
      REAL, intent(in):: paprs(klon,klev+1)
      real, intent(in):: pplay(klon,klev)
      REAL t(klon,klev)
c
      REAL pclc(klon,klev)
      REAL pqlwp(klon,klev)
      REAL pcltau(klon,klev), pclemi(klon,klev)
c
      REAL pct(klon), pctlwp(klon), pch(klon), pcl(klon), pcm(klon)
c
      LOGICAL lo
c
      REAL cetahb, cetamb
      PARAMETER (cetahb = 0.45, cetamb = 0.80)
C
      INTEGER i, k
      REAL zflwp, zradef, zfice, zmsac
c
      REAL radius, rad_froid, rad_chaud, rad_chau1, rad_chau2
      PARAMETER (rad_chau1=13.0, rad_chau2=9.0, rad_froid=35.0)
ccc      PARAMETER (rad_chaud=15.0, rad_froid=35.0)
c sintex initial      PARAMETER (rad_chaud=10.0, rad_froid=30.0)
      REAL coef, coef_froi, coef_chau
      PARAMETER (coef_chau=0.13, coef_froi=0.09)
      REAL seuil_neb, t_glace
      PARAMETER (seuil_neb=0.001, t_glace=273.0-15.0)
      INTEGER nexpo ! exponentiel pour glace/eau
      PARAMETER (nexpo=6)
      
cjq for the aerosol indirect effect
cjq introduced by Johannes Quaas (quaas@lmd.jussieu.fr), 27/11/2003
cjq      
      LOGICAL ok_aie            ! Apply AIE or not?
      
      REAL sulfate(klon, klev)  ! sulfate aerosol mass concentration [ug m-3]
      REAL cdnc(klon, klev)     ! cloud droplet number concentration [m-3]
      REAL re(klon, klev)       ! cloud droplet effective radius [um]
      REAL sulfate_pi(klon, klev)  ! sulfate aerosol mass concentration [ug m-3] (pre-industrial value)
      REAL cdnc_pi(klon, klev)     ! cloud droplet number concentration [m-3] (pi value)
      REAL re_pi(klon, klev)       ! cloud droplet effective radius [um] (pi value)
      
      REAL fl(klon, klev)  ! xliq * rneb (denominator to re; fraction of liquid water clouds within the grid cell)
      
      REAL bl95_b0, bl95_b1     ! Parameter in B&L 95-Formula
      
      REAL cldtaupi(klon, klev) ! pre-industrial cloud opt thickness for diag
cjq-end      
      
ccc      PARAMETER (nexpo=1)
c
c Calculer l'epaisseur optique et l'emmissivite des nuages
c
      DO k = 1, klev
      DO i = 1, klon
         rad_chaud = rad_chau1
         IF (k.LE.3) rad_chaud = rad_chau2
            
         pclc(i,k) = MAX(pclc(i,k), seuil_neb)
         zflwp = 1000.*pqlwp(i,k)/RG/pclc(i,k)
     .          *(paprs(i,k)-paprs(i,k+1))
         zfice = 1.0 - (t(i,k)-t_glace) / (273.13-t_glace)
         zfice = MIN(MAX(zfice,0.0),1.0)
         zfice = zfice**nexpo
         
         IF (ok_aie) THEN
            ! Formula "D" of Boucher and Lohmann, Tellus, 1995
            !             
            cdnc(i,k) = 10.**(bl95_b0+bl95_b1*
     .           log(MAX(sulfate(i,k),1.e-4))/log(10.))*1.e6 !-m-3
            ! Cloud droplet number concentration (CDNC) is restricted
            ! to be within [20, 1000 cm^3]
            ! 
            cdnc(i,k)=MIN(1000.e6,MAX(20.e6,cdnc(i,k)))
            cdnc_pi(i,k) = 10.**(bl95_b0+bl95_b1*
     .           log(MAX(sulfate_pi(i,k),1.e-4))/log(10.))*1.e6 !-m-3
            cdnc_pi(i,k)=MIN(1000.e6,MAX(20.e6,cdnc_pi(i,k)))
            !            
            !
            ! air density: pplay(i,k) / (RD * zT(i,k)) 
            ! factor 1.1: derive effective radius from volume-mean radius
            ! factor 1000 is the water density
            ! _chaud means that this is the CDR for liquid water clouds
            !
            rad_chaud = 
     .           1.1 * ( (pqlwp(i,k) * pplay(i,k) / (RD * T(i,k)) )  
     .               / (4./3. * RPI * 1000. * cdnc(i,k)) )**(1./3.)
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
            radius = MAX(1.1e6 * ( (pqlwp(i,k)*pplay(i,k)/(RD*T(i,k)))  
     .                / (4./3.*RPI*1000.*cdnc_pi(i,k)) )**(1./3.), 
     .               3.) * (1.-zfice) + rad_froid * zfice           
            cldtaupi(i,k) = 3.0/2.0 * zflwp / radius
     .           
         ENDIF                  ! ok_aie
         
         radius = rad_chaud * (1.-zfice) + rad_froid * zfice
         coef = coef_chau * (1.-zfice) + coef_froi * zfice
         pcltau(i,k) = 3.0/2.0 * zflwp / radius
         pclemi(i,k) = 1.0 - EXP( - coef * zflwp)
         lo = (pclc(i,k) .LE. seuil_neb)
         IF (lo) pclc(i,k) = 0.0
         IF (lo) pcltau(i,k) = 0.0
         IF (lo) pclemi(i,k) = 0.0
         
         IF (.NOT.ok_aie) cldtaupi(i,k)=pcltau(i,k)            
      ENDDO
      ENDDO
ccc      DO k = 1, klev
ccc      DO i = 1, klon
ccc         t(i,k) = t(i,k)
ccc         pclc(i,k) = MAX( 1.e-5 , pclc(i,k) )
ccc         lo = pclc(i,k) .GT. (2.*1.e-5)
ccc         zflwp = pqlwp(i,k)*1000.*(paprs(i,k)-paprs(i,k+1))
ccc     .          /(rg*pclc(i,k))
ccc         zradef = 10.0 + (1.-sigs(k))*45.0
ccc         pcltau(i,k) = 1.5 * zflwp / zradef
ccc         zfice=1.0-MIN(MAX((t(i,k)-263.)/(273.-263.),0.0),1.0)
ccc         zmsac = 0.13*(1.0-zfice) + 0.08*zfice
ccc         pclemi(i,k) = 1.-EXP(-zmsac*zflwp)
ccc         if (.NOT.lo) pclc(i,k) = 0.0
ccc         if (.NOT.lo) pcltau(i,k) = 0.0
ccc         if (.NOT.lo) pclemi(i,k) = 0.0
ccc      ENDDO
ccc      ENDDO
cccccc      print*, 'pas de nuage dans le rayonnement'
cccccc      DO k = 1, klev
cccccc      DO i = 1, klon
cccccc         pclc(i,k) = 0.0
cccccc         pcltau(i,k) = 0.0
cccccc         pclemi(i,k) = 0.0
cccccc      ENDDO
cccccc      ENDDO
C
C COMPUTE CLOUD LIQUID PATH AND TOTAL CLOUDINESS
C
      DO i = 1, klon
         pct(i)=1.0
         pch(i)=1.0
         pcm(i) = 1.0
         pcl(i) = 1.0
         pctlwp(i) = 0.0
      ENDDO
C
      DO k = klev, 1, -1
      DO i = 1, klon
         pctlwp(i) = pctlwp(i) 
     .             + pqlwp(i,k)*(paprs(i,k)-paprs(i,k+1))/RG
         pct(i) = pct(i)*(1.0-pclc(i,k))
         if (pplay(i,k).LE.cetahb*paprs(i,1))
     .      pch(i) = pch(i)*(1.0-pclc(i,k))
         if (pplay(i,k).GT.cetahb*paprs(i,1) .AND.
     .       pplay(i,k).LE.cetamb*paprs(i,1)) 
     .      pcm(i) = pcm(i)*(1.0-pclc(i,k))
         if (pplay(i,k).GT.cetamb*paprs(i,1))
     .      pcl(i) = pcl(i)*(1.0-pclc(i,k))
      ENDDO
      ENDDO
C
      DO i = 1, klon
         pct(i)=1.-pct(i)
         pch(i)=1.-pch(i)
         pcm(i)=1.-pcm(i)
         pcl(i)=1.-pcl(i)
      ENDDO
C
      RETURN
      END
      SUBROUTINE diagcld1(paprs,pplay,rain,snow,kbot,ktop,
     .                   diafra,dialiq)
      use dimens_m
      use dimphy
      use YOMCST
      IMPLICIT none
c
c Laurent Li (LMD/CNRS), le 12 octobre 1998
c                        (adaptation du code ECMWF)
c
c Dans certains cas, le schema pronostique des nuages n'est
c pas suffisament performant. On a donc besoin de diagnostiquer
c ces nuages. Je dois avouer que c'est une frustration.
c
c
c Arguments d'entree:
      REAL, intent(in):: paprs(klon,klev+1) ! pression (Pa) a inter-couche
      REAL, intent(in):: pplay(klon,klev) ! pression (Pa) au milieu de couche
      REAL t(klon,klev) ! temperature (K)
      REAL q(klon,klev) ! humidite specifique (Kg/Kg)
      REAL rain(klon) ! pluie convective (kg/m2/s)
      REAL snow(klon) ! neige convective (kg/m2/s)
      INTEGER ktop(klon) ! sommet de la convection
      INTEGER kbot(klon) ! bas de la convection
c
c Arguments de sortie:
      REAL diafra(klon,klev) ! fraction nuageuse diagnostiquee
      REAL dialiq(klon,klev) ! eau liquide nuageuse
c
c Constantes ajustables:
      REAL CANVA, CANVB, CANVH
      PARAMETER (CANVA=2.0, CANVB=0.3, CANVH=0.4)
      REAL CCA, CCB, CCC
      PARAMETER (CCA=0.125, CCB=1.5, CCC=0.8)
      REAL CCFCT, CCSCAL
      PARAMETER (CCFCT=0.400)
      PARAMETER (CCSCAL=1.0E+11)
      REAL CETAHB, CETAMB
      PARAMETER (CETAHB=0.45, CETAMB=0.80)
      REAL CCLWMR
      PARAMETER (CCLWMR=1.E-04)
      REAL ZEPSCR
      PARAMETER (ZEPSCR=1.0E-10)
c
c Variables locales:
      INTEGER i, k
      REAL zcc(klon)
c
c Initialisation:
c
      DO k = 1, klev
      DO i = 1, klon
         diafra(i,k) = 0.0
         dialiq(i,k) = 0.0
      ENDDO
      ENDDO
c
      DO i = 1, klon ! Calculer la fraction nuageuse
      zcc(i) = 0.0
      IF((rain(i)+snow(i)).GT.0.) THEN
         zcc(i)= CCA * LOG(MAX(ZEPSCR,(rain(i)+snow(i))*CCSCAL))-CCB
         zcc(i)= MIN(CCC,MAX(0.0,zcc(i)))
      ENDIF
      ENDDO
c
      DO i = 1, klon ! pour traiter les enclumes
      diafra(i,ktop(i)) = MAX(diafra(i,ktop(i)),zcc(i)*CCFCT)
      IF ((zcc(i).GE.CANVH) .AND.
     .    (pplay(i,ktop(i)).LE.CETAHB*paprs(i,1)))
     . diafra(i,ktop(i)) = MAX(diafra(i,ktop(i)),
     .                         MAX(zcc(i)*CCFCT,CANVA*(zcc(i)-CANVB)))
      dialiq(i,ktop(i))=CCLWMR*diafra(i,ktop(i))
      ENDDO
c
      DO k = 1, klev ! nuages convectifs (sauf enclumes)
      DO i = 1, klon
      IF (k.LT.ktop(i) .AND. k.GE.kbot(i)) THEN
         diafra(i,k)=MAX(diafra(i,k),zcc(i)*CCFCT)
         dialiq(i,k)=CCLWMR*diafra(i,k)
      ENDIF
      ENDDO
      ENDDO
c
      RETURN
      END
      SUBROUTINE diagcld2(paprs,pplay,t,q, diafra,dialiq)
      use dimens_m
      use dimphy
      use YOMCST
      use yoethf
c Fonctions thermodynamiques:
      use fcttre
      IMPLICIT none
c
c
c Arguments d'entree:
      REAL, intent(in):: paprs(klon,klev+1) ! pression (Pa) a inter-couche
      REAL, intent(in):: pplay(klon,klev) ! pression (Pa) au milieu de couche
      REAL t(klon,klev) ! temperature (K)
      REAL q(klon,klev) ! humidite specifique (Kg/Kg)
c
c Arguments de sortie:
      REAL diafra(klon,klev) ! fraction nuageuse diagnostiquee
      REAL dialiq(klon,klev) ! eau liquide nuageuse
c
      REAL CETAMB
      PARAMETER (CETAMB=0.80)
      REAL CLOIA, CLOIB, CLOIC, CLOID
      PARAMETER (CLOIA=1.0E+02, CLOIB=-10.00, CLOIC=-0.6, CLOID=5.0)
ccc      PARAMETER (CLOIA=1.0E+02, CLOIB=-10.00, CLOIC=-0.9, CLOID=5.0)
      REAL RGAMMAS
      PARAMETER (RGAMMAS=0.05)
      REAL CRHL
      PARAMETER (CRHL=0.15)
ccc      PARAMETER (CRHL=0.70)
      REAL t_coup
      PARAMETER (t_coup=234.0)
c
c Variables locales:
      INTEGER i, k, kb, invb(klon)
      REAL zqs, zrhb, zcll, zdthmin(klon), zdthdp
      REAL zdelta, zcor
c
c
c Initialisation:
c
      DO k = 1, klev
      DO i = 1, klon
         diafra(i,k) = 0.0
         dialiq(i,k) = 0.0
      ENDDO
      ENDDO
c
      DO i = 1, klon
         invb(i) = klev
         zdthmin(i)=0.0
      ENDDO

      DO k = 2, klev/2-1
      DO i = 1, klon
         zdthdp = (t(i,k)-t(i,k+1))/(pplay(i,k)-pplay(i,k+1))
     .          - RD * 0.5*(t(i,k)+t(i,k+1))/RCPD/paprs(i,k+1)
         zdthdp = zdthdp * CLOIA
         IF (pplay(i,k).GT.CETAMB*paprs(i,1) .AND.
     .       zdthdp.LT.zdthmin(i) ) THEN
            zdthmin(i) = zdthdp
            invb(i) = k
         ENDIF
      ENDDO
      ENDDO

      DO i = 1, klon
         kb=invb(i)
         IF (thermcep) THEN
           zdelta=MAX(0.,SIGN(1.,RTT-t(i,kb)))
           zqs= R2ES*FOEEW(t(i,kb),zdelta)/pplay(i,kb)
           zqs=MIN(0.5,zqs)
           zcor=1./(1.-RETV*zqs)
           zqs=zqs*zcor
         ELSE
           IF (t(i,kb) .LT. t_coup) THEN
              zqs = qsats(t(i,kb)) / pplay(i,kb)
           ELSE
              zqs = qsatl(t(i,kb)) / pplay(i,kb)
           ENDIF
         ENDIF
         zcll = CLOIB * zdthmin(i) + CLOIC
         zcll = MIN(1.0,MAX(0.0,zcll))
         zrhb= q(i,kb)/zqs
         IF (zcll.GT.0.0.AND.zrhb.LT.CRHL)
     .   zcll=zcll*(1.-(CRHL-zrhb)*CLOID)
         zcll=MIN(1.0,MAX(0.0,zcll))
         diafra(i,kb) = MAX(diafra(i,kb),zcll)
         dialiq(i,kb)= diafra(i,kb) * RGAMMAS*zqs
      ENDDO
c
      RETURN
      END
