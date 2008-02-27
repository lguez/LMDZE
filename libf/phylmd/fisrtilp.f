!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/fisrtilp.F,v 1.2 2004/11/09 16:55:40 lmdzadmin Exp $
!
c
      SUBROUTINE fisrtilp(dtime,paprs,pplay,t,q,ptconv,ratqs,
     s                   d_t, d_q, d_ql, rneb, radliq, rain, snow,
     s                   pfrac_impa, pfrac_nucl, pfrac_1nucl,
     s                   frac_impa, frac_nucl,
     s                   prfl, psfl, rhcl)

c
      use dimens_m
      use dimphy
      use tracstoke
      use YOMCST
      use yoethf
      use fcttre
      use comfisrtilp
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS)
c Date: le 20 mars 1995
c Objet: condensation et precipitation stratiforme.
c        schema de nuage
c======================================================================
c======================================================================
c
c Arguments:
c
      REAL dtime ! intervalle du temps (s)
      REAL, intent(in):: paprs(klon,klev+1) ! pression a inter-couche
      REAL pplay(klon,klev) ! pression au milieu de couche
      REAL t(klon,klev) ! temperature (K)
      REAL q(klon,klev) ! humidite specifique (kg/kg)
      REAL d_t(klon,klev) ! incrementation de la temperature (K)
      REAL d_q(klon,klev) ! incrementation de la vapeur d'eau
      REAL d_ql(klon,klev) ! incrementation de l'eau liquide
      REAL rneb(klon,klev) ! fraction nuageuse
      REAL radliq(klon,klev) ! eau liquide utilisee dans rayonnements
      REAL rhcl(klon,klev) ! humidite relative en ciel clair
      REAL rain(klon) ! pluies (mm/s)
      REAL snow(klon) ! neige (mm/s)
      REAL prfl(klon,klev+1) ! flux d'eau precipitante aux interfaces (kg/m2/s)
      REAL psfl(klon,klev+1) ! flux d'eau precipitante aux interfaces (kg/m2/s)
cAA
c Coeffients de fraction lessivee : pour OFF-LINE
c
      REAL pfrac_nucl(klon,klev)
      REAL pfrac_1nucl(klon,klev)
      REAL pfrac_impa(klon,klev)
c
c Fraction d'aerosols lessivee par impaction et par nucleation
c POur ON-LINE
c
      REAL frac_impa(klon,klev)
      REAL frac_nucl(klon,klev)
      real zct(klon),zcl(klon)
cAA
c
c Options du programme:
c
      REAL seuil_neb ! un nuage existe vraiment au-dela
      PARAMETER (seuil_neb=0.001)

      INTEGER ninter ! sous-intervals pour la precipitation
      PARAMETER (ninter=5)
      LOGICAL evap_prec ! evaporation de la pluie
      PARAMETER (evap_prec=.TRUE.)
      REAL ratqs(klon,klev) ! determine la largeur de distribution de vapeur
      logical ptconv(klon,klev) ! determine la largeur de distribution de vapeur

      real zpdf_sig(klon),zpdf_k(klon),zpdf_delta(klon)
      real Zpdf_a(klon),zpdf_b(klon),zpdf_e1(klon),zpdf_e2(klon)
      real erf
c
      LOGICAL cpartiel ! condensation partielle
      PARAMETER (cpartiel=.TRUE.)
      REAL t_coup
      PARAMETER (t_coup=234.0)
c
c Variables locales:
c
      INTEGER i, k, n, kk
      REAL zqs(klon), zdqs(klon), zdelta, zcor, zcvm5
      REAL zrfl(klon), zrfln(klon), zqev, zqevt
      REAL zoliq(klon), zcond(klon), zq(klon), zqn(klon), zdelq
      REAL ztglace, zt(klon)
      INTEGER nexpo ! exponentiel pour glace/eau
      REAL zdz(klon),zrho(klon),ztot(klon), zrhol(klon)
      REAL zchau(klon),zfroi(klon),zfice(klon),zneb(klon)
c
      LOGICAL appel1er
      SAVE appel1er
c
c---------------------------------------------------------------
c
cAA Variables traceurs:
cAA  Provisoire !!! Parametres alpha du lessivage
cAA  A priori on a 4 scavenging # possibles
c
      REAL a_tr_sca(4)
      save a_tr_sca
c
c Variables intermediaires
c
      REAL zalpha_tr
      REAL zfrac_lessi
      REAL zprec_cond(klon)
cAA
      REAL zmair, zcpair, zcpeau
C     Pour la conversion eau-neige
      REAL zlh_solid(klon), zm_solid
cIM 
      INTEGER klevm1
c---------------------------------------------------------------
c
c Fonctions en ligne:
c
      REAL fallvs,fallvc ! vitesse de chute pour crystaux de glace
      REAL zzz
      fallvc (zzz) = 3.29/2.0 * ((zzz)**0.16) * ffallv_con
      fallvs (zzz) = 3.29/2.0 * ((zzz)**0.16) * ffallv_lsc
c
      DATA appel1er /.TRUE./
cym
      zdelq=0.0
      
      IF (appel1er) THEN
c
         PRINT*, 'fisrtilp, ninter:', ninter
         PRINT*, 'fisrtilp, evap_prec:', evap_prec
         PRINT*, 'fisrtilp, cpartiel:', cpartiel
         IF (ABS(dtime/FLOAT(ninter)-360.0).GT.0.001) THEN
          PRINT*, 'fisrtilp: Ce n est pas prevu, voir Z.X.Li', dtime
          PRINT*, 'Je prefere un sous-intervalle de 6 minutes'
c         stop 1
         ENDIF
         appel1er = .FALSE.
c
cAA initialiation provisoire
       a_tr_sca(1) = -0.5
       a_tr_sca(2) = -0.5
       a_tr_sca(3) = -0.5
       a_tr_sca(4) = -0.5
c
cAA Initialisation a 1 des coefs des fractions lessivees 
c
      DO k = 1, klev
       DO i = 1, klon
          pfrac_nucl(i,k)=1.
          pfrac_1nucl(i,k)=1.
          pfrac_impa(i,k)=1.
       ENDDO 
      ENDDO 

      ENDIF          !  test sur appel1er
c
cMAf Initialisation a 0 de zoliq
       DO i = 1, klon
          zoliq(i)=0.
       ENDDO 
c Determiner les nuages froids par leur temperature
c  nexpo regle la raideur de la transition eau liquide / eau glace.
c
      ztglace = RTT - 15.0
      nexpo = 6
ccc      nexpo = 1
c
c Initialiser les sorties:
c
      DO k = 1, klev+1
      DO i = 1, klon
         prfl(i,k) = 0.0
         psfl(i,k) = 0.0
      ENDDO
      ENDDO

      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = 0.0
         d_q(i,k) = 0.0
         d_ql(i,k) = 0.0
         rneb(i,k) = 0.0
         radliq(i,k) = 0.0
         frac_nucl(i,k) = 1. 
         frac_impa(i,k) = 1. 
      ENDDO
      ENDDO
      DO i = 1, klon
         rain(i) = 0.0
         snow(i) = 0.0
      ENDDO
c
c Initialiser le flux de precipitation a zero
c
      DO i = 1, klon
         zrfl(i) = 0.0
         zneb(i) = seuil_neb
      ENDDO
c
c
cAA Pour plus de securite 

      zalpha_tr   = 0.
      zfrac_lessi = 0.

cAA----------------------------------------------------------
c
c Boucle verticale (du haut vers le bas)
c
cIM : klevm1
      klevm1=klev-1
      DO 9999 k = klev, 1, -1
c
cAA----------------------------------------------------------
c
      DO i = 1, klon
         zt(i)=t(i,k)
         zq(i)=q(i,k)
      ENDDO
c
c Calculer la varition de temp. de l'air du a la chaleur sensible
C transporter par la pluie.
C Il resterait a rajouter cet effet de la chaleur sensible sur les
C flux de surface, du a la diff. de temp. entre le 1er niveau et la
C surface.
C
      DO i = 1, klon
cIM
       IF(k.LE.klevm1) THEN         
        zmair=(paprs(i,k)-paprs(i,k+1))/RG
        zcpair=RCPD*(1.0+RVTMP2*zq(i))
        zcpeau=RCPD*RVTMP2
        zt(i) = ( (t(i,k+1)+d_t(i,k+1))*zrfl(i)*dtime*zcpeau
     $      + zmair*zcpair*zt(i) )
     $      / (zmair*zcpair + zrfl(i)*dtime*zcpeau)
CC        WRITE (6,*) 'cppluie ', zt(i)-(t(i,k+1)+d_t(i,k+1))
       ENDIF
      ENDDO
c
c
c Calculer l'evaporation de la precipitation
c


      IF (evap_prec) THEN
      DO i = 1, klon
      IF (zrfl(i) .GT.0.) THEN
         IF (thermcep) THEN
           zdelta=MAX(0.,SIGN(1.,RTT-zt(i)))
           zqs(i)= R2ES*FOEEW(zt(i),zdelta)/pplay(i,k)
           zqs(i)=MIN(0.5,zqs(i))
           zcor=1./(1.-RETV*zqs(i))
           zqs(i)=zqs(i)*zcor
         ELSE
           IF (zt(i) .LT. t_coup) THEN
              zqs(i) = qsats(zt(i)) / pplay(i,k)
           ELSE
              zqs(i) = qsatl(zt(i)) / pplay(i,k)
           ENDIF
         ENDIF
         zqev = MAX (0.0, (zqs(i)-zq(i))*zneb(i) )
         zqevt = coef_eva * (1.0-zq(i)/zqs(i)) * SQRT(zrfl(i))
     .         * (paprs(i,k)-paprs(i,k+1))/pplay(i,k)*zt(i)*RD/RG
         zqevt = MAX(0.0,MIN(zqevt,zrfl(i)))
     .         * RG*dtime/(paprs(i,k)-paprs(i,k+1))
         zqev = MIN (zqev, zqevt)
         zrfln(i) = zrfl(i) - zqev*(paprs(i,k)-paprs(i,k+1))
     .                            /RG/dtime

c pour la glace, on réévapore toute la précip dans la couche du dessous
c la glace venant de la couche du dessus est simplement dans la couche
c du dessous.

         IF (zt(i) .LT. t_coup.and.reevap_ice) zrfln(i)=0.

         zq(i) = zq(i) - (zrfln(i)-zrfl(i))
     .             * (RG/(paprs(i,k)-paprs(i,k+1)))*dtime
         zt(i) = zt(i) + (zrfln(i)-zrfl(i))
     .             * (RG/(paprs(i,k)-paprs(i,k+1)))*dtime
     .             * RLVTT/RCPD/(1.0+RVTMP2*zq(i))
         zrfl(i) = zrfln(i)
      ENDIF
      ENDDO
      ENDIF
c
c Calculer Qs et L/Cp*dQs/dT:
c
      IF (thermcep) THEN
         DO i = 1, klon
           zdelta = MAX(0.,SIGN(1.,RTT-zt(i)))
           zcvm5 = R5LES*RLVTT*(1.-zdelta) + R5IES*RLSTT*zdelta
           zcvm5 = zcvm5 /RCPD/(1.0+RVTMP2*zq(i))
           zqs(i) = R2ES*FOEEW(zt(i),zdelta)/pplay(i,k)
           zqs(i) = MIN(0.5,zqs(i))
           zcor = 1./(1.-RETV*zqs(i))
           zqs(i) = zqs(i)*zcor
           zdqs(i) = FOEDE(zt(i),zdelta,zcvm5,zqs(i),zcor)
         ENDDO
      ELSE
         DO i = 1, klon
            IF (zt(i).LT.t_coup) THEN
               zqs(i) = qsats(zt(i))/pplay(i,k)
               zdqs(i) = dqsats(zt(i),zqs(i))
            ELSE
               zqs(i) = qsatl(zt(i))/pplay(i,k)
               zdqs(i) = dqsatl(zt(i),zqs(i))
            ENDIF
         ENDDO
      ENDIF
c
c Determiner la condensation partielle et calculer la quantite
c de l'eau condensee:
c
      IF (cpartiel) THEN

c        print*,'Dans partiel k=',k
c
c   Calcul de l'eau condensee et de la fraction nuageuse et de l'eau
c   nuageuse a partir des PDF de Sandrine Bony.
c   rneb  : fraction nuageuse
c   zqn   : eau totale dans le nuage
c   zcond : eau condensee moyenne dans la maille.
c           on prend en compte le réchauffement qui diminue la partie condensee
c
c   Version avec les raqts

         if (iflag_pdf.eq.0) then

           do i=1,klon
            zdelq = min(ratqs(i,k),0.99) * zq(i)
            rneb(i,k) = (zq(i)+zdelq-zqs(i)) / (2.0*zdelq)
            zqn(i) = (zq(i)+zdelq+zqs(i))/2.0
           enddo

         else
c
c   Version avec les nouvelles PDFs.
           do i=1,klon
              if(zq(i).lt.1.e-15) then
CC Lionel GUEZ                print*,'ZQ(',i,',',k,')=',zq(i)
                zq(i)=1.e-15
              endif
           enddo
           do i=1,klon
            zpdf_sig(i)=ratqs(i,k)*zq(i)
            zpdf_k(i)=-sqrt(log(1.+(zpdf_sig(i)/zq(i))**2))
            zpdf_delta(i)=log(zq(i)/zqs(i))
            zpdf_a(i)=zpdf_delta(i)/(zpdf_k(i)*sqrt(2.))
            zpdf_b(i)=zpdf_k(i)/(2.*sqrt(2.))
            zpdf_e1(i)=zpdf_a(i)-zpdf_b(i)
            zpdf_e1(i)=sign(min(abs(zpdf_e1(i)),5.),zpdf_e1(i))
            zpdf_e1(i)=1.-erf(zpdf_e1(i))
            zpdf_e2(i)=zpdf_a(i)+zpdf_b(i)
            zpdf_e2(i)=sign(min(abs(zpdf_e2(i)),5.),zpdf_e2(i))
            zpdf_e2(i)=1.-erf(zpdf_e2(i))
            if (zpdf_e1(i).lt.1.e-10) then
               rneb(i,k)=0.
               zqn(i)=zqs(i)
            else
               rneb(i,k)=0.5*zpdf_e1(i)
               zqn(i)=zq(i)*zpdf_e2(i)/zpdf_e1(i)
            endif
            
           enddo

        endif ! iflag_pdf

         do i=1,klon
            IF (rneb(i,k) .LE. 0.0) zqn(i) = 0.0
            IF (rneb(i,k) .GE. 1.0) zqn(i) = zq(i)
            rneb(i,k) = MAX(0.0,MIN(1.0,rneb(i,k)))
c           zcond(i) = MAX(0.0,zqn(i)-zqs(i))*rneb(i,k)/(1.+zdqs(i))
c  On ne divise pas par 1+zdqs pour forcer a avoir l'eau predite par
c  la convection.
c  ATTENTION !!! Il va falloir verifier tout ca.
            zcond(i) = MAX(0.0,zqn(i)-zqs(i))*rneb(i,k)
c           print*,'ZDQS ',zdqs(i)
c--Olivier
            rhcl(i,k)=(zqs(i)+zq(i)-zdelq)/2./zqs(i)
            IF (rneb(i,k) .LE. 0.0) rhcl(i,k)=zq(i)/zqs(i)
            IF (rneb(i,k) .GE. 1.0) rhcl(i,k)=1.0
c--fin
           ENDDO
      ELSE
         DO i = 1, klon
            IF (zq(i).GT.zqs(i)) THEN
               rneb(i,k) = 1.0
            ELSE
               rneb(i,k) = 0.0
            ENDIF
            zcond(i) = MAX(0.0,zq(i)-zqs(i))/(1.+zdqs(i))
         ENDDO
      ENDIF
c
      DO i = 1, klon
         zq(i) = zq(i) - zcond(i)
c         zt(i) = zt(i) + zcond(i) * RLVTT/RCPD
         zt(i) = zt(i) + zcond(i) * RLVTT/RCPD/(1.0+RVTMP2*zq(i))
      ENDDO
c
c Partager l'eau condensee en precipitation et eau liquide nuageuse
c
      DO i = 1, klon
      IF (rneb(i,k).GT.0.0) THEN
         zoliq(i) = zcond(i)
         zrho(i) = pplay(i,k) / zt(i) / RD
         zdz(i) = (paprs(i,k)-paprs(i,k+1)) / (zrho(i)*RG)
         zfice(i) = 1.0 - (zt(i)-ztglace) / (273.13-ztglace)
         zfice(i) = MIN(MAX(zfice(i),0.0),1.0)
         zfice(i) = zfice(i)**nexpo
         zneb(i) = MAX(rneb(i,k), seuil_neb)
         radliq(i,k) = zoliq(i)/FLOAT(ninter+1)
      ENDIF
      ENDDO
c
      DO n = 1, ninter
      DO i = 1, klon
      IF (rneb(i,k).GT.0.0) THEN
         zrhol(i) = zrho(i) * zoliq(i) / zneb(i)

         if (ptconv(i,k)) then
            zcl(i)=cld_lc_con
            zct(i)=1./cld_tau_con
         else
            zcl(i)=cld_lc_lsc
            zct(i)=1./cld_tau_lsc
         endif
c  quantité d'eau à élminier.
         zchau(i) = zct(i)*dtime/FLOAT(ninter) * zoliq(i)
     .         *(1.0-EXP(-(zoliq(i)/zneb(i)/zcl(i))**2)) *(1.-zfice(i))
c  meme chose pour la glace.
         if (ptconv(i,k)) then
            zfroi(i) = dtime/FLOAT(ninter)/zdz(i)*zoliq(i)
     .              *fallvc(zrhol(i)) * zfice(i)
         else
            zfroi(i) = dtime/FLOAT(ninter)/zdz(i)*zoliq(i)
     .              *fallvs(zrhol(i)) * zfice(i)
         endif
         ztot(i) = zchau(i) + zfroi(i)
         IF (zneb(i).EQ.seuil_neb) ztot(i) = 0.0
         ztot(i) = MIN(MAX(ztot(i),0.0),zoliq(i))
         zoliq(i) = MAX(zoliq(i)-ztot(i), 0.0)
         radliq(i,k) = radliq(i,k) + zoliq(i)/FLOAT(ninter+1)
      ENDIF
      ENDDO
      ENDDO
c
      DO i = 1, klon
      IF (rneb(i,k).GT.0.0) THEN
         d_ql(i,k) = zoliq(i)
         zrfl(i) = zrfl(i)+ MAX(zcond(i)-zoliq(i),0.0)
     .                    * (paprs(i,k)-paprs(i,k+1))/(RG*dtime)
      ENDIF
      IF (zt(i).LT.RTT) THEN
        psfl(i,k)=zrfl(i)
      ELSE
        prfl(i,k)=zrfl(i)
      ENDIF
      ENDDO
c
c Calculer les tendances de q et de t:
c
      DO i = 1, klon
         d_q(i,k) = zq(i) - q(i,k)
         d_t(i,k) = zt(i) - t(i,k)
      ENDDO
c
cAA--------------- Calcul du lessivage stratiforme  -------------

      DO i = 1,klon
c
         zprec_cond(i) = MAX(zcond(i)-zoliq(i),0.0)
     .                * (paprs(i,k)-paprs(i,k+1))/RG
         IF (rneb(i,k).GT.0.0.and.zprec_cond(i).gt.0.) THEN
cAA lessivage nucleation LMD5 dans la couche elle-meme
            if (t(i,k) .GE. ztglace) THEN
               zalpha_tr = a_tr_sca(3)
            else
               zalpha_tr = a_tr_sca(4)
            endif
            zfrac_lessi = 1. - EXP(zalpha_tr*zprec_cond(i)/zneb(i))
            pfrac_nucl(i,k)=pfrac_nucl(i,k)*(1.-zneb(i)*zfrac_lessi)
            frac_nucl(i,k)= 1.-zneb(i)*zfrac_lessi 
c
c nucleation avec un facteur -1 au lieu de -0.5
            zfrac_lessi = 1. - EXP(-zprec_cond(i)/zneb(i))
            pfrac_1nucl(i,k)=pfrac_1nucl(i,k)*(1.-zneb(i)*zfrac_lessi)
         ENDIF
c
      ENDDO      ! boucle sur i
c
cAA Lessivage par impaction dans les couches en-dessous
      DO kk = k-1, 1, -1
        DO i = 1, klon
          IF (rneb(i,k).GT.0.0.and.zprec_cond(i).gt.0.) THEN
            if (t(i,kk) .GE. ztglace) THEN
              zalpha_tr = a_tr_sca(1)
            else
              zalpha_tr = a_tr_sca(2)
            endif
            zfrac_lessi = 1. - EXP(zalpha_tr*zprec_cond(i)/zneb(i))
            pfrac_impa(i,kk)=pfrac_impa(i,kk)*(1.-zneb(i)*zfrac_lessi)
            frac_impa(i,kk)= 1.-zneb(i)*zfrac_lessi
          ENDIF
        ENDDO
      ENDDO
c
cAA----------------------------------------------------------
c                     FIN DE BOUCLE SUR K   
 9999 CONTINUE
c
cAA-----------------------------------------------------------
c
c Pluie ou neige au sol selon la temperature de la 1ere couche
c
      DO i = 1, klon
      IF ((t(i,1)+d_t(i,1)) .LT. RTT) THEN
         snow(i) = zrfl(i)
         zlh_solid(i) = RLSTT-RLVTT
      ELSE
         rain(i) = zrfl(i)
         zlh_solid(i) = 0.
      ENDIF
      ENDDO
C
C For energy conservation : when snow is present, the solification
c latent heat is considered.
      DO k = 1, klev
        DO i = 1, klon
          zcpair=RCPD*(1.0+RVTMP2*(q(i,k)+d_q(i,k)))
          zmair=(paprs(i,k)-paprs(i,k+1))/RG
          zm_solid = (prfl(i,k)-prfl(i,k+1)+psfl(i,k)-psfl(i,k+1))*dtime
          d_t(i,k) = d_t(i,k) + zlh_solid(i) *zm_solid / (zcpair*zmair)
        END DO 
      END DO
c
      RETURN
      END
