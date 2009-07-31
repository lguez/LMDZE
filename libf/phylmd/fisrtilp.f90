SUBROUTINE fisrtilp(dtime,paprs,pplay,t,q,ptconv,ratqs,d_t,d_q,d_ql,rneb, &
     radliq,rain,snow,pfrac_impa,pfrac_nucl,pfrac_1nucl,frac_impa, &
     frac_nucl,prfl,psfl,rhcl)

  ! From phylmd/fisrtilp.F,v 1.2 2004/11/09 16:55:40
  ! Auteur(s): Z.X. Li (LMD/CNRS)
  ! Date: le 20 mars 1995
  ! Objet: condensation et precipitation stratiforme.
  !        schema de nuage

  USE dimens_m
  USE dimphy
  USE tracstoke
  USE yomcst
  USE yoethf
  USE fcttre
  USE comfisrtilp
  use numer_rec, only: erf

  IMPLICIT NONE

  ! Arguments:

  REAL, INTENT (IN) :: & ! intervalle du temps (s)                
       dtime
  REAL, INTENT (IN) :: paprs(klon,klev+1) ! pression a inter-couche   
  REAL, INTENT (IN) :: pplay(klon,klev) ! pression au milieu de couche
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
  ! Coeffients de fraction lessivee : pour OFF-LINE

  REAL pfrac_nucl(klon,klev)
  REAL pfrac_1nucl(klon,klev)
  REAL pfrac_impa(klon,klev)

  ! Fraction d'aerosols lessivee par impaction et par nucleation
  ! POur ON-LINE

  REAL frac_impa(klon,klev)
  REAL frac_nucl(klon,klev)
  REAL zct(klon), zcl(klon)
  !AA

  ! Options du programme:

  REAL seuil_neb ! un nuage existe vraiment au-dela
  PARAMETER (seuil_neb=0.001)

  INTEGER ninter ! sous-intervals pour la precipitation
  PARAMETER (ninter=5)
  LOGICAL evap_prec ! evaporation de la pluie                       
  PARAMETER (evap_prec=.TRUE.)
  REAL ratqs(klon,klev) ! determine la largeur de distribution de vapeur
  LOGICAL ptconv(klon,klev) ! determine la largeur de distribution de vapeur
  REAL zpdf_sig(klon), zpdf_k(klon), zpdf_delta(klon)
  REAL zpdf_a(klon), zpdf_b(klon), zpdf_e1(klon), zpdf_e2(klon)

  LOGICAL cpartiel ! condensation partielle                         
  PARAMETER (cpartiel=.TRUE.)
  REAL t_coup
  PARAMETER (t_coup=234.0)

  ! Variables locales:

  INTEGER i, k, n, kk
  REAL zqs(klon), zdqs(klon), zdelta, zcor, zcvm5
  REAL zrfl(klon), zrfln(klon), zqev, zqevt
  REAL zoliq(klon), zcond(klon), zq(klon), zqn(klon), zdelq
  REAL ztglace, zt(klon)
  INTEGER nexpo ! exponentiel pour glace/eau                        
  REAL zdz(klon), zrho(klon), ztot(klon), zrhol(klon)
  REAL zchau(klon), zfroi(klon), zfice(klon), zneb(klon)

  LOGICAL appel1er
  SAVE appel1er

  !---------------------------------------------------------------

  !AA Variables traceurs:
  !AA  Provisoire !!! Parametres alpha du lessivage
  !AA  A priori on a 4 scavenging numbers possibles

  REAL a_tr_sca(4)
  SAVE a_tr_sca

  ! Variables intermediaires

  REAL zalpha_tr
  REAL zfrac_lessi
  REAL zprec_cond(klon)
  !AA
  REAL zmair, zcpair, zcpeau
  !     Pour la conversion eau-neige
  REAL zlh_solid(klon), zm_solid
  !IM
  INTEGER klevm1
  !---------------------------------------------------------------

  ! Fonctions en ligne:

  REAL fallvs, fallvc ! vitesse de chute pour crystaux de glace      
  REAL zzz

  fallvc(zzz) = 3.29/2.0*((zzz)**0.16)*ffallv_con
  fallvs(zzz) = 3.29/2.0*((zzz)**0.16)*ffallv_lsc

  DATA appel1er/ .TRUE./
  !ym
  zdelq = 0.0

  IF (appel1er) THEN

     PRINT *, 'fisrtilp, ninter:', ninter
     PRINT *, 'fisrtilp, evap_prec:', evap_prec
     PRINT *, 'fisrtilp, cpartiel:', cpartiel
     IF (abs(dtime/float(ninter)-360.0)>0.001) THEN
        PRINT *, 'fisrtilp: Ce n est pas prevu, voir Z.X.Li', dtime
        PRINT *, 'Je prefere un sous-intervalle de 6 minutes'
        !         stop 1
     END IF
     appel1er = .FALSE.

     !AA initialiation provisoire
     a_tr_sca(1) = -0.5
     a_tr_sca(2) = -0.5
     a_tr_sca(3) = -0.5
     a_tr_sca(4) = -0.5

     !AA Initialisation a 1 des coefs des fractions lessivees

     DO k = 1, klev
        DO i = 1, klon
           pfrac_nucl(i,k) = 1.
           pfrac_1nucl(i,k) = 1.
           pfrac_impa(i,k) = 1.
        END DO
     END DO


  END IF !  test sur appel1er
  !MAf Initialisation a 0 de zoliq
  DO i = 1, klon
     zoliq(i) = 0.
  END DO
  ! Determiner les nuages froids par leur temperature
  !  nexpo regle la raideur de la transition eau liquide / eau glace.

  ztglace = rtt - 15.0
  nexpo = 6
  !cc      nexpo = 1

  ! Initialiser les sorties:

  DO k = 1, klev + 1
     DO i = 1, klon
        prfl(i,k) = 0.0
        psfl(i,k) = 0.0
     END DO
  END DO

  DO k = 1, klev
     DO i = 1, klon
        d_t(i,k) = 0.0
        d_q(i,k) = 0.0
        d_ql(i,k) = 0.0
        rneb(i,k) = 0.0
        radliq(i,k) = 0.0
        frac_nucl(i,k) = 1.
        frac_impa(i,k) = 1.
     END DO
  END DO
  DO i = 1, klon
     rain(i) = 0.0
     snow(i) = 0.0
  END DO

  ! Initialiser le flux de precipitation a zero

  DO i = 1, klon
     zrfl(i) = 0.0
     zneb(i) = seuil_neb
  END DO


  !AA Pour plus de securite

  zalpha_tr = 0.
  zfrac_lessi = 0.

  !AA----------------------------------------------------------

  ! Boucle verticale (du haut vers le bas)

  !IM : klevm1
  klevm1 = klev - 1
  DO  k = klev, 1, -1

     !AA----------------------------------------------------------

     DO i = 1, klon
        zt(i) = t(i,k)
        zq(i) = q(i,k)
     END DO

     ! Calculer la varition de temp. de l'air du a la chaleur sensible
     ! transporter par la pluie.
     ! Il resterait a rajouter cet effet de la chaleur sensible sur les
     ! flux de surface, du a la diff. de temp. entre le 1er niveau et la
     ! surface.

     DO i = 1, klon
        IF (k<=klevm1) THEN
           zmair = (paprs(i,k)-paprs(i,k+1))/rg
           zcpair = rcpd*(1.0+rvtmp2*zq(i))
           zcpeau = rcpd*rvtmp2
           zt(i) = ((t(i,k+1)+d_t(i,k+1))*zrfl(i)*dtime*zcpeau+zmair*zcpair* &
                zt(i))/(zmair*zcpair+zrfl(i)*dtime*zcpeau)
           !C        WRITE (6,*) 'cppluie ', zt(i)-(t(i,k+1)+d_t(i,k+1))
        END IF
     END DO

     ! Calculer l'evaporation de la precipitation



     IF (evap_prec) THEN
        DO i = 1, klon
           IF (zrfl(i)>0.) THEN
              IF (thermcep) THEN
                 zdelta = max(0.,sign(1.,rtt-zt(i)))
                 zqs(i) = r2es*foeew(zt(i),zdelta)/pplay(i,k)
                 zqs(i) = min(0.5,zqs(i))
                 zcor = 1./(1.-retv*zqs(i))
                 zqs(i) = zqs(i)*zcor
              ELSE
                 IF (zt(i)<t_coup) THEN
                    zqs(i) = qsats(zt(i))/pplay(i,k)
                 ELSE
                    zqs(i) = qsatl(zt(i))/pplay(i,k)
                 END IF
              END IF
              zqev = max(0.0,(zqs(i)-zq(i))*zneb(i))
              zqevt = coef_eva*(1.0-zq(i)/zqs(i))*sqrt(zrfl(i))* &
                   (paprs(i,k)-paprs(i,k+1))/pplay(i,k)*zt(i)*rd/rg
              zqevt = max(0.0,min(zqevt,zrfl(i)))*rg*dtime/ &
                   (paprs(i,k)-paprs(i,k+1))
              zqev = min(zqev,zqevt)
              zrfln(i) = zrfl(i) - zqev*(paprs(i,k)-paprs(i,k+1))/rg/dtime

              ! pour la glace, on réévapore toute la précip dans la couche du dessous
              ! la glace venant de la couche du dessus est simplement dans la couche
              ! du dessous.

              IF (zt(i)<t_coup .AND. reevap_ice) zrfln(i) = 0.

              zq(i) = zq(i) - (zrfln(i)-zrfl(i))*(rg/(paprs(i,k)-paprs(i, &
                   k+1)))*dtime
              zt(i) = zt(i) + (zrfln(i)-zrfl(i))*(rg/(paprs(i,k)-paprs(i, &
                   k+1)))*dtime*rlvtt/rcpd/(1.0+rvtmp2*zq(i))
              zrfl(i) = zrfln(i)
           END IF
        END DO
     END IF

     ! Calculer Qs et L/Cp*dQs/dT:

     IF (thermcep) THEN
        DO i = 1, klon
           zdelta = max(0.,sign(1.,rtt-zt(i)))
           zcvm5 = r5les*rlvtt*(1.-zdelta) + r5ies*rlstt*zdelta
           zcvm5 = zcvm5/rcpd/(1.0+rvtmp2*zq(i))
           zqs(i) = r2es*foeew(zt(i),zdelta)/pplay(i,k)
           zqs(i) = min(0.5,zqs(i))
           zcor = 1./(1.-retv*zqs(i))
           zqs(i) = zqs(i)*zcor
           zdqs(i) = foede(zt(i),zdelta,zcvm5,zqs(i),zcor)
        END DO
     ELSE
        DO i = 1, klon
           IF (zt(i)<t_coup) THEN
              zqs(i) = qsats(zt(i))/pplay(i,k)
              zdqs(i) = dqsats(zt(i),zqs(i))
           ELSE
              zqs(i) = qsatl(zt(i))/pplay(i,k)
              zdqs(i) = dqsatl(zt(i),zqs(i))
           END IF
        END DO
     END IF

     ! Determiner la condensation partielle et calculer la quantite
     ! de l'eau condensee:

     IF (cpartiel) THEN

        !        print*,'Dans partiel k=',k

        !   Calcul de l'eau condensee et de la fraction nuageuse et de l'eau
        !   nuageuse a partir des PDF de Sandrine Bony.
        !   rneb  : fraction nuageuse
        !   zqn   : eau totale dans le nuage
        !   zcond : eau condensee moyenne dans la maille.
        !           on prend en compte le réchauffement qui diminue la partie condensee

        !   Version avec les raqts

        IF (iflag_pdf==0) THEN

           DO i = 1, klon
              zdelq = min(ratqs(i,k),0.99)*zq(i)
              rneb(i,k) = (zq(i)+zdelq-zqs(i))/(2.0*zdelq)
              zqn(i) = (zq(i)+zdelq+zqs(i))/2.0
           END DO

        ELSE

           !   Version avec les nouvelles PDFs.
           DO i = 1, klon
              IF (zq(i)<1.E-15) THEN
                 !C Lionel GUEZ                print*,'ZQ(',i,',',k,')=',zq(i)
                 zq(i) = 1.E-15
              END IF
           END DO
           DO i = 1, klon
              zpdf_sig(i) = ratqs(i,k)*zq(i)
              zpdf_k(i) = -sqrt(log(1.+(zpdf_sig(i)/zq(i))**2))
              zpdf_delta(i) = log(zq(i)/zqs(i))
              zpdf_a(i) = zpdf_delta(i)/(zpdf_k(i)*sqrt(2.))
              zpdf_b(i) = zpdf_k(i)/(2.*sqrt(2.))
              zpdf_e1(i) = zpdf_a(i) - zpdf_b(i)
              zpdf_e1(i) = sign(min(abs(zpdf_e1(i)),5.),zpdf_e1(i))
              zpdf_e1(i) = 1. - erf(zpdf_e1(i))
              zpdf_e2(i) = zpdf_a(i) + zpdf_b(i)
              zpdf_e2(i) = sign(min(abs(zpdf_e2(i)),5.),zpdf_e2(i))
              zpdf_e2(i) = 1. - erf(zpdf_e2(i))
              IF (zpdf_e1(i)<1.E-10) THEN
                 rneb(i,k) = 0.
                 zqn(i) = zqs(i)
              ELSE
                 rneb(i,k) = 0.5*zpdf_e1(i)
                 zqn(i) = zq(i)*zpdf_e2(i)/zpdf_e1(i)
              END IF

           END DO


        END IF
        ! iflag_pdf                                               
        DO i = 1, klon
           IF (rneb(i,k)<=0.0) zqn(i) = 0.0
           IF (rneb(i,k)>=1.0) zqn(i) = zq(i)
           rneb(i,k) = max(0.0,min(1.0,rneb(i,k)))
           !           zcond(i) = MAX(0.0,zqn(i)-zqs(i))*rneb(i,k)/(1.+zdqs(i))
           !  On ne divise pas par 1+zdqs pour forcer a avoir l'eau predite par
           !  la convection.
           !  ATTENTION !!! Il va falloir verifier tout ca.
           zcond(i) = max(0.0,zqn(i)-zqs(i))*rneb(i,k)
           !           print*,'ZDQS ',zdqs(i)
           !--Olivier
           rhcl(i,k) = (zqs(i)+zq(i)-zdelq)/2./zqs(i)
           IF (rneb(i,k)<=0.0) rhcl(i,k) = zq(i)/zqs(i)
           IF (rneb(i,k)>=1.0) rhcl(i,k) = 1.0
           !--fin
        END DO
     ELSE
        DO i = 1, klon
           IF (zq(i)>zqs(i)) THEN
              rneb(i,k) = 1.0
           ELSE
              rneb(i,k) = 0.0
           END IF
           zcond(i) = max(0.0,zq(i)-zqs(i))/(1.+zdqs(i))
        END DO
     END IF

     DO i = 1, klon
        zq(i) = zq(i) - zcond(i)
        !         zt(i) = zt(i) + zcond(i) * RLVTT/RCPD
        zt(i) = zt(i) + zcond(i)*rlvtt/rcpd/(1.0+rvtmp2*zq(i))
     END DO

     ! Partager l'eau condensee en precipitation et eau liquide nuageuse

     DO i = 1, klon
        IF (rneb(i,k)>0.0) THEN
           zoliq(i) = zcond(i)
           zrho(i) = pplay(i,k)/zt(i)/rd
           zdz(i) = (paprs(i,k)-paprs(i,k+1))/(zrho(i)*rg)
           zfice(i) = 1.0 - (zt(i)-ztglace)/(273.13-ztglace)
           zfice(i) = min(max(zfice(i),0.0),1.0)
           zfice(i) = zfice(i)**nexpo
           zneb(i) = max(rneb(i,k),seuil_neb)
           radliq(i,k) = zoliq(i)/float(ninter+1)
        END IF
     END DO

     DO n = 1, ninter
        DO i = 1, klon
           IF (rneb(i,k)>0.0) THEN
              zrhol(i) = zrho(i)*zoliq(i)/zneb(i)

              IF (ptconv(i,k)) THEN
                 zcl(i) = cld_lc_con
                 zct(i) = 1./cld_tau_con
              ELSE
                 zcl(i) = cld_lc_lsc
                 zct(i) = 1./cld_tau_lsc
              END IF
              !  quantité d'eau à élminier.
              zchau(i) = zct(i)*dtime/float(ninter)*zoliq(i)* &
                   (1.0-exp(-(zoliq(i)/zneb(i)/zcl(i))**2))*(1.-zfice(i))
              !  meme chose pour la glace.
              IF (ptconv(i,k)) THEN
                 zfroi(i) = dtime/float(ninter)/zdz(i)*zoliq(i)* &
                      fallvc(zrhol(i))*zfice(i)
              ELSE
                 zfroi(i) = dtime/float(ninter)/zdz(i)*zoliq(i)* &
                      fallvs(zrhol(i))*zfice(i)
              END IF
              ztot(i) = zchau(i) + zfroi(i)
              IF (zneb(i)==seuil_neb) ztot(i) = 0.0
              ztot(i) = min(max(ztot(i),0.0),zoliq(i))
              zoliq(i) = max(zoliq(i)-ztot(i),0.0)
              radliq(i,k) = radliq(i,k) + zoliq(i)/float(ninter+1)
           END IF
        END DO
     END DO

     DO i = 1, klon
        IF (rneb(i,k)>0.0) THEN
           d_ql(i,k) = zoliq(i)
           zrfl(i) = zrfl(i) + max(zcond(i)-zoliq(i),0.0)*(paprs(i,k)-paprs(i &
                ,k+1))/(rg*dtime)
        END IF
        IF (zt(i)<rtt) THEN
           psfl(i,k) = zrfl(i)
        ELSE
           prfl(i,k) = zrfl(i)
        END IF
     END DO

     ! Calculer les tendances de q et de t:

     DO i = 1, klon
        d_q(i,k) = zq(i) - q(i,k)
        d_t(i,k) = zt(i) - t(i,k)
     END DO

     !AA--------------- Calcul du lessivage stratiforme  -------------

     DO i = 1, klon
        zprec_cond(i) = max(zcond(i)-zoliq(i),0.0)* &
             (paprs(i,k)-paprs(i,k+1))/rg
        IF (rneb(i,k)>0.0 .AND. zprec_cond(i)>0.) THEN
           !AA lessivage nucleation LMD5 dans la couche elle-meme
           IF (t(i,k)>=ztglace) THEN
              zalpha_tr = a_tr_sca(3)
           ELSE
              zalpha_tr = a_tr_sca(4)
           END IF
           zfrac_lessi = 1. - exp(zalpha_tr*zprec_cond(i)/zneb(i))
           pfrac_nucl(i,k) = pfrac_nucl(i,k)*(1.-zneb(i)*zfrac_lessi)
           frac_nucl(i,k) = 1. - zneb(i)*zfrac_lessi

           ! nucleation avec un facteur -1 au lieu de -0.5
           zfrac_lessi = 1. - exp(-zprec_cond(i)/zneb(i))
           pfrac_1nucl(i,k) = pfrac_1nucl(i,k)*(1.-zneb(i)*zfrac_lessi)
        END IF


     END DO
     !AA Lessivage par impaction dans les couches en-dessous
     ! boucle sur i                                         
     DO kk = k - 1, 1, -1
        DO i = 1, klon
           IF (rneb(i,k)>0.0 .AND. zprec_cond(i)>0.) THEN
              IF (t(i,kk)>=ztglace) THEN
                 zalpha_tr = a_tr_sca(1)
              ELSE
                 zalpha_tr = a_tr_sca(2)
              END IF
              zfrac_lessi = 1. - exp(zalpha_tr*zprec_cond(i)/zneb(i))
              pfrac_impa(i,kk) = pfrac_impa(i,kk)*(1.-zneb(i)*zfrac_lessi)
              frac_impa(i,kk) = 1. - zneb(i)*zfrac_lessi
           END IF
        END DO
     END DO

     !AA----------------------------------------------------------
     !                     FIN DE BOUCLE SUR K
  end DO

  !AA-----------------------------------------------------------

  ! Pluie ou neige au sol selon la temperature de la 1ere couche

  DO i = 1, klon
     IF ((t(i,1)+d_t(i,1))<rtt) THEN
        snow(i) = zrfl(i)
        zlh_solid(i) = rlstt - rlvtt
     ELSE
        rain(i) = zrfl(i)
        zlh_solid(i) = 0.
     END IF
  END DO

  ! For energy conservation : when snow is present, the solification
  ! latent heat is considered.
  DO k = 1, klev
     DO i = 1, klon
        zcpair = rcpd*(1.0+rvtmp2*(q(i,k)+d_q(i,k)))
        zmair = (paprs(i,k)-paprs(i,k+1))/rg
        zm_solid = (prfl(i,k)-prfl(i,k+1)+psfl(i,k)-psfl(i,k+1))*dtime
        d_t(i,k) = d_t(i,k) + zlh_solid(i)*zm_solid/(zcpair*zmair)
     END DO
  END DO

END SUBROUTINE fisrtilp
