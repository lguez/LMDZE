!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/albedo.F,v 1.2 2005/02/07 15:00:52 fairhead Exp $
!
c
c
      SUBROUTINE alboc(rjour,rlat,albedo)
      use dimens_m
      use dimphy
      use YOMCST
      use orbite_m, only: orbite
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) (adaptation du GCM du LMD)
c Date: le 16 mars 1995
c Objet: Calculer l'albedo sur l'ocean
c Methode: Integrer numeriquement l'albedo pendant une journee
c
c Arguments;
c rjour (in,R)  : jour dans l'annee (a compter du 1 janvier)
c rlat (in,R)   : latitude en degre
c albedo (out,R): albedo obtenu (de 0 a 1)
c======================================================================
c
      REAL fmagic ! un facteur magique pour regler l'albedo
ccc      PARAMETER (fmagic=0.7)
cccIM => a remplacer  
c       PARAMETER (fmagic=1.32)
        PARAMETER (fmagic=1.0)
c       PARAMETER (fmagic=0.7)
      INTEGER npts ! il controle la precision de l'integration
      PARAMETER (npts=120) ! 120 correspond a l'interval 6 minutes
c
      REAL rlat(klon), rjour, albedo(klon)
      REAL zdist, zlonsun, zpi, zdeclin
      REAL rmu,alb, srmu, salb, fauxo, aa, bb
      INTEGER i, k
cccIM
      LOGICAL ancien_albedo
      PARAMETER(ancien_albedo=.FALSE.) 
c     SAVE albedo
c
      IF ( ancien_albedo ) THEN
c
      zpi = 4. * ATAN(1.)
c
c Calculer la longitude vraie de l'orbite terrestre:
      CALL orbite(rjour,zlonsun,zdist)
c
c Calculer la declinaison du soleil (qui varie entre + et - R_incl):
      zdeclin = ASIN(SIN(zlonsun*zpi/180.0)*SIN(R_incl*zpi/180.0))
c
      DO 999 i=1,klon
      aa = SIN(rlat(i)*zpi/180.0) * SIN(zdeclin)
      bb = COS(rlat(i)*zpi/180.0) * COS(zdeclin)
c
c Midi local (angle du temps = 0.0):
      rmu = aa + bb * COS(0.0)
      rmu = MAX(0.0, rmu)
      fauxo = (1.47-ACOS(rmu))/.15
      alb = 0.03+0.630/(1.+fauxo*fauxo)
      srmu = rmu
      salb = alb * rmu
c
c Faire l'integration numerique de midi a minuit (le facteur 2
c prend en compte l'autre moitie de la journee):
      DO k = 1, npts
         rmu = aa + bb * COS(FLOAT(k)/FLOAT(npts)*zpi)
         rmu = MAX(0.0, rmu)
         fauxo = (1.47-ACOS(rmu))/.15
         alb = 0.03+0.630/(1.+fauxo*fauxo)
         srmu = srmu + rmu * 2.0
         salb = salb + alb*rmu * 2.0
      ENDDO
      IF (srmu .NE. 0.0) THEN
         albedo(i) = salb / srmu * fmagic
      ELSE ! nuit polaire (on peut prendre une valeur quelconque)
         albedo(i) = fmagic
      ENDIF
  999 CONTINUE
c
c nouvel albedo 
c
      ELSE
c
      zpi = 4. * ATAN(1.)
c
c Calculer la longitude vraie de l'orbite terrestre:
      CALL orbite(rjour,zlonsun,zdist)
c
c Calculer la declinaison du soleil (qui varie entre + et - R_incl):
      zdeclin = ASIN(SIN(zlonsun*zpi/180.0)*SIN(R_incl*zpi/180.0))
c
      DO 1999 i=1,klon
      aa = SIN(rlat(i)*zpi/180.0) * SIN(zdeclin)
      bb = COS(rlat(i)*zpi/180.0) * COS(zdeclin)
c
c Midi local (angle du temps = 0.0):
      rmu = aa + bb * COS(0.0)
      rmu = MAX(0.0, rmu)
cIM cf. PB  alb = 0.058/(rmu + 0.30)
c     alb = 0.058/(rmu + 0.30) * 1.5
      alb = 0.058/(rmu + 0.30) * 1.2
c     alb = 0.058/(rmu + 0.30) * 1.3
      srmu = rmu
      salb = alb * rmu
c
c Faire l'integration numerique de midi a minuit (le facteur 2
c prend en compte l'autre moitie de la journee):
      DO k = 1, npts
         rmu = aa + bb * COS(FLOAT(k)/FLOAT(npts)*zpi)
         rmu = MAX(0.0, rmu)
cIM cf. PB      alb = 0.058/(rmu + 0.30)
c        alb = 0.058/(rmu + 0.30) * 1.5
         alb = 0.058/(rmu + 0.30) * 1.2
c        alb = 0.058/(rmu + 0.30) * 1.3
         srmu = srmu + rmu * 2.0
         salb = salb + alb*rmu * 2.0
      ENDDO
      IF (srmu .NE. 0.0) THEN
         albedo(i) = salb / srmu * fmagic
      ELSE ! nuit polaire (on peut prendre une valeur quelconque)
         albedo(i) = fmagic
      ENDIF
1999  CONTINUE
      ENDIF
      RETURN
      END
c=====================================================================
      SUBROUTINE alboc_cd(rmu0,albedo)
      use dimens_m
      use dimphy
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS)
c date: 19940624
c Calculer l'albedo sur l'ocean en fonction de l'angle zenithal moyen
c Formule due a Larson and Barkstrom (1977) Proc. of the symposium
C on radiation in the atmosphere, 19-28 August 1976, science Press,
C 1977 pp 451-453, ou These de 3eme cycle de Sylvie Joussaume.
c
c Arguments
c rmu0    (in): cosinus de l'angle solaire zenithal
c albedo (out): albedo de surface de l'ocean
c======================================================================
      REAL rmu0(klon), albedo(klon)
c
      REAL fmagic ! un facteur magique pour regler l'albedo
ccc      PARAMETER (fmagic=0.7)
cccIM => a remplacer  
c       PARAMETER (fmagic=1.32)
        PARAMETER (fmagic=1.0)
c       PARAMETER (fmagic=0.7) 
c
      REAL fauxo
      INTEGER i
cccIM
      LOGICAL ancien_albedo
      PARAMETER(ancien_albedo=.FALSE.) 
c     SAVE albedo
c
      IF ( ancien_albedo ) THEN
c
      DO i = 1, klon
c
         rmu0(i) = MAX(rmu0(i),0.0)
c
         fauxo = ( 1.47 - ACOS( rmu0(i) ) )/0.15
         albedo(i) = fmagic*( .03 + .630/( 1. + fauxo*fauxo))
         albedo(i) = MAX(MIN(albedo(i),0.60),0.04)
      ENDDO
c
c nouvel albedo 
c
      ELSE
c
      DO i = 1, klon
         rmu0(i) = MAX(rmu0(i),0.0)
cIM:orig albedo(i) = 0.058/(rmu0(i) + 0.30)
         albedo(i) = fmagic * 0.058/(rmu0(i) + 0.30)
         albedo(i) = MAX(MIN(albedo(i),0.60),0.04)
      ENDDO
c
      ENDIF
c
      RETURN
      END
c========================================================================
