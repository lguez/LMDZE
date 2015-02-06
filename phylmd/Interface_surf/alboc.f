module alboc_m

  IMPLICIT NONE

contains

  SUBROUTINE alboc(rjour, rlat, albedo)

    ! From LMDZ4/libf/phylmd/albedo.F, version 1.2 2005/02/07 15:00:52

    ! Auteur(s): Z.X. Li (LMD/CNRS) (adaptation du GCM du LMD)
    ! Date: le 16 mars 1995
    ! Objet: Calculer l'albedo sur l'ocean
    ! Methode: Integrer numeriquement l'albedo pendant une journee

    USE dimphy, only: klon
    USE yomcst, only: r_incl
    USE orbite_m, ONLY: orbite

    ! Arguments;
    ! rjour (in, R) : jour dans l'annee (a compter du 1 janvier)
    ! rlat (in, R) : latitude en degre
    ! albedo (out, R): albedo obtenu (de 0 a 1)

    REAL fmagic ! un facteur magique pour regler l'albedo
    ! cc PARAMETER (fmagic=0.7)
    ! ccIM => a remplacer
    ! PARAMETER (fmagic=1.32)
    PARAMETER (fmagic=1.0)
    ! PARAMETER (fmagic=0.7)
    INTEGER npts ! il controle la precision de l'integration
    PARAMETER (npts=120) ! 120 correspond a l'interval 6 minutes

    REAL rlat(klon), rjour, albedo(klon)
    REAL zdist, zlonsun, zpi, zdeclin
    REAL rmu, alb, srmu, salb, fauxo, aa, bb
    INTEGER i, k
    ! ccIM
    LOGICAL ancien_albedo
    PARAMETER (ancien_albedo=.FALSE.)
    ! SAVE albedo

    IF (ancien_albedo) THEN

       zpi = 4.*atan(1.)

       ! Calculer la longitude vraie de l'orbite terrestre:
       CALL orbite(rjour, zlonsun, zdist)

       ! Calculer la declinaison du soleil (qui varie entre + et - R_incl):
       zdeclin = asin(sin(zlonsun*zpi/180.0)*sin(r_incl*zpi/180.0))

       DO i = 1, klon
          aa = sin(rlat(i)*zpi/180.0)*sin(zdeclin)
          bb = cos(rlat(i)*zpi/180.0)*cos(zdeclin)

          ! Midi local (angle du temps = 0.0):
          rmu = aa + bb*cos(0.0)
          rmu = max(0.0, rmu)
          fauxo = (1.47-acos(rmu))/.15
          alb = 0.03 + 0.630/(1.+fauxo*fauxo)
          srmu = rmu
          salb = alb*rmu

          ! Faire l'integration numerique de midi a minuit (le facteur 2
          ! prend en compte l'autre moitie de la journee):
          DO k = 1, npts
             rmu = aa + bb*cos(float(k)/float(npts)*zpi)
             rmu = max(0.0, rmu)
             fauxo = (1.47-acos(rmu))/.15
             alb = 0.03 + 0.630/(1.+fauxo*fauxo)
             srmu = srmu + rmu*2.0
             salb = salb + alb*rmu*2.0
          END DO
          IF (srmu/=0.0) THEN
             albedo(i) = salb/srmu*fmagic
          ELSE ! nuit polaire (on peut prendre une valeur quelconque)
             albedo(i) = fmagic
          END IF
       END DO

       ! nouvel albedo

    ELSE

       zpi = 4.*atan(1.)

       ! Calculer la longitude vraie de l'orbite terrestre:
       CALL orbite(rjour, zlonsun, zdist)

       ! Calculer la declinaison du soleil (qui varie entre + et - R_incl):
       zdeclin = asin(sin(zlonsun*zpi/180.0)*sin(r_incl*zpi/180.0))

       DO i = 1, klon
          aa = sin(rlat(i)*zpi/180.0)*sin(zdeclin)
          bb = cos(rlat(i)*zpi/180.0)*cos(zdeclin)

          ! Midi local (angle du temps = 0.0):
          rmu = aa + bb*cos(0.0)
          rmu = max(0.0, rmu)
          ! IM cf. PB alb = 0.058/(rmu + 0.30)
          ! alb = 0.058/(rmu + 0.30) * 1.5
          alb = 0.058/(rmu+0.30)*1.2
          ! alb = 0.058/(rmu + 0.30) * 1.3
          srmu = rmu
          salb = alb*rmu

          ! Faire l'integration numerique de midi a minuit (le facteur 2
          ! prend en compte l'autre moitie de la journee):
          DO k = 1, npts
             rmu = aa + bb*cos(float(k)/float(npts)*zpi)
             rmu = max(0.0, rmu)
             ! IM cf. PB alb = 0.058/(rmu + 0.30)
             ! alb = 0.058/(rmu + 0.30) * 1.5
             alb = 0.058/(rmu+0.30)*1.2
             ! alb = 0.058/(rmu + 0.30) * 1.3
             srmu = srmu + rmu*2.0
             salb = salb + alb*rmu*2.0
          END DO
          IF (srmu/=0.0) THEN
             albedo(i) = salb/srmu*fmagic
          ELSE ! nuit polaire (on peut prendre une valeur quelconque)
             albedo(i) = fmagic
          END IF
       END DO
    END IF

  END SUBROUTINE alboc

end module alboc_m
