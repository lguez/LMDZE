module alboc_m

  IMPLICIT NONE

contains

  SUBROUTINE alboc(jour, rlat, albedo)

    ! From LMDZ4/libf/phylmd/albedo.F, version 1.2 2005/02/07 15:00:52

    ! Author: Z. X. Li (LMD/CNRS) (adaptation du GCM du LMD)
    ! Date : 16 mars 1995
    ! Objet : Calculer l'alb\'edo sur l'oc\'ean
    ! M\'ethode: int\'egrer num\'eriquement l'alb\'edo pendant une journ\'ee

    use nr_util, only: pi
    USE orbite_m, ONLY: orbite
    USE yomcst, only: r_incl

    integer, intent(in):: jour ! jour dans l'annee (a compter du 1 janvier)
    REAL, intent(in):: rlat(:) ! latitude en degre
    real, intent(out):: albedo(:) ! albedo obtenu (de 0 a 1)

    ! Local:

    REAL, PARAMETER:: fmagic=1. ! un facteur magique pour regler l'albedo

    INTEGER, PARAMETER:: npts =120
    ! Contr\^ole la pr\'ecision de l'int\'egration. 120 correspond \`a
    ! l'intervalle 6 minutes.

    REAL zdist, zlonsun, zdeclin
    REAL rmu, alb, srmu, salb, aa, bb
    INTEGER i, k

    !----------------------------------------------------------------------

    ! Calculer la longitude vraie de l'orbite terrestre:
    CALL orbite(real(jour), zlonsun, zdist)

    ! Calculer la declinaison du soleil (qui varie entre + et - R_incl):
    zdeclin = asin(sin(zlonsun*pi/180.0)*sin(r_incl*pi/180.0))

    DO i = 1, size(rlat)
       aa = sin(rlat(i)*pi/180.0)*sin(zdeclin)
       bb = cos(rlat(i)*pi/180.0)*cos(zdeclin)

       ! Midi local (angle du temps = 0.0):
       rmu = aa + bb*cos(0.0)
       rmu = max(0.0, rmu)
       alb = 0.058/(rmu+0.30)*1.2
       srmu = rmu
       salb = alb*rmu

       ! Faire l'integration numerique de midi a minuit (le facteur 2
       ! prend en compte l'autre moitie de la journee):
       DO k = 1, npts
          rmu = aa + bb*cos(float(k)/float(npts)*pi)
          rmu = max(0.0, rmu)
          alb = 0.058/(rmu+0.30)*1.2
          srmu = srmu + rmu*2.0
          salb = salb + alb*rmu*2.0
       END DO
       IF (srmu/=0.0) THEN
          albedo(i) = salb/srmu*fmagic
       ELSE ! nuit polaire (on peut prendre une valeur quelconque)
          albedo(i) = fmagic
       END IF
    END DO

  END SUBROUTINE alboc

end module alboc_m
