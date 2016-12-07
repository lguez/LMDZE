module alboc_m

  IMPLICIT NONE

contains

  pure function alboc(jour, rlat)

    ! From LMDZ4/libf/phylmd/albedo.F, version 1.2, 2005/02/07 15:00:52

    ! Author: Z. X. Li (LMD/CNRS)
    ! Date : 16 mars 1995
    ! Objet : Calculer l'alb\'edo sur l'oc\'ean
    ! M\'ethode: int\'egrer num\'eriquement l'alb\'edo pendant une journ\'ee

    use nr_util, only: pi
    USE orbite_m, ONLY: orbite
    USE yomcst, only: r_incl

    integer, intent(in):: jour ! jour dans l'annee (a compter du 1 janvier)
    REAL, intent(in):: rlat(:) ! latitude en degre
    real alboc(size(rlat)) ! albedo obtenu (de 0 a 1)

    ! Local:

    INTEGER, PARAMETER:: npts =120
    ! Contr\^ole la pr\'ecision de l'int\'egration. 120 correspond \`a
    ! l'intervalle 6 minutes.

    REAL lonsun, declin
    REAL rmu, srmu, salb, aa, bb
    INTEGER i, k

    !----------------------------------------------------------------------

    ! Calculer la longitude vraie de l'orbite terrestre:
    CALL orbite(real(jour), lonsun)

    ! Calculer la declinaison du soleil (qui varie entre + R_incl et - R_incl) :
    declin = asin(sin(lonsun * pi / 180.) * sin(r_incl * pi / 180.))

    DO i = 1, size(rlat)
       aa = sin(rlat(i) * pi / 180.) * sin(declin)
       bb = cos(rlat(i) * pi / 180.) * cos(declin)

       ! Midi local (angle du temps = 0.):
       rmu = max(0., aa + bb)
       srmu = rmu
       salb = 0.058 / (rmu + 0.30) * 1.2 * rmu

       ! Faire l'integration numerique de midi a minuit (le facteur 2
       ! prend en compte l'autre moitie de la journee):
       DO k = 1, npts
          rmu = max(0., aa + bb * cos(real(k) / real(npts) * pi))
          srmu = srmu + rmu * 2.
          salb = salb + 0.058 / (rmu + 0.30) * 1.2 * rmu * 2.
       END DO
       IF (srmu /= 0.) THEN
          alboc(i) = salb / srmu
       ELSE
          ! nuit polaire (on peut prendre une valeur quelconque)
          alboc(i) = 1.
       END IF
    END DO

  END function alboc

end module alboc_m
