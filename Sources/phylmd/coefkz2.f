module coefkz2_m

  IMPLICIT none

contains

  SUBROUTINE coefkz2(nsrf, paprs, pplay, t, pcfm, pcfh)

    ! J'introduit un peu de diffusion sauf dans les endroits o\`u une
    ! forte inversion est pr\'esente. On peut dire que la diffusion
    ! repr\'esente la convection peu profonde.

    use indicesol, only: is_oce
    use dimphy, only: klev
    use SUPHEC_M, only: RCPD, rd

    integer, intent(in):: nsrf ! indicateur de la nature du sol

    REAL, intent(in):: paprs(:, :) ! (knon, klev+1)
    ! pression a chaque intercouche (en Pa)

    REAL, intent(in):: pplay(:, :) ! (knon, klev)
    ! pression au milieu de chaque couche (en Pa)
    
    REAL, intent(in):: t(:, :) ! (knon, klev) temperature (K)

    REAL, intent(out):: pcfm(:, 2:) ! (knon, 2:klev) coefficient vitesse

    REAL, intent(out):: pcfh(:, 2:) ! (knon, 2:klev)
    ! coefficient chaleur et humidite)

    ! Local:
    
    ! Quelques constantes et options:

    REAL prandtl
    PARAMETER (prandtl=0.4)
    REAL kstable
    PARAMETER (kstable=0.002)
    REAL mixlen ! constante controlant longueur de melange
    PARAMETER (mixlen=35.0)
    REAL seuil ! au-dela l'inversion est consideree trop faible
    PARAMETER (seuil=-0.02)

    INTEGER knon ! nombre de points a traiter
    INTEGER i, k
    INTEGER invb(size(paprs, 1)) ! (knon)
    REAL zl2(size(paprs, 1)) ! (knon)
    REAL zdthmin(size(paprs, 1)) ! (knon)
    real zdthdp

    !----------------------------------------------------------

    knon = size(paprs, 1)
    
    ! Initialiser les sorties
    DO k = 2, klev
       DO i = 1, knon
          pcfm(i, k) = 0.0
          pcfh(i, k) = 0.0
       ENDDO
    ENDDO

    ! Chercher la zone d'inversion forte

    DO i = 1, knon
       invb(i) = klev
       zdthmin(i)=0.0
    ENDDO
    DO k = 2, klev/ 2 - 1
       DO i = 1, knon
          zdthdp = (t(i, k) - t(i, k + 1)) / (pplay(i, k) - pplay(i, k + 1)) &
               - RD * 0.5 * (t(i, k) + t(i, k + 1)) / RCPD / paprs(i, k + 1)
          zdthdp = zdthdp * 100.
          IF (pplay(i, k) > 0.8 * paprs(i, 1) .AND. zdthdp < zdthmin(i)) THEN
             zdthmin(i) = zdthdp
             invb(i) = k
          ENDIF
       ENDDO
    ENDDO

    ! Introduire une diffusion:
    DO k = 2, klev
       DO i = 1, knon
          ! si on est sur ocean et s'il n'y a pas d'inversion ou si
          ! l'inversion est trop faible:
          IF ((nsrf.EQ.is_oce) .AND.  &  
               ((invb(i).EQ.klev) .OR. (zdthmin(i) > seuil))) THEN
             zl2(i)=(mixlen*MAX(0.0, (paprs(i, k)-paprs(i, klev+1)) &
                  /(paprs(i, 2)-paprs(i, klev+1))))**2
             pcfm(i, k)= zl2(i)* kstable
             pcfh(i, k) = pcfm(i, k) /prandtl ! h et m different
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE coefkz2

end module coefkz2_m
