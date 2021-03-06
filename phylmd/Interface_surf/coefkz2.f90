module coefkz2_m

  IMPLICIT none

contains

  SUBROUTINE coefkz2(paprs, pplay, t, coefm, coefh)

    ! J'introduit un peu de diffusion sauf dans les endroits o\`u une
    ! forte inversion est pr\'esente. On peut dire que la diffusion
    ! repr\'esente la convection peu profonde.

    use dimphy, only: klev
    use SUPHEC_M, only: RCPD, rd

    REAL, intent(in):: paprs(:, :) ! (knon, klev + 1)
    ! pression a chaque intercouche (en Pa)

    REAL, intent(in):: pplay(:, :) ! (knon, klev)
    ! pression au milieu de chaque couche (en Pa)

    REAL, intent(in):: t(:, :) ! (knon, klev) temperature (K)

    REAL, intent(out):: coefm(:, 2:) ! (knon, 2:klev) coefficient vitesse

    REAL, intent(out):: coefh(:, 2:) ! (knon, 2:klev)
    ! coefficient chaleur et humidite)

    ! Local:

    ! Quelques constantes et options:

    REAL, PARAMETER:: prandtl = 0.4
    REAL, PARAMETER:: kstable = 0.002 ! in s-1
    REAL, PARAMETER:: mixlen = 35.0 ! constante controlant longueur de melange

    REAL, PARAMETER:: seuil = - 0.02
    ! au-dela l'inversion est consideree trop faible

    INTEGER knon ! nombre de points a traiter
    INTEGER i, k
    INTEGER invb(size(paprs, 1)) ! (knon)
    REAL zdthmin(size(paprs, 1)) ! (knon)
    real zdthdp

    !----------------------------------------------------------

    knon = size(paprs, 1)

    ! Chercher la zone d'inversion forte :

    DO i = 1, knon
       invb(i) = klev
       zdthmin(i) = 0.0
    ENDDO

    DO k = 2, klev / 2 - 1
       DO i = 1, knon
          zdthdp = ((t(i, k) - t(i, k + 1)) / (pplay(i, k) - pplay(i, k + 1)) &
               - RD * 0.5 * (t(i, k) + t(i, k + 1)) / RCPD / paprs(i, k + 1)) &
               * 100.
          IF (pplay(i, k) > 0.8 * paprs(i, 1) .AND. zdthdp < zdthmin(i)) THEN
             zdthmin(i) = zdthdp
             invb(i) = k
          ENDIF
       ENDDO
    ENDDO

    ! Introduire une diffusion :
    DO i = 1, knon
       IF (invb(i) == klev .OR. zdthmin(i) > seuil) THEN
          ! S'il n'y a pas d'inversion ou si l'inversion est trop
          ! faible :
          coefm(i, :) = (mixlen * MAX(0.0, (paprs(i, 2:klev) &
               - paprs(i, klev + 1)) / (paprs(i, 2) - paprs(i, klev + 1))))**2 &
               * kstable
          coefh(i, :) = coefm(i, :) / prandtl ! h et m different
       else
          coefm(i, :) = 0. 
          coefh(i, :) = 0. 
       ENDIF
    ENDDO

  END SUBROUTINE coefkz2

end module coefkz2_m
