module qcheck_m

  IMPLICIT none

contains

  pure FUNCTION qcheck(paprs, q, ql)

    ! From phylmd/physiq.F, v 1.22 2006/02/20 09:38:28

    ! Calculer et imprimer l'eau totale. A utiliser pour vérifier
    ! la conservation de l'eau.

    use comgeomphy, only: airephy
    use dimphy, only: klon, klev
    use SUPHEC_M, ONLY: rg

    REAL, intent(in):: paprs(:, :) ! (klon, klev + 1)
    real, intent(in):: q(:, :), ql(:, :) ! (klon, klev)

    ! Local:
    REAL qtotal, zx, qcheck
    INTEGER i, k

    !---------------------------------------------------------

    zx = 0.0
    DO i = 1, klon
       zx = zx + airephy(i)
    ENDDO
    qtotal = 0.0
    DO k = 1, klev
       DO i = 1, klon
          qtotal = qtotal + (q(i, k)+ql(i, k)) * airephy(i) &
               *(paprs(i, k)-paprs(i, k+1))/RG
       ENDDO
    ENDDO

    qcheck = qtotal / zx

  END FUNCTION qcheck

end module qcheck_m
