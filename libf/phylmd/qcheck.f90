module qcheck_m

  IMPLICIT none

contains

  FUNCTION qcheck(klon, klev, paprs, q, ql, aire)

    ! From phylmd/physiq.F, v 1.22 2006/02/20 09:38:28

    use YOMCST

    ! Calculer et imprimer l'eau totale. A utiliser pour verifier
    ! la conservation de l'eau

    INTEGER klon, klev
    REAL, intent(in):: paprs(klon, klev+1)
    real q(klon, klev), ql(klon, klev)
    REAL aire(klon)
    REAL qtotal, zx, qcheck
    INTEGER i, k

    zx = 0.0
    DO i = 1, klon
       zx = zx + aire(i)
    ENDDO
    qtotal = 0.0
    DO k = 1, klev
       DO i = 1, klon
          qtotal = qtotal + (q(i, k)+ql(i, k)) * aire(i) &
               *(paprs(i, k)-paprs(i, k+1))/RG
       ENDDO
    ENDDO

    qcheck = qtotal/zx

  END FUNCTION qcheck

end module qcheck_m
