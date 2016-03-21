module gr_fi_ecrit_m

  IMPLICIT none

contains

  SUBROUTINE gr_fi_ecrit(nfield, nlon, iim, jjmp1, fi, ecrit)

    ! From phylmd/physiq.F, version 1.22 2006/02/20 09:38:28

    ! Transforme une variable de la grille physique \`a la grille
    ! d'\'ecriture.  Cf. version moderne "gr_phy_write_2d", dans le
    ! cas o\`u "nfield" vaut 1.

    INTEGER, intent(in):: nfield, nlon, iim, jjmp1
    REAL, intent(in):: fi(nlon, nfield)
    real ecrit(iim*jjmp1, nfield)

    ! Variables local to the procedure:

    integer jjm
    INTEGER i, n, ig

    !---------------

    jjm = jjmp1 - 1
    DO n = 1, nfield
       DO i=1, iim
          ecrit(i, n) = fi(1, n)
          ecrit(i+jjm*iim, n) = fi(nlon, n)
       ENDDO
       DO ig = 1, nlon - 2
          ecrit(iim+ig, n) = fi(1+ig, n)
       ENDDO
    ENDDO

  END SUBROUTINE gr_fi_ecrit

end module gr_fi_ecrit_m
