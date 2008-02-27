SUBROUTINE gr_fi_ecrit(nfield, nlon, iim, jjmp1, fi, ecrit)

  ! From phylmd/physiq.F, version 1.22 2006/02/20 09:38:28
  ! This procedure is clean: no C preprocessor directive, no include line.

  IMPLICIT none

  ! Transforme une variable de la grille physique à la grille d'écriture.
  ! Cf. version moderne "gr_phy_write", dans le cas où "nfield" vaut 1.

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
