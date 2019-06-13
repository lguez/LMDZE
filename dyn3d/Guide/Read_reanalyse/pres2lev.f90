module pres2lev_m

  IMPLICIT NONE

contains

  SUBROUTINE pres2lev(varo, varn, po, pn)

    ! From LMDZ4/libf/dyn3d/pres2lev.F, version 1.1.1.1 2004/05/19 12:53:07

    ! Interpolation lin\'eaire pour passer \`a une nouvelle
    ! discr\'etisation verticale pour les variables de GCM.

    ! Francois Forget (January 1995)

    REAL, intent(in):: varo(:, :, :) ! (ni, nj, lmo) var in the old grid
    REAL, intent(out):: varn(:, :, :) ! (ni, nj, lmn)! var in the new grid

    REAL, intent(in):: po(:) ! (lmo) 
    ! pressure levels, old  (in monotonic order), in hPa

    REAL, intent(in):: pn(:, :, :) ! (ni, nj, lmn) pressure levels, new, in Pa

    ! Local:
    INTEGER lmn ! dimensions nouvelle couches
    INTEGER ni, nj
    INTEGER lmo ! dimensions ancienne couches
    INTEGER i, j
    REAL zvaro(size(po))
    real zpo(size(po)) ! pressure levels, old, in descending order, in hPa
    INTEGER ln, lo

    !--------------------------------------------------------------

    lmo = size(po)
    ni = size(varn, 1)
    nj = size(varn, 2)
    lmn = size(varn, 3)

    DO i = 1, ni
       DO j = 1, nj
          if (po(1) < po(2)) then
             ! Inversion de l'ordre des niveaux verticaux :
             zpo = po(lmo:1:- 1)
             zvaro = varo(i, j, lmo:1:- 1)
          else
             zpo = po
             zvaro = varo(i, j, :)
          end if

          DO ln = 1, lmn
             IF (pn(i, j, ln) >= zpo(1)) THEN
                varn(i, j, ln) = zvaro(1)
             ELSE IF (pn(i, j, ln) <= zpo(lmo)) THEN
                varn(i, j, ln) = zvaro(lmo)
             ELSE
                DO lo = 1, lmo - 1
                   IF ((pn(i, j, ln) <= zpo(lo)) &
                        .AND. (pn(i, j, ln) > zpo(lo + 1))) THEN
                      varn(i, j, ln) = zvaro(lo) + (pn(i, j, ln) - zpo(lo)) &
                           / (zpo(lo + 1) - zpo(lo)) * (zvaro(lo + 1) &
                           - zvaro(lo))
                   END IF
                END DO
             END IF
          END DO
       END DO
    END DO

  END SUBROUTINE pres2lev

end module pres2lev_m
