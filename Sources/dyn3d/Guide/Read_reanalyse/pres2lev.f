module pres2lev_m

  IMPLICIT NONE

contains

  SUBROUTINE pres2lev(varo, varn, lmo, lmn, po, pn, ni, nj)

    ! From LMDZ4/libf/dyn3d/pres2lev.F, v 1.1.1.1 2004/05/19 12:53:07

    ! interpolation lineaire pour passer
    ! a une nouvelle discretisation verticale pour
    ! les variables de GCM
    ! Francois Forget (01/1995)

    ! MOdif remy roca 12/97 pour passer de pres2sig

    ! Declarations:

    ! ARGUMENTS
    ! """""""""

    INTEGER, intent(in):: lmo ! dimensions ancienne couches (input)
    INTEGER lmn ! dimensions nouvelle couches (input)
    INTEGER lmomx ! dimensions ancienne couches (input)
    INTEGER lmnmx ! dimensions nouvelle couches (input)

    PARAMETER (lmomx=10000, lmnmx=10000)

    REAL, intent(in):: po(lmo) ! niveau de pression en millibars
    INTEGER ni, nj
    REAL pn(ni, nj, lmn) ! niveau de pression en pascals

    INTEGER i, j ! nombre de point horizontale (input)

    REAL varo(ni, nj, lmo) ! var dans l'ancienne grille (input)
    REAL varn(ni, nj, lmn) ! var dans la nouvelle grille (output)

    REAL zvaro(lmomx), zpo(lmomx)

    ! Autres variables
    ! """"""""""""""""
    INTEGER ln, lo
    REAL coef

    ! run

    DO i = 1, ni
       DO j = 1, nj
          ! a chaque point de grille correspond un nouveau sigma old
          ! qui vaut pres(l)/ps(i, j)
          DO lo = 1, lmo
             zpo(lo) = po(lmo+1-lo)
             zvaro(lo) = varo(i, j, lmo+1-lo)
          END DO

          DO ln = 1, lmn
             IF (pn(i, j, ln)>=zpo(1)) THEN
                varn(i, j, ln) = zvaro(1)
             ELSE IF (pn(i, j, ln)<=zpo(lmo)) THEN
                varn(i, j, ln) = zvaro(lmo)
             ELSE
                DO lo = 1, lmo - 1
                   IF ((pn(i, j, ln)<=zpo(lo)) .AND. (pn(i, j, ln)>zpo(lo+1))) THEN
                      coef = (pn(i, j, ln)-zpo(lo))/(zpo(lo+1)-zpo(lo))
                      varn(i, j, ln) = zvaro(lo) + coef*(zvaro(lo+1)-zvaro(lo))
                      ! print*, 'pn(', ln, ')=', pn(i, j, ln), varn(i, j, ln)
                   END IF
                END DO
             END IF
          END DO

       END DO
    END DO
    RETURN
  END SUBROUTINE pres2lev

end module pres2lev_m
