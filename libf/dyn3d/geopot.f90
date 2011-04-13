module geopot_m

  IMPLICIT NONE

contains

  SUBROUTINE geopot(ngrid, teta, pk, pks, phis, phi)

    ! From libf/dyn3d/geopot.F, version 1.1.1.1 2004/05/19
    ! Author: P. Le Van
    ! Objet : calcul du g�opotentiel aux milieux des couches
    ! L'int�gration se fait de bas en haut.

    USE dimens_m, ONLY: llm

    INTEGER, INTENT (IN):: ngrid
    REAL, INTENT (IN):: teta(ngrid, llm), pks(ngrid)
    REAL, INTENT (IN) :: phis(ngrid)
    REAL, INTENT (IN) :: pk(ngrid, llm)
    REAL, INTENT (out)::  phi(ngrid, llm)

    ! Local:
    INTEGER l

    ! -----------------------------------------------------------------------

    ! Calcul de phi au niveau 1 pr�s du sol :
    phi(:, 1) = phis + teta(:, 1) * (pks - pk(:, 1))

    ! Calcul de phi aux niveaux sup�rieurs :
    DO l = 2, llm
       phi(:, l) = phi(:, l-1) + 0.5 * (teta(:, l) + teta(:, l-1)) &
            * (pk(:, l-1) - pk(:, l))
    END DO

  END SUBROUTINE geopot

end module geopot_m
