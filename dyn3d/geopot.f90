module geopot_m

  IMPLICIT NONE

contains

  SUBROUTINE geopot(teta, pk, pks, phi)

    ! From libf/dyn3d/geopot.F, version 1.1.1.1 2004/05/19
    ! Author: P. Le Van
    ! Objet : calcul du géopotentiel aux milieux des couches
    ! L'intégration se fait de bas en haut.

    use jumble, only: assert

    USE dimensions, ONLY: iim, jjm, llm
    use grid_noro_m, only: phis

    REAL, INTENT(IN):: teta(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: pk(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: pks(:, :) ! (iim + 1, jjm + 1)
    REAL, INTENT(out)::  phi(:, :, :) ! (iim + 1, jjm + 1, llm)

    ! Local:
    INTEGER l

    ! -----------------------------------------------------------------------

    call assert((/size(teta, 1), size(pk, 1), size(pks, 1), &
         size(phi, 1)/) == iim + 1, "geopot iim")
    call assert((/size(teta, 2), size(pk, 2), size(pks, 2), &
         size(phi, 2)/) == jjm + 1, "geopot jjm")
    call assert((/size(teta, 3), size(pk, 3), size(phi, 3)/) == llm, &
         "geopot llm")

    ! Calcul de phi au niveau 1 près du sol :
    phi(:, :, 1) = phis + teta(:, :, 1) * (pks - pk(:, :, 1))

    ! Calcul de phi aux niveaux supérieurs :
    DO l = 2, llm
       phi(:, :, l) = phi(:, :, l-1) + 0.5 * (teta(:, :, l) + teta(:, :, l-1)) &
            * (pk(:, :, l-1) - pk(:, :, l))
    END DO

  END SUBROUTINE geopot

end module geopot_m
