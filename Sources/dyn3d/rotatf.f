module rotatf_m

  IMPLICIT NONE

contains

  SUBROUTINE rotatf(klevel, x, y, rot)

    ! From LMDZ4/libf/dyn3d/rotatf.F, version 1.1.1.1 2004/05/19 12:53:05

    ! Author: P.Le Van

    ! Calcule le rotationnel à tous les niveaux d'un vecteur de
    ! composantes covariantes x et y.

    USE dimens_m
    USE paramet_m
    USE comgeom
    USE filtreg_m, ONLY: filtreg


    ! .....  variables en arguments  ......

    INTEGER, INTENT (IN) :: klevel
    ! klevel, x  et y   sont des arguments d'entree pour le s-prog
    ! rot          est  un argument  de sortie pour le s-prog

    REAL rot(iim + 1, jjm, klevel)
    REAL, INTENT (IN) :: x(iim + 1, jjm + 1, klevel), y(iim + 1, jjm, klevel)

    ! ...   variables  locales  ...

    INTEGER l, i, j


    DO l = 1, klevel
       forall (i = 1:iim, j = 1:jjm) rot(i, j, l) = y(i + 1, j, l) &
            - y(i, j, l) + x(i, j + 1, l) - x(i, j, l)
       rot(iim + 1, :, l) = rot(1, :, l)
    END DO

    CALL filtreg(rot, direct = .true., intensive = .false.)

    DO l = 1, klevel
       rot(:, :, l) = rot(:, :, l)*unsairez_2d
    END DO

  END SUBROUTINE rotatf

end module rotatf_m
