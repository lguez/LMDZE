module grad_m

  IMPLICIT NONE

contains

  SUBROUTINE grad(klevel, g, gx, gy)

    ! From LMDZ4/libf/dyn3d/grad.F, version 1.1.1.1 2004/05/19 12:53:05
    ! P. Le Van

    ! Calcul des composantes covariantes en x et y du gradient de g.

    USE dimens_m, ONLY: iim, jjm

    INTEGER, intent(in):: klevel
    REAL, intent(in):: g(iim + 1, jjm + 1, klevel)
    REAL, intent(out):: gx(iim + 1, jjm + 1, klevel) , gy(iim + 1, jjm, klevel)

    ! Local:
    INTEGER i, j

    !----------------------------------------------------------------

    forall (i = 1:iim) gx(i, :, :) = g(i + 1, :, :) - g(i, :, :)
    gx(iim + 1, :, :)= gx(1, :, :)

    forall (j = 1:jjm) gy(:, j, :) = g(:, j, :) - g(:, j + 1, :)

  END SUBROUTINE grad

end module grad_m
