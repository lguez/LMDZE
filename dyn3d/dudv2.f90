module dudv2_m

  IMPLICIT NONE

contains

  SUBROUTINE dudv2(teta, pkf, bern, du, dv)

    ! From LMDZ4/libf/dyn3d/dudv2.F, version 1.1.1.1, 2004/05/19 12:53:06

    ! Author: P. Le Van

    ! Objet : calcul du terme de pression (gradient de pression /
    ! densité) et du terme "- gradient de la fonction de
    ! Bernoulli". Ces termes sont ajoutés à d(ucov)/dt et à
    ! d(vcov)/dt.

    USE dimensions, ONLY: iim, jjm

    REAL, INTENT(IN):: teta(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: pkf(:, :, :) ! (iim + 1, jjm + 1, llm)
    real, INTENT(IN):: bern(:, :, :) ! (iim + 1, jjm + 1, llm)
    real, intent(inout):: du(:, 2:, :) ! (iim + 1, 2:jjm, llm)
    real, intent(inout):: dv(:, :, :) ! (iim + 1, jjm, llm)

    ! Local:
    INTEGER i, j

    !-----------------------------------------------------------------

    forall (i = 1:iim) du(i, :, :) = du(i, :, :) + 0.5 &
         * (teta(i, 2:jjm, :) + teta(i + 1, 2:jjm, :)) * (pkf(i, 2:jjm, :) &
         - pkf(i + 1, 2:jjm, :)) + bern(i, 2:jjm, :) - bern(i + 1, 2:jjm, :)

    du(iim + 1, :, :) = du(1, :, :)

    forall (j = 1:jjm) dv(:, j, :) = dv(:, j, :) + 0.5 * (teta(:, j, :) &
         + teta(:, j + 1, :)) * (pkf(:, j + 1, :) - pkf(:, j, :)) &
         + bern(:, j + 1, :) - bern(:, j, :)

  END SUBROUTINE dudv2

end module dudv2_m
