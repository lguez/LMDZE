module covcont_m

  IMPLICIT NONE

contains

  SUBROUTINE covcont(ucov, vcov, ucont, vcont)

    ! From LMDZ4/libf/dyn3d/covcont.F, version 1.1.1.1, 2004/05/19 12:53:07

    ! Author:  P. Le Van

    ! Objet : calcul des composantes contravariantes \`a partir des
    ! composantes covariantes

    use dimensions, only: jjm, llm
    USE comgeom, only: unscu2, unscv2

    REAL, INTENT(IN):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm) vent covariant
    REAL, INTENT(IN):: vcov(:, :, :) ! (iim + 1, jjm, llm) vent covariant
    REAL, INTENT(out):: ucont(:, :, :) ! (iim + 1, jjm + 1, llm)
    real, INTENT(out):: vcont(:, :, :) ! (iim + 1, jjm, llm)

    ! Local:
    INTEGER l

    !-------------------------------------------------------------------

    forall (l = 1:llm)
       ucont(:, 2:jjm, l) = ucov(:, 2:jjm, l) * unscu2(:, 2:jjm)
       vcont(:, :, l) = vcov(:, :, l) * unscv2
    END forall

    ucont(:, 1, :) = 0.
    ucont(:, jjm + 1, :) = 0.

  END SUBROUTINE covcont

end module covcont_m
