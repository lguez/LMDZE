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

    USE dimensions, ONLY: iim, jjm, llm

    REAL, INTENT(IN):: teta(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: pkf(:, :, :) ! (iim + 1, jjm + 1, llm)
    real, INTENT(IN):: bern(:, :, :) ! (iim + 1, jjm + 1, llm)
    real, intent(inout):: du(:, :, :) ! (iim + 1, jjm + 1, llm)
    real, intent(inout):: dv(:, :, :) ! (iim + 1, jjm, llm)

    ! Local:
    INTEGER l, i, j

    !-----------------------------------------------------------------

    DO l = 1, llm
       do j = 2, jjm
          do i = 1, iim
             du(i, j, l) = du(i, j, l) + 0.5 * (teta(i, j, l) &
                  + teta(i + 1, j, l)) * (pkf(i, j, l) - pkf(i + 1, j, l)) &
                  + bern(i, j, l) - bern(i + 1, j, l)
          end do
       end do

       do j = 2, jjm
          du(iim + 1, j, l) = du(1, j, l)
       end do

       do j = 1, jjm
          do i = 1, iim + 1
             dv(i, j, l) = dv(i, j, l) + 0.5 * (teta(i, j, l) &
                  + teta(i, j + 1, l)) * (pkf(i, j + 1, l) - pkf(i, j, l)) &
                  + bern(i, j + 1, l) - bern(i, j, l)
          end do
       END DO
    END DO

  END SUBROUTINE dudv2

end module dudv2_m
