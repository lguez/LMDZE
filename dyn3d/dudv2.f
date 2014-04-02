module dudv2_m

  IMPLICIT NONE

contains

  SUBROUTINE dudv2(teta, pkf, bern, du, dv)

    ! From LMDZ4/libf/dyn3d/dudv2.F, version 1.1.1.1, 2004/05/19 12:53:06

    ! Author: P. Le Van

    ! Objet : calcul du terme de pression (gradient de p / densité) et
    ! du terme "- gradient de la fonction de Bernouilli". Ces termes
    ! sont ajoutés à d(ucov)/dt et à d(vcov)/dt.

    USE dimens_m, ONLY: iim, llm
    USE paramet_m, ONLY: iip1, iip2, ip1jm, ip1jmp1

    REAL, INTENT(IN):: teta(ip1jmp1, llm)
    REAL, INTENT(IN):: pkf(ip1jmp1, llm)
    real, INTENT(IN):: bern(ip1jmp1, llm)
    real du(ip1jmp1, llm), dv(ip1jm, llm)
    ! du et dv sont des arguments de sortie pour le s-pg

    ! Local:
    INTEGER l, ij

    !-----------------------------------------------------------------

    DO l = 1, llm
       DO ij = iip2, ip1jm - 1
          du(ij, l) = du(ij, l) + 0.5 * (teta(ij, l) + teta(ij + 1, l)) &
               * (pkf(ij, l) - pkf(ij + 1, l)) + bern(ij, l) - bern(ij + 1, l)
       END DO

       ! correction pour du(iip1, j, l), j=2, jjm
       ! du(iip1, j, l) = du(1, j, l)
       DO ij = iip1 + iip1, ip1jm, iip1
          du(ij, l) = du(ij - iim, l)
       END DO

       DO ij = 1, ip1jm
          dv(ij, l) = dv(ij, l) + 0.5 * (teta(ij, l) + teta(ij + iip1, l)) &
               * (pkf(ij + iip1, l) - pkf(ij, l)) + bern(ij + iip1, l) &
               - bern(ij, l)
       END DO
    END DO

  END SUBROUTINE dudv2

end module dudv2_m
