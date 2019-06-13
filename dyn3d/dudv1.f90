module dudv1_m

  IMPLICIT NONE

contains

  SUBROUTINE dudv1(vorpot, pbaru, pbarv, du, dv)

    ! From LMDZ4/libf/dyn3d/dudv1.F, version 1.1.1.1, 2004/05/19 12:53:06 

    ! Author: P. Le Van

    ! Objet: calcul du terme de rotation. Ce terme est ajouté à
    ! d(ucov)/dt et à d(vcov)/dt.

    USE dimensions, ONLY: iim, jjm, llm
    USE paramet_m, ONLY: iip1, iip2, ip1jm, ip1jmp1

    REAL, intent(in):: vorpot(ip1jm, llm)
    REAL, intent(in):: pbaru(ip1jmp1, llm), pbarv(ip1jm, llm)
    real, intent(out):: du(iim + 2: (iim + 1) * jjm, llm), dv(ip1jm, llm)

    ! Local:
    INTEGER l, ij

    !----------------------------------------------------------------------

    DO l = 1, llm
       DO ij = iip2, ip1jm - 1
          du(ij, l) = 0.125 * (vorpot(ij - iip1, l) + vorpot(ij, l)) &
               * (pbarv(ij - iip1, l) + pbarv(ij - iim, l) + pbarv(ij, l) &
               + pbarv(ij + 1, l))
       END DO

       DO ij = 1, ip1jm - 1
          dv(ij + 1, l) = - 0.125 * (vorpot(ij, l) + vorpot(ij + 1, l)) &
               * (pbaru(ij, l) + pbaru(ij + 1, l) + pbaru(ij + iip1, l) &
               + pbaru(ij + iip2, l))
       END DO

       ! correction pour dv(1, j, l) 
       ! dv(1, j, l) = dv(iip1, j, l) 
       DO ij = 1, ip1jm, iip1
          dv(ij, l) = dv(ij + iim, l)
       END DO
    END DO

  END SUBROUTINE dudv1

end module dudv1_m
