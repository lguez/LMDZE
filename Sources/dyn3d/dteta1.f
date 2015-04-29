module dteta1_m

  IMPLICIT NONE

contains

  SUBROUTINE dteta1(teta, pbaru, pbarv, dteta)

    ! From LMDZ4/libf/dyn3d/dteta1.F, version 1.1.1.1, 2004/05/19 12:53:06
    ! Authors: P. Le Van, F. Forget

    ! Calcul du terme de convergence horizontale du flux d'enthalpie
    ! potentielle.

    USE dimens_m, ONLY: iim, jjm, llm
    USE paramet_m, ONLY: iip1, iip2, ip1jm, ip1jmp1, jjp1
    USE filtreg_m, ONLY: filtreg

    REAL, intent(in):: teta(ip1jmp1, llm)
    REAL, intent(in):: pbaru(ip1jmp1, llm), pbarv(ip1jm, llm)
    REAL, intent(out):: dteta(iim + 1, jjm + 1, llm)

    ! Local:
    INTEGER l, ij
    REAL hbyv(ip1jm, llm), hbxu(ip1jmp1, llm)

    !----------------------------------------------------------------

    DO l = 1, llm
       DO ij = iip2, ip1jm - 1
          hbxu(ij, l) = pbaru(ij, l) * 0.5 * (teta(ij, l) + teta(ij + 1, l))
       end DO

       DO ij = iip1+ iip1, ip1jm, iip1
          hbxu(ij, l) = hbxu(ij - iim, l)
       end DO

       DO ij = 1, ip1jm
          hbyv(ij, l)= pbarv(ij, l) * 0.5 * (teta(ij, l) + teta(ij + iip1, l))
       end DO
    end DO

    CALL convflu(hbxu, hbyv, llm, dteta)

    ! stockage dans dh de la convergence horizontale filtrée du flux
    ! d'enthalpie potentielle
    CALL filtreg(dteta, direct = .true., intensive = .false.)

  END SUBROUTINE dteta1

end module dteta1_m
