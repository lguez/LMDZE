module convmas_m

  IMPLICIT NONE

contains

  SUBROUTINE convmas(pbaru, pbarv, convm)

    ! From LMDZ4/libf/dyn3d/convmas.F, version 1.1.1.1, 2004/05/19 12:53:07

    USE dimens_m, ONLY: llm
    USE paramet_m, ONLY: ip1jm, ip1jmp1, jjp1, llmm1
    USE filtreg_m, ONLY: filtreg

    ! Authors: P. Le Van, F. Hourdin
    ! Objet: calcul de la convergence du flux de masse aux niveaux p

    ! Le calcul se fait de haut en bas, la convergence de masse au
    ! niveau p(llm+1) est égale à 0 et n'est pas stockée dans le
    ! tableau convm.

    REAL, INTENT(IN):: pbaru(ip1jmp1, llm), pbarv(ip1jm, llm)
    REAL, INTENT(OUT):: convm(ip1jmp1, llm)

    ! Local:
    INTEGER l, ij

    !-----------------------------------------------------------------------

    ! Calcul de - (d(pbaru)/dx + d(pbarv)/dy) :
    CALL convflu(pbaru, pbarv, llm, convm)

    ! Filtrage :
    CALL filtreg(convm, jjp1, llm, 2, 2, .TRUE.)

    ! Intégration de la convergence de masse de haut en bas :
    DO l = llmm1, 1, -1
       DO ij = 1, ip1jmp1
          convm(ij, l) = convm(ij, l) + convm(ij, l+1)
       END DO
    END DO

  END SUBROUTINE convmas

end module convmas_m
