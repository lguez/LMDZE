module convmas_m

  IMPLICIT NONE

contains

  SUBROUTINE convmas(pbaru, pbarv, convm)

    ! From LMDZ4/libf/dyn3d/convmas.F, version 1.1.1.1, 2004/05/19 12:53:07

    USE dimens_m, ONLY: iim, jjm, llm
    USE paramet_m, ONLY: ip1jm, ip1jmp1, llmm1
    USE filtreg_scal_m, ONLY: filtreg_scal

    ! Authors: P. Le Van, F. Hourdin
    ! Objet: calcul de la convergence du flux de masse aux niveaux p

    ! Le calcul se fait de haut en bas, la convergence de masse au
    ! niveau p(llm + 1) est égale à 0 et n'est pas stockée dans le
    ! tableau convm.

    REAL, INTENT(IN):: pbaru(ip1jmp1, llm), pbarv(ip1jm, llm)
    REAL, INTENT(OUT):: convm(iim + 1, jjm + 1, llm)

    ! Local:
    INTEGER l

    !-----------------------------------------------------------------------

    ! Calcul de - (d(pbaru)/dx + d(pbarv)/dy) :
    CALL convflu(pbaru, pbarv, llm, convm)

    ! Filtrage :
    CALL filtreg_scal(convm, direct = .true., intensive = .false.)

    ! Intégration de la convergence de masse de haut en bas :
    DO l = llmm1, 1, -1
       convm(:, :, l) = convm(:, :, l) + convm(:, :, l + 1)
    END DO

  END SUBROUTINE convmas

end module convmas_m
