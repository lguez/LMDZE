module convmas_m

  IMPLICIT NONE

contains

  function convmas(pbaru, pbarv)

    ! From LMDZ4/libf/dyn3d/convmas.F, version 1.1.1.1, 2004/05/19 12:53:07

    use convflu_m, only: convflu
    USE dimensions, ONLY: iim, jjm, llm
    USE paramet_m, ONLY: ip1jm, ip1jmp1, llmm1
    USE filtreg_scal_m, ONLY: filtreg_scal

    ! Authors: P. Le Van, F. Hourdin
    ! Objet: calcul de la convergence du flux de masse aux niveaux p

    ! Le calcul se fait de haut en bas, la convergence de masse au
    ! niveau p(llm + 1) est égale à 0 et n'est pas stockée dans le
    ! tableau convmas.

    REAL, INTENT(IN):: pbaru(ip1jmp1, llm), pbarv(ip1jm, llm)
    REAL convmas(iim + 1, jjm + 1, llm)

    ! Local:
    INTEGER l

    !-----------------------------------------------------------------------

    ! Calcul de - (d(pbaru)/dx + d(pbarv)/dy) :
    CALL convflu(pbaru, pbarv, llm, convmas)

    ! Filtrage :
    CALL filtreg_scal(convmas, direct = .true., intensive = .false.)

    ! Intégration de la convergence de masse de haut en bas :
    DO l = llmm1, 1, -1
       convmas(:, :, l) = convmas(:, :, l) + convmas(:, :, l + 1)
    END DO

  END function convmas

end module convmas_m
