module bernoui_m

  IMPLICIT NONE

contains

  SUBROUTINE bernoui(phi, ecin, bern)

    ! From LMDZ4/libf/dyn3d/bernoui.F, version 1.1.1.1 2004/05/19 12:53:06

    ! Author: P. Le Van

    ! Objet: calcul de la fonction de Bernouilli aux niveaux s.
    ! fonction de Bernouilli = bern = filtre de(geopotentiel + energ.cinet.)

    USE dimens_m, ONLY: jjm, llm
    USE filtreg_m, ONLY: filtreg

    REAL, INTENT(IN):: phi(:, :, :), ecin (:, :, :)! (iim + 1, jjm + 1, llm)
    REAL, intent(out):: bern(:, :, :) ! (iim + 1, jjm + 1, llm)

    !-----------------------------------------------------------------------

    bern = phi + ecin
    CALL filtreg(bern, direct = .true., intensive = .true.)

  END SUBROUTINE bernoui

end module bernoui_m
