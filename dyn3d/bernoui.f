module bernoui_m

  IMPLICIT NONE

contains

  function bernoui(phi, ecin)

    ! From LMDZ4/libf/dyn3d/bernoui.F, version 1.1.1.1 2004/05/19 12:53:06

    ! Author: P. Le Van

    ! Objet : calcul de la fonction de Bernouilli aux niveaux s.
 
    USE dimens_m, ONLY: iim, jjm, llm
    USE filtreg_m, ONLY: filtreg

    REAL, INTENT(IN):: phi(:, :, :), ecin (:, :, :) ! (iim + 1, jjm + 1, llm)

    REAL bernoui(iim + 1, jjm + 1, llm)
    ! fonction de Bernouilli = filtre de (géopotentiel + énergie cinétique)

    !-----------------------------------------------------------------------

    bernoui = phi + ecin
    CALL filtreg(bernoui, direct = .true., intensive = .true.)

  END function bernoui

end module bernoui_m
