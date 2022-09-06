module massdair_m

  IMPLICIT NONE

contains

  pure function massdair(p)

    ! From LMDZ4/libf/dyn3d/massdair.F, version 1.1.1.1 2004/05/19 12:53:07

    ! Calcule la masse d'air dans chaque maille.
    ! Authors: P. Le Van, F. Hourdin

    USE comgeom, ONLY: airesurg
    USE dimensions, ONLY: iim, jjm, llm

    REAL, intent(in):: p(:, :, :) ! (iim + 1, jjm + 1, llm + 1)
    ! pressure at layer interface, in Pa

    real massdair(iim + 1, jjm + 1, llm) ! in kg

    ! Local:
    INTEGER l

    !----------------------------------------------------------

    forall (l = 1: llm) massdair(:iim, :, l) = airesurg(:iim, :) &
         * (p(:iim, :, l) - p(:iim, :, l + 1))
    massdair(iim + 1, :, :) = massdair(1, :, :)

  END function massdair

end module massdair_m
