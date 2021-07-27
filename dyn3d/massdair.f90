module massdair_m

  IMPLICIT NONE

contains

  SUBROUTINE massdair(p, masse)

    ! From LMDZ4/libf/dyn3d/massdair.F, version 1.1.1.1 2004/05/19 12:53:07

    ! Calcule la masse d'air dans chaque maille.
    ! Authors: P. Le Van, F. Hourdin

    USE comgeom, ONLY: airesurg
    USE dimensions, ONLY: iim, jjm, llm

    REAL, intent(in):: p(iim + 1, jjm + 1, llm + 1)
    ! aux interfaces des llm couches

    real, intent(out):: masse(iim + 1, jjm + 1, llm)

    ! Local:
    INTEGER l

    !----------------------------------------------------------

    forall (l = 1: llm) masse(:iim, :, l) = airesurg(:iim, :) &
         * (p(:iim, :, l) - p(:iim, :, l + 1))
    masse(iim + 1, :, :) = masse(1, :, :)

  END SUBROUTINE massdair

end module massdair_m
