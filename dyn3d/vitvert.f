module vitvert_m

  IMPLICIT NONE

contains

  pure function vitvert(convm)

    ! From libf/dyn3d/vitvert.F, version 1.1.1.1, 2004/05/19 12:53:05
    ! Authors: P. Le Van, F. Hourdin

    ! Purpose: Compute vertical speed at sigma levels.

    ! Vertical speed is oriented from bottom to top. At ground-level
    ! sigma(1): vitvert(i, j, 1) = 0. At top-level sigma(llm + 1), vertical
    ! speed is 0 too and is not stored in vitvert.
    
    USE dimens_m, ONLY : llm
    USE disvert_m, ONLY : bp

    real, intent(in):: convm(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL vitvert(size(convm, 1), size(convm, 2), size(convm, 3))

    ! Local:
    INTEGER l

    !------------------------------------------------------

    forall (l = 2: llm) &
         vitvert(:, :, l) = convm(:, :, l) - bp(l) * convm(:, :, 1)
    vitvert(:, :, 1) = 0.

  END function vitvert

end module vitvert_m
