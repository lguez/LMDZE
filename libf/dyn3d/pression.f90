module pression_m

  IMPLICIT NONE

contains

  SUBROUTINE pression(ngrid, ap, bp, ps, p)

    ! From dyn3d/pression.F, version 1.1.1.1 2004/05/19 12:53:07

    ! Authors : P. Le Van, F. Hourdin

    ! Calcule la pression "p(l)" aux différents niveaux, de l = 1 (niveau du
    ! sol) à l = llm +1.
    ! Ces niveaux correspondent aux interfaces des "llm" couches, avec :
    ! p(:, llm +1) = 0
    ! et :
    ! p(:, 1) = ps

    use dimens_m, only: llm

    INTEGER, intent(in):: ngrid
    REAL, intent(in):: ap(llm + 1), bp(llm + 1)
    real, intent(in):: ps(ngrid)
    real, intent(out):: p(ngrid, llm + 1) 

    INTEGER l

    !---------------------

    forall(l = 1: llm + 1) p(:, l) = ap(l) + bp(l) * ps(:)

  END SUBROUTINE pression

end module pression_m
