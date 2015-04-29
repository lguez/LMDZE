module vitvert_m

  IMPLICIT NONE

contains

  SUBROUTINE vitvert(convm, w)

    ! From libf/dyn3d/vitvert.F, version 1.1.1.1, 2004/05/19 12:53:05
    ! Authors: P. Le Van, F. Hourdin

    ! Objet : calcul de la vitesse verticale aux niveaux sigma

    ! La vitesse verticale est orientée de haut en bas. Au sol, au
    ! niveau sigma(1), w(i, j, 1) = 0. Au sommet, au niveau
    ! sigma(llm+1), la vitesse verticale est aussi égale à 0 et n'est
    ! pas stockée dans le tableau w.

    USE dimens_m, ONLY : llm
    USE disvert_m, ONLY : bp
    USE paramet_m, ONLY : ip1jmp1

    real, intent(in):: convm(ip1jmp1, llm)
    REAL, intent(out):: w(ip1jmp1, llm)

    ! Local:
    INTEGER l

    !------------------------------------------------------

    forall (l = 2: llm) w(:, l) = convm(:, l) - bp(l) * convm(:, 1)
    w(:, 1) = 0.

  END SUBROUTINE vitvert

end module vitvert_m
