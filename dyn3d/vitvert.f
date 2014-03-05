SUBROUTINE vitvert(convm , w)

  ! From libf/dyn3d/vitvert.F, version 1.1.1.1 2004/05/19 12:53:05
  ! Authors: P. Le Van , F. Hourdin

  ! Objet : calcul de la vitesse verticale aux niveaux sigma

  ! La vitesse verticale est orientee de haut en bas .
  ! au sol, au niveau sigma(1), w(i, j, 1) = 0.
  ! au sommet, au niveau sigma(llm+1) , la vit.verticale est aussi
  ! egale a 0. et n'est pas stockee dans le tableau w .

  USE dimens_m, ONLY : llm
  USE paramet_m, ONLY : ip1jmp1
  USE disvert_m, ONLY : bp

  IMPLICIT NONE

  real, intent(in):: convm(ip1jmp1, llm)
  REAL, intent(out):: w(ip1jmp1, llm)

  ! Local:
  INTEGER l

  !------------------------------------------------------

  forall (l = 2: llm) w(:, l) = convm(:, l) - bp(l) * convm(:, 1)
  w(:, 1) = 0.

END SUBROUTINE vitvert
