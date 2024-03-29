module dissip_m

  IMPLICIT NONE

contains

  SUBROUTINE dissip(vcov, ucov, teta, p, dv, du, dh)

    ! From dyn3d/dissip.F, version 1.1.1.1 2004/05/19 12:53:05
    ! Author: P. Le Van
    
    ! Objet : calcul de la dissipation horizontale. Avec op\'erateurs
    ! star : gradiv2, divgrad2, nxgraro2.

    use jumble, only: assert

    USE comdissnew, ONLY: nitergdiv, nitergrot, niterh
    USE dimensions, ONLY: iim, jjm, llm
    use divgrad2_m, only: divgrad2
    use gradiv2_m, only: gradiv2
    USE inidissip_m, ONLY: dtdiss, tetah, tetaudiv, tetaurot, cdivu, crot, cdivh
    use nxgraro2_m, only: nxgraro2

    REAL, intent(in):: vcov(:, :, :) ! (iim + 1, jjm, llm)
    REAL, intent(in):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(in):: teta(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: p(:, :, :) ! (iim + 1, jjm + 1, llm + 1)
    REAL, intent(out):: dv(:, :, :) ! (iim + 1, jjm, llm)
    REAL, intent(out):: du(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(out):: dh(:, :, :) ! (iim + 1, jjm + 1, llm)

    ! Local:
    REAL gdx(iim + 1, jjm + 1, llm), gdy(iim + 1, jjm, llm)
    REAL tedt(llm)
    REAL deltapres(iim + 1, jjm + 1, llm)
    INTEGER l

    !-----------------------------------------------------------------------

    call assert((/size(vcov, 1), size(ucov, 1), size(teta, 1), size(p, 1), &
         size(dv, 1), size(du, 1), size(dh, 1)/) == iim + 1, "dissip iim")
    call assert((/size(vcov, 2), size(ucov, 2) - 1, size(teta, 2) - 1, &
         size(p, 2) - 1, size(dv, 2), size(du, 2) - 1, size(dh, 2) - 1/) &
         == jjm, "dissip jjm")
    call assert((/size(vcov, 3), size(ucov, 3), size(teta, 3), size(p, 3) - 1, &
         size(dv, 3), size(du, 3), size(dh, 3)/) == llm, "dissip llm")

    du(:, 1, :) = 0.
    du(:, jjm + 1, :) = 0.

    ! Calcul de la partie grad(div) :
    CALL gradiv2(ucov, vcov, nitergdiv, gdx, gdy, cdivu)
    tedt = tetaudiv * dtdiss
    forall (l = 1: llm)
       du(:, 2: jjm, l) = - tedt(l) * gdx(:, 2: jjm, l)
       dv(:, :, l) = - tedt(l) * gdy(:, :, l)
    END forall

    ! Calcul de la partie n \wedge grad(rot) :
    CALL nxgraro2(ucov, vcov, nitergrot, gdx, gdy, crot)
    tedt = tetaurot * dtdiss
    forall (l = 1: llm)
       du(:, 2: jjm, l) = du(:, 2: jjm, l) - tedt(l) * gdx(:, 2: jjm, l)
       dv(:, :, l) = dv(:, :, l) - tedt(l) * gdy(:, :, l)
    END forall

    ! calcul de la partie div(grad) :
    forall (l = 1: llm) &
         deltapres(:, :, l) = max(0., p(:, :, l) - p(:, :, l + 1))
    CALL divgrad2(llm, teta, deltapres, niterh, gdx, cdivh)
    forall (l = 1: llm) dh(:, :, l) = - tetah(l) * dtdiss * gdx(:, :, l)

  END SUBROUTINE dissip

end module dissip_m
