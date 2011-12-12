module dissip_m

  IMPLICIT NONE

contains

  SUBROUTINE dissip(vcov, ucov, teta, p, dv, du, dh)

    ! From dyn3d/dissip.F, version 1.1.1.1 2004/05/19 12:53:05
    ! Avec nouveaux opérateurs star : gradiv2, divgrad2, nxgraro2
    ! Author: P. Le Van
    ! Objet : dissipation horizontale

    USE dimens_m, ONLY: iim, jjm, llm
    USE paramet_m, ONLY: iip1, iip2, ip1jmp1, llmp1
    USE comdissnew, ONLY: lstardis, nitergdiv, nitergrot, niterh
    USE inidissip_m, ONLY: dtdiss, tetah, tetaudiv, tetaurot, cdivu, crot, cdivh
    use gradiv2_m, only: gradiv2

    REAL, intent(in):: vcov(:, :, :) ! (iim + 1, jjm, llm)
    REAL, intent(in):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(in):: teta((iim + 1) * (jjm + 1), llm)
    REAL, INTENT(IN):: p((iim + 1) * (jjm + 1), llmp1)
    REAL, intent(out):: dv(:, :, :) ! (iim + 1, jjm, llm)
    REAL, intent(out):: du(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(out):: dh(:, :, :) ! (iim + 1, jjm + 1, llm)

    ! Local:
    REAL gdx(iim + 1, jjm + 1, llm), gdy(iim + 1, jjm, llm)
    REAL grx(iim + 1, jjm + 1, llm), gry(iim + 1, jjm, llm)
    REAL te1dt(llm), te2dt(llm), te3dt(llm)
    REAL deltapres((iim + 1) * (jjm + 1), llm)
    INTEGER l

    !-----------------------------------------------------------------------

    ! Initializations:
    te1dt = tetaudiv * dtdiss
    te2dt = tetaurot * dtdiss
    te3dt = tetah * dtdiss
    du = 0.
    dv = 0.
    dh = 0.

    ! Calcul de la dissipation:

    ! Calcul de la partie grad (div) :

    IF (lstardis) THEN
       CALL gradiv2(llm, ucov, vcov, nitergdiv, gdx, gdy, cdivu)
    ELSE
       CALL gradiv(llm, ucov, vcov, nitergdiv, gdx, gdy, cdivu)
    END IF

    gdx(:, 1, :) = 0.
    gdx(:, jjm + 1, :) = 0.
    forall (l = 1: llm)
       du(:, 2: jjm, l) = du(:, 2: jjm, l) - te1dt(l) * gdx(:, 2: jjm, l)
       dv(:, :, l) = dv(:, :, l) - te1dt(l) * gdy(:, :, l)
    END forall

    ! calcul de la partie n X grad (rot) :

    IF (lstardis) THEN
       CALL nxgraro2(llm, ucov, vcov, nitergrot, grx, gry, crot)
    ELSE
       CALL nxgrarot(llm, ucov, vcov, nitergrot, grx, gry, crot)
    END IF


    grx(:, 1, :) = 0.
    forall (l = 1: llm)
       du(:, 2: jjm, l) = du(:, 2: jjm, l) - te2dt(l) * grx(:, 2: jjm, l)
       dv(:, :, l) = dv(:, :, l) - te2dt(l) * gry(:, :, l)
    END forall

    ! calcul de la partie div (grad) :

    IF (lstardis) THEN
       forall (l = 1: llm) deltapres(:, l) = max(0., p(:, l) - p(:, l + 1))
       CALL divgrad2(llm, teta, deltapres, niterh, gdx, cdivh)
    ELSE
       CALL divgrad(llm, teta, niterh, gdx, cdivh)
    END IF

    forall (l = 1: llm) dh(:, :, l) = dh(:, :, l) - te3dt(l) * gdx(:, :, l)

  END SUBROUTINE dissip

end module dissip_m
