module dissip_m

  IMPLICIT NONE

contains

  SUBROUTINE dissip(vcov, ucov, teta, p, dv, du, dh)

    ! From dyn3d/dissip.F, version 1.1.1.1 2004/05/19 12:53:05
    ! Avec nouveaux opérateurs star : gradiv2, divgrad2, nxgraro2
    ! Author: P. Le Van
    ! Objet : dissipation horizontale

    USE dimens_m, ONLY : iim, jjm, llm
    USE paramet_m, ONLY : iip1, iip2, ip1jmp1, llmp1
    USE comdissnew, ONLY : lstardis, nitergdiv, nitergrot, niterh
    USE inidissip_m, ONLY : dtdiss, tetah, tetaudiv, tetaurot

    ! Arguments:
    REAL vcov((iim + 1) * jjm, llm), ucov(ip1jmp1, llm), teta(ip1jmp1, llm)
    REAL, INTENT (IN) :: p(ip1jmp1, llmp1)
    REAL dv((iim + 1) * jjm, llm), du(ip1jmp1, llm), dh(ip1jmp1, llm)

    ! Local:
    REAL gdx(ip1jmp1, llm), gdy((iim + 1) * jjm, llm)
    REAL grx(ip1jmp1, llm), gry((iim + 1) * jjm, llm)
    REAL te1dt(llm), te2dt(llm), te3dt(llm)
    REAL deltapres(ip1jmp1, llm)

    INTEGER l, ij

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
       CALL gradiv2(llm, ucov, vcov, nitergdiv, gdx, gdy)
    ELSE
       CALL gradiv(llm, ucov, vcov, nitergdiv, gdx, gdy)
    END IF

    DO l = 1, llm
       DO ij = 1, iip1
          gdx(ij, l) = 0.
          gdx(ij+(iim + 1) * jjm, l) = 0.
       END DO

       DO ij = iip2, (iim + 1) * jjm
          du(ij, l) = du(ij, l) - te1dt(l) * gdx(ij, l)
       END DO
       DO ij = 1, (iim + 1) * jjm
          dv(ij, l) = dv(ij, l) - te1dt(l) * gdy(ij, l)
       END DO
    END DO

    ! calcul de la partie n X grad (rot) :

    IF (lstardis) THEN
       CALL nxgraro2(llm, ucov, vcov, nitergrot, grx, gry)
    ELSE
       CALL nxgrarot(llm, ucov, vcov, nitergrot, grx, gry)
    END IF


    DO l = 1, llm
       DO ij = 1, iip1
          grx(ij, l) = 0.
       END DO

       DO ij = iip2, (iim + 1) * jjm
          du(ij, l) = du(ij, l) - te2dt(l) * grx(ij, l)
       END DO
       DO ij = 1, (iim + 1) * jjm
          dv(ij, l) = dv(ij, l) - te2dt(l) * gry(ij, l)
       END DO
    END DO

    ! calcul de la partie div (grad) :

    IF (lstardis) THEN
       DO l = 1, llm
          DO ij = 1, ip1jmp1
             deltapres(ij, l) = max(0., p(ij, l) - p(ij, l + 1))
          END DO
       END DO

       CALL divgrad2(llm, teta, deltapres, niterh, gdx)
    ELSE
       CALL divgrad(llm, teta, niterh, gdx)
    END IF

    forall (l = 1: llm) dh(:, l) = dh(:, l) - te3dt(l) * gdx(:, l)

  END SUBROUTINE dissip

end module dissip_m
