module caldyn_m

  IMPLICIT NONE

contains

  SUBROUTINE caldyn(itau, ucov, vcov, teta, ps, masse, pk, pkf, phis, phi, &
       dudyn, dv, dteta, dp, w, pbaru, pbarv, time_0, conser)

    ! From dyn3d/caldyn.F, version 1.1.1.1 2004/05/19 12:53:06
    ! Auteur : P. Le Van
    ! Objet : calcul des tendances dynamiques

    use advect_m, only: advect
    USE comgeom, ONLY: airesurg, constang_2d
    USE dimens_m, ONLY: iim, jjm, llm
    USE disvert_m, ONLY: ap, bp
    use flumass_m, only: flumass
    use massbarxy_m, only: massbarxy
    use massdair_m, only: massdair
    USE paramet_m, ONLY: iip1, ip1jm, ip1jmp1, jjp1, llmp1
    use sortvarc_m, only: sortvarc
    use tourpot_m, only: tourpot

    ! Arguments:
    INTEGER, INTENT(IN):: itau
    REAL, INTENT(IN):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm) vent covariant
    REAL, INTENT(IN):: vcov(:, :, :) ! (iim + 1, jjm, llm) ! vent covariant
    real, intent(in):: teta(ip1jmp1, llm)
    REAL, INTENT(IN):: ps(ip1jmp1)
    real, intent(out):: masse(ip1jmp1, llm)
    REAL, INTENT(IN):: pk(iip1, jjp1, llm)
    REAL pkf(ip1jmp1, llm)
    REAL, INTENT(IN):: phis(ip1jmp1)
    REAL, INTENT(IN):: phi(ip1jmp1, llm)
    REAL dudyn(ip1jmp1, llm), dv(ip1jm, llm)
    REAL dteta(ip1jmp1, llm)
    real, INTENT(out):: dp(ip1jmp1)
    REAL, INTENT(out):: w(ip1jmp1, llm)
    REAL pbaru(ip1jmp1, llm), pbarv(ip1jm, llm)
    REAL, intent(in):: time_0
    LOGICAL, INTENT(IN):: conser

    ! Local:

    REAL vcont(ip1jm, llm), ucont(ip1jmp1, llm)
    REAL ang(iim + 1, jjm + 1, llm), p(ip1jmp1, llmp1)
    REAL massebx(ip1jmp1, llm), masseby(ip1jm, llm)
    REAL vorpot(ip1jm, llm)
    real ecin(ip1jmp1, llm), convm(ip1jmp1, llm)
    REAL bern(ip1jmp1, llm)
    REAL massebxy(ip1jm, llm)

    INTEGER ij, l

    !-----------------------------------------------------------------------

    CALL covcont(llm, ucov, vcov, ucont, vcont)
    forall (l = 1: llm + 1) p(:, l) = ap(l) + bp(l) * ps
    CALL massdair(p, masse)
    CALL massbar(masse, massebx, masseby)
    CALL massbarxy(masse, massebxy)
    CALL flumass(massebx, masseby, vcont, ucont, pbaru, pbarv)
    CALL dteta1(teta, pbaru, pbarv, dteta)
    CALL convmas(pbaru, pbarv, convm)
    dp = convm(:, 1) / airesurg
    CALL vitvert(convm, w)
    CALL tourpot(vcov, ucov, massebxy, vorpot)
    CALL dudv1(vorpot, pbaru, pbarv, dudyn, dv)
    CALL enercin(vcov, ucov, vcont, ucont, ecin)
    CALL bernoui(ip1jmp1, llm, phi, ecin, bern)
    CALL dudv2(teta, pkf, bern, dudyn, dv)

    forall (l = 1: llm) ang(:, :, l) = ucov(:, :, l) + constang_2d
    CALL advect(ang, vcov, teta, w, massebx, masseby, dudyn, dv, dteta)

    ! WARNING probleme de peridocite de dv sur les PC/linux. Pb d'arrondi
    ! probablement. Observe sur le code compile avec pgf90 3.0-1
    DO l = 1, llm
       DO ij = 1, ip1jm, iip1
          IF (dv(ij, l)/=dv(ij+iim, l)) THEN
             dv(ij+iim, l) = dv(ij, l)
          END IF
       END DO
    END DO

    ! Sorties eventuelles des variables de controle :
    IF (conser) CALL sortvarc(itau, ucov, teta, ps, masse, pk, phis, vorpot, &
         phi, bern, dp, time_0)

  END SUBROUTINE caldyn

end module caldyn_m
