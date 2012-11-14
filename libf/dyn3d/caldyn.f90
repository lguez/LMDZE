module caldyn_m

  IMPLICIT NONE

contains

  SUBROUTINE caldyn(itau, ucov, vcov, teta, ps, masse, pk, pkf, phis, phi, &
       dudyn, dv, dteta, dp, w, pbaru, pbarv, time_0, conser)

    ! From dyn3d/caldyn.F, v 1.1.1.1 2004/05/19 12:53:06
    ! Auteur : P. Le Van
    ! Objet : calcul des tendances dynamiques

    use advect_m, only: advect
    USE comgeom, ONLY : airesurg, constang
    USE dimens_m, ONLY : iim, llm
    USE disvert_m, ONLY : ap, bp
    use massdair_m, only: massdair
    USE paramet_m, ONLY : iip1, ip1jm, ip1jmp1, jjp1, llmp1
    use sortvarc_m, only: sortvarc

    ! Arguments:

    LOGICAL, INTENT(IN):: conser
    INTEGER, INTENT(IN):: itau
    REAL vcov(ip1jm, llm), ucov(ip1jmp1, llm)
    real, intent(in):: teta(ip1jmp1, llm)
    REAL, INTENT(IN):: ps(ip1jmp1), phis(ip1jmp1)
    REAL, INTENT(IN):: pk(iip1, jjp1, llm)
    REAL pkf(ip1jmp1, llm)
    REAL vcont(ip1jm, llm), ucont(ip1jmp1, llm)
    REAL, INTENT(IN):: phi(ip1jmp1, llm)
    real masse(ip1jmp1, llm)
    REAL dv(ip1jm, llm), dudyn(ip1jmp1, llm)
    REAL dteta(ip1jmp1, llm)
    real, INTENT(out):: dp(ip1jmp1)
    REAL pbaru(ip1jmp1, llm), pbarv(ip1jm, llm)
    REAL, intent(in):: time_0
    REAL, INTENT(out):: w(ip1jmp1, llm)

    ! Local:

    REAL ang(ip1jmp1, llm), p(ip1jmp1, llmp1)
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

    DO l = 1, llm
       DO ij = 1, ip1jmp1
          ang(ij, l) = ucov(ij, l) + constang(ij)
       END DO
    END DO

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
