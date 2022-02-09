module caldyn_m

  IMPLICIT NONE

contains


  SUBROUTINE caldyn(itau, ucov, vcov, teta, ps, masse, pk, pkf, phi, du, dv, &
       dteta, dp, w, pbaru, pbarv, conser)

    ! From dyn3d/caldyn.F, version 1.1.1.1, 2004/05/19 12:53:06
    ! Author: P. Le Van
    ! Objet : calcul des tendances dynamiques

    use advect_m, only: advect
    use bernoui_m, only: bernoui
    use caldyn0_m, only: ang0, etot0, ptot0, stot0, ztot0
    USE comconst, ONLY: daysec
    USE comgeom, ONLY: airesurg, constang_2d
    USE conf_gcm_m, ONLY: day_step, dtvr
    use convmas_m, only: convmas
    use covcont_m, only: covcont
    USE dimensions, ONLY: iim, jjm, llm
    USE disvert_m, ONLY: ap, bp
    use dteta1_m, only: dteta1
    use dudv1_m, only: dudv1
    use dudv2_m, only: dudv2
    USE dynetat0_m, ONLY: day_ini
    use enercin_m, only: enercin
    use flumass_m, only: flumass
    use massbar_m, only: massbar
    use massbarxy_m, only: massbarxy
    use massdair_m, only: massdair
    USE paramet_m, ONLY: iip1, ip1jmp1, jjp1, llmp1
    use sortvarc_m, only: sortvarc
    use tourpot_m, only: tourpot
    use vitvert_m, only: vitvert

    INTEGER, INTENT(IN):: itau
    REAL, INTENT(IN):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm) vent covariant
    REAL, INTENT(IN):: vcov(:, :, :) ! (iim + 1, jjm, llm) vent covariant
    REAL, INTENT(IN):: teta(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, INTENT (IN):: ps(ip1jmp1)
    real, intent(out):: masse(ip1jmp1, llm)
    REAL, INTENT(IN):: pk(iip1, jjp1, llm)
    REAL, INTENT(IN):: pkf(:, :, :) ! (iim + 1, jjm + 1, llm)

    REAL, INTENT(IN):: phi(iim + 1, jjm + 1, llm)
    ! geopotential at mid-layer, in m2 s-2

    REAL, INTENT(out):: du(:, :, :) ! (iim + 1, jjm + 1, llm)
    real, INTENT(out):: dv(:, :, :) ! (iim + 1, jjm, llm)
    REAL, INTENT(out):: dteta(:, :, :) ! (iim + 1, jjm + 1, llm)
    real, INTENT(out):: dp(:, :) ! (iim + 1, jjm + 1)
    REAL, INTENT(out):: w(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(out):: pbaru(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(out):: pbarv(:, :, :) ! (iim + 1, jjm, llm)
    LOGICAL, INTENT(IN):: conser

    ! Local:
    REAL vcont(iim + 1, jjm, llm), ucont(iim + 1, jjm + 1, llm)
    REAL ang_3d(iim + 1, jjm + 1, llm), p(ip1jmp1, llmp1)
    REAL massebx(ip1jmp1, llm), masseby((iim + 1) * jjm, llm)
    REAL vorpot(iim + 1, jjm, llm)
    real ecin(iim + 1, jjm + 1, llm), convm(iim + 1, jjm + 1, llm)
    REAL bern(iim + 1, jjm + 1, llm)
    REAL massebxy(iim + 1, jjm, llm)
    INTEGER l
    real heure, time
    real ang, etot, ptot, ztot, stot, rmsdpdt, rmsv

    !-----------------------------------------------------------------------

    CALL covcont(ucov, vcov, ucont, vcont)
    forall (l = 1: llm + 1) p(:, l) = ap(l) + bp(l) * ps
    CALL massdair(p, masse)
    CALL massbar(masse, massebx, masseby)
    CALL massbarxy(masse, massebxy)
    CALL flumass(massebx, masseby, vcont, ucont, pbaru, pbarv)
    CALL dteta1(teta, pbaru, pbarv, dteta)
    CALL convmas(pbaru, pbarv, convm)
    dp = convm(:, :, 1) / airesurg
    w = vitvert(convm)
    CALL tourpot(vcov, ucov, massebxy, vorpot)
    CALL dudv1(vorpot, pbaru, pbarv, du(:, 2: jjm, :), dv)
    CALL enercin(vcov, ucov, vcont, ucont, ecin)
    bern = bernoui(phi, ecin)
    CALL dudv2(teta, pkf, bern, du(:, 2:jjm, :), dv)
    forall (l = 1: llm) ang_3d(:, :, l) = ucov(:, :, l) + constang_2d
    CALL advect(ang_3d, vcov, teta, w, massebx, masseby, du, dv, dteta)

    ! Problème de périodicité de dv sur les PC Linux. Problème
    ! d'arrondi probablement. Observé sur le code compilé avec pgf90
    ! 3.0-1.
    where (dv(1, :, :) /= dv(iim + 1, :, :)) dv(iim + 1, :, :) = dv(1, :, :)

    ! Sortie éventuelle des variables de contrôle :
    IF (conser) then
       CALL sortvarc(ucov, teta, ps, masse, pk, vorpot, phi, bern, dp, ang, &
            etot, ptot, ztot, stot, rmsdpdt, rmsv)
       time = real(itau) / day_step
       heure = mod(itau * dtvr / daysec, 1.) * 24.
       IF (abs(heure - 24.) <= 1e-4) heure = 0.
       PRINT 3500, itau, int(day_ini + time), heure, time
       PRINT 4000, ptot / ptot0, rmsdpdt, etot / etot0, ztot / ztot0, &
            stot / stot0, sqrt(rmsv / ptot), ang / ang0
    end IF

3500 FORMAT (4X, 'pas', I7, 5X, 'jour', i5, 1X, 'heure', F5.1, 4X, 'date', &
         F10.5)
4000 FORMAT (10X, 'masse', 4X, 'rmsdpdt', 7X, 'energie', 2X, 'enstrophie', &
         2X, 'entropie', 3X, 'rmsv', 4X, 'mt.ang', /, 'GLOB  ', F10.6, &
         E13.6, 5F10.3/)

  END SUBROUTINE caldyn

end module caldyn_m
