module calcul_fluxs_m

  implicit none

contains

  SUBROUTINE calcul_fluxs(dtime, tsurf, p1lay, cal, beta, coef1lay, ps, &
       qsurf, radsol, dif_grnd, t1lay, q1lay, u1lay, v1lay, petAcoef, &
       peqAcoef, petBcoef, peqBcoef, tsurf_new, evap, fluxlat, fluxsens, &
       dflux_s, dflux_l)

    ! Cette routine calcule les fluxs en h et q à l'interface et une
    ! température de surface.

    ! L. Fairhead April 2000

    USE abort_gcm_m, ONLY: abort_gcm
    USE fcttre, ONLY: dqsatl, dqsats, foede, foeew, qsatl, qsats, thermcep
    USE indicesol, ONLY: is_ter
    use nr_util, only: assert_eq
    USE suphec_m, ONLY: rcpd, rd, retv, rkappa, rlstt, rlvtt, rtt
    USE yoethf_m, ONLY: r2es, r5ies, r5les, rvtmp2

    real, intent(IN):: dtime
    real, intent(IN):: tsurf(:) ! (knon) temperature de surface
    real, intent(IN):: p1lay(:) ! (knon) pression 1er niveau (milieu de couche)
    real, intent(IN):: cal(:) ! (knon) capacité calorifique du sol
    real, intent(IN):: beta(:) ! (knon) evap reelle
    real, intent(IN):: coef1lay(:) ! (knon) coefficient d'échange
    real, intent(IN):: ps(:) ! (knon) pression au sol
    real, intent(OUT):: qsurf(:) ! (knon) humidite de l'air au dessus du sol
    real, intent(IN):: radsol(:) ! (knon) rayonnement net au sol (LW + SW)

    real, intent(IN):: dif_grnd(:) ! (knon)
    ! coefficient diffusion vers le sol profond

    real, intent(IN):: t1lay(:), q1lay(:), u1lay(:), v1lay(:) ! (knon)

    real, intent(IN):: petAcoef(:), peqAcoef(:) ! (knon)
    ! coefficients A de la résolution de la couche limite pour t et q

    real, intent(IN):: petBcoef(:), peqBcoef(:) ! (knon)
    ! coeff. B de la resolution de la CL pour t et q

    real, intent(OUT):: tsurf_new(:) ! (knon) température au sol
    real, intent(OUT):: evap(:) ! (knon)

    real, intent(OUT):: fluxlat(:), fluxsens(:) ! (knon)
    ! flux de chaleur latente et sensible

    real, intent(OUT):: dflux_s(:), dflux_l(:) ! (knon)
    ! dérivées des flux de chaleurs sensible et latente par rapport à
    ! Ts (W m-2 K-1)

    ! Local:
    integer i
    integer knon ! nombre de points a traiter
    real, dimension(size(ps)):: mh, oh, mq, nq, oq
    real, dimension(size(ps)):: dq_s_dt, coef
    real qsat(size(ps)) ! qsat en kg/kg
    real sl(size(ps)) ! chaleur latente d'evaporation ou de sublimation
    logical delta
    real zcor
    real, parameter:: t_grnd = 271.35, t_coup = 273.15

    !---------------------------------------------------------------------

    knon = assert_eq((/size(tsurf), size(p1lay), size(cal), size(beta), &
         size(coef1lay), size(ps), size(qsurf), size(radsol), size(dif_grnd), &
         size(t1lay), size(q1lay), size(u1lay), size(v1lay), size(petAcoef), &
         size(peqAcoef), size(petBcoef), size(peqBcoef), size(tsurf_new), &
         size(evap), size(fluxlat), size(fluxsens), size(dflux_s), &
         size(dflux_l)/), "calcul_fluxs knon")

    ! Traitement humidite du sol

    IF (thermcep) THEN
       DO i = 1, knon
          delta = rtt >= tsurf(i)
          qsat(i) = MIN(0.5, r2es * FOEEW(tsurf(i), delta) / ps(i))
          zcor = 1. / (1. - retv * qsat(i))
          qsat(i) = qsat(i) * zcor
          dq_s_dt(i) = RCPD * FOEDE(tsurf(i), delta, merge(R5IES * RLSTT, &
               R5LES * RLVTT, delta) / RCPD / (1. + RVTMP2 * q1lay(i)), &
               qsat(i), zcor) / RLVTT
       ENDDO
    ELSE
       DO i = 1, knon
          IF (tsurf(i) < t_coup) THEN
             qsat(i) = qsats(tsurf(i)) / ps(i)
             dq_s_dt(i) = RCPD * dqsats(tsurf(i), qsat(i)) / RLVTT
          ELSE
             qsat(i) = qsatl(tsurf(i)) / ps(i)
             dq_s_dt(i) = RCPD * dqsatl(tsurf(i), qsat(i)) / RLVTT
          ENDIF
       ENDDO
    ENDIF

    coef = coef1lay * (1. + SQRT(u1lay**2 + v1lay**2)) * p1lay / (RD * t1lay)
    sl = merge(RLSTT, RLVTT, tsurf < RTT)

    ! Q
    oq = 1. - (beta * coef * peqBcoef * dtime)
    mq = beta * coef * (peqAcoef - qsat + dq_s_dt * tsurf) / oq
    nq = beta * coef * (- 1. * dq_s_dt) / oq

    ! H
    oh = 1. - (coef * petBcoef * dtime)
    mh = coef * petAcoef / oh
    dflux_s = - (coef * RCPD)/ oh

    ! Tsurface
    tsurf_new = (tsurf + cal / RCPD * dtime * (radsol + mh + sl * mq) &
         + dif_grnd * t_grnd * dtime) / (1. - dtime * cal / RCPD * (dflux_s &
         + sl * nq) + dtime * dif_grnd)

    evap = - mq - nq * tsurf_new
    fluxlat = - evap * sl
    fluxsens = mh + dflux_s * tsurf_new
    dflux_l = sl * nq

    ! Nouvelle valeur de l'humidité au dessus du sol :
    qsurf = (peqAcoef - peqBcoef * evap * dtime) * (1. - beta) + beta * (qsat &
         + dq_s_dt * (tsurf_new - tsurf))

  END SUBROUTINE calcul_fluxs

end module calcul_fluxs_m
