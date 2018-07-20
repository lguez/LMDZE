module calcul_fluxs_m

  implicit none

contains

  SUBROUTINE calcul_fluxs(dtime, tsurf, p1lay, cal, beta, coef1lay, ps, qsurf, &
       radsol, dif_grnd, t1lay, q1lay, u1lay, v1lay, petAcoef, peqAcoef, &
       petBcoef, peqBcoef, tsurf_new, evap, fluxlat, flux_t, dflux_s, dflux_l)

    ! Cette routine calcule les flux en h et q à l'interface et une
    ! température de surface.

    ! L. Fairhead, April 2000

    ! Note that, if cal = 0, beta = 1 and dif_grnd = 0, then tsurf_new
    ! = tsurf and qsurf = qsat.

    use nr_util, only: assert_eq

    USE fcttre, ONLY: foede, foeew
    USE suphec_m, ONLY: rcpd, rd, retv, rlstt, rlvtt, rtt
    USE yoethf_m, ONLY: r2es, r5ies, r5les, rvtmp2

    real, intent(IN):: dtime
    real, intent(IN):: tsurf(:) ! (knon) température de surface

    real, intent(IN):: p1lay(:) ! (knon)
    ! pression première couche (milieu de couche)

    real, intent(IN):: cal(:) ! (knon) capacité calorifique du sol
    real, intent(IN):: beta(:) ! (knon) évaporation réelle
    real, intent(IN):: coef1lay(:) ! (knon) coefficient d'échange
    real, intent(IN):: ps(:) ! (knon) pression au sol
    real, intent(OUT):: qsurf(:) ! (knon) humidité de l'air au-dessus du sol

    real, intent(IN):: radsol(:) ! (knon)
    ! rayonnement net au sol (longwave + shortwave)

    real, intent(IN):: dif_grnd(:) ! (knon)
    ! coefficient de diffusion vers le sol profond

    real, intent(IN):: t1lay(:), q1lay(:), u1lay(:), v1lay(:) ! (knon)

    real, intent(IN):: petAcoef(:), peqAcoef(:) ! (knon)
    ! coefficients A de la résolution de la couche limite pour T et q

    real, intent(IN):: petBcoef(:), peqBcoef(:) ! (knon)
    ! coefficients B de la résolution de la couche limite pour t et q

    real, intent(OUT):: tsurf_new(:) ! (knon) température au sol
    real, intent(OUT):: evap(:) ! (knon)

    real, intent(OUT):: fluxlat(:), flux_t(:) ! (knon)
    ! flux de chaleurs latente et sensible

    real, intent(OUT):: dflux_s(:), dflux_l(:) ! (knon)
    ! dérivées des flux de chaleurs sensible et latente par rapport à
    ! Ts (W m-2 K-1)

    ! Local:
    integer i
    integer knon ! nombre de points \`a traiter
    real, dimension(size(ps)):: mh, oh, mq, nq, oq, dq_s_dt, coef ! (knon)
    real qsat(size(ps)) ! (knon) mass fraction
    real sl(size(ps)) ! (knon) chaleur latente d'évaporation ou de sublimation
    logical delta
    real zcor
    real, parameter:: t_grnd = 271.35

    !---------------------------------------------------------------------

    knon = assert_eq([size(tsurf), size(p1lay), size(cal), size(beta), &
         size(coef1lay), size(ps), size(qsurf), size(radsol), size(dif_grnd), &
         size(t1lay), size(q1lay), size(u1lay), size(v1lay), size(petAcoef), &
         size(peqAcoef), size(petBcoef), size(peqBcoef), size(tsurf_new), &
         size(evap), size(fluxlat), size(flux_t), size(dflux_s), &
         size(dflux_l)], "calcul_fluxs knon")

    ! Traitement de l'humidité du sol

    DO i = 1, knon
       delta = rtt >= tsurf(i)
       qsat(i) = MIN(0.5, r2es * FOEEW(tsurf(i), delta) / ps(i))
       zcor = 1. / (1. - retv * qsat(i))
       qsat(i) = qsat(i) * zcor
       dq_s_dt(i) = RCPD * FOEDE(tsurf(i), delta, merge(R5IES * RLSTT, &
            R5LES * RLVTT, delta) / RCPD / (1. + RVTMP2 * q1lay(i)), &
            qsat(i), zcor) / RLVTT
    ENDDO

    coef = coef1lay * (1. + SQRT(u1lay**2 + v1lay**2)) * p1lay / (RD * t1lay)
    sl = merge(RLSTT, RLVTT, tsurf < RTT)

    ! Q
    oq = 1. - beta * coef * peqBcoef * dtime
    mq = beta * coef * (peqAcoef - qsat + dq_s_dt * tsurf) / oq
    nq = - beta * coef * dq_s_dt / oq

    ! H
    oh = 1. - coef * petBcoef * dtime
    mh = coef * petAcoef / oh
    dflux_s = - coef * RCPD / oh

    tsurf_new = (tsurf + cal / RCPD * dtime * (radsol + mh + sl * mq) &
         + dif_grnd * t_grnd * dtime) / (1. - dtime * cal / RCPD * (dflux_s &
         + sl * nq) + dtime * dif_grnd)
    evap = - mq - nq * tsurf_new
    fluxlat = - evap * sl
    flux_t = mh + dflux_s * tsurf_new
    dflux_l = sl * nq
    qsurf = (peqAcoef - peqBcoef * evap * dtime) * (1. - beta) + beta * (qsat &
         + dq_s_dt * (tsurf_new - tsurf))

  END SUBROUTINE calcul_fluxs

end module calcul_fluxs_m
