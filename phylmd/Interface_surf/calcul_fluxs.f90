module calcul_fluxs_m

  implicit none

contains

  SUBROUTINE calcul_fluxs(qsurf, tsurf_new, flux_q, fluxlat, flux_t, dflux_s, &
       dflux_l, tsurf, p1lay, cdragh, ps, radsol, t1lay, q1lay, u1lay, v1lay, &
       ah, aq, bh, bq, soilflux, cal, beta, dif_grnd)

    ! Cette routine calcule les flux en h et q à l'interface et une
    ! température de surface.

    ! L. Fairhead, April 2000

    ! Note that if cal = 0, beta = 1 and dif_grnd = 0, then tsurf_new
    ! = tsurf and qsurf = qsat.

    ! Libraries:
    use nr_util, only: assert_eq

    use comconst, only: dtphys
    USE fcttre, ONLY: foede, foeew
    USE suphec_m, ONLY: rcpd, rd, retv, rlstt, rlvtt, rtt
    USE yoethf_m, ONLY: r2es, r5ies, r5les, rvtmp2

    real, intent(OUT):: qsurf(:) ! (knon)
    ! humidité de l'air au-dessus de la surface
    
    real, intent(OUT):: tsurf_new(:) ! (knon) new surface temperature, in K

    real, intent(OUT):: flux_q(:) ! (knon)
    ! downward water vapor flux at the surface, in kg m-2 s-1

    real, intent(OUT):: fluxlat(:), flux_t(:) ! (knon)
    ! downward flux of latent and sensible heat (c_p T) at the surface, in W m-2

    real, intent(OUT):: dflux_s(:), dflux_l(:) ! (knon)
    ! dérivées des flux de chaleurs sensible et latente par rapport à
    ! T_surf (W m-2 K-1)

    real, intent(IN):: tsurf(:) ! (knon)
    ! surface temperature at previous time step (K)

    real, intent(IN):: p1lay(:) ! (knon)
    ! pression première couche (milieu de couche)

    real, intent(IN):: cdragh(:) ! (knon)
    ! coefficient d'\'echange, sans dimension
    
    real, intent(IN):: ps(:) ! (knon) pression au sol, en Pa

    real, intent(IN):: radsol(:) ! (knon)
    ! net downward radiative (longwave + shortwave) flux at the
    ! surface, in W / m2

    real, intent(IN):: t1lay(:) ! (knon) temp\'erature de l'air 1\`ere couche
    real, intent(IN):: q1lay(:), u1lay(:), v1lay(:) ! (knon)

    real, intent(IN):: ah(:), aq(:) ! (knon)
    ! coefficients A de la résolution de la couche limite pour T et q

    real, intent(IN):: bh(:), bq(:) ! (knon)
    ! coefficients B de la résolution de la couche limite pour t et q

    REAL, intent(IN):: soilflux(:) ! (knon)
    real, intent(IN):: cal(:) ! (knon) 1 / soilcap, in J-1 K m2
    real, intent(IN):: beta(:) ! (knon) évaporation réelle
    real, intent(IN):: dif_grnd ! coefficient de diffusion vers le sol profond

    ! Local:
    integer i
    integer knon ! nombre de points \`a traiter
    real, dimension(size(ps)):: mh, oh, mq, nq, oq, dq_s_dt, coef ! (knon)
    real qsat(size(ps)) ! (knon) mass fraction
    real sl(size(ps)) ! (knon) chaleur latente d'évaporation ou de sublimation
    logical ice
    real zcor
    real, parameter:: t_grnd = 271.35 ! in K
    real, parameter:: min_wind_speed = 1. ! in m s-1

    !---------------------------------------------------------------------

    knon = assert_eq([size(tsurf), size(p1lay), size(cal), size(beta), &
         size(cdragh), size(ps), size(qsurf), size(radsol), size(t1lay), &
         size(q1lay), size(u1lay), size(v1lay), size(ah), size(aq), &
         size(bh), size(bq), size(tsurf_new), size(flux_q), &
         size(fluxlat), size(flux_t), size(dflux_s), size(dflux_l)], &
         "calcul_fluxs knon")

    ! Traitement de l'humidité du sol

    DO i = 1, knon
       ice = tsurf(i) <= rtt
       qsat(i) = MIN(0.5, r2es * FOEEW(tsurf(i), ice) / ps(i))
       zcor = 1. / (1. - retv * qsat(i))
       qsat(i) = qsat(i) * zcor
       dq_s_dt(i) = RCPD * FOEDE(tsurf(i), ice, merge(R5IES * RLSTT, &
            R5LES * RLVTT, ice) / RCPD / (1. + RVTMP2 * q1lay(i)), qsat(i), &
            zcor) / RLVTT
    ENDDO

    coef = cdragh * (min_wind_speed + SQRT(u1lay**2 + v1lay**2)) * p1lay &
         / (RD * t1lay)
    sl = merge(RLSTT, RLVTT, tsurf < RTT)

    ! Q
    oq = 1. - beta * coef * bq * dtphys
    mq = beta * coef * (aq - qsat + dq_s_dt * tsurf) / oq
    nq = - beta * coef * dq_s_dt / oq

    ! H
    oh = 1. - coef * bh * dtphys
    mh = coef * ah / oh
    dflux_s = - coef * RCPD / oh

    tsurf_new = (tsurf + cal * dtphys * (radsol + soilflux + mh + sl * mq) &
         + dif_grnd * t_grnd * dtphys) / (1. - dtphys * cal * (dflux_s + sl &
         * nq) + dtphys * dif_grnd)
    flux_q = mq + nq * tsurf_new
    fluxlat = flux_q * sl
    flux_t = mh + dflux_s * tsurf_new
    dflux_l = sl * nq
    qsurf = (aq + bq * flux_q * dtphys) * (1. - beta) + beta * (qsat + dq_s_dt &
         * (tsurf_new - tsurf))

  END SUBROUTINE calcul_fluxs

end module calcul_fluxs_m
