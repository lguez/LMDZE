module calcul_fluxs_m

  implicit none

contains

  SUBROUTINE calcul_fluxs(nisurf, dtime, tsurf, p1lay, cal, beta, coef1lay, &
       ps, qsurf, radsol, dif_grnd, t1lay, q1lay, u1lay, v1lay, petAcoef, &
       peqAcoef, petBcoef, peqBcoef, tsurf_new, evap, fluxlat, fluxsens, &
       dflux_s, dflux_l)

    ! Cette routine calcule les fluxs en h et q à l'interface et une
    ! température de surface.

    ! L. Fairhead April 2000

    USE abort_gcm_m, ONLY: abort_gcm
    USE indicesol, ONLY: is_ter
    USE fcttre, ONLY: dqsatl, dqsats, foede, foeew, qsatl, qsats, thermcep
    USE interface_surf, ONLY: run_off
    use nr_util, only: assert_eq
    USE suphec_m, ONLY: rcpd, rd, retv, rkappa, rlstt, rlvtt, rtt
    USE yoethf_m, ONLY: r2es, r5ies, r5les, rvtmp2

    integer, intent(IN):: nisurf ! surface a traiter
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
    ! petBcoef coeff. B de la resolution de la CL pour t
    ! peqBcoef coeff. B de la resolution de la CL pour q

    real, intent(OUT):: tsurf_new(:) ! (knon) température au sol
    real, intent(OUT):: evap(:), fluxlat(:), fluxsens(:) ! (knon)
    ! fluxlat flux de chaleur latente
    ! fluxsens flux de chaleur sensible
    real, intent(OUT):: dflux_s(:), dflux_l(:) ! (knon)
    ! Dérivées des flux dF/dTs (W m-2 K-1)
    ! dflux_s derivee du flux de chaleur sensible / Ts
    ! dflux_l derivee du flux de chaleur latente / Ts

    ! Local:
    integer i
    real, dimension(size(ps)) :: zx_mh, zx_nh, zx_oh
    real, dimension(size(ps)) :: zx_mq, zx_nq, zx_oq
    real, dimension(size(ps)) :: zx_pkh, zx_dq_s_dt, zx_qsat, zx_coef
    real, dimension(size(ps)) :: zx_sl, zx_k1
    real, dimension(size(ps)) :: zx_q_0 , d_ts
    logical zdelta
    real zcvm5, zx_qs, zcor, zx_dq_s_dh
    real :: bilan_f, fq_fonte
    REAL :: subli, fsno
    REAL :: qsat_new, q1_new
    integer knon ! nombre de points a traiter
    real, parameter:: t_grnd = 271.35, t_coup = 273.15

    !---------------------------------------------------------------------

    knon = assert_eq((/size(tsurf), size(p1lay), size(cal), size(beta), &
         size(coef1lay), size(ps), size(qsurf), size(radsol), size(dif_grnd), &
         size(t1lay), size(q1lay), size(u1lay), size(v1lay), size(petAcoef), &
         size(peqAcoef), size(petBcoef), size(peqBcoef), size(tsurf_new), &
         size(evap), size(fluxlat), size(fluxsens), size(dflux_s), &
         size(dflux_l)/), "calcul_fluxs knon")

    if (size(run_off) /= knon .AND. nisurf == is_ter) then
       print *, 'Bizarre, le nombre de points continentaux'
       print *, 'a change entre deux appels. J''arrete.'
       call abort_gcm('calcul_fluxs', 'Pb run_off', 1)
    endif

    ! Traitement humidite du sol

    evap = 0.
    fluxsens=0.
    fluxlat=0.
    dflux_s = 0.
    dflux_l = 0.

    ! zx_qs = qsat en kg/kg

    DO i = 1, knon
       zx_pkh(i) = (ps(i)/ps(i))**RKAPPA
       IF (thermcep) THEN
          zdelta= rtt >= tsurf(i)
          zcvm5 = merge(R5IES*RLSTT, R5LES*RLVTT, zdelta)
          zcvm5 = zcvm5 / RCPD / (1.0+RVTMP2*q1lay(i))
          zx_qs= r2es * FOEEW(tsurf(i), zdelta)/ps(i)
          zx_qs=MIN(0.5, zx_qs)
          zcor=1./(1.-retv*zx_qs)
          zx_qs=zx_qs*zcor
          zx_dq_s_dh = FOEDE(tsurf(i), zdelta, zcvm5, zx_qs, zcor) &
               /RLVTT / zx_pkh(i)
       ELSE
          IF (tsurf(i).LT.t_coup) THEN
             zx_qs = qsats(tsurf(i)) / ps(i)
             zx_dq_s_dh = dqsats(tsurf(i), zx_qs)/RLVTT &
                  / zx_pkh(i)
          ELSE
             zx_qs = qsatl(tsurf(i)) / ps(i)
             zx_dq_s_dh = dqsatl(tsurf(i), zx_qs)/RLVTT &
                  / zx_pkh(i)
          ENDIF
       ENDIF
       zx_dq_s_dt(i) = RCPD * zx_pkh(i) * zx_dq_s_dh
       zx_qsat(i) = zx_qs
       zx_coef(i) = coef1lay(i) &
            * (1.0+SQRT(u1lay(i)**2+v1lay(i)**2)) &
            * p1lay(i)/(RD*t1lay(i))

    ENDDO

    ! === Calcul de la temperature de surface ===

    ! zx_sl = chaleur latente d'evaporation ou de sublimation

    do i = 1, knon
       zx_sl(i) = RLVTT
       if (tsurf(i) .LT. RTT) zx_sl(i) = RLSTT
       zx_k1(i) = zx_coef(i)
    enddo

    do i = 1, knon
       ! Q
       zx_oq(i) = 1. - (beta(i) * zx_k1(i) * peqBcoef(i) * dtime)
       zx_mq(i) = beta(i) * zx_k1(i) * &
            (peqAcoef(i) - zx_qsat(i) &
            + zx_dq_s_dt(i) * tsurf(i)) &
            / zx_oq(i)
       zx_nq(i) = beta(i) * zx_k1(i) * (-1. * zx_dq_s_dt(i)) &
            / zx_oq(i)

       ! H
       zx_oh(i) = 1. - (zx_k1(i) * petBcoef(i) * dtime)
       zx_mh(i) = zx_k1(i) * petAcoef(i) / zx_oh(i)
       zx_nh(i) = - (zx_k1(i) * RCPD * zx_pkh(i))/ zx_oh(i)

       ! Tsurface
       tsurf_new(i) = (tsurf(i) + cal(i)/(RCPD * zx_pkh(i)) * dtime * &
            (radsol(i) + zx_mh(i) + zx_sl(i) * zx_mq(i)) &
            + dif_grnd(i) * t_grnd * dtime)/ &
            ( 1. - dtime * cal(i)/(RCPD * zx_pkh(i)) * ( &
            zx_nh(i) + zx_sl(i) * zx_nq(i)) &
            + dtime * dif_grnd(i))


       ! Y'a-t-il fonte de neige?

       ! fonte_neige = (nisurf /= is_oce) .AND. &
       ! & (snow(i) > epsfra .OR. nisurf == is_sic .OR. nisurf == is_lic) &
       ! & .AND. (tsurf_new(i) >= RTT)
       ! if (fonte_neige) tsurf_new(i) = RTT
       d_ts(i) = tsurf_new(i) - tsurf(i)
       ! zx_h_ts(i) = tsurf_new(i) * RCPD * zx_pkh(i)
       ! zx_q_0(i) = zx_qsat(i) + zx_dq_s_dt(i) * d_ts(i)
       !== flux_q est le flux de vapeur d'eau: kg/(m**2 s) positive vers bas
       !== flux_t est le flux de cpt (energie sensible): j/(m**2 s)
       evap(i) = - zx_mq(i) - zx_nq(i) * tsurf_new(i)
       fluxlat(i) = - evap(i) * zx_sl(i)
       fluxsens(i) = zx_mh(i) + zx_nh(i) * tsurf_new(i)
       ! Derives des flux dF/dTs (W m-2 K-1):
       dflux_s(i) = zx_nh(i)
       dflux_l(i) = (zx_sl(i) * zx_nq(i))
       ! Nouvelle valeure de l'humidite au dessus du sol
       qsat_new=zx_qsat(i) + zx_dq_s_dt(i) * d_ts(i)
       q1_new = peqAcoef(i) - peqBcoef(i)*evap(i)*dtime
       qsurf(i)=q1_new*(1.-beta(i)) + beta(i)*qsat_new
    ENDDO

  END SUBROUTINE calcul_fluxs

end module calcul_fluxs_m
