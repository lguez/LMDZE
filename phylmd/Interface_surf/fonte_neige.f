module fonte_neige_m

  implicit none

contains

  SUBROUTINE fonte_neige( klon, knon, nisurf, dtime,  &
       tsurf, p1lay, cal, beta, coef1lay, ps,  &
       precip_rain, precip_snow, snow, qsol,  &
       radsol, dif_grnd, t1lay, q1lay, u1lay, v1lay,  &
       petAcoef, peqAcoef, petBcoef, peqBcoef,  &
       tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l,  &
       fqcalving, ffonte, run_off_lic_0)

    ! Routine de traitement de la fonte de la neige dans le cas du traitement
    ! de sol simplifié

    ! LF 03/2001
    ! input:
    ! knon nombre de points a traiter
    ! nisurf surface a traiter
    ! tsurf temperature de surface
    ! p1lay pression 1er niveau (milieu de couche)
    ! cal capacite calorifique du sol
    ! beta evap reelle
    ! coef1lay coefficient d'echange
    ! ps pression au sol
    ! precip_rain precipitations liquides
    ! precip_snow precipitations solides
    ! snow champs hauteur de neige
    ! qsol hauteur d'eau contenu dans le sol
    ! runoff runoff en cas de trop plein
    ! petAcoef coeff. A de la resolution de la CL pour t
    ! peqAcoef coeff. A de la resolution de la CL pour q
    ! petBcoef coeff. B de la resolution de la CL pour t
    ! peqBcoef coeff. B de la resolution de la CL pour q
    ! radsol rayonnement net aus sol (LW + SW)
    ! dif_grnd coeff. diffusion vers le sol profond

    ! output:
    ! tsurf_new temperature au sol
    ! fluxsens flux de chaleur sensible
    ! fluxlat flux de chaleur latente
    ! dflux_s derivee du flux de chaleur sensible / Ts
    ! dflux_l derivee du flux de chaleur latente / Ts
    ! in/out:
    ! run_off_lic_0 run off glacier du pas de temps précedent


    use indicesol
    use SUPHEC_M
    use yoethf_m
    use fcttre
    use interface_surf
    !IM cf JLD

    ! Parametres d'entree
    integer, intent(IN) :: knon, nisurf, klon
    real , intent(IN) :: dtime
    real, dimension(klon), intent(IN) :: petAcoef, peqAcoef
    real, dimension(klon), intent(IN) :: petBcoef, peqBcoef
    real, dimension(klon), intent(IN) :: ps, q1lay
    real, dimension(klon), intent(IN) :: tsurf, p1lay, cal, beta, coef1lay
    real, dimension(klon), intent(IN) :: precip_rain, precip_snow
    real, dimension(klon), intent(IN) :: radsol, dif_grnd
    real, dimension(klon), intent(IN) :: t1lay, u1lay, v1lay
    real, dimension(klon), intent(INOUT) :: snow, qsol

    ! Parametres sorties
    real, dimension(klon), intent(INOUT):: tsurf_new, evap, fluxsens, fluxlat
    real, dimension(klon), intent(INOUT):: dflux_s, dflux_l
    ! Flux thermique utiliser pour fondre la neige
    real, dimension(klon), intent(INOUT):: ffonte
    ! Flux d'eau "perdue" par la surface et necessaire pour que limiter la
    ! hauteur de neige, en kg/m2/s
    real, dimension(klon), intent(INOUT):: fqcalving
    real, dimension(klon), intent(INOUT):: run_off_lic_0
    ! Variables locales
    ! Masse maximum de neige (kg/m2). Au dessus de ce seuil, la neige
    ! en exces "s'ecoule" (calving)
    ! real, parameter :: snow_max=1.
    !IM cf JLD/GK
    real, parameter :: snow_max=3000.
    integer :: i
    real, dimension(klon) :: zx_mh, zx_nh, zx_oh
    real, dimension(klon) :: zx_mq, zx_nq, zx_oq
    real, dimension(klon) :: zx_pkh, zx_dq_s_dt, zx_qsat, zx_coef
    real, dimension(klon) :: zx_sl, zx_k1
    real, dimension(klon) :: zx_q_0 , d_ts
    real :: zdelta, zcvm5, zx_qs, zcor, zx_dq_s_dh
    real :: bilan_f, fq_fonte
    REAL :: subli, fsno
    REAL, DIMENSION(klon) :: bil_eau_s, snow_evap
    real, parameter :: t_grnd = 271.35, t_coup = 273.15
    !! PB temporaire en attendant mieux pour le modele de neige
    ! REAL, parameter :: chasno = RLMLT/(2.3867E+06*0.15)
    REAL, parameter :: chasno = 3.334E+05/(2.3867E+06*0.15)
    !IM cf JLD/ GKtest
    REAL, parameter :: chaice = 3.334E+05/(2.3867E+06*0.15)
    ! fin GKtest

    logical, save :: check = .FALSE.
    character (len = 20) :: modname = 'fonte_neige'
    logical, save :: neige_fond = .false.
    real, save :: max_eau_sol = 150.0
    character (len = 80) :: abort_message
    logical, save :: first = .true., second=.false.
    real :: coeff_rel

    if (check) write(*, *)'Entree ', modname, ' surface = ', nisurf

    ! Initialisations
    coeff_rel = dtime/(tau_calv * rday)
    bil_eau_s = 0.
    DO i = 1, knon
       zx_pkh(i) = (ps(i)/ps(i))**RKAPPA
       IF (thermcep) THEN
          zdelta=MAX(0., SIGN(1., rtt-tsurf(i)))
          zcvm5 = R5LES*RLVTT*(1.-zdelta) + R5IES*RLSTT*zdelta
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
    enddo

    WHERE (precip_snow > 0.) snow = snow + (precip_snow * dtime)
    snow_evap = 0.
    WHERE (evap > 0. )
       snow_evap = MIN (snow / dtime, evap)
       snow = snow - snow_evap * dtime
       snow = MAX(0.0, snow)
    end where

    ! bil_eau_s = bil_eau_s + (precip_rain * dtime) - (evap - snow_evap) * dtime
    bil_eau_s = (precip_rain * dtime) - (evap - snow_evap) * dtime


    ! Y'a-t-il fonte de neige?

    ffonte=0.
    do i = 1, knon
       neige_fond = ((snow(i) > epsfra .OR. nisurf == is_sic .OR. nisurf == is_lic) &
            .AND. tsurf_new(i) >= RTT)
       if (neige_fond) then
          fq_fonte = MIN( MAX((tsurf_new(i)-RTT )/chasno, 0.0), snow(i))
          ffonte(i) = fq_fonte * RLMLT/dtime
          snow(i) = max(0., snow(i) - fq_fonte)
          bil_eau_s(i) = bil_eau_s(i) + fq_fonte
          tsurf_new(i) = tsurf_new(i) - fq_fonte * chasno
          !IM cf JLD OK
          !IM cf JLD/ GKtest fonte aussi pour la glace
          IF (nisurf == is_sic .OR. nisurf == is_lic ) THEN
             fq_fonte = MAX((tsurf_new(i)-RTT )/chaice, 0.0)
             ffonte(i) = ffonte(i) + fq_fonte * RLMLT/dtime
             bil_eau_s(i) = bil_eau_s(i) + fq_fonte
             tsurf_new(i) = RTT
          ENDIF
          d_ts(i) = tsurf_new(i) - tsurf(i)
       endif

       ! s'il y a une hauteur trop importante de neige, elle s'coule
       fqcalving(i) = max(0., snow(i) - snow_max)/dtime
       snow(i)=min(snow(i), snow_max)

       IF (nisurf == is_ter) then
          qsol(i) = qsol(i) + bil_eau_s(i)
          run_off(i) = run_off(i) + MAX(qsol(i) - max_eau_sol, 0.0)
          qsol(i) = MIN(qsol(i), max_eau_sol)
       else if (nisurf == is_lic) then
          run_off_lic(i) = (coeff_rel * fqcalving(i)) + &
               (1. - coeff_rel) * run_off_lic_0(i)
          run_off_lic_0(i) = run_off_lic(i)
          run_off_lic(i) = run_off_lic(i) + bil_eau_s(i)/dtime
       endif
    enddo

  END SUBROUTINE fonte_neige

end module fonte_neige_m
