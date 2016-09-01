module fonte_neige_m

  implicit none

contains

  SUBROUTINE fonte_neige(nisurf, dtime, tsurf, p1lay, beta, coef1lay, ps, &
       precip_rain, precip_snow, snow, qsol, t1lay, q1lay, u1lay, v1lay, &
       petAcoef, peqAcoef, petBcoef, peqBcoef, tsurf_new, evap, fqcalving, &
       ffonte, run_off_lic_0)

    ! Routine de traitement de la fonte de la neige dans le cas du traitement
    ! de sol simplifi\'e

    ! Laurent Fairhead, March, 2001

    USE fcttre, ONLY: foeew, qsatl, qsats
    USE indicesol, ONLY: epsfra, is_lic, is_sic, is_ter
    USE interface_surf, ONLY: run_off_lic, tau_calv
    use nr_util, only: assert_eq
    USE suphec_m, ONLY: rcpd, rday, retv, rlmlt, rlstt, rlvtt, rtt
    USE yoethf_m, ONLY: r2es, r5ies, r5les, rvtmp2

    integer, intent(IN):: nisurf ! surface \`a traiter
    real, intent(IN):: dtime ! pas de temps de la physique (en s)
    real, intent(IN):: tsurf(:) ! (knon) temperature de surface
    real, intent(IN):: p1lay(:) ! (knon) pression 1er niveau (milieu de couche)
    real, intent(IN):: beta(:) ! (knon) evap reelle
    real, intent(IN):: coef1lay(:) ! (knon) coefficient d'echange
    real, intent(IN):: ps(:) ! (knon) pression au sol

    real, intent(IN):: precip_rain(:) ! (knon)
    ! precipitation, liquid water mass flux (kg / m2 / s), positive down

    real, intent(IN):: precip_snow(:) ! (knon)
    ! precipitation, solid water mass flux (kg / m2 / s), positive down

    real, intent(INOUT):: snow(:) ! (knon)
    ! column-density of mass of snow, in kg m-2

    real, intent(INOUT):: qsol(:) ! (knon)
    ! column-density of water in soil, in kg m-2

    real, intent(IN):: t1lay(:) ! (knon)
    real, intent(IN):: q1lay(:) ! (knon)
    real, intent(IN):: u1lay(:), v1lay(:) ! (knon)

    real, intent(IN):: petAcoef(:), peqAcoef(:) ! (knon)
    ! coefficients A de la r\'esolution de la couche limite pour t et q

    real, intent(IN):: petBcoef(:), peqBcoef(:) ! (knon)
    ! coefficients B de la r\'esolution de la couche limite pour t et q

    real, intent(INOUT):: tsurf_new(:) ! (knon) temp\'erature au sol
    real, intent(IN):: evap(:) ! (knon)

    real, intent(OUT):: fqcalving(:) ! (knon)
    ! flux d'eau "perdue" par la surface et n\'ecessaire pour limiter la
    ! hauteur de neige, en kg / m2 / s

    real, intent(OUT):: ffonte(:) ! (knon)
    ! flux thermique utilis\'é pour fondre la neige

    real, intent(INOUT):: run_off_lic_0(:) ! (knon)
    ! run off glacier du pas de temps pr\'ecedent

    ! Local:

    integer knon ! nombre de points \`a traiter
    real, parameter:: snow_max=3000.
    ! Masse maximum de neige (kg / m2). Au dessus de ce seuil, la neige
    ! en exces "s'ecoule" (calving)

    integer i
    logical zdelta
    real zcvm5, zx_qs, zcor
    real fq_fonte
    REAL bil_eau_s(size(ps)) ! in kg m-2
    real snow_evap(size(ps)) ! in kg m-2 s-1
    REAL, parameter:: chasno = 3.334E5 / (2.3867E6*0.15)
    REAL, parameter:: chaice = 3.334E5 / (2.3867E6*0.15)
    real, parameter:: max_eau_sol = 150. ! in kg m-2
    real coeff_rel

    !--------------------------------------------------------------------

    knon = assert_eq((/size(tsurf), size(p1lay), size(beta), size(coef1lay), &
         size(ps), size(precip_rain), size(precip_snow), size(snow), &
         size(qsol), size(t1lay), size(q1lay), size(u1lay), size(v1lay), &
         size(petAcoef), size(peqAcoef), size(petBcoef), size(peqBcoef), &
         size(tsurf_new), size(evap), size(fqcalving), size(ffonte), &
         size(run_off_lic_0)/), "fonte_neige knon")

    ! Initialisations
    coeff_rel = dtime / (tau_calv * rday)
    bil_eau_s = 0.
    DO i = 1, knon
       zdelta= rtt >= tsurf(i)
       zcvm5 = merge(R5IES*RLSTT, R5LES*RLVTT, zdelta)
       zcvm5 = zcvm5 / RCPD / (1. + RVTMP2*q1lay(i))
       zx_qs= r2es * FOEEW(tsurf(i), zdelta) / ps(i)
       zx_qs=MIN(0.5, zx_qs)
       zcor=1. / (1. - retv*zx_qs)
       zx_qs=zx_qs*zcor
    ENDDO

    ! Calcul de la temperature de surface

    WHERE (precip_snow > 0.) snow = snow + precip_snow * dtime

    WHERE (evap > 0.)
       snow_evap = MIN(snow / dtime, evap)
       snow = snow - snow_evap * dtime
       snow = MAX(0., snow)
    elsewhere
       snow_evap = 0.
    end where

    bil_eau_s = precip_rain * dtime - (evap(:knon) - snow_evap(:knon)) * dtime

    ! Y a-t-il fonte de neige ?

    ffonte=0.
    do i = 1, knon
       if ((snow(i) > epsfra .OR. nisurf == is_sic &
            .OR. nisurf == is_lic) .AND. tsurf_new(i) >= RTT) then
          fq_fonte = MIN(MAX((tsurf_new(i) - RTT) / chasno, 0.), snow(i))
          ffonte(i) = fq_fonte * RLMLT / dtime
          snow(i) = max(0., snow(i) - fq_fonte)
          bil_eau_s(i) = bil_eau_s(i) + fq_fonte
          tsurf_new(i) = tsurf_new(i) - fq_fonte * chasno
          !IM cf JLD/ GKtest fonte aussi pour la glace
          IF (nisurf == is_sic .OR. nisurf == is_lic) THEN
             fq_fonte = MAX((tsurf_new(i) - RTT) / chaice, 0.)
             ffonte(i) = ffonte(i) + fq_fonte * RLMLT / dtime
             bil_eau_s(i) = bil_eau_s(i) + fq_fonte
             tsurf_new(i) = RTT
          ENDIF
       endif

       ! S'il y a une hauteur trop importante de neige, elle s'\'ecoule
       fqcalving(i) = max(0., snow(i) - snow_max) / dtime
       snow(i)=min(snow(i), snow_max)

       IF (nisurf == is_ter) then
          qsol(i) = qsol(i) + bil_eau_s(i)
          qsol(i) = MIN(qsol(i), max_eau_sol)
       else if (nisurf == is_lic) then
          run_off_lic(i) = (coeff_rel * fqcalving(i)) + &
               (1. - coeff_rel) * run_off_lic_0(i)
          run_off_lic_0(i) = run_off_lic(i)
          run_off_lic(i) = run_off_lic(i) + bil_eau_s(i) / dtime
       endif
    enddo

  END SUBROUTINE fonte_neige

end module fonte_neige_m
