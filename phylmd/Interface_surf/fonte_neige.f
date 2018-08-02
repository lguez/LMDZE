module fonte_neige_m

  implicit none

contains

  SUBROUTINE fonte_neige(nisurf, precip_rain, precip_snow, snow, qsol, &
       tsurf_new, evap, fqcalving, ffonte, run_off_lic_0)

    ! Routine de traitement de la fonte de la neige dans le cas du traitement
    ! de sol simplifi\'e

    ! Laurent Fairhead, March, 2001

    use comconst, only: dtphys
    USE indicesol, ONLY: epsfra, is_lic, is_sic, is_ter
    USE conf_interface_m, ONLY: tau_calv
    use nr_util, only: assert_eq
    USE suphec_m, ONLY: rday, rlmlt, rtt

    integer, intent(IN):: nisurf ! surface \`a traiter

    real, intent(IN):: precip_rain(:) ! (knon)
    ! precipitation, liquid water mass flux (kg / m2 / s), positive down

    real, intent(IN):: precip_snow(:) ! (knon)
    ! precipitation, solid water mass flux (kg / m2 / s), positive down

    real, intent(INOUT):: snow(:) ! (knon)
    ! column-density of mass of snow, in kg m-2

    real, intent(INOUT):: qsol(:) ! (knon)
    ! column-density of water in soil, in kg m-2

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
    ! en exces "s'\'ecoule" (calving).

    integer i
    real fq_fonte
    REAL bil_eau_s(size(precip_rain)) ! (knon) in kg m-2
    real snow_evap(size(precip_rain)) ! (knon) in kg m-2 s-1
    REAL, parameter:: chasno = 3.334E5 / (2.3867E6 * 0.15)
    REAL, parameter:: chaice = 3.334E5 / (2.3867E6 * 0.15)
    real, parameter:: max_eau_sol = 150. ! in kg m-2
    real coeff_rel
    REAL, ALLOCATABLE, SAVE:: run_off_lic(:) ! ruissellement total

    !--------------------------------------------------------------------

    knon = assert_eq((/size(precip_rain), size(precip_snow), size(snow), &
         size(qsol), size(tsurf_new), size(evap), size(fqcalving), &
         size(ffonte), size(run_off_lic_0)/), "fonte_neige knon")

    coeff_rel = dtphys / (tau_calv * rday)
    WHERE (precip_snow > 0.) snow = snow + precip_snow * dtphys

    WHERE (evap > 0.)
       snow_evap = MIN(snow / dtphys, evap)
       snow = snow - snow_evap * dtphys
       snow = MAX(0., snow)
    elsewhere
       snow_evap = 0.
    end where

    bil_eau_s = (precip_rain - evap + snow_evap) * dtphys

    ! Y a-t-il fonte de neige ?

    do i = 1, knon
       if ((snow(i) > epsfra .OR. nisurf == is_sic &
            .OR. nisurf == is_lic) .AND. tsurf_new(i) >= RTT) then
          fq_fonte = MIN(MAX((tsurf_new(i) - RTT) / chasno, 0.), snow(i))
          ffonte(i) = fq_fonte * RLMLT / dtphys
          snow(i) = max(0., snow(i) - fq_fonte)
          bil_eau_s(i) = bil_eau_s(i) + fq_fonte
          tsurf_new(i) = tsurf_new(i) - fq_fonte * chasno

          !IM cf. JLD/ GKtest fonte aussi pour la glace
          IF (nisurf == is_sic .OR. nisurf == is_lic) THEN
             fq_fonte = MAX((tsurf_new(i) - RTT) / chaice, 0.)
             ffonte(i) = ffonte(i) + fq_fonte * RLMLT / dtphys
             bil_eau_s(i) = bil_eau_s(i) + fq_fonte
             tsurf_new(i) = RTT
          ENDIF
       else
          ffonte(i) = 0.
       endif

       ! S'il y a une hauteur trop importante de neige, elle s'\'ecoule
       fqcalving(i) = max(0., snow(i) - snow_max) / dtphys
       snow(i) = min(snow(i), snow_max)
    enddo

    IF (nisurf == is_ter) then
       qsol = MIN(qsol + bil_eau_s, max_eau_sol)
    else if (nisurf == is_lic) then
       if (.not. allocated(run_off_lic)) allocate(run_off_lic(knon))
       ! assumes that the fraction of land-ice does not change during the run

       do i = 1, knon
          run_off_lic(i) = (coeff_rel * fqcalving(i)) + &
               (1. - coeff_rel) * run_off_lic_0(i)
          run_off_lic_0(i) = run_off_lic(i)
          run_off_lic(i) = run_off_lic(i) + bil_eau_s(i) / dtphys
       enddo
    endif

  END SUBROUTINE fonte_neige

end module fonte_neige_m
