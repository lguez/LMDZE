module clqh_m

  IMPLICIT none

contains

  SUBROUTINE clqh(julien, nisurf, knindex, tsoil, qsol, mu0, rugos, rugoro, &
       u1lay, v1lay, coef, cdragh, t, q, ts, paprs, pplay, delp, radsol, &
       albedo, snow, qsurf, rain_fall, snow_fall, fluxlat, pctsrf_new_sic, &
       agesno, d_t, d_q, tsurf_new, z0_new, flux_t, flux_q, dflux_s, dflux_l, &
       fqcalving, ffonte, run_off_lic_0, run_off_lic)

    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1993 Aug. 18th
    ! Objet : diffusion verticale de "q" et de "h"

    use climb_hq_down_m, only: climb_hq_down
    use climb_hq_up_m, only: climb_hq_up
    USE dimphy, ONLY: klev
    USE interfsurf_hq_m, ONLY: interfsurf_hq
    USE suphec_m, ONLY: rkappa

    integer, intent(in):: julien ! jour de l'annee en cours
    integer, intent(in):: nisurf
    integer, intent(in):: knindex(:) ! (knon)

    REAL, intent(inout):: tsoil(:, :) ! (knon, nsoilmx)
    ! temperature inside the ground, in K, layer 1 nearest to the surface

    REAL, intent(inout):: qsol(:) ! (knon)
    ! column-density of water in soil, in kg m-2

    real, intent(in):: mu0(:) ! (knon) cosinus de l'angle solaire zenithal
    real, intent(in):: rugos(:) ! (knon) rugosite
    REAL, intent(in):: rugoro(:) ! (knon)

    REAL, intent(in):: u1lay(:), v1lay(:) ! (knon)
    ! vitesse de la 1ere couche (m / s)

    REAL, intent(in):: coef(:, 2:) ! (knon, 2:klev)
    ! Le coefficient d'echange (m^2 / s) multiplie par le cisaillement
    ! du vent (dV / dz)

    REAL, intent(in):: cdragh(:) ! (knon) sans unite

    REAL, intent(in):: t(:, :) ! (knon, klev) air temperature, in K
    REAL, intent(in):: q(:, :) ! (knon, klev) humidit\'e sp\'ecifique
    REAL, intent(in):: ts(:) ! (knon) temperature du sol (K)

    REAL, intent(in):: paprs(:, :) ! (knon, klev + 1)
    ! pression \`a l'inter-couche (Pa)

    REAL, intent(in):: pplay(:, :) ! (knon, klev)
    ! pression au milieu de couche (Pa)

    REAL, intent(in):: delp(:, :) ! (knon, klev)
    ! epaisseur de couche en pression (Pa)

    REAL, intent(in):: radsol(:) ! (knon)
    ! surface net downward radiative flux, in W / m2

    REAL, intent(inout):: albedo(:) ! (knon) albedo de la surface

    REAL, intent(inout):: snow(:) ! (knon)
    ! column-density of mass of snow at the surface, in kg m-2

    REAL, intent(out):: qsurf(:) ! (knon)
    ! humidite de l'air au dessus de la surface

    real, intent(in):: rain_fall(:) ! (knon)
    ! liquid water mass flux (kg / m2 / s), positive down

    real, intent(in):: snow_fall(:) ! (knon)
    ! solid water mass flux (kg / m2 / s), positive down

    real, intent(out):: fluxlat(:) ! (knon) flux de chaleur latente, en W m-2
    real, intent(in):: pctsrf_new_sic(:) ! (knon)
    REAL, intent(inout):: agesno(:) ! (knon)
    REAL, intent(out):: d_t(:, :) ! (knon, klev) variation of air temperature t
    REAL, intent(out):: d_q(:, :) ! (knon, klev) incrementation de "q"
    REAL, intent(out):: tsurf_new(:) ! (knon) new surface temperature, in K
    real, intent(out):: z0_new(:) ! (knon)

    REAL, intent(out):: flux_t(:) ! (knon)
    ! (diagnostic) flux de chaleur sensible (Cp T) à la surface,
    ! positif vers le bas, W / m2

    REAL, intent(out):: flux_q(:) ! (knon)
    ! flux de la vapeur d'eau à la surface, en kg / (m**2 s)

    REAL, intent(out):: dflux_s(:) ! (knon) derivee du flux sensible dF / dTs
    REAL, intent(out):: dflux_l(:) ! (knon) derivee du flux latent dF / dTs

    REAL, intent(out):: fqcalving(:) ! (knon)
    ! Flux d'eau "perdue" par la surface et n\'ecessaire pour que limiter la
    ! hauteur de neige, en kg / m2 / s

    REAL, intent(out):: ffonte(:) ! (knon)
    ! flux thermique utilis\'e pour fondre la neige
    
    REAL, intent(inout):: run_off_lic_0(:) ! (knon)
    ! run-off glacier au pas de temps precedent

    REAL, intent(OUT):: run_off_lic(:) ! (knon) ruissellement total

    ! Local:

    INTEGER k
    REAL evap(size(knindex)) ! (knon) evaporation au sol
    REAL, dimension(size(knindex), klev):: cq, dq, ch, dh ! (knon, klev)
    REAL pkf(size(knindex), klev) ! (knon, klev)

    !----------------------------------------------------------------

    forall (k = 1:klev) pkf(:, k) = (paprs(:, 1) / pplay(:, k))**RKAPPA
    ! (La pression de r\'ef\'erence est celle au sol.)

    call climb_hq_down(pkf, cq, dq, ch, dh, paprs, pplay, t, coef, delp, q)
    CALL interfsurf_hq(julien, mu0, nisurf, knindex, tsoil, qsol, u1lay, &
         v1lay, t(:, 1), q(:, 1), cdragh, ch(:, 1), cq(:, 1), dh(:, 1), &
         dq(:, 1), rain_fall, snow_fall, rugos, rugoro, snow, qsurf, ts, &
         pplay(:, 1), paprs(:, 1), radsol, evap, flux_t, fluxlat, dflux_l, &
         dflux_s, tsurf_new, albedo, z0_new, pctsrf_new_sic, agesno, &
         fqcalving, ffonte, run_off_lic_0, run_off_lic)
    flux_q = - evap
    call climb_hq_up(d_t, d_q, cq, dq, ch, dh, flux_t, flux_q, pkf, t, q)

  END SUBROUTINE clqh

end module clqh_m
