module clqh_m

  IMPLICIT none

contains

  SUBROUTINE clqh(julien, nisurf, knindex, tsoil, qsol, mu0, rugos, rugoro, &
       u1lay, v1lay, coefh, cdragh, t, q, ts, paprs, pplay, delp, radsol, &
       albedo, snow, qsurf, rain_fall, snow_fall, fluxlat, pctsrf_new_sic, &
       agesno, d_t, d_q, tsurf_new, z0_new, flux_t, flux_q, dflux_s, dflux_l, &
       fqcalving, ffonte, run_off_lic_0, run_off_lic)

    ! Authors: Z. X. Li (LMD/CNRS), Laurent Fairhead
    ! Dates: 1993 Aug. 18th, February 2000
    ! Objet : diffusion verticale de "q" et de "h"

    USE abort_gcm_m, ONLY: abort_gcm
    use alboc_cd_m, only: alboc_cd
    USE albsno_m, ONLY: albsno
    USE calcul_fluxs_m, ONLY: calcul_fluxs
    use climb_hq_down_m, only: climb_hq_down
    use climb_hq_up_m, only: climb_hq_up
    USE dimphy, ONLY: klev
    USE fonte_neige_m, ONLY: fonte_neige
    USE indicesol, ONLY: epsfra, is_lic, is_oce, is_sic, is_ter
    USE interfsur_lim_m, ONLY: interfsur_lim
    use limit_read_sst_m, only: limit_read_sst
    use soil_m, only: soil
    USE suphec_m, ONLY: rcpd, rtt, rkappa

    integer, intent(in):: julien ! jour de l'annee en cours
    integer, intent(in):: nisurf ! index de la surface a traiter

    integer, intent(in):: knindex(:) ! (knon)
    ! index des points de la surface a traiter

    REAL, intent(inout):: tsoil(:, :) ! (knon, nsoilmx)
    ! temperature inside the ground, in K, layer 1 nearest to the surface

    REAL, intent(inout):: qsol(:) ! (knon)
    ! column-density of water in soil, in kg m-2

    real, intent(in):: mu0(:) ! (knon) cosinus de l'angle solaire zenithal
    real, intent(in):: rugos(:) ! (knon) rugosite
    REAL, intent(in):: rugoro(:) ! (knon) rugosite orographique

    REAL, intent(in):: u1lay(:), v1lay(:) ! (knon)
    ! vitesse de la 1ere couche (m / s)

    REAL, intent(in):: coefh(:, 2:) ! (knon, 2:klev)
    ! diffusion coefficient at layer interface, for heat and humidity, in m2 s-1

    REAL, intent(in):: cdragh(:) ! (knon) coefficient d'\'echange, sans unite

    REAL, intent(in):: t(:, :) ! (knon, klev) air temperature, in K
    REAL, intent(in):: q(:, :) ! (knon, klev) humidit\'e sp\'ecifique
    REAL, intent(in):: ts(:) ! (knon) temperature de surface (K)

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
    ! precipitation, liquid water mass flux (kg / m2 / s), positive down

    real, intent(in):: snow_fall(:) ! (knon)
    ! precipitation, solid water mass flux (kg / m2 / s), positive down

    real, intent(out):: fluxlat(:) ! (knon) flux de chaleur latente, en W m-2

    real, intent(in):: pctsrf_new_sic(:) ! (knon)
    ! nouvelle repartition des surfaces

    REAL, intent(inout):: agesno(:) ! (knon)
    REAL, intent(out):: d_t(:, :) ! (knon, klev) variation of air temperature t
    REAL, intent(out):: d_q(:, :) ! (knon, klev) variation of q
    REAL, intent(out):: tsurf_new(:) ! (knon) new surface temperature, in K
    real, intent(out):: z0_new(:) ! (knon) surface roughness

    REAL, intent(out):: flux_t(:) ! (knon)
    ! (diagnostic) flux de chaleur sensible (Cp T) à la surface,
    ! positif vers le bas, W / m2

    REAL, intent(out):: flux_q(:) ! (knon)
    ! flux de la vapeur d'eau à la surface, en kg / (m**2 s)

    REAL, intent(out):: dflux_s(:) ! (knon) derivee du flux sensible dF / dTs
    REAL, intent(out):: dflux_l(:) ! (knon) derivee du flux latent dF / dTs

    REAL, intent(out):: fqcalving(:) ! (knon)
    ! Flux d'eau "perdue" par la surface et n\'ecessaire pour limiter la
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
    ! coefficients de la r\'esolution de la couche limite pour t et q

    REAL pkf(size(knindex), klev) ! (knon, klev)
    REAL soilcap(size(knindex)) ! (knon)
    REAL soilflux(size(knindex)) ! (knon)
    integer i, knon
    real tsurf(size(knindex)) ! (knon)
    real alb_neig(size(knindex)) ! (knon)
    real zfra(size(knindex)) ! (knon) fraction of surface covered by snow
    REAL, PARAMETER:: fmagic = 1. ! facteur magique pour r\'egler l'alb\'edo
    REAL, PARAMETER:: max_eau_sol = 150. ! in kg m-2

    REAL, PARAMETER:: tau_gl = 86400. * 5.
    ! constante de rappel de la temp\'erature \`a la surface de la glace, en s

    !----------------------------------------------------------------

    forall (k = 1:klev) pkf(:, k) = (paprs(:, 1) / pplay(:, k))**RKAPPA
    ! (La pression de r\'ef\'erence est celle au sol.)

    call climb_hq_down(pkf, cq, dq, ch, dh, paprs, pplay, t, coefh, delp, q)
    knon  = size(knindex)

    select case (nisurf)
    case (is_ter)
       ! Surface "terre", appel \`a l'interface avec les sols continentaux

       CALL soil(is_ter, snow, ts, tsoil, soilcap, soilflux)
       CALL calcul_fluxs(qsurf, tsurf_new, evap, fluxlat, flux_t, dflux_s, &
            dflux_l, ts, pplay(:, 1), cdragh, paprs(:, 1), radsol + soilflux, &
            t(:, 1), q(:, 1), u1lay, v1lay, ch(:, 1), cq(:, 1), dh(:, 1), &
            dq(:, 1), cal = RCPD / soilcap, &
            beta = min(2. * qsol / max_eau_sol, 1.), dif_grnd = 0.)
       CALL fonte_neige(is_ter, rain_fall, snow_fall, snow, qsol, tsurf_new, &
            evap, fqcalving, ffonte, run_off_lic_0, run_off_lic)
       call albsno(agesno, alb_neig, snow_fall)
       where (snow < 1e-4) agesno = 0.
       zfra = snow / (snow + 10.)

       ! Read albedo from the file containing boundary conditions then
       ! add the albedo of snow:
       call interfsur_lim(julien, knindex, albedo, z0_new)
       albedo = alb_neig * zfra + albedo * (1. - zfra)

       z0_new = sqrt(z0_new**2 + rugoro**2)
    case (is_oce)
       ! Surface oc\'ean, appel \`a l'interface avec l'oc\'ean

       ffonte = 0.
       call limit_read_sst(julien, knindex, tsurf)
       call calcul_fluxs(qsurf, tsurf_new, evap, fluxlat, flux_t, dflux_s, &
            dflux_l, tsurf, pplay(:, 1), cdragh, paprs(:, 1), radsol, t(:, 1), &
            q(:, 1), u1lay, v1lay, ch(:, 1), cq(:, 1), dh(:, 1), dq(:, 1), &
            cal = [(0., i = 1, knon)], beta = [(1., i = 1, knon)], &
            dif_grnd = 0.)
       agesno = 0.
       albedo = alboc_cd(mu0) * fmagic
       z0_new = sqrt(rugos**2 + rugoro**2)
       fqcalving = 0.
    case (is_sic)
       ! Surface glace de mer

       DO i = 1, knon
          IF (pctsrf_new_sic(i) < EPSFRA) then
             snow(i) = 0.
             tsurf(i) = RTT - 1.8
             tsoil(i, :) = tsurf(i)
          else
             tsurf(i) = ts(i)
          endif
       enddo

       CALL soil(is_sic, snow, tsurf, tsoil, soilcap, soilflux)
       CALL calcul_fluxs(qsurf, tsurf_new, evap, fluxlat, flux_t, dflux_s, &
            dflux_l, tsurf, pplay(:, 1), cdragh, paprs(:, 1), &
            radsol + soilflux, t(:, 1), q(:, 1), u1lay, v1lay, ch(:, 1), &
            cq(:, 1), dh(:, 1), dq(:, 1), cal = RCPD / soilcap, &
            beta = [(1., i = 1, knon)], dif_grnd = 1. / tau_gl)
       CALL fonte_neige(is_sic, rain_fall, snow_fall, snow, qsol, &
            tsurf_new, evap, fqcalving, ffonte, run_off_lic_0, run_off_lic)

       ! Compute the albedo:

       CALL albsno(agesno, alb_neig, snow_fall)
       WHERE (snow < 1e-4) agesno = 0.
       zfra = snow / (snow + 10.)
       albedo = alb_neig * zfra + 0.6 * (1. - zfra)

       z0_new = SQRT(0.002**2 + rugoro**2)
    case (is_lic)
       ! Surface "glaciers continentaux" appel \`a l'interface avec le sol

       CALL soil(is_lic, snow, ts, tsoil, soilcap, soilflux)
       call calcul_fluxs(qsurf, tsurf_new, evap, fluxlat, flux_t, dflux_s, &
            dflux_l, ts, pplay(:, 1), cdragh, paprs(:, 1), radsol + soilflux, &
            t(:, 1), q(:, 1), u1lay, v1lay, ch(:, 1), cq(:, 1), dh(:, 1), &
            dq(:, 1), cal = RCPD / soilcap, beta = [(1., i = 1, knon)], &
            dif_grnd = 0.)
       call fonte_neige(is_lic, rain_fall, snow_fall, snow, qsol, &
            tsurf_new, evap, fqcalving, ffonte, run_off_lic_0, run_off_lic)

       ! calcul albedo
       CALL albsno(agesno, alb_neig, snow_fall)
       WHERE (snow < 1e-4) agesno = 0.
       albedo = 0.77

       ! Rugosite
       z0_new = rugoro
    case default
       print *, 'Index of surface = ', nisurf
       call abort_gcm("clqh", 'Index surface non valable')
    end select

    flux_q = - evap
    call climb_hq_up(d_t, d_q, cq, dq, ch, dh, flux_t, flux_q, pkf, t, q)

  END SUBROUTINE clqh

end module clqh_m
