module interfsurf_hq_m

  implicit none

contains

  SUBROUTINE interfsurf_hq(julien, mu0, nisurf, knindex, tsoil, qsol, u1lay, &
       v1lay, t1lay, q1lay, cdragh, tAcoef, qAcoef, tBcoef, qBcoef, &
       rain_fall, snow_fall, rugos, rugoro, snow, qsurf, ts, p1lay, ps, &
       radsol, evap, flux_t, fluxlat, dflux_l, dflux_s, tsurf_new, albedo, &
       z0_new, pctsrf_new_sic, agesno, fqcalving, ffonte, run_off_lic_0, &
       run_off_lic)

    ! Cette routine sert d'aiguillage entre l'atmosph\`ere et la surface
    ! en g\'en\'eral (sols continentaux, oc\'eans, glaces) pour les flux de
    ! chaleur et d'humidit\'e.

    ! Laurent Fairhead, February 2000

    USE abort_gcm_m, ONLY: abort_gcm
    use alboc_cd_m, only: alboc_cd
    USE albsno_m, ONLY: albsno
    USE calcul_fluxs_m, ONLY: calcul_fluxs
    USE fonte_neige_m, ONLY: fonte_neige
    USE indicesol, ONLY: epsfra, is_lic, is_oce, is_sic, is_ter
    USE interfsur_lim_m, ONLY: interfsur_lim
    use limit_read_sst_m, only: limit_read_sst
    use soil_m, only: soil
    USE suphec_m, ONLY: rcpd, rtt

    integer, intent(IN):: julien ! jour dans l'annee en cours
    real, intent(IN):: mu0(:) ! (knon) cosinus de l'angle solaire zenithal
    integer, intent(IN):: nisurf ! index de la surface a traiter

    integer, intent(in):: knindex(:) ! (knon)
    ! index des points de la surface a traiter

    REAL, intent(inout):: tsoil(:, :) ! (knon, nsoilmx)
    ! temperature inside the ground, in K, layer 1 nearest to the surface

    REAL, intent(INOUT):: qsol(:) ! (knon)
    ! column-density of water in soil, in kg m-2

    real, intent(IN):: u1lay(:), v1lay(:) ! (knon) vitesse 1ere couche
    real, intent(IN):: t1lay(:) ! (knon) temp\'erature de l'air 1\`ere couche

    real, intent(IN):: q1lay(:) ! (knon)
    ! humidit\'e sp\'ecifique de la premi\`ere couche

    real, intent(IN):: cdragh(:) ! (knon) coefficient d'\'echange

    real, intent(IN):: tAcoef(:), qAcoef(:) ! (knon)
    ! coefficients A de la r\'esolution de la couche limite pour t et q

    real, intent(IN):: tBcoef(:), qBcoef(:) ! (knon)
    ! coefficients B de la r\'esolution de la couche limite pour t et q

    real, intent(IN):: rain_fall(:) ! (knon)
    ! precipitation, liquid water mass flux (kg / m2 / s), positive down

    real, intent(IN):: snow_fall(:) ! (knon)
    ! precipitation, solid water mass flux (kg / m2 / s), positive down

    real, intent(IN):: rugos(:) ! (knon) rugosite
    real, intent(IN):: rugoro(:) ! (knon) rugosite orographique

    real, intent(INOUT):: snow(:) ! (knon)
    ! column-density of mass of snow at the surface, in kg m-2
    
    real, intent(OUT):: qsurf(:) ! (knon)
    real, intent(IN):: ts(:) ! (knon) temp\'erature de surface
    real, intent(IN):: p1lay(:) ! (knon) pression 1er niveau (milieu de couche)
    real, intent(IN):: ps(:) ! (knon) pression au sol, en Pa

    REAL, INTENT(IN):: radsol(:) ! (knon)
    ! surface net downward radiative flux, in W / m2

    real, intent(OUT):: evap(:) ! (knon) evaporation totale

    real, intent(OUT):: flux_t(:) ! (knon) flux de chaleur sensible
    ! (Cp T) à la surface, positif vers le bas, W / m2

    real, intent(OUT):: fluxlat(:) ! (knon) flux de chaleur latente, en W m-2
    real, intent(OUT):: dflux_l(:), dflux_s(:) ! (knon)
    real, intent(OUT):: tsurf_new(:) ! (knon) temp\'erature au sol
    real, intent(OUT):: albedo(:) ! (knon) albedo
    real, intent(OUT):: z0_new(:) ! (knon) surface roughness

    real, intent(in):: pctsrf_new_sic(:) ! (knon) 
    ! nouvelle repartition des surfaces

    real, intent(INOUT):: agesno(:) ! (knon)

    real, intent(OUT):: fqcalving(:) ! (knon)
    ! Flux d'eau "perdue" par la surface et n\'ecessaire pour limiter la
    ! hauteur de neige, en kg / m2 / s

    real, intent(OUT):: ffonte(:) ! (knon)
    ! flux thermique utilis\'e pour fondre la neige

    real, intent(INOUT):: run_off_lic_0(:) ! (knon)
    ! run_off_lic_0 runoff glacier du pas de temps precedent

    REAL, intent(OUT):: run_off_lic(:) ! (knon) ruissellement total

    ! Local:

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

    !-------------------------------------------------------------

    knon  = size(knindex)

    select case (nisurf)
    case (is_ter)
       ! Surface "terre", appel \`a l'interface avec les sols continentaux

       CALL soil(is_ter, snow, ts, tsoil, soilcap, soilflux)
       CALL calcul_fluxs(qsurf, tsurf_new, evap, fluxlat, flux_t, dflux_s, &
            dflux_l, ts, p1lay, cdragh, ps, radsol + soilflux, t1lay, q1lay, &
            u1lay, v1lay, tAcoef, qAcoef, tBcoef, qBcoef, &
            cal = RCPD / soilcap, beta = min(2. * qsol / max_eau_sol, 1.), &
            dif_grnd = 0.)
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
            dflux_l, tsurf, p1lay, cdragh, ps, radsol, t1lay, q1lay, u1lay, &
            v1lay, tAcoef, qAcoef, tBcoef, qBcoef, cal = [(0., i = 1, knon)], &
            beta = [(1., i = 1, knon)], dif_grnd = 0.)
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
            dflux_l, tsurf, p1lay, cdragh, ps, radsol + soilflux, t1lay, &
            q1lay, u1lay, v1lay, tAcoef, qAcoef, tBcoef, qBcoef, &
            cal = RCPD / soilcap, beta = [(1., i = 1, knon)], &
            dif_grnd = 1. / tau_gl)
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
            dflux_l, ts, p1lay, cdragh, ps, radsol + soilflux, t1lay, q1lay, &
            u1lay, v1lay, tAcoef, qAcoef, tBcoef, qBcoef, &
            cal = RCPD / soilcap, beta = [(1., i = 1, knon)], dif_grnd = 0.)
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
       call abort_gcm("interfsurf_hq", 'Index surface non valable')
    end select

  END SUBROUTINE interfsurf_hq

end module interfsurf_hq_m
