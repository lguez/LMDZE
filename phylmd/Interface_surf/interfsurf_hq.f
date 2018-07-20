module interfsurf_hq_m

  implicit none

contains

  SUBROUTINE interfsurf_hq(dtime, julien, rmu0, nisurf, knindex, debut, &
       tsoil, qsol, u1_lay, v1_lay, temp_air, spechum, tq_cdrag, petAcoef, &
       peqAcoef, petBcoef, peqBcoef, precip_rain, precip_snow, rugos, rugoro, &
       snow, qsurf, ts, p1lay, ps, radsol, evap, flux_t, fluxlat, dflux_l, &
       dflux_s, tsurf_new, albedo, z0_new, pctsrf_new_sic, agesno, fqcalving, &
       ffonte, run_off_lic_0)

    ! Cette routine sert d'aiguillage entre l'atmosph\`ere et la surface
    ! en g\'en\'eral (sols continentaux, oc\'eans, glaces) pour les flux de
    ! chaleur et d'humidit\'e.

    ! Laurent Fairhead, February 2000

    USE abort_gcm_m, ONLY: abort_gcm
    use alboc_cd_m, only: alboc_cd
    USE albsno_m, ONLY: albsno
    USE calcul_fluxs_m, ONLY: calcul_fluxs
    USE dimphy, ONLY: klon
    USE fonte_neige_m, ONLY: fonte_neige
    USE indicesol, ONLY: epsfra, is_lic, is_oce, is_sic, is_ter
    USE interface_surf, ONLY: conf_interface
    USE interfsur_lim_m, ONLY: interfsur_lim
    use read_sst_m, only: read_sst
    use soil_m, only: soil
    USE suphec_m, ONLY: rcpd, rtt

    real, intent(IN):: dtime ! pas de temps de la physique (en s)
    integer, intent(IN):: julien ! jour dans l'annee en cours
    real, intent(IN):: rmu0(klon) ! cosinus de l'angle solaire zenithal
    integer, intent(IN):: nisurf ! index de la surface a traiter

    integer, intent(in):: knindex(:) ! (knon)
    ! index des points de la surface a traiter

    logical, intent(IN):: debut ! 1er appel a la physique
    ! (si false calcul simplifie des fluxs sur les continents)

    REAL, intent(inout):: tsoil(:, :) ! (knon, nsoilmx)

    REAL, intent(INOUT):: qsol(:) ! (knon)
    ! column-density of water in soil, in kg m-2

    real, intent(IN):: u1_lay(:), v1_lay(:) ! (knon) vitesse 1ere couche

    real, dimension(klon), intent(IN):: temp_air, spechum
    ! temp_air temperature de l'air 1ere couche
    ! spechum humidite specifique 1ere couche
    real, intent(IN):: tq_cdrag(:) ! (knon) coefficient d'echange

    real, dimension(klon), intent(IN):: petAcoef, peqAcoef
    ! coefficients A de la r\'esolution de la couche limite pour t et q

    real, dimension(klon), intent(IN):: petBcoef, peqBcoef
    ! coefficients B de la r\'esolution de la couche limite pour t et q

    real, intent(IN):: precip_rain(klon)
    ! precipitation, liquid water mass flux (kg / m2 / s), positive down

    real, intent(IN):: precip_snow(klon)
    ! precipitation, solid water mass flux (kg / m2 / s), positive down

    real, intent(IN):: rugos(klon) ! rugosite
    real, intent(IN):: rugoro(klon) ! rugosite orographique
    real, intent(INOUT):: snow(:) ! (knon)
    real, intent(INOUT):: qsurf(klon)
    real, intent(IN):: ts(:) ! (knon) temp\'erature de surface
    real, intent(IN):: p1lay(klon) ! pression 1er niveau (milieu de couche)
    real, intent(IN):: ps(:) ! (knon) pression au sol
    REAL, INTENT(IN):: radsol(:) ! (knon) rayonnement net au sol (LW + SW)
    real, intent(OUT):: evap(:) ! (knon) evaporation totale

    real, intent(OUT):: flux_t(:) ! (knon) flux de chaleur sensible
    ! (Cp T) à la surface, positif vers le bas, W / m2

    real, intent(OUT):: fluxlat(:) ! (knon) flux de chaleur latente
    real, intent(OUT):: dflux_l(:), dflux_s(:) ! (knon)
    real, intent(OUT):: tsurf_new(:) ! (knon) temp\'erature au sol
    real, intent(OUT):: albedo(:) ! (knon) albedo
    real, intent(OUT):: z0_new(klon) ! surface roughness

    real, intent(in):: pctsrf_new_sic(:) ! (klon) 
    ! nouvelle repartition des surfaces

    real, intent(INOUT):: agesno(:) ! (knon)

    ! Flux d'eau "perdue" par la surface et n\'ecessaire pour limiter la
    ! hauteur de neige, en kg / m2 / s
    real, intent(OUT):: fqcalving(:) ! (knon)

    ! Flux thermique utiliser pour fondre la neige
    real, dimension(klon), intent(INOUT):: ffonte

    real, dimension(klon), intent(INOUT):: run_off_lic_0
    ! run_off_lic_0 runoff glacier du pas de temps precedent

    ! Local:
    integer knon ! nombre de points de la surface a traiter
    REAL soilcap(size(knindex)) ! (knon)
    REAL soilflux(size(knindex)) ! (knon)
    logical:: first_call = .true.
    integer ii
    real cal(size(knindex)) ! (knon)
    real beta(size(knindex)) ! (knon) evap reelle
    real dif_grnd(klon)
    real tsurf(size(knindex)) ! (knon)
    real alb_neig(size(knindex)) ! (knon)
    real zfra(size(knindex)) ! (knon)
    REAL, PARAMETER:: fmagic = 1. ! facteur magique pour r\'egler l'alb\'edo
    REAL, PARAMETER:: max_eau_sol = 150. ! in kg m-2
    REAL, PARAMETER:: tau_gl = 86400. * 5.

    !-------------------------------------------------------------

    knon = size(knindex)

    ! On doit commencer par appeler les schemas de surfaces continentales
    ! car l'ocean a besoin du ruissellement qui est y calcule

    if (first_call) then
       call conf_interface

       if (nisurf /= is_ter .and. klon > 1) then
          print *, ' nisurf = ', nisurf, ' /= is_ter = ', is_ter
          print *, 'or on doit commencer par les surfaces continentales'
          call abort_gcm("interfsurf_hq", &
               'On doit commencer par les surfaces continentales')
       endif

       if (is_oce > is_sic) then
          print *, 'is_oce = ', is_oce, '> is_sic = ', is_sic
          call abort_gcm("interfsurf_hq", &
               "L'ocean doit etre traite avant la banquise")
       endif

       first_call = .false.
    endif

    ! Initialisations diverses

    ffonte(1:knon) = 0.
    dif_grnd = 999999.
    z0_new = 999999.

    ! Aiguillage vers les differents schemas de surface

    select case (nisurf)
    case (is_ter)
       ! Surface "terre", appel \`a l'interface avec les sols continentaux

       ! Calcul age de la neige

       ! Read albedo from the file containing boundary conditions then
       ! add the albedo of snow:

       call interfsur_lim(dtime, julien, knindex, debut, albedo, z0_new)

       beta = min(2. * qsol / max_eau_sol, 1.)
       dif_grnd(:knon) = 0.
       CALL soil(dtime, is_ter, snow, ts, tsoil, soilcap, soilflux)
       cal = RCPD / soilcap

       CALL calcul_fluxs(dtime, ts, p1lay(:knon), cal, beta, tq_cdrag, ps, &
            qsurf(:knon), radsol + soilflux, dif_grnd(:knon), temp_air(:knon), &
            spechum(:knon), u1_lay, v1_lay, petAcoef(:knon), peqAcoef(:knon), &
            petBcoef(:knon), peqBcoef(:knon), tsurf_new, evap, fluxlat, &
            flux_t, dflux_s, dflux_l)
       CALL fonte_neige(is_ter, dtime, precip_rain(:knon), precip_snow(:knon), &
            snow, qsol, tsurf_new, evap, fqcalving, ffonte(:knon), &
            run_off_lic_0(:knon))

       call albsno(dtime, agesno, alb_neig, precip_snow(:knon))
       where (snow < 0.0001) agesno = 0.
       zfra = max(0., min(1., snow / (snow + 10.)))
       albedo = alb_neig * zfra + albedo * (1. - zfra)
       z0_new = sqrt(z0_new**2 + rugoro**2)
    case (is_oce)
       ! Surface "oc\'ean", appel \`a l'interface avec l'oc\'ean

       call read_sst(julien, knindex, tsurf)
       cal = 0.
       beta = 1.
       dif_grnd = 0.
       call calcul_fluxs(dtime, tsurf, p1lay(:knon), cal, beta, tq_cdrag, ps, &
            qsurf(:knon), radsol, dif_grnd(:knon), temp_air(:knon), &
            spechum(:knon), u1_lay, v1_lay, petAcoef(:knon), peqAcoef(:knon), &
            petBcoef(:knon), peqBcoef(:knon), tsurf_new, evap, fluxlat, &
            flux_t, dflux_s, dflux_l)
       agesno = 0.
       albedo = alboc_cd(rmu0(knindex)) * fmagic
       z0_new = sqrt(rugos**2 + rugoro**2)
       fqcalving = 0.
    case (is_sic)
       ! Surface "glace de mer" appel a l'interface avec l'ocean

       DO ii = 1, knon
          IF (pctsrf_new_sic(knindex(ii)) < EPSFRA) then
             snow(ii) = 0.
             tsurf_new(ii) = RTT - 1.8
             tsoil(ii, :) = RTT - 1.8
          else
             tsurf_new(ii) = ts(ii)
          endif
       enddo

       CALL soil(dtime, is_sic, snow, tsurf_new, tsoil, soilcap, soilflux)
       cal = RCPD / soilcap
       dif_grnd = 1. / tau_gl
       tsurf = tsurf_new
       beta = 1.

       CALL calcul_fluxs(dtime, tsurf, p1lay(:knon), cal, beta, tq_cdrag, ps, &
            qsurf(:knon), radsol + soilflux, dif_grnd(:knon), temp_air(:knon), &
            spechum(:knon), u1_lay, v1_lay, petAcoef(:knon), peqAcoef(:knon), &
            petBcoef(:knon), peqBcoef(:knon), tsurf_new, evap, fluxlat, &
            flux_t, dflux_s, dflux_l)
       CALL fonte_neige(is_sic, dtime, precip_rain(:knon), &
            precip_snow(:knon), snow, qsol, tsurf_new, evap, &
            fqcalving, ffonte(:knon), run_off_lic_0(:knon))

       ! Compute the albedo:

       CALL albsno(dtime, agesno, alb_neig, precip_snow(:knon))
       WHERE (snow < 0.0001) agesno = 0.
       zfra = MAX(0., MIN(1., snow / (snow + 10.)))
       albedo = alb_neig * zfra + 0.6 * (1. - zfra)

       z0_new = SQRT(0.002**2 + rugoro**2)
    case (is_lic)
       ! Surface "glacier continentaux" appel a l'interface avec le sol

       CALL soil(dtime, is_lic, snow, ts, tsoil, soilcap, soilflux)
       cal = RCPD / soilcap
       beta = 1.
       dif_grnd = 0.

       call calcul_fluxs(dtime, ts, p1lay(:knon), cal, beta, tq_cdrag, ps, &
            qsurf(:knon), radsol + soilflux, dif_grnd(:knon), temp_air(:knon), &
            spechum(:knon), u1_lay, v1_lay, petAcoef(:knon), peqAcoef(:knon), &
            petBcoef(:knon), peqBcoef(:knon), tsurf_new, evap, fluxlat, &
            flux_t, dflux_s, dflux_l)
       call fonte_neige(is_lic, dtime, precip_rain(:knon), &
            precip_snow(:knon), snow, qsol, tsurf_new, evap, &
            fqcalving, ffonte(:knon), run_off_lic_0(:knon))

       ! calcul albedo
       CALL albsno(dtime, agesno, alb_neig, precip_snow(:knon))
       WHERE (snow < 0.0001) agesno = 0.
       albedo = 0.77

       ! Rugosite
       z0_new = rugoro
    case default
       print *, 'Index surface = ', nisurf
       call abort_gcm("interfsurf_hq", 'Index surface non valable')
    end select

  END SUBROUTINE interfsurf_hq

end module interfsurf_hq_m
