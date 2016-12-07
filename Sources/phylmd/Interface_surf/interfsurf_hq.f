module interfsurf_hq_m

  implicit none

contains

  SUBROUTINE interfsurf_hq(dtime, jour, rmu0, nisurf, knon, knindex, rlat, &
       debut, tsoil, qsol, u1_lay, v1_lay, temp_air, spechum, tq_cdrag, &
       petAcoef, peqAcoef, petBcoef, peqBcoef, precip_rain, precip_snow, &
       fder, rugos, rugoro, snow, qsurf, tsurf, p1lay, ps, radsol, evap, &
       flux_t, fluxlat, dflux_l, dflux_s, tsurf_new, albedo, z0_new, &
       pctsrf_new_sic, agesno, fqcalving, ffonte, run_off_lic_0)

    ! Cette routine sert d'aiguillage entre l'atmosph\`ere et la surface
    ! en g\'en\'eral (sols continentaux, oc\'eans, glaces) pour les flux de
    ! chaleur et d'humidit\'e.

    ! Laurent Fairhead, February 2000

    USE abort_gcm_m, ONLY: abort_gcm
    use alboc_cd_m, only: alboc_cd
    use alboc_m, only: alboc
    USE albsno_m, ONLY: albsno
    use calbeta_m, only: calbeta
    USE calcul_fluxs_m, ONLY: calcul_fluxs
    use clesphys2, only: soil_model, cycle_diurne
    USE dimphy, ONLY: klon
    USE fonte_neige_m, ONLY: fonte_neige
    USE indicesol, ONLY: epsfra, is_lic, is_oce, is_sic, is_ter
    USE interface_surf, ONLY: run_off_lic, conf_interface
    USE interfsur_lim_m, ONLY: interfsur_lim
    use read_sst_m, only: read_sst
    use soil_m, only: soil
    USE suphec_m, ONLY: rcpd, rtt

    real, intent(IN):: dtime ! pas de temps de la physique (en s)
    integer, intent(IN):: jour ! jour dans l'annee en cours
    real, intent(IN):: rmu0(klon) ! cosinus de l'angle solaire zenithal
    integer, intent(IN):: nisurf ! index de la surface a traiter
    integer, intent(IN):: knon ! nombre de points de la surface a traiter

    integer, intent(in):: knindex(:) ! (knon)
    ! index des points de la surface a traiter

    real, intent(IN):: rlat(klon) ! latitudes

    logical, intent(IN):: debut ! 1er appel a la physique
    ! (si false calcul simplifie des fluxs sur les continents)

    REAL, intent(inout):: tsoil(:, :) ! (knon, nsoilmx)

    REAL, intent(INOUT):: qsol(klon)
    ! column-density of water in soil, in kg m-2

    real, dimension(klon), intent(IN):: u1_lay, v1_lay
    ! u1_lay vitesse u 1ere couche
    ! v1_lay vitesse v 1ere couche
    real, dimension(klon), intent(IN):: temp_air, spechum
    ! temp_air temperature de l'air 1ere couche
    ! spechum humidite specifique 1ere couche
    real, dimension(klon), intent(INOUT):: tq_cdrag
    ! tq_cdrag cdrag
    real, dimension(klon), intent(IN):: petAcoef, peqAcoef
    ! petAcoef coeff. A de la resolution de la CL pour t
    ! peqAcoef coeff. A de la resolution de la CL pour q
    real, dimension(klon), intent(IN):: petBcoef, peqBcoef
    ! petBcoef coeff. B de la resolution de la CL pour t
    ! peqBcoef coeff. B de la resolution de la CL pour q

    real, intent(IN):: precip_rain(klon)
    ! precipitation, liquid water mass flux (kg / m2 / s), positive down

    real, intent(IN):: precip_snow(klon)
    ! precipitation, solid water mass flux (kg / m2 / s), positive down

    REAL, INTENT(INOUT):: fder(klon) ! derivee des flux (pour le couplage)
    real, intent(IN):: rugos(klon) ! rugosite
    real, intent(IN):: rugoro(klon) ! rugosite orographique
    real, intent(INOUT):: snow(klon), qsurf(klon)
    real, intent(IN):: tsurf(:) ! (knon) temp\'erature de surface
    real, dimension(klon), intent(IN):: p1lay
    ! p1lay pression 1er niveau (milieu de couche)
    real, dimension(klon), intent(IN):: ps
    ! ps pression au sol

    REAL, DIMENSION(klon), INTENT(INOUT):: radsol
    ! rayonnement net au sol (LW + SW)

    real, intent(OUT):: evap(:) ! (knon) evaporation totale
    real, intent(OUT):: flux_t(:) ! (knon) flux de chaleur sensible
    real, dimension(klon), intent(OUT):: fluxlat ! flux de chaleur latente
    real, dimension(klon), intent(OUT):: dflux_l, dflux_s
    real, intent(OUT):: tsurf_new(:) ! (knon) temp\'erature au sol
    real, intent(OUT):: albedo(:) ! (knon) albedo
    real, intent(OUT):: z0_new(klon) ! surface roughness

    real, intent(in):: pctsrf_new_sic(:) ! (klon) 
    ! nouvelle repartition des surfaces

    real, intent(INOUT):: agesno(:) ! (knon)

    ! Flux d'eau "perdue" par la surface et n\'ecessaire pour que limiter la
    ! hauteur de neige, en kg / m2 / s
    !jld a rajouter real, dimension(klon), intent(INOUT):: fqcalving
    real, dimension(klon), intent(INOUT):: fqcalving

    ! Flux thermique utiliser pour fondre la neige
    !jld a rajouter real, dimension(klon), intent(INOUT):: ffonte
    real, dimension(klon), intent(INOUT):: ffonte

    real, dimension(klon), intent(INOUT):: run_off_lic_0
    ! run_off_lic_0 runoff glacier du pas de temps precedent

    ! Local:
    REAL soilcap(knon)
    REAL soilflux(knon)
    logical:: first_call = .true.
    integer ii
    real, dimension(klon):: cal, beta, dif_grnd, capsol
    real, parameter:: calice = 1. / (5.1444e6 * 0.15), tau_gl = 86400. * 5.
    real, parameter:: calsno = 1. / (2.3867e6 * 0.15)
    real tsurf_temp(knon)
    real alb_neig(knon)
    real zfra(knon)
    REAL, PARAMETER:: fmagic = 1. ! facteur magique pour r\'egler l'alb\'edo

    !-------------------------------------------------------------

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
    fqcalving(1:knon) = 0.
    cal = 999999.
    beta = 999999.
    dif_grnd = 999999.
    capsol = 999999.
    z0_new = 999999.
    tsurf_new = 999999.

    ! Aiguillage vers les differents schemas de surface

    select case (nisurf)
    case (is_ter)
       ! Surface "terre", appel \`a l'interface avec les sols continentaux

       ! Calcul age de la neige

       ! Read albedo from the file containing boundary conditions then
       ! add the albedo of snow:

       call interfsur_lim(dtime, jour, knindex, debut, albedo, z0_new)

       ! Calcul snow et qsurf, hydrologie adapt\'ee
       CALL calbeta(is_ter, snow(:knon), qsol(:knon), beta(:knon), &
            capsol(:knon), dif_grnd(:knon))

       IF (soil_model) THEN
          CALL soil(dtime, is_ter, snow(:knon), tsurf, tsoil, soilcap, soilflux)
          cal(1:knon) = RCPD / soilcap
          radsol(1:knon) = radsol(1:knon) + soilflux
       ELSE
          cal = RCPD * capsol
       ENDIF

       CALL calcul_fluxs(dtime, tsurf, p1lay(:knon), cal(:knon), &
            beta(:knon), tq_cdrag(:knon), ps(:knon), qsurf(:knon), &
            radsol(:knon), dif_grnd(:knon), temp_air(:knon), spechum(:knon), &
            u1_lay(:knon), v1_lay(:knon), petAcoef(:knon), peqAcoef(:knon), &
            petBcoef(:knon), peqBcoef(:knon), tsurf_new, evap, &
            fluxlat(:knon), flux_t, dflux_s(:knon), dflux_l(:knon))

       CALL fonte_neige(is_ter, dtime, tsurf, p1lay(:knon), beta(:knon), &
            tq_cdrag(:knon), ps(:knon), precip_rain(:knon), &
            precip_snow(:knon), snow(:knon), qsol(:knon), temp_air(:knon), &
            spechum(:knon), u1_lay(:knon), v1_lay(:knon), petAcoef(:knon), &
            peqAcoef(:knon), petBcoef(:knon), peqBcoef(:knon), tsurf_new, &
            evap, fqcalving(:knon), ffonte(:knon), run_off_lic_0(:knon))

       call albsno(dtime, agesno, alb_neig, precip_snow(:knon))
       where (snow(:knon) < 0.0001) agesno = 0.
       zfra = max(0., min(1., snow(:knon) / (snow(:knon) + 10.)))
       albedo = alb_neig * zfra + albedo * (1. - zfra)
       z0_new = sqrt(z0_new**2 + rugoro**2)
    case (is_oce)
       ! Surface "oc\'ean", appel \`a l'interface avec l'oc\'ean

       call read_sst(dtime, jour, knindex, debut, tsurf_temp)

       cal = 0.
       beta = 1.
       dif_grnd = 0.
       agesno = 0.
       call calcul_fluxs(dtime, tsurf_temp, p1lay(:knon), cal(:knon), &
            beta(:knon), tq_cdrag(:knon), ps(:knon), qsurf(:knon), &
            radsol(:knon), dif_grnd(:knon), temp_air(:knon), spechum(:knon), &
            u1_lay(:knon), v1_lay(:knon), petAcoef(:knon), peqAcoef(:knon), &
            petBcoef(:knon), peqBcoef(:knon), tsurf_new, evap, &
            fluxlat(:knon), flux_t, dflux_s(:knon), dflux_l(:knon))
       fder = fder + dflux_s + dflux_l

       ! Compute the albedo:

       if (cycle_diurne) then
          albedo = alboc_cd(rmu0(knindex))
       else
          albedo = alboc(jour, rlat(knindex))
       endif

       albedo = albedo * fmagic

       z0_new = sqrt(rugos**2 + rugoro**2)
    case (is_sic)
       ! Surface "glace de mer" appel a l'interface avec l'ocean

       DO ii = 1, knon
          tsurf_new(ii) = tsurf(ii)
          IF (pctsrf_new_sic(knindex(ii)) < EPSFRA) then
             snow(ii) = 0.
             tsurf_new(ii) = RTT - 1.8
             IF (soil_model) tsoil(ii, :) = RTT - 1.8
          endif
       enddo

       CALL calbeta(is_sic, snow(:knon), qsol(:knon), beta(:knon), &
            capsol(:knon), dif_grnd(:knon))

       IF (soil_model) THEN
          CALL soil(dtime, is_sic, snow(:knon), tsurf_new, tsoil, soilcap, &
               soilflux)
          cal(1:knon) = RCPD / soilcap
          radsol(1:knon) = radsol(1:knon) + soilflux
          dif_grnd = 0.
       ELSE
          dif_grnd = 1. / tau_gl
          cal = RCPD * calice
          WHERE (snow > 0.) cal = RCPD * calsno
       ENDIF
       tsurf_temp = tsurf_new
       beta = 1.

       CALL calcul_fluxs(dtime, tsurf_temp, p1lay(:knon), cal(:knon), &
            beta(:knon), tq_cdrag(:knon), ps(:knon), qsurf(:knon), &
            radsol(:knon), dif_grnd(:knon), temp_air(:knon), spechum(:knon), &
            u1_lay(:knon), v1_lay(:knon), petAcoef(:knon), peqAcoef(:knon), &
            petBcoef(:knon), peqBcoef(:knon), tsurf_new, evap, &
            fluxlat(:knon), flux_t, dflux_s(:knon), dflux_l(:knon))

       CALL fonte_neige(is_sic, dtime, tsurf_temp, p1lay(:knon), beta(:knon), &
            tq_cdrag(:knon), ps(:knon), precip_rain(:knon), &
            precip_snow(:knon), snow(:knon), qsol(:knon), temp_air(:knon), &
            spechum(:knon), u1_lay(:knon), v1_lay(:knon), petAcoef(:knon), &
            peqAcoef(:knon), petBcoef(:knon), peqBcoef(:knon), tsurf_new, &
            evap, fqcalving(:knon), ffonte(:knon), run_off_lic_0(:knon))

       ! Compute the albedo:

       CALL albsno(dtime, agesno, alb_neig, precip_snow(:knon))
       WHERE (snow(:knon) < 0.0001) agesno = 0.
       zfra = MAX(0., MIN(1., snow(:knon) / (snow(:knon) + 10.)))
       albedo = alb_neig * zfra + 0.6 * (1. - zfra)

       fder = fder + dflux_s + dflux_l
       z0_new = SQRT(0.002**2 + rugoro**2)
    case (is_lic)
       if (.not. allocated(run_off_lic)) then
          allocate(run_off_lic(knon))
          run_off_lic = 0.
       endif

       ! Surface "glacier continentaux" appel a l'interface avec le sol

       IF (soil_model) THEN
          CALL soil(dtime, is_lic, snow(:knon), tsurf, tsoil, soilcap, soilflux)
          cal(1:knon) = RCPD / soilcap
          radsol(1:knon) = radsol(1:knon) + soilflux
       ELSE
          cal = RCPD * calice
          WHERE (snow > 0.) cal = RCPD * calsno
       ENDIF
       beta = 1.
       dif_grnd = 0.

       call calcul_fluxs(dtime, tsurf, p1lay(:knon), cal(:knon), &
            beta(:knon), tq_cdrag(:knon), ps(:knon), qsurf(:knon), &
            radsol(:knon), dif_grnd(:knon), temp_air(:knon), spechum(:knon), &
            u1_lay(:knon), v1_lay(:knon), petAcoef(:knon), peqAcoef(:knon), &
            petBcoef(:knon), peqBcoef(:knon), tsurf_new, evap, &
            fluxlat(:knon), flux_t, dflux_s(:knon), dflux_l(:knon))

       call fonte_neige(is_lic, dtime, tsurf, p1lay(:knon), beta(:knon), &
            tq_cdrag(:knon), ps(:knon), precip_rain(:knon), &
            precip_snow(:knon), snow(:knon), qsol(:knon), temp_air(:knon), &
            spechum(:knon), u1_lay(:knon), v1_lay(:knon), petAcoef(:knon), &
            peqAcoef(:knon), petBcoef(:knon), peqBcoef(:knon), tsurf_new, &
            evap, fqcalving(:knon), ffonte(:knon), run_off_lic_0(:knon))

       ! calcul albedo
       CALL albsno(dtime, agesno, alb_neig, precip_snow(:knon))
       WHERE (snow(:knon) < 0.0001) agesno = 0.
       albedo = 0.77

       ! Rugosite
       z0_new = rugoro
    case default
       print *, 'Index surface = ', nisurf
       call abort_gcm("interfsurf_hq", 'Index surface non valable')
    end select

  END SUBROUTINE interfsurf_hq

end module interfsurf_hq_m
