module interfsurf_hq_m

  implicit none

contains

  SUBROUTINE interfsurf_hq(itime, dtime, jour, rmu0, nisurf, knon, knindex, &
       pctsrf, rlat, debut, nsoilmx, tsoil, qsol, u1_lay, v1_lay, temp_air, &
       spechum, tq_cdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
       precip_rain, precip_snow, fder, rugos, rugoro, snow, qsurf, tsurf, &
       p1lay, ps, radsol, evap, fluxsens, fluxlat, dflux_l, dflux_s, &
       tsurf_new, albedo, z0_new, pctsrf_new, agesno, fqcalving, ffonte, &
       run_off_lic_0, flux_o, flux_g)

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
    USE indicesol, ONLY: epsfra, is_lic, is_oce, is_sic, is_ter, nbsrf
    USE interface_surf, ONLY: run_off, run_off_lic, conf_interface
    USE interfoce_lim_m, ONLY: interfoce_lim
    USE interfsur_lim_m, ONLY: interfsur_lim
    use soil_m, only: soil
    USE suphec_m, ONLY: rcpd, rlstt, rlvtt, rtt

    integer, intent(IN):: itime ! numero du pas de temps
    real, intent(IN):: dtime ! pas de temps de la physique (en s)
    integer, intent(IN):: jour ! jour dans l'annee en cours
    real, intent(IN):: rmu0(klon) ! cosinus de l'angle solaire zenithal
    integer, intent(IN):: nisurf ! index de la surface a traiter
    integer, intent(IN):: knon ! nombre de points de la surface a traiter

    integer, intent(in):: knindex(:) ! (knon)
    ! index des points de la surface a traiter

    real, intent(IN):: pctsrf(klon, nbsrf)
    ! tableau des pourcentages de surface de chaque maille

    real, intent(IN):: rlat(klon) ! latitudes

    logical, intent(IN):: debut ! 1er appel a la physique
    ! (si false calcul simplifie des fluxs sur les continents)

    integer, intent(in):: nsoilmx
    REAL tsoil(klon, nsoilmx)

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
    ! precipitation, liquid water mass flux (kg/m2/s), positive down

    real, intent(IN):: precip_snow(klon)
    ! precipitation, solid water mass flux (kg/m2/s), positive down

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
    ! radsol rayonnement net aus sol (LW + SW)
    real, intent(INOUT):: evap(klon) ! evaporation totale
    real, dimension(klon), intent(OUT):: fluxsens, fluxlat
    ! fluxsens flux de chaleur sensible
    ! fluxlat flux de chaleur latente
    real, dimension(klon), intent(OUT):: dflux_l, dflux_s
    real, intent(OUT):: tsurf_new(knon) ! temp\'erature au sol
    real, intent(OUT):: albedo(:) ! (knon) albedo
    real, intent(OUT):: z0_new(klon) ! surface roughness
    real, dimension(klon, nbsrf), intent(OUT):: pctsrf_new
    ! pctsrf_new nouvelle repartition des surfaces
    real, dimension(klon), intent(INOUT):: agesno

    ! Flux d'eau "perdue" par la surface et n\'ecessaire pour que limiter la
    ! hauteur de neige, en kg/m2/s
    !jld a rajouter real, dimension(klon), intent(INOUT):: fqcalving
    real, dimension(klon), intent(INOUT):: fqcalving

    ! Flux thermique utiliser pour fondre la neige
    !jld a rajouter real, dimension(klon), intent(INOUT):: ffonte
    real, dimension(klon), intent(INOUT):: ffonte

    real, dimension(klon), intent(INOUT):: run_off_lic_0
    ! run_off_lic_0 runoff glacier du pas de temps precedent

    !IM: "slab" ocean
    real, dimension(klon), intent(OUT):: flux_o, flux_g

    ! Local:

    REAL, dimension(klon):: soilcap
    REAL, dimension(klon):: soilflux

    !IM: "slab" ocean
    real, parameter:: t_grnd=271.35
    integer i

    character (len = 20), save:: modname = 'interfsurf_hq'
    character (len = 80):: abort_message
    logical, save:: first_call = .true.
    integer:: ii
    real, dimension(klon):: cal, beta, dif_grnd, capsol
    real, parameter:: calice=1.0/(5.1444e6 * 0.15), tau_gl=86400. * 5.
    real, parameter:: calsno=1./(2.3867e6 * 0.15)
    real tsurf_temp(knon)
    real alb_neig(klon)
    real zfra(klon)

    !-------------------------------------------------------------

    ! On doit commencer par appeler les schemas de surfaces continentales
    ! car l'ocean a besoin du ruissellement qui est y calcule

    if (first_call) then
       call conf_interface
       if (nisurf /= is_ter .and. klon > 1) then
          print *, ' Warning:'
          print *, ' nisurf = ', nisurf, ' /= is_ter = ', is_ter
          print *, 'or on doit commencer par les surfaces continentales'
          abort_message='voir ci-dessus'
          call abort_gcm(modname, abort_message)
       endif
       if (is_oce > is_sic) then
          print *, 'Warning:'
          print *, ' Pour des raisons de sequencement dans le code'
          print *, ' l''ocean doit etre traite avant la banquise'
          print *, ' or is_oce = ', is_oce, '> is_sic = ', is_sic
          abort_message='voir ci-dessus'
          call abort_gcm(modname, abort_message)
       endif
    endif
    first_call = .false.

    ! Initialisations diverses

    ffonte(1:knon)=0.
    fqcalving(1:knon)=0.
    cal = 999999.
    beta = 999999.
    dif_grnd = 999999.
    capsol = 999999.
    z0_new = 999999.
    alb_neig = 999999.
    tsurf_new = 999999.

    !IM: "slab" ocean; initialisations
    flux_o = 0.
    flux_g = 0.

    ! Aiguillage vers les differents schemas de surface

    select case (nisurf)
    case (is_ter)
       ! Surface "terre", appel \`a l'interface avec les sols continentaux

       ! allocation du run-off
       if (.not. allocated(run_off)) then
          allocate(run_off(knon))
          run_off = 0.
       else if (size(run_off) /= knon) then
          call abort_gcm(modname, 'Something is wrong: the number of ' &
               // 'continental points has changed since last call.')
       endif

       ! Calcul age de la neige

       ! calcul albedo: lecture albedo fichier boundary conditions
       ! puis ajout albedo neige
       call interfsur_lim(itime, dtime, jour, knindex, debut, albedo, z0_new)

       ! calcul snow et qsurf, hydrol adapt\'e
       CALL calbeta(nisurf, snow(:knon), qsol(:knon), beta(:knon), &
            capsol(:knon), dif_grnd(:knon))

       IF (soil_model) THEN
          CALL soil(dtime, nisurf, knon, snow, tsurf, tsoil, soilcap, soilflux)
          cal(1:knon) = RCPD / soilcap(1:knon)
          radsol(1:knon) = radsol(1:knon) + soilflux(:knon)
       ELSE
          cal = RCPD * capsol
       ENDIF

       CALL calcul_fluxs(dtime, tsurf, p1lay(:knon), cal(:knon), &
            beta(:knon), tq_cdrag(:knon), ps(:knon), qsurf(:knon), &
            radsol(:knon), dif_grnd(:knon), temp_air(:knon), spechum(:knon), &
            u1_lay(:knon), v1_lay(:knon), petAcoef(:knon), peqAcoef(:knon), &
            petBcoef(:knon), peqBcoef(:knon), tsurf_new, evap(:knon), &
            fluxlat(:knon), fluxsens(:knon), dflux_s(:knon), dflux_l(:knon))

       CALL fonte_neige(nisurf, dtime, tsurf, p1lay(:knon), beta(:knon), &
            tq_cdrag(:knon), ps(:knon), precip_rain(:knon), &
            precip_snow(:knon), snow(:knon), qsol(:knon), temp_air(:knon), &
            spechum(:knon), u1_lay(:knon), v1_lay(:knon), petAcoef(:knon), &
            peqAcoef(:knon), petBcoef(:knon), peqBcoef(:knon), tsurf_new, &
            evap(:knon), fqcalving(:knon), ffonte(:knon), run_off_lic_0(:knon))

       call albsno(klon, knon, dtime, agesno, alb_neig, precip_snow)
       where (snow(1 : knon) < 0.0001) agesno(1 : knon) = 0.
       zfra(:knon) = max(0.0, min(1.0, snow(1:knon)/(snow(1:knon) + 10.0)))
       albedo = alb_neig(:knon) * zfra(:knon) &
            + albedo * (1. - zfra(:knon))
       z0_new = sqrt(z0_new**2 + rugoro**2)

       ! Remplissage des pourcentages de surface
       pctsrf_new(:, nisurf) = pctsrf(:, nisurf)
    case (is_oce)
       ! Surface "ocean" appel \`a l'interface avec l'oc\'ean
       ! lecture conditions limites
       call interfoce_lim(itime, dtime, jour, knindex, debut, tsurf_temp, &
            pctsrf_new)

       cal = 0.
       beta = 1.
       dif_grnd = 0.
       alb_neig = 0.
       agesno = 0.
       call calcul_fluxs(dtime, tsurf_temp, p1lay(:knon), &
            cal(:knon), beta(:knon), tq_cdrag(:knon), ps(:knon), &
            qsurf(:knon), radsol(:knon), dif_grnd(:knon), temp_air(:knon), &
            spechum(:knon), u1_lay(:knon), v1_lay(:knon), petAcoef(:knon), &
            peqAcoef(:knon), petBcoef(:knon), peqBcoef(:knon), &
            tsurf_new, evap(:knon), fluxlat(:knon), fluxsens(:knon), &
            dflux_s(:knon), dflux_l(:knon))
       fder = fder + dflux_s + dflux_l

       !IM: flux ocean-atmosphere utile pour le "slab" ocean
       flux_o(:knon) = fluxsens(:knon) - evap(:knon) &
            * merge(RLSTT, RLVTT, tsurf_new < RTT)

       ! calcul albedo
       if (cycle_diurne) then
          CALL alboc_cd(rmu0(knindex), albedo)
       else
          CALL alboc(jour, rlat(knindex), albedo)
       endif

       z0_new = sqrt(rugos**2 + rugoro**2)
    case (is_sic)
       ! Surface "glace de mer" appel a l'interface avec l'ocean

       ! ! lecture conditions limites
       CALL interfoce_lim(itime, dtime, jour, knindex, debut, tsurf_new, &
            pctsrf_new)

       DO ii = 1, knon
          tsurf_new(ii) = tsurf(ii)
          IF (pctsrf_new(knindex(ii), nisurf) < EPSFRA) then
             snow(ii) = 0.0
             tsurf_new(ii) = RTT - 1.8
             IF (soil_model) tsoil(ii, :) = RTT - 1.8
          endif
       enddo

       CALL calbeta(nisurf, snow(:knon), qsol(:knon), beta(:knon), &
            capsol(:knon), dif_grnd(:knon))

       IF (soil_model) THEN
          CALL soil(dtime, nisurf, knon, snow, tsurf_new, tsoil, soilcap, &
               soilflux)
          cal(1:knon) = RCPD / soilcap(1:knon)
          radsol(1:knon) = radsol(1:knon) + soilflux(1:knon)
          dif_grnd = 0.
       ELSE
          dif_grnd = 1.0 / tau_gl
          cal = RCPD * calice
          WHERE (snow > 0.0) cal = RCPD * calsno
       ENDIF
       tsurf_temp = tsurf_new
       beta = 1.0

       CALL calcul_fluxs(dtime, tsurf_temp, p1lay(:knon), cal(:knon), &
            beta(:knon), tq_cdrag(:knon), ps(:knon), qsurf(:knon), &
            radsol(:knon), dif_grnd(:knon), temp_air(:knon), spechum(:knon), &
            u1_lay(:knon), v1_lay(:knon), petAcoef(:knon), peqAcoef(:knon), &
            petBcoef(:knon), peqBcoef(:knon), tsurf_new, evap(:knon), &
            fluxlat(:knon), fluxsens(:knon), dflux_s(:knon), dflux_l(:knon))

       !IM: flux entre l'ocean et la glace de mer pour le "slab" ocean
       DO i = 1, knon
          flux_g(i) = 0.0
          IF (cal(i) > 1e-15) flux_g(i) = (tsurf_new(i) - t_grnd) &
               * dif_grnd(i) * RCPD / cal(i)
       ENDDO

       CALL fonte_neige(nisurf, dtime, tsurf_temp, p1lay(:knon), beta(:knon), &
            tq_cdrag(:knon), ps(:knon), precip_rain(:knon), &
            precip_snow(:knon), snow(:knon), qsol(:knon), temp_air(:knon), &
            spechum(:knon), u1_lay(:knon), v1_lay(:knon), petAcoef(:knon), &
            peqAcoef(:knon), petBcoef(:knon), peqBcoef(:knon), tsurf_new, &
            evap(:knon), fqcalving(:knon), ffonte(:knon), run_off_lic_0(:knon))

       ! calcul albedo

       CALL albsno(klon, knon, dtime, agesno, alb_neig, precip_snow)
       WHERE (snow(1 : knon) < 0.0001) agesno(1 : knon) = 0.
       zfra(:knon) = MAX(0.0, MIN(1.0, snow(1:knon)/(snow(1:knon) + 10.0)))
       albedo = alb_neig(:knon) * zfra(:knon) + 0.6 * (1.0 - zfra(:knon))

       fder = fder + dflux_s + dflux_l

       ! 2eme appel a interfoce pour le cumul et le passage des flux a l'ocean

       z0_new = 0.002
       z0_new = SQRT(z0_new**2 + rugoro**2)
    case (is_lic)
       if (.not. allocated(run_off_lic)) then
          allocate(run_off_lic(knon))
          run_off_lic = 0.
       endif

       ! Surface "glacier continentaux" appel a l'interface avec le sol

       IF (soil_model) THEN
          CALL soil(dtime, nisurf, knon, snow, tsurf, tsoil, soilcap, soilflux)
          cal(1:knon) = RCPD / soilcap(1:knon)
          radsol(1:knon) = radsol(1:knon) + soilflux(1:knon)
       ELSE
          cal = RCPD * calice
          WHERE (snow > 0.0) cal = RCPD * calsno
       ENDIF
       beta = 1.0
       dif_grnd = 0.0

       call calcul_fluxs(dtime, tsurf, p1lay(:knon), cal(:knon), &
            beta(:knon), tq_cdrag(:knon), ps(:knon), qsurf(:knon), &
            radsol(:knon), dif_grnd(:knon), temp_air(:knon), spechum(:knon), &
            u1_lay(:knon), v1_lay(:knon), petAcoef(:knon), peqAcoef(:knon), &
            petBcoef(:knon), peqBcoef(:knon), tsurf_new, evap(:knon), &
            fluxlat(:knon), fluxsens(:knon), dflux_s(:knon), dflux_l(:knon))

       call fonte_neige(nisurf, dtime, tsurf, p1lay(:knon), beta(:knon), &
            tq_cdrag(:knon), ps(:knon), precip_rain(:knon), &
            precip_snow(:knon), snow(:knon), qsol(:knon), temp_air(:knon), &
            spechum(:knon), u1_lay(:knon), v1_lay(:knon), petAcoef(:knon), &
            peqAcoef(:knon), petBcoef(:knon), peqBcoef(:knon), tsurf_new, &
            evap(:knon), fqcalving(:knon), ffonte(:knon), run_off_lic_0(:knon))

       ! calcul albedo
       CALL albsno(klon, knon, dtime, agesno, alb_neig, precip_snow)
       WHERE (snow(1 : knon) < 0.0001) agesno(1 : knon) = 0.
       zfra(:knon) = MAX(0.0, MIN(1.0, snow(1:knon)/(snow(1:knon) + 10.0)))
       albedo = 0.77

       ! Rugosite
       z0_new = rugoro

       ! Remplissage des pourcentages de surface
       pctsrf_new(:, nisurf) = pctsrf(:, nisurf)

    case default
       print *, 'Index surface = ', nisurf
       abort_message = 'Index surface non valable'
       call abort_gcm(modname, abort_message)
    end select

  END SUBROUTINE interfsurf_hq

end module interfsurf_hq_m
