module interfsurf_hq_m

  implicit none

contains

  SUBROUTINE interfsurf_hq(itime, dtime, jour, rmu0, iim, jjm, nisurf, knon, &
       knindex, pctsrf, rlat, debut, soil_model, nsoilmx, tsoil, qsol, &
       u1_lay, v1_lay, temp_air, spechum, tq_cdrag, petAcoef, peqAcoef, &
       petBcoef, peqBcoef, precip_rain, precip_snow, fder, rugos, rugoro, &
       snow, qsurf, tsurf, p1lay, ps, radsol, evap, fluxsens, fluxlat, &
       dflux_l, dflux_s, tsurf_new, alb_new, alblw, z0_new, pctsrf_new, &
       agesno, fqcalving, ffonte, run_off_lic_0, flux_o, flux_g, tslab, seaice)

    ! Cette routine sert d'aiguillage entre l'atmosph�re et la surface
    ! en g�n�ral (sols continentaux, oc�ans, glaces) pour les flux de
    ! chaleur et d'humidit�.

    ! Laurent Fairhead, 02/2000

    USE abort_gcm_m, ONLY: abort_gcm
    USE albsno_m, ONLY: albsno
    USE calcul_fluxs_m, ONLY: calcul_fluxs
    USE dimphy, ONLY: klon
    USE fonte_neige_m, ONLY: fonte_neige
    USE indicesol, ONLY: epsfra, is_lic, is_oce, is_sic, is_ter, nbsrf
    USE interface_surf, ONLY: coastalflow, riverflow, run_off, &
         run_off_lic, conf_interface
    USE interfoce_lim_m, ONLY: interfoce_lim
    USE interfoce_slab_m, ONLY: interfoce_slab
    USE interfsur_lim_m, ONLY: interfsur_lim
    USE suphec_m, ONLY: rcpd, rlstt, rlvtt, rtt

    integer, intent(IN):: itime ! numero du pas de temps
    real, intent(IN):: dtime ! pas de temps de la physique (en s)
    integer, intent(IN):: jour ! jour dans l'annee en cours
    real, intent(IN):: rmu0(klon) ! cosinus de l'angle solaire zenithal
    integer, intent(IN):: iim, jjm
    ! iim, jjm nbres de pts de grille
    integer, intent(IN):: nisurf
    ! nisurf index de la surface a traiter (1 = sol continental)
    integer, intent(IN):: knon
    ! knon nombre de points de la surface a traiter
    integer, intent(in):: knindex(klon)
    ! knindex index des points de la surface a traiter
    real, intent(IN):: pctsrf(klon, nbsrf)
    ! pctsrf tableau des pourcentages de surface de chaque maille
    real, dimension(klon), intent(IN):: rlat
    ! rlat latitudes
    logical, intent(IN):: debut
    ! debut logical: 1er appel a la physique
    ! (si false calcul simplifie des fluxs sur les continents)
    !! PB ajout pour soil
    logical, intent(in):: soil_model
    integer:: nsoilmx
    REAL, DIMENSION(klon, nsoilmx):: tsoil
    REAL, dimension(klon), intent(INOUT):: qsol
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
    real, dimension(klon), intent(IN):: precip_rain, precip_snow
    ! precip_rain precipitation liquide
    ! precip_snow precipitation solide
    REAL, DIMENSION(klon), INTENT(INOUT):: fder
    ! fder derivee des flux (pour le couplage)
    real, dimension(klon), intent(IN):: rugos, rugoro
    ! rugos rugosite
    ! rugoro rugosite orographique
    real, dimension(klon), intent(INOUT):: snow, qsurf
    real, dimension(klon), intent(IN):: tsurf, p1lay
    ! tsurf temperature de surface
    ! p1lay pression 1er niveau (milieu de couche)
    real, dimension(klon), intent(IN):: ps
    ! ps pression au sol
    REAL, DIMENSION(klon), INTENT(INOUT):: radsol
    ! radsol rayonnement net aus sol (LW + SW)
    real, dimension(klon), intent(INOUT):: evap
    ! evap evaporation totale
    real, dimension(klon), intent(OUT):: fluxsens, fluxlat
    ! fluxsens flux de chaleur sensible
    ! fluxlat flux de chaleur latente
    real, dimension(klon), intent(OUT):: dflux_l, dflux_s
    real, dimension(klon), intent(OUT):: tsurf_new, alb_new
    ! tsurf_new temperature au sol
    ! alb_new albedo
    real, dimension(klon), intent(OUT):: alblw
    real, dimension(klon), intent(OUT):: z0_new
    ! z0_new surface roughness
    real, dimension(klon, nbsrf), intent(OUT):: pctsrf_new
    ! pctsrf_new nouvelle repartition des surfaces
    real, dimension(klon), intent(INOUT):: agesno

    ! Flux d'eau "perdue" par la surface et n�cessaire pour que limiter la
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
    real, dimension(klon), intent(INOUT):: tslab
    ! tslab temperature slab ocean
    real, dimension(klon), intent(INOUT):: seaice ! glace de mer (kg/m2)

    ! Local:

    real, allocatable, dimension(:), save:: tmp_tslab
    REAL, dimension(klon):: soilcap
    REAL, dimension(klon):: soilflux

    !IM: "slab" ocean
    real, parameter:: t_grnd=271.35
    real, dimension(klon):: zx_sl
    integer i
    real, allocatable, dimension(:), save:: tmp_flux_o, tmp_flux_g
    real, allocatable, dimension(:), save:: tmp_radsol
    real, allocatable, dimension(:, :), save:: tmp_pctsrf_slab
    ! pctsrf_slab pourcentages (0-1) des sous-surfaces dans le slab
    ! tmp_pctsrf_slab = pctsrf_slab
    real, allocatable, dimension(:), save:: tmp_seaice

    character (len = 20), save:: modname = 'interfsurf_hq'
    character (len = 80):: abort_message
    logical, save:: first_call = .true.
    integer, save:: error
    integer:: ii
    logical, save:: check = .false.
    real, dimension(klon):: cal, beta, dif_grnd, capsol
    real, parameter:: calice=1.0/(5.1444e+06*0.15), tau_gl=86400.*5.
    real, parameter:: calsno=1./(2.3867e+06*.15)
    real, dimension(klon):: tsurf_temp
    real, dimension(klon):: alb_neig, alb_eau
    real, DIMENSION(klon):: zfra
    INTEGER, dimension(1):: iloc
    real, dimension(klon):: fder_prev
    REAL, dimension(klon):: bidule

    !-------------------------------------------------------------

    if (check) write(*, *) 'Entree ', modname

    ! On doit commencer par appeler les schemas de surfaces continentales
    ! car l'ocean a besoin du ruissellement qui est y calcule

    if (first_call) then
       call conf_interface
       if (nisurf /= is_ter .and. klon > 1) then
          write(*, *)' *** Warning ***'
          write(*, *)' nisurf = ', nisurf, ' /= is_ter = ', is_ter
          write(*, *)'or on doit commencer par les surfaces continentales'
          abort_message='voir ci-dessus'
          call abort_gcm(modname, abort_message, 1)
       endif
       if ( is_oce > is_sic ) then
          write(*, *)' *** Warning ***'
          write(*, *)' Pour des raisons de sequencement dans le code'
          write(*, *)' l''ocean doit etre traite avant la banquise'
          write(*, *)' or is_oce = ', is_oce, '> is_sic = ', is_sic
          abort_message='voir ci-dessus'
          call abort_gcm(modname, abort_message, 1)
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
    alb_new = 999999.
    z0_new = 999999.
    alb_neig = 999999.
    tsurf_new = 999999.
    alblw = 999999.

    !IM: "slab" ocean; initialisations
    flux_o = 0.
    flux_g = 0.

    if (.not. allocated(tmp_flux_o)) then
       allocate(tmp_flux_o(klon), stat = error)
       DO i=1, knon
          tmp_flux_o(knindex(i))=flux_o(i)
       ENDDO
       if (error /= 0) then
          abort_message='Pb allocation tmp_flux_o'
          call abort_gcm(modname, abort_message, 1)
       endif
    endif
    if (.not. allocated(tmp_flux_g)) then
       allocate(tmp_flux_g(klon), stat = error)
       DO i=1, knon
          tmp_flux_g(knindex(i))=flux_g(i)
       ENDDO
       if (error /= 0) then
          abort_message='Pb allocation tmp_flux_g'
          call abort_gcm(modname, abort_message, 1)
       endif
    endif
    if (.not. allocated(tmp_radsol)) then
       allocate(tmp_radsol(klon), stat = error)
       if (error /= 0) then
          abort_message='Pb allocation tmp_radsol'
          call abort_gcm(modname, abort_message, 1)
       endif
    endif
    DO i=1, knon
       tmp_radsol(knindex(i))=radsol(i)
    ENDDO
    if (.not. allocated(tmp_pctsrf_slab)) then
       allocate(tmp_pctsrf_slab(klon, nbsrf), stat = error)
       if (error /= 0) then
          abort_message='Pb allocation tmp_pctsrf_slab'
          call abort_gcm(modname, abort_message, 1)
       endif
       DO i=1, klon
          tmp_pctsrf_slab(i, 1:nbsrf)=pctsrf(i, 1:nbsrf)
       ENDDO
    endif

    if (.not. allocated(tmp_seaice)) then
       allocate(tmp_seaice(klon), stat = error)
       if (error /= 0) then
          abort_message='Pb allocation tmp_seaice'
          call abort_gcm(modname, abort_message, 1)
       endif
       DO i=1, klon
          tmp_seaice(i)=seaice(i)
       ENDDO
    endif

    if (.not. allocated(tmp_tslab)) then
       allocate(tmp_tslab(klon), stat = error)
       if (error /= 0) then
          abort_message='Pb allocation tmp_tslab'
          call abort_gcm(modname, abort_message, 1)
       endif
    endif
    DO i=1, klon
       tmp_tslab(i)=tslab(i)
    ENDDO

    ! Aiguillage vers les differents schemas de surface

    if (nisurf == is_ter) then
       ! Surface "terre" appel a l'interface avec les sols continentaux

       ! allocation du run-off
       if (.not. allocated(coastalflow)) then
          allocate(coastalflow(knon), stat = error)
          if (error /= 0) then
             abort_message='Pb allocation coastalflow'
             call abort_gcm(modname, abort_message, 1)
          endif
          allocate(riverflow(knon), stat = error)
          if (error /= 0) then
             abort_message='Pb allocation riverflow'
             call abort_gcm(modname, abort_message, 1)
          endif
          allocate(run_off(knon), stat = error)
          if (error /= 0) then
             abort_message='Pb allocation run_off'
             call abort_gcm(modname, abort_message, 1)
          endif

          run_off=0.0
       else if (size(coastalflow) /= knon) then
          write(*, *)'Bizarre, le nombre de points continentaux'
          write(*, *)'a change entre deux appels. J''arrete ...'
          abort_message='voir ci-dessus'
          call abort_gcm(modname, abort_message, 1)
       endif
       coastalflow = 0.
       riverflow = 0.

       ! Calcul age de la neige

       ! calcul albedo: lecture albedo fichier boundary conditions
       ! puis ajout albedo neige
       call interfsur_lim(itime, dtime, jour, nisurf, knon, knindex, &
            debut, alb_new, z0_new)

       ! calcul snow et qsurf, hydrol adapt�
       CALL calbeta(dtime, nisurf, knon, snow, qsol, beta, capsol, dif_grnd)

       IF (soil_model) THEN
          CALL soil(dtime, nisurf, knon, snow, tsurf, tsoil, soilcap, &
               soilflux)
          cal(1:knon) = RCPD / soilcap(1:knon)
          radsol(1:knon) = radsol(1:knon) + soilflux(1:knon)
       ELSE
          cal = RCPD * capsol
       ENDIF
       CALL calcul_fluxs( klon, knon, nisurf, dtime, &
            tsurf, p1lay, cal, beta, tq_cdrag, ps, &
            precip_rain, precip_snow, snow, qsurf, &
            radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, &
            petAcoef, peqAcoef, petBcoef, peqBcoef, &
            tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l)

       CALL fonte_neige( klon, knon, nisurf, dtime, &
            tsurf, p1lay, cal, beta, tq_cdrag, ps, &
            precip_rain, precip_snow, snow, qsol, &
            radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, &
            petAcoef, peqAcoef, petBcoef, peqBcoef, &
            tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l, &
            fqcalving, ffonte, run_off_lic_0)

       call albsno(klon, knon, dtime, agesno, alb_neig, precip_snow)
       where (snow(1 : knon) .LT. 0.0001) agesno(1 : knon) = 0.
       zfra(1:knon) = max(0.0, min(1.0, snow(1:knon)/(snow(1:knon)+10.0)))
       alb_new(1 : knon) = alb_neig(1 : knon) *zfra(1:knon) + &
            alb_new(1 : knon)*(1.0-zfra(1:knon))
       z0_new = sqrt(z0_new**2+rugoro**2)
       alblw(1 : knon) = alb_new(1 : knon)

       ! Remplissage des pourcentages de surface
       pctsrf_new(:, nisurf) = pctsrf(:, nisurf)
    else if (nisurf == is_oce) then
       ! Surface "ocean" appel a l'interface avec l'ocean
       ! lecture conditions limites
       call interfoce_lim(itime, dtime, jour, klon, nisurf, knon, knindex, &
            debut, tsurf_new, pctsrf_new)

       tsurf_temp = tsurf_new
       cal = 0.
       beta = 1.
       dif_grnd = 0.
       alb_neig = 0.
       agesno = 0.

       call calcul_fluxs( klon, knon, nisurf, dtime, &
            tsurf_temp, p1lay, cal, beta, tq_cdrag, ps, &
            precip_rain, precip_snow, snow, qsurf, &
            radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, &
            petAcoef, peqAcoef, petBcoef, peqBcoef, &
            tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l)

       fder_prev = fder
       fder = fder_prev + dflux_s + dflux_l

       iloc = maxloc(fder(1:klon))
       if (check.and.fder(iloc(1))> 0.) then
          WRITE(*, *)'**** Debug fder****'
          WRITE(*, *)'max fder(', iloc(1), ') = ', fder(iloc(1))
          WRITE(*, *)'fder_prev, dflux_s, dflux_l', fder_prev(iloc(1)), &
               dflux_s(iloc(1)), dflux_l(iloc(1))
       endif

       !IM: flux ocean-atmosphere utile pour le "slab" ocean
       DO i=1, knon
          zx_sl(i) = RLVTT
          if (tsurf_new(i) .LT. RTT) zx_sl(i) = RLSTT
          flux_o(i) = fluxsens(i)-evap(i)*zx_sl(i)
          tmp_flux_o(knindex(i)) = flux_o(i)
          tmp_radsol(knindex(i))=radsol(i)
       ENDDO

       ! calcul albedo
       if ( minval(rmu0) == maxval(rmu0) .and. minval(rmu0) == -999.999 ) then
          CALL alboc(FLOAT(jour), rlat, alb_eau)
       else ! cycle diurne
          CALL alboc_cd(rmu0, alb_eau)
       endif
       DO ii =1, knon
          alb_new(ii) = alb_eau(knindex(ii))
       enddo

       z0_new = sqrt(rugos**2 + rugoro**2)
       alblw(1:knon) = alb_new(1:knon)
    else if (nisurf == is_sic) then
       if (check) write(*, *)'sea ice, nisurf = ', nisurf

       ! Surface "glace de mer" appel a l'interface avec l'ocean

       ! ! lecture conditions limites
       CALL interfoce_lim(itime, dtime, jour, &
            klon, nisurf, knon, knindex, &
            debut, &
            tsurf_new, pctsrf_new)

       !IM cf LF
       DO ii = 1, knon
          tsurf_new(ii) = tsurf(ii)
          !IMbad IF (pctsrf_new(ii, nisurf) < EPSFRA) then
          IF (pctsrf_new(knindex(ii), nisurf) < EPSFRA) then
             snow(ii) = 0.0
             !IM cf LF/JLD tsurf(ii) = RTT - 1.8
             tsurf_new(ii) = RTT - 1.8
             IF (soil_model) tsoil(ii, :) = RTT -1.8
          endif
       enddo

       CALL calbeta(dtime, nisurf, knon, snow, qsol, beta, capsol, dif_grnd)

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
       !IMbadtsurf_temp = tsurf
       tsurf_temp = tsurf_new
       beta = 1.0

       CALL calcul_fluxs( klon, knon, nisurf, dtime, &
            tsurf_temp, p1lay, cal, beta, tq_cdrag, ps, &
            precip_rain, precip_snow, snow, qsurf, &
            radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, &
            petAcoef, peqAcoef, petBcoef, peqBcoef, &
            tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l)

       !IM: flux entre l'ocean et la glace de mer pour le "slab" ocean
       DO i = 1, knon
          flux_g(i) = 0.0

          !IM: faire dependre le coefficient de conduction de la glace de mer
          ! de l'epaisseur de la glace de mer, dans l'hypothese ou le coeff.
          ! actuel correspond a 3m de glace de mer, cf. L.Li

          ! IF(1.EQ.0) THEN
          ! IF(siceh(i).GT.0.) THEN
          ! new_dif_grnd(i) = dif_grnd(i)*3./siceh(i)
          ! ELSE
          ! new_dif_grnd(i) = 0.
          ! ENDIF
          ! ENDIF !(1.EQ.0) THEN

          IF (cal(i).GT.1.0e-15) flux_g(i)=(tsurf_new(i)-t_grnd) &
               * dif_grnd(i) *RCPD/cal(i)
          ! & * new_dif_grnd(i) *RCPD/cal(i)
          tmp_flux_g(knindex(i))=flux_g(i)
          tmp_radsol(knindex(i))=radsol(i)
       ENDDO

       CALL fonte_neige( klon, knon, nisurf, dtime, &
            tsurf_temp, p1lay, cal, beta, tq_cdrag, ps, &
            precip_rain, precip_snow, snow, qsol, &
            radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, &
            petAcoef, peqAcoef, petBcoef, peqBcoef, &
            tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l, &
            fqcalving, ffonte, run_off_lic_0)

       ! calcul albedo

       CALL albsno(klon, knon, dtime, agesno, alb_neig, precip_snow)
       WHERE (snow(1 : knon) .LT. 0.0001) agesno(1 : knon) = 0.
       zfra(1:knon) = MAX(0.0, MIN(1.0, snow(1:knon)/(snow(1:knon)+10.0)))
       alb_new(1 : knon) = alb_neig(1 : knon) *zfra(1:knon) + &
            0.6 * (1.0-zfra(1:knon))

       fder_prev = fder
       fder = fder_prev + dflux_s + dflux_l

       iloc = maxloc(fder(1:klon))
       if (check.and.fder(iloc(1))> 0.) then
          WRITE(*, *)'**** Debug fder ****'
          WRITE(*, *)'max fder(', iloc(1), ') = ', fder(iloc(1))
          WRITE(*, *)'fder_prev, dflux_s, dflux_l', fder_prev(iloc(1)), &
               dflux_s(iloc(1)), dflux_l(iloc(1))
       endif


       ! 2eme appel a interfoce pour le cumul et le passage des flux a l'ocean

       z0_new = 0.002
       z0_new = SQRT(z0_new**2+rugoro**2)
       alblw(1:knon) = alb_new(1:knon)

    else if (nisurf == is_lic) then

       if (check) write(*, *)'glacier, nisurf = ', nisurf

       if (.not. allocated(run_off_lic)) then
          allocate(run_off_lic(knon), stat = error)
          if (error /= 0) then
             abort_message='Pb allocation run_off_lic'
             call abort_gcm(modname, abort_message, 1)
          endif
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

       call calcul_fluxs( klon, knon, nisurf, dtime, &
            tsurf, p1lay, cal, beta, tq_cdrag, ps, &
            precip_rain, precip_snow, snow, qsurf, &
            radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, &
            petAcoef, peqAcoef, petBcoef, peqBcoef, &
            tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l)

       call fonte_neige( klon, knon, nisurf, dtime, &
            tsurf, p1lay, cal, beta, tq_cdrag, ps, &
            precip_rain, precip_snow, snow, qsol, &
            radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, &
            petAcoef, peqAcoef, petBcoef, peqBcoef, &
            tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l, &
            fqcalving, ffonte, run_off_lic_0)

       ! passage du run-off des glaciers calcule dans fonte_neige au coupleur
       bidule=0.
       bidule(1:knon)= run_off_lic(1:knon)

       ! calcul albedo

       CALL albsno(klon, knon, dtime, agesno, alb_neig, precip_snow)
       WHERE (snow(1 : knon) .LT. 0.0001) agesno(1 : knon) = 0.
       zfra(1:knon) = MAX(0.0, MIN(1.0, snow(1:knon)/(snow(1:knon)+10.0)))
       alb_new(1 : knon) = alb_neig(1 : knon)*zfra(1:knon) + &
            0.6 * (1.0-zfra(1:knon))

       !IM: plusieurs choix/tests sur l'albedo des "glaciers continentaux"
       ! alb_new(1 : knon) = 0.6 !IM cf FH/GK
       ! alb_new(1 : knon) = 0.82
       ! alb_new(1 : knon) = 0.77 !211003 Ksta0.77
       ! alb_new(1 : knon) = 0.8 !KstaTER0.8 & LMD_ARMIP5
       !IM: KstaTER0.77 & LMD_ARMIP6
       alb_new(1 : knon) = 0.77


       ! Rugosite

       z0_new = rugoro

       ! Remplissage des pourcentages de surface

       pctsrf_new(:, nisurf) = pctsrf(:, nisurf)

       alblw(1:knon) = alb_new(1:knon)
    else
       write(*, *)'Index surface = ', nisurf
       abort_message = 'Index surface non valable'
       call abort_gcm(modname, abort_message, 1)
    endif

  END SUBROUTINE interfsurf_hq

end module interfsurf_hq_m
