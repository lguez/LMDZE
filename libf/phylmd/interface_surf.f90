MODULE interface_surf

  ! From phylmd/interface_surf.F90,v 1.8 2005/05/25 13:10:09

  ! Ce module regroupe toutes les routines gérant l'interface entre le modèle 
  ! atmosphérique et les modèles de surface (sols continentaux,
  ! océans, glaces).
  ! Les routines sont les suivantes:
  ! interfsurf_hq : routine d'aiguillage vers les interfaces avec les
  !    différents modèles de surface
  ! interfoce_* : routines d'interface proprement dites

  ! L. Fairhead, LMD, 02/2000

  IMPLICIT none

  PRIVATE
  PUBLIC :: interfsurf_hq

  ! run_off      ruissellement total
  REAL, ALLOCATABLE, DIMENSION(:),SAVE    :: run_off, run_off_lic
  real, allocatable, dimension(:),save    :: coastalflow, riverflow

  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE :: tmp_rriv, tmp_rcoa,tmp_rlic
  !! pour simuler la fonte des glaciers antarctiques
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE :: coeff_iceberg
  real, save                              :: surf_maille 
  real, save                              :: cte_flux_iceberg = 6.3e7
  integer, save                           :: num_antarctic = 1
  REAL, save                              :: tau_calv

CONTAINS

  SUBROUTINE interfsurf_hq(itime, dtime, date0, jour, rmu0, &
       & klon, iim, jjm, nisurf, knon, knindex, pctsrf, &
       & rlon, rlat, cufi, cvfi,&
       & debut, lafin, ok_veget, soil_model, nsoilmx, tsoil, qsol,&
       & zlev,  u1_lay, v1_lay, temp_air, spechum, epot_air, ccanopy, & 
       & tq_cdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
       & precip_rain, precip_snow, sollw, sollwdown, swnet, swdown, &
       & fder, taux, tauy, &
       & windsp, &
       & rugos, rugoro, &
       & albedo, snow, qsurf, &
       & tsurf, p1lay, ps, radsol, &
       & ocean, npas, nexca, zmasq, &
       & evap, fluxsens, fluxlat, dflux_l, dflux_s, &              
       & tsol_rad, tsurf_new, alb_new, alblw, emis_new, &
       & z0_new, pctsrf_new, agesno,fqcalving,ffonte, run_off_lic_0,&
       !IM "slab" ocean
       & flux_o, flux_g, tslab, seaice)

    ! Cette routine sert d'aiguillage entre l'atmosphère et la surface
    ! en général (sols continentaux, océans, glaces) pour les fluxs de
    ! chaleur et d'humidité.
    ! En pratique l'interface se fait entre la couche limite du modèle 
    ! atmosphérique ("clmain.F") et les routines de surface
    ! ("sechiba", "oasis"...).

    ! L.Fairhead 02/2000

    ! input:
    !   klon         nombre total de points de grille
    !   iim, jjm     nbres de pts de grille
    !   dtime        pas de temps de la physique (en s)
    !   date0        jour initial 
    !   jour         jour dans l'annee en cours,
    !   rmu0         cosinus de l'angle solaire zenithal
    !   nexca        pas de temps couplage
    !   nisurf       index de la surface a traiter (1 = sol continental)
    !   knon         nombre de points de la surface a traiter
    !   knindex      index des points de la surface a traiter
    !   pctsrf       tableau des pourcentages de surface de chaque maille
    !   rlon         longitudes
    !   rlat         latitudes
    !   cufi,cvfi    resolution des mailles en x et y (m)
    !   debut        logical: 1er appel a la physique
    !   lafin        logical: dernier appel a la physique
    !   ok_veget     logical: appel ou non au schema de surface continental
    !                  (si false calcul simplifie des fluxs sur les continents)
    !   zlev         hauteur de la premiere couche
    !   u1_lay       vitesse u 1ere couche
    !   v1_lay       vitesse v 1ere couche
    !   temp_air     temperature de l'air 1ere couche
    !   spechum      humidite specifique 1ere couche
    !   epot_air     temp potentielle de l'air
    !   ccanopy      concentration CO2 canopee
    !   tq_cdrag     cdrag
    !   petAcoef     coeff. A de la resolution de la CL pour t
    !   peqAcoef     coeff. A de la resolution de la CL pour q
    !   petBcoef     coeff. B de la resolution de la CL pour t
    !   peqBcoef     coeff. B de la resolution de la CL pour q
    !   precip_rain  precipitation liquide
    !   precip_snow  precipitation solide
    !   sollw        flux IR net a la surface
    !   sollwdown    flux IR descendant a la surface
    !   swnet        flux solaire net
    !   swdown       flux solaire entrant a la surface
    !   albedo       albedo de la surface
    !   tsurf        temperature de surface
    !   tslab        temperature slab ocean
    !   pctsrf_slab  pourcentages (0-1) des sous-surfaces dans le slab
    !   tmp_pctsrf_slab = pctsrf_slab
    !   p1lay        pression 1er niveau (milieu de couche)
    !   ps           pression au sol
    !   radsol       rayonnement net aus sol (LW + SW)
    !   ocean        type d'ocean utilise (force, slab, couple)
    !   fder         derivee des flux (pour le couplage)
    !   taux, tauy   tension de vents
    !   windsp       module du vent a 10m
    !   rugos        rugosite
    !   zmasq        masque terre/ocean
    !   rugoro       rugosite orographique
    !   run_off_lic_0 runoff glacier du pas de temps precedent

    ! output:
    !   evap         evaporation totale
    !   fluxsens     flux de chaleur sensible
    !   fluxlat      flux de chaleur latente
    !   tsol_rad     
    !   tsurf_new    temperature au sol
    !   alb_new      albedo
    !   emis_new     emissivite
    !   z0_new       surface roughness
    !   pctsrf_new   nouvelle repartition des surfaces

    use abort_gcm_m, only: abort_gcm
    use gath_cpl, only: gath2cpl
    use indicesol
    use YOMCST
    use albsno_m, only: albsno

    ! Parametres d'entree
    integer, intent(IN) :: itime     !  numero du pas de temps
    integer, intent(IN) :: iim, jjm
    integer, intent(IN) :: klon
    real, intent(IN) :: dtime
    real, intent(IN) :: date0
    integer, intent(IN) :: jour
    real, intent(IN)    :: rmu0(klon)
    integer, intent(IN) :: nisurf
    integer, intent(IN) :: knon
    integer, dimension(klon), intent(in) :: knindex
    real, dimension(klon,nbsrf), intent(IN) :: pctsrf
    logical, intent(IN) :: debut, lafin, ok_veget
    real, dimension(klon), intent(IN) :: rlon, rlat
    real, dimension(klon), intent(IN) :: cufi, cvfi
    real, dimension(klon), intent(INOUT) :: tq_cdrag
    real, dimension(klon), intent(IN) :: zlev
    real, dimension(klon), intent(IN) :: u1_lay, v1_lay
    real, dimension(klon), intent(IN) :: temp_air, spechum
    real, dimension(klon), intent(IN) :: epot_air, ccanopy
    real, dimension(klon), intent(IN) :: petAcoef, peqAcoef
    real, dimension(klon), intent(IN) :: petBcoef, peqBcoef
    real, dimension(klon), intent(IN) :: precip_rain, precip_snow
    real, dimension(klon), intent(IN) :: sollw, sollwdown, swnet, swdown
    real, dimension(klon), intent(IN) :: ps, albedo
    real, dimension(klon), intent(IN) :: tsurf, p1lay
    !IM: "slab" ocean
    real, dimension(klon), intent(INOUT) :: tslab
    real, allocatable, dimension(:), save :: tmp_tslab
    real, dimension(klon), intent(OUT) :: flux_o, flux_g
    real, dimension(klon), intent(INOUT)      :: seaice ! glace de mer (kg/m2)
    REAL, DIMENSION(klon), INTENT(INOUT) :: radsol,fder
    real, dimension(klon), intent(IN) :: zmasq
    real, dimension(klon), intent(IN) :: taux, tauy, rugos, rugoro
    real, dimension(klon), intent(IN) :: windsp
    character (len = 6)  :: ocean
    integer              :: npas, nexca ! nombre et pas de temps couplage
    real, dimension(klon), intent(INOUT) :: evap, snow, qsurf
    !! PB ajout pour soil
    logical          :: soil_model
    integer          :: nsoilmx
    REAL, DIMENSION(klon, nsoilmx) :: tsoil
    REAL, dimension(klon), intent(INOUT) :: qsol
    REAL, dimension(klon)          :: soilcap
    REAL, dimension(klon)          :: soilflux
    ! Parametres de sortie
    real, dimension(klon), intent(OUT):: fluxsens, fluxlat
    real, dimension(klon), intent(OUT):: tsol_rad, tsurf_new, alb_new
    real, dimension(klon), intent(OUT):: alblw
    real, dimension(klon), intent(OUT):: emis_new, z0_new
    real, dimension(klon), intent(OUT):: dflux_l, dflux_s
    real, dimension(klon,nbsrf), intent(OUT) :: pctsrf_new
    real, dimension(klon), intent(INOUT):: agesno
    real, dimension(klon), intent(INOUT):: run_off_lic_0

    ! Flux thermique utiliser pour fondre la neige
    !jld a rajouter   real, dimension(klon), intent(INOUT):: ffonte
    real, dimension(klon), intent(INOUT):: ffonte
    ! Flux d'eau "perdue" par la surface et nécessaire pour que limiter la
    ! hauteur de neige, en kg/m2/s
    !jld a rajouter   real, dimension(klon), intent(INOUT):: fqcalving
    real, dimension(klon), intent(INOUT):: fqcalving
    !IM: "slab" ocean - Local
    real, parameter :: t_grnd=271.35
    real, dimension(klon) :: zx_sl
    integer i
    real, allocatable, dimension(:), save :: tmp_flux_o, tmp_flux_g
    real, allocatable, dimension(:), save :: tmp_radsol
    real, allocatable, dimension(:,:), save :: tmp_pctsrf_slab
    real, allocatable, dimension(:), save :: tmp_seaice

    ! Local
    character (len = 20),save :: modname = 'interfsurf_hq'
    character (len = 80) :: abort_message 
    logical, save        :: first_call = .true.
    integer, save        :: error
    integer              :: ii
    logical,save              :: check = .false.
    real, dimension(klon):: cal, beta, dif_grnd, capsol
!!$PB  real, parameter      :: calice=1.0/(5.1444e+06*0.15), tau_gl=86400.*5.
    real, parameter      :: calice=1.0/(5.1444e+06*0.15), tau_gl=86400.*5.
    real, parameter      :: calsno=1./(2.3867e+06*.15)
    real, dimension(klon):: tsurf_temp
    real, dimension(klon):: alb_neig, alb_eau
    real, DIMENSION(klon):: zfra
    logical              :: cumul = .false.
    INTEGER,dimension(1) :: iloc
    real, dimension(klon):: fder_prev
    REAL, dimension(klon) :: bidule

    !-------------------------------------------------------------

    if (check) write(*,*) 'Entree ', modname

    ! On doit commencer par appeler les schemas de surfaces continentales
    ! car l'ocean a besoin du ruissellement qui est y calcule

    if (first_call) then
       call conf_interface(tau_calv)
       if (nisurf /= is_ter .and. klon > 1) then 
          write(*,*)' *** Warning ***'
          write(*,*)' nisurf = ',nisurf,' /= is_ter = ',is_ter
          write(*,*)'or on doit commencer par les surfaces continentales'
          abort_message='voir ci-dessus'
          call abort_gcm(modname,abort_message,1)
       endif
       if (ocean /= 'slab' .and. ocean /= 'force' .and. ocean /= 'couple') then
          write(*,*)' *** Warning ***'
          write(*,*)'Option couplage pour l''ocean = ', ocean
          abort_message='option pour l''ocean non valable'
          call abort_gcm(modname,abort_message,1)
       endif
       if ( is_oce > is_sic ) then
          write(*,*)' *** Warning ***'
          write(*,*)' Pour des raisons de sequencement dans le code'
          write(*,*)' l''ocean doit etre traite avant la banquise'
          write(*,*)' or is_oce = ',is_oce, '> is_sic = ',is_sic
          abort_message='voir ci-dessus'
          call abort_gcm(modname,abort_message,1)
       endif
    endif
    first_call = .false.

    ! Initialisations diverses
    !
    ffonte(1:knon)=0.
    fqcalving(1:knon)=0.

    cal = 999999. ; beta = 999999. ; dif_grnd = 999999. ; capsol = 999999.
    alb_new = 999999. ; z0_new = 999999. ; alb_neig = 999999.
    tsurf_new = 999999.
    alblw = 999999.

    !IM: "slab" ocean; initialisations
    flux_o = 0.
    flux_g = 0.
    !
    if (.not. allocated(tmp_flux_o)) then
       allocate(tmp_flux_o(klon), stat = error)
       DO i=1, knon
          tmp_flux_o(knindex(i))=flux_o(i)
       ENDDO
       if (error /= 0) then
          abort_message='Pb allocation tmp_flux_o'
          call abort_gcm(modname,abort_message,1)
       endif
    endif
    if (.not. allocated(tmp_flux_g)) then  
       allocate(tmp_flux_g(klon), stat = error)
       DO i=1, knon
          tmp_flux_g(knindex(i))=flux_g(i)
       ENDDO
       if (error /= 0) then
          abort_message='Pb allocation tmp_flux_g'
          call abort_gcm(modname,abort_message,1)
       endif
    endif
    if (.not. allocated(tmp_radsol)) then  
       allocate(tmp_radsol(klon), stat = error)
       if (error /= 0) then
          abort_message='Pb allocation tmp_radsol'
          call abort_gcm(modname,abort_message,1)
       endif
    endif
    DO i=1, knon
       tmp_radsol(knindex(i))=radsol(i)
    ENDDO
    if (.not. allocated(tmp_pctsrf_slab)) then
       allocate(tmp_pctsrf_slab(klon,nbsrf), stat = error)
       if (error /= 0) then
          abort_message='Pb allocation tmp_pctsrf_slab'
          call abort_gcm(modname,abort_message,1)
       endif
       DO i=1, klon
          tmp_pctsrf_slab(i,1:nbsrf)=pctsrf(i,1:nbsrf)
       ENDDO
    endif
    !
    if (.not. allocated(tmp_seaice)) then
       allocate(tmp_seaice(klon), stat = error)
       if (error /= 0) then
          abort_message='Pb allocation tmp_seaice'
          call abort_gcm(modname,abort_message,1)
       endif
       DO i=1, klon
          tmp_seaice(i)=seaice(i)
       ENDDO
    endif
    !
    if (.not. allocated(tmp_tslab)) then
       allocate(tmp_tslab(klon), stat = error)
       if (error /= 0) then
          abort_message='Pb allocation tmp_tslab'
          call abort_gcm(modname,abort_message,1)
       endif
    endif
    DO i=1, klon
       tmp_tslab(i)=tslab(i)
    ENDDO
    !
    ! Aiguillage vers les differents schemas de surface

    if (nisurf == is_ter) then
       !
       ! Surface "terre" appel a l'interface avec les sols continentaux
       !
       ! allocation du run-off
       if (.not. allocated(coastalflow)) then
          allocate(coastalflow(knon), stat = error)
          if (error /= 0) then
             abort_message='Pb allocation coastalflow'
             call abort_gcm(modname,abort_message,1)
          endif
          allocate(riverflow(knon), stat = error)
          if (error /= 0) then
             abort_message='Pb allocation riverflow'
             call abort_gcm(modname,abort_message,1)
          endif
          allocate(run_off(knon), stat = error)
          if (error /= 0) then
             abort_message='Pb allocation run_off'
             call abort_gcm(modname,abort_message,1)
          endif
          !cym      
          run_off=0.0
          !cym

!!$PB
          ALLOCATE (tmp_rriv(iim,jjm+1), stat=error)
          if (error /= 0) then
             abort_message='Pb allocation tmp_rriv'
             call abort_gcm(modname,abort_message,1)
          endif
          ALLOCATE (tmp_rcoa(iim,jjm+1), stat=error)
          if (error /= 0) then
             abort_message='Pb allocation tmp_rcoa'
             call abort_gcm(modname,abort_message,1)
          endif
          ALLOCATE (tmp_rlic(iim,jjm+1), stat=error)
          if (error /= 0) then
             abort_message='Pb allocation tmp_rlic'
             call abort_gcm(modname,abort_message,1)
          endif
          tmp_rriv = 0.0
          tmp_rcoa = 0.0
          tmp_rlic = 0.0

!!$
       else if (size(coastalflow) /= knon) then
          write(*,*)'Bizarre, le nombre de points continentaux'
          write(*,*)'a change entre deux appels. J''arrete ...'
          abort_message='voir ci-dessus'
          call abort_gcm(modname,abort_message,1)
       endif
       coastalflow = 0.
       riverflow = 0.
       !
       ! Calcul age de la neige

       if (.not. ok_veget) then
          !
          ! calcul albedo: lecture albedo fichier CL puis ajout albedo neige 
          ! 
          call interfsur_lim(itime, dtime, jour, &
               & klon, nisurf, knon, knindex, debut,  &
               & alb_new, z0_new)
          !  
          ! calcul snow et qsurf, hydrol adapté
          !
          CALL calbeta(dtime, nisurf, knon, snow, qsol, beta, capsol, dif_grnd)

          IF (soil_model) THEN 
             CALL soil(dtime, nisurf, knon,snow, tsurf, tsoil,soilcap, soilflux)
             cal(1:knon) = RCPD / soilcap(1:knon)
             radsol(1:knon) = radsol(1:knon)  + soilflux(1:knon)
          ELSE 
             cal = RCPD * capsol
!!$      cal = capsol
          ENDIF
          CALL calcul_fluxs( klon, knon, nisurf, dtime, &
               &   tsurf, p1lay, cal, beta, tq_cdrag, ps, &
               &   precip_rain, precip_snow, snow, qsurf,  &
               &   radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, &
               &   petAcoef, peqAcoef, petBcoef, peqBcoef, &
               &   tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l)

          CALL fonte_neige( klon, knon, nisurf, dtime, &
               &   tsurf, p1lay, cal, beta, tq_cdrag, ps, &
               &   precip_rain, precip_snow, snow, qsol,  &
               &   radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, &
               &   petAcoef, peqAcoef, petBcoef, peqBcoef, &
               &   tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l, &
               &   fqcalving,ffonte, run_off_lic_0)

          call albsno(klon,knon,dtime,agesno(:),alb_neig(:), precip_snow(:))  
          where (snow(1 : knon) .LT. 0.0001) agesno(1 : knon) = 0.
          zfra(1:knon) = max(0.0,min(1.0,snow(1:knon)/(snow(1:knon)+10.0)))
          alb_new(1 : knon)  = alb_neig(1 : knon) *zfra(1:knon) + &
               &                     alb_new(1 : knon)*(1.0-zfra(1:knon))
          z0_new = sqrt(z0_new**2+rugoro**2)
          alblw(1 : knon) = alb_new(1 : knon)

       else
       endif
       !
       ! Remplissage des pourcentages de surface
       !
       pctsrf_new(:,nisurf) = pctsrf(:,nisurf)

    else if (nisurf == is_oce) then

       if (check) write(*,*)'ocean, nisurf = ',nisurf 

       !
       ! Surface "ocean" appel a l'interface avec l'ocean
       !
       if (ocean == 'couple') then
          if (nexca == 0) then
             abort_message='nexca = 0 dans interfoce_cpl'
             call abort_gcm(modname,abort_message,1)
          endif

          cumul = .false.

          iloc = maxloc(fder(1:klon))
          if (check) then
             if (fder(iloc(1))> 0.) then
                WRITE(*,*)'**** Debug fder ****'
                WRITE(*,*)'max fder(',iloc(1),') = ',fder(iloc(1))
             endif
          endif
!!$
!!$      where(fder.gt.0.) 
!!$        fder = 0.
!!$      endwhere

          call interfoce_cpl(itime, dtime, cumul, &
               & klon, iim, jjm, nisurf, pctsrf, knon, knindex, rlon, rlat, &
               & ocean, npas, nexca, debut, lafin, &
               & swdown, sollw, precip_rain, precip_snow, evap, tsurf, &
               & fluxlat, fluxsens, fder, albedo, taux, tauy, &
               & windsp, &
               & zmasq, &
               & tsurf_new, alb_new, &
               & pctsrf_new)

          !IM: "slab" ocean
       else if (ocean == 'slab  ') then
          tsurf_new = tsurf
          pctsrf_new = tmp_pctsrf_slab
          !
       else                              ! lecture conditions limites
          call interfoce_lim(itime, dtime, jour, & 
               &  klon, nisurf, knon, knindex, &
               &  debut, &
               &  tsurf_new, pctsrf_new)

       endif

       tsurf_temp = tsurf_new
       cal = 0.
       beta = 1.
       dif_grnd = 0.
       alb_neig(:) = 0.
       agesno(:) = 0.

       call calcul_fluxs( klon, knon, nisurf, dtime, &
            &   tsurf_temp, p1lay, cal, beta, tq_cdrag, ps, &
            &   precip_rain, precip_snow, snow, qsurf,  &
            &   radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, &
            &   petAcoef, peqAcoef, petBcoef, peqBcoef, &
            &   tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l)

       fder_prev = fder    
       fder = fder_prev + dflux_s + dflux_l

       iloc = maxloc(fder(1:klon))
       if (check.and.fder(iloc(1))> 0.) then
          WRITE(*,*)'**** Debug fder****'
          WRITE(*,*)'max fder(',iloc(1),') = ',fder(iloc(1))
          WRITE(*,*)'fder_prev, dflux_s, dflux_l',fder_prev(iloc(1)), &
               &                        dflux_s(iloc(1)), dflux_l(iloc(1))
       endif
!!$
!!$      where(fder.gt.0.) 
!!$        fder = 0.
!!$      endwhere

       !IM: flux ocean-atmosphere utile pour le "slab" ocean
       DO i=1, knon
          zx_sl(i) = RLVTT
          if (tsurf_new(i) .LT. RTT) zx_sl(i) = RLSTT
          flux_o(i) = fluxsens(i)-evap(i)*zx_sl(i)
          tmp_flux_o(knindex(i)) = flux_o(i)
          tmp_radsol(knindex(i))=radsol(i)
       ENDDO
       !
       ! 2eme appel a interfoce pour le cumul des champs (en particulier
       ! fluxsens et fluxlat calcules dans calcul_fluxs)
       !
       if (ocean == 'couple') then

          cumul = .true.

          call interfoce_cpl(itime, dtime, cumul, &
               & klon, iim, jjm, nisurf, pctsrf, knon, knindex, rlon, rlat, &
               & ocean, npas, nexca, debut, lafin, &
               & swdown, sollw, precip_rain, precip_snow, evap, tsurf, &
               & fluxlat, fluxsens, fder, albedo, taux, tauy, &
               & windsp, &
               & zmasq, &
               & tsurf_new, alb_new, &
               & pctsrf_new)

          !IM: "slab" ocean 
       else if (ocean == 'slab  ') then
          !
          seaice=tmp_seaice
          cumul = .true.
          call interfoce_slab(klon, debut, itime, dtime, jour, &
               & tmp_radsol, tmp_flux_o, tmp_flux_g, pctsrf, &
               & tslab, seaice, pctsrf_new)
          !
          tmp_pctsrf_slab=pctsrf_new
          DO i=1, knon
             tsurf_new(i)=tslab(knindex(i))
          ENDDO !i
          !
       endif

       !
       ! calcul albedo
       !

       if ( minval(rmu0) == maxval(rmu0) .and. minval(rmu0) == -999.999 ) then
          CALL alboc(FLOAT(jour),rlat,alb_eau)
       else  ! cycle diurne
          CALL alboc_cd(rmu0,alb_eau)
       endif
       DO ii =1, knon
          alb_new(ii) = alb_eau(knindex(ii))
       enddo

       z0_new = sqrt(rugos**2 + rugoro**2)
       alblw(1:knon) = alb_new(1:knon)

       !
    else if (nisurf == is_sic) then

       if (check) write(*,*)'sea ice, nisurf = ',nisurf 

       !
       ! Surface "glace de mer" appel a l'interface avec l'ocean
       !
       !
       if (ocean == 'couple') then

          cumul =.false.

          iloc = maxloc(fder(1:klon))
          if (check.and.fder(iloc(1))> 0.) then
             WRITE(*,*)'**** Debug fder ****'
             WRITE(*,*)'max fder(',iloc(1),') = ',fder(iloc(1))
          endif
!!$
!!$      where(fder.gt.0.) 
!!$        fder = 0.
!!$      endwhere

          call interfoce_cpl(itime, dtime, cumul, &
               & klon, iim, jjm, nisurf, pctsrf, knon, knindex, rlon, rlat, &
               & ocean, npas, nexca, debut, lafin, &
               & swdown, sollw, precip_rain, precip_snow, evap, tsurf, &
               & fluxlat, fluxsens, fder, albedo, taux, tauy, &
               & windsp, &
               & zmasq, &
               & tsurf_new, alb_new, &
               & pctsrf_new)

          tsurf_temp = tsurf_new
          cal = 0.
          dif_grnd = 0.
          beta = 1.0

          !IM: "slab" ocean
       else if (ocean == 'slab  ') then
          pctsrf_new=tmp_pctsrf_slab
          !
          DO ii = 1, knon
             tsurf_new(ii) = tsurf(ii)
             IF (pctsrf_new(knindex(ii),nisurf) < EPSFRA) then
                snow(ii) = 0.0
                tsurf_new(ii) = RTT - 1.8
                IF (soil_model) tsoil(ii,:) = RTT -1.8
             ENDIF
          ENDDO

          CALL calbeta(dtime, nisurf, knon, snow, qsol, beta, capsol, dif_grnd)

          IF (soil_model) THEN 
             CALL soil(dtime, nisurf, knon,snow, tsurf_new, tsoil,soilcap, soilflux)
             cal(1:knon) = RCPD / soilcap(1:knon)
             radsol(1:knon) = radsol(1:knon)  + soilflux(1:knon)
          ELSE 
             dif_grnd = 1.0 / tau_gl
             cal = RCPD * calice
             WHERE (snow > 0.0) cal = RCPD * calsno 
          ENDIF
          tsurf_temp = tsurf_new
          beta = 1.0
          !
       ELSE
          !                              ! lecture conditions limites
          CALL interfoce_lim(itime, dtime, jour, & 
               &  klon, nisurf, knon, knindex, &
               &  debut, &
               &  tsurf_new, pctsrf_new)

          !IM cf LF
          DO ii = 1, knon
             tsurf_new(ii) = tsurf(ii)
             !IMbad IF (pctsrf_new(ii,nisurf) < EPSFRA) then
             IF (pctsrf_new(knindex(ii),nisurf) < EPSFRA) then
                snow(ii) = 0.0
                !IM cf LF/JLD         tsurf(ii) = RTT - 1.8
                tsurf_new(ii) = RTT - 1.8
                IF (soil_model) tsoil(ii,:) = RTT -1.8
             endif
          enddo

          CALL calbeta(dtime, nisurf, knon, snow, qsol, beta, capsol, dif_grnd)

          IF (soil_model) THEN 
             !IM cf LF/JLD        CALL soil(dtime, nisurf, knon,snow, tsurf, tsoil,soilcap, soilflux)
             CALL soil(dtime, nisurf, knon,snow, tsurf_new, tsoil,soilcap, soilflux)
             cal(1:knon) = RCPD / soilcap(1:knon)
             radsol(1:knon) = radsol(1:knon)  + soilflux(1:knon)
             dif_grnd = 0.
          ELSE 
             dif_grnd = 1.0 / tau_gl
             cal = RCPD * calice
             WHERE (snow > 0.0) cal = RCPD * calsno 
          ENDIF
          !IMbadtsurf_temp = tsurf
          tsurf_temp = tsurf_new
          beta = 1.0
       ENDIF

       CALL calcul_fluxs( klon, knon, nisurf, dtime, &
            &   tsurf_temp, p1lay, cal, beta, tq_cdrag, ps, &
            &   precip_rain, precip_snow, snow, qsurf,  &
            &   radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, &
            &   petAcoef, peqAcoef, petBcoef, peqBcoef, &
            &   tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l)
       !
       !IM: flux entre l'ocean et la glace de mer pour le "slab" ocean
       DO i = 1, knon
          flux_g(i) = 0.0
          !
          !IM: faire dependre le coefficient de conduction de la glace de mer
          !    de l'epaisseur de la glace de mer, dans l'hypothese ou le coeff.
          !    actuel correspond a 3m de glace de mer, cf. L.Li
          !
          !      IF(1.EQ.0) THEN
          !       IF(siceh(i).GT.0.) THEN
          !        new_dif_grnd(i) = dif_grnd(i)*3./siceh(i)
          !       ELSE
          !        new_dif_grnd(i) = 0.
          !       ENDIF
          !      ENDIF !(1.EQ.0) THEN
          !
          IF (cal(i).GT.1.0e-15) flux_g(i)=(tsurf_new(i)-t_grnd) &
               &                          * dif_grnd(i) *RCPD/cal(i)
          !   &                          * new_dif_grnd(i) *RCPD/cal(i)
          tmp_flux_g(knindex(i))=flux_g(i)
          tmp_radsol(knindex(i))=radsol(i)
       ENDDO

       IF (ocean /= 'couple') THEN
          CALL fonte_neige( klon, knon, nisurf, dtime, &
               &   tsurf_temp, p1lay, cal, beta, tq_cdrag, ps, &
               &   precip_rain, precip_snow, snow, qsol,  &
               &   radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, &
               &   petAcoef, peqAcoef, petBcoef, peqBcoef, &
               &   tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l, &
               &   fqcalving,ffonte, run_off_lic_0)

          !     calcul albedo

          CALL albsno(klon,knon,dtime,agesno(:),alb_neig(:), precip_snow(:))  
          WHERE (snow(1 : knon) .LT. 0.0001) agesno(1 : knon) = 0.
          zfra(1:knon) = MAX(0.0,MIN(1.0,snow(1:knon)/(snow(1:knon)+10.0)))
          alb_new(1 : knon) = alb_neig(1 : knon) *zfra(1:knon) + & 
               &                    0.6 * (1.0-zfra(1:knon))
          !!      alb_new(1 : knon) = 0.6
       ENDIF

       fder_prev = fder    
       fder = fder_prev + dflux_s + dflux_l

       iloc = maxloc(fder(1:klon))
       if (check.and.fder(iloc(1))> 0.) then
          WRITE(*,*)'**** Debug fder ****'
          WRITE(*,*)'max fder(',iloc(1),') = ',fder(iloc(1))
          WRITE(*,*)'fder_prev, dflux_s, dflux_l',fder_prev(iloc(1)), &
               &                        dflux_s(iloc(1)), dflux_l(iloc(1))
       endif
!!$      where(fder.gt.0.) 
!!$        fder = 0.
!!$      endwhere

       !
       ! 2eme appel a interfoce pour le cumul et le passage des flux a l'ocean
       !
       if (ocean == 'couple') then

          cumul =.true.

          call interfoce_cpl(itime, dtime, cumul, &
               & klon, iim, jjm, nisurf, pctsrf, knon, knindex, rlon, rlat, &
               & ocean, npas, nexca, debut, lafin, &
               & swdown, sollw, precip_rain, precip_snow, evap, tsurf, &
               & fluxlat, fluxsens, fder, albedo, taux, tauy, &
               & windsp, &
               & zmasq, &
               & tsurf_new, alb_new, &
               & pctsrf_new)
       endif

       z0_new = 0.002
       z0_new = SQRT(z0_new**2+rugoro**2)
       alblw(1:knon) = alb_new(1:knon)

    else if (nisurf == is_lic) then

       if (check) write(*,*)'glacier, nisurf = ',nisurf 

       if (.not. allocated(run_off_lic)) then
          allocate(run_off_lic(knon), stat = error)
          if (error /= 0) then
             abort_message='Pb allocation run_off_lic'
             call abort_gcm(modname,abort_message,1)
          endif
          run_off_lic = 0.
       endif
       !
       ! Surface "glacier continentaux" appel a l'interface avec le sol
       !
       IF (soil_model) THEN 
          CALL soil(dtime, nisurf, knon, snow, tsurf, tsoil,soilcap, soilflux)
          cal(1:knon) = RCPD / soilcap(1:knon)
          radsol(1:knon)  = radsol(1:knon) + soilflux(1:knon)
       ELSE 
          cal = RCPD * calice
          WHERE (snow > 0.0) cal = RCPD * calsno
       ENDIF
       beta = 1.0
       dif_grnd = 0.0

       call calcul_fluxs( klon, knon, nisurf, dtime, &
            &   tsurf, p1lay, cal, beta, tq_cdrag, ps, &
            &   precip_rain, precip_snow, snow, qsurf,  &
            &   radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, &
            &   petAcoef, peqAcoef, petBcoef, peqBcoef, &
            &   tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l)

       call fonte_neige( klon, knon, nisurf, dtime, &
            &   tsurf, p1lay, cal, beta, tq_cdrag, ps, &
            &   precip_rain, precip_snow, snow, qsol,  &
            &   radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, &
            &   petAcoef, peqAcoef, petBcoef, peqBcoef, &
            &   tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l, &
            &   fqcalving,ffonte, run_off_lic_0)

       ! passage du run-off des glaciers calcule dans fonte_neige au coupleur
       bidule=0.
       bidule(1:knon)= run_off_lic(1:knon)    
       call gath2cpl(bidule, tmp_rlic, klon, knon,iim,jjm,knindex)
       !
       ! calcul albedo
       !
       CALL albsno(klon,knon,dtime,agesno(:),alb_neig(:), precip_snow(:))  
       WHERE (snow(1 : knon) .LT. 0.0001) agesno(1 : knon) = 0.
       zfra(1:knon) = MAX(0.0,MIN(1.0,snow(1:knon)/(snow(1:knon)+10.0)))
       alb_new(1 : knon)  = alb_neig(1 : knon)*zfra(1:knon) + &
            &                     0.6 * (1.0-zfra(1:knon))
       !
       !IM: plusieurs choix/tests sur l'albedo des "glaciers continentaux"
       !       alb_new(1 : knon)  = 0.6 !IM cf FH/GK 
       !       alb_new(1 : knon)  = 0.82
       !       alb_new(1 : knon)  = 0.77 !211003 Ksta0.77
       !       alb_new(1 : knon)  = 0.8 !KstaTER0.8 & LMD_ARMIP5
       !IM: KstaTER0.77 & LMD_ARMIP6    
       alb_new(1 : knon)  = 0.77

       !
       ! Rugosite
       !
       z0_new = rugoro
       !
       ! Remplissage des pourcentages de surface
       !
       pctsrf_new(:,nisurf) = pctsrf(:,nisurf)

       alblw(1:knon) = alb_new(1:knon)
    else
       write(*,*)'Index surface = ',nisurf
       abort_message = 'Index surface non valable'
       call abort_gcm(modname,abort_message,1)
    endif

  END SUBROUTINE interfsurf_hq

  !************************

  SUBROUTINE interfoce_cpl(itime, dtime, cumul, &
       & klon, iim, jjm, nisurf, pctsrf, knon, knindex, rlon, rlat, &
       & ocean, npas, nexca, debut, lafin, &
       & swdown, lwdown, precip_rain, precip_snow, evap, tsurf, &
       & fluxlat, fluxsens, fder, albsol, taux, tauy, &
       & windsp, &
       & zmasq, &
       & tsurf_new, alb_new, &
       & pctsrf_new)

    ! Cette routine sert d'interface entre le modele atmospherique et un 
    ! coupleur avec un modele d'ocean 'complet' derriere
    !
    ! Le modele de glace qu'il est prevu d'utiliser etant couple directement a 
    ! l'ocean presentement, on va passer deux fois dans cette routine par pas de 
    ! temps physique, une fois avec les points oceans et l'autre avec les points
    ! glace. A chaque pas de temps de couplage, la lecture des champs provenant
    ! du coupleur se fera "dans" l'ocean et l'ecriture des champs a envoyer
    ! au coupleur "dans" la glace. Il faut donc des tableaux de travail "tampons"
    ! dimensionnes sur toute la grille qui remplissent les champs sur les
    ! domaines ocean/glace quand il le faut. Il est aussi necessaire que l'index
    ! ocean soit traiter avant l'index glace (sinon tout intervertir)
    !
    !
    ! L. Fairhead 02/2000
    !
    ! input:
    !   itime        numero du pas de temps
    !   iim, jjm     nbres de pts de grille
    !   dtime        pas de temps de la physique
    !   klon         nombre total de points de grille
    !   nisurf       index de la surface a traiter (1 = sol continental)
    !   pctsrf       tableau des fractions de surface de chaque maille
    !   knon         nombre de points de la surface a traiter
    !   knindex      index des points de la surface a traiter
    !   rlon         longitudes
    !   rlat         latitudes
    !   debut        logical: 1er appel a la physique
    !   lafin        logical: dernier appel a la physique
    !   ocean        type d'ocean
    !   nexca        frequence de couplage
    !   swdown       flux solaire entrant a la surface
    !   lwdown       flux IR net a la surface
    !   precip_rain  precipitation liquide
    !   precip_snow  precipitation solide
    !   evap         evaporation
    !   tsurf        temperature de surface
    !   fder         derivee dF/dT
    !   albsol       albedo du sol (coherent avec swdown)
    !   taux         tension de vent en x
    !   tauy         tension de vent en y
    !    windsp       module du vent a 10m
    !   nexca        frequence de couplage
    !   zmasq        masque terre/ocean
    !
    !
    ! output:
    !   tsurf_new    temperature au sol
    !   alb_new      albedo
    !   pctsrf_new   nouvelle repartition des surfaces
    !   alb_ice      albedo de la glace
    !
    use temps
    use iniprint
    use abort_gcm_m, only: abort_gcm
    use gath_cpl, only: gath2cpl, cpl2gath
    use ioipsl
    use indicesol
    use YOMCST

    ! Parametres d'entree
    integer, intent(IN) :: itime
    integer, intent(IN) :: iim, jjm
    real, intent(IN) :: dtime
    integer, intent(IN) :: klon
    integer, intent(IN) :: nisurf
    integer, intent(IN) :: knon
    real, dimension(klon,nbsrf), intent(IN) :: pctsrf
    integer, dimension(klon), intent(in) :: knindex
    logical, intent(IN) :: debut, lafin
    real, dimension(klon), intent(IN) :: rlon, rlat
    character (len = 6)  :: ocean
    real, dimension(klon), intent(IN) :: lwdown, swdown
    real, dimension(klon), intent(IN) :: precip_rain, precip_snow
    real, dimension(klon), intent(IN) :: tsurf, fder, albsol, taux, tauy
    real, dimension(klon), intent(IN) :: windsp
    INTEGER              :: nexca, npas
    real, dimension(klon), intent(IN) :: zmasq
    real, dimension(klon), intent(IN) :: fluxlat, fluxsens
    logical, intent(IN)               :: cumul
    real, dimension(klon), intent(INOUT) :: evap

    ! Parametres de sortie
    real, dimension(klon), intent(OUT):: tsurf_new, alb_new
    real, dimension(klon,nbsrf), intent(OUT) :: pctsrf_new

    ! Variables locales
    integer                    :: j, error, sum_error, ig, cpl_index,i
    INTEGER :: nsrf
    character (len = 20) :: modname = 'interfoce_cpl'
    character (len = 80) :: abort_message
    logical,save              :: check = .FALSE.
    ! variables pour moyenner les variables de couplage
    real, allocatable, dimension(:,:),save :: cpl_sols, cpl_nsol, cpl_rain
    real, allocatable, dimension(:,:),save :: cpl_snow, cpl_evap, cpl_tsol
    real, allocatable, dimension(:,:),save :: cpl_fder, cpl_albe, cpl_taux
    real, allocatable, dimension(:,:),save :: cpl_windsp
    real, allocatable, dimension(:,:),save :: cpl_tauy
    REAL, ALLOCATABLE, DIMENSION(:,:),SAVE :: cpl_rriv, cpl_rcoa, cpl_rlic
!!$
    ! variables tampons avant le passage au coupleur
    real, allocatable, dimension(:,:,:),save :: tmp_sols, tmp_nsol, tmp_rain
    real, allocatable, dimension(:,:,:),save :: tmp_snow, tmp_evap, tmp_tsol
    real, allocatable, dimension(:,:,:),save :: tmp_fder, tmp_albe, tmp_taux
    real, allocatable, dimension(:,:,:),save :: tmp_windsp
    REAL, ALLOCATABLE, DIMENSION(:,:,:),SAVE :: tmp_tauy
    ! variables a passer au coupleur
    real, dimension(iim, jjm+1) :: wri_sol_ice, wri_sol_sea, wri_nsol_ice 
    real, dimension(iim, jjm+1) :: wri_nsol_sea, wri_fder_ice, wri_evap_ice
    REAL, DIMENSION(iim, jjm+1) :: wri_evap_sea, wri_rcoa, wri_rriv
    REAL, DIMENSION(iim, jjm+1) :: wri_rain, wri_snow, wri_taux, wri_tauy
    REAL, DIMENSION(iim, jjm+1) :: wri_windsp
    REAL, DIMENSION(iim, jjm+1) :: wri_calv
    REAL, DIMENSION(iim, jjm+1) :: wri_tauxx, wri_tauyy, wri_tauzz
    REAL, DIMENSION(iim, jjm+1) :: tmp_lon, tmp_lat
    ! variables relues par le coupleur
    ! read_sic = fraction de glace
    ! read_sit = temperature de glace
    real, allocatable, dimension(:,:),save :: read_sst, read_sic, read_sit
    real, allocatable, dimension(:,:),save :: read_alb_sic
    ! variable tampon
    real, dimension(klon)       :: tamp_sic
    ! sauvegarde des fractions de surface d'un pas de temps a l'autre apres 
    ! l'avoir lu
    real, allocatable,dimension(:,:),save :: pctsrf_sav
    real, dimension(iim, jjm+1, 2) :: tamp_srf
    integer, allocatable, dimension(:), save :: tamp_ind
    real, allocatable, dimension(:,:),save :: tamp_zmasq
    real, dimension(iim, jjm+1) :: deno
    integer                     :: idtime
    integer, allocatable,dimension(:),save :: unity
    ! 
    logical, save    :: first_appel = .true.
    logical,save          :: print
    !maf
    ! variables pour avoir une sortie IOIPSL des champs echanges
    CHARACTER(len=80),SAVE :: clintocplnam, clfromcplnam
    INTEGER, SAVE :: jf,nhoridct,nidct
    INTEGER, SAVE :: nhoridcs,nidcs
    INTEGER :: ndexct(iim*(jjm+1)),ndexcs(iim*(jjm+1))
    REAL :: zx_lon(iim,jjm+1), zx_lat(iim,jjm+1), zjulian
    INTEGER,save :: idayref
    !med  integer :: itau_w
    integer,save :: itau_w
    integer :: nb_interf_cpl
    include "param_cou.h"
    include "inc_cpl.h"
    !
    ! Initialisation
    !
    if (check) write(*,*)'Entree ',modname,'nisurf = ',nisurf

    if (first_appel) then
       error = 0
       allocate(unity(klon), stat = error)
       if ( error  /=0) then
          abort_message='Pb allocation variable unity'
          call abort_gcm(modname,abort_message,1)
       endif
       allocate(pctsrf_sav(klon,nbsrf), stat = error)
       if ( error  /=0) then
          abort_message='Pb allocation variable pctsrf_sav'
          call abort_gcm(modname,abort_message,1)
       endif
       pctsrf_sav = 0.

       do ig = 1, klon
          unity(ig) = ig
       enddo
       sum_error = 0
       allocate(cpl_sols(klon,2), stat = error); sum_error = sum_error + error
       allocate(cpl_nsol(klon,2), stat = error); sum_error = sum_error + error
       allocate(cpl_rain(klon,2), stat = error); sum_error = sum_error + error
       allocate(cpl_snow(klon,2), stat = error); sum_error = sum_error + error
       allocate(cpl_evap(klon,2), stat = error); sum_error = sum_error + error
       allocate(cpl_tsol(klon,2), stat = error); sum_error = sum_error + error
       allocate(cpl_fder(klon,2), stat = error); sum_error = sum_error + error
       allocate(cpl_albe(klon,2), stat = error); sum_error = sum_error + error
       allocate(cpl_taux(klon,2), stat = error); sum_error = sum_error + error
       allocate(cpl_windsp(klon,2), stat = error); sum_error = sum_error + error
       allocate(cpl_tauy(klon,2), stat = error); sum_error = sum_error + error
       ALLOCATE(cpl_rriv(iim,jjm+1), stat=error); sum_error = sum_error + error
       ALLOCATE(cpl_rcoa(iim,jjm+1), stat=error); sum_error = sum_error + error
       ALLOCATE(cpl_rlic(iim,jjm+1), stat=error); sum_error = sum_error + error
       !!
       allocate(read_sst(iim, jjm+1), stat = error); sum_error = sum_error + error
       allocate(read_sic(iim, jjm+1), stat = error); sum_error = sum_error + error
       allocate(read_sit(iim, jjm+1), stat = error); sum_error = sum_error + error
       allocate(read_alb_sic(iim, jjm+1), stat = error); sum_error = sum_error + error

       if (sum_error /= 0) then
          abort_message='Pb allocation variables couplees'
          call abort_gcm(modname,abort_message,1)
       endif
       cpl_sols = 0.; cpl_nsol = 0.; cpl_rain = 0.; cpl_snow = 0.
       cpl_evap = 0.; cpl_tsol = 0.; cpl_fder = 0.; cpl_albe = 0.
       cpl_taux = 0.; cpl_tauy = 0.; cpl_rriv = 0.; cpl_rcoa = 0.; cpl_rlic = 0.
       cpl_windsp = 0.

       sum_error = 0
       allocate(tamp_ind(klon), stat = error); sum_error = sum_error + error
       allocate(tamp_zmasq(iim, jjm+1), stat = error); sum_error = sum_error + error    
       do ig = 1, klon
          tamp_ind(ig) = ig
       enddo
       call gath2cpl(zmasq, tamp_zmasq, klon, klon, iim, jjm, tamp_ind)
       !
       ! initialisation couplage
       !
       idtime = int(dtime)
       !
       ! initialisation sorties netcdf
       !
       idayref = day_ini
       CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)
       CALL gr_fi_ecrit(1,klon,iim,jjm+1,rlon,zx_lon)
       DO i = 1, iim
          zx_lon(i,1) = rlon(i+1)
          zx_lon(i,jjm+1) = rlon(i+1)
       ENDDO
       CALL gr_fi_ecrit(1,klon,iim,jjm+1,rlat,zx_lat)
       clintocplnam="cpl_atm_tauflx"
       CALL histbeg_totreg(clintocplnam, iim,zx_lon(:,1),jjm+1,zx_lat(1,:),1,iim,1,jjm+1, &
            & itau_phy,zjulian,dtime,nhoridct,nidct) 
       ! no vertical axis
       CALL histdef(nidct, 'tauxe','tauxe', &
            & "-",iim, jjm+1, nhoridct, 1, 1, 1, -99, 32, "inst", dtime,dtime)
       CALL histdef(nidct, 'tauyn','tauyn', &
            & "-",iim, jjm+1, nhoridct, 1, 1, 1, -99, 32, "inst", dtime,dtime)
       CALL histdef(nidct, 'tmp_lon','tmp_lon', &
            & "-",iim, jjm+1, nhoridct, 1, 1, 1, -99, 32, "inst", dtime,dtime)
       CALL histdef(nidct, 'tmp_lat','tmp_lat', &
            & "-",iim, jjm+1, nhoridct, 1, 1, 1, -99, 32, "inst", dtime,dtime)
       DO jf=1,jpflda2o1 + jpflda2o2
          CALL histdef(nidct, cl_writ(jf),cl_writ(jf), &
               & "-",iim, jjm+1, nhoridct, 1, 1, 1, -99, 32, "inst", dtime,dtime)
       END DO
       CALL histend(nidct)
       CALL histsync(nidct)

       clfromcplnam="cpl_atm_sst"
       CALL histbeg_totreg(clfromcplnam, iim,zx_lon(:,1),jjm+1,zx_lat(1,:),1,iim,1,jjm+1, &
            & 0,zjulian,dtime,nhoridcs,nidcs) 
       ! no vertical axis
       DO jf=1,jpfldo2a
          CALL histdef(nidcs, cl_read(jf),cl_read(jf), &
               & "-",iim, jjm+1, nhoridcs, 1, 1, 1, -99, 32, "inst", dtime,dtime)
       END DO
       CALL histend(nidcs)
       CALL histsync(nidcs)

       ! pour simuler la fonte des glaciers antarctiques
       ! 
       surf_maille = (4. * rpi * ra**2) / (iim * (jjm +1))
       ALLOCATE(coeff_iceberg(iim,jjm+1), stat=error)
       if (error /= 0) then
          abort_message='Pb allocation variable coeff_iceberg'
          call abort_gcm(modname,abort_message,1)
       endif
       open (12,file='flux_iceberg',form='formatted',status='old')
       read (12,*) coeff_iceberg
       close (12)
       num_antarctic = max(1, count(coeff_iceberg > 0))

       first_appel = .false.
    endif ! fin if (first_appel)

    ! Initialisations

    ! calcul des fluxs a passer
    nb_interf_cpl = nb_interf_cpl + 1  
    if (check) write(lunout,*)'passage dans interface_surf.F90 :  ',nb_interf_cpl 
    cpl_index = 1
    if (nisurf == is_sic) cpl_index = 2
    if (cumul) then
       if (check) write(lunout,*)'passage dans cumul '
       if (check) write(lunout,*)'valeur de cpl_index ', cpl_index 
       ! -- LOOP  
       if (check) write(*,*) modname, 'cumul des champs'
       do ig = 1, knon
          cpl_sols(ig,cpl_index) = cpl_sols(ig,cpl_index) &
               &                          + swdown(ig)      / FLOAT(nexca)
          cpl_nsol(ig,cpl_index) = cpl_nsol(ig,cpl_index) &
               &                          + (lwdown(ig) + fluxlat(ig) +fluxsens(ig))&
               &                                / FLOAT(nexca)
          cpl_rain(ig,cpl_index) = cpl_rain(ig,cpl_index) &
               &                          + precip_rain(ig) / FLOAT(nexca)
          cpl_snow(ig,cpl_index) = cpl_snow(ig,cpl_index) &
               &                          + precip_snow(ig) / FLOAT(nexca)
          cpl_evap(ig,cpl_index) = cpl_evap(ig,cpl_index) &
               &                          + evap(ig)        / FLOAT(nexca)
          cpl_tsol(ig,cpl_index) = cpl_tsol(ig,cpl_index) &
               &                          + tsurf(ig)       / FLOAT(nexca)
          cpl_fder(ig,cpl_index) = cpl_fder(ig,cpl_index) &
               &                          + fder(ig)        / FLOAT(nexca)
          cpl_albe(ig,cpl_index) = cpl_albe(ig,cpl_index) &
               &                          + albsol(ig)      / FLOAT(nexca)
          cpl_taux(ig,cpl_index) = cpl_taux(ig,cpl_index) &
               &                          + taux(ig)        / FLOAT(nexca)
          cpl_tauy(ig,cpl_index) = cpl_tauy(ig,cpl_index) &
               &                          + tauy(ig)        / FLOAT(nexca)
          IF (cpl_index .EQ. 1) THEN
             cpl_windsp(ig,cpl_index) = cpl_windsp(ig,cpl_index) &
                  &                          + windsp(ig)      / FLOAT(nexca) 
          ENDIF
       enddo
       IF (cpl_index .EQ. 1) THEN 
          cpl_rriv(:,:) = cpl_rriv(:,:) + tmp_rriv(:,:) / FLOAT(nexca)
          cpl_rcoa(:,:) = cpl_rcoa(:,:) + tmp_rcoa(:,:) / FLOAT(nexca)
          cpl_rlic(:,:) = cpl_rlic(:,:) + tmp_rlic(:,:) / FLOAT(nexca)
       ENDIF
    endif

    if (mod(itime, nexca) == 1) then
       !
       ! Demande des champs au coupleur
       !
       ! Si le domaine considere est l'ocean, on lit les champs venant du coupleur
       !
       if (nisurf == is_oce .and. .not. cumul) then
          if (check) write(*,*)'rentree fromcpl, itime-1 = ',itime-1
          !
          ! sorties NETCDF des champs recus
          !
          ndexcs(:)=0
          itau_w = itau_phy + itime
          CALL histwrite(nidcs,cl_read(1),itau_w,read_sst,iim*(jjm+1),ndexcs)
          CALL histwrite(nidcs,cl_read(2),itau_w,read_sic,iim*(jjm+1),ndexcs)
          CALL histwrite(nidcs,cl_read(3),itau_w,read_alb_sic,iim*(jjm+1),ndexcs)
          CALL histwrite(nidcs,cl_read(4),itau_w,read_sit,iim*(jjm+1),ndexcs)
          CALL histsync(nidcs)
          ! pas utile      IF (npas-itime.LT.nexca )CALL histclo(nidcs)

          do j = 1, jjm + 1
             do ig = 1, iim
                if (abs(1. - read_sic(ig,j)) < 0.00001) then
                   read_sst(ig,j) = RTT - 1.8
                   read_sit(ig,j) = read_sit(ig,j) / read_sic(ig,j)
                   read_alb_sic(ig,j) = read_alb_sic(ig,j) / read_sic(ig,j)
                else if (abs(read_sic(ig,j)) < 0.00001) then
                   read_sst(ig,j) = read_sst(ig,j) / (1. - read_sic(ig,j))
                   read_sit(ig,j) = read_sst(ig,j)
                   read_alb_sic(ig,j) =  0.6
                else
                   read_sst(ig,j) = read_sst(ig,j) / (1. - read_sic(ig,j))
                   read_sit(ig,j) = read_sit(ig,j) / read_sic(ig,j)
                   read_alb_sic(ig,j) = read_alb_sic(ig,j) / read_sic(ig,j)
                endif
             enddo
          enddo
          !
          ! transformer read_sic en pctsrf_sav
          !
          call cpl2gath(read_sic, tamp_sic , klon, klon,iim,jjm, unity)
          do ig = 1, klon
             IF (pctsrf(ig,is_oce) > epsfra .OR.            &
                  &             pctsrf(ig,is_sic) > epsfra) THEN
                pctsrf_sav(ig,is_sic) = (pctsrf(ig,is_oce) + pctsrf(ig,is_sic)) &
                     &                               * tamp_sic(ig)
                pctsrf_sav(ig,is_oce) = (pctsrf(ig,is_oce) + pctsrf(ig,is_sic)) &
                     &                        - pctsrf_sav(ig,is_sic)
             endif
          enddo
          !
          ! Pour rattraper des erreurs d'arrondis
          !
          where (abs(pctsrf_sav(:,is_sic)) .le. 2.*epsilon(pctsrf_sav(1,is_sic)))
             pctsrf_sav(:,is_sic) = 0.
             pctsrf_sav(:,is_oce) = pctsrf(:,is_oce) + pctsrf(:,is_sic)
          endwhere
          where (abs(pctsrf_sav(:,is_oce)) .le. 2.*epsilon(pctsrf_sav(1,is_oce)))
             pctsrf_sav(:,is_sic) = pctsrf(:,is_oce) + pctsrf(:,is_sic)
             pctsrf_sav(:,is_oce) = 0.
          endwhere
          if (minval(pctsrf_sav(:,is_oce)) < 0.) then
             write(*,*)'Pb fraction ocean inferieure a 0'
             write(*,*)'au point ',minloc(pctsrf_sav(:,is_oce))
             write(*,*)'valeur = ',minval(pctsrf_sav(:,is_oce)) 
             abort_message = 'voir ci-dessus'
             call abort_gcm(modname,abort_message,1)
          endif
          if (minval(pctsrf_sav(:,is_sic)) < 0.) then
             write(*,*)'Pb fraction glace inferieure a 0'
             write(*,*)'au point ',minloc(pctsrf_sav(:,is_sic))
             write(*,*)'valeur = ',minval(pctsrf_sav(:,is_sic)) 
             abort_message = 'voir ci-dessus'
             call abort_gcm(modname,abort_message,1)
          endif
       endif
    endif                         ! fin mod(itime, nexca) == 1

    if (mod(itime, nexca) == 0) then
       !
       ! allocation memoire
       if (nisurf == is_oce .and. (.not. cumul) ) then
          sum_error = 0
          allocate(tmp_sols(iim,jjm+1,2), stat=error); sum_error = sum_error + error
          allocate(tmp_nsol(iim,jjm+1,2), stat=error); sum_error = sum_error + error
          allocate(tmp_rain(iim,jjm+1,2), stat=error); sum_error = sum_error + error
          allocate(tmp_snow(iim,jjm+1,2), stat=error); sum_error = sum_error + error
          allocate(tmp_evap(iim,jjm+1,2), stat=error); sum_error = sum_error + error
          allocate(tmp_tsol(iim,jjm+1,2), stat=error); sum_error = sum_error + error
          allocate(tmp_fder(iim,jjm+1,2), stat=error); sum_error = sum_error + error
          allocate(tmp_albe(iim,jjm+1,2), stat=error); sum_error = sum_error + error
          allocate(tmp_taux(iim,jjm+1,2), stat=error); sum_error = sum_error + error
          allocate(tmp_tauy(iim,jjm+1,2), stat=error); sum_error = sum_error + error
          allocate(tmp_windsp(iim,jjm+1,2), stat=error); sum_error = sum_error + error
!!$      allocate(tmp_rriv(iim,jjm+1,2), stat=error); sum_error = sum_error + error
!!$      allocate(tmp_rcoa(iim,jjm+1,2), stat=error); sum_error = sum_error + error
          if (sum_error /= 0) then
             abort_message='Pb allocation variables couplees pour l''ecriture'
             call abort_gcm(modname,abort_message,1)
          endif
       endif

       !
       ! Mise sur la bonne grille des champs a passer au coupleur
       !
       cpl_index = 1
       if (nisurf == is_sic) cpl_index = 2
       call gath2cpl(cpl_sols(1,cpl_index), tmp_sols(1,1,cpl_index), klon, knon,iim,jjm,                  knindex)
       call gath2cpl(cpl_nsol(1,cpl_index), tmp_nsol(1,1,cpl_index), klon, knon,iim,jjm,                  knindex)
       call gath2cpl(cpl_rain(1,cpl_index), tmp_rain(1,1,cpl_index), klon, knon,iim,jjm,                  knindex)
       call gath2cpl(cpl_snow(1,cpl_index), tmp_snow(1,1,cpl_index), klon, knon,iim,jjm,                  knindex)
       call gath2cpl(cpl_evap(1,cpl_index), tmp_evap(1,1,cpl_index), klon, knon,iim,jjm,                  knindex)
       call gath2cpl(cpl_tsol(1,cpl_index), tmp_tsol(1,1,cpl_index), klon, knon,iim,jjm,                  knindex)
       call gath2cpl(cpl_fder(1,cpl_index), tmp_fder(1,1,cpl_index), klon, knon,iim,jjm,                  knindex)
       call gath2cpl(cpl_albe(1,cpl_index), tmp_albe(1,1,cpl_index), klon, knon,iim,jjm,                  knindex)
       call gath2cpl(cpl_taux(1,cpl_index), tmp_taux(1,1,cpl_index), klon, knon,iim,jjm,                  knindex)
       call gath2cpl(cpl_windsp(1,cpl_index), tmp_windsp(1,1,cpl_index), klon, knon,iim,jjm,             knindex)
       call gath2cpl(cpl_tauy(1,cpl_index), tmp_tauy(1,1,cpl_index), klon, knon,iim,jjm,                  knindex)

       !
       ! Si le domaine considere est la banquise, on envoie les champs au coupleur
       !
       if (nisurf == is_sic .and. cumul) then
          wri_rain = 0.; wri_snow = 0.; wri_rcoa = 0.; wri_rriv = 0. 
          wri_taux = 0.; wri_tauy = 0.
          wri_windsp = 0.
          ! -- LOOP      
          call gath2cpl(pctsrf(1,is_oce), tamp_srf(1,1,1), klon, klon, iim, jjm, tamp_ind)
          call gath2cpl(pctsrf(1,is_sic), tamp_srf(1,1,2), klon, klon, iim, jjm, tamp_ind)

          wri_sol_ice = tmp_sols(:,:,2)
          wri_sol_sea = tmp_sols(:,:,1)
          wri_nsol_ice = tmp_nsol(:,:,2)
          wri_nsol_sea = tmp_nsol(:,:,1)
          wri_fder_ice = tmp_fder(:,:,2)
          wri_evap_ice = tmp_evap(:,:,2)
          wri_evap_sea = tmp_evap(:,:,1)
          wri_windsp = tmp_windsp(:,:,1)

!!$PB
          wri_rriv = cpl_rriv(:,:)
          wri_rcoa = cpl_rcoa(:,:)
          DO j = 1, jjm + 1
             wri_calv(:,j) = sum(cpl_rlic(:,j)) / iim
          enddo

          where (tamp_zmasq /= 1.)
             deno =  tamp_srf(:,:,1) + tamp_srf(:,:,2)
             wri_rain = tmp_rain(:,:,1) * tamp_srf(:,:,1) / deno +    &
                  &            tmp_rain(:,:,2) * tamp_srf(:,:,2) / deno
             wri_snow = tmp_snow(:,:,1) * tamp_srf(:,:,1) / deno +    &
                  &            tmp_snow(:,:,2) * tamp_srf(:,:,2) / deno
             wri_taux = tmp_taux(:,:,1) * tamp_srf(:,:,1) / deno +    &
                  &            tmp_taux(:,:,2) * tamp_srf(:,:,2) / deno
             wri_tauy = tmp_tauy(:,:,1) * tamp_srf(:,:,1) / deno +    &
                  &            tmp_tauy(:,:,2) * tamp_srf(:,:,2) / deno
          endwhere
          !
          ! pour simuler la fonte des glaciers antarctiques
          !
          !$$$        wri_rain = wri_rain      &
          !$$$      &     + coeff_iceberg * cte_flux_iceberg / (num_antarctic * surf_maille)
          !      wri_calv = coeff_iceberg * cte_flux_iceberg / (num_antarctic * surf_maille)
          !
          ! on passe les coordonnées de la grille
          !

          CALL gr_fi_ecrit(1,klon,iim,jjm+1,rlon,tmp_lon)
          CALL gr_fi_ecrit(1,klon,iim,jjm+1,rlat,tmp_lat)

          DO i = 1, iim
             tmp_lon(i,1) = rlon(i+1)
             tmp_lon(i,jjm + 1) = rlon(i+1)
          ENDDO
          !
          ! sortie netcdf des champs pour le changement de repere
          !
          ndexct(:)=0
          CALL histwrite(nidct,'tauxe',itau_w,wri_taux,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,'tauyn',itau_w,wri_tauy,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,'tmp_lon',itau_w,tmp_lon,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,'tmp_lat',itau_w,tmp_lat,iim*(jjm+1),ndexct)

          !
          ! calcul 3 coordonnées du vent
          !
          CALL atm2geo (iim , jjm + 1, wri_taux, wri_tauy, tmp_lon, tmp_lat, &
               & wri_tauxx, wri_tauyy, wri_tauzz )
          !
          ! sortie netcdf des champs apres changement de repere et juste avant
          ! envoi au coupleur
          !
          CALL histwrite(nidct,cl_writ(8),itau_w,wri_sol_ice,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,cl_writ(9),itau_w,wri_sol_sea,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,cl_writ(10),itau_w,wri_nsol_ice,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,cl_writ(11),itau_w,wri_nsol_sea,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,cl_writ(12),itau_w,wri_fder_ice,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,cl_writ(13),itau_w,wri_evap_ice,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,cl_writ(14),itau_w,wri_evap_sea,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,cl_writ(15),itau_w,wri_rain,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,cl_writ(16),itau_w,wri_snow,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,cl_writ(17),itau_w,wri_rcoa,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,cl_writ(18),itau_w,wri_rriv,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,cl_writ(19),itau_w,wri_calv,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,cl_writ(1),itau_w,wri_tauxx,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,cl_writ(2),itau_w,wri_tauyy,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,cl_writ(3),itau_w,wri_tauzz,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,cl_writ(4),itau_w,wri_tauxx,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,cl_writ(5),itau_w,wri_tauyy,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,cl_writ(6),itau_w,wri_tauzz,iim*(jjm+1),ndexct)
          CALL histwrite(nidct,cl_writ(7),itau_w,wri_windsp,iim*(jjm+1),ndexct)
          CALL histsync(nidct)
          ! pas utile      IF (lafin) CALL histclo(nidct)
          ! 
          cpl_sols = 0.; cpl_nsol = 0.; cpl_rain = 0.; cpl_snow = 0.
          cpl_evap = 0.; cpl_tsol = 0.; cpl_fder = 0.; cpl_albe = 0.
          cpl_taux = 0.; cpl_tauy = 0.; cpl_rriv = 0.; cpl_rcoa = 0.; cpl_rlic = 0.
          cpl_windsp = 0.
          !
          ! deallocation memoire variables temporaires
          !
          sum_error = 0
          deallocate(tmp_sols, stat=error); sum_error = sum_error + error
          deallocate(tmp_nsol, stat=error); sum_error = sum_error + error
          deallocate(tmp_rain, stat=error); sum_error = sum_error + error
          deallocate(tmp_snow, stat=error); sum_error = sum_error + error
          deallocate(tmp_evap, stat=error); sum_error = sum_error + error
          deallocate(tmp_fder, stat=error); sum_error = sum_error + error
          deallocate(tmp_tsol, stat=error); sum_error = sum_error + error
          deallocate(tmp_albe, stat=error); sum_error = sum_error + error
          deallocate(tmp_taux, stat=error); sum_error = sum_error + error
          deallocate(tmp_tauy, stat=error); sum_error = sum_error + error
          deallocate(tmp_windsp, stat=error); sum_error = sum_error + error
          if (sum_error /= 0) then
             abort_message='Pb deallocation variables couplees'
             call abort_gcm(modname,abort_message,1)
          endif

       endif

    endif            ! fin (mod(itime, nexca) == 0)
    !
    ! on range les variables lues/sauvegardees dans les bonnes variables de sortie
    !
    if (nisurf == is_oce) then
       call cpl2gath(read_sst, tsurf_new, klon, knon,iim,jjm, knindex)
    else if (nisurf == is_sic) then
       call cpl2gath(read_sit, tsurf_new, klon, knon,iim,jjm, knindex)
       call cpl2gath(read_alb_sic, alb_new, klon, knon,iim,jjm, knindex)
    endif
    pctsrf_new(:,nisurf) = pctsrf_sav(:,nisurf)

    !  if (lafin) call quitcpl

  END SUBROUTINE interfoce_cpl

  !************************

  SUBROUTINE interfoce_slab(klon, debut, itap, dtime, ijour, &
       & radsol, fluxo, fluxg, pctsrf, &
       & tslab, seaice, pctsrf_slab)
    !
    ! Cette routine calcule la temperature d'un slab ocean, la glace de mer 
    ! et les pourcentages de la maille couverte par l'ocean libre et/ou 
    ! la glace de mer pour un "slab" ocean de 50m
    !
    ! I. Musat 04.02.2005
    !
    ! input:
    !   klon         nombre total de points de grille
    !   debut        logical: 1er appel a la physique
    !   itap         numero du pas de temps
    !   dtime        pas de temps de la physique (en s)
    !   ijour        jour dans l'annee en cours
    !   radsol       rayonnement net au sol (LW + SW)
    !   fluxo        flux turbulent (sensible + latent) sur les mailles oceaniques 
    !   fluxg        flux de conduction entre la surface de la glace de mer et l'ocean
    !   pctsrf       tableau des pourcentages de surface de chaque maille
    ! output: 
    !   tslab        temperature de l'ocean libre
    !   seaice       glace de mer (kg/m2)
    !   pctsrf_slab  "pourcentages" (valeurs entre 0. et 1.) surfaces issus du slab
    !
    use indicesol
    use clesphys
    use abort_gcm_m, only: abort_gcm
    use YOMCST

    ! Parametres d'entree
    integer, intent(IN) :: klon
    logical, intent(IN) :: debut
    INTEGER, intent(IN) :: itap
    REAL, intent(IN) :: dtime
    INTEGER, intent(IN) :: ijour
    REAL, dimension(klon), intent(IN) :: radsol
    REAL, dimension(klon), intent(IN) :: fluxo
    REAL, dimension(klon), intent(IN) :: fluxg
    real, dimension(klon, nbsrf), intent(IN) :: pctsrf
    ! Parametres de sortie
    real, dimension(klon), intent(INOUT) :: tslab
    real, dimension(klon), intent(INOUT)        :: seaice ! glace de mer (kg/m2)
    real, dimension(klon, nbsrf), intent(OUT) :: pctsrf_slab
    !
    ! Variables locales :
    INTEGER, save :: lmt_pas, julien, idayvrai
    REAL, parameter :: unjour=86400.
    real, allocatable, dimension(:), save :: tmp_tslab, tmp_seaice
    REAL, allocatable, dimension(:), save :: slab_bils
    REAL, allocatable, dimension(:), save :: lmt_bils
    logical,save              :: check = .false.
    !
    REAL, parameter :: cyang=50.0 * 4.228e+06 ! capacite calorifique volumetrique de l'eau J/(m2 K)
    REAL, parameter :: cbing=0.334e+05        ! J/kg
    real, dimension(klon)                 :: siceh !hauteur de la glace de mer (m)
    INTEGER :: i
    integer :: sum_error, error
    REAL :: zz, za, zb
    !
    character (len = 80) :: abort_message
    character (len = 20) :: modname = 'interfoce_slab'
    !
    julien = MOD(ijour,360)
    sum_error = 0
    IF (debut) THEN
       allocate(slab_bils(klon), stat = error); sum_error = sum_error + error
       allocate(lmt_bils(klon), stat = error); sum_error = sum_error + error
       allocate(tmp_tslab(klon), stat = error); sum_error = sum_error + error
       allocate(tmp_seaice(klon), stat = error); sum_error = sum_error + error
       if (sum_error /= 0) then
          abort_message='Pb allocation var. slab_bils,lmt_bils,tmp_tslab,tmp_seaice'
          call abort_gcm(modname,abort_message,1)
       endif
       tmp_tslab=tslab
       tmp_seaice=seaice
       lmt_pas = nint(86400./dtime * 1.0) ! pour une lecture une fois par jour
       !
       IF (check) THEN
          PRINT*,'interfoce_slab klon, debut, itap, dtime, ijour, &
               &          lmt_pas ', klon, debut, itap, dtime, ijour, &
               &          lmt_pas
       ENDIF !check
       !
       PRINT*, '************************'
       PRINT*, 'SLAB OCEAN est actif, prenez precautions !'
       PRINT*, '************************'
       !
       ! a mettre un slab_bils aussi en force !!!
       !
       DO i = 1, klon
          slab_bils(i) = 0.0
       ENDDO
       !
    ENDIF !debut
    pctsrf_slab(1:klon,1:nbsrf) = pctsrf(1:klon,1:nbsrf)
    !
    ! lecture du bilan au sol lmt_bils issu d'une simulation forcee en debut de journee
    ! 
    IF (MOD(itap,lmt_pas) .EQ. 1) THEN !1er pas de temps de la journee
       idayvrai = ijour
       CALL condsurf(julien,idayvrai, lmt_bils)
    ENDIF !(MOD(itap-1,lmt_pas) .EQ. 0) THEN

    DO i = 1, klon
       IF((pctsrf_slab(i,is_oce).GT.epsfra).OR. &
            &  (pctsrf_slab(i,is_sic).GT.epsfra)) THEN
          !
          ! fabriquer de la glace si congelation atteinte:
          !
          IF (tmp_tslab(i).LT.(RTT-1.8)) THEN
             zz =  (RTT-1.8)-tmp_tslab(i)
             tmp_seaice(i) = tmp_seaice(i) + cyang/cbing * zz
             seaice(i) = tmp_seaice(i)
             tmp_tslab(i) = RTT-1.8
          ENDIF
          !
          ! faire fondre de la glace si temperature est superieure a 0:
          !
          IF ((tmp_tslab(i).GT.RTT) .AND. (tmp_seaice(i).GT.0.0)) THEN
             zz = cyang/cbing * (tmp_tslab(i)-RTT)
             zz = MIN(zz,tmp_seaice(i))
             tmp_seaice(i) = tmp_seaice(i) - zz
             seaice(i) = tmp_seaice(i)
             tmp_tslab(i) = tmp_tslab(i) - zz*cbing/cyang
          ENDIF
          !
          ! limiter la glace de mer a 10 metres (10000 kg/m2)
          !
          IF(tmp_seaice(i).GT.45.) THEN
             tmp_seaice(i) = MIN(tmp_seaice(i),10000.0)
          ELSE
             tmp_seaice(i) = 0. 
          ENDIF
          seaice(i) = tmp_seaice(i)
          siceh(i)=tmp_seaice(i)/1000. !en metres
          !
          ! determiner la nature du sol (glace de mer ou ocean libre):
          !
          ! on fait dependre la fraction de seaice "pctsrf(i,is_sic)" 
          ! de l'epaisseur de seaice :
          ! pctsrf(i,is_sic)=1. si l'epaisseur de la glace de mer est >= a 20cm
          ! et pctsrf(i,is_sic) croit lineairement avec seaice de 0. a 20cm d'epaisseur
          !
          pctsrf_slab(i,is_sic)=MIN(siceh(i)/0.20, &
               &                      1.-(pctsrf_slab(i,is_ter)+pctsrf_slab(i,is_lic)))
          pctsrf_slab(i,is_oce)=1.0 - &
               &      (pctsrf_slab(i,is_ter)+pctsrf_slab(i,is_lic)+pctsrf_slab(i,is_sic))
       ENDIF !pctsrf
    ENDDO
    !
    ! Calculer le bilan du flux de chaleur au sol :
    !
    DO i = 1, klon
       za = radsol(i) + fluxo(i)
       zb = fluxg(i)
       IF((pctsrf_slab(i,is_oce).GT.epsfra).OR. &
            &   (pctsrf_slab(i,is_sic).GT.epsfra)) THEN
          slab_bils(i)=slab_bils(i)+(za*pctsrf_slab(i,is_oce) &
               &             +zb*pctsrf_slab(i,is_sic))/ FLOAT(lmt_pas)
       ENDIF
    ENDDO !klon
    !
    ! calcul tslab 
    !
    IF (MOD(itap,lmt_pas).EQ.0) THEN !fin de journee
       DO i = 1, klon
          IF ((pctsrf_slab(i,is_oce).GT.epsfra).OR. &
               &    (pctsrf_slab(i,is_sic).GT.epsfra)) THEN
             tmp_tslab(i) = tmp_tslab(i) + &
                  & (slab_bils(i)-lmt_bils(i)) &
                  &                         /cyang*unjour
             ! on remet l'accumulation a 0
             slab_bils(i) = 0.
          ENDIF !pctsrf
       ENDDO !klon
    ENDIF !(MOD(itap,lmt_pas).EQ.0) THEN
    !
    tslab = tmp_tslab
    seaice = tmp_seaice
  END SUBROUTINE interfoce_slab

  !************************

  SUBROUTINE interfoce_lim(itime, dtime, jour, &
       & klon, nisurf, knon, knindex, &
       & debut,  &
       & lmt_sst, pctsrf_new)

    ! Cette routine sert d'interface entre le modele atmospherique et un fichier
    ! de conditions aux limites
    !
    ! L. Fairhead 02/2000
    !
    ! input:
    !   itime        numero du pas de temps courant
    !   dtime        pas de temps de la physique (en s)
    !   jour         jour a lire dans l'annee
    !   nisurf       index de la surface a traiter (1 = sol continental)
    !   knon         nombre de points dans le domaine a traiter
    !   knindex      index des points de la surface a traiter
    !   klon         taille de la grille
    !   debut        logical: 1er appel a la physique (initialisation)
    !
    ! output:
    !   lmt_sst      SST lues dans le fichier de CL
    !   pctsrf_new   sous-maille fractionnelle
    !

    use abort_gcm_m, only: abort_gcm
    use indicesol

    ! Parametres d'entree
    integer, intent(IN) :: itime
    real   , intent(IN) :: dtime
    integer, intent(IN) :: jour
    integer, intent(IN) :: nisurf
    integer, intent(IN) :: knon
    integer, intent(IN) :: klon
    integer, dimension(klon), intent(in) :: knindex
    logical, intent(IN) :: debut

    ! Parametres de sortie
    real, intent(out), dimension(klon) :: lmt_sst
    real, intent(out), dimension(klon,nbsrf) :: pctsrf_new

    ! Variables locales
    integer     :: ii
    INTEGER,save :: lmt_pas     ! frequence de lecture des conditions limites 
    ! (en pas de physique)
    logical,save :: deja_lu    ! pour indiquer que le jour a lire a deja
    ! lu pour une surface precedente
    integer,save :: jour_lu 
    integer      :: ierr
    character (len = 20) :: modname = 'interfoce_lim'
    character (len = 80) :: abort_message
    logical, save     :: newlmt = .TRUE.
    logical, save     :: check = .FALSE.
    ! Champs lus dans le fichier de CL
    real, allocatable , save, dimension(:) :: sst_lu, rug_lu, nat_lu
    real, allocatable , save, dimension(:,:) :: pct_tmp
    !
    ! quelques variables pour netcdf
    !
    include "netcdf.inc"
    integer              :: nid, nvarid
    integer, dimension(2) :: start, epais
    !
    ! Fin déclaration
    !

    if (debut .and. .not. allocated(sst_lu)) then
       lmt_pas = nint(86400./dtime * 1.0) ! pour une lecture une fois par jour
       jour_lu = jour - 1
       allocate(sst_lu(klon))
       allocate(nat_lu(klon))
       allocate(pct_tmp(klon,nbsrf))
    endif

    if ((jour - jour_lu) /= 0) deja_lu = .false.

    if (check) write(*,*)modname,' :: jour, jour_lu, deja_lu', jour, jour_lu, &
         deja_lu 
    if (check) write(*,*)modname,' :: itime, lmt_pas ', itime, lmt_pas,dtime

    ! Tester d'abord si c'est le moment de lire le fichier
    if (mod(itime-1, lmt_pas) == 0 .and. .not. deja_lu) then
       !
       ! Ouverture du fichier
       !
       ierr = NF_OPEN ('limit.nc', NF_NOWRITE,nid)
       if (ierr.NE.NF_NOERR) then
          abort_message &
               = 'Pb d''ouverture du fichier de conditions aux limites'
          call abort_gcm(modname,abort_message,1)
       endif
       ! 
       ! La tranche de donnees a lire:
       !
       start(1) = 1
       start(2) = jour
       epais(1) = klon
       epais(2) = 1
       !
       if (newlmt) then
          !
          ! Fraction "ocean" 
          !
          ierr = NF_INQ_VARID(nid, 'FOCE', nvarid)
          if (ierr /= NF_NOERR) then
             abort_message = 'Le champ <FOCE> est absent'
             call abort_gcm(modname,abort_message,1)
          endif
          ierr = NF_GET_VARA_REAL(nid,nvarid,start,epais,pct_tmp(1,is_oce))
          if (ierr /= NF_NOERR) then
             abort_message = 'Lecture echouee pour <FOCE>'
             call abort_gcm(modname,abort_message,1)
          endif
          !
          ! Fraction "glace de mer" 
          !
          ierr = NF_INQ_VARID(nid, 'FSIC', nvarid)
          if (ierr /= NF_NOERR) then
             abort_message = 'Le champ <FSIC> est absent'
             call abort_gcm(modname,abort_message,1)
          endif
          ierr = NF_GET_VARA_REAL(nid,nvarid,start,epais,pct_tmp(1,is_sic))
          if (ierr /= NF_NOERR) then
             abort_message = 'Lecture echouee pour <FSIC>'
             call abort_gcm(modname,abort_message,1)
          endif
          !
          ! Fraction "terre" 
          !
          ierr = NF_INQ_VARID(nid, 'FTER', nvarid)
          if (ierr /= NF_NOERR) then
             abort_message = 'Le champ <FTER> est absent'
             call abort_gcm(modname,abort_message,1)
          endif
          ierr = NF_GET_VARA_REAL(nid,nvarid,start,epais,pct_tmp(1,is_ter))
          if (ierr /= NF_NOERR) then
             abort_message = 'Lecture echouee pour <FTER>'
             call abort_gcm(modname,abort_message,1)
          endif
          !
          ! Fraction "glacier terre" 
          !
          ierr = NF_INQ_VARID(nid, 'FLIC', nvarid)
          if (ierr /= NF_NOERR) then
             abort_message = 'Le champ <FLIC> est absent'
             call abort_gcm(modname,abort_message,1)
          endif
          ierr = NF_GET_VARA_REAL(nid,nvarid,start,epais,pct_tmp(1,is_lic))
          if (ierr /= NF_NOERR) then
             abort_message = 'Lecture echouee pour <FLIC>'
             call abort_gcm(modname,abort_message,1)
          endif
          !
       else  ! on en est toujours a rnatur
          ! 
          ierr = NF_INQ_VARID(nid, 'NAT', nvarid)
          if (ierr /= NF_NOERR) then
             abort_message = 'Le champ <NAT> est absent'
             call abort_gcm(modname,abort_message,1)
          endif
          ierr = NF_GET_VARA_REAL(nid,nvarid,start,epais, nat_lu)
          if (ierr /= NF_NOERR) then
             abort_message = 'Lecture echouee pour <NAT>'
             call abort_gcm(modname,abort_message,1)
          endif
          !
          ! Remplissage des fractions de surface
          ! nat = 0, 1, 2, 3 pour ocean, terre, glacier, seaice
          ! 
          pct_tmp = 0.0
          do ii = 1, klon
             pct_tmp(ii,nint(nat_lu(ii)) + 1) = 1.
          enddo

          !
          !  On se retrouve avec ocean en 1 et terre en 2 alors qu'on veut le contraire
          !
          pctsrf_new = pct_tmp
          pctsrf_new (:,2)= pct_tmp (:,1)
          pctsrf_new (:,1)= pct_tmp (:,2)
          pct_tmp = pctsrf_new 
       endif ! fin test sur newlmt
       !
       ! Lecture SST
       !
       ierr = NF_INQ_VARID(nid, 'SST', nvarid)
       if (ierr /= NF_NOERR) then
          abort_message = 'Le champ <SST> est absent'
          call abort_gcm(modname,abort_message,1)
       endif
       ierr = NF_GET_VARA_REAL(nid,nvarid,start,epais, sst_lu)
       if (ierr /= NF_NOERR) then
          abort_message = 'Lecture echouee pour <SST>'
          call abort_gcm(modname,abort_message,1)
       endif

       !
       ! Fin de lecture
       !
       ierr = NF_CLOSE(nid)
       deja_lu = .true.
       jour_lu = jour
    endif
    !
    ! Recopie des variables dans les champs de sortie
    !
    lmt_sst = 999999999.
    do ii = 1, knon
       lmt_sst(ii) = sst_lu(knindex(ii))
    enddo

    pctsrf_new(:,is_oce) = pct_tmp(:,is_oce)
    pctsrf_new(:,is_sic) = pct_tmp(:,is_sic)

  END SUBROUTINE interfoce_lim

  !************************

  SUBROUTINE interfsur_lim(itime, dtime, jour, &
       & klon, nisurf, knon, knindex, &
       & debut,  &
       & lmt_alb, lmt_rug)

    ! Cette routine sert d'interface entre le modèle atmosphérique et
    ! un fichier de conditions aux limites.
    !
    ! L. Fairhead 02/2000
    !
    ! input:
    !   itime        numero du pas de temps courant
    !   dtime        pas de temps de la physique (en s)
    !   jour         jour a lire dans l'annee
    !   nisurf       index de la surface a traiter (1 = sol continental)
    !   knon         nombre de points dans le domaine a traiter
    !   knindex      index des points de la surface a traiter
    !   klon         taille de la grille
    !   debut        logical: 1er appel a la physique (initialisation)
    !
    ! output:
    !   lmt_sst      SST lues dans le fichier de CL
    !   lmt_alb      Albedo lu 
    !   lmt_rug      longueur de rugosité lue
    !   pctsrf_new   sous-maille fractionnelle
    !

    use abort_gcm_m, only: abort_gcm

    ! Parametres d'entree
    integer, intent(IN) :: itime
    real   , intent(IN) :: dtime
    integer, intent(IN) :: jour
    integer, intent(IN) :: nisurf
    integer, intent(IN) :: knon
    integer, intent(IN) :: klon
    integer, dimension(klon), intent(in) :: knindex
    logical, intent(IN) :: debut

    ! Parametres de sortie
    real, intent(out), dimension(klon) :: lmt_alb
    real, intent(out), dimension(klon) :: lmt_rug

    ! Variables locales
    integer     :: ii
    integer,save :: lmt_pas     ! frequence de lecture des conditions limites 
    ! (en pas de physique)
    logical,save :: deja_lu_sur! pour indiquer que le jour a lire a deja
    ! lu pour une surface precedente
    integer,save :: jour_lu_sur 
    integer      :: ierr
    character (len = 20) :: modname = 'interfsur_lim'
    character (len = 80) :: abort_message
    logical,save     :: newlmt = .false.
    logical,save     :: check = .false.
    ! Champs lus dans le fichier de CL
    real, allocatable , save, dimension(:) :: alb_lu, rug_lu
    !
    ! quelques variables pour netcdf
    !
    include "netcdf.inc"
    integer ,save             :: nid, nvarid
    integer, dimension(2),save :: start, epais
    !
    ! Fin déclaration
    !

    if (debut) then
       lmt_pas = nint(86400./dtime * 1.0) ! pour une lecture une fois par jour
       jour_lu_sur = jour - 1
       allocate(alb_lu(klon))
       allocate(rug_lu(klon))
    endif

    if ((jour - jour_lu_sur) /= 0) deja_lu_sur = .false.

    if (check) write(*,*)modname,':: jour_lu_sur, deja_lu_sur', jour_lu_sur, &
         deja_lu_sur 
    if (check) write(*,*)modname,':: itime, lmt_pas', itime, lmt_pas

    ! Tester d'abord si c'est le moment de lire le fichier
    if (mod(itime-1, lmt_pas) == 0 .and. .not. deja_lu_sur) then
       !
       ! Ouverture du fichier
       !
       ierr = NF_OPEN ('limit.nc', NF_NOWRITE,nid)
       if (ierr.NE.NF_NOERR) then
          abort_message &
               = 'Pb d''ouverture du fichier de conditions aux limites'
          call abort_gcm(modname,abort_message,1)
       endif
       ! 
       ! La tranche de donnees a lire:

       start(1) = 1
       start(2) = jour
       epais(1) = klon
       epais(2) = 1
       !
       ! Lecture Albedo
       !
       ierr = NF_INQ_VARID(nid, 'ALB', nvarid)
       if (ierr /= NF_NOERR) then
          abort_message = 'Le champ <ALB> est absent'
          call abort_gcm(modname,abort_message,1)
       endif
       ierr = NF_GET_VARA_REAL(nid,nvarid,start,epais, alb_lu)
       if (ierr /= NF_NOERR) then
          abort_message = 'Lecture echouee pour <ALB>'
          call abort_gcm(modname,abort_message,1)
       endif
       !
       ! Lecture rugosité
       !
       ierr = NF_INQ_VARID(nid, 'RUG', nvarid)
       if (ierr /= NF_NOERR) then
          abort_message = 'Le champ <RUG> est absent'
          call abort_gcm(modname,abort_message,1)
       endif
       ierr = NF_GET_VARA_REAL(nid,nvarid,start,epais, rug_lu)
       if (ierr /= NF_NOERR) then
          abort_message = 'Lecture echouee pour <RUG>'
          call abort_gcm(modname,abort_message,1)
       endif

       !
       ! Fin de lecture
       !
       ierr = NF_CLOSE(nid)
       deja_lu_sur = .true.
       jour_lu_sur = jour
    endif
    !
    ! Recopie des variables dans les champs de sortie
    !
!!$  lmt_alb(:) = 0.0
!!$  lmt_rug(:) = 0.0
    lmt_alb(:) = 999999.
    lmt_rug(:) = 999999.
    DO ii = 1, knon
       lmt_alb(ii) = alb_lu(knindex(ii))
       lmt_rug(ii) = rug_lu(knindex(ii))
    enddo

  END SUBROUTINE interfsur_lim

  !************************

  SUBROUTINE calcul_fluxs( klon, knon, nisurf, dtime, &
       & tsurf, p1lay, cal, beta, coef1lay, ps, &
       & precip_rain, precip_snow, snow, qsurf, &
       & radsol, dif_grnd, t1lay, q1lay, u1lay, v1lay, &
       & petAcoef, peqAcoef, petBcoef, peqBcoef, &
       & tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l)

    ! Cette routine calcule les fluxs en h et q a l'interface et eventuellement
    ! une temperature de surface (au cas ou ok_veget = false)
    !
    ! L. Fairhead 4/2000
    !
    ! input:
    !   knon         nombre de points a traiter
    !   nisurf       surface a traiter
    !   tsurf        temperature de surface
    !   p1lay        pression 1er niveau (milieu de couche)
    !   cal          capacite calorifique du sol
    !   beta         evap reelle
    !   coef1lay     coefficient d'echange
    !   ps           pression au sol
    !   precip_rain  precipitations liquides
    !   precip_snow  precipitations solides
    !   snow         champs hauteur de neige
    !   runoff       runoff en cas de trop plein
    !   petAcoef     coeff. A de la resolution de la CL pour t
    !   peqAcoef     coeff. A de la resolution de la CL pour q
    !   petBcoef     coeff. B de la resolution de la CL pour t
    !   peqBcoef     coeff. B de la resolution de la CL pour q
    !   radsol       rayonnement net aus sol (LW + SW)
    !   dif_grnd     coeff. diffusion vers le sol profond
    !
    ! output:
    !   tsurf_new    temperature au sol
    !   qsurf        humidite de l'air au dessus du sol
    !   fluxsens     flux de chaleur sensible
    !   fluxlat      flux de chaleur latente
    !   dflux_s      derivee du flux de chaleur sensible / Ts
    !   dflux_l      derivee du flux de chaleur latente  / Ts
    !

    use indicesol
    use abort_gcm_m, only: abort_gcm
    use yoethf
    use fcttre, only: thermcep, foeew, qsats, qsatl, foede, dqsats, dqsatl
    use YOMCST

    ! Parametres d'entree
    integer, intent(IN) :: knon, nisurf, klon
    real   , intent(IN) :: dtime
    real, dimension(klon), intent(IN) :: petAcoef, peqAcoef
    real, dimension(klon), intent(IN) :: petBcoef, peqBcoef
    real, dimension(klon), intent(IN) :: ps, q1lay
    real, dimension(klon), intent(IN) :: tsurf, p1lay, cal, beta, coef1lay
    real, dimension(klon), intent(IN) :: precip_rain, precip_snow
    real, dimension(klon), intent(IN) :: radsol, dif_grnd
    real, dimension(klon), intent(IN) :: t1lay, u1lay, v1lay
    real, dimension(klon), intent(INOUT) :: snow, qsurf

    ! Parametres sorties
    real, dimension(klon), intent(OUT):: tsurf_new, evap, fluxsens, fluxlat
    real, dimension(klon), intent(OUT):: dflux_s, dflux_l

    ! Variables locales
    integer :: i
    real, dimension(klon) :: zx_mh, zx_nh, zx_oh
    real, dimension(klon) :: zx_mq, zx_nq, zx_oq
    real, dimension(klon) :: zx_pkh, zx_dq_s_dt, zx_qsat, zx_coef
    real, dimension(klon) :: zx_sl, zx_k1
    real, dimension(klon) :: zx_q_0 , d_ts
    real                  :: zdelta, zcvm5, zx_qs, zcor, zx_dq_s_dh
    real                  :: bilan_f, fq_fonte
    REAL                  :: subli, fsno
    REAL                  :: qsat_new, q1_new
    real, parameter :: t_grnd = 271.35, t_coup = 273.15
    !! PB temporaire en attendant mieux pour le modele de neige
    REAL, parameter :: chasno = 3.334E+05/(2.3867E+06*0.15)
    !
    logical, save         :: check = .false.
    character (len = 20)  :: modname = 'calcul_fluxs'
    logical, save         :: fonte_neige = .false.
    real, save            :: max_eau_sol = 150.0
    character (len = 80) :: abort_message 
    logical,save         :: first = .true.,second=.false.

    if (check) write(*,*)'Entree ', modname,' surface = ',nisurf

    IF (check) THEN
       WRITE(*,*)' radsol (min, max)' &
            &     , MINVAL(radsol(1:knon)), MAXVAL(radsol(1:knon))
       !!CALL flush(6)
    ENDIF

    if (size(coastalflow) /= knon .AND. nisurf == is_ter) then
       write(*,*)'Bizarre, le nombre de points continentaux'
       write(*,*)'a change entre deux appels. J''arrete ...'
       abort_message='Pb run_off'
       call abort_gcm(modname,abort_message,1)
    endif
    !
    ! Traitement neige et humidite du sol
    !
!!$  WRITE(*,*)'test calcul_flux, surface ', nisurf
    !!PB test
!!$    if (nisurf == is_oce) then
!!$      snow = 0.
!!$      qsol = max_eau_sol
!!$    else
!!$      where (precip_snow > 0.) snow = snow + (precip_snow * dtime)
!!$      where (snow > epsilon(snow)) snow = max(0.0, snow - (evap * dtime))
!!$!      snow = max(0.0, snow + (precip_snow - evap) * dtime)
!!$      where (precip_rain > 0.) qsol = qsol + (precip_rain - evap) * dtime
!!$    endif 
!!$    IF (nisurf /= is_ter) qsol = max_eau_sol

    ! 
    ! Initialisation
    !
    evap = 0.
    fluxsens=0.
    fluxlat=0.
    dflux_s = 0.
    dflux_l = 0.
    !
    ! zx_qs = qsat en kg/kg
    !
    DO i = 1, knon
       zx_pkh(i) = (ps(i)/ps(i))**RKAPPA
       IF (thermcep) THEN
          zdelta=MAX(0.,SIGN(1.,rtt-tsurf(i)))
          zcvm5 = R5LES*RLVTT*(1.-zdelta) + R5IES*RLSTT*zdelta
          zcvm5 = zcvm5 / RCPD / (1.0+RVTMP2*q1lay(i))
          zx_qs= r2es * FOEEW(tsurf(i),zdelta)/ps(i)
          zx_qs=MIN(0.5,zx_qs)
          zcor=1./(1.-retv*zx_qs)
          zx_qs=zx_qs*zcor
          zx_dq_s_dh = FOEDE(tsurf(i),zdelta,zcvm5,zx_qs,zcor) &
               &                 /RLVTT / zx_pkh(i)
       ELSE
          IF (tsurf(i).LT.t_coup) THEN
             zx_qs = qsats(tsurf(i)) / ps(i)
             zx_dq_s_dh = dqsats(tsurf(i),zx_qs)/RLVTT &
                  &                    / zx_pkh(i)
          ELSE
             zx_qs = qsatl(tsurf(i)) / ps(i)
             zx_dq_s_dh = dqsatl(tsurf(i),zx_qs)/RLVTT &
                  &               / zx_pkh(i)
          ENDIF
       ENDIF
       zx_dq_s_dt(i) = RCPD * zx_pkh(i) * zx_dq_s_dh
       zx_qsat(i) = zx_qs
       zx_coef(i) = coef1lay(i) &
            & * (1.0+SQRT(u1lay(i)**2+v1lay(i)**2)) &
            & * p1lay(i)/(RD*t1lay(i))

    ENDDO

    ! === Calcul de la temperature de surface ===
    ! 
    ! zx_sl = chaleur latente d'evaporation ou de sublimation
    !
    do i = 1, knon
       zx_sl(i) = RLVTT
       if (tsurf(i) .LT. RTT) zx_sl(i) = RLSTT
       zx_k1(i) = zx_coef(i)
    enddo

    do i = 1, knon
       ! Q
       zx_oq(i) = 1. - (beta(i) * zx_k1(i) * peqBcoef(i) * dtime)
       zx_mq(i) = beta(i) * zx_k1(i) * &
            &             (peqAcoef(i) - zx_qsat(i) &
            &                          + zx_dq_s_dt(i) * tsurf(i)) &
            &             / zx_oq(i)
       zx_nq(i) = beta(i) * zx_k1(i) * (-1. * zx_dq_s_dt(i)) &
            &                              / zx_oq(i)

       ! H
       zx_oh(i) = 1. - (zx_k1(i) * petBcoef(i) * dtime)
       zx_mh(i) = zx_k1(i) * petAcoef(i) / zx_oh(i)
       zx_nh(i) = - (zx_k1(i) * RCPD * zx_pkh(i))/ zx_oh(i)

       ! Tsurface
       tsurf_new(i) = (tsurf(i) + cal(i)/(RCPD * zx_pkh(i)) * dtime * &
            &             (radsol(i) + zx_mh(i) + zx_sl(i) * zx_mq(i)) & 
            &                 + dif_grnd(i) * t_grnd * dtime)/ &
            &          ( 1. - dtime * cal(i)/(RCPD * zx_pkh(i)) * ( &
            &                       zx_nh(i) + zx_sl(i) * zx_nq(i)) &  
            &                     + dtime * dif_grnd(i))

       !
       ! Y'a-t-il fonte de neige?
       !
       !    fonte_neige = (nisurf /= is_oce) .AND. &
       !     & (snow(i) > epsfra .OR. nisurf == is_sic .OR. nisurf == is_lic) &
       !     & .AND. (tsurf_new(i) >= RTT)
       !    if (fonte_neige) tsurf_new(i) = RTT  
       d_ts(i) = tsurf_new(i) - tsurf(i)
       !    zx_h_ts(i) = tsurf_new(i) * RCPD * zx_pkh(i)
       !    zx_q_0(i) = zx_qsat(i) + zx_dq_s_dt(i) * d_ts(i)
       !== flux_q est le flux de vapeur d'eau: kg/(m**2 s)  positive vers bas
       !== flux_t est le flux de cpt (energie sensible): j/(m**2 s)
       evap(i) = - zx_mq(i) - zx_nq(i) * tsurf_new(i) 
       fluxlat(i) = - evap(i) * zx_sl(i)
       fluxsens(i) = zx_mh(i) + zx_nh(i) * tsurf_new(i)
       ! Derives des flux dF/dTs (W m-2 K-1):
       dflux_s(i) = zx_nh(i)
       dflux_l(i) = (zx_sl(i) * zx_nq(i))
       ! Nouvelle valeure de l'humidite au dessus du sol
       qsat_new=zx_qsat(i) + zx_dq_s_dt(i) * d_ts(i)
       q1_new = peqAcoef(i) - peqBcoef(i)*evap(i)*dtime
       qsurf(i)=q1_new*(1.-beta(i)) + beta(i)*qsat_new
    ENDDO

  END SUBROUTINE calcul_fluxs

  !************************

  SUBROUTINE fonte_neige( klon, knon, nisurf, dtime, &
       & tsurf, p1lay, cal, beta, coef1lay, ps, &
       & precip_rain, precip_snow, snow, qsol, &
       & radsol, dif_grnd, t1lay, q1lay, u1lay, v1lay, &
       & petAcoef, peqAcoef, petBcoef, peqBcoef, &
       & tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l, &
       & fqcalving,ffonte,run_off_lic_0)

    ! Routine de traitement de la fonte de la neige dans le cas du traitement
    ! de sol simplifié
    !
    ! LF 03/2001
    ! input:
    !   knon         nombre de points a traiter
    !   nisurf       surface a traiter
    !   tsurf        temperature de surface
    !   p1lay        pression 1er niveau (milieu de couche)
    !   cal          capacite calorifique du sol
    !   beta         evap reelle
    !   coef1lay     coefficient d'echange
    !   ps           pression au sol
    !   precip_rain  precipitations liquides
    !   precip_snow  precipitations solides
    !   snow         champs hauteur de neige
    !   qsol         hauteur d'eau contenu dans le sol
    !   runoff       runoff en cas de trop plein
    !   petAcoef     coeff. A de la resolution de la CL pour t
    !   peqAcoef     coeff. A de la resolution de la CL pour q
    !   petBcoef     coeff. B de la resolution de la CL pour t
    !   peqBcoef     coeff. B de la resolution de la CL pour q
    !   radsol       rayonnement net aus sol (LW + SW)
    !   dif_grnd     coeff. diffusion vers le sol profond
    !
    ! output:
    !   tsurf_new    temperature au sol
    !   fluxsens     flux de chaleur sensible
    !   fluxlat      flux de chaleur latente
    !   dflux_s      derivee du flux de chaleur sensible / Ts
    !   dflux_l      derivee du flux de chaleur latente  / Ts
    ! in/out:
    !   run_off_lic_0 run off glacier du pas de temps précedent
    !

    use indicesol
    use YOMCST
    use yoethf
    use fcttre
    !IM cf JLD

    ! Parametres d'entree
    integer, intent(IN) :: knon, nisurf, klon
    real   , intent(IN) :: dtime
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
    !  real, parameter :: snow_max=1.
    !IM cf JLD/GK
    real, parameter :: snow_max=3000.
    integer :: i
    real, dimension(klon) :: zx_mh, zx_nh, zx_oh
    real, dimension(klon) :: zx_mq, zx_nq, zx_oq
    real, dimension(klon) :: zx_pkh, zx_dq_s_dt, zx_qsat, zx_coef
    real, dimension(klon) :: zx_sl, zx_k1
    real, dimension(klon) :: zx_q_0 , d_ts
    real                  :: zdelta, zcvm5, zx_qs, zcor, zx_dq_s_dh
    real                  :: bilan_f, fq_fonte
    REAL                  :: subli, fsno
    REAL, DIMENSION(klon) :: bil_eau_s, snow_evap
    real, parameter :: t_grnd = 271.35, t_coup = 273.15
    !! PB temporaire en attendant mieux pour le modele de neige
    ! REAL, parameter :: chasno = RLMLT/(2.3867E+06*0.15)
    REAL, parameter :: chasno = 3.334E+05/(2.3867E+06*0.15)
    !IM cf JLD/ GKtest
    REAL, parameter :: chaice = 3.334E+05/(2.3867E+06*0.15)
    ! fin GKtest
    !
    logical, save         :: check = .FALSE.
    character (len = 20)  :: modname = 'fonte_neige'
    logical, save         :: neige_fond = .false.
    real, save            :: max_eau_sol = 150.0
    character (len = 80) :: abort_message 
    logical,save         :: first = .true.,second=.false.
    real                 :: coeff_rel

    if (check) write(*,*)'Entree ', modname,' surface = ',nisurf

    ! Initialisations
    coeff_rel = dtime/(tau_calv * rday)
    bil_eau_s(:) = 0.
    DO i = 1, knon
       zx_pkh(i) = (ps(i)/ps(i))**RKAPPA
       IF (thermcep) THEN
          zdelta=MAX(0.,SIGN(1.,rtt-tsurf(i)))
          zcvm5 = R5LES*RLVTT*(1.-zdelta) + R5IES*RLSTT*zdelta
          zcvm5 = zcvm5 / RCPD / (1.0+RVTMP2*q1lay(i))
          zx_qs= r2es * FOEEW(tsurf(i),zdelta)/ps(i)
          zx_qs=MIN(0.5,zx_qs)
          zcor=1./(1.-retv*zx_qs)
          zx_qs=zx_qs*zcor
          zx_dq_s_dh = FOEDE(tsurf(i),zdelta,zcvm5,zx_qs,zcor) &
               &                 /RLVTT / zx_pkh(i)
       ELSE
          IF (tsurf(i).LT.t_coup) THEN
             zx_qs = qsats(tsurf(i)) / ps(i)
             zx_dq_s_dh = dqsats(tsurf(i),zx_qs)/RLVTT &
                  &                    / zx_pkh(i)
          ELSE
             zx_qs = qsatl(tsurf(i)) / ps(i)
             zx_dq_s_dh = dqsatl(tsurf(i),zx_qs)/RLVTT &
                  &               / zx_pkh(i)
          ENDIF
       ENDIF
       zx_dq_s_dt(i) = RCPD * zx_pkh(i) * zx_dq_s_dh
       zx_qsat(i) = zx_qs
       zx_coef(i) = coef1lay(i) &
            & * (1.0+SQRT(u1lay(i)**2+v1lay(i)**2)) &
            & * p1lay(i)/(RD*t1lay(i))
    ENDDO

    ! === Calcul de la temperature de surface ===
    ! 
    ! zx_sl = chaleur latente d'evaporation ou de sublimation
    !
    do i = 1, knon
       zx_sl(i) = RLVTT
       if (tsurf(i) .LT. RTT) zx_sl(i) = RLSTT
       zx_k1(i) = zx_coef(i)
    enddo

    do i = 1, knon
       ! Q
       zx_oq(i) = 1. - (beta(i) * zx_k1(i) * peqBcoef(i) * dtime)
       zx_mq(i) = beta(i) * zx_k1(i) * &
            &             (peqAcoef(i) - zx_qsat(i) &
            &                          + zx_dq_s_dt(i) * tsurf(i)) &
            &             / zx_oq(i)
       zx_nq(i) = beta(i) * zx_k1(i) * (-1. * zx_dq_s_dt(i)) &
            &                              / zx_oq(i)

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

    !  bil_eau_s = bil_eau_s + (precip_rain * dtime) - (evap - snow_evap) * dtime
    bil_eau_s = (precip_rain * dtime) - (evap - snow_evap) * dtime

    !
    ! Y'a-t-il fonte de neige?
    !
    ffonte=0.
    do i = 1, knon
       neige_fond = ((snow(i) > epsfra .OR. nisurf == is_sic .OR. nisurf == is_lic) &
            & .AND. tsurf_new(i) >= RTT)
       if (neige_fond) then
          fq_fonte = MIN( MAX((tsurf_new(i)-RTT )/chasno,0.0),snow(i))
          ffonte(i) = fq_fonte * RLMLT/dtime
          snow(i) = max(0., snow(i) - fq_fonte)
          bil_eau_s(i) = bil_eau_s(i) + fq_fonte 
          tsurf_new(i) = tsurf_new(i) - fq_fonte * chasno  
          !IM cf JLD OK     
          !IM cf JLD/ GKtest fonte aussi pour la glace
          IF (nisurf == is_sic .OR. nisurf == is_lic ) THEN
             fq_fonte = MAX((tsurf_new(i)-RTT )/chaice,0.0)
             ffonte(i) = ffonte(i) + fq_fonte * RLMLT/dtime
             bil_eau_s(i) = bil_eau_s(i) + fq_fonte
             tsurf_new(i) = RTT
          ENDIF
          d_ts(i) = tsurf_new(i) - tsurf(i)
       endif
       !
       !   s'il y a une hauteur trop importante de neige, elle s'coule
       fqcalving(i) = max(0., snow(i) - snow_max)/dtime
       snow(i)=min(snow(i),snow_max)
       !
       IF (nisurf == is_ter) then
          qsol(i) = qsol(i) + bil_eau_s(i)
          run_off(i) = run_off(i) + MAX(qsol(i) - max_eau_sol, 0.0)
          qsol(i) = MIN(qsol(i), max_eau_sol) 
       else if (nisurf == is_lic) then
          run_off_lic(i) = (coeff_rel *  fqcalving(i)) + &
               &                        (1. - coeff_rel) * run_off_lic_0(i)
          run_off_lic_0(i) = run_off_lic(i)
          run_off_lic(i) = run_off_lic(i) + bil_eau_s(i)/dtime
       endif
    enddo

  END SUBROUTINE fonte_neige

END MODULE interface_surf
