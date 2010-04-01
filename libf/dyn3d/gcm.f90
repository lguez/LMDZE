PROGRAM gcm

  ! General circulation model of LMD

  ! From "gcm.F", version 1.4, 2006/04/04 15:05:16

  ! Avec coordonnées verticales hybrides, avec nouveaux opérateurs de
  ! dissipation * (gradiv2, divgrad2, nxgraro2)

  ! Authors: P. Le Van, L. Fairhead, F. Hourdin

  ! Possibilité de choisir le schéma pour l'advection de "q", en
  ! modifiant "iadv" dans "traceur.def".

  ! Pour Van-Leer plus vapeur d'eau saturée : iadv(1)=4
  ! Pour Van-Leer : iadv=10

  use abort_gcm_m, only: abort_gcm
  use clesphys2, only: read_clesphys2
  use com_io_dyn, only: histid, histvid, histaveid
  use comconst, only: daysec, cpp, dtvr, g, rad, r, initialize
  use comgeom, only: rlatu, aire_2d, cu_2d, cv_2d, rlonv
  use conf_gcm_m, only: day_step, iperiod, anneeref, dayref, iecri, iphysiq, &
       nday, raz_date, periodav, conf_gcm
  use dimens_m, only: iim, jjm, llm, nqmx
  use dimphy, only: klon
  use dynetat0_m, only: dynetat0, day_ini
  use dynredem0_m, only: dynredem0
  use grid_change, only: dyn_phy, init_dyn_phy
  use iniadvtrac_m, only: iniadvtrac
  use inidissip_m, only: inidissip
  use inigeom_m, only: inigeom
  use initdynav_m, only: initdynav
  use inithist_m, only: inithist
  USE calendar, only: ioconf_calendar
  use histcom, only: histclo
  use leapfrog_m, only: leapfrog
  use logic, only: iflag_phys
  use paramet_m, only: ip1jm, ip1jmp1
  use temps, only: day_ref, annee_ref, day_end, itau_dyn
  use tracstoke, only: istdyn, istphy

  IMPLICIT NONE

  REAL zdtvr ! time step for dynamics, in s

  ! Variables dynamiques :
  REAL vcov(ip1jm, llm), ucov(ip1jmp1, llm) ! vents covariants
  REAL teta(ip1jmp1, llm) ! température potentielle 
  REAL q(ip1jmp1, llm, nqmx) ! champs advectés
  REAL ps(ip1jmp1) ! pression au sol (Pa)

  REAL masse(ip1jmp1, llm) ! masse d'air
  REAL phis(iim + 1, jjm + 1) ! géopotentiel au sol

  ! Variables pour le fichier histoire :
  REAL time_0 ! time in day, as a fraction of day, in [0, 1[

  ! Calendrier :
  LOGICAL:: true_calendar = .false. ! default value

  ! Variables pour l'initialisation de la physique :
  REAL zcufi(klon), zcvfi(klon) ! "cu" and "cv" values on the scalar grid
  REAL latfi(klon), lonfi(klon)
  REAL airefi(klon)

  logical mask_v(iim + 1, jjm) 
  ! (mask for points in the "v" grid, first index is for longitude,
  ! second index is for latitude)

  namelist /main_nml/true_calendar

  !------------------------------------------------------------

  print *, "Enter namelist 'main_nml'."
  read (unit=*, nml=main_nml)
  write(unit=*, nml=main_nml)

  ! Initialisations:
  call initialize

  ! Choix du calendrier :
  if (true_calendar) then
     call ioconf_calendar('gregorian')
  else
     call ioconf_calendar('360d')
  endif

  ! Lecture des fichiers "gcm.def" ou "run.def" :
  call read_clesphys2
  CALL conf_gcm

  ! Initialisation des traceurs
  ! Choix du schéma pour l'advection dans le fichier "traceur.def" ou via INCA
  call iniadvtrac

  ! Lecture du fichier "start.nc" :
  CALL dynetat0(vcov, ucov, teta, q, masse, ps, phis, time_0)

  ! Lecture des paramètres de contrôle pour la simulation :
  ! on recalcule éventuellement le pas de temps
  IF (MOD(day_step, iperiod) /= 0) THEN
     call abort_gcm(modname = "gcm", message = &
          'Il faut choisir un nombre de pas par jour multiple de "iperiod".', &
          ierr = 1)
  ENDIF

  IF (MOD(day_step,iphysiq)/=0) THEN
     call abort_gcm(modname = "gcm", message = &
          'Il faut choisir un nombre de pas par jour multiple de "iphysiq".', &
          ierr = 1)
  ENDIF

  ! On remet le calendrier à zero si demandé:
  if (annee_ref /= anneeref .or. day_ref /= dayref) then
     print *, 'Attention : les dates initiales lues dans le fichier ' // &
          '"start" ne correspondent pas à celles lues dans "gcm.def".'
     if (raz_date) then
        print *, 'On réinitialise à la date lue dans "gcm.def".'
        annee_ref = anneeref
        day_ref = dayref
        day_ini = dayref
        itau_dyn = 0
        time_0 = 0.
     else
        print *, 'On garde les dates du fichier "start".'
     endif
  ELSE
     raz_date = .false.
  endif

  ! Initialisation des constantes dynamiques :
  zdtvr = daysec / REAL(day_step)
  IF (dtvr /= zdtvr) THEN
     print *, 'Warning: the time steps in the ".def" file and in ' // &
          '"start.nc" are different'
     print *, 'dtvr (from "start.nc") = ', dtvr
     print *, 'zdtvr (from ".def") = ', zdtvr
     print *, 'Using the value from the ".def" file.'
     dtvr = zdtvr
  ENDIF
  CALL iniconst

  ! Initialisation de la géometrie :
  CALL inigeom

  ! Initialisation du filtre :
  CALL inifilr

  CALL inidissip
  call init_dyn_phy

  ! Initialisation de la physique :
  IF (iflag_phys == 1) THEN
     latfi(1)=rlatu(1)
     latfi(2:klon-1) = pack(spread(rlatu(2:jjm), 1, iim), .true.)
     latfi(klon)= rlatu(jjm + 1)

     lonfi(1)=0.
     lonfi(2:klon-1) = pack(spread(rlonv(:iim), 2, jjm - 1), .true.)
     lonfi(klon)= 0.

     zcufi = pack(cu_2d, dyn_phy)

     ! Construct a mask for points in the "v" grid:
     mask_v = .true.
     mask_v(2:, 1) = .false.
     mask_v(iim + 1, 2:) = .false.

     zcvfi(:klon - 1) = pack(cv_2d, mask_v)
     zcvfi(klon) = cv_2d(1, jjm)
     ! (that value of "cv_2d" is used twice in "zcvfi")

     airefi = pack(aire_2d, dyn_phy)
     print *, 'Attention : vitesse verticale nulle dans la physique.'
     CALL iniphysiq(klon, llm, latfi, lonfi, airefi, &
          zcufi, zcvfi, rad, g, r, cpp)
  ENDIF

  ! Initialisation des entrées-sorties :
  day_end = day_ini + nday
  print *, "day_ini = ", day_ini
  print *, "day_end = ", day_end

  CALL dynredem0("restart.nc", day_end, phis)
  CALL inithist(day_ref, annee_ref, zdtvr, nqmx, histid, histvid, &
       t_ops = iecri * daysec, t_wrt = iecri * daysec)
  CALL initdynav(day_ref, annee_ref, zdtvr, nqmx, histaveid, &
       t_ops = iperiod * zdtvr, t_wrt = periodav * daysec)

  ! Choix des fréquences de stockage pour le hors-ligne :
  istdyn = day_step / 4     ! stockage toutes les 6 h = 1 jour / 4
  istphy = istdyn / iphysiq     

  ! Intégration temporelle du modèle :
  CALL leapfrog(ucov, vcov, teta, ps, masse, phis, q, time_0)

  call histclo
  print *, 'Simulation finished'
  print *, 'Everything is cool'

END PROGRAM gcm
