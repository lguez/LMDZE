PROGRAM gcm

  ! Authors: P. Le Van, L. Fairhead, F. Hourdin
  ! From "gcm.F", version 1.4, 2006/04/04 15:05:16

  ! General circulation model of LMD. Avec coordonnée verticale
  ! hybride, avec nouveaux opérateurs de dissipation "*" ("gradiv2",
  ! "divgrad2", "nxgraro2"). Possibilité de choisir le schéma pour
  ! l'advection de "q", en modifiant "iadv" dans "traceur.def".

  use comconst, only: daysec, dtvr, iniconst
  use comgeom, only: rlatu, aire_2d, cu_2d, cv_2d, rlonv, inigeom
  use comgeomphy, only: airephy, cuphy, cvphy, rlatd, rlond
  use conf_gcm_m, only: day_step, iperiod, anneeref, dayref, iecri, iphysiq, &
       nday, raz_date, periodav, conf_gcm, iflag_phys
  use conf_guide_m, only: conf_guide
  use dimens_m, only: iim, jjm, llm, nqmx
  use dimphy, only: klon
  USE disvert_m, ONLY : disvert
  use dynetat0_m, only: dynetat0, day_ini
  use dynredem0_m, only: dynredem0
  use grid_change, only: dyn_phy, init_dyn_phy
  use histclo_m, only: histclo
  use iniadvtrac_m, only: iniadvtrac
  use inidissip_m, only: inidissip
  use inifilr_m, only: inifilr
  use initdynav_m, only: initdynav
  use inithist_m, only: inithist
  use init_dynzon_m, only: init_dynzon
  USE ioconf_calendar_m, only: ioconf_calendar
  use jumble, only: new_unit
  use leapfrog_m, only: leapfrog
  use netcdf95, only: nf95_close
  use suphec_m, only: suphec
  use temps, only: day_ref, annee_ref, day_end, itau_dyn
  use tracstoke, only: istdyn, istphy
  use unit_nml_m, only: unit_nml
  use yoethf_m, only: yoethf
  use createnewfield_m, only: NbField, Ncid

  IMPLICIT NONE

  ! Variables dynamiques :
  REAL ucov(iim + 1, jjm + 1, llm), vcov(iim + 1, jjm, llm)  ! vent covariant
  REAL teta(iim + 1, jjm + 1, llm) ! température potentielle 
  REAL q(iim + 1, jjm + 1, llm, nqmx) ! champs advectés
  REAL ps(iim + 1, jjm + 1) ! pression au sol (Pa)
  REAL masse(iim + 1, jjm + 1, llm) ! masse d'air
  REAL phis(iim + 1, jjm + 1) ! géopotentiel au sol

  ! Variables pour le fichier histoire :
  REAL time_0 ! time in day, as a fraction of day, in [0, 1[

  ! Calendrier :
  LOGICAL:: true_calendar = .false. ! default value

  logical mask_v(iim + 1, jjm) 
  ! (mask for points in the "v" grid, first index is for longitude,
  ! second index is for latitude)

  integer i

  namelist /main_nml/true_calendar

  !------------------------------------------------------------

  call new_unit(unit_nml)
  open(unit_nml, file="used_namelists.txt", status="replace", action="write")

  CALL conf_gcm

  print *, "Enter namelist 'main_nml'."
  read (unit=*, nml=main_nml)
  write(unit_nml, nml=main_nml)

  ! Choix du calendrier :
  if (true_calendar) then
     call ioconf_calendar('gregorian')
  else
     call ioconf_calendar('360d')
  endif

  ! Initialisation des traceurs
  ! Choix du schéma pour l'advection dans le fichier "traceur.def" ou via INCA
  call iniadvtrac

  CALL iniconst

  ! Lecture du fichier "start.nc" :
  CALL dynetat0(vcov, ucov, teta, q, masse, ps, phis, time_0)

  ! On remet le calendrier à zéro si demandé :
  if (raz_date) then
     print *, 'On réinitialise à la date lue dans la namelist.'
     annee_ref = anneeref
     day_ref = dayref
     day_ini = dayref
     itau_dyn = 0
     time_0 = 0.
  else
     print *, 'On garde les dates du fichier "start".'
  endif

  CALL disvert
  CALL inigeom ! initialisation de la géometrie
  CALL inifilr ! initialisation du filtre
  CALL inidissip
  call init_dyn_phy

  ! Initialisation de la physique :
  IF (iflag_phys == 1) THEN
     rlatd(1)=rlatu(1)
     rlatd(2:klon-1) = pack(spread(rlatu(2:jjm), 1, iim), .true.)
     rlatd(klon)= rlatu(jjm + 1)

     rlond(1)=0.
     rlond(2:klon-1) = pack(spread(rlonv(:iim), 2, jjm - 1), .true.)
     rlond(klon)= 0.

     cuphy = pack(cu_2d, dyn_phy)

     ! Construct a mask for points in the "v" grid:
     mask_v = .true.
     mask_v(2:, 1) = .false.
     mask_v(iim + 1, 2:) = .false.

     cvphy(:klon - 1) = pack(cv_2d, mask_v)
     cvphy(klon) = cv_2d(1, jjm)
     ! (that value of "cv_2d" is used twice in "cvphy")

     airephy = pack(aire_2d, dyn_phy)
     CALL suphec
     call yoethf
  ENDIF

  ! Initialisation des entrées-sorties :
  day_end = day_ini + nday
  print *, "day_ini = ", day_ini
  print *, "day_end = ", day_end

  CALL dynredem0("restart.nc", day_end, phis)
  CALL inithist(day_ref, annee_ref, dtvr, nqmx, t_ops = iecri * daysec, &
       t_wrt = iecri * daysec)
  CALL initdynav(day_ref, annee_ref, dtvr, nqmx, t_ops = iperiod * dtvr, &
       t_wrt = periodav * daysec)
  call init_dynzon(dt_app = dtvr * iperiod)

  ! Choix des fréquences de stockage pour le hors-ligne :
  istdyn = day_step / 4 ! stockage toutes les 6 h = 1 jour / 4
  istphy = istdyn / iphysiq     

  CALL conf_guide

  ! Intégration temporelle du modèle :
  CALL leapfrog(ucov, vcov, teta, ps, masse, phis, q, time_0)

  close(unit_nml)
  call histclo

  do i = 1, nbfield
     call nf95_close(Ncid(i))
  end do

  print *, 'Simulation finished'
  print *, 'Everything is cool'

END PROGRAM gcm
