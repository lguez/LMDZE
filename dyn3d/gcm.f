PROGRAM gcm

  ! Authors: P. Le Van, L. Fairhead, F. Hourdin
  ! From "gcm.F", version 1.4, 2006/04/04 15:05:16

  ! General circulation model of LMD. Avec coordonn\'ee verticale
  ! hybride, avec nouveaux op\'erateurs de dissipation "*" ("gradiv2",
  ! "divgrad2", "nxgraro2"). Possibilit\'e de choisir le sch\'ema pour
  ! l'advection de "q", en modifiant "iadv" dans "traceur.def".

  use comconst, only: dtvr, iniconst
  use comdissnew, only: read_comdissnew
  use comgeom, only:  aire_2d, inigeom
  use comgeomphy, only: airephy
  use conf_gcm_m, only: day_step, iperiod, iecri, iphysiq, nday, conf_gcm, &
       iflag_phys
  use conf_guide_m, only: conf_guide
  use dimensions, only: iim, jjm, llm, nqmx
  USE disvert_m, ONLY : disvert
  use dynetat0_m, only: dynetat0, day_ini
  use dynredem0_m, only: dynredem0
  use grid_change, only: dyn_phy, init_dyn_phy
  use histclo_m, only: histclo
  use iniadvtrac_m, only: iniadvtrac
  use inidissip_m, only: inidissip
  use inifilr_m, only: inifilr
  use inithist_m, only: inithist
  use init_dynzon_m, only: init_dynzon
  USE ioconf_calendar_m, only: ioconf_calendar
  use jumble, only: new_unit
  use leapfrog_m, only: leapfrog
  use netcdf95, only: nf95_close
  use suphec_m, only: suphec
  use unit_nml_m, only: unit_nml
  use yoethf_m, only: yoethf
  use createnewfield_m, only: NbField, Ncid

  IMPLICIT NONE

  ! Variables dynamiques :
  REAL ucov(iim + 1, jjm + 1, llm), vcov(iim + 1, jjm, llm)  ! vent covariant
  REAL teta(iim + 1, jjm + 1, llm) ! temp\'erature potentielle 
  REAL q(iim + 1, jjm + 1, llm, nqmx) ! champs advect\'es
  REAL ps(iim + 1, jjm + 1) ! pression au sol (Pa)
  REAL masse(iim + 1, jjm + 1, llm) ! masse d'air
  REAL phis(iim + 1, jjm + 1) ! g\'eopotentiel au sol

  LOGICAL:: true_calendar = .false. ! default value
  integer i

  namelist /main_nml/true_calendar

  !------------------------------------------------------------

  call new_unit(unit_nml)
  open(unit_nml, file="used_namelists.txt", status="replace", action="write")

  CALL conf_gcm
  call read_comdissnew

  print *, "Enter namelist 'main_nml'."
  read (unit=*, nml=main_nml)
  write(unit_nml, nml=main_nml)

  ! Choix du calendrier :
  if (true_calendar) then
     call ioconf_calendar('gregorian')
  else
     call ioconf_calendar('360d')
  endif

  call iniadvtrac
  CALL iniconst
  CALL dynetat0(vcov, ucov, teta, q, masse, ps, phis)
  CALL disvert
  CALL inigeom ! initialisation de la g\'eometrie
  CALL inifilr ! initialisation du filtre
  CALL inidissip
  call init_dyn_phy

  ! Initialisation de la physique :
  IF (iflag_phys) THEN
     airephy = pack(aire_2d, dyn_phy)
     CALL suphec
     call yoethf
  ENDIF

  ! Initialisation des entr\'ees-sorties :
  CALL dynredem0(day_ini + nday, phis)
  CALL inithist(t_ops = iecri * dtvr, t_wrt = iecri * dtvr)
  call init_dynzon(dt_app = dtvr * iperiod)

  CALL conf_guide
  CALL leapfrog(ucov, vcov, teta, ps, masse, phis, q)

  close(unit_nml)
  call histclo

  do i = 1, nbfield
     call nf95_close(Ncid(i))
  end do

  print *, 'Simulation finished'
  print *, 'Everything is cool'

END PROGRAM gcm
