PROGRAM gcm

  ! Authors: P. Le Van, L. Fairhead, F. Hourdin
  ! From "gcm.F", version 1.4, 2006/04/04

  ! General circulation model of LMD. Avec coordonn\'ee verticale
  ! hybride, avec nouveaux op\'erateurs de dissipation "*" ("gradiv2",
  ! "divgrad2", "nxgraro2"). Possibilit\'e de choisir le sch\'ema pour
  ! l'advection de "q", en modifiant "iadv" dans "traceur.def".

  ! Libraries:
  use mpi_f08, only: mpi_init, mpi_finalize, mpi_comm_size, mpi_comm_world, &
       mpi_abort
  use netcdf, only: NF90_NOWRITE
  use netcdf95, only: nf95_close, nf95_open
  use xios, only: xios_initialize, xios_finalize, xios_context_initialize, &
       xios_context_finalize, xios_close_context_definition, &
       xios_set_time_origin, xios_date

  use caldyn0_m, only: read_caldyn0
  use comdissnew, only: read_comdissnew
  use comgeom, only:  inigeom
  use conf_gcm_m, only: day_step, iperiod, iphysiq, nday, conf_gcm, iflag_phys
  use conf_guide_m, only: conf_guide
  use dimensions, only: iim, jjm, llm, nqmx, set_dimensions
  use dimphy, only: init_dimphy
  USE disvert_m, ONLY : disvert
  use dynetat0_m, only: dynetat0, day_ini
  use dynetat0_chosen_m, only: dynetat0_chosen, annee_ref, day_ref
  use dynredem0_m, only: dynredem0
  use grid_change, only: init_dyn_phy
  use grid_noro_m, only: read_phis
  use histclo_m, only: histclo
  use infotrac_init_m, only: infotrac_init
  use inidissip_m, only: inidissip
  use inifilr_m, only: inifilr
  use inithist_m, only: inithist
  use init_dynzon_m, only: init_dynzon
  USE ioconf_calendar_m, only: ioconf_calendar
  use leapfrog_m, only: leapfrog
  use paramet_m, only: paramet
  use suphec_m, only: suphec
  use unit_nml_m, only: unit_nml, set_unit_nml
  use createnewfield_m, only: NbField, Ncid

  IMPLICIT NONE

  REAL, ALLOCATABLE:: ucov(:, :, :) ! (iim + 1, jjm + 1, llm) ! vent covariant
  REAL, ALLOCATABLE:: vcov(:, :, :) ! (iim + 1, jjm, llm) ! vent covariant

  REAL, ALLOCATABLE:: teta(:, :, :) ! (iim + 1, jjm + 1, llm)
  ! temp\'erature potentielle 

  REAL, ALLOCATABLE:: q(:, :, :, :) ! (iim + 1, jjm + 1, llm, nqmx)
  ! mass fraction of advected species

  REAL, ALLOCATABLE:: ps(:, :) ! (iim + 1, jjm + 1) ! pression au sol (Pa)

  REAL, ALLOCATABLE:: masse(:, :, :) ! (iim + 1, jjm + 1, llm)
  ! mass in a grid cell, in kg

  integer i, n_proc, return_comm, ncid_start

  !------------------------------------------------------------

  call mpi_init
  call mpi_comm_size(mpi_comm_world, n_proc)
  if (n_proc /= 1) call mpi_abort(mpi_comm_world, errorcode = 1)
  call xios_initialize("LMDZE", return_comm = return_comm)
  CALL xios_context_initialize("LMDZE_context", return_comm)
  call set_unit_nml
  open(unit_nml, file="used_namelists.txt", status="replace", action="write")
  call set_dimensions
  ALLOCATE(ucov(iim + 1, jjm + 1, llm), vcov(iim + 1, jjm, llm))
  ALLOCATE(teta(iim + 1, jjm + 1, llm))
  ALLOCATE(q(iim + 1, jjm + 1, llm, nqmx))
  ALLOCATE(ps(iim + 1, jjm + 1))
  ALLOCATE(masse(iim + 1, jjm + 1, llm))
  call paramet
  call init_dimphy
  CALL conf_gcm
  call read_comdissnew
  call ioconf_calendar('360_day')
  call infotrac_init
  call nf95_open("start.nc", NF90_NOWRITE, ncid_start) ! fichier \'etat initial
  CALL dynetat0_chosen(ncid_start)
  CALL dynetat0(vcov, ucov, teta, q, masse, ps, ncid_start)
  call xios_set_time_origin(xios_date(annee_ref, 1,day_ref, 0, 0, 0))
  call xios_close_context_definition
  call read_phis(ncid_start)
  call read_caldyn0(ncid_start)
  call NF95_CLOSE(ncid_start)
  CALL disvert
  CALL inigeom ! initialisation de la g\'eometrie
  CALL inifilr ! initialisation du filtre
  CALL inidissip
  call init_dyn_phy

  ! Initialisation de la physique :
  IF (iflag_phys) CALL suphec

  ! Initialisation des entr\'ees-sorties :
  CALL dynredem0(iday_end = day_ini + nday)
  CALL inithist
  call init_dynzon

  CALL conf_guide
  CALL leapfrog(ucov, vcov, teta, ps, masse, q)
  close(unit_nml)
  call histclo

  do i = 1, nbfield
     call nf95_close(Ncid(i))
  end do

  print *, 'Simulation finished'
  print *, 'Everything is cool'
  CALL xios_context_finalize
  call xios_finalize
  call mpi_finalize

END PROGRAM gcm
