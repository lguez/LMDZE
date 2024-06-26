module conf_gcm_m

  IMPLICIT NONE

  INTEGER, protected:: nday = 1 ! nombre de jours d'int\'egration

  integer, protected:: day_step = 240
  ! nombre de pas de temps de la dynamique par jour

  integer, protected:: iperiod = 5
  ! periode pour le pas Matsuno (en pas de temps)

  integer, protected:: iapp_tracvl = 5
  ! Should normally be equal to "iperiod"
  ! frequence du groupement des flux (en pas de temps) 

  integer, protected:: iconser = 240
  ! number of time steps of dynamics between output of control variables

  integer, protected:: iphysiq = 5
  ! Number of time steps of dynamics between calls to physics. Must be
  ! >= 1. 1 means one call to physics at each time step of dynamics.

  logical, protected:: raz_date = .false.
  ! prise en compte de la date initiale de la namelist et remise \`a
  ! z\'ero des compteurs de pas de temps (sinon on garde la date du
  ! fichier restart)

  integer, protected:: periodav = 1 
  ! time interval between outputs in the dynamical part, in days

  LOGICAL, protected:: purmats= .FALSE.
  ! Help = Choix du schema d'integration temporel.
  ! y = pure Matsuno sinon c'est du Matsuno-leapfrog

  logical, protected:: iflag_phys = .true. ! call parameterizations of physics
  INTEGER, SAVE, protected:: lmt_pas ! number of time steps of "physics" per day
  INTEGER, protected:: ngroup = 3
  REAL, protected:: dtvr ! time step for dynamics, in s
  REAL, protected:: dtphys ! time step for physics, in s

contains

  SUBROUTINE conf_gcm

    ! Auteurs : L. Fairhead, P. Le Van

    use jumble, only: assert

    use abort_gcm_m, only: abort_gcm
    use comconst, only: daysec
    use dimensions, only: iim, jjm
    use unit_nml_m, only: unit_nml

    ! Local:
    namelist /conf_gcm_nml/ raz_date, nday, day_step, iperiod, iapp_tracvl, &
         iconser, periodav, iphysiq, ngroup
    namelist /logic_nml/ purmats, iflag_phys

    !------------------------------------

    print *, "Enter namelist 'logic_nml'."
    read(unit=*, nml=logic_nml)
    write(unit_nml, nml=logic_nml)

    print *, "Enter namelist 'conf_gcm_nml'."
    read(unit=*, nml=conf_gcm_nml)
    write(unit_nml, nml=conf_gcm_nml)

    if (iphysiq <= 0 .or. iperiod <= 0 .or. day_step <= 0) &
         call abort_gcm("conf_gcm", &
         "iphysiq <= 0 or iperiod <= 0 or day_step <= 0")
    IF (MOD(day_step, iperiod) /= 0) call abort_gcm("conf_gcm", &
         'day_step must be multiple of iperiod.')
    IF (MOD(day_step, iphysiq)/= 0) call abort_gcm("conf_gcm", &
         'day_step must be multiple of iphysiq.')
    call assert(mod(iphysiq, iperiod) == 0, &
         "conf_gcm -- iphysiq must be multiple of iperiod")

    lmt_pas = day_step / iphysiq
    print *, 'Number of time steps of "physics" per day: ', lmt_pas

    call assert(mod(iim, 2**ngroup) == 0 .and. 2**ngroup <= jjm + 1, &
         "conf_gcm: ngroup")

    dtvr = daysec / real(day_step)
    dtphys  = iphysiq * dtvr
    print *, 'dtvr = ', dtvr
    print *, 'dtphys = ', dtphys

  END SUBROUTINE conf_gcm

end module conf_gcm_m
