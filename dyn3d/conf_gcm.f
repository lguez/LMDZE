module conf_gcm_m

  IMPLICIT NONE

  INTEGER:: nday = 1 ! nombre de jours d'int\'egration

  integer, protected:: day_step = 240
  ! nombre de pas de temps de la dynamique par jour

  integer, protected:: iperiod = 5
  ! periode pour le pas Matsuno (en pas de temps)

  integer, protected:: iapp_tracvl = 5
  ! Should normally be equal to "iperiod"
  ! frequence du groupement des flux (en pas de temps) 

  integer, protected:: iconser = 240
  ! number of time steps between output of control variables

  integer, protected:: iphysiq = 5
  ! number of time steps of dynamics between calls to physics

  logical, protected:: raz_date = .false.
  ! prise en compte de la date initiale de la namelist et remise \`a
  ! z\'ero des compteurs de pas de temps (sinon on garde la date du
  ! fichier restart)

  integer, protected:: periodav = 1 
  ! time interval between outputs in the dynamical part, in days

  integer, protected:: prt_level = 0
  ! niveau d'impression souhait\'e (0 = minimum)

  LOGICAL, protected:: purmats= .FALSE.
  ! Help = Choix du schema d'integration temporel.
  ! y = pure Matsuno sinon c'est du Matsuno-leapfrog

  logical, protected:: iflag_phys = .true. ! call parameterizations of physics
  INTEGER, SAVE, protected:: lmt_pas ! number of time steps of "physics" per day

contains

  SUBROUTINE conf_gcm

    ! Auteurs : L. Fairhead, P. Le Van
    ! Version du 29/04/97

    use abort_gcm_m, only: abort_gcm
    use nr_util, only: assert
    use unit_nml_m, only: unit_nml

    namelist /conf_gcm_nml/ raz_date, nday, day_step, iperiod, iapp_tracvl, &
         iconser, periodav, iphysiq

    namelist /iniprint_nml/ prt_level

    namelist /logic_nml/ purmats, iflag_phys

    !------------------------------------

    print *, "Call sequence information: conf_gcm"

    print *, "Enter namelist 'iniprint_nml'."
    read(unit=*, nml=iniprint_nml)
    write(unit_nml, nml=iniprint_nml)

    print *, "Enter namelist 'logic_nml'."
    read(unit=*, nml=logic_nml)
    write(unit_nml, nml=logic_nml)

    print *, "Enter namelist 'conf_gcm_nml'."
    read(unit=*, nml=conf_gcm_nml)
    write(unit_nml, nml=conf_gcm_nml)

    IF (MOD(day_step, iperiod) /= 0) call abort_gcm("conf_gcm", &
         'Il faut choisir un nombre de pas par jour multiple de "iperiod".')

    IF (MOD(day_step, iphysiq)/= 0) call abort_gcm("conf_gcm", &
         'Il faut choisir un nombre de pas par jour multiple de "iphysiq".')

    call assert(mod(iphysiq, iperiod) == 0, &
         "conf_gcm -- iphysiq must be multiple of iperiod")

    lmt_pas = day_step / iphysiq
    print *, 'Number of time steps of "physics" per day: ', lmt_pas

  END SUBROUTINE conf_gcm

end module conf_gcm_m
