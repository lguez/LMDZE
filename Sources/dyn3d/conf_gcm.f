module conf_gcm_m

  IMPLICIT NONE

  INTEGER:: nday = 10 ! nombre de jours d'int\'egration
  integer:: day_step = 240 ! nombre de pas de temps de la dynamique par jour
  integer:: iperiod = 5 ! periode pour le pas Matsuno (en pas de temps)

  integer:: iapp_tracvl = 5
  ! Should normally be equal to "iperiod"
  ! frequence du groupement des flux (en pas de temps) 

  integer:: iconser = 240
  ! number of time steps between output of control variables

  integer:: iecri = 1 ! time interval between outputs to "dyn_hist.nc" (in days)

  integer:: iphysiq = 5
  ! number of time steps of dynamics between calls to physics

  logical:: raz_date = .false.
  ! prise en compte de la date initiale de la namelist et remise \`a
  ! z\'ero des compteurs de pas de temps (sinon on garde la date du
  ! fichier restart)

  integer:: periodav = 1 
  ! time interval between outputs in the dynamical part, in days

  logical:: offline = .FALSE.
  ! permet de mettre en route la nouvelle parametrisation de l'eau liquide

  integer:: prt_level = 0 ! niveau d'impression souhait\'e (0 = minimum)

  LOGICAL:: purmats= .FALSE.
  ! Help = Choix du schema d'integration temporel.
  ! y = pure Matsuno sinon c'est du Matsuno-leapfrog

  INTEGER:: iflag_phys = 1
  ! contr\^ole l'appel \`a la physique :
  ! 0 : pas de physique
  ! 1 : physique normale (appel \`a phylmd, phymars...) (default)
  ! 2 : rappel Newtonien pour la temp\'erature + friction au sol

contains

  SUBROUTINE conf_gcm

    ! Auteurs : L. Fairhead, P. Le Van
    ! Version du 29/04/97

    use abort_gcm_m, only: abort_gcm
    use comdissnew, only: read_comdissnew
    use unit_nml_m, only: unit_nml

    namelist /conf_gcm_nml/ raz_date, nday, day_step, iperiod, iapp_tracvl, &
         iconser, iecri, periodav, iphysiq, offline

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

    call read_comdissnew

    print *, "Enter namelist 'conf_gcm_nml'."
    read(unit=*, nml=conf_gcm_nml)
    write(unit_nml, nml=conf_gcm_nml)

    IF (MOD(day_step, iperiod) /= 0) call abort_gcm(modname = "conf_gcm", &
         message = &
         'Il faut choisir un nombre de pas par jour multiple de "iperiod".', &
         ierr = 1)

    IF (MOD(day_step, iphysiq)/= 0) call abort_gcm(modname = "conf_gcm", &
         message = &
         'Il faut choisir un nombre de pas par jour multiple de "iphysiq".', &
         ierr = 1)

  END SUBROUTINE conf_gcm

end module conf_gcm_m
