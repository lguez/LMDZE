module conf_gcm_m

  IMPLICIT NONE

  INTEGER:: nday = 10 ! nombre de jours d'intégration
  integer:: day_step = 240 ! nombre de pas par jour, multiple de iperiod

  integer:: iperiod = 5
  ! periode pour le pas Matsuno (en pas de temps)

  integer:: iapp_tracvl = 5
  ! Should normally be equal to "iperiod"
  ! frequence du groupement des flux (en pas de temps) 

  integer:: iconser = 240
  ! number of time steps between output of control variables

  integer:: iecri = 1 ! time interval between outputs to "dyn_hist.nc" (in days)

  integer:: iphysiq = 5
  ! number of time steps of dynamics between calls to physics

  integer:: dayref = 1 ! jour de l'année de l'état initial
  ! (= 350 si 20 décembre par exemple)

  integer:: anneeref = 1998 ! Annee de l'etat initial (avec 4 chiffres)

  logical:: raz_date = .false.
  ! (remise a zero de la date initiale, prise en compte de la date de
  ! gcm.def avec remise a zero des compteurs de pas de temps)
  ! (pas de remise a zero: on garde la date du fichier restart)

  integer:: periodav = 1 
  ! time interval between outputs in the dynamical part, in days

  logical:: offline = .FALSE.
  ! permet de mettre en route la nouvelle parametrisation de l'eau liquide

  integer:: prt_level = 0 ! niveau d'impression souhaité (0 = minimum)

  LOGICAL:: purmats= .FALSE.
  ! Help = Choix du schema d'integration temporel.
  ! y = pure Matsuno sinon c'est du Matsuno-leapfrog

  logical:: fxyhypb = .TRUE.
  ! (fonction f(y) à dérivée tangente hyperbolique, sinon à dérivée
  ! sinusoïdale)

  logical:: ysinus = .TRUE.
  ! (Fonction f(y) avec y = Sin(latit.) si = .true. sinon y = latit.)

  logical:: ok_guide= .FALSE. ! guidage

  INTEGER:: iflag_phys = 1
  ! contrôle l'appel à la physique :
  ! 0 : pas de physique
  ! 1 : physique normale (appel à phylmd, phymars...) (default)
  ! 2 : rappel Newtonien pour la température + friction au sol

contains

  SUBROUTINE conf_gcm

    ! Auteurs : L. Fairhead, P. Le Van
    ! Version du 29/04/97

    ! On ne compare pas les valeurs des paramètres du zoom, grossismx, 
    ! grossismy, clon, clat, fxyhypb lues sur le fichier start avec
    ! celles passées par run.def, au début du gcm. 
    ! Ces paramètres définissent entre autres la grille et doivent être
    ! cohérents, sinon il y aura divergence du gcm.

    use abort_gcm_m, only: abort_gcm
    use comdissnew, only: read_comdissnew
    use serre, only: clon, clat, grossismx, grossismy, alphax, alphay, &
         dzoomx, dzoomy, taux, tauy
    use unit_nml_m, only: unit_nml

    namelist /conf_gcm_nml/dayref, anneeref, raz_date, nday, day_step, &
         iperiod, iapp_tracvl, iconser, iecri, periodav, &
         iphysiq, clon, clat, grossismx, grossismy, dzoomx, dzoomy, taux, &
         tauy, offline

    namelist /iniprint_nml/prt_level

    namelist /logic_nml/ purmats, fxyhypb, ysinus, ok_guide, iflag_phys

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

    IF (grossismx < 1.) THEN
       PRINT *, 'Error: grossismx < 1'
       STOP 1
    ELSE
       alphax = 1. - 1. / grossismx
    ENDIF
    IF (grossismy < 1.) THEN
       PRINT *, 'Error: grossismy < 1'
       STOP 1
    ELSE
       alphay = 1. - 1. / grossismy
    ENDIF
    PRINT *, 'alphax = ', alphax
    PRINT *, 'alphay = ', alphay

  END SUBROUTINE conf_gcm

end module conf_gcm_m
