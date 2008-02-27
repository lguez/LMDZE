module conf_gcm_m

  ! This module is clean: no C preprocessor directive, no include line

  IMPLICIT NONE

  INTEGER:: nday= 10
  ! Nombre de jours d'integration
  ! On pourait aussi permettre des mois ou des annees !

  integer:: day_step= 240 ! nombre de pas par jour, multiple de iperiod

  integer:: iperiod= 5
  ! periode pour le pas Matsuno (en pas de temps)

  integer:: iapp_tracvl= 5
  ! Should normally be equal to "iperiod"
  ! frequence du groupement des flux (en pas de temps) 

  integer:: iconser= 240
  ! periode de sortie des variables de controle
  ! (En pas de temps)

  integer:: iecri= 1 ! période d'écriture du fichier "dyn_hist.nc" (en jours)

  integer:: idissip= 10 ! periode de la dissipation (en pas)

  integer:: iphysiq= 5
  ! Help = Periode de la physique en pas de temps de la dynamique.

  integer:: dayref = 1 ! jour de l'année de l'état initial
  ! (= 350 si 20 décembre par exemple)

  integer:: anneeref = 1998 ! Annee de l'etat initial (avec 4 chiffres)

  integer:: raz_date = 0 ! Remise a zero de la date initiale
  ! 0 pas de remise a zero, on garde la date du fichier restart
  ! 1 prise en compte de la date de gcm.def avec remise a zero
  ! des compteurs de pas de temps

  REAL:: periodav= 1.
  ! periode de stockage fichier histmoy (en jour) 

  logical, save:: offline

contains

  SUBROUTINE conf_gcm(clesphy0)

    ! Auteurs : L. Fairhead, P. Le Van
    ! Version du 29/04/97

    ! Nouveaux paramètres nitergdiv, nitergrot, niterh, tetagdiv, tetagrot, 
    ! tetatemp ajoutés pour la dissipation.

    ! Autre paramètre ajouté en fin de liste de tapedef : fxyhypb
    ! Si fxyhypb = .TRUE., choix de la fonction à dérivée tangente
    ! hyperbolique
    ! Sinon, choix de fxynew, à dérivée sinusoïdale

    ! On ne compare pas les valeurs des paramètres du zoom, grossismx, 
    ! grossismy, clon, clat, fxyhypb lues sur le fichier start avec
    ! celles passées par run.def, au début du gcm. 
    ! Ces paramètres définissent entre autres la grille et doivent être
    ! cohérents, sinon il y aura divergence du gcm.

    use IOIPSL, only: getin
    use dimens_m
    use paramet_m
    use comdissnew, only: read_comdissnew
    use logic, only: iflag_phys, fxyhypb, ysinus, purmats, ok_guide
    use serre, only: clon, clat, grossismx, grossismy, alphax, alphay, &
         dzoomx, dzoomy, taux, tauy
    use clesphys, only: ok_limitvrai, ok_orolf, iflag_con, nbapp_rad, &
         ok_orodr, cycle_diurne, new_oliq, soil_model
    use iniprint, only: lunout, read_iniprint

    INTEGER, PARAMETER:: longcles = 20

    REAL, intent(out), optional:: clesphy0(longcles)

    namelist /conf_gcm_nml/dayref, anneeref, raz_date, nday, day_step, &
         iperiod, iapp_tracvl, iconser, iecri, periodav, idissip, purmats, &
         ok_guide, iflag_phys, iphysiq, cycle_diurne, soil_model, new_oliq, &
         ok_orodr, ok_orolf, ok_limitvrai, nbapp_rad, iflag_con, clon, clat, &
         grossismx, grossismy, dzoomx, dzoomy, taux, tauy

    !------------------------------------

    print *, "Call sequence information: conf_gcm"
    call read_iniprint

    print *, "Enter namelist 'conf_gcm_nml'."
    read(unit=*, nml=conf_gcm_nml)
    write(unit=*, nml=conf_gcm_nml)

    call read_comdissnew

    IF (lunout /= 5 .and. lunout /= 6) THEN
       OPEN(lunout, FILE='lmdz.out')
    ENDIF

    if (present(clesphy0)) then
       clesphy0(:) = 0.
       clesphy0(1) = REAL(iflag_con)
       clesphy0(2) = REAL(nbapp_rad)

       IF (cycle_diurne) clesphy0(3) = 1.
       IF (soil_model) clesphy0(4) = 1.
       IF (new_oliq) clesphy0(5) = 1.
       IF (ok_orodr) clesphy0(6) = 1.
       IF (ok_orolf) clesphy0(7) = 1.
       IF (ok_limitvrai) clesphy0(8) = 1.
    end if

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

    ! Key = fxyhypb
    ! Desc = Fonction hyperbolique
    ! Def = y
    ! Help = Fonction f(y) hyperbolique si = .true. 
    ! sinon sinusoidale
    fxyhypb = .TRUE.
    CALL getin('fxyhypb', fxyhypb)

    ! Key = ysinus
    ! IF = !fxyhypb
    ! Desc = Fonction en Sinus
    ! Def = y
    ! Help = Fonction f(y) avec y = Sin(latit.) si = .true. 
    ! sinon y = latit.
    ysinus = .TRUE.
    CALL getin('ysinus', ysinus)

    ! Key = offline
    ! Desc = Nouvelle eau liquide
    ! Def = n
    ! Help = Permet de mettre en route la
    ! nouvelle parametrisation de l'eau liquide
    offline = .FALSE.
    CALL getin('offline', offline)
    write(lunout, *)' offline = ', offline

  END SUBROUTINE conf_gcm

end module conf_gcm_m
