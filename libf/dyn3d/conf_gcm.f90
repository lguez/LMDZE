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

  integer:: iecri= 1 ! p�riode d'�criture du fichier "dyn_hist.nc" (en jours)

  integer:: idissip= 10 ! periode de la dissipation (en pas)

  integer:: iphysiq= 5
  ! Help = Periode de la physique en pas de temps de la dynamique.

  integer:: dayref = 1 ! jour de l'ann�e de l'�tat initial
  ! (= 350 si 20 d�cembre par exemple)

  integer:: anneeref = 1998 ! Annee de l'etat initial (avec 4 chiffres)

  integer:: raz_date = 0 ! Remise a zero de la date initiale
  ! 0 pas de remise a zero, on garde la date du fichier restart
  ! 1 prise en compte de la date de gcm.def avec remise a zero
  ! des compteurs de pas de temps

  REAL:: periodav= 1.
  ! periode de stockage fichier histmoy (en jour) 

  logical:: offline = .FALSE.
  ! Nouvelle eau liquide
  ! Permet de mettre en route la nouvelle parametrisation de l'eau liquide

contains

  SUBROUTINE conf_gcm

    ! Auteurs : L. Fairhead, P. Le Van
    ! Version du 29/04/97

    ! Nouveaux param�tres nitergdiv, nitergrot, niterh, tetagdiv, tetagrot, 
    ! tetatemp ajout�s pour la dissipation.

    ! On ne compare pas les valeurs des param�tres du zoom, grossismx, 
    ! grossismy, clon, clat, fxyhypb lues sur le fichier start avec
    ! celles pass�es par run.def, au d�but du gcm. 
    ! Ces param�tres d�finissent entre autres la grille et doivent �tre
    ! coh�rents, sinon il y aura divergence du gcm.

    use comdissnew, only: read_comdissnew
    use logic, only: read_logic
    use serre, only: clon, clat, grossismx, grossismy, alphax, alphay, &
         dzoomx, dzoomy, taux, tauy
    use iniprint, only: read_iniprint

    namelist /conf_gcm_nml/dayref, anneeref, raz_date, nday, day_step, &
         iperiod, iapp_tracvl, iconser, iecri, periodav, idissip, &
         iphysiq, clon, clat, grossismx, grossismy, dzoomx, dzoomy, taux, &
         tauy, offline

    !------------------------------------

    print *, "Call sequence information: conf_gcm"

    call read_iniprint
    call read_logic
    call read_comdissnew

    print *, "Enter namelist 'conf_gcm_nml'."
    read(unit=*, nml=conf_gcm_nml)
    write(unit=*, nml=conf_gcm_nml)

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
