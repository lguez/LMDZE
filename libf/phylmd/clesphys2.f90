module clesphys2

  ! From version 1.3 2005/06/06 13:16:33

  implicit none

  LOGICAL:: cycle_diurne = .TRUE.
  ! Cette option permet d'éteindre le cycle diurne.
  ! Peut être utile pour accélérer le code.

  LOGICAL:: soil_model = .TRUE.
  ! Choix du modele de sol (Thermique ?)
  ! Option qui pourait un string afin de pouvoir
  ! plus de choix ! Ou meme une liste d'options

  LOGICAL:: new_oliq = .TRUE.
  ! Permet de mettre en route la nouvelle parametrisation de l'eau liquide

  LOGICAL:: ok_orodr = .TRUE. ! pour l'orographie
  LOGICAL:: ok_orolf = .TRUE. ! pour l'orographie

  LOGICAL:: ok_limitvrai = .FALSE.
  ! On peut forcer le modele a lire le fichier SST de la bonne
  ! annee. C'est une tres bonne idee, pourquoi ne pas mettre toujours
  ! a y ?

  INTEGER:: nbapp_rad = 12
  ! (nombre d'appels des routines de rayonnements par jour)

  INTEGER:: iflag_con = 2
  ! Flag pour la convection :
  ! 1 LMD,
  ! 2 Tiedtke,
  ! 3 CCM(NCAR)  
  ! 3 KE
  ! 4 KE vect

contains

  subroutine read_clesphys2

    namelist /clesphys2_nml/cycle_diurne, soil_model, new_oliq, &
         ok_orodr, ok_orolf, ok_limitvrai, nbapp_rad, iflag_con

    !------------------------------------

    print *, "Enter namelist 'clesphys2_nml'."
    read(unit=*, nml=clesphys2_nml)
    write(unit=*, nml=clesphys2_nml)

    select case (iflag_con)
    case (1)
       PRINT *, 'Schéma convection LMD'
    case (2)
       PRINT *, 'Schéma convection Tiedtke'
    case (3)
       PRINT *, 'Schéma convection CCM'
    END select

  end subroutine read_clesphys2

end module clesphys2
