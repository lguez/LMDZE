module clesphys2

  ! v 1.3 2005/06/06 13:16:33 fairhead

  implicit none

  LOGICAL:: cycle_diurne= .TRUE.
  ! Cette option permet d'éteindre le cycle diurne.
  ! Peut être utile pour accélérer le code.

  LOGICAL:: soil_model= .TRUE.
  ! Help = Choix du modele de sol (Thermique ?)
  ! Option qui pourait un string afin de pouvoir
  ! plus de choix ! Ou meme une liste d'options

  LOGICAL:: new_oliq= .TRUE.
  ! Help = Permet de mettre en route la
  ! nouvelle parametrisation de l'eau liquide

  LOGICAL:: ok_orodr= .TRUE.
  ! Desc = Orodr  ou  non   pour l orographie

  LOGICAL:: ok_orolf = .TRUE.
  ! Desc = Orolf  ou  non   pour l orographie

  LOGICAL:: ok_limitvrai= .FALSE.
  ! Help = On peut forcer le modele a lire le
  ! fichier SST de la bonne annee. C'est une tres bonne
  ! idee, pourquoi ne pas mettre toujours a y ?

  INTEGER:: nbapp_rad= 12
  ! (nombre d'appels des routines de rayonnements par jour)

  INTEGER:: iflag_con = 2
  ! Help = Flag  pour la convection les options suivantes existent :
  ! 1 pour LMD,
  ! 2 pour Tiedtke,
  ! 3 pour CCM(NCAR)  
  ! Flag  pour la convection (1 pour LMD, 2 pour Tiedtke, 3 KE, 4 KE vect)

end module clesphys2
