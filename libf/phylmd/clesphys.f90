module clesphys

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

  LOGICAL bug_ozone

  INTEGER:: nbapp_rad= 12
  ! (nombre d'appels des routines de rayonnements par jour)

  INTEGER:: iflag_con = 2
  ! Help = Flag  pour la convection les options suivantes existent :
  ! 1 pour LMD,
  ! 2 pour Tiedtke,
  ! 3 pour CCM(NCAR)  
  ! Flag  pour la convection (1 pour LMD, 2 pour Tiedtke, 3 KE, 4 KE vect)

  REAL co2_ppm, solaire
  double precision RCO2, RCH4, RN2O, RCFC11, RCFC12  
  double precision CH4_ppb, N2O_ppb, CFC11_ppt, CFC12_ppt
  !IM simulateur ISCCP 
  INTEGER top_height, overlap
  !IM seuils cdrm, cdrh
  REAL cdmmax, cdhmax
  !IM param. stabilite s/ terres et en dehors
  REAL ksta, ksta_ter
  !IM ok_kzmin : clef calcul Kzmin dans la CL de surface cf FH
  LOGICAL ok_kzmin
  INTEGER lev_histhf ! niveau sorties 6h
  integer lev_histday ! niveau sorties journalieres
  integer lev_histmth ! niveau sorties mensuelles
  CHARACTER(len=4) type_run
  LOGICAL ok_isccp, ok_regdyn
  REAL lonmin_ins, lonmax_ins, latmin_ins, latmax_ins
  INTEGER ecrit_ins, ecrit_hf, ecrit_hf2mth, ecrit_day
  INTEGER ecrit_mth, ecrit_tra, ecrit_reg

  save

end module clesphys
