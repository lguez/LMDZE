module clesphys

  ! From version 1.3 2005/06/06 13:16:33

  implicit none

  REAL co2_ppm ! concentration du gaz carbonique en ppmv
  real solaire ! Constante solaire en W/m2
  double precision RCO2 ! Concentration du CO2
  double precision RCH4 ! Concentration du CH4
  double precision RN2O ! Concentration du N2O
  double precision RCFC11 ! Concentration du CFC11
  double precision RCFC12 ! Concentration du CFC12
  double precision CH4_ppb, N2O_ppb, CFC11_ppt, CFC12_ppt

  ! Simulateur ISCCP:
  INTEGER top_height
  INTEGER overlap ! 1, 2 or 3

  REAL cdmmax, cdhmax ! seuils cdrm, cdrh
  REAL ksta, ksta_ter ! paramètres stabilité sur terres et en dehors
  LOGICAL ok_kzmin ! calcul Kzmin dans la couche limite de surface

  INTEGER lev_histhf ! niveau sorties 6h
  ! 4: histhf3d.nc champs 3d niveaux modele

  integer lev_histday ! niveau sorties journalieres
  ! 3: champs 3D => F. Lott
  ! 4: + champs sous-surfaces

  integer lev_histmth ! niveau sorties mensuelles
  ! 3: albedo, rugosite sous-surfaces
  ! 4: champs tendances 3d

  CHARACTER(len=4) type_run ! 'AMIP' ou 'ENSP'
  LOGICAL ok_isccp, ok_regdyn

  REAL lonmin_ins, lonmax_ins, latmin_ins, latmax_ins
  ! longitude et latitude minimales et maximales pour la zone avec
  ! sorties instantanées tous les pas de temps de la physique,
  ! fichier "histbilKP_ins.nc"

  INTEGER ecrit_ins, ecrit_hf, ecrit_hf2mth, ecrit_day
  INTEGER ecrit_mth, ecrit_tra, ecrit_reg

  save

end module clesphys
