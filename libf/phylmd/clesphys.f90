module clesphys

  ! v 1.3 2005/06/06 13:16:33 fairhead

  implicit none

  LOGICAL bug_ozone
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
