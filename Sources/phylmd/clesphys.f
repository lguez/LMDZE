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
  logical:: ok_instan = .false. ! sorties instantanees dans le fichier histins

  save

contains

  subroutine read_clesphys

    use unit_nml_m, only: unit_nml

    namelist /clesphys_nml/ solaire, co2_ppm, CH4_ppb, N2O_ppb, CFC11_ppt, &
         CFC12_ppt, top_height, overlap, cdmmax, cdhmax, ksta, &
         ksta_ter, ok_kzmin, lev_histhf, lev_histday, lev_histmth, &
         type_run, ok_isccp, ok_regdyn, lonmin_ins, lonmax_ins, latmin_ins, &
         latmax_ins, ecrit_ins, ecrit_hf, ecrit_hf2mth, ecrit_day, ecrit_mth, &
         ecrit_tra, ecrit_reg, ok_instan

    !---------------------------------------------------------

    solaire = 1365. ! AMIP II
    co2_ppm = 348. ! AMIP II
    CH4_ppb = 1650.
    N2O_ppb = 306.
    CFC11_ppt = 280.
    CFC12_ppt = 484.
    top_height = 3
    overlap = 3
    cdmmax = 1.3E-3
    cdhmax = 1.1E-3
    ksta = 1.0e-10
    ksta_ter = 1.0e-10
    ok_kzmin = .true.
    lev_histhf = 0
    lev_histday = 1
    lev_histmth = 2
    type_run = 'AMIP'
    ok_isccp = .false.
    ok_regdyn = .false.
    lonmin_ins = 100.
    lonmax_ins = 130.
    latmin_ins = -20.
    latmax_ins = 20.
    ecrit_ins = NINT(86400./48.)
    ecrit_hf = NINT(86400. *0.25)
    ecrit_hf2mth = 4*30 ! ecriture mens. a partir de val. inst. toutes les 6h
    ecrit_day = 86400
    ecrit_mth = 86400
    ecrit_tra = 1
    ecrit_reg = NINT(86400. *0.25) ! 4 fois par jour

    print *, "Enter namelist 'clesphys_nml'."
    read(unit=*, nml=clesphys_nml)
    write(unit_nml, nml=clesphys_nml)

    RCO2 = co2_ppm * 1.0e-06 * 44.011/28.97
    RCH4 = CH4_ppb * 1.0E-09 * 16.043/28.97
    RN2O = N2O_ppb * 1.0E-09 * 44.013/28.97
    RCFC11=CFC11_ppt* 1.0E-12 * 137.3686/28.97
    RCFC12 = CFC12_ppt * 1.0E-12 * 120.9140/28.97

    print *, ' RCO2 = ', RCO2
    print *, ' RCH4 = ', RCH4
    print *, ' RN2O = ', RN2O
    print *, ' RCFC11 = ', RCFC11
    print *, ' RCFC12 = ', RCFC12

  end subroutine read_clesphys

end module clesphys
