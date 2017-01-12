module clesphys

  ! From version 1.3 2005/06/06 13:16:33

  implicit none

  REAL:: co2_ppm = 348. ! AMIP II, concentration du gaz carbonique en ppmv
  real:: solaire = 1365. ! AMIP II, constante solaire en W/m2
  double precision, save:: RCO2 ! Concentration du CO2
  double precision, save:: RCH4 ! Concentration du CH4
  double precision, save:: RN2O ! Concentration du N2O
  double precision, save:: RCFC11 ! Concentration du CFC11
  double precision, save:: RCFC12 ! Concentration du CFC12
  double precision:: CH4_ppb = 1650., N2O_ppb = 306.
  double precision:: CFC11_ppt = 280., CFC12_ppt = 484.
  REAL:: cdmmax = 1.3E-3, cdhmax = 1.1E-3 ! seuils cdrm, cdrh
  
  REAL:: ksta = 1.0e-10, ksta_ter = 1.0e-10
  ! param\`etres stabilit\'e sur terre et en dehors

  LOGICAL:: ok_kzmin  = .true.! calcul Kzmin dans la couche limite de surface
  INTEGER, save:: ecrit_ins
  logical:: ok_instan = .false. ! sorties instantanees dans le fichier histins

contains

  subroutine read_clesphys

    use unit_nml_m, only: unit_nml

    namelist /clesphys_nml/ solaire, co2_ppm, CH4_ppb, N2O_ppb, CFC11_ppt, &
         CFC12_ppt, cdmmax, cdhmax, ksta, ksta_ter, ok_kzmin, ecrit_ins, &
         ok_instan

    !---------------------------------------------------------

    ecrit_ins = NINT(86400./48.)

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
