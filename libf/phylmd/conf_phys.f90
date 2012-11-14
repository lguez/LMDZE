module conf_phys_m

  implicit none

  integer iflag_pbl ! parameter for the planetary boundary layer

contains

  subroutine conf_phys

    ! From phylmd/conf_phys.F90, version 1.7 2005/07/05 07:21:23

    ! Configuration de la "physique" de LMDZ.

    USE clesphys, ONLY: bug_ozone, cdhmax, cdmmax, cfc11_ppt, cfc12_ppt, &
         ch4_ppb, co2_ppm, ecrit_day, ecrit_hf, ecrit_hf2mth, ecrit_ins, &
         ecrit_mth, ecrit_reg, ecrit_tra, ksta, ksta_ter, latmax_ins, &
         latmin_ins, lev_histday, lev_histhf, lev_histmth, lonmax_ins, &
         lonmin_ins, n2o_ppb, ok_isccp, ok_kzmin, ok_regdyn, overlap, rcfc11, &
         rcfc12, rch4, rco2, rn2o, solaire, top_height, type_run
    USE comfisrtilp, ONLY: cld_lc_con, cld_lc_lsc, cld_tau_con, &
         cld_tau_lsc, coef_eva, ffallv_con, ffallv_lsc, iflag_pdf, reevap_ice
    USE conema3_m, ONLY: epmax, iflag_clw, ok_adj_ema
    USE nuagecom, ONLY: rad_chau1, rad_chau2, rad_froid
    use unit_nml_m, only: unit_nml
    USE yomcst, ONLY: r_ecc, r_incl, r_peri

    namelist /conf_phys_nml/ R_ecc, R_peri, R_incl, solaire, co2_ppm, &
         CH4_ppb, N2O_ppb, CFC11_ppt, CFC12_ppt, epmax, ok_adj_ema, &
         iflag_clw, cld_lc_lsc, cld_lc_con, cld_tau_lsc, cld_tau_con, &
         ffallv_lsc, ffallv_con, coef_eva, reevap_ice, iflag_pdf, &
         rad_froid, rad_chau1, rad_chau2, top_height, overlap, cdmmax, &
         cdhmax, ksta, ksta_ter, ok_kzmin, iflag_pbl, lev_histhf, &
         lev_histday, lev_histmth, type_run, ok_isccp, ok_regdyn, &
         lonmin_ins, lonmax_ins, latmin_ins, latmax_ins, ecrit_ins, &
         ecrit_hf, ecrit_hf2mth, ecrit_day, ecrit_mth, ecrit_tra, &
         ecrit_reg, bug_ozone

    !-----------------------------------------------------------

    print *, "Call sequence information: conf_phys"

    R_ecc = 0.016715 ! AMIP II
    R_peri = 102.7 ! AMIP II
    R_incl = 23.441 ! AMIP II
    solaire = 1365. ! AMIP II
    co2_ppm = 348. ! AMIP II
    CH4_ppb = 1650.
    N2O_ppb = 306.
    CFC11_ppt = 280.
    CFC12_ppt = 484.
    epmax = .993
    ok_adj_ema = .false.
    iflag_clw = 0
    cld_lc_lsc = 2.6e-4
    cld_lc_con = 2.6e-4
    cld_tau_lsc = 3600.
    cld_tau_con = 3600.
    ffallv_lsc = 1.
    ffallv_con = 1.
    coef_eva = 2.e-5
    reevap_ice = .false.
    iflag_pdf = 0
    rad_froid = 35.0
    rad_chau1 = 13.0
    rad_chau2 = 9.0
    top_height = 3
    overlap = 3
    cdmmax = 1.3E-3
    cdhmax = 1.1E-3
    ksta = 1.0e-10
    ksta_ter = 1.0e-10
    ok_kzmin = .true.
    iflag_pbl = 1
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
    bug_ozone = .false.

    print *, "Enter namelist 'conf_phys_nml'."
    read(unit=*, nml=conf_phys_nml)
    write(unit_nml, nml=conf_phys_nml)

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

  end subroutine conf_phys

end module conf_phys_m
