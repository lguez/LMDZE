module conf_phys_m

  ! This module is clean: no C preprocessor directive, no include line.

  implicit none

  integer iflag_pbl

contains

  subroutine conf_phys(ocean, ok_veget, ok_journe, ok_mensuel, ok_instan, &
       fact_cldcon, facttemps,ok_newmicro,iflag_cldcon, &
       ratqsbas,ratqshaut,if_ebil, &
       ok_ade, ok_aie, &
       bl95_b0, bl95_b1,&
       iflag_thermals,nsplit_thermals)

    ! From phylmd/conf_phys.F90,v 1.7 2005/07/05 07:21:23

    use getincom, only: getin
    use clesphys, only: solaire, co2_ppm, rco2, ch4_ppb, rch4, n2o_ppb, rn2o, &
         cfc11_ppt, rcfc11, cfc12_ppt, rcfc12, top_height, overlap, cdmmax, &
         cdhmax, ksta, ksta_ter, ok_kzmin, lev_histhf, lev_histday, &
         lev_histmth, type_run, ok_isccp, ok_regdyn, lonmin_ins, lonmax_ins, &
         latmin_ins, latmax_ins, ecrit_ins, ecrit_hf, ecrit_hf2mth, &
         ecrit_day, ecrit_mth, ecrit_tra, ecrit_reg, bug_ozone
    use YOMCST
    use conema3_m
    use comfisrtilp
    use nuagecom

    !IM : on inclut/initialise les taux de CH4, N2O, CFC11 et CFC12

    ! Configuration de la "physique" de LMDZ a l'aide de la fonction
    ! GETIN de IOIPSL

    ! ocean:      type d'ocean (force, slab, couple)
    ! ok_veget:   type de modele de vegetation
    ! ok_journe:  sorties journalieres
    ! ok_mensuel: sorties mensuelles
    ! ok_instan:  sorties instantanees
    ! ok_ade, ok_aie: apply or not aerosol direct and indirect effects
    ! bl95_b*: parameters in the formula to link CDNC to aerosol mass conc 

    ! Sortie:
    character(len=*), intent(out):: ocean
    logical              :: ok_veget, ok_newmicro
    logical              :: ok_journe, ok_mensuel, ok_instan        
    LOGICAL              :: ok_ade, ok_aie
    REAL                 :: bl95_b0, bl95_b1
    real, intent(out):: fact_cldcon
    real, intent(out):: facttemps
    real ratqsbas,ratqshaut
    integer              :: iflag_cldcon, if_ebil

    ! Local
    real                 :: zzz

    integer :: iflag_thermals,nsplit_thermals

    !-----------------------------------------------------------

    print *, "Call sequence information: conf_phys"

    !Config Key  = OCEAN 
    !Config Desc = Type d'ocean
    !Config Def  = force
    !Config Help = Type d'ocean utilise: force ou slab mais pas couple
    !
    ocean = 'force '
    call getin('OCEAN', ocean)
    !
    !Config Key  = VEGET 
    !Config Desc = Type de modele de vegetation
    !Config Def  = .false.
    !Config Help = Type de modele de vegetation utilise
    !
    ok_veget = .false.
    call getin('VEGET', ok_veget)
    !
    !Config Key  = OK_journe
    !Config Desc = Pour des sorties journalieres 
    !Config Def  = .false.
    !Config Help = Pour creer le fichier histday contenant les sorties
    !              journalieres 
    !
    ok_journe = .false.
    call getin('OK_journe', ok_journe)
    !
    !Config Key  = OK_mensuel
    !Config Desc = Pour des sorties mensuelles 
    !Config Def  = .true.
    !Config Help = Pour creer le fichier histmth contenant les sorties
    !              mensuelles 
    !
    ok_mensuel = .true.
    call getin('OK_mensuel', ok_mensuel)
    !
    !Config Key  = OK_instan
    !Config Desc = Pour des sorties instantanees 
    !Config Def  = .false.
    !Config Help = Pour creer le fichier histins contenant les sorties
    !              instantanees 
    !
    ok_instan = .false.
    call getin('OK_instan', ok_instan)
    !
    !Config Key  = ok_ade
    !Config Desc = Aerosol direct effect or not?
    !Config Def  = .false.
    !Config Help = Used in radlwsw.F
    !
    ok_ade = .false.
    call getin('ok_ade', ok_ade)

    !
    !Config Key  = ok_aie
    !Config Desc = Aerosol indirect effect or not?
    !Config Def  = .false.
    !Config Help = Used in nuage.F and radlwsw.F
    !
    ok_aie = .false.
    call getin('ok_aie', ok_aie)

    !
    !Config Key  = bl95_b0
    !Config Desc = Parameter in CDNC-maer link (Boucher&Lohmann 1995)
    !Config Def  = .false.
    !Config Help = Used in nuage.F
    !
    bl95_b0 = 2.
    call getin('bl95_b0', bl95_b0)

    !Config Key  = bl95_b1
    !Config Desc = Parameter in CDNC-maer link (Boucher&Lohmann 1995)
    !Config Def  = .false.
    !Config Help = Used in nuage.F
    !
    bl95_b1 = 0.2
    call getin('bl95_b1', bl95_b1)

    !
    !
    !Config Key  = if_ebil
    !Config Desc = Niveau de sortie pour les diags bilan d'energie 
    !Config Def  = 0
    !Config Help = 
    !               
    !
    if_ebil = 0
    call getin('if_ebil', if_ebil)
    !!
    !! Constante solaire & Parametres orbitaux & taux gaz effet de serre BEG
    !!
    !Config Key  = R_ecc
    !Config Desc = Excentricite
    !Config Def  = 0.016715
    !Config Help = 
    !               
    !valeur AMIP II
    R_ecc = 0.016715
    call getin('R_ecc', R_ecc)
    !!
    !Config Key  = R_peri
    !Config Desc = Equinoxe
    !Config Def  = 
    !Config Help = 
    !               
    !
    !valeur AMIP II
    R_peri = 102.7
    call getin('R_peri', R_peri)
    !!
    !Config Key  = R_incl
    !Config Desc = Inclinaison
    !Config Def  = 
    !Config Help = 
    !               
    !
    !valeur AMIP II
    R_incl = 23.441
    call getin('R_incl', R_incl)
    !!
    !Config Key  = solaire
    !Config Desc = Constante solaire en W/m2
    !Config Def  = 1365.
    !Config Help = 
    !               
    !
    !valeur AMIP II
    solaire = 1365.
    call getin('solaire', solaire)
    !!
    !Config Key  = co2_ppm
    !Config Desc = concentration du gaz carbonique en ppmv
    !Config Def  = 348.
    !Config Help = 
    !               
    !
    !valeur AMIP II
    co2_ppm = 348.
    call getin('co2_ppm', co2_ppm)
    !!
    !Config Key  = RCO2
    !Config Desc = Concentration du CO2
    !Config Def  = co2_ppm * 1.0e-06  * 44.011/28.97
    !Config Def  = 348. * 1.0e-06  * 44.011/28.97
    !Config Help = 
    !               
    ! RCO2 = 5.286789092164308E-04
    !ancienne valeur
    RCO2 = co2_ppm * 1.0e-06  * 44.011/28.97 ! pour co2_ppm=348.

    !!  call getin('RCO2', RCO2)
    !!
    !Config Key  = RCH4
    !Config Desc = Concentration du CH4
    !Config Def  = 1.65E-06* 16.043/28.97
    !Config Help = 
    !               
    !
    !valeur AMIP II
    !OK  RCH4 = 1.65E-06* 16.043/28.97
    ! RCH4 = 9.137366240938903E-07
    !
    !ancienne valeur
    ! RCH4 = 1.72E-06* 16.043/28.97
    !OK call getin('RCH4', RCH4)
    zzz = 1650.
    call getin('CH4_ppb', zzz)
    CH4_ppb = zzz
    RCH4 = CH4_ppb * 1.0E-09 * 16.043/28.97
    !!
    !Config Key  = RN2O
    !Config Desc = Concentration du N2O
    !Config Def  = 306.E-09* 44.013/28.97
    !Config Help = 
    !               
    !
    !valeur AMIP II
    !OK  RN2O = 306.E-09* 44.013/28.97
    ! RN2O = 4.648939592682085E-07
    !
    !ancienne valeur
    ! RN2O = 310.E-09* 44.013/28.97
    !OK  call getin('RN2O', RN2O)
    zzz=306.
    call getin('N2O_ppb', zzz)
    N2O_ppb = zzz
    RN2O = N2O_ppb * 1.0E-09 * 44.013/28.97
    !!
    !Config Key  = RCFC11
    !Config Desc = Concentration du CFC11
    !Config Def  = 280.E-12* 137.3686/28.97
    !Config Help = 
    !               
    !
    !OK RCFC11 = 280.E-12* 137.3686/28.97
    zzz = 280.
    call getin('CFC11_ppt',zzz)
    CFC11_ppt = zzz
    RCFC11=CFC11_ppt* 1.0E-12 * 137.3686/28.97
    ! RCFC11 = 1.327690990680013E-09
    !OK call getin('RCFC11', RCFC11)
    !!
    !Config Key  = RCFC12
    !Config Desc = Concentration du CFC12
    !Config Def  = 484.E-12* 120.9140/28.97
    !Config Help = 
    !               
    !
    !OK RCFC12 = 484.E-12* 120.9140/28.97
    zzz = 484.
    call getin('CFC12_ppt',zzz)
    CFC12_ppt = zzz
    RCFC12 = CFC12_ppt * 1.0E-12 * 120.9140/28.97
    ! RCFC12 = 2.020102726958923E-09
    !OK call getin('RCFC12', RCFC12)
    !!
    !! Constante solaire & Parametres orbitaux & taux gaz effet de serre END
    !!
    !! KE
    !
    !Config Key  = epmax
    !Config Desc = Efficacite precip
    !Config Def  = 0.993
    !Config Help = 
    !
    epmax = .993
    call getin('epmax', epmax)
    !
    !Config Key  = ok_adj_ema
    !Config Desc =  
    !Config Def  = false
    !Config Help = 
    !
    ok_adj_ema = .false.
    call getin('ok_adj_ema',ok_adj_ema)
    !
    !Config Key  = iflag_clw
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    iflag_clw = 0
    call getin('iflag_clw',iflag_clw)
    !
    !Config Key  = cld_lc_lsc 
    !Config Desc =  
    !Config Def  = 2.6e-4
    !Config Help = 
    !
    cld_lc_lsc = 2.6e-4
    call getin('cld_lc_lsc',cld_lc_lsc)
    !
    !Config Key  = cld_lc_con
    !Config Desc =  
    !Config Def  = 2.6e-4
    !Config Help = 
    !
    cld_lc_con = 2.6e-4
    call getin('cld_lc_con',cld_lc_con)
    !
    !Config Key  = cld_tau_lsc
    !Config Desc =  
    !Config Def  = 3600.
    !Config Help = 
    !
    cld_tau_lsc = 3600.
    call getin('cld_tau_lsc',cld_tau_lsc)
    !
    !Config Key  = cld_tau_con
    !Config Desc =  
    !Config Def  = 3600.
    !Config Help = 
    !
    cld_tau_con = 3600.
    call getin('cld_tau_con',cld_tau_con)
    !
    !Config Key  = ffallv_lsc
    !Config Desc =  
    !Config Def  = 1.
    !Config Help = 
    !
    ffallv_lsc = 1.
    call getin('ffallv_lsc',ffallv_lsc)
    !
    !Config Key  = ffallv_con
    !Config Desc =  
    !Config Def  = 1.
    !Config Help = 
    !
    ffallv_con = 1.
    call getin('ffallv_con',ffallv_con)
    !
    !Config Key  = coef_eva
    !Config Desc =  
    !Config Def  = 2.e-5
    !Config Help = 
    !
    coef_eva = 2.e-5
    call getin('coef_eva',coef_eva)
    !
    !Config Key  = reevap_ice
    !Config Desc =  
    !Config Def  = .false.
    !Config Help = 
    !
    reevap_ice = .false.
    call getin('reevap_ice',reevap_ice)
    !
    !Config Key  = iflag_cldcon 
    !Config Desc =  
    !Config Def  = 1
    !Config Help = 
    !
    iflag_cldcon = 1
    call getin('iflag_cldcon',iflag_cldcon)

    !
    !Config Key  = iflag_pdf 
    !Config Desc =  
    !Config Def  = 0
    !Config Help = 
    !
    iflag_pdf = 0
    call getin('iflag_pdf',iflag_pdf)
    !
    !Config Key  = fact_cldcon
    !Config Desc =  
    !Config Def  = 0.375
    !Config Help = 
    !
    fact_cldcon = 0.375
    call getin('fact_cldcon',fact_cldcon)

    !
    !Config Key  = facttemps
    !Config Desc =  
    !Config Def  = 1.e-4
    !Config Help = 
    !
    facttemps = 1.e-4
    call getin('facttemps',facttemps)

    !
    !Config Key  = ok_newmicro
    !Config Desc =  
    !Config Def  = .true.
    !Config Help = 
    !
    ok_newmicro = .true.
    call getin('ok_newmicro',ok_newmicro)
    !
    !Config Key  = ratqsbas
    !Config Desc =  
    !Config Def  = 0.01
    !Config Help = 
    !
    ratqsbas = 0.01
    call getin('ratqsbas',ratqsbas)
    !
    !Config Key  = ratqshaut
    !Config Desc =  
    !Config Def  = 0.3
    !Config Help = 
    !
    ratqshaut = 0.3
    call getin('ratqshaut',ratqshaut)

    !
    !Config Key  = rad_froid
    !Config Desc =  
    !Config Def  = 35.0
    !Config Help = 
    !
    rad_froid = 35.0
    call getin('rad_froid',rad_froid)

    !
    !Config Key  = rad_chau1
    !Config Desc =  
    !Config Def  = 13.0
    !Config Help = 
    !
    rad_chau1 = 13.0
    call getin('rad_chau1',rad_chau1)

    !
    !Config Key  = rad_chau2
    !Config Desc =  
    !Config Def  = 9.0
    !Config Help = 
    !
    rad_chau2 = 9.0
    call getin('rad_chau2',rad_chau2)

    !
    !Config Key  = top_height
    !Config Desc =
    !Config Def  = 3
    !Config Help =
    !
    top_height = 3
    call getin('top_height',top_height)

    !
    !Config Key  = overlap
    !Config Desc =
    !Config Def  = 3
    !Config Help =
    !
    overlap = 3
    call getin('overlap',overlap)


    !
    !
    !Config Key  = cdmmax
    !Config Desc =
    !Config Def  = 1.3E-3
    !Config Help =
    !
    cdmmax = 1.3E-3
    call getin('cdmmax',cdmmax)

    !
    !Config Key  = cdhmax
    !Config Desc =
    !Config Def  = 1.1E-3
    !Config Help =
    !
    cdhmax = 1.1E-3
    call getin('cdhmax',cdhmax)

    !261103
    !
    !Config Key  = ksta
    !Config Desc =
    !Config Def  = 1.0e-10
    !Config Help =
    !
    ksta = 1.0e-10
    call getin('ksta',ksta)

    !
    !Config Key  = ksta_ter
    !Config Desc =
    !Config Def  = 1.0e-10
    !Config Help =
    !
    ksta_ter = 1.0e-10
    call getin('ksta_ter',ksta_ter)

    !
    !Config Key  = ok_kzmin
    !Config Desc =
    !Config Def  = .true.
    !Config Help =
    !
    ok_kzmin = .true.
    call getin('ok_kzmin',ok_kzmin)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! PARAMETER FOR THE PLANETARY BOUNDARY LAYER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Config Key  = iflag_pbl
    !Config Desc =
    !Config Def  = 1
    !Config Help =
    !
    iflag_pbl = 1
    call getin('iflag_pbl',iflag_pbl)
    !
    !Config Key  = iflag_thermals
    !Config Desc =
    !Config Def  = 0
    !Config Help =
    !
    iflag_thermals = 0
    call getin('iflag_thermals',iflag_thermals)
    !
    !
    !Config Key  = nsplit_thermals
    !Config Desc =
    !Config Def  = 1
    !Config Help =
    !
    nsplit_thermals = 1
    call getin('nsplit_thermals',nsplit_thermals)



    !
    !Config Key  = lev_histhf
    !Config Desc =
    !Config Def  = 0
    !Config Help =
    !
    lev_histhf = 0
    call getin('lev_histhf',lev_histhf)

    !
    !Config Key  = lev_histday
    !Config Desc =
    !Config Def  = 1
    !Config Help =
    !
    lev_histday = 1
    call getin('lev_histday',lev_histday)

    !
    !Config Key  = lev_histmth
    !Config Desc =
    !Config Def  = 2
    !Config Help =
    !
    lev_histmth = 2
    call getin('lev_histmth',lev_histmth)

    !
    !Config Key  = type_run
    !Config Desc =
    !Config Def  = 'AMIP' ou 'ENSP'
    !Config Help =
    !
    type_run = 'AMIP'
    call getin('type_run',type_run)

    !
    !Config Key  = ok_isccp
    !Config Desc =
    !Config Def  = .true.
    !Config Help =
    !
    ! ok_isccp = .true.
    ok_isccp = .false.
    call getin('ok_isccp',ok_isccp)

    !
    !
    !Config Key  = ok_regdyn
    !Config Desc =
    !Config Def  = 'AMIP'
    !Config Help =
    !
    ! ok_regdyn = .true.
    ok_regdyn = .false.
    call getin('ok_regdyn',ok_regdyn)
    !
    ! coordonnees (lonmin_ins, lonmax_ins, latmin_ins, latmax_ins) pour la zone 
    ! avec sorties instantannees tous les pas de temps de la physique => "histbilKP_ins.nc"
    !
    !Config Key  = lonmin_ins
    !Config Desc = 100.  
    !Config Def  = longitude minimale sorties "bilKP_ins"
    !Config Help = 
    !
    lonmin_ins = 100.
    call getin('lonmin_ins',lonmin_ins)
    !
    !Config Key  = lonmax_ins
    !Config Desc = 130. 
    !Config Def  = longitude maximale sorties "bilKP_ins"
    !Config Help =
    !
    lonmax_ins = 130.
    call getin('lonmax_ins',lonmax_ins)
    !
    !Config Key  = latmin_ins
    !Config Desc = -20.  
    !Config Def  = latitude minimale sorties "bilKP_ins"
    !Config Help = 
    !
    latmin_ins = -20.
    call getin('latmin_ins',latmin_ins)
    !
    !Config Key  = latmax_ins
    !Config Desc = 20. 
    !Config Def  = latitude maximale sorties "bilKP_ins"
    !Config Help =
    !
    latmax_ins = 20.
    call getin('latmax_ins',latmax_ins)
    !
    !Config Key  = ecrit_ins
    !Config Desc =
    !Config Def  = NINT(86400./dtime/48.) !a chaque pas de temps physique
    !Config Help =
    !
    ! ecrit_ins = NINT(86400./dtime/48.)
    ecrit_ins = NINT(86400./48.)
    call getin('ecrit_ins',ecrit_ins)
    !
    !Config Key  = ecrit_hf
    !Config Desc =
    !Config Def  = NINT(86400./dtime *0.25) !toutes les 6h
    !Config Help =
    !
    ! ecrit_hf = NINT(86400./dtime *0.25)
    ecrit_hf = NINT(86400. *0.25)
    call getin('ecrit_hf',ecrit_hf)
    !
    !Config Key  = ecrit_hf2mth
    !Config Desc =
    !Config Def  = 4*30 !ecriture mens. a partir de val. inst. toutes les 6h
    !Config Help =
    !
    ecrit_hf2mth = 4*30
    call getin('ecrit_hf2mth',ecrit_hf2mth)
    !
    !Config Key  = ecrit_day
    ecrit_day = 86400
    call getin('ecrit_day',ecrit_day)
    !
    ecrit_mth = 86400
    call getin('ecrit_mth',ecrit_mth)
    !
    ecrit_tra = 1
    call getin('ecrit_tra',ecrit_tra)
    !
    !Config Key  = ecrit_reg
    !Config Desc =
    !Config Def  = NINT(86400./dtime *0.25)  !4 fois par jour
    !Config Help =
    !
    ! ecrit_reg = NINT(86400./dtime *0.25)  !4 fois par jour
    ecrit_reg = NINT(86400. *0.25)  !4 fois par jour
    call getin('ecrit_reg',ecrit_reg)
    !
    !
    !
    !Config Key  = bug_ozone
    !Config Desc = Pour retrouver le bug de l'ozone (IPCC), mettre a true
    !Config Def  = false
    !Config Help =
    !
    bug_ozone = .false.
    call getin('bug_ozone',bug_ozone)

    print *, ' *********'
    print *, ' Configuration des parametres de la physique: '
    print *, ' Config ocean = ', ocean
    print *, ' Config veget = ', ok_veget
    print *, ' Sortie journaliere = ', ok_journe
    print *, ' Sortie mensuelle = ', ok_mensuel
    print *, ' Sortie instantanee = ', ok_instan
    print *, ' Sortie bilan d''energie, if_ebil =', if_ebil
    print *, ' Excentricite = ',R_ecc
    print *, ' Equinoxe = ',R_peri
    print *, ' Inclinaison =',R_incl
    print *, ' Constante solaire =',solaire
    print *, ' co2_ppm =',co2_ppm
    print *, ' RCO2 = ',RCO2
    print *, ' CH4_ppb =',CH4_ppb,' RCH4 = ',RCH4
    print *, ' N2O_ppb =',N2O_ppb,' RN2O =  ',RN2O
    print *, ' CFC11_ppt=',CFC11_ppt,' RCFC11 =  ',RCFC11
    print *, ' CFC12_ppt=',CFC12_ppt,' RCFC12 =  ',RCFC12
    print *, ' epmax = ', epmax
    print *, ' ok_adj_ema = ', ok_adj_ema
    print *, ' iflag_clw = ', iflag_clw
    print *, ' cld_lc_lsc = ', cld_lc_lsc
    print *, ' cld_lc_con = ', cld_lc_con
    print *, ' cld_tau_lsc = ', cld_tau_lsc
    print *, ' cld_tau_con = ', cld_tau_con
    print *, ' ffallv_lsc = ', ffallv_lsc
    print *, ' ffallv_con = ', ffallv_con
    print *, ' coef_eva = ', coef_eva
    print *, ' reevap_ice = ', reevap_ice
    print *, ' iflag_pdf = ', iflag_pdf
    print *, ' iflag_cldcon = ', iflag_cldcon
    print *, ' fact_cldcon = ', fact_cldcon
    print *, ' facttemps = ', facttemps
    print *, ' ok_newmicro = ',ok_newmicro 
    print *, ' ratqsbas = ',ratqsbas 
    print *, ' ratqshaut = ',ratqshaut 
    print *, ' top_height = ',top_height 
    print *, ' overlap = ',overlap 
    print *, ' cdmmax = ',cdmmax 
    print *, ' cdhmax = ',cdhmax 
    print *, ' ksta = ',ksta 
    print *, ' ksta_ter = ',ksta_ter 
    print *, ' ok_kzmin = ',ok_kzmin 
    print *, ' ok_ade = ',ok_ade
    print *, ' ok_aie = ',ok_aie
    print *, ' bl95_b0 = ',bl95_b0
    print *, ' bl95_b1 = ',bl95_b1
    print *, ' lev_histhf = ',lev_histhf 
    print *, ' lev_histday = ',lev_histday 
    print *, ' lev_histmth = ',lev_histmth 
    print *, ' iflag_pbl = ', iflag_pbl
    print *, ' iflag_thermals = ', iflag_thermals
    print *, ' type_run = ',type_run 
    print *, ' ok_isccp = ',ok_isccp 
    print *, ' ok_regdyn = ',ok_regdyn
    print *, ' lonmin lonmax latmin latmax bilKP_ins =',&
         lonmin_ins, lonmax_ins, latmin_ins, latmax_ins
    print *,  'ecrit_ ins, hf, hf2mth, day, mth, reg, tra', ecrit_ins, &
         ecrit_hf, ecrit_hf2mth, ecrit_day, ecrit_mth, ecrit_reg, ecrit_tra
    print *, ' bug_ozone = ', bug_ozone

  end subroutine conf_phys

end module conf_phys_m
