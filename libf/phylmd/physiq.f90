module physiq_m

  ! This module is clean: no C preprocessor directive, no include line.

  IMPLICIT none

  private
  public physiq

contains

  SUBROUTINE physiq(nq, firstcal, lafin, rdayvrai, gmtime, pdtphys, paprs, &
       pplay, pphi, pphis, presnivs, u, v, t, qx, omega, d_u, d_v, &
       d_t, d_qx, d_ps, dudyn, PVteta)

    ! From phylmd/physiq.F, v 1.22 2006/02/20 09:38:28

    ! Author : Z.X. Li (LMD/CNRS), date: 1993/08/18

    ! Objet: Moniteur general de la physique du modele
    !AA      Modifications quant aux traceurs :
    !AA                  -  uniformisation des parametrisations ds phytrac
    !AA                  -  stockage des moyennes des champs necessaires
    !AA                     en mode traceur off-line 

    USE ioipsl, only: ymds2ju, histwrite, histsync
    use dimens_m, only: jjm, iim, llm
    use indicesol, only: nbsrf, is_ter, is_lic, is_sic, is_oce, &
         clnsurf, epsfra
    use dimphy, only: klon, nbtr
    use conf_gcm_m, only: raz_date, offline, iphysiq
    use dimsoil, only: nsoilmx
    use temps, only: itau_phy, day_ref, annee_ref, itaufin
    use clesphys, only: ecrit_hf, ecrit_ins, ecrit_mth, &
         cdmmax, cdhmax, &
         co2_ppm, ecrit_reg, ecrit_tra, ksta, ksta_ter, &
         ok_kzmin
    use clesphys2, only: iflag_con, ok_orolf, ok_orodr, nbapp_rad, &
         cycle_diurne, new_oliq, soil_model
    use iniprint, only: prt_level
    use abort_gcm_m, only: abort_gcm
    use YOMCST, only: rcpd, rtt, rlvtt, rg, ra, rsigma, retv, romega
    use comgeomphy
    use ctherm
    use phytrac_m, only: phytrac
    use oasis_m
    use radepsi
    use radopt
    use yoethf
    use ini_hist, only: ini_histhf, ini_histday, ini_histins
    use orbite_m, only: orbite, zenang
    use phyetat0_m, only: phyetat0, rlat, rlon
    use hgardfou_m, only: hgardfou
    use conf_phys_m, only: conf_phys
    use phyredem_m, only: phyredem
    use qcheck_m, only: qcheck

    ! Declaration des constantes et des fonctions thermodynamiques :
    use fcttre, only: thermcep, foeew, qsats, qsatl

    ! Variables argument:

    INTEGER, intent(in):: nq ! nombre de traceurs (y compris vapeur d'eau)

    REAL, intent(in):: rdayvrai
    ! (elapsed time since January 1st 0h of the starting year, in days)

    REAL, intent(in):: gmtime ! heure de la journée en fraction de jour
    REAL, intent(in):: pdtphys ! pas d'integration pour la physique (seconde)
    LOGICAL, intent(in):: firstcal ! first call to "calfis"
    logical, intent(in):: lafin ! dernier passage

    REAL, intent(in):: paprs(klon, llm+1)
    ! (pression pour chaque inter-couche, en Pa)

    REAL, intent(in):: pplay(klon, llm)
    ! (input pression pour le mileu de chaque couche (en Pa))

    REAL pphi(klon, llm)  
    ! (input geopotentiel de chaque couche (g z) (reference sol))

    REAL pphis(klon) ! input geopotentiel du sol

    REAL presnivs(llm)
    ! (input pressions approximat. des milieux couches ( en PA))

    REAL u(klon, llm)  ! input vitesse dans la direction X (de O a E) en m/s
    REAL v(klon, llm)  ! input vitesse Y (de S a N) en m/s
    REAL t(klon, llm)  ! input temperature (K)

    REAL, intent(in):: qx(klon, llm, nq)
    ! (humidite specifique (kg/kg) et fractions massiques des autres traceurs)

    REAL omega(klon, llm)  ! input vitesse verticale en Pa/s
    REAL d_u(klon, llm)  ! output tendance physique de "u" (m/s/s)
    REAL d_v(klon, llm)  ! output tendance physique de "v" (m/s/s)
    REAL d_t(klon, llm)  ! output tendance physique de "t" (K/s)
    REAL d_qx(klon, llm, nq)  ! output tendance physique de "qx" (kg/kg/s)
    REAL d_ps(klon)  ! output tendance physique de la pression au sol

    INTEGER nbteta
    PARAMETER(nbteta=3)

    REAL PVteta(klon, nbteta) 
    ! (output vorticite potentielle a des thetas constantes)

    LOGICAL ok_cvl  ! pour activer le nouveau driver pour convection KE
    PARAMETER (ok_cvl=.TRUE.)
    LOGICAL ok_gust ! pour activer l'effet des gust sur flux surface
    PARAMETER (ok_gust=.FALSE.)

    LOGICAL check ! Verifier la conservation du modele en eau
    PARAMETER (check=.FALSE.)
    LOGICAL ok_stratus ! Ajouter artificiellement les stratus
    PARAMETER (ok_stratus=.FALSE.)

    ! Parametres lies au coupleur OASIS:
    INTEGER, SAVE :: npas, nexca
    logical rnpb
    parameter(rnpb=.true.)

    character(len=6), save:: ocean
    ! (type de modèle océan à utiliser: "force" ou "slab" mais pas "couple")

    logical ok_ocean
    SAVE ok_ocean

    !IM "slab" ocean
    REAL tslab(klon)    !Temperature du slab-ocean
    SAVE tslab
    REAL seaice(klon)   !glace de mer (kg/m2) 
    SAVE seaice
    REAL fluxo(klon)    !flux turbulents ocean-glace de mer 
    REAL fluxg(klon)    !flux turbulents ocean-atmosphere

    ! Modele thermique du sol, a activer pour le cycle diurne:
    logical, save:: ok_veget
    LOGICAL, save:: ok_journe ! sortir le fichier journalier

    LOGICAL ok_mensuel ! sortir le fichier mensuel

    LOGICAL ok_instan ! sortir le fichier instantane
    save ok_instan

    LOGICAL ok_region ! sortir le fichier regional
    PARAMETER (ok_region=.FALSE.)

    !     pour phsystoke avec thermiques
    REAL fm_therm(klon, llm+1)
    REAL entr_therm(klon, llm)
    real q2(klon, llm+1, nbsrf)
    save q2

    INTEGER ivap          ! indice de traceurs pour vapeur d'eau
    PARAMETER (ivap=1)
    INTEGER iliq          ! indice de traceurs pour eau liquide
    PARAMETER (iliq=2)

    REAL t_ancien(klon, llm), q_ancien(klon, llm)
    SAVE t_ancien, q_ancien
    LOGICAL ancien_ok
    SAVE ancien_ok

    REAL d_t_dyn(klon, llm) ! tendance dynamique pour "t" (K/s)
    REAL d_q_dyn(klon, llm)  ! tendance dynamique pour "q" (kg/kg/s)

    real da(klon, llm), phi(klon, llm, llm), mp(klon, llm)

    !IM Amip2 PV a theta constante 

    CHARACTER(LEN=3) ctetaSTD(nbteta)
    DATA ctetaSTD/'350', '380', '405'/
    REAL rtetaSTD(nbteta)
    DATA rtetaSTD/350., 380., 405./

    !MI Amip2 PV a theta constante

    INTEGER klevp1
    PARAMETER(klevp1=llm+1)

    REAL swdn0(klon, klevp1), swdn(klon, klevp1)
    REAL swup0(klon, klevp1), swup(klon, klevp1)
    SAVE swdn0, swdn, swup0, swup

    REAL SWdn200clr(klon), SWdn200(klon)
    REAL SWup200clr(klon), SWup200(klon)
    SAVE SWdn200clr, SWdn200, SWup200clr, SWup200

    REAL lwdn0(klon, klevp1), lwdn(klon, klevp1)
    REAL lwup0(klon, klevp1), lwup(klon, klevp1)
    SAVE lwdn0, lwdn, lwup0, lwup 

    REAL LWdn200clr(klon), LWdn200(klon)
    REAL LWup200clr(klon), LWup200(klon)
    SAVE LWdn200clr, LWdn200, LWup200clr, LWup200

    !IM Amip2
    ! variables a une pression donnee

    integer nlevSTD
    PARAMETER(nlevSTD=17)
    real rlevSTD(nlevSTD)
    DATA rlevSTD/100000., 92500., 85000., 70000., &
         60000., 50000., 40000., 30000., 25000., 20000., &
         15000., 10000., 7000., 5000., 3000., 2000., 1000./
    CHARACTER(LEN=4) clevSTD(nlevSTD)
    DATA clevSTD/'1000', '925 ', '850 ', '700 ', '600 ', &
         '500 ', '400 ', '300 ', '250 ', '200 ', '150 ', '100 ', &
         '70  ', '50  ', '30  ', '20  ', '10  '/

    real tlevSTD(klon, nlevSTD), qlevSTD(klon, nlevSTD)
    real rhlevSTD(klon, nlevSTD), philevSTD(klon, nlevSTD)
    real ulevSTD(klon, nlevSTD), vlevSTD(klon, nlevSTD)
    real wlevSTD(klon, nlevSTD) 

    ! nout : niveau de output des variables a une pression donnee
    INTEGER nout
    PARAMETER(nout=3) !nout=1 : day; =2 : mth; =3 : NMC

    logical oknondef(klon, nlevSTD, nout)
    real tnondef(klon, nlevSTD, nout) 
    save tnondef

    ! les produits uvSTD, vqSTD, .., T2STD sont calcules
    ! a partir des valeurs instantannees toutes les 6 h
    ! qui sont moyennees sur le mois

    real uvSTD(klon, nlevSTD)
    real vqSTD(klon, nlevSTD)
    real vTSTD(klon, nlevSTD)
    real wqSTD(klon, nlevSTD)

    real vphiSTD(klon, nlevSTD)
    real wTSTD(klon, nlevSTD)
    real u2STD(klon, nlevSTD)
    real v2STD(klon, nlevSTD)
    real T2STD(klon, nlevSTD)

    ! prw: precipitable water
    real prw(klon)

    ! flwp, fiwp = Liquid Water Path & Ice Water Path (kg/m2)
    ! flwc, fiwc = Liquid Water Content & Ice Water Content (kg/kg)
    REAL flwp(klon), fiwp(klon)
    REAL flwc(klon, llm), fiwc(klon, llm)

    INTEGER l, kmax, lmax
    PARAMETER(kmax=8, lmax=8)
    INTEGER kmaxm1, lmaxm1
    PARAMETER(kmaxm1=kmax-1, lmaxm1=lmax-1)

    REAL zx_tau(kmaxm1), zx_pc(lmaxm1)
    DATA zx_tau/0.0, 0.3, 1.3, 3.6, 9.4, 23., 60./
    DATA zx_pc/50., 180., 310., 440., 560., 680., 800./

    ! cldtopres pression au sommet des nuages
    REAL cldtopres(lmaxm1)
    DATA cldtopres/50., 180., 310., 440., 560., 680., 800./

    ! taulev: numero du niveau de tau dans les sorties ISCCP
    CHARACTER(LEN=4) taulev(kmaxm1)

    DATA taulev/'tau0', 'tau1', 'tau2', 'tau3', 'tau4', 'tau5', 'tau6'/
    CHARACTER(LEN=3) pclev(lmaxm1)
    DATA pclev/'pc1', 'pc2', 'pc3', 'pc4', 'pc5', 'pc6', 'pc7'/

    CHARACTER(LEN=28) cnameisccp(lmaxm1, kmaxm1)
    DATA cnameisccp/'pc< 50hPa, tau< 0.3', 'pc= 50-180hPa, tau< 0.3', &
         'pc= 180-310hPa, tau< 0.3', 'pc= 310-440hPa, tau< 0.3', &
         'pc= 440-560hPa, tau< 0.3', 'pc= 560-680hPa, tau< 0.3', &
         'pc= 680-800hPa, tau< 0.3', 'pc< 50hPa, tau= 0.3-1.3', &
         'pc= 50-180hPa, tau= 0.3-1.3', 'pc= 180-310hPa, tau= 0.3-1.3', &
         'pc= 310-440hPa, tau= 0.3-1.3', 'pc= 440-560hPa, tau= 0.3-1.3', &
         'pc= 560-680hPa, tau= 0.3-1.3', 'pc= 680-800hPa, tau= 0.3-1.3', &
         'pc< 50hPa, tau= 1.3-3.6', 'pc= 50-180hPa, tau= 1.3-3.6', &
         'pc= 180-310hPa, tau= 1.3-3.6', 'pc= 310-440hPa, tau= 1.3-3.6', &
         'pc= 440-560hPa, tau= 1.3-3.6', 'pc= 560-680hPa, tau= 1.3-3.6', &
         'pc= 680-800hPa, tau= 1.3-3.6', 'pc< 50hPa, tau= 3.6-9.4', &
         'pc= 50-180hPa, tau= 3.6-9.4', 'pc= 180-310hPa, tau= 3.6-9.4', &
         'pc= 310-440hPa, tau= 3.6-9.4', 'pc= 440-560hPa, tau= 3.6-9.4', &
         'pc= 560-680hPa, tau= 3.6-9.4', 'pc= 680-800hPa, tau= 3.6-9.4', &
         'pc< 50hPa, tau= 9.4-23', 'pc= 50-180hPa, tau= 9.4-23', &
         'pc= 180-310hPa, tau= 9.4-23', 'pc= 310-440hPa, tau= 9.4-23', &
         'pc= 440-560hPa, tau= 9.4-23', 'pc= 560-680hPa, tau= 9.4-23', &
         'pc= 680-800hPa, tau= 9.4-23', 'pc< 50hPa, tau= 23-60', &
         'pc= 50-180hPa, tau= 23-60', 'pc= 180-310hPa, tau= 23-60', &
         'pc= 310-440hPa, tau= 23-60', 'pc= 440-560hPa, tau= 23-60', &
         'pc= 560-680hPa, tau= 23-60', 'pc= 680-800hPa, tau= 23-60', &
         'pc< 50hPa, tau> 60.', 'pc= 50-180hPa, tau> 60.', &
         'pc= 180-310hPa, tau> 60.', 'pc= 310-440hPa, tau> 60.', &
         'pc= 440-560hPa, tau> 60.', 'pc= 560-680hPa, tau> 60.', &
         'pc= 680-800hPa, tau> 60.'/

    !IM ISCCP simulator v3.4

    integer nid_hf, nid_hf3d
    save nid_hf, nid_hf3d

    INTEGER        longcles
    PARAMETER    ( longcles = 20 )

    ! Variables propres a la physique

    INTEGER, save:: radpas
    ! (Radiative transfer computations are made every "radpas" call to
    ! "physiq".)

    REAL radsol(klon)
    SAVE radsol               ! bilan radiatif au sol calcule par code radiatif

    INTEGER, SAVE:: itap ! number of calls to "physiq"

    REAL ftsol(klon, nbsrf)
    SAVE ftsol                  ! temperature du sol

    REAL ftsoil(klon, nsoilmx, nbsrf)
    SAVE ftsoil                 ! temperature dans le sol

    REAL fevap(klon, nbsrf)
    SAVE fevap                 ! evaporation
    REAL fluxlat(klon, nbsrf)
    SAVE fluxlat

    REAL fqsurf(klon, nbsrf)
    SAVE fqsurf                 ! humidite de l'air au contact de la surface

    REAL qsol(klon)
    SAVE qsol                  ! hauteur d'eau dans le sol

    REAL fsnow(klon, nbsrf)
    SAVE fsnow                  ! epaisseur neigeuse

    REAL falbe(klon, nbsrf)
    SAVE falbe                  ! albedo par type de surface
    REAL falblw(klon, nbsrf)
    SAVE falblw                 ! albedo par type de surface

    ! Paramètres de l'orographie à l'échelle sous-maille (OESM) :
    REAL, save:: zmea(klon) ! orographie moyenne
    REAL, save:: zstd(klon) ! deviation standard de l'OESM
    REAL, save:: zsig(klon) ! pente de l'OESM
    REAL, save:: zgam(klon) ! anisotropie de l'OESM
    REAL, save:: zthe(klon) ! orientation de l'OESM
    REAL, save:: zpic(klon) ! Maximum de l'OESM
    REAL, save:: zval(klon) ! Minimum de l'OESM
    REAL, save:: rugoro(klon) ! longueur de rugosite de l'OESM

    REAL zulow(klon), zvlow(klon)

    INTEGER igwd, idx(klon), itest(klon)

    REAL agesno(klon, nbsrf)
    SAVE agesno                 ! age de la neige

    REAL run_off_lic_0(klon)
    SAVE run_off_lic_0
    !KE43
    ! Variables liees a la convection de K. Emanuel (sb):

    REAL bas, top             ! cloud base and top levels
    SAVE bas
    SAVE top

    REAL Ma(klon, llm)        ! undilute upward mass flux
    SAVE Ma
    REAL qcondc(klon, llm)    ! in-cld water content from convect
    SAVE qcondc 
    REAL ema_work1(klon, llm), ema_work2(klon, llm)
    SAVE ema_work1, ema_work2

    REAL wd(klon) ! sb
    SAVE wd       ! sb

    ! Variables locales pour la couche limite (al1):

    ! Variables locales:

    REAL cdragh(klon) ! drag coefficient pour T and Q
    REAL cdragm(klon) ! drag coefficient pour vent

    !AA  Pour phytrac 
    REAL ycoefh(klon, llm)    ! coef d'echange pour phytrac
    REAL yu1(klon)            ! vents dans la premiere couche U
    REAL yv1(klon)            ! vents dans la premiere couche V
    REAL ffonte(klon, nbsrf)    !Flux thermique utilise pour fondre la neige
    REAL fqcalving(klon, nbsrf) !Flux d'eau "perdue" par la surface 
    !                               !et necessaire pour limiter la
    !                               !hauteur de neige, en kg/m2/s
    REAL zxffonte(klon), zxfqcalving(klon)

    REAL pfrac_impa(klon, llm)! Produits des coefs lessivage impaction
    save pfrac_impa
    REAL pfrac_nucl(klon, llm)! Produits des coefs lessivage nucleation
    save pfrac_nucl
    REAL pfrac_1nucl(klon, llm)! Produits des coefs lessi nucl (alpha = 1)
    save pfrac_1nucl
    REAL frac_impa(klon, llm) ! fractions d'aerosols lessivees (impaction)
    REAL frac_nucl(klon, llm) ! idem (nucleation)

    !AA
    REAL rain_fall(klon) ! pluie
    REAL snow_fall(klon) ! neige
    save snow_fall, rain_fall
    !IM cf FH pour Tiedtke 080604
    REAL rain_tiedtke(klon), snow_tiedtke(klon)

    REAL evap(klon), devap(klon) ! evaporation et sa derivee
    REAL sens(klon), dsens(klon) ! chaleur sensible et sa derivee
    REAL dlw(klon)    ! derivee infra rouge
    SAVE dlw
    REAL bils(klon) ! bilan de chaleur au sol
    REAL fder(klon) ! Derive de flux (sensible et latente) 
    save fder
    REAL ve(klon) ! integr. verticale du transport meri. de l'energie
    REAL vq(klon) ! integr. verticale du transport meri. de l'eau
    REAL ue(klon) ! integr. verticale du transport zonal de l'energie
    REAL uq(klon) ! integr. verticale du transport zonal de l'eau

    REAL frugs(klon, nbsrf) ! longueur de rugosite
    save frugs
    REAL zxrugs(klon) ! longueur de rugosite

    ! Conditions aux limites

    INTEGER julien

    INTEGER, SAVE:: lmt_pas ! number of time steps of "physics" per day
    REAL pctsrf(klon, nbsrf)
    !IM
    REAL pctsrf_new(klon, nbsrf) !pourcentage surfaces issus d'ORCHIDEE

    SAVE pctsrf                 ! sous-fraction du sol
    REAL albsol(klon)
    SAVE albsol                 ! albedo du sol total
    REAL albsollw(klon)
    SAVE albsollw                 ! albedo du sol total

    REAL, SAVE:: wo(klon, llm) ! column density of ozone in a cell, in kDU

    ! Declaration des procedures appelees

    EXTERNAL alboc     ! calculer l'albedo sur ocean
    EXTERNAL ajsec     ! ajustement sec
    EXTERNAL clmain    ! couche limite 
    !KE43
    EXTERNAL conema3  ! convect4.3
    EXTERNAL fisrtilp  ! schema de condensation a grande echelle (pluie)
    EXTERNAL nuage     ! calculer les proprietes radiatives
    EXTERNAL ozonecm   ! prescrire l'ozone
    EXTERNAL radlwsw   ! rayonnements solaire et infrarouge
    EXTERNAL transp    ! transport total de l'eau et de l'energie

    ! Variables locales

    real clwcon(klon, llm), rnebcon(klon, llm)
    real clwcon0(klon, llm), rnebcon0(klon, llm)

    save rnebcon, clwcon

    REAL rhcl(klon, llm)    ! humiditi relative ciel clair
    REAL dialiq(klon, llm)  ! eau liquide nuageuse
    REAL diafra(klon, llm)  ! fraction nuageuse
    REAL cldliq(klon, llm)  ! eau liquide nuageuse
    REAL cldfra(klon, llm)  ! fraction nuageuse
    REAL cldtau(klon, llm)  ! epaisseur optique
    REAL cldemi(klon, llm)  ! emissivite infrarouge

    REAL fluxq(klon, llm, nbsrf)   ! flux turbulent d'humidite
    REAL fluxt(klon, llm, nbsrf)   ! flux turbulent de chaleur
    REAL fluxu(klon, llm, nbsrf)   ! flux turbulent de vitesse u
    REAL fluxv(klon, llm, nbsrf)   ! flux turbulent de vitesse v

    REAL zxfluxt(klon, llm)
    REAL zxfluxq(klon, llm)
    REAL zxfluxu(klon, llm)
    REAL zxfluxv(klon, llm)

    REAL heat(klon, llm)    ! chauffage solaire
    REAL heat0(klon, llm)   ! chauffage solaire ciel clair
    REAL cool(klon, llm)    ! refroidissement infrarouge
    REAL cool0(klon, llm)   ! refroidissement infrarouge ciel clair
    REAL topsw(klon), toplw(klon), solsw(klon), sollw(klon)
    real sollwdown(klon)    ! downward LW flux at surface
    REAL topsw0(klon), toplw0(klon), solsw0(klon), sollw0(klon)
    REAL albpla(klon)
    REAL fsollw(klon, nbsrf)   ! bilan flux IR pour chaque sous surface
    REAL fsolsw(klon, nbsrf)   ! flux solaire absorb. pour chaque sous surface
    ! Le rayonnement n'est pas calcule tous les pas, il faut donc
    !                      sauvegarder les sorties du rayonnement
    SAVE  heat, cool, albpla, topsw, toplw, solsw, sollw, sollwdown
    SAVE  topsw0, toplw0, solsw0, sollw0, heat0, cool0

    INTEGER itaprad
    SAVE itaprad

    REAL conv_q(klon, llm) ! convergence de l'humidite (kg/kg/s)
    REAL conv_t(klon, llm) ! convergence de la temperature(K/s)

    REAL cldl(klon), cldm(klon), cldh(klon) !nuages bas, moyen et haut
    REAL cldt(klon), cldq(klon) !nuage total, eau liquide integree

    REAL zxtsol(klon), zxqsurf(klon), zxsnow(klon), zxfluxlat(klon)

    REAL dist, rmu0(klon), fract(klon)
    REAL zdtime ! pas de temps du rayonnement (s)
    real zlongi

    REAL z_avant(klon), z_apres(klon), z_factor(klon)
    LOGICAL zx_ajustq

    REAL za, zb
    REAL zx_t, zx_qs, zdelta, zcor, zlvdcp, zlsdcp
    real zqsat(klon, llm)
    INTEGER i, k, iq, nsrf
    REAL t_coup
    PARAMETER (t_coup=234.0)

    REAL zphi(klon, llm)

    !IM cf. AM Variables locales pour la CLA (hbtm2)

    REAL pblh(klon, nbsrf)           ! Hauteur de couche limite
    REAL plcl(klon, nbsrf)           ! Niveau de condensation de la CLA
    REAL capCL(klon, nbsrf)          ! CAPE de couche limite
    REAL oliqCL(klon, nbsrf)          ! eau_liqu integree de couche limite
    REAL cteiCL(klon, nbsrf)          ! cloud top instab. crit. couche limite
    REAL pblt(klon, nbsrf)          ! T a la Hauteur de couche limite
    REAL therm(klon, nbsrf)
    REAL trmb1(klon, nbsrf)          ! deep_cape
    REAL trmb2(klon, nbsrf)          ! inhibition 
    REAL trmb3(klon, nbsrf)          ! Point Omega
    ! Grdeurs de sorties
    REAL s_pblh(klon), s_lcl(klon), s_capCL(klon)
    REAL s_oliqCL(klon), s_cteiCL(klon), s_pblt(klon)
    REAL s_therm(klon), s_trmb1(klon), s_trmb2(klon)
    REAL s_trmb3(klon)

    ! Variables locales pour la convection de K. Emanuel (sb):

    REAL upwd(klon, llm)      ! saturated updraft mass flux
    REAL dnwd(klon, llm)      ! saturated downdraft mass flux
    REAL dnwd0(klon, llm)     ! unsaturated downdraft mass flux
    REAL tvp(klon, llm)       ! virtual temp of lifted parcel
    REAL cape(klon)           ! CAPE
    SAVE cape

    REAL pbase(klon)          ! cloud base pressure
    SAVE pbase
    REAL bbase(klon)          ! cloud base buoyancy
    SAVE bbase
    REAL rflag(klon)          ! flag fonctionnement de convect
    INTEGER iflagctrl(klon)          ! flag fonctionnement de convect
    ! -- convect43:
    INTEGER ntra              ! nb traceurs pour convect4.3
    REAL dtvpdt1(klon, llm), dtvpdq1(klon, llm)
    REAL dplcldt(klon), dplcldr(klon)

    ! Variables du changement

    ! con: convection
    ! lsc: condensation a grande echelle (Large-Scale-Condensation)
    ! ajs: ajustement sec
    ! eva: evaporation de l'eau liquide nuageuse
    ! vdf: couche limite (Vertical DiFfusion)
    REAL d_t_con(klon, llm), d_q_con(klon, llm)
    REAL d_u_con(klon, llm), d_v_con(klon, llm)
    REAL d_t_lsc(klon, llm), d_q_lsc(klon, llm), d_ql_lsc(klon, llm)
    REAL d_t_ajs(klon, llm), d_q_ajs(klon, llm)
    REAL d_u_ajs(klon, llm), d_v_ajs(klon, llm)
    REAL rneb(klon, llm)

    REAL pmfu(klon, llm), pmfd(klon, llm)
    REAL pen_u(klon, llm), pen_d(klon, llm)
    REAL pde_u(klon, llm), pde_d(klon, llm)
    INTEGER kcbot(klon), kctop(klon), kdtop(klon)
    REAL pmflxr(klon, llm+1), pmflxs(klon, llm+1)
    REAL prfl(klon, llm+1), psfl(klon, llm+1)

    INTEGER ibas_con(klon), itop_con(klon)

    SAVE ibas_con, itop_con

    REAL rain_con(klon), rain_lsc(klon)
    REAL snow_con(klon), snow_lsc(klon)
    REAL d_ts(klon, nbsrf)

    REAL d_u_vdf(klon, llm), d_v_vdf(klon, llm)
    REAL d_t_vdf(klon, llm), d_q_vdf(klon, llm)

    REAL d_u_oro(klon, llm), d_v_oro(klon, llm)
    REAL d_t_oro(klon, llm)
    REAL d_u_lif(klon, llm), d_v_lif(klon, llm)
    REAL d_t_lif(klon, llm)

    REAL ratqs(klon, llm), ratqss(klon, llm), ratqsc(klon, llm)
    real ratqsbas, ratqshaut
    save ratqsbas, ratqshaut, ratqs

    ! Parametres lies au nouveau schema de nuages (SB, PDF)
    real, save:: fact_cldcon
    real, save:: facttemps
    logical ok_newmicro
    save ok_newmicro
    real facteur

    integer iflag_cldcon
    save iflag_cldcon

    logical ptconv(klon, llm)

    ! Variables locales pour effectuer les appels en serie

    REAL t_seri(klon, llm), q_seri(klon, llm)
    REAL ql_seri(klon, llm), qs_seri(klon, llm)
    REAL u_seri(klon, llm), v_seri(klon, llm)

    REAL tr_seri(klon, llm, nbtr)
    REAL d_tr(klon, llm, nbtr)

    REAL zx_rh(klon, llm)

    REAL zustrdr(klon), zvstrdr(klon)
    REAL zustrli(klon), zvstrli(klon)
    REAL zustrph(klon), zvstrph(klon)
    REAL aam, torsfc

    REAL dudyn(iim+1, jjm + 1, llm)

    REAL zx_tmp_fi2d(klon)      ! variable temporaire grille physique
    REAL zx_tmp_fi3d(klon, llm) ! variable temporaire pour champs 3D 

    REAL zx_tmp_2d(iim, jjm + 1), zx_tmp_3d(iim, jjm + 1, llm)

    INTEGER, SAVE:: nid_day, nid_ins

    REAL ve_lay(klon, llm) ! transport meri. de l'energie a chaque niveau vert.
    REAL vq_lay(klon, llm) ! transport meri. de l'eau a chaque niveau vert.
    REAL ue_lay(klon, llm) ! transport zonal de l'energie a chaque niveau vert.
    REAL uq_lay(klon, llm) ! transport zonal de l'eau a chaque niveau vert.

    REAL zsto

    character(len=20) modname
    character(len=80) abort_message
    logical ok_sync
    real date0

    !     Variables liees au bilan d'energie et d'enthalpi
    REAL ztsol(klon)
    REAL      d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec
    REAL      d_h_vcol_phy
    REAL      fs_bound, fq_bound
    SAVE      d_h_vcol_phy
    REAL      zero_v(klon)
    CHARACTER(LEN=15) ztit
    INTEGER   ip_ebil  ! PRINT level for energy conserv. diag.
    SAVE      ip_ebil
    DATA      ip_ebil/0/
    INTEGER   if_ebil ! level for energy conserv. dignostics
    SAVE      if_ebil
    !+jld ec_conser
    REAL d_t_ec(klon, llm)    ! tendance du a la conersion Ec -> E thermique
    REAL ZRCPD
    !-jld ec_conser
    !IM: t2m, q2m, u10m, v10m
    REAL t2m(klon, nbsrf), q2m(klon, nbsrf)   !temperature, humidite a 2m
    REAL u10m(klon, nbsrf), v10m(klon, nbsrf) !vents a 10m
    REAL zt2m(klon), zq2m(klon)             !temp., hum. 2m moyenne s/ 1 maille
    REAL zu10m(klon), zv10m(klon)           !vents a 10m moyennes s/1 maille
    !jq   Aerosol effects (Johannes Quaas, 27/11/2003)
    REAL sulfate(klon, llm) ! SO4 aerosol concentration [ug/m3]

    REAL sulfate_pi(klon, llm)
    ! (SO4 aerosol concentration [ug/m3] (pre-industrial value))
    SAVE sulfate_pi

    REAL cldtaupi(klon, llm)
    ! (Cloud optical thickness for pre-industrial (pi) aerosols)

    REAL re(klon, llm)       ! Cloud droplet effective radius
    REAL fl(klon, llm)  ! denominator of re

    ! Aerosol optical properties
    REAL tau_ae(klon, llm, 2), piz_ae(klon, llm, 2)
    REAL cg_ae(klon, llm, 2)

    REAL topswad(klon), solswad(klon) ! Aerosol direct effect.
    ! ok_ade=T -ADE=topswad-topsw

    REAL topswai(klon), solswai(klon) ! Aerosol indirect effect.
    ! ok_aie=T ->
    !        ok_ade=T -AIE=topswai-topswad
    !        ok_ade=F -AIE=topswai-topsw

    REAL aerindex(klon)       ! POLDER aerosol index

    ! Parameters
    LOGICAL ok_ade, ok_aie    ! Apply aerosol (in)direct effects or not
    REAL bl95_b0, bl95_b1   ! Parameter in Boucher and Lohmann (1995)

    SAVE ok_ade, ok_aie, bl95_b0, bl95_b1
    SAVE u10m
    SAVE v10m
    SAVE t2m
    SAVE q2m
    SAVE ffonte
    SAVE fqcalving
    SAVE piz_ae
    SAVE tau_ae
    SAVE cg_ae
    SAVE rain_con
    SAVE snow_con
    SAVE topswai
    SAVE topswad
    SAVE solswai
    SAVE solswad
    SAVE d_u_con
    SAVE d_v_con
    SAVE rnebcon0
    SAVE clwcon0
    SAVE pblh
    SAVE plcl
    SAVE capCL
    SAVE oliqCL
    SAVE cteiCL
    SAVE pblt
    SAVE therm
    SAVE trmb1
    SAVE trmb2
    SAVE trmb3

    real zmasse(klon, llm) 
    ! (column-density of mass of air in a cell, in kg m-2)

    real, parameter:: dobson_u = 2.1415e-05 ! Dobson unit, in kg m-2

    !----------------------------------------------------------------

    modname = 'physiq'
    IF (if_ebil >= 1) THEN
       DO i=1, klon
          zero_v(i)=0.
       END DO
    END IF
    ok_sync=.TRUE.
    IF (nq  <  2) THEN
       abort_message = 'eaux vapeur et liquide sont indispensables'
       CALL abort_gcm(modname, abort_message, 1)
    ENDIF

    test_firstcal: IF (firstcal) THEN
       !  initialiser
       u10m=0.
       v10m=0.
       t2m=0.
       q2m=0.
       ffonte=0.
       fqcalving=0.
       piz_ae(:, :, :)=0.
       tau_ae(:, :, :)=0.
       cg_ae(:, :, :)=0.
       rain_con(:)=0.
       snow_con(:)=0.
       bl95_b0=0.
       bl95_b1=0.
       topswai(:)=0.
       topswad(:)=0.
       solswai(:)=0.
       solswad(:)=0.

       d_u_con = 0.0
       d_v_con = 0.0
       rnebcon0 = 0.0
       clwcon0 = 0.0
       rnebcon = 0.0
       clwcon = 0.0

       pblh   =0.        ! Hauteur de couche limite
       plcl   =0.        ! Niveau de condensation de la CLA
       capCL  =0.        ! CAPE de couche limite
       oliqCL =0.        ! eau_liqu integree de couche limite
       cteiCL =0.        ! cloud top instab. crit. couche limite
       pblt   =0.        ! T a la Hauteur de couche limite
       therm  =0.
       trmb1  =0.        ! deep_cape
       trmb2  =0.        ! inhibition 
       trmb3  =0.        ! Point Omega

       IF (if_ebil >= 1) d_h_vcol_phy=0.

       ! appel a la lecture du run.def physique

       call conf_phys(ocean, ok_veget, ok_journe, ok_mensuel, &
            ok_instan, fact_cldcon, facttemps, ok_newmicro, &
            iflag_cldcon, ratqsbas, ratqshaut, if_ebil, &
            ok_ade, ok_aie,  &
            bl95_b0, bl95_b1, &
            iflag_thermals, nsplit_thermals)

       ! Initialiser les compteurs:

       frugs = 0.
       itap = 0
       itaprad = 0
       CALL phyetat0("startphy.nc", pctsrf, ftsol, ftsoil, ocean, tslab, &
            seaice, fqsurf, qsol, fsnow, &
            falbe, falblw, fevap, rain_fall, snow_fall, solsw, sollwdown, &
            dlw, radsol, frugs, agesno, &
            zmea, zstd, zsig, zgam, zthe, zpic, zval, &
            t_ancien, q_ancien, ancien_ok, rnebcon, ratqs, clwcon,  &
            run_off_lic_0)

       !   ATTENTION : il faudra a terme relire q2 dans l'etat initial
       q2(:, :, :)=1.e-8

       radpas = NINT( 86400. / pdtphys / nbapp_rad)

       ! on remet le calendrier a zero
       IF (raz_date) itau_phy = 0

       PRINT *, 'cycle_diurne = ', cycle_diurne

       IF(ocean.NE.'force ') THEN
          ok_ocean=.TRUE.
       ENDIF

       CALL printflag(radpas, ok_ocean, ok_oasis, ok_journe, ok_instan, &
            ok_region)

       IF (pdtphys*REAL(radpas).GT.21600..AND.cycle_diurne) THEN 
          print *,'Nbre d appels au rayonnement insuffisant'
          print *,"Au minimum 4 appels par jour si cycle diurne"
          abort_message='Nbre d appels au rayonnement insuffisant'
          call abort_gcm(modname, abort_message, 1)
       ENDIF
       print *,"Clef pour la convection, iflag_con=", iflag_con
       print *,"Clef pour le driver de la convection, ok_cvl=", &
            ok_cvl

       ! Initialisation pour la convection de K.E. (sb):
       IF (iflag_con >= 3) THEN

          print *,"*** Convection de Kerry Emanuel 4.3  "

          !IM15/11/02 rajout initialisation ibas_con, itop_con cf. SB =>BEG
          DO i = 1, klon
             ibas_con(i) = 1
             itop_con(i) = 1
          ENDDO
          !IM15/11/02 rajout initialisation ibas_con, itop_con cf. SB =>END

       ENDIF

       IF (ok_orodr) THEN
          rugoro = MAX(1e-5, zstd * zsig / 2)
          CALL SUGWD(klon, llm, paprs, pplay)
       else
          rugoro = 0.
       ENDIF

       lmt_pas = NINT(86400. / pdtphys)  ! tous les jours
       print *, 'Number of time steps of "physics" per day: ', lmt_pas

       ecrit_ins = NINT(ecrit_ins/pdtphys)
       ecrit_hf = NINT(ecrit_hf/pdtphys)
       ecrit_mth = NINT(ecrit_mth/pdtphys)
       ecrit_tra = NINT(86400.*ecrit_tra/pdtphys)
       ecrit_reg = NINT(ecrit_reg/pdtphys)

       ! Initialiser le couplage si necessaire

       npas = 0
       nexca = 0

       print *,'AVANT HIST IFLAG_CON=', iflag_con

       !   Initialisation des sorties

       call ini_histhf(pdtphys, presnivs, nid_hf, nid_hf3d)
       call ini_histday(pdtphys, presnivs, ok_journe, nid_day, nq)
       call ini_histins(pdtphys, presnivs, ok_instan, nid_ins)
       CALL ymds2ju(annee_ref, 1, int(day_ref), 0., date0)
       !XXXPB Positionner date0 pour initialisation de ORCHIDEE
       WRITE(*, *) 'physiq date0 : ', date0
    ENDIF test_firstcal

    ! Mettre a zero des variables de sortie (pour securite)

    DO i = 1, klon
       d_ps(i) = 0.0
    ENDDO
    DO k = 1, llm
       DO i = 1, klon
          d_t(i, k) = 0.0
          d_u(i, k) = 0.0
          d_v(i, k) = 0.0
       ENDDO
    ENDDO
    DO iq = 1, nq
       DO k = 1, llm
          DO i = 1, klon
             d_qx(i, k, iq) = 0.0
          ENDDO
       ENDDO
    ENDDO
    da=0.
    mp=0.
    phi(:, :, :)=0.

    ! Ne pas affecter les valeurs entrees de u, v, h, et q

    DO k = 1, llm
       DO i = 1, klon
          t_seri(i, k)  = t(i, k)
          u_seri(i, k)  = u(i, k)
          v_seri(i, k)  = v(i, k)
          q_seri(i, k)  = qx(i, k, ivap)
          ql_seri(i, k) = qx(i, k, iliq)
          qs_seri(i, k) = 0.
       ENDDO
    ENDDO
    IF (nq >= 3) THEN
       tr_seri(:, :, :nq-2) = qx(:, :, 3:nq)
    ELSE
       tr_seri(:, :, 1) = 0.
    ENDIF

    DO i = 1, klon
       ztsol(i) = 0.
    ENDDO
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          ztsol(i) = ztsol(i) + ftsol(i, nsrf)*pctsrf(i, nsrf)
       ENDDO
    ENDDO

    IF (if_ebil >= 1) THEN 
       ztit='after dynamic'
       CALL diagetpq(airephy, ztit, ip_ebil, 1, 1, pdtphys &
            , t_seri, q_seri, ql_seri, qs_seri, u_seri, v_seri, paprs &
            , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
       !     Comme les tendances de la physique sont ajoute dans la dynamique, 
       !     on devrait avoir que la variation d'entalpie par la dynamique
       !     est egale a la variation de la physique au pas de temps precedent.
       !     Donc la somme de ces 2 variations devrait etre nulle.
       call diagphy(airephy, ztit, ip_ebil &
            , zero_v, zero_v, zero_v, zero_v, zero_v &
            , zero_v, zero_v, zero_v, ztsol &
            , d_h_vcol+d_h_vcol_phy, d_qt, 0. &
            , fs_bound, fq_bound )
    END IF

    ! Diagnostiquer la tendance dynamique

    IF (ancien_ok) THEN
       DO k = 1, llm
          DO i = 1, klon
             d_t_dyn(i, k) = (t_seri(i, k)-t_ancien(i, k))/pdtphys
             d_q_dyn(i, k) = (q_seri(i, k)-q_ancien(i, k))/pdtphys
          ENDDO
       ENDDO
    ELSE
       DO k = 1, llm
          DO i = 1, klon
             d_t_dyn(i, k) = 0.0
             d_q_dyn(i, k) = 0.0
          ENDDO
       ENDDO
       ancien_ok = .TRUE.
    ENDIF

    ! Ajouter le geopotentiel du sol:

    DO k = 1, llm
       DO i = 1, klon
          zphi(i, k) = pphi(i, k) + pphis(i)
       ENDDO
    ENDDO

    ! Verifier les temperatures

    CALL hgardfou(t_seri, ftsol)

    ! Incrementer le compteur de la physique

    itap = itap + 1
    julien = MOD(NINT(rdayvrai), 360)
    if (julien == 0) julien = 360

    forall (k = 1: llm) zmasse(:, k) = (paprs(:, k)-paprs(:, k+1)) / rg

    ! Mettre en action les conditions aux limites (albedo, sst, etc.).
    ! Prescrire l'ozone et calculer l'albedo sur l'ocean.

!!$    if (nq >= 5) then
!!$       wo = qx(:, :, 5) * zmasse / dobson_u / 1e3
!!$    else IF (MOD(itap - 1, lmt_pas) == 0) THEN
    IF (MOD(itap - 1, lmt_pas) == 0) THEN
       CALL ozonecm(REAL(julien), rlat, paprs, wo)
    ENDIF

    ! Re-evaporer l'eau liquide nuageuse

    DO k = 1, llm  ! re-evaporation de l'eau liquide nuageuse
       DO i = 1, klon
          zlvdcp=RLVTT/RCPD/(1.0+RVTMP2*q_seri(i, k))
          zlsdcp=RLVTT/RCPD/(1.0+RVTMP2*q_seri(i, k))
          zdelta = MAX(0., SIGN(1., RTT-t_seri(i, k)))
          zb = MAX(0.0, ql_seri(i, k))
          za = - MAX(0.0, ql_seri(i, k)) &
               * (zlvdcp*(1.-zdelta)+zlsdcp*zdelta)
          t_seri(i, k) = t_seri(i, k) + za
          q_seri(i, k) = q_seri(i, k) + zb
          ql_seri(i, k) = 0.0
       ENDDO
    ENDDO

    IF (if_ebil >= 2) THEN 
       ztit='after reevap'
       CALL diagetpq(airephy, ztit, ip_ebil, 2, 1, pdtphys &
            , t_seri, q_seri, ql_seri, qs_seri, u_seri, v_seri, paprs &
            , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
       call diagphy(airephy, ztit, ip_ebil &
            , zero_v, zero_v, zero_v, zero_v, zero_v &
            , zero_v, zero_v, zero_v, ztsol &
            , d_h_vcol, d_qt, d_ec &
            , fs_bound, fq_bound )

    END IF

    ! Appeler la diffusion verticale (programme de couche limite)

    DO i = 1, klon
       zxrugs(i) = 0.0
    ENDDO
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          frugs(i, nsrf) = MAX(frugs(i, nsrf), 0.000015)
       ENDDO
    ENDDO
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          zxrugs(i) = zxrugs(i) + frugs(i, nsrf)*pctsrf(i, nsrf)
       ENDDO
    ENDDO

    ! calculs necessaires au calcul de l'albedo dans l'interface

    CALL orbite(REAL(julien), zlongi, dist)
    IF (cycle_diurne) THEN
       zdtime = pdtphys * REAL(radpas)
       CALL zenang(zlongi, gmtime, zdtime, rmu0, fract)
    ELSE
       rmu0 = -999.999
    ENDIF

    !     Calcul de l'abedo moyen par maille
    albsol(:)=0.
    albsollw(:)=0.
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          albsol(i) = albsol(i) + falbe(i, nsrf) * pctsrf(i, nsrf)
          albsollw(i) = albsollw(i) + falblw(i, nsrf) * pctsrf(i, nsrf)
       ENDDO
    ENDDO

    !     Repartition sous maille des flux LW et SW
    ! Repartition du longwave par sous-surface linearisee

    DO nsrf = 1, nbsrf
       DO i = 1, klon
          fsollw(i, nsrf) = sollw(i) &
               + 4.0*RSIGMA*ztsol(i)**3 * (ztsol(i)-ftsol(i, nsrf))
          fsolsw(i, nsrf) = solsw(i)*(1.-falbe(i, nsrf))/(1.-albsol(i))
       ENDDO
    ENDDO

    fder = dlw

    CALL clmain(pdtphys, itap, date0, pctsrf, pctsrf_new, &
         t_seri, q_seri, u_seri, v_seri, &
         julien, rmu0, co2_ppm,  &
         ok_veget, ocean, npas, nexca, ftsol, &
         soil_model, cdmmax, cdhmax, &
         ksta, ksta_ter, ok_kzmin, ftsoil, qsol,  &
         paprs, pplay, fsnow, fqsurf, fevap, falbe, falblw, &
         fluxlat, rain_fall, snow_fall, &
         fsolsw, fsollw, sollwdown, fder, &
         rlon, rlat, cuphy, cvphy, frugs, &
         firstcal, lafin, agesno, rugoro, &
         d_t_vdf, d_q_vdf, d_u_vdf, d_v_vdf, d_ts, &
         fluxt, fluxq, fluxu, fluxv, cdragh, cdragm, &
         q2, dsens, devap, &
         ycoefh, yu1, yv1, t2m, q2m, u10m, v10m, &
         pblh, capCL, oliqCL, cteiCL, pblT, &
         therm, trmb1, trmb2, trmb3, plcl, &
         fqcalving, ffonte, run_off_lic_0, &
         fluxo, fluxg, tslab, seaice)

    !XXX Incrementation des flux

    zxfluxt=0.
    zxfluxq=0.
    zxfluxu=0.
    zxfluxv=0.
    DO nsrf = 1, nbsrf
       DO k = 1, llm
          DO i = 1, klon
             zxfluxt(i, k) = zxfluxt(i, k) +  &
                  fluxt(i, k, nsrf) * pctsrf( i, nsrf)
             zxfluxq(i, k) = zxfluxq(i, k) +  &
                  fluxq(i, k, nsrf) * pctsrf( i, nsrf)
             zxfluxu(i, k) = zxfluxu(i, k) +  &
                  fluxu(i, k, nsrf) * pctsrf( i, nsrf)
             zxfluxv(i, k) = zxfluxv(i, k) +  &
                  fluxv(i, k, nsrf) * pctsrf( i, nsrf)
          END DO
       END DO
    END DO
    DO i = 1, klon
       sens(i) = - zxfluxt(i, 1) ! flux de chaleur sensible au sol
       evap(i) = - zxfluxq(i, 1) ! flux d'evaporation au sol
       fder(i) = dlw(i) + dsens(i) + devap(i)
    ENDDO

    DO k = 1, llm
       DO i = 1, klon
          t_seri(i, k) = t_seri(i, k) + d_t_vdf(i, k)
          q_seri(i, k) = q_seri(i, k) + d_q_vdf(i, k)
          u_seri(i, k) = u_seri(i, k) + d_u_vdf(i, k)
          v_seri(i, k) = v_seri(i, k) + d_v_vdf(i, k)
       ENDDO
    ENDDO

    IF (if_ebil >= 2) THEN 
       ztit='after clmain'
       CALL diagetpq(airephy, ztit, ip_ebil, 2, 2, pdtphys &
            , t_seri, q_seri, ql_seri, qs_seri, u_seri, v_seri, paprs &
            , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
       call diagphy(airephy, ztit, ip_ebil &
            , zero_v, zero_v, zero_v, zero_v, sens &
            , evap, zero_v, zero_v, ztsol &
            , d_h_vcol, d_qt, d_ec &
            , fs_bound, fq_bound )
    END IF

    ! Incrementer la temperature du sol

    DO i = 1, klon
       zxtsol(i) = 0.0
       zxfluxlat(i) = 0.0

       zt2m(i) = 0.0
       zq2m(i) = 0.0
       zu10m(i) = 0.0
       zv10m(i) = 0.0
       zxffonte(i) = 0.0
       zxfqcalving(i) = 0.0

       s_pblh(i) = 0.0 
       s_lcl(i) = 0.0 
       s_capCL(i) = 0.0
       s_oliqCL(i) = 0.0
       s_cteiCL(i) = 0.0
       s_pblT(i) = 0.0
       s_therm(i) = 0.0
       s_trmb1(i) = 0.0
       s_trmb2(i) = 0.0
       s_trmb3(i) = 0.0

       IF ( abs( pctsrf(i, is_ter) + pctsrf(i, is_lic) +  &
            pctsrf(i, is_oce) + pctsrf(i, is_sic)  - 1.) .GT. EPSFRA)  &
            THEN 
          WRITE(*, *) 'physiq : pb sous surface au point ', i,  &
               pctsrf(i, 1 : nbsrf)
       ENDIF
    ENDDO
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          ftsol(i, nsrf) = ftsol(i, nsrf) + d_ts(i, nsrf)
          zxtsol(i) = zxtsol(i) + ftsol(i, nsrf)*pctsrf(i, nsrf)
          zxfluxlat(i) = zxfluxlat(i) + fluxlat(i, nsrf)*pctsrf(i, nsrf)

          zt2m(i) = zt2m(i) + t2m(i, nsrf)*pctsrf(i, nsrf)
          zq2m(i) = zq2m(i) + q2m(i, nsrf)*pctsrf(i, nsrf)
          zu10m(i) = zu10m(i) + u10m(i, nsrf)*pctsrf(i, nsrf)
          zv10m(i) = zv10m(i) + v10m(i, nsrf)*pctsrf(i, nsrf)
          zxffonte(i) = zxffonte(i) + ffonte(i, nsrf)*pctsrf(i, nsrf)
          zxfqcalving(i) = zxfqcalving(i) +  &
               fqcalving(i, nsrf)*pctsrf(i, nsrf)
          s_pblh(i) = s_pblh(i) + pblh(i, nsrf)*pctsrf(i, nsrf)
          s_lcl(i) = s_lcl(i) + plcl(i, nsrf)*pctsrf(i, nsrf)
          s_capCL(i) = s_capCL(i) + capCL(i, nsrf) *pctsrf(i, nsrf)
          s_oliqCL(i) = s_oliqCL(i) + oliqCL(i, nsrf) *pctsrf(i, nsrf)
          s_cteiCL(i) = s_cteiCL(i) + cteiCL(i, nsrf) *pctsrf(i, nsrf)
          s_pblT(i) = s_pblT(i) + pblT(i, nsrf) *pctsrf(i, nsrf)
          s_therm(i) = s_therm(i) + therm(i, nsrf) *pctsrf(i, nsrf)
          s_trmb1(i) = s_trmb1(i) + trmb1(i, nsrf) *pctsrf(i, nsrf)
          s_trmb2(i) = s_trmb2(i) + trmb2(i, nsrf) *pctsrf(i, nsrf)
          s_trmb3(i) = s_trmb3(i) + trmb3(i, nsrf) *pctsrf(i, nsrf)
       ENDDO
    ENDDO

    ! Si une sous-fraction n'existe pas, elle prend la temp. moyenne

    DO nsrf = 1, nbsrf
       DO i = 1, klon
          IF (pctsrf(i, nsrf)  <  epsfra) ftsol(i, nsrf) = zxtsol(i)

          IF (pctsrf(i, nsrf)  <  epsfra) t2m(i, nsrf) = zt2m(i)
          IF (pctsrf(i, nsrf)  <  epsfra) q2m(i, nsrf) = zq2m(i)
          IF (pctsrf(i, nsrf)  <  epsfra) u10m(i, nsrf) = zu10m(i)
          IF (pctsrf(i, nsrf)  <  epsfra) v10m(i, nsrf) = zv10m(i)
          IF (pctsrf(i, nsrf)  <  epsfra) ffonte(i, nsrf) = zxffonte(i)
          IF (pctsrf(i, nsrf)  <  epsfra)  &
               fqcalving(i, nsrf) = zxfqcalving(i)
          IF (pctsrf(i, nsrf)  <  epsfra) pblh(i, nsrf)=s_pblh(i)
          IF (pctsrf(i, nsrf)  <  epsfra) plcl(i, nsrf)=s_lcl(i)
          IF (pctsrf(i, nsrf)  <  epsfra) capCL(i, nsrf)=s_capCL(i)
          IF (pctsrf(i, nsrf)  <  epsfra) oliqCL(i, nsrf)=s_oliqCL(i)
          IF (pctsrf(i, nsrf)  <  epsfra) cteiCL(i, nsrf)=s_cteiCL(i)
          IF (pctsrf(i, nsrf)  <  epsfra) pblT(i, nsrf)=s_pblT(i)
          IF (pctsrf(i, nsrf)  <  epsfra) therm(i, nsrf)=s_therm(i)
          IF (pctsrf(i, nsrf)  <  epsfra) trmb1(i, nsrf)=s_trmb1(i)
          IF (pctsrf(i, nsrf)  <  epsfra) trmb2(i, nsrf)=s_trmb2(i)
          IF (pctsrf(i, nsrf)  <  epsfra) trmb3(i, nsrf)=s_trmb3(i)
       ENDDO
    ENDDO

    ! Calculer la derive du flux infrarouge

    DO i = 1, klon
       dlw(i) = - 4.0*RSIGMA*zxtsol(i)**3 
    ENDDO

    ! Appeler la convection (au choix)

    DO k = 1, llm
       DO i = 1, klon
          conv_q(i, k) = d_q_dyn(i, k)  &
               + d_q_vdf(i, k)/pdtphys
          conv_t(i, k) = d_t_dyn(i, k)  &
               + d_t_vdf(i, k)/pdtphys
       ENDDO
    ENDDO
    IF (check) THEN
       za = qcheck(klon, llm, paprs, q_seri, ql_seri, airephy)
       print *, "avantcon=", za
    ENDIF
    zx_ajustq = .FALSE.
    IF (iflag_con == 2) zx_ajustq=.TRUE.
    IF (zx_ajustq) THEN
       DO i = 1, klon
          z_avant(i) = 0.0
       ENDDO
       DO k = 1, llm
          DO i = 1, klon
             z_avant(i) = z_avant(i) + (q_seri(i, k)+ql_seri(i, k)) &
                  *zmasse(i, k)
          ENDDO
       ENDDO
    ENDIF
    IF (iflag_con == 1) THEN
       stop 'reactiver le call conlmd dans physiq.F'
    ELSE IF (iflag_con == 2) THEN
       CALL conflx(pdtphys, paprs, pplay, t_seri, q_seri, &
            conv_t, conv_q, zxfluxq(1, 1), omega, &
            d_t_con, d_q_con, rain_con, snow_con, &
            pmfu, pmfd, pen_u, pde_u, pen_d, pde_d, &
            kcbot, kctop, kdtop, pmflxr, pmflxs)
       WHERE (rain_con < 0.) rain_con = 0.
       WHERE (snow_con < 0.) snow_con = 0.
       DO i = 1, klon
          ibas_con(i) = llm+1 - kcbot(i)
          itop_con(i) = llm+1 - kctop(i)
       ENDDO
    ELSE IF (iflag_con >= 3) THEN
       ! nb of tracers for the KE convection:
       ! MAF la partie traceurs est faite dans phytrac
       ! on met ntra=1 pour limiter les appels mais on peut
       ! supprimer les calculs / ftra.
       ntra = 1
       ! Schema de convection modularise et vectorise:
       ! (driver commun aux versions 3 et 4)

       IF (ok_cvl) THEN ! new driver for convectL
          CALL concvl (iflag_con, pdtphys, paprs, pplay, t_seri, q_seri, &
               u_seri, v_seri, tr_seri, ntra, &
               ema_work1, ema_work2, &
               d_t_con, d_q_con, d_u_con, d_v_con, d_tr, &
               rain_con, snow_con, ibas_con, itop_con, &
               upwd, dnwd, dnwd0, &
               Ma, cape, tvp, iflagctrl, &
               pbase, bbase, dtvpdt1, dtvpdq1, dplcldt, dplcldr, qcondc, wd, &
               pmflxr, pmflxs, &
               da, phi, mp)

          clwcon0=qcondc
          pmfu=upwd+dnwd
       ELSE ! ok_cvl
          ! MAF conema3 ne contient pas les traceurs
          CALL conema3 (pdtphys, paprs, pplay, t_seri, q_seri, &
               u_seri, v_seri, tr_seri, ntra, &
               ema_work1, ema_work2, &
               d_t_con, d_q_con, d_u_con, d_v_con, d_tr, &
               rain_con, snow_con, ibas_con, itop_con, &
               upwd, dnwd, dnwd0, bas, top, &
               Ma, cape, tvp, rflag, &
               pbase &
               , bbase, dtvpdt1, dtvpdq1, dplcldt, dplcldr &
               , clwcon0)
       ENDIF ! ok_cvl

       IF (.NOT. ok_gust) THEN
          do i = 1, klon
             wd(i)=0.0
          enddo
       ENDIF

       ! Calcul des proprietes des nuages convectifs

       DO k = 1, llm
          DO i = 1, klon
             zx_t = t_seri(i, k)
             IF (thermcep) THEN
                zdelta = MAX(0., SIGN(1., rtt-zx_t))
                zx_qs  = r2es * FOEEW(zx_t, zdelta)/pplay(i, k)
                zx_qs  = MIN(0.5, zx_qs)
                zcor   = 1./(1.-retv*zx_qs)
                zx_qs  = zx_qs*zcor
             ELSE
                IF (zx_t < t_coup) THEN
                   zx_qs = qsats(zx_t)/pplay(i, k)
                ELSE
                   zx_qs = qsatl(zx_t)/pplay(i, k)
                ENDIF
             ENDIF
             zqsat(i, k)=zx_qs
          ENDDO
       ENDDO

       !   calcul des proprietes des nuages convectifs
       clwcon0=fact_cldcon*clwcon0
       call clouds_gno &
            (klon, llm, q_seri, zqsat, clwcon0, ptconv, ratqsc, rnebcon0)
    ELSE
       print *, "iflag_con non-prevu", iflag_con
       stop 1
    ENDIF

    DO k = 1, llm
       DO i = 1, klon
          t_seri(i, k) = t_seri(i, k) + d_t_con(i, k)
          q_seri(i, k) = q_seri(i, k) + d_q_con(i, k)
          u_seri(i, k) = u_seri(i, k) + d_u_con(i, k)
          v_seri(i, k) = v_seri(i, k) + d_v_con(i, k)
       ENDDO
    ENDDO

    IF (if_ebil >= 2) THEN 
       ztit='after convect'
       CALL diagetpq(airephy, ztit, ip_ebil, 2, 2, pdtphys &
            , t_seri, q_seri, ql_seri, qs_seri, u_seri, v_seri, paprs &
            , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
       call diagphy(airephy, ztit, ip_ebil &
            , zero_v, zero_v, zero_v, zero_v, zero_v &
            , zero_v, rain_con, snow_con, ztsol &
            , d_h_vcol, d_qt, d_ec &
            , fs_bound, fq_bound )
    END IF

    IF (check) THEN
       za = qcheck(klon, llm, paprs, q_seri, ql_seri, airephy)
       print *,"aprescon=", za
       zx_t = 0.0
       za = 0.0
       DO i = 1, klon
          za = za + airephy(i)/REAL(klon)
          zx_t = zx_t + (rain_con(i)+ &
               snow_con(i))*airephy(i)/REAL(klon)
       ENDDO
       zx_t = zx_t/za*pdtphys
       print *,"Precip=", zx_t
    ENDIF
    IF (zx_ajustq) THEN
       DO i = 1, klon
          z_apres(i) = 0.0
       ENDDO
       DO k = 1, llm
          DO i = 1, klon
             z_apres(i) = z_apres(i) + (q_seri(i, k)+ql_seri(i, k)) &
                  *zmasse(i, k)
          ENDDO
       ENDDO
       DO i = 1, klon
          z_factor(i) = (z_avant(i)-(rain_con(i)+snow_con(i))*pdtphys) &
               /z_apres(i)
       ENDDO
       DO k = 1, llm
          DO i = 1, klon
             IF (z_factor(i).GT.(1.0+1.0E-08) .OR. &
                  z_factor(i) < (1.0-1.0E-08)) THEN
                q_seri(i, k) = q_seri(i, k) * z_factor(i)
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    zx_ajustq=.FALSE.

    ! Convection seche (thermiques ou ajustement)

    d_t_ajs=0.
    d_u_ajs=0.
    d_v_ajs=0.
    d_q_ajs=0.
    fm_therm=0.
    entr_therm=0.

    IF(prt_level>9)print *, &
         'AVANT LA CONVECTION SECHE, iflag_thermals=' &
         , iflag_thermals, '   nsplit_thermals=', nsplit_thermals
    if(iflag_thermals < 0) then
       !  Rien
       IF(prt_level>9)print *,'pas de convection'
    else if(iflag_thermals == 0) then
       !  Ajustement sec
       IF(prt_level>9)print *,'ajsec'
       CALL ajsec(paprs, pplay, t_seri, q_seri, d_t_ajs, d_q_ajs)
       t_seri = t_seri + d_t_ajs
       q_seri = q_seri + d_q_ajs
    else
       !  Thermiques
       IF(prt_level>9)print *,'JUSTE AVANT, iflag_thermals=' &
            , iflag_thermals, '   nsplit_thermals=', nsplit_thermals
       call calltherm(pdtphys &
            , pplay, paprs, pphi &
            , u_seri, v_seri, t_seri, q_seri &
            , d_u_ajs, d_v_ajs, d_t_ajs, d_q_ajs &
            , fm_therm, entr_therm)
    endif

    IF (if_ebil >= 2) THEN 
       ztit='after dry_adjust'
       CALL diagetpq(airephy, ztit, ip_ebil, 2, 2, pdtphys &
            , t_seri, q_seri, ql_seri, qs_seri, u_seri, v_seri, paprs &
            , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
    END IF

    !  Caclul des ratqs

    !   ratqs convectifs a l'ancienne en fonction de q(z=0)-q / q
    !   on ecrase le tableau ratqsc calcule par clouds_gno
    if (iflag_cldcon == 1) then
       do k=1, llm
          do i=1, klon
             if(ptconv(i, k)) then
                ratqsc(i, k)=ratqsbas &
                     +fact_cldcon*(q_seri(i, 1)-q_seri(i, k))/q_seri(i, k)
             else
                ratqsc(i, k)=0.
             endif
          enddo
       enddo
    endif

    !   ratqs stables
    do k=1, llm
       do i=1, klon
          ratqss(i, k)=ratqsbas+(ratqshaut-ratqsbas)* &
               min((paprs(i, 1)-pplay(i, k))/(paprs(i, 1)-30000.), 1.) 
       enddo
    enddo

    !  ratqs final
    if (iflag_cldcon == 1 .or.iflag_cldcon == 2) then
       !   les ratqs sont une conbinaison de ratqss et ratqsc
       !   ratqs final
       !   1e4 (en gros 3 heures), en dur pour le moment, est le temps de
       !   relaxation des ratqs
       facteur=exp(-pdtphys*facttemps)
       ratqs=max(ratqs*facteur, ratqss)
       ratqs=max(ratqs, ratqsc)
    else
       !   on ne prend que le ratqs stable pour fisrtilp
       ratqs=ratqss
    endif

    ! Appeler le processus de condensation a grande echelle
    ! et le processus de precipitation
    CALL fisrtilp(pdtphys, paprs, pplay, &
         t_seri, q_seri, ptconv, ratqs, &
         d_t_lsc, d_q_lsc, d_ql_lsc, rneb, cldliq, &
         rain_lsc, snow_lsc, &
         pfrac_impa, pfrac_nucl, pfrac_1nucl, &
         frac_impa, frac_nucl, &
         prfl, psfl, rhcl)

    WHERE (rain_lsc < 0) rain_lsc = 0.
    WHERE (snow_lsc < 0) snow_lsc = 0.
    DO k = 1, llm
       DO i = 1, klon
          t_seri(i, k) = t_seri(i, k) + d_t_lsc(i, k)
          q_seri(i, k) = q_seri(i, k) + d_q_lsc(i, k)
          ql_seri(i, k) = ql_seri(i, k) + d_ql_lsc(i, k)
          cldfra(i, k) = rneb(i, k)
          IF (.NOT.new_oliq) cldliq(i, k) = ql_seri(i, k)
       ENDDO
    ENDDO
    IF (check) THEN
       za = qcheck(klon, llm, paprs, q_seri, ql_seri, airephy)
       print *,"apresilp=", za
       zx_t = 0.0
       za = 0.0
       DO i = 1, klon
          za = za + airephy(i)/REAL(klon)
          zx_t = zx_t + (rain_lsc(i) &
               + snow_lsc(i))*airephy(i)/REAL(klon)
       ENDDO
       zx_t = zx_t/za*pdtphys
       print *,"Precip=", zx_t
    ENDIF

    IF (if_ebil >= 2) THEN 
       ztit='after fisrt'
       CALL diagetpq(airephy, ztit, ip_ebil, 2, 2, pdtphys &
            , t_seri, q_seri, ql_seri, qs_seri, u_seri, v_seri, paprs &
            , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
       call diagphy(airephy, ztit, ip_ebil &
            , zero_v, zero_v, zero_v, zero_v, zero_v &
            , zero_v, rain_lsc, snow_lsc, ztsol &
            , d_h_vcol, d_qt, d_ec &
            , fs_bound, fq_bound )
    END IF

    !  PRESCRIPTION DES NUAGES POUR LE RAYONNEMENT

    ! 1. NUAGES CONVECTIFS

    IF (iflag_cldcon.le.-1) THEN ! seulement pour Tiedtke
       snow_tiedtke=0.
       if (iflag_cldcon == -1) then
          rain_tiedtke=rain_con
       else
          rain_tiedtke=0.
          do k=1, llm
             do i=1, klon
                if (d_q_con(i, k) < 0.) then
                   rain_tiedtke(i)=rain_tiedtke(i)-d_q_con(i, k)/pdtphys &
                        *zmasse(i, k)
                endif
             enddo
          enddo
       endif

       ! Nuages diagnostiques pour Tiedtke
       CALL diagcld1(paprs, pplay, &
            rain_tiedtke, snow_tiedtke, ibas_con, itop_con, &
            diafra, dialiq)
       DO k = 1, llm
          DO i = 1, klon
             IF (diafra(i, k).GT.cldfra(i, k)) THEN
                cldliq(i, k) = dialiq(i, k)
                cldfra(i, k) = diafra(i, k)
             ENDIF
          ENDDO
       ENDDO

    ELSE IF (iflag_cldcon == 3) THEN
       ! On prend pour les nuages convectifs le max du calcul de la
       ! convection et du calcul du pas de temps précédent diminué d'un facteur
       ! facttemps
       facteur = pdtphys *facttemps
       do k=1, llm
          do i=1, klon
             rnebcon(i, k)=rnebcon(i, k)*facteur
             if (rnebcon0(i, k)*clwcon0(i, k).gt.rnebcon(i, k)*clwcon(i, k)) &
                  then
                rnebcon(i, k)=rnebcon0(i, k)
                clwcon(i, k)=clwcon0(i, k)
             endif
          enddo
       enddo

       !   On prend la somme des fractions nuageuses et des contenus en eau
       cldfra=min(max(cldfra, rnebcon), 1.)
       cldliq=cldliq+rnebcon*clwcon

    ENDIF

    ! 2. NUAGES STARTIFORMES

    IF (ok_stratus) THEN
       CALL diagcld2(paprs, pplay, t_seri, q_seri, diafra, dialiq)
       DO k = 1, llm
          DO i = 1, klon
             IF (diafra(i, k).GT.cldfra(i, k)) THEN
                cldliq(i, k) = dialiq(i, k)
                cldfra(i, k) = diafra(i, k)
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    ! Precipitation totale

    DO i = 1, klon
       rain_fall(i) = rain_con(i) + rain_lsc(i)
       snow_fall(i) = snow_con(i) + snow_lsc(i)
    ENDDO

    IF (if_ebil >= 2) THEN 
       ztit="after diagcld"
       CALL diagetpq(airephy, ztit, ip_ebil, 2, 2, pdtphys &
            , t_seri, q_seri, ql_seri, qs_seri, u_seri, v_seri, paprs &
            , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
    END IF

    ! Calculer l'humidite relative pour diagnostique

    DO k = 1, llm
       DO i = 1, klon
          zx_t = t_seri(i, k)
          IF (thermcep) THEN
             zdelta = MAX(0., SIGN(1., rtt-zx_t))
             zx_qs  = r2es * FOEEW(zx_t, zdelta)/pplay(i, k)
             zx_qs  = MIN(0.5, zx_qs)
             zcor   = 1./(1.-retv*zx_qs)
             zx_qs  = zx_qs*zcor
          ELSE
             IF (zx_t < t_coup) THEN
                zx_qs = qsats(zx_t)/pplay(i, k)
             ELSE
                zx_qs = qsatl(zx_t)/pplay(i, k)
             ENDIF
          ENDIF
          zx_rh(i, k) = q_seri(i, k)/zx_qs
          zqsat(i, k)=zx_qs
       ENDDO
    ENDDO
    !jq - introduce the aerosol direct and first indirect radiative forcings
    !jq - Johannes Quaas, 27/11/2003 (quaas@lmd.jussieu.fr)
    IF (ok_ade.OR.ok_aie) THEN
       ! Get sulfate aerosol distribution
       CALL readsulfate(rdayvrai, firstcal, sulfate)
       CALL readsulfate_preind(rdayvrai, firstcal, sulfate_pi)

       ! Calculate aerosol optical properties (Olivier Boucher)
       CALL aeropt(pplay, paprs, t_seri, sulfate, rhcl, &
            tau_ae, piz_ae, cg_ae, aerindex)
    ELSE
       tau_ae(:, :, :)=0.0
       piz_ae(:, :, :)=0.0
       cg_ae(:, :, :)=0.0
    ENDIF

    ! Calculer les parametres optiques des nuages et quelques
    ! parametres pour diagnostiques:

    if (ok_newmicro) then
       CALL newmicro (paprs, pplay, ok_newmicro, &
            t_seri, cldliq, cldfra, cldtau, cldemi, &
            cldh, cldl, cldm, cldt, cldq, &
            flwp, fiwp, flwc, fiwc, &
            ok_aie, &
            sulfate, sulfate_pi, &
            bl95_b0, bl95_b1, &
            cldtaupi, re, fl)
    else
       CALL nuage (paprs, pplay, &
            t_seri, cldliq, cldfra, cldtau, cldemi, &
            cldh, cldl, cldm, cldt, cldq, &
            ok_aie, &
            sulfate, sulfate_pi, &
            bl95_b0, bl95_b1, &
            cldtaupi, re, fl)

    endif

    ! Appeler le rayonnement mais calculer tout d'abord l'albedo du sol.

    IF (MOD(itaprad, radpas) == 0) THEN
       DO i = 1, klon
          albsol(i) = falbe(i, is_oce) * pctsrf(i, is_oce) &
               + falbe(i, is_lic) * pctsrf(i, is_lic) &
               + falbe(i, is_ter) * pctsrf(i, is_ter) &
               + falbe(i, is_sic) * pctsrf(i, is_sic)
          albsollw(i) = falblw(i, is_oce) * pctsrf(i, is_oce) &
               + falblw(i, is_lic) * pctsrf(i, is_lic) &
               + falblw(i, is_ter) * pctsrf(i, is_ter) &
               + falblw(i, is_sic) * pctsrf(i, is_sic)
       ENDDO
       ! nouveau rayonnement (compatible Arpege-IFS):
       CALL radlwsw(dist, rmu0, fract,  &
            paprs, pplay, zxtsol, albsol, albsollw, t_seri, q_seri, &
            wo, &
            cldfra, cldemi, cldtau, &
            heat, heat0, cool, cool0, radsol, albpla, &
            topsw, toplw, solsw, sollw, &
            sollwdown, &
            topsw0, toplw0, solsw0, sollw0, &
            lwdn0, lwdn, lwup0, lwup,  &
            swdn0, swdn, swup0, swup, &
            ok_ade, ok_aie, & ! new for aerosol radiative effects
            tau_ae, piz_ae, cg_ae, &
            topswad, solswad, &
            cldtaupi, &
            topswai, solswai)
       itaprad = 0
    ENDIF
    itaprad = itaprad + 1

    ! Ajouter la tendance des rayonnements (tous les pas)

    DO k = 1, llm
       DO i = 1, klon
          t_seri(i, k) = t_seri(i, k) &
               + (heat(i, k)-cool(i, k)) * pdtphys/86400.
       ENDDO
    ENDDO

    IF (if_ebil >= 2) THEN 
       ztit='after rad'
       CALL diagetpq(airephy, ztit, ip_ebil, 2, 2, pdtphys &
            , t_seri, q_seri, ql_seri, qs_seri, u_seri, v_seri, paprs &
            , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
       call diagphy(airephy, ztit, ip_ebil &
            , topsw, toplw, solsw, sollw, zero_v &
            , zero_v, zero_v, zero_v, ztsol &
            , d_h_vcol, d_qt, d_ec &
            , fs_bound, fq_bound )
    END IF

    ! Calculer l'hydrologie de la surface

    DO i = 1, klon
       zxqsurf(i) = 0.0
       zxsnow(i) = 0.0
    ENDDO
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          zxqsurf(i) = zxqsurf(i) + fqsurf(i, nsrf)*pctsrf(i, nsrf)
          zxsnow(i) = zxsnow(i) + fsnow(i, nsrf)*pctsrf(i, nsrf)
       ENDDO
    ENDDO

    ! Calculer le bilan du sol et la derive de temperature (couplage)

    DO i = 1, klon
       bils(i) = radsol(i) - sens(i) + zxfluxlat(i)
    ENDDO

    !mod deb lott(jan95)
    ! Appeler le programme de parametrisation de l'orographie
    ! a l'echelle sous-maille:

    IF (ok_orodr) THEN
       !  selection des points pour lesquels le shema est actif:
       igwd=0
       DO i=1, klon
          itest(i)=0
          IF (((zpic(i)-zmea(i)).GT.100.).AND.(zstd(i).GT.10.0)) THEN
             itest(i)=1
             igwd=igwd+1
             idx(igwd)=i
          ENDIF
       ENDDO

       CALL drag_noro(klon, llm, pdtphys, paprs, pplay, &
            zmea, zstd, zsig, zgam, zthe, zpic, zval, &
            igwd, idx, itest, &
            t_seri, u_seri, v_seri, &
            zulow, zvlow, zustrdr, zvstrdr, &
            d_t_oro, d_u_oro, d_v_oro)

       !  ajout des tendances
       DO k = 1, llm
          DO i = 1, klon
             t_seri(i, k) = t_seri(i, k) + d_t_oro(i, k)
             u_seri(i, k) = u_seri(i, k) + d_u_oro(i, k)
             v_seri(i, k) = v_seri(i, k) + d_v_oro(i, k)
          ENDDO
       ENDDO
    ENDIF

    IF (ok_orolf) THEN

       !  selection des points pour lesquels le shema est actif:
       igwd=0
       DO i=1, klon
          itest(i)=0
          IF ((zpic(i)-zmea(i)).GT.100.) THEN
             itest(i)=1
             igwd=igwd+1
             idx(igwd)=i
          ENDIF
       ENDDO

       CALL lift_noro(klon, llm, pdtphys, paprs, pplay, &
            rlat, zmea, zstd, zpic, &
            itest, &
            t_seri, u_seri, v_seri, &
            zulow, zvlow, zustrli, zvstrli, &
            d_t_lif, d_u_lif, d_v_lif)

       !  ajout des tendances
       DO k = 1, llm
          DO i = 1, klon
             t_seri(i, k) = t_seri(i, k) + d_t_lif(i, k)
             u_seri(i, k) = u_seri(i, k) + d_u_lif(i, k)
             v_seri(i, k) = v_seri(i, k) + d_v_lif(i, k)
          ENDDO
       ENDDO

    ENDIF ! fin de test sur ok_orolf

    ! STRESS NECESSAIRES: TOUTE LA PHYSIQUE

    DO i = 1, klon
       zustrph(i)=0.
       zvstrph(i)=0.
    ENDDO
    DO k = 1, llm
       DO i = 1, klon
          zustrph(i)=zustrph(i)+(u_seri(i, k)-u(i, k))/pdtphys* zmasse(i, k)
          zvstrph(i)=zvstrph(i)+(v_seri(i, k)-v(i, k))/pdtphys* zmasse(i, k)
       ENDDO
    ENDDO

    !IM calcul composantes axiales du moment angulaire et couple des montagnes

    CALL aaam_bud(27, klon, llm, gmtime, &
         ra, rg, romega, &
         rlat, rlon, pphis, &
         zustrdr, zustrli, zustrph, &
         zvstrdr, zvstrli, zvstrph, &
         paprs, u, v, &
         aam, torsfc)

    IF (if_ebil >= 2) THEN 
       ztit='after orography'
       CALL diagetpq(airephy, ztit, ip_ebil, 2, 2, pdtphys &
            , t_seri, q_seri, ql_seri, qs_seri, u_seri, v_seri, paprs &
            , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
    END IF

    !AA Installation de l'interface online-offline pour traceurs

    !   Calcul  des tendances traceurs

    call phytrac(rnpb, itap, lmt_pas, julien,  gmtime, firstcal, lafin, nq-2, &
         pdtphys, u, v, t, paprs, pplay, pmfu,  pmfd,  pen_u,  pde_u,  pen_d, &
         pde_d, ycoefh, fm_therm, entr_therm, yu1, yv1, ftsol, pctsrf, &
         frac_impa,  frac_nucl, presnivs, pphis, pphi, albsol, rhcl, cldfra, &
         rneb,  diafra,  cldliq, itop_con, ibas_con, pmflxr, pmflxs, prfl, &
         psfl, da, phi, mp, upwd, dnwd, tr_seri, zmasse)

    IF (offline) THEN

       print*, 'Attention on met a 0 les thermiques pour phystoke'
       call phystokenc(pdtphys, rlon, rlat, &
            t, pmfu, pmfd, pen_u, pde_u, pen_d, pde_d, &
            fm_therm, entr_therm, &
            ycoefh, yu1, yv1, ftsol, pctsrf, &
            frac_impa, frac_nucl, &
            pphis, airephy, pdtphys, itap)

    ENDIF

    ! Calculer le transport de l'eau et de l'energie (diagnostique)

    CALL transp (paprs, zxtsol, &
         t_seri, q_seri, u_seri, v_seri, zphi, &
         ve, vq, ue, uq)

    !IM diag. bilKP

    CALL transp_lay (paprs, zxtsol, &
         t_seri, q_seri, u_seri, v_seri, zphi, &
         ve_lay, vq_lay, ue_lay, uq_lay)

    ! Accumuler les variables a stocker dans les fichiers histoire:

    !+jld ec_conser
    DO k = 1, llm
       DO i = 1, klon
          ZRCPD = RCPD*(1.0+RVTMP2*q_seri(i, k))
          d_t_ec(i, k)=0.5/ZRCPD &
               *(u(i, k)**2+v(i, k)**2-u_seri(i, k)**2-v_seri(i, k)**2)
          t_seri(i, k)=t_seri(i, k)+d_t_ec(i, k)
          d_t_ec(i, k) = d_t_ec(i, k)/pdtphys
       END DO
    END DO
    !-jld ec_conser
    IF (if_ebil >= 1) THEN 
       ztit='after physic'
       CALL diagetpq(airephy, ztit, ip_ebil, 1, 1, pdtphys &
            , t_seri, q_seri, ql_seri, qs_seri, u_seri, v_seri, paprs &
            , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
       !     Comme les tendances de la physique sont ajoute dans la dynamique, 
       !     on devrait avoir que la variation d'entalpie par la dynamique
       !     est egale a la variation de la physique au pas de temps precedent.
       !     Donc la somme de ces 2 variations devrait etre nulle.
       call diagphy(airephy, ztit, ip_ebil &
            , topsw, toplw, solsw, sollw, sens &
            , evap, rain_fall, snow_fall, ztsol &
            , d_h_vcol, d_qt, d_ec &
            , fs_bound, fq_bound )

       d_h_vcol_phy=d_h_vcol

    END IF

    !   SORTIES

    !cc prw = eau precipitable
    DO i = 1, klon
       prw(i) = 0.
       DO k = 1, llm
          prw(i) = prw(i) + q_seri(i, k)*zmasse(i, k)
       ENDDO
    ENDDO

    ! Convertir les incrementations en tendances

    DO k = 1, llm
       DO i = 1, klon
          d_u(i, k) = ( u_seri(i, k) - u(i, k) ) / pdtphys
          d_v(i, k) = ( v_seri(i, k) - v(i, k) ) / pdtphys
          d_t(i, k) = ( t_seri(i, k)-t(i, k) ) / pdtphys
          d_qx(i, k, ivap) = ( q_seri(i, k) - qx(i, k, ivap) ) / pdtphys
          d_qx(i, k, iliq) = ( ql_seri(i, k) - qx(i, k, iliq) ) / pdtphys
       ENDDO
    ENDDO

    IF (nq >= 3) THEN
       DO iq = 3, nq
          DO  k = 1, llm
             DO  i = 1, klon
                d_qx(i, k, iq) = (tr_seri(i, k, iq-2) - qx(i, k, iq)) / pdtphys
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    ! Sauvegarder les valeurs de t et q a la fin de la physique:
    DO k = 1, llm
       DO i = 1, klon
          t_ancien(i, k) = t_seri(i, k)
          q_ancien(i, k) = q_seri(i, k)
       ENDDO
    ENDDO

    !   Ecriture des sorties
    call write_histhf
    call write_histday
    call write_histins

    ! Si c'est la fin, il faut conserver l'etat de redemarrage
    IF (lafin) THEN
       itau_phy = itau_phy + itap
       CALL phyredem("restartphy.nc", rlat, rlon, pctsrf, ftsol, &
            ftsoil, tslab, seaice, fqsurf, qsol, &
            fsnow, falbe, falblw, fevap, rain_fall, snow_fall, &
            solsw, sollwdown, dlw, &
            radsol, frugs, agesno, &
            zmea, zstd, zsig, zgam, zthe, zpic, zval, &
            t_ancien, q_ancien, rnebcon, ratqs, clwcon, run_off_lic_0)
    ENDIF

  contains

    subroutine write_histday

      use grid_change, only: gr_phy_write_3d
      integer itau_w  ! pas de temps ecriture

      !------------------------------------------------

      if (ok_journe) THEN
         itau_w = itau_phy + itap
         if (nq <= 4) then
            call histwrite(nid_day, "Sigma_O3_Royer", itau_w, &
                 gr_phy_write_3d(wo) * 1e3)
            ! (convert "wo" from kDU to DU)
         end if
         if (ok_sync) then
            call histsync(nid_day)
         endif
      ENDIF

    End subroutine write_histday

    !****************************

    subroutine write_histhf

      ! From phylmd/write_histhf.h, v 1.5 2005/05/25 13:10:09

      !------------------------------------------------

      call write_histhf3d

      IF (ok_sync) THEN
         call histsync(nid_hf)
      ENDIF

    end subroutine write_histhf

    !***************************************************************

    subroutine write_histins

      ! From phylmd/write_histins.h, v 1.2 2005/05/25 13:10:09

      real zout
      integer itau_w  ! pas de temps ecriture

      !--------------------------------------------------

      IF (ok_instan) THEN
         ! Champs 2D:

         zsto = pdtphys * ecrit_ins
         zout = pdtphys * ecrit_ins
         itau_w = itau_phy + itap

         i = NINT(zout/zsto)
         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), pphis, zx_tmp_2d)
         CALL histwrite(nid_ins, "phis", itau_w, zx_tmp_2d)

         i = NINT(zout/zsto)
         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), airephy, zx_tmp_2d)
         CALL histwrite(nid_ins, "aire", itau_w, zx_tmp_2d)

         DO i = 1, klon
            zx_tmp_fi2d(i) = paprs(i, 1)
         ENDDO
         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zx_tmp_fi2d, zx_tmp_2d)
         CALL histwrite(nid_ins, "psol", itau_w, zx_tmp_2d)

         DO i = 1, klon
            zx_tmp_fi2d(i) = rain_fall(i) + snow_fall(i)
         ENDDO
         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zx_tmp_fi2d, zx_tmp_2d)
         CALL histwrite(nid_ins, "precip", itau_w, zx_tmp_2d)

         DO i = 1, klon
            zx_tmp_fi2d(i) = rain_lsc(i) + snow_lsc(i)
         ENDDO
         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zx_tmp_fi2d, zx_tmp_2d)
         CALL histwrite(nid_ins, "plul", itau_w, zx_tmp_2d)

         DO i = 1, klon
            zx_tmp_fi2d(i) = rain_con(i) + snow_con(i)
         ENDDO
         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zx_tmp_fi2d, zx_tmp_2d)
         CALL histwrite(nid_ins, "pluc", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zxtsol, zx_tmp_2d)
         CALL histwrite(nid_ins, "tsol", itau_w, zx_tmp_2d)
         !ccIM
         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zt2m, zx_tmp_2d)
         CALL histwrite(nid_ins, "t2m", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zq2m, zx_tmp_2d)
         CALL histwrite(nid_ins, "q2m", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zu10m, zx_tmp_2d)
         CALL histwrite(nid_ins, "u10m", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zv10m, zx_tmp_2d)
         CALL histwrite(nid_ins, "v10m", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), snow_fall, zx_tmp_2d)
         CALL histwrite(nid_ins, "snow", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), cdragm, zx_tmp_2d)
         CALL histwrite(nid_ins, "cdrm", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), cdragh, zx_tmp_2d)
         CALL histwrite(nid_ins, "cdrh", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), toplw, zx_tmp_2d)
         CALL histwrite(nid_ins, "topl", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), evap, zx_tmp_2d)
         CALL histwrite(nid_ins, "evap", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), solsw, zx_tmp_2d)
         CALL histwrite(nid_ins, "sols", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), sollw, zx_tmp_2d)
         CALL histwrite(nid_ins, "soll", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), sollwdown, zx_tmp_2d)
         CALL histwrite(nid_ins, "solldown", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), bils, zx_tmp_2d)
         CALL histwrite(nid_ins, "bils", itau_w, zx_tmp_2d)

         zx_tmp_fi2d(1:klon)=-1*sens(1:klon)
         !     CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), sens, zx_tmp_2d)
         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zx_tmp_fi2d, zx_tmp_2d)
         CALL histwrite(nid_ins, "sens", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), fder, zx_tmp_2d)
         CALL histwrite(nid_ins, "fder", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), d_ts(1, is_oce), zx_tmp_2d)
         CALL histwrite(nid_ins, "dtsvdfo", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), d_ts(1, is_ter), zx_tmp_2d)
         CALL histwrite(nid_ins, "dtsvdft", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), d_ts(1, is_lic), zx_tmp_2d)
         CALL histwrite(nid_ins, "dtsvdfg", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), d_ts(1, is_sic), zx_tmp_2d)
         CALL histwrite(nid_ins, "dtsvdfi", itau_w, zx_tmp_2d)

         DO nsrf = 1, nbsrf
            !XXX
            zx_tmp_fi2d(1 : klon) = pctsrf( 1 : klon, nsrf)*100.
            CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zx_tmp_fi2d, zx_tmp_2d)
            CALL histwrite(nid_ins, "pourc_"//clnsurf(nsrf), itau_w, &
                 zx_tmp_2d) 

            zx_tmp_fi2d(1 : klon) = pctsrf( 1 : klon, nsrf)
            CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zx_tmp_fi2d, zx_tmp_2d)
            CALL histwrite(nid_ins, "fract_"//clnsurf(nsrf), itau_w, &
                 zx_tmp_2d) 

            zx_tmp_fi2d(1 : klon) = fluxt( 1 : klon, 1, nsrf)
            CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zx_tmp_fi2d, zx_tmp_2d)
            CALL histwrite(nid_ins, "sens_"//clnsurf(nsrf), itau_w, &
                 zx_tmp_2d) 

            zx_tmp_fi2d(1 : klon) = fluxlat( 1 : klon, nsrf)
            CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zx_tmp_fi2d, zx_tmp_2d)
            CALL histwrite(nid_ins, "lat_"//clnsurf(nsrf), itau_w, &
                 zx_tmp_2d) 

            zx_tmp_fi2d(1 : klon) = ftsol( 1 : klon, nsrf)
            CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zx_tmp_fi2d, zx_tmp_2d)
            CALL histwrite(nid_ins, "tsol_"//clnsurf(nsrf), itau_w, &
                 zx_tmp_2d) 

            zx_tmp_fi2d(1 : klon) = fluxu( 1 : klon, 1, nsrf)
            CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zx_tmp_fi2d, zx_tmp_2d)
            CALL histwrite(nid_ins, "taux_"//clnsurf(nsrf), itau_w, &
                 zx_tmp_2d) 

            zx_tmp_fi2d(1 : klon) = fluxv( 1 : klon, 1, nsrf)
            CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zx_tmp_fi2d, zx_tmp_2d)
            CALL histwrite(nid_ins, "tauy_"//clnsurf(nsrf), itau_w, &
                 zx_tmp_2d)

            zx_tmp_fi2d(1 : klon) = frugs( 1 : klon, nsrf)
            CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zx_tmp_fi2d, zx_tmp_2d)
            CALL histwrite(nid_ins, "rugs_"//clnsurf(nsrf), itau_w, &
                 zx_tmp_2d) 

            zx_tmp_fi2d(1 : klon) = falbe( 1 : klon, nsrf)
            CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zx_tmp_fi2d, zx_tmp_2d)
            CALL histwrite(nid_ins, "albe_"//clnsurf(nsrf), itau_w, &
                 zx_tmp_2d) 

         END DO
         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), albsol, zx_tmp_2d)
         CALL histwrite(nid_ins, "albs", itau_w, zx_tmp_2d)
         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), albsollw, zx_tmp_2d)
         CALL histwrite(nid_ins, "albslw", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), zxrugs, zx_tmp_2d)
         CALL histwrite(nid_ins, "rugs", itau_w, zx_tmp_2d)

         !IM cf. AM 081204 BEG

         !HBTM2

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), s_pblh, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_pblh", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), s_pblt, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_pblt", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), s_lcl, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_lcl", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), s_capCL, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_capCL", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), s_oliqCL, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_oliqCL", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), s_cteiCL, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_cteiCL", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), s_therm, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_therm", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), s_trmb1, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_trmb1", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), s_trmb2, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_trmb2", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), s_trmb3, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_trmb3", itau_w, zx_tmp_2d)

         !IM cf. AM 081204 END

         ! Champs 3D:

         CALL gr_fi_ecrit(llm, klon, iim, (jjm + 1), t_seri, zx_tmp_3d)
         CALL histwrite(nid_ins, "temp", itau_w, zx_tmp_3d)

         CALL gr_fi_ecrit(llm, klon, iim, (jjm + 1), u_seri, zx_tmp_3d)
         CALL histwrite(nid_ins, "vitu", itau_w, zx_tmp_3d)

         CALL gr_fi_ecrit(llm, klon, iim, (jjm + 1), v_seri, zx_tmp_3d)
         CALL histwrite(nid_ins, "vitv", itau_w, zx_tmp_3d)

         CALL gr_fi_ecrit(llm, klon, iim, (jjm + 1), zphi, zx_tmp_3d)
         CALL histwrite(nid_ins, "geop", itau_w, zx_tmp_3d)

         CALL gr_fi_ecrit(llm, klon, iim, (jjm + 1), pplay, zx_tmp_3d)
         CALL histwrite(nid_ins, "pres", itau_w, zx_tmp_3d)

         CALL gr_fi_ecrit(llm, klon, iim, (jjm + 1), d_t_vdf, zx_tmp_3d)
         CALL histwrite(nid_ins, "dtvdf", itau_w, zx_tmp_3d)

         CALL gr_fi_ecrit(llm, klon, iim, (jjm + 1), d_q_vdf, zx_tmp_3d)
         CALL histwrite(nid_ins, "dqvdf", itau_w, zx_tmp_3d)

         if (ok_sync) then
            call histsync(nid_ins)
         endif
      ENDIF

    end subroutine write_histins

    !****************************************************

    subroutine write_histhf3d

      ! From phylmd/write_histhf3d.h, v 1.2 2005/05/25 13:10:09

      integer itau_w  ! pas de temps ecriture

      !-------------------------------------------------------

      itau_w = itau_phy + itap

      ! Champs 3D:

      CALL gr_fi_ecrit(llm, klon, iim, (jjm + 1), t_seri, zx_tmp_3d)
      CALL histwrite(nid_hf3d, "temp", itau_w, zx_tmp_3d)

      CALL gr_fi_ecrit(llm, klon, iim, (jjm + 1), qx(1, 1, ivap), zx_tmp_3d)
      CALL histwrite(nid_hf3d, "ovap", itau_w, zx_tmp_3d)

      CALL gr_fi_ecrit(llm, klon, iim, (jjm + 1), u_seri, zx_tmp_3d)
      CALL histwrite(nid_hf3d, "vitu", itau_w, zx_tmp_3d)

      CALL gr_fi_ecrit(llm, klon, iim, (jjm + 1), v_seri, zx_tmp_3d)
      CALL histwrite(nid_hf3d, "vitv", itau_w, zx_tmp_3d)

      if (nbtr >= 3) then
         CALL gr_fi_ecrit(llm, klon, iim, (jjm + 1), tr_seri(1, 1, 3), &
              zx_tmp_3d)
         CALL histwrite(nid_hf3d, "O3", itau_w, zx_tmp_3d)
      end if

      if (ok_sync) then
         call histsync(nid_hf3d)
      endif

    end subroutine write_histhf3d

  END SUBROUTINE physiq

end module physiq_m
