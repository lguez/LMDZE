module physiq_m

  IMPLICIT none

contains

  SUBROUTINE physiq(lafin, dayvrai, time, paprs, play, pphi, pphis, u, v, t, &
       qx, omega, d_u, d_v, d_t, d_qx)

    ! From phylmd/physiq.F, version 1.22 2006/02/20 09:38:28
    ! (subversion revision 678)

    ! Author: Z. X. Li (LMD/CNRS) 1993

    ! This is the main procedure for the "physics" part of the program.

    use aaam_bud_m, only: aaam_bud
    USE abort_gcm_m, ONLY: abort_gcm
    use aeropt_m, only: aeropt
    use ajsec_m, only: ajsec
    use calltherm_m, only: calltherm
    USE clesphys, ONLY: cdhmax, cdmmax, co2_ppm, ecrit_hf, ecrit_ins, &
         ecrit_mth, ecrit_reg, ecrit_tra, ksta, ksta_ter, ok_kzmin
    USE clesphys2, ONLY: cycle_diurne, iflag_con, nbapp_rad, new_oliq, &
         ok_orodr, ok_orolf
    USE clmain_m, ONLY: clmain
    use clouds_gno_m, only: clouds_gno
    use comconst, only: dtphys
    USE comgeomphy, ONLY: airephy
    USE concvl_m, ONLY: concvl
    USE conf_gcm_m, ONLY: offline, raz_date, day_step, iphysiq
    USE conf_phys_m, ONLY: conf_phys
    use conflx_m, only: conflx
    USE ctherm, ONLY: iflag_thermals, nsplit_thermals
    use diagcld2_m, only: diagcld2
    use diagetpq_m, only: diagetpq
    use diagphy_m, only: diagphy
    USE dimens_m, ONLY: llm, nqmx
    USE dimphy, ONLY: klon
    USE dimsoil, ONLY: nsoilmx
    use drag_noro_m, only: drag_noro
    use dynetat0_m, only: day_ref, annee_ref
    USE fcttre, ONLY: foeew, qsatl, qsats, thermcep
    use fisrtilp_m, only: fisrtilp
    USE hgardfou_m, ONLY: hgardfou
    USE indicesol, ONLY: clnsurf, epsfra, is_lic, is_oce, is_sic, is_ter, &
         nbsrf
    USE ini_histins_m, ONLY: ini_histins
    use netcdf95, only: NF95_CLOSE
    use newmicro_m, only: newmicro
    USE orbite_m, ONLY: orbite
    USE ozonecm_m, ONLY: ozonecm
    USE phyetat0_m, ONLY: phyetat0, rlat, rlon
    USE phyredem_m, ONLY: phyredem
    USE phyredem0_m, ONLY: phyredem0
    USE phystokenc_m, ONLY: phystokenc
    USE phytrac_m, ONLY: phytrac
    USE qcheck_m, ONLY: qcheck
    use radlwsw_m, only: radlwsw
    use readsulfate_m, only: readsulfate
    use readsulfate_preind_m, only: readsulfate_preind
    use yoegwd, only: sugwd
    USE suphec_m, ONLY: ra, rcpd, retv, rg, rlvtt, romega, rsigma, rtt
    USE temps, ONLY: itau_phy
    use unit_nml_m, only: unit_nml
    USE ymds2ju_m, ONLY: ymds2ju
    USE yoethf_m, ONLY: r2es, rvtmp2
    use zenang_m, only: zenang

    logical, intent(in):: lafin ! dernier passage

    integer, intent(in):: dayvrai
    ! current day number, based at value 1 on January 1st of annee_ref

    REAL, intent(in):: time ! heure de la journ\'ee en fraction de jour

    REAL, intent(in):: paprs(:, :) ! (klon, llm + 1)
    ! pression pour chaque inter-couche, en Pa

    REAL, intent(in):: play(:, :) ! (klon, llm)
    ! pression pour le mileu de chaque couche (en Pa)

    REAL, intent(in):: pphi(:, :) ! (klon, llm) 
    ! géopotentiel de chaque couche (référence sol)

    REAL, intent(in):: pphis(:) ! (klon) géopotentiel du sol

    REAL, intent(in):: u(:, :) ! (klon, llm)
    ! vitesse dans la direction X (de O a E) en m/s

    REAL, intent(in):: v(:, :) ! (klon, llm) vitesse Y (de S a N) en m/s
    REAL, intent(in):: t(:, :) ! (klon, llm) temperature (K)

    REAL, intent(in):: qx(:, :, :) ! (klon, llm, nqmx)
    ! (humidit\'e sp\'ecifique et fractions massiques des autres traceurs)

    REAL, intent(in):: omega(:, :) ! (klon, llm) vitesse verticale en Pa/s
    REAL, intent(out):: d_u(:, :) ! (klon, llm) tendance physique de "u" (m s-2)
    REAL, intent(out):: d_v(:, :) ! (klon, llm) tendance physique de "v" (m s-2)
    REAL, intent(out):: d_t(:, :) ! (klon, llm) tendance physique de "t" (K/s)

    REAL, intent(out):: d_qx(:, :, :) ! (klon, llm, nqmx)
    ! tendance physique de "qx" (s-1)

    ! Local:

    LOGICAL:: firstcal = .true.

    LOGICAL ok_gust ! pour activer l'effet des gust sur flux surface
    PARAMETER (ok_gust = .FALSE.)

    LOGICAL, PARAMETER:: check = .FALSE. 
    ! Verifier la conservation du modele en eau

    LOGICAL, PARAMETER:: ok_stratus = .FALSE.
    ! Ajouter artificiellement les stratus

    ! "slab" ocean
    REAL, save:: tslab(klon) ! temperature of ocean slab
    REAL, save:: seaice(klon) ! glace de mer (kg/m2) 
    REAL fluxo(klon) ! flux turbulents ocean-glace de mer 
    REAL fluxg(klon) ! flux turbulents ocean-atmosphere

    logical:: ok_journe = .false., ok_mensuel = .true., ok_instan = .false.
    ! sorties journalieres, mensuelles et instantanees dans les
    ! fichiers histday, histmth et histins

    LOGICAL ok_region ! sortir le fichier regional
    PARAMETER (ok_region = .FALSE.)

    ! pour phsystoke avec thermiques
    REAL fm_therm(klon, llm + 1)
    REAL entr_therm(klon, llm)
    real, save:: q2(klon, llm + 1, nbsrf)

    INTEGER, PARAMETER:: ivap = 1 ! indice de traceur pour vapeur d'eau
    INTEGER, PARAMETER:: iliq = 2 ! indice de traceur pour eau liquide

    REAL, save:: t_ancien(klon, llm), q_ancien(klon, llm)
    LOGICAL, save:: ancien_ok

    REAL d_t_dyn(klon, llm) ! tendance dynamique pour "t" (K/s)
    REAL d_q_dyn(klon, llm) ! tendance dynamique pour "q" (kg/kg/s)

    real da(klon, llm), phi(klon, llm, llm), mp(klon, llm)

    REAL swdn0(klon, llm + 1), swdn(klon, llm + 1)
    REAL swup0(klon, llm + 1), swup(klon, llm + 1)
    SAVE swdn0, swdn, swup0, swup

    REAL lwdn0(klon, llm + 1), lwdn(klon, llm + 1)
    REAL lwup0(klon, llm + 1), lwup(klon, llm + 1)
    SAVE lwdn0, lwdn, lwup0, lwup 

    ! Amip2
    ! variables a une pression donnee

    integer nlevSTD
    PARAMETER(nlevSTD = 17)

    ! prw: precipitable water
    real prw(klon)

    ! flwp, fiwp = Liquid Water Path & Ice Water Path (kg/m2)
    ! flwc, fiwc = Liquid Water Content & Ice Water Content (kg/kg)
    REAL flwp(klon), fiwp(klon)
    REAL flwc(klon, llm), fiwc(klon, llm)

    INTEGER kmax, lmax
    PARAMETER(kmax = 8, lmax = 8)
    INTEGER kmaxm1, lmaxm1
    PARAMETER(kmaxm1 = kmax-1, lmaxm1 = lmax-1)

    ! Variables propres a la physique

    INTEGER, save:: radpas
    ! Radiative transfer computations are made every "radpas" call to
    ! "physiq".

    REAL radsol(klon)
    SAVE radsol ! bilan radiatif au sol calcule par code radiatif

    INTEGER:: itap = 0 ! number of calls to "physiq"

    REAL, save:: ftsol(klon, nbsrf) ! skin temperature of surface fraction

    REAL, save:: ftsoil(klon, nsoilmx, nbsrf)
    ! soil temperature of surface fraction

    REAL, save:: fevap(klon, nbsrf) ! evaporation
    REAL fluxlat(klon, nbsrf)
    SAVE fluxlat

    REAL, save:: fqsurf(klon, nbsrf)
    ! humidite de l'air au contact de la surface

    REAL, save:: qsol(klon)
    ! column-density of water in soil, in kg m-2

    REAL, save:: fsnow(klon, nbsrf) ! epaisseur neigeuse
    REAL, save:: falbe(klon, nbsrf) ! albedo visible par type de surface

    ! Param\`etres de l'orographie \`a l'\'echelle sous-maille (OESM) :
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
    SAVE agesno ! age de la neige

    REAL run_off_lic_0(klon)
    SAVE run_off_lic_0
    !KE43
    ! Variables liees a la convection de K. Emanuel (sb):

    REAL Ma(klon, llm) ! undilute upward mass flux
    SAVE Ma
    REAL qcondc(klon, llm) ! in-cld water content from convect
    SAVE qcondc 
    REAL, save:: sig1(klon, llm), w01(klon, llm)
    REAL, save:: wd(klon)

    ! Variables locales pour la couche limite (al1):

    ! Variables locales:

    REAL cdragh(klon) ! drag coefficient pour T and Q
    REAL cdragm(klon) ! drag coefficient pour vent

    ! Pour phytrac :
    REAL ycoefh(klon, llm) ! coef d'echange pour phytrac
    REAL yu1(klon) ! vents dans la premiere couche U
    REAL yv1(klon) ! vents dans la premiere couche V
    REAL ffonte(klon, nbsrf) !Flux thermique utilise pour fondre la neige
    REAL fqcalving(klon, nbsrf) !Flux d'eau "perdue" par la surface 
    ! !et necessaire pour limiter la
    ! !hauteur de neige, en kg/m2/s
    REAL zxffonte(klon), zxfqcalving(klon)

    REAL pfrac_impa(klon, llm)! Produits des coefs lessivage impaction
    save pfrac_impa
    REAL pfrac_nucl(klon, llm)! Produits des coefs lessivage nucleation
    save pfrac_nucl
    REAL pfrac_1nucl(klon, llm)! Produits des coefs lessi nucl (alpha = 1)
    save pfrac_1nucl
    REAL frac_impa(klon, llm) ! fractions d'aerosols lessivees (impaction)
    REAL frac_nucl(klon, llm) ! idem (nucleation)

    REAL, save:: rain_fall(klon)
    ! liquid water mass flux (kg/m2/s), positive down

    REAL, save:: snow_fall(klon)
    ! solid water mass flux (kg/m2/s), positive down

    REAL rain_tiedtke(klon), snow_tiedtke(klon)

    REAL evap(klon), devap(klon) ! evaporation and its derivative
    REAL sens(klon), dsens(klon) ! chaleur sensible et sa derivee
    REAL dlw(klon) ! derivee infra rouge
    SAVE dlw
    REAL bils(klon) ! bilan de chaleur au sol
    REAL, save:: fder(klon) ! Derive de flux (sensible et latente) 
    REAL ve(klon) ! integr. verticale du transport meri. de l'energie
    REAL vq(klon) ! integr. verticale du transport meri. de l'eau
    REAL ue(klon) ! integr. verticale du transport zonal de l'energie
    REAL uq(klon) ! integr. verticale du transport zonal de l'eau

    REAL, save:: frugs(klon, nbsrf) ! longueur de rugosite
    REAL zxrugs(klon) ! longueur de rugosite

    ! Conditions aux limites

    INTEGER julien
    INTEGER, SAVE:: lmt_pas ! number of time steps of "physics" per day
    REAL, save:: pctsrf(klon, nbsrf) ! percentage of surface
    REAL pctsrf_new(klon, nbsrf) ! pourcentage surfaces issus d'ORCHIDEE
    REAL, save:: albsol(klon) ! albedo du sol total visible
    REAL, SAVE:: wo(klon, llm) ! column density of ozone in a cell, in kDU

    ! Declaration des procedures appelees

    EXTERNAL nuage ! calculer les proprietes radiatives
    EXTERNAL transp ! transport total de l'eau et de l'energie

    ! Variables locales

    real, save:: clwcon(klon, llm), rnebcon(klon, llm)
    real, save:: clwcon0(klon, llm), rnebcon0(klon, llm)

    REAL rhcl(klon, llm) ! humiditi relative ciel clair
    REAL dialiq(klon, llm) ! eau liquide nuageuse
    REAL diafra(klon, llm) ! fraction nuageuse
    REAL cldliq(klon, llm) ! eau liquide nuageuse
    REAL cldfra(klon, llm) ! fraction nuageuse
    REAL cldtau(klon, llm) ! epaisseur optique
    REAL cldemi(klon, llm) ! emissivite infrarouge

    REAL fluxq(klon, llm, nbsrf) ! flux turbulent d'humidite
    REAL fluxt(klon, llm, nbsrf) ! flux turbulent de chaleur
    REAL fluxu(klon, llm, nbsrf) ! flux turbulent de vitesse u
    REAL fluxv(klon, llm, nbsrf) ! flux turbulent de vitesse v

    REAL zxfluxt(klon, llm)
    REAL zxfluxq(klon, llm)
    REAL zxfluxu(klon, llm)
    REAL zxfluxv(klon, llm)

    ! Le rayonnement n'est pas calcul\'e tous les pas, il faut donc que
    ! les variables soient r\'emanentes.
    REAL, save:: heat(klon, llm) ! chauffage solaire
    REAL, save:: heat0(klon, llm) ! chauffage solaire ciel clair
    REAL, save:: cool(klon, llm) ! refroidissement infrarouge
    REAL, save:: cool0(klon, llm) ! refroidissement infrarouge ciel clair
    REAL, save:: topsw(klon), toplw(klon), solsw(klon)
    REAL, save:: sollw(klon) ! rayonnement infrarouge montant \`a la surface
    real, save:: sollwdown(klon) ! downward LW flux at surface
    REAL, save:: topsw0(klon), toplw0(klon), solsw0(klon), sollw0(klon)
    REAL, save:: albpla(klon)
    REAL fsollw(klon, nbsrf) ! bilan flux IR pour chaque sous surface
    REAL fsolsw(klon, nbsrf) ! flux solaire absorb. pour chaque sous surface

    REAL conv_q(klon, llm) ! convergence de l'humidite (kg/kg/s)
    REAL conv_t(klon, llm) ! convergence of temperature (K/s)

    REAL cldl(klon), cldm(klon), cldh(klon) !nuages bas, moyen et haut
    REAL cldt(klon), cldq(klon) !nuage total, eau liquide integree

    REAL zxtsol(klon), zxqsurf(klon), zxsnow(klon), zxfluxlat(klon)

    REAL dist, mu0(klon), fract(klon)
    real longi
    REAL z_avant(klon), z_apres(klon), z_factor(klon)
    REAL za, zb
    REAL zx_t, zx_qs, zcor
    real zqsat(klon, llm)
    INTEGER i, k, iq, nsrf
    REAL, PARAMETER:: t_coup = 234.
    REAL zphi(klon, llm)

    ! cf. AM Variables locales pour la CLA (hbtm2)

    REAL, SAVE:: pblh(klon, nbsrf) ! Hauteur de couche limite
    REAL, SAVE:: plcl(klon, nbsrf) ! Niveau de condensation de la CLA
    REAL, SAVE:: capCL(klon, nbsrf) ! CAPE de couche limite
    REAL, SAVE:: oliqCL(klon, nbsrf) ! eau_liqu integree de couche limite
    REAL, SAVE:: cteiCL(klon, nbsrf) ! cloud top instab. crit. couche limite
    REAL, SAVE:: pblt(klon, nbsrf) ! T a la Hauteur de couche limite
    REAL, SAVE:: therm(klon, nbsrf)
    REAL, SAVE:: trmb1(klon, nbsrf) ! deep_cape
    REAL, SAVE:: trmb2(klon, nbsrf) ! inhibition 
    REAL, SAVE:: trmb3(klon, nbsrf) ! Point Omega
    ! Grdeurs de sorties
    REAL s_pblh(klon), s_lcl(klon), s_capCL(klon)
    REAL s_oliqCL(klon), s_cteiCL(klon), s_pblt(klon)
    REAL s_therm(klon), s_trmb1(klon), s_trmb2(klon)
    REAL s_trmb3(klon)

    ! Variables locales pour la convection de K. Emanuel :

    REAL upwd(klon, llm) ! saturated updraft mass flux
    REAL dnwd(klon, llm) ! saturated downdraft mass flux
    REAL dnwd0(klon, llm) ! unsaturated downdraft mass flux
    REAL cape(klon) ! CAPE
    SAVE cape

    INTEGER iflagctrl(klon) ! flag fonctionnement de convect

    ! Variables du changement

    ! con: convection
    ! lsc: large scale condensation
    ! ajs: ajustement sec
    ! eva: \'evaporation de l'eau liquide nuageuse
    ! vdf: vertical diffusion in boundary layer
    REAL d_t_con(klon, llm), d_q_con(klon, llm)
    REAL d_u_con(klon, llm), d_v_con(klon, llm)
    REAL d_t_lsc(klon, llm), d_q_lsc(klon, llm), d_ql_lsc(klon, llm)
    REAL d_t_ajs(klon, llm), d_q_ajs(klon, llm)
    REAL d_u_ajs(klon, llm), d_v_ajs(klon, llm)
    REAL rneb(klon, llm)

    REAL mfu(klon, llm), mfd(klon, llm)
    REAL pen_u(klon, llm), pen_d(klon, llm)
    REAL pde_u(klon, llm), pde_d(klon, llm)
    INTEGER kcbot(klon), kctop(klon), kdtop(klon)
    REAL pmflxr(klon, llm + 1), pmflxs(klon, llm + 1)
    REAL prfl(klon, llm + 1), psfl(klon, llm + 1)

    INTEGER, save:: ibas_con(klon), itop_con(klon)

    REAL rain_con(klon), rain_lsc(klon)
    REAL snow_con(klon), snow_lsc(klon)
    REAL d_ts(klon, nbsrf)

    REAL d_u_vdf(klon, llm), d_v_vdf(klon, llm)
    REAL d_t_vdf(klon, llm), d_q_vdf(klon, llm)

    REAL d_u_oro(klon, llm), d_v_oro(klon, llm)
    REAL d_t_oro(klon, llm)
    REAL d_u_lif(klon, llm), d_v_lif(klon, llm)
    REAL d_t_lif(klon, llm)

    REAL, save:: ratqs(klon, llm)
    real ratqss(klon, llm), ratqsc(klon, llm)
    real:: ratqsbas = 0.01, ratqshaut = 0.3

    ! Parametres lies au nouveau schema de nuages (SB, PDF)
    real:: fact_cldcon = 0.375
    real:: facttemps = 1.e-4
    logical:: ok_newmicro = .true.
    real facteur

    integer:: iflag_cldcon = 1
    logical ptconv(klon, llm)

    ! Variables locales pour effectuer les appels en s\'erie :

    REAL t_seri(klon, llm), q_seri(klon, llm)
    REAL ql_seri(klon, llm)
    REAL u_seri(klon, llm), v_seri(klon, llm)
    REAL tr_seri(klon, llm, nqmx - 2)

    REAL zx_rh(klon, llm)

    REAL zustrdr(klon), zvstrdr(klon)
    REAL zustrli(klon), zvstrli(klon)
    REAL zustrph(klon), zvstrph(klon)
    REAL aam, torsfc

    REAL zx_tmp_fi2d(klon) ! variable temporaire grille physique

    INTEGER, SAVE:: nid_ins

    REAL ve_lay(klon, llm) ! transport meri. de l'energie a chaque niveau vert.
    REAL vq_lay(klon, llm) ! transport meri. de l'eau a chaque niveau vert.
    REAL ue_lay(klon, llm) ! transport zonal de l'energie a chaque niveau vert.
    REAL uq_lay(klon, llm) ! transport zonal de l'eau a chaque niveau vert.

    real date0

    ! Variables li\'ees au bilan d'\'energie et d'enthalpie :
    REAL ztsol(klon)
    REAL d_h_vcol, d_qt, d_ec
    REAL, SAVE:: d_h_vcol_phy
    REAL zero_v(klon)
    CHARACTER(LEN = 20) tit
    INTEGER:: ip_ebil = 0 ! print level for energy conservation diagnostics
    INTEGER:: if_ebil = 0 ! verbosity for diagnostics of energy conservation 

    REAL d_t_ec(klon, llm) ! tendance due \`a la conversion Ec -> E thermique
    REAL ZRCPD

    REAL t2m(klon, nbsrf), q2m(klon, nbsrf) ! temperature and humidity at 2 m
    REAL u10m(klon, nbsrf), v10m(klon, nbsrf) ! vents a 10 m
    REAL zt2m(klon), zq2m(klon) ! temp., hum. 2 m moyenne s/ 1 maille
    REAL zu10m(klon), zv10m(klon) ! vents a 10 m moyennes s/1 maille

    ! Aerosol effects:

    REAL sulfate(klon, llm) ! SO4 aerosol concentration (micro g/m3)

    REAL, save:: sulfate_pi(klon, llm)
    ! SO4 aerosol concentration, in micro g/m3, pre-industrial value

    REAL cldtaupi(klon, llm)
    ! cloud optical thickness for pre-industrial (pi) aerosols

    REAL re(klon, llm) ! Cloud droplet effective radius
    REAL fl(klon, llm) ! denominator of re

    ! Aerosol optical properties
    REAL, save:: tau_ae(klon, llm, 2), piz_ae(klon, llm, 2)
    REAL, save:: cg_ae(klon, llm, 2)

    REAL topswad(klon), solswad(klon) ! aerosol direct effect
    REAL topswai(klon), solswai(klon) ! aerosol indirect effect

    REAL aerindex(klon) ! POLDER aerosol index

    LOGICAL:: ok_ade = .false. ! apply aerosol direct effect
    LOGICAL:: ok_aie = .false. ! apply aerosol indirect effect

    REAL:: bl95_b0 = 2., bl95_b1 = 0.2
    ! Parameters in equation (D) of Boucher and Lohmann (1995, Tellus
    ! B). They link cloud droplet number concentration to aerosol mass
    ! concentration.

    SAVE u10m
    SAVE v10m
    SAVE t2m
    SAVE q2m
    SAVE ffonte
    SAVE fqcalving
    SAVE rain_con
    SAVE snow_con
    SAVE topswai
    SAVE topswad
    SAVE solswai
    SAVE solswad
    SAVE d_u_con
    SAVE d_v_con

    real zmasse(klon, llm) 
    ! (column-density of mass of air in a cell, in kg m-2)

    real, parameter:: dobson_u = 2.1415e-05 ! Dobson unit, in kg m-2
    integer, save:: ncid_startphy

    namelist /physiq_nml/ ok_journe, ok_mensuel, ok_instan, fact_cldcon, &
         facttemps, ok_newmicro, iflag_cldcon, ratqsbas, ratqshaut, if_ebil, &
         ok_ade, ok_aie, bl95_b0, bl95_b1, iflag_thermals, nsplit_thermals

    !----------------------------------------------------------------

    IF (if_ebil >= 1) zero_v = 0.
    IF (nqmx < 2) CALL abort_gcm('physiq', &
         'eaux vapeur et liquide sont indispensables', 1)

    test_firstcal: IF (firstcal) THEN
       ! initialiser
       u10m = 0.
       v10m = 0.
       t2m = 0.
       q2m = 0.
       ffonte = 0.
       fqcalving = 0.
       piz_ae = 0.
       tau_ae = 0.
       cg_ae = 0.
       rain_con = 0.
       snow_con = 0.
       topswai = 0.
       topswad = 0.
       solswai = 0.
       solswad = 0.

       d_u_con = 0.
       d_v_con = 0.
       rnebcon0 = 0.
       clwcon0 = 0.
       rnebcon = 0.
       clwcon = 0.

       pblh =0. ! Hauteur de couche limite
       plcl =0. ! Niveau de condensation de la CLA
       capCL =0. ! CAPE de couche limite
       oliqCL =0. ! eau_liqu integree de couche limite
       cteiCL =0. ! cloud top instab. crit. couche limite
       pblt =0. ! T a la Hauteur de couche limite
       therm =0.
       trmb1 =0. ! deep_cape
       trmb2 =0. ! inhibition 
       trmb3 =0. ! Point Omega

       IF (if_ebil >= 1) d_h_vcol_phy = 0.

       iflag_thermals = 0
       nsplit_thermals = 1
       print *, "Enter namelist 'physiq_nml'."
       read(unit=*, nml=physiq_nml)
       write(unit_nml, nml=physiq_nml)

       call conf_phys

       ! Initialiser les compteurs:

       frugs = 0.
       CALL phyetat0(pctsrf, ftsol, ftsoil, tslab, seaice, fqsurf, qsol, &
            fsnow, falbe, fevap, rain_fall, snow_fall, solsw, sollw, dlw, &
            radsol, frugs, agesno, zmea, zstd, zsig, zgam, zthe, zpic, zval, &
            t_ancien, q_ancien, ancien_ok, rnebcon, ratqs, clwcon, &
            run_off_lic_0, sig1, w01, ncid_startphy)

       ! ATTENTION : il faudra a terme relire q2 dans l'etat initial
       q2 = 1e-8

       lmt_pas = day_step / iphysiq
       print *, 'Number of time steps of "physics" per day: ', lmt_pas

       radpas = lmt_pas / nbapp_rad

       ! On remet le calendrier a zero
       IF (raz_date) itau_phy = 0

       CALL printflag(radpas, ok_journe, ok_instan, ok_region)

       ! Initialisation pour le sch\'ema de convection d'Emanuel :
       IF (iflag_con >= 3) THEN
          ibas_con = 1
          itop_con = 1
       ENDIF

       IF (ok_orodr) THEN
          rugoro = MAX(1e-5, zstd * zsig / 2)
          CALL SUGWD(paprs, play)
       else
          rugoro = 0.
       ENDIF

       ecrit_ins = NINT(ecrit_ins/dtphys)
       ecrit_hf = NINT(ecrit_hf/dtphys)
       ecrit_mth = NINT(ecrit_mth/dtphys)
       ecrit_tra = NINT(86400.*ecrit_tra/dtphys)
       ecrit_reg = NINT(ecrit_reg/dtphys)

       ! Initialisation des sorties

       call ini_histins(dtphys, ok_instan, nid_ins)
       CALL ymds2ju(annee_ref, 1, day_ref, 0., date0)
       ! Positionner date0 pour initialisation de ORCHIDEE
       print *, 'physiq date0: ', date0
       CALL phyredem0(lmt_pas)
    ENDIF test_firstcal

    ! We will modify variables *_seri and we will not touch variables
    ! u, v, t, qx:
    t_seri = t
    u_seri = u
    v_seri = v
    q_seri = qx(:, :, ivap)
    ql_seri = qx(:, :, iliq)
    tr_seri = qx(:, :, 3:nqmx)

    ztsol = sum(ftsol * pctsrf, dim = 2)

    IF (if_ebil >= 1) THEN 
       tit = 'after dynamics'
       CALL diagetpq(airephy, tit, ip_ebil, 1, 1, dtphys, t_seri, q_seri, &
            ql_seri, u_seri, v_seri, paprs, d_h_vcol, d_qt, d_ec)
       ! Comme les tendances de la physique sont ajout\'es dans la
       !  dynamique, la variation d'enthalpie par la dynamique devrait
       !  \^etre \'egale \`a la variation de la physique au pas de temps
       !  pr\'ec\'edent.  Donc la somme de ces 2 variations devrait \^etre
       !  nulle.
       call diagphy(airephy, tit, ip_ebil, zero_v, zero_v, zero_v, zero_v, &
            zero_v, zero_v, zero_v, zero_v, ztsol, d_h_vcol + d_h_vcol_phy, &
            d_qt, 0.)
    END IF

    ! Diagnostic de la tendance dynamique :
    IF (ancien_ok) THEN
       DO k = 1, llm
          DO i = 1, klon
             d_t_dyn(i, k) = (t_seri(i, k) - t_ancien(i, k)) / dtphys
             d_q_dyn(i, k) = (q_seri(i, k) - q_ancien(i, k)) / dtphys
          ENDDO
       ENDDO
    ELSE
       DO k = 1, llm
          DO i = 1, klon
             d_t_dyn(i, k) = 0.
             d_q_dyn(i, k) = 0.
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

    ! Check temperatures:
    CALL hgardfou(t_seri, ftsol)

    ! Incrémenter le compteur de la physique
    itap = itap + 1
    julien = MOD(dayvrai, 360)
    if (julien == 0) julien = 360

    forall (k = 1: llm) zmasse(:, k) = (paprs(:, k) - paprs(:, k + 1)) / rg

    ! Prescrire l'ozone :
    wo = ozonecm(REAL(julien), paprs)

    ! \'Evaporation de l'eau liquide nuageuse :
    DO k = 1, llm
       DO i = 1, klon
          zb = MAX(0., ql_seri(i, k))
          t_seri(i, k) = t_seri(i, k) &
               - zb * RLVTT / RCPD / (1. + RVTMP2 * q_seri(i, k))
          q_seri(i, k) = q_seri(i, k) + zb
       ENDDO
    ENDDO
    ql_seri = 0.

    IF (if_ebil >= 2) THEN 
       tit = 'after reevap'
       CALL diagetpq(airephy, tit, ip_ebil, 2, 1, dtphys, t_seri, q_seri, &
            ql_seri, u_seri, v_seri, paprs, d_h_vcol, d_qt, d_ec)
       call diagphy(airephy, tit, ip_ebil, zero_v, zero_v, zero_v, zero_v, &
            zero_v, zero_v, zero_v, zero_v, ztsol, d_h_vcol, d_qt, d_ec)
    END IF

    frugs = MAX(frugs, 0.000015)
    zxrugs = sum(frugs * pctsrf, dim = 2)

    ! Calculs nécessaires au calcul de l'albedo dans l'interface avec
    ! la surface.

    CALL orbite(REAL(julien), longi, dist)
    IF (cycle_diurne) THEN
       CALL zenang(longi, time, dtphys * radpas, mu0, fract)
    ELSE
       mu0 = -999.999
    ENDIF

    ! Calcul de l'abedo moyen par maille
    albsol = sum(falbe * pctsrf, dim = 2)

    ! R\'epartition sous maille des flux longwave et shortwave
    ! R\'epartition du longwave par sous-surface lin\'earis\'ee

    forall (nsrf = 1: nbsrf)
       fsollw(:, nsrf) = sollw + 4. * RSIGMA * ztsol**3 &
            * (ztsol - ftsol(:, nsrf))
       fsolsw(:, nsrf) = solsw * (1. - falbe(:, nsrf)) / (1. - albsol)
    END forall

    fder = dlw

    ! Couche limite:

    CALL clmain(dtphys, itap, pctsrf, pctsrf_new, t_seri, q_seri, u_seri, &
         v_seri, julien, mu0, co2_ppm, ftsol, cdmmax, cdhmax, ksta, ksta_ter, &
         ok_kzmin, ftsoil, qsol, paprs, play, fsnow, fqsurf, fevap, falbe, &
         fluxlat, rain_fall, snow_fall, fsolsw, fsollw, fder, rlat, frugs, &
         firstcal, agesno, rugoro, d_t_vdf, d_q_vdf, d_u_vdf, d_v_vdf, d_ts, &
         fluxt, fluxq, fluxu, fluxv, cdragh, cdragm, q2, dsens, devap, &
         ycoefh, yu1, yv1, t2m, q2m, u10m, v10m, pblh, capCL, oliqCL, cteiCL, &
         pblT, therm, trmb1, trmb2, trmb3, plcl, fqcalving, ffonte, &
         run_off_lic_0, fluxo, fluxg, tslab)

    ! Incr\'ementation des flux

    zxfluxt = 0.
    zxfluxq = 0.
    zxfluxu = 0.
    zxfluxv = 0.
    DO nsrf = 1, nbsrf
       DO k = 1, llm
          DO i = 1, klon
             zxfluxt(i, k) = zxfluxt(i, k) + fluxt(i, k, nsrf) * pctsrf(i, nsrf)
             zxfluxq(i, k) = zxfluxq(i, k) + fluxq(i, k, nsrf) * pctsrf(i, nsrf)
             zxfluxu(i, k) = zxfluxu(i, k) + fluxu(i, k, nsrf) * pctsrf(i, nsrf)
             zxfluxv(i, k) = zxfluxv(i, k) + fluxv(i, k, nsrf) * pctsrf(i, nsrf)
          END DO
       END DO
    END DO
    DO i = 1, klon
       sens(i) = - zxfluxt(i, 1) ! flux de chaleur sensible au sol
       evap(i) = - zxfluxq(i, 1) ! flux d'\'evaporation au sol
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
       tit = 'after clmain'
       CALL diagetpq(airephy, tit, ip_ebil, 2, 2, dtphys, t_seri, q_seri, &
            ql_seri, u_seri, v_seri, paprs, d_h_vcol, d_qt, d_ec)
       call diagphy(airephy, tit, ip_ebil, zero_v, zero_v, zero_v, zero_v, &
            sens, evap, zero_v, zero_v, ztsol, d_h_vcol, d_qt, d_ec)
    END IF

    ! Update surface temperature:

    DO i = 1, klon
       zxtsol(i) = 0.
       zxfluxlat(i) = 0.

       zt2m(i) = 0.
       zq2m(i) = 0.
       zu10m(i) = 0.
       zv10m(i) = 0.
       zxffonte(i) = 0.
       zxfqcalving(i) = 0.

       s_pblh(i) = 0. 
       s_lcl(i) = 0. 
       s_capCL(i) = 0.
       s_oliqCL(i) = 0.
       s_cteiCL(i) = 0.
       s_pblT(i) = 0.
       s_therm(i) = 0.
       s_trmb1(i) = 0.
       s_trmb2(i) = 0.
       s_trmb3(i) = 0.

       IF (abs(pctsrf(i, is_ter) + pctsrf(i, is_lic) + pctsrf(i, is_oce) &
            + pctsrf(i, is_sic) - 1.)  >  EPSFRA) print *, &
            'physiq : probl\`eme sous surface au point ', i, &
            pctsrf(i, 1 : nbsrf)
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
          zxfqcalving(i) = zxfqcalving(i) + &
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

    ! Si une sous-fraction n'existe pas, elle prend la température moyenne :
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          IF (pctsrf(i, nsrf) < epsfra) ftsol(i, nsrf) = zxtsol(i)

          IF (pctsrf(i, nsrf) < epsfra) t2m(i, nsrf) = zt2m(i)
          IF (pctsrf(i, nsrf) < epsfra) q2m(i, nsrf) = zq2m(i)
          IF (pctsrf(i, nsrf) < epsfra) u10m(i, nsrf) = zu10m(i)
          IF (pctsrf(i, nsrf) < epsfra) v10m(i, nsrf) = zv10m(i)
          IF (pctsrf(i, nsrf) < epsfra) ffonte(i, nsrf) = zxffonte(i)
          IF (pctsrf(i, nsrf) < epsfra) &
               fqcalving(i, nsrf) = zxfqcalving(i)
          IF (pctsrf(i, nsrf) < epsfra) pblh(i, nsrf) = s_pblh(i)
          IF (pctsrf(i, nsrf) < epsfra) plcl(i, nsrf) = s_lcl(i)
          IF (pctsrf(i, nsrf) < epsfra) capCL(i, nsrf) = s_capCL(i)
          IF (pctsrf(i, nsrf) < epsfra) oliqCL(i, nsrf) = s_oliqCL(i)
          IF (pctsrf(i, nsrf) < epsfra) cteiCL(i, nsrf) = s_cteiCL(i)
          IF (pctsrf(i, nsrf) < epsfra) pblT(i, nsrf) = s_pblT(i)
          IF (pctsrf(i, nsrf) < epsfra) therm(i, nsrf) = s_therm(i)
          IF (pctsrf(i, nsrf) < epsfra) trmb1(i, nsrf) = s_trmb1(i)
          IF (pctsrf(i, nsrf) < epsfra) trmb2(i, nsrf) = s_trmb2(i)
          IF (pctsrf(i, nsrf) < epsfra) trmb3(i, nsrf) = s_trmb3(i)
       ENDDO
    ENDDO

    ! Calculer la dérive du flux infrarouge

    DO i = 1, klon
       dlw(i) = - 4. * RSIGMA * zxtsol(i)**3 
    ENDDO

    IF (check) print *, "avantcon = ", qcheck(paprs, q_seri, ql_seri)

    ! Appeler la convection (au choix)

    if (iflag_con == 2) then
       conv_q = d_q_dyn + d_q_vdf / dtphys
       conv_t = d_t_dyn + d_t_vdf / dtphys
       z_avant = sum((q_seri + ql_seri) * zmasse, dim=2)
       CALL conflx(dtphys, paprs, play, t_seri(:, llm:1:-1), &
            q_seri(:, llm:1:-1), conv_t, conv_q, zxfluxq(:, 1), omega, &
            d_t_con, d_q_con, rain_con, snow_con, mfu(:, llm:1:-1), &
            mfd(:, llm:1:-1), pen_u, pde_u, pen_d, pde_d, kcbot, kctop, &
            kdtop, pmflxr, pmflxs)
       WHERE (rain_con < 0.) rain_con = 0.
       WHERE (snow_con < 0.) snow_con = 0.
       ibas_con = llm + 1 - kcbot
       itop_con = llm + 1 - kctop
    else
       ! iflag_con >= 3

       da = 0.
       mp = 0.
       phi = 0.
       CALL concvl(dtphys, paprs, play, t_seri, q_seri, u_seri, v_seri, sig1, &
            w01, d_t_con, d_q_con, d_u_con, d_v_con, rain_con, snow_con, &
            ibas_con, itop_con, upwd, dnwd, dnwd0, Ma, cape, iflagctrl, &
            qcondc, wd, pmflxr, pmflxs, da, phi, mp)
       clwcon0 = qcondc
       mfu = upwd + dnwd
       IF (.NOT. ok_gust) wd = 0.

       IF (thermcep) THEN
          zqsat = MIN(0.5, r2es * FOEEW(t_seri, rtt >= t_seri) / play)
          zqsat = zqsat / (1. - retv * zqsat)
       ELSE
          zqsat = merge(qsats(t_seri), qsatl(t_seri), t_seri < t_coup) / play
       ENDIF

       ! Properties of convective clouds
       clwcon0 = fact_cldcon * clwcon0
       call clouds_gno(klon, llm, q_seri, zqsat, clwcon0, ptconv, ratqsc, &
            rnebcon0)

       mfd = 0.
       pen_u = 0.
       pen_d = 0.
       pde_d = 0.
       pde_u = 0.
    END if

    DO k = 1, llm
       DO i = 1, klon
          t_seri(i, k) = t_seri(i, k) + d_t_con(i, k)
          q_seri(i, k) = q_seri(i, k) + d_q_con(i, k)
          u_seri(i, k) = u_seri(i, k) + d_u_con(i, k)
          v_seri(i, k) = v_seri(i, k) + d_v_con(i, k)
       ENDDO
    ENDDO

    IF (if_ebil >= 2) THEN 
       tit = 'after convect'
       CALL diagetpq(airephy, tit, ip_ebil, 2, 2, dtphys, t_seri, q_seri, &
            ql_seri, u_seri, v_seri, paprs, d_h_vcol, d_qt, d_ec)
       call diagphy(airephy, tit, ip_ebil, zero_v, zero_v, zero_v, zero_v, &
            zero_v, zero_v, rain_con, snow_con, ztsol, d_h_vcol, d_qt, d_ec)
    END IF

    IF (check) THEN
       za = qcheck(paprs, q_seri, ql_seri)
       print *, "aprescon = ", za
       zx_t = 0.
       za = 0.
       DO i = 1, klon
          za = za + airephy(i)/REAL(klon)
          zx_t = zx_t + (rain_con(i)+ &
               snow_con(i))*airephy(i)/REAL(klon)
       ENDDO
       zx_t = zx_t/za*dtphys
       print *, "Precip = ", zx_t
    ENDIF

    IF (iflag_con == 2) THEN
       z_apres = sum((q_seri + ql_seri) * zmasse, dim=2)
       z_factor = (z_avant - (rain_con + snow_con) * dtphys) / z_apres
       DO k = 1, llm
          DO i = 1, klon
             IF (z_factor(i) > 1. + 1E-8 .OR. z_factor(i) < 1. - 1E-8) THEN
                q_seri(i, k) = q_seri(i, k) * z_factor(i)
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    ! Convection s\`eche (thermiques ou ajustement)

    d_t_ajs = 0.
    d_u_ajs = 0.
    d_v_ajs = 0.
    d_q_ajs = 0.
    fm_therm = 0.
    entr_therm = 0.

    if (iflag_thermals == 0) then
       ! Ajustement sec
       CALL ajsec(paprs, play, t_seri, q_seri, d_t_ajs, d_q_ajs)
       t_seri = t_seri + d_t_ajs
       q_seri = q_seri + d_q_ajs
    else
       ! Thermiques
       call calltherm(dtphys, play, paprs, pphi, u_seri, v_seri, t_seri, &
            q_seri, d_u_ajs, d_v_ajs, d_t_ajs, d_q_ajs, fm_therm, entr_therm)
    endif

    IF (if_ebil >= 2) THEN 
       tit = 'after dry_adjust'
       CALL diagetpq(airephy, tit, ip_ebil, 2, 2, dtphys, t_seri, q_seri, &
            ql_seri, u_seri, v_seri, paprs, d_h_vcol, d_qt, d_ec)
    END IF

    ! Caclul des ratqs

    ! ratqs convectifs \`a l'ancienne en fonction de (q(z = 0) - q) / q
    ! on \'ecrase le tableau ratqsc calcul\'e par clouds_gno
    if (iflag_cldcon == 1) then
       do k = 1, llm
          do i = 1, klon
             if(ptconv(i, k)) then
                ratqsc(i, k) = ratqsbas + fact_cldcon &
                     * (q_seri(i, 1) - q_seri(i, k)) / q_seri(i, k)
             else
                ratqsc(i, k) = 0.
             endif
          enddo
       enddo
    endif

    ! ratqs stables
    do k = 1, llm
       do i = 1, klon
          ratqss(i, k) = ratqsbas + (ratqshaut - ratqsbas) &
               * min((paprs(i, 1) - play(i, k)) / (paprs(i, 1) - 3e4), 1.) 
       enddo
    enddo

    ! ratqs final
    if (iflag_cldcon == 1 .or. iflag_cldcon == 2) then
       ! les ratqs sont une conbinaison de ratqss et ratqsc
       ! ratqs final
       ! 1e4 (en gros 3 heures), en dur pour le moment, est le temps de
       ! relaxation des ratqs
       ratqs = max(ratqs * exp(- dtphys * facttemps), ratqss)
       ratqs = max(ratqs, ratqsc)
    else
       ! on ne prend que le ratqs stable pour fisrtilp
       ratqs = ratqss
    endif

    CALL fisrtilp(dtphys, paprs, play, t_seri, q_seri, ptconv, ratqs, &
         d_t_lsc, d_q_lsc, d_ql_lsc, rneb, cldliq, rain_lsc, snow_lsc, &
         pfrac_impa, pfrac_nucl, pfrac_1nucl, frac_impa, frac_nucl, prfl, &
         psfl, rhcl)

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
       za = qcheck(paprs, q_seri, ql_seri)
       print *, "apresilp = ", za
       zx_t = 0.
       za = 0.
       DO i = 1, klon
          za = za + airephy(i)/REAL(klon)
          zx_t = zx_t + (rain_lsc(i) &
               + snow_lsc(i))*airephy(i)/REAL(klon)
       ENDDO
       zx_t = zx_t/za*dtphys
       print *, "Precip = ", zx_t
    ENDIF

    IF (if_ebil >= 2) THEN 
       tit = 'after fisrt'
       CALL diagetpq(airephy, tit, ip_ebil, 2, 2, dtphys, t_seri, q_seri, &
            ql_seri, u_seri, v_seri, paprs, d_h_vcol, d_qt, d_ec)
       call diagphy(airephy, tit, ip_ebil, zero_v, zero_v, zero_v, zero_v, &
            zero_v, zero_v, rain_lsc, snow_lsc, ztsol, d_h_vcol, d_qt, d_ec)
    END IF

    ! PRESCRIPTION DES NUAGES POUR LE RAYONNEMENT

    ! 1. NUAGES CONVECTIFS

    IF (iflag_cldcon <= -1) THEN
       ! seulement pour Tiedtke
       snow_tiedtke = 0.
       if (iflag_cldcon == -1) then
          rain_tiedtke = rain_con
       else
          rain_tiedtke = 0.
          do k = 1, llm
             do i = 1, klon
                if (d_q_con(i, k) < 0.) then
                   rain_tiedtke(i) = rain_tiedtke(i)-d_q_con(i, k)/dtphys &
                        *zmasse(i, k)
                endif
             enddo
          enddo
       endif

       ! Nuages diagnostiques pour Tiedtke
       CALL diagcld1(paprs, play, rain_tiedtke, snow_tiedtke, ibas_con, &
            itop_con, diafra, dialiq)
       DO k = 1, llm
          DO i = 1, klon
             IF (diafra(i, k) > cldfra(i, k)) THEN
                cldliq(i, k) = dialiq(i, k)
                cldfra(i, k) = diafra(i, k)
             ENDIF
          ENDDO
       ENDDO
    ELSE IF (iflag_cldcon == 3) THEN
       ! On prend pour les nuages convectifs le maximum du calcul de
       ! la convection et du calcul du pas de temps pr\'ec\'edent diminu\'e
       ! d'un facteur facttemps.
       facteur = dtphys * facttemps
       do k = 1, llm
          do i = 1, klon
             rnebcon(i, k) = rnebcon(i, k) * facteur
             if (rnebcon0(i, k) * clwcon0(i, k) &
                  > rnebcon(i, k) * clwcon(i, k)) then
                rnebcon(i, k) = rnebcon0(i, k)
                clwcon(i, k) = clwcon0(i, k)
             endif
          enddo
       enddo

       ! On prend la somme des fractions nuageuses et des contenus en eau
       cldfra = min(max(cldfra, rnebcon), 1.)
       cldliq = cldliq + rnebcon*clwcon
    ENDIF

    ! 2. Nuages stratiformes

    IF (ok_stratus) THEN
       CALL diagcld2(paprs, play, t_seri, q_seri, diafra, dialiq)
       DO k = 1, llm
          DO i = 1, klon
             IF (diafra(i, k) > cldfra(i, k)) THEN
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

    IF (if_ebil >= 2) CALL diagetpq(airephy, "after diagcld", ip_ebil, 2, 2, &
         dtphys, t_seri, q_seri, ql_seri, u_seri, v_seri, paprs, d_h_vcol, &
         d_qt, d_ec)

    ! Humidit\'e relative pour diagnostic :
    DO k = 1, llm
       DO i = 1, klon
          zx_t = t_seri(i, k)
          IF (thermcep) THEN
             zx_qs = r2es * FOEEW(zx_t, rtt >= zx_t)/play(i, k)
             zx_qs = MIN(0.5, zx_qs)
             zcor = 1./(1.-retv*zx_qs)
             zx_qs = zx_qs*zcor
          ELSE
             IF (zx_t < t_coup) THEN
                zx_qs = qsats(zx_t)/play(i, k)
             ELSE
                zx_qs = qsatl(zx_t)/play(i, k)
             ENDIF
          ENDIF
          zx_rh(i, k) = q_seri(i, k)/zx_qs
          zqsat(i, k) = zx_qs
       ENDDO
    ENDDO

    ! Introduce the aerosol direct and first indirect radiative forcings:
    IF (ok_ade .OR. ok_aie) THEN
       ! Get sulfate aerosol distribution :
       CALL readsulfate(dayvrai, time, firstcal, sulfate)
       CALL readsulfate_preind(dayvrai, time, firstcal, sulfate_pi)

       CALL aeropt(play, paprs, t_seri, sulfate, rhcl, tau_ae, piz_ae, cg_ae, &
            aerindex)
    ELSE
       tau_ae = 0.
       piz_ae = 0.
       cg_ae = 0.
    ENDIF

    ! Param\`etres optiques des nuages et quelques param\`etres pour
    ! diagnostics :
    if (ok_newmicro) then
       CALL newmicro(paprs, play, t_seri, cldliq, cldfra, cldtau, cldemi, &
            cldh, cldl, cldm, cldt, cldq, flwp, fiwp, flwc, fiwc, ok_aie, &
            sulfate, sulfate_pi, bl95_b0, bl95_b1, cldtaupi, re, fl)
    else
       CALL nuage(paprs, play, t_seri, cldliq, cldfra, cldtau, cldemi, cldh, &
            cldl, cldm, cldt, cldq, ok_aie, sulfate, sulfate_pi, bl95_b0, &
            bl95_b1, cldtaupi, re, fl)
    endif

    IF (MOD(itap - 1, radpas) == 0) THEN
       ! Appeler le rayonnement mais calculer tout d'abord l'albedo du sol.
       ! Calcul de l'abedo moyen par maille
       albsol = sum(falbe * pctsrf, dim = 2)

       ! Rayonnement (compatible Arpege-IFS) :
       CALL radlwsw(dist, mu0, fract, paprs, play, zxtsol, albsol, t_seri, &
            q_seri, wo, cldfra, cldemi, cldtau, heat, heat0, cool, cool0, &
            radsol, albpla, topsw, toplw, solsw, sollw, sollwdown, topsw0, &
            toplw0, solsw0, sollw0, lwdn0, lwdn, lwup0, lwup, swdn0, swdn, &
            swup0, swup, ok_ade, ok_aie, tau_ae, piz_ae, cg_ae, topswad, &
            solswad, cldtaupi, topswai, solswai)
    ENDIF

    ! Ajouter la tendance des rayonnements (tous les pas)

    DO k = 1, llm
       DO i = 1, klon
          t_seri(i, k) = t_seri(i, k) + (heat(i, k)-cool(i, k)) * dtphys/86400.
       ENDDO
    ENDDO

    IF (if_ebil >= 2) THEN 
       tit = 'after rad'
       CALL diagetpq(airephy, tit, ip_ebil, 2, 2, dtphys, t_seri, q_seri, &
            ql_seri, u_seri, v_seri, paprs, d_h_vcol, d_qt, d_ec)
       call diagphy(airephy, tit, ip_ebil, topsw, toplw, solsw, sollw, &
            zero_v, zero_v, zero_v, zero_v, ztsol, d_h_vcol, d_qt, d_ec)
    END IF

    ! Calculer l'hydrologie de la surface
    DO i = 1, klon
       zxqsurf(i) = 0.
       zxsnow(i) = 0.
    ENDDO
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          zxqsurf(i) = zxqsurf(i) + fqsurf(i, nsrf)*pctsrf(i, nsrf)
          zxsnow(i) = zxsnow(i) + fsnow(i, nsrf)*pctsrf(i, nsrf)
       ENDDO
    ENDDO

    ! Calculer le bilan du sol et la d\'erive de temp\'erature (couplage)

    DO i = 1, klon
       bils(i) = radsol(i) - sens(i) + zxfluxlat(i)
    ENDDO

    ! Param\'etrisation de l'orographie \`a l'\'echelle sous-maille :

    IF (ok_orodr) THEN
       ! selection des points pour lesquels le shema est actif:
       igwd = 0
       DO i = 1, klon
          itest(i) = 0
          IF (((zpic(i)-zmea(i)) > 100.).AND.(zstd(i) > 10.)) THEN
             itest(i) = 1
             igwd = igwd + 1
             idx(igwd) = i
          ENDIF
       ENDDO

       CALL drag_noro(klon, llm, dtphys, paprs, play, zmea, zstd, zsig, zgam, &
            zthe, zpic, zval, itest, t_seri, u_seri, v_seri, zulow, zvlow, &
            zustrdr, zvstrdr, d_t_oro, d_u_oro, d_v_oro)

       ! ajout des tendances
       DO k = 1, llm
          DO i = 1, klon
             t_seri(i, k) = t_seri(i, k) + d_t_oro(i, k)
             u_seri(i, k) = u_seri(i, k) + d_u_oro(i, k)
             v_seri(i, k) = v_seri(i, k) + d_v_oro(i, k)
          ENDDO
       ENDDO
    ENDIF

    IF (ok_orolf) THEN
       ! S\'election des points pour lesquels le sch\'ema est actif :
       igwd = 0
       DO i = 1, klon
          itest(i) = 0
          IF ((zpic(i) - zmea(i)) > 100.) THEN
             itest(i) = 1
             igwd = igwd + 1
             idx(igwd) = i
          ENDIF
       ENDDO

       CALL lift_noro(klon, llm, dtphys, paprs, play, rlat, zmea, zstd, zpic, &
            itest, t_seri, u_seri, v_seri, zulow, zvlow, zustrli, zvstrli, &
            d_t_lif, d_u_lif, d_v_lif)

       ! Ajout des tendances :
       DO k = 1, llm
          DO i = 1, klon
             t_seri(i, k) = t_seri(i, k) + d_t_lif(i, k)
             u_seri(i, k) = u_seri(i, k) + d_u_lif(i, k)
             v_seri(i, k) = v_seri(i, k) + d_v_lif(i, k)
          ENDDO
       ENDDO
    ENDIF

    ! Stress n\'ecessaires : toute la physique

    DO i = 1, klon
       zustrph(i) = 0.
       zvstrph(i) = 0.
    ENDDO
    DO k = 1, llm
       DO i = 1, klon
          zustrph(i) = zustrph(i) + (u_seri(i, k) - u(i, k)) / dtphys &
               * zmasse(i, k)
          zvstrph(i) = zvstrph(i) + (v_seri(i, k) - v(i, k)) / dtphys &
               * zmasse(i, k)
       ENDDO
    ENDDO

    CALL aaam_bud(ra, rg, romega, rlat, rlon, pphis, zustrdr, zustrli, &
         zustrph, zvstrdr, zvstrli, zvstrph, paprs, u, v, aam, torsfc)

    IF (if_ebil >= 2) CALL diagetpq(airephy, 'after orography', ip_ebil, 2, &
         2, dtphys, t_seri, q_seri, ql_seri, u_seri, v_seri, paprs, d_h_vcol, &
         d_qt, d_ec)

    ! Calcul des tendances traceurs
    call phytrac(itap, lmt_pas, julien, time, firstcal, lafin, dtphys, t, &
         paprs, play, mfu, mfd, pde_u, pen_d, ycoefh, fm_therm, entr_therm, &
         yu1, yv1, ftsol, pctsrf, frac_impa, frac_nucl, da, phi, mp, upwd, &
         dnwd, tr_seri, zmasse, ncid_startphy, nid_ins)

    IF (offline) call phystokenc(dtphys, rlon, rlat, t, mfu, mfd, pen_u, &
         pde_u, pen_d, pde_d, fm_therm, entr_therm, ycoefh, yu1, yv1, ftsol, &
         pctsrf, frac_impa, frac_nucl, pphis, airephy, dtphys, itap)

    ! Calculer le transport de l'eau et de l'energie (diagnostique)
    CALL transp(paprs, zxtsol, t_seri, q_seri, u_seri, v_seri, zphi, ve, vq, &
         ue, uq)

    ! diag. bilKP

    CALL transp_lay(paprs, zxtsol, t_seri, q_seri, u_seri, v_seri, zphi, &
         ve_lay, vq_lay, ue_lay, uq_lay)

    ! Accumuler les variables a stocker dans les fichiers histoire:

    ! conversion Ec -> E thermique
    DO k = 1, llm
       DO i = 1, klon
          ZRCPD = RCPD * (1. + RVTMP2 * q_seri(i, k))
          d_t_ec(i, k) = 0.5 / ZRCPD &
               * (u(i, k)**2 + v(i, k)**2 - u_seri(i, k)**2 - v_seri(i, k)**2)
          t_seri(i, k) = t_seri(i, k) + d_t_ec(i, k)
          d_t_ec(i, k) = d_t_ec(i, k) / dtphys
       END DO
    END DO

    IF (if_ebil >= 1) THEN 
       tit = 'after physic'
       CALL diagetpq(airephy, tit, ip_ebil, 1, 1, dtphys, t_seri, q_seri, &
            ql_seri, u_seri, v_seri, paprs, d_h_vcol, d_qt, d_ec)
       ! Comme les tendances de la physique sont ajoute dans la dynamique, 
       ! on devrait avoir que la variation d'entalpie par la dynamique
       ! est egale a la variation de la physique au pas de temps precedent.
       ! Donc la somme de ces 2 variations devrait etre nulle.
       call diagphy(airephy, tit, ip_ebil, topsw, toplw, solsw, sollw, sens, &
            evap, rain_fall, snow_fall, ztsol, d_h_vcol, d_qt, d_ec)
       d_h_vcol_phy = d_h_vcol
    END IF

    ! SORTIES

    ! prw = eau precipitable
    DO i = 1, klon
       prw(i) = 0.
       DO k = 1, llm
          prw(i) = prw(i) + q_seri(i, k)*zmasse(i, k)
       ENDDO
    ENDDO

    ! Convertir les incrementations en tendances

    DO k = 1, llm
       DO i = 1, klon
          d_u(i, k) = (u_seri(i, k) - u(i, k)) / dtphys
          d_v(i, k) = (v_seri(i, k) - v(i, k)) / dtphys
          d_t(i, k) = (t_seri(i, k) - t(i, k)) / dtphys
          d_qx(i, k, ivap) = (q_seri(i, k) - qx(i, k, ivap)) / dtphys
          d_qx(i, k, iliq) = (ql_seri(i, k) - qx(i, k, iliq)) / dtphys
       ENDDO
    ENDDO

    DO iq = 3, nqmx
       DO k = 1, llm
          DO i = 1, klon
             d_qx(i, k, iq) = (tr_seri(i, k, iq-2) - qx(i, k, iq)) / dtphys
          ENDDO
       ENDDO
    ENDDO

    ! Sauvegarder les valeurs de t et q a la fin de la physique:
    DO k = 1, llm
       DO i = 1, klon
          t_ancien(i, k) = t_seri(i, k)
          q_ancien(i, k) = q_seri(i, k)
       ENDDO
    ENDDO

    call write_histins

    IF (lafin) then
       call NF95_CLOSE(ncid_startphy)
       CALL phyredem(pctsrf, ftsol, ftsoil, tslab, seaice, fqsurf, qsol, &
            fsnow, falbe, fevap, rain_fall, snow_fall, solsw, sollw, dlw, &
            radsol, frugs, agesno, zmea, zstd, zsig, zgam, zthe, zpic, zval, &
            t_ancien, q_ancien, rnebcon, ratqs, clwcon, run_off_lic_0, sig1, &
            w01)
    end IF

    firstcal = .FALSE.

  contains

    subroutine write_histins

      ! From phylmd/write_histins.h, version 1.2 2005/05/25 13:10:09

      ! Ecriture des sorties

      use dimens_m, only: iim, jjm
      USE histsync_m, ONLY: histsync
      USE histwrite_m, ONLY: histwrite

      integer i, itau_w ! pas de temps ecriture
      REAL zx_tmp_2d(iim, jjm + 1), zx_tmp_3d(iim, jjm + 1, llm)

      !--------------------------------------------------

      IF (ok_instan) THEN
         ! Champs 2D:

         itau_w = itau_phy + itap

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, pphis, zx_tmp_2d)
         CALL histwrite(nid_ins, "phis", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, airephy, zx_tmp_2d)
         CALL histwrite(nid_ins, "aire", itau_w, zx_tmp_2d)

         DO i = 1, klon
            zx_tmp_fi2d(i) = paprs(i, 1)
         ENDDO
         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zx_tmp_fi2d, zx_tmp_2d)
         CALL histwrite(nid_ins, "psol", itau_w, zx_tmp_2d)

         DO i = 1, klon
            zx_tmp_fi2d(i) = rain_fall(i) + snow_fall(i)
         ENDDO
         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zx_tmp_fi2d, zx_tmp_2d)
         CALL histwrite(nid_ins, "precip", itau_w, zx_tmp_2d)

         DO i = 1, klon
            zx_tmp_fi2d(i) = rain_lsc(i) + snow_lsc(i)
         ENDDO
         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zx_tmp_fi2d, zx_tmp_2d)
         CALL histwrite(nid_ins, "plul", itau_w, zx_tmp_2d)

         DO i = 1, klon
            zx_tmp_fi2d(i) = rain_con(i) + snow_con(i)
         ENDDO
         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zx_tmp_fi2d, zx_tmp_2d)
         CALL histwrite(nid_ins, "pluc", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zxtsol, zx_tmp_2d)
         CALL histwrite(nid_ins, "tsol", itau_w, zx_tmp_2d)
         !ccIM
         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zt2m, zx_tmp_2d)
         CALL histwrite(nid_ins, "t2m", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zq2m, zx_tmp_2d)
         CALL histwrite(nid_ins, "q2m", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zu10m, zx_tmp_2d)
         CALL histwrite(nid_ins, "u10m", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zv10m, zx_tmp_2d)
         CALL histwrite(nid_ins, "v10m", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, snow_fall, zx_tmp_2d)
         CALL histwrite(nid_ins, "snow", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, cdragm, zx_tmp_2d)
         CALL histwrite(nid_ins, "cdrm", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, cdragh, zx_tmp_2d)
         CALL histwrite(nid_ins, "cdrh", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, toplw, zx_tmp_2d)
         CALL histwrite(nid_ins, "topl", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, evap, zx_tmp_2d)
         CALL histwrite(nid_ins, "evap", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, solsw, zx_tmp_2d)
         CALL histwrite(nid_ins, "sols", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, sollw, zx_tmp_2d)
         CALL histwrite(nid_ins, "soll", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, sollwdown, zx_tmp_2d)
         CALL histwrite(nid_ins, "solldown", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, bils, zx_tmp_2d)
         CALL histwrite(nid_ins, "bils", itau_w, zx_tmp_2d)

         zx_tmp_fi2d(1:klon) = -1*sens(1:klon)
         ! CALL gr_fi_ecrit(1, klon, iim, jjm + 1, sens, zx_tmp_2d)
         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zx_tmp_fi2d, zx_tmp_2d)
         CALL histwrite(nid_ins, "sens", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, fder, zx_tmp_2d)
         CALL histwrite(nid_ins, "fder", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, d_ts(1, is_oce), zx_tmp_2d)
         CALL histwrite(nid_ins, "dtsvdfo", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, d_ts(1, is_ter), zx_tmp_2d)
         CALL histwrite(nid_ins, "dtsvdft", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, d_ts(1, is_lic), zx_tmp_2d)
         CALL histwrite(nid_ins, "dtsvdfg", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, d_ts(1, is_sic), zx_tmp_2d)
         CALL histwrite(nid_ins, "dtsvdfi", itau_w, zx_tmp_2d)

         DO nsrf = 1, nbsrf
            !XXX
            zx_tmp_fi2d(1 : klon) = pctsrf(1 : klon, nsrf)*100.
            CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zx_tmp_fi2d, zx_tmp_2d)
            CALL histwrite(nid_ins, "pourc_"//clnsurf(nsrf), itau_w, &
                 zx_tmp_2d) 

            zx_tmp_fi2d(1 : klon) = pctsrf(1 : klon, nsrf)
            CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zx_tmp_fi2d, zx_tmp_2d)
            CALL histwrite(nid_ins, "fract_"//clnsurf(nsrf), itau_w, &
                 zx_tmp_2d) 

            zx_tmp_fi2d(1 : klon) = fluxt(1 : klon, 1, nsrf)
            CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zx_tmp_fi2d, zx_tmp_2d)
            CALL histwrite(nid_ins, "sens_"//clnsurf(nsrf), itau_w, &
                 zx_tmp_2d) 

            zx_tmp_fi2d(1 : klon) = fluxlat(1 : klon, nsrf)
            CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zx_tmp_fi2d, zx_tmp_2d)
            CALL histwrite(nid_ins, "lat_"//clnsurf(nsrf), itau_w, &
                 zx_tmp_2d) 

            zx_tmp_fi2d(1 : klon) = ftsol(1 : klon, nsrf)
            CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zx_tmp_fi2d, zx_tmp_2d)
            CALL histwrite(nid_ins, "tsol_"//clnsurf(nsrf), itau_w, &
                 zx_tmp_2d) 

            zx_tmp_fi2d(1 : klon) = fluxu(1 : klon, 1, nsrf)
            CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zx_tmp_fi2d, zx_tmp_2d)
            CALL histwrite(nid_ins, "taux_"//clnsurf(nsrf), itau_w, &
                 zx_tmp_2d) 

            zx_tmp_fi2d(1 : klon) = fluxv(1 : klon, 1, nsrf)
            CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zx_tmp_fi2d, zx_tmp_2d)
            CALL histwrite(nid_ins, "tauy_"//clnsurf(nsrf), itau_w, &
                 zx_tmp_2d)

            zx_tmp_fi2d(1 : klon) = frugs(1 : klon, nsrf)
            CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zx_tmp_fi2d, zx_tmp_2d)
            CALL histwrite(nid_ins, "rugs_"//clnsurf(nsrf), itau_w, &
                 zx_tmp_2d) 

            zx_tmp_fi2d(1 : klon) = falbe(:, nsrf)
            CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zx_tmp_fi2d, zx_tmp_2d)
            CALL histwrite(nid_ins, "albe_"//clnsurf(nsrf), itau_w, &
                 zx_tmp_2d) 

         END DO
         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, albsol, zx_tmp_2d)
         CALL histwrite(nid_ins, "albs", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, zxrugs, zx_tmp_2d)
         CALL histwrite(nid_ins, "rugs", itau_w, zx_tmp_2d)

         !HBTM2

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, s_pblh, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_pblh", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, s_pblt, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_pblt", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, s_lcl, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_lcl", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, s_capCL, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_capCL", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, s_oliqCL, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_oliqCL", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, s_cteiCL, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_cteiCL", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, s_therm, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_therm", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, s_trmb1, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_trmb1", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, s_trmb2, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_trmb2", itau_w, zx_tmp_2d)

         CALL gr_fi_ecrit(1, klon, iim, jjm + 1, s_trmb3, zx_tmp_2d)
         CALL histwrite(nid_ins, "s_trmb3", itau_w, zx_tmp_2d)

         ! Champs 3D:

         CALL gr_fi_ecrit(llm, klon, iim, jjm + 1, t_seri, zx_tmp_3d)
         CALL histwrite(nid_ins, "temp", itau_w, zx_tmp_3d)

         CALL gr_fi_ecrit(llm, klon, iim, jjm + 1, u_seri, zx_tmp_3d)
         CALL histwrite(nid_ins, "vitu", itau_w, zx_tmp_3d)

         CALL gr_fi_ecrit(llm, klon, iim, jjm + 1, v_seri, zx_tmp_3d)
         CALL histwrite(nid_ins, "vitv", itau_w, zx_tmp_3d)

         CALL gr_fi_ecrit(llm, klon, iim, jjm + 1, zphi, zx_tmp_3d)
         CALL histwrite(nid_ins, "geop", itau_w, zx_tmp_3d)

         CALL gr_fi_ecrit(llm, klon, iim, jjm + 1, play, zx_tmp_3d)
         CALL histwrite(nid_ins, "pres", itau_w, zx_tmp_3d)

         CALL gr_fi_ecrit(llm, klon, iim, jjm + 1, d_t_vdf, zx_tmp_3d)
         CALL histwrite(nid_ins, "dtvdf", itau_w, zx_tmp_3d)

         CALL gr_fi_ecrit(llm, klon, iim, jjm + 1, d_q_vdf, zx_tmp_3d)
         CALL histwrite(nid_ins, "dqvdf", itau_w, zx_tmp_3d)

         call histsync(nid_ins)
      ENDIF

    end subroutine write_histins

  END SUBROUTINE physiq

end module physiq_m
