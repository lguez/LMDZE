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
    use ajsec_m, only: ajsec
    use calltherm_m, only: calltherm
    USE clesphys, ONLY: cdhmax, cdmmax, ecrit_ins, ksta, ksta_ter, ok_kzmin, &
         ok_instan
    USE clesphys2, ONLY: conv_emanuel, nbapp_rad, new_oliq, ok_orodr, ok_orolf
    USE clmain_m, ONLY: clmain
    use clouds_gno_m, only: clouds_gno
    use comconst, only: dtphys
    USE comgeomphy, ONLY: airephy
    USE concvl_m, ONLY: concvl
    USE conf_gcm_m, ONLY: lmt_pas
    USE conf_phys_m, ONLY: conf_phys
    use conflx_m, only: conflx
    USE ctherm, ONLY: iflag_thermals, nsplit_thermals
    use diagcld2_m, only: diagcld2
    USE dimens_m, ONLY: llm, nqmx
    USE dimphy, ONLY: klon
    USE dimsoil, ONLY: nsoilmx
    use drag_noro_m, only: drag_noro
    use dynetat0_m, only: day_ref, annee_ref
    USE fcttre, ONLY: foeew
    use fisrtilp_m, only: fisrtilp
    USE hgardfou_m, ONLY: hgardfou
    USE histsync_m, ONLY: histsync
    USE histwrite_phy_m, ONLY: histwrite_phy
    USE indicesol, ONLY: clnsurf, epsfra, is_lic, is_oce, is_sic, is_ter, &
         nbsrf
    USE ini_histins_m, ONLY: ini_histins, nid_ins
    use lift_noro_m, only: lift_noro
    use netcdf95, only: NF95_CLOSE
    use newmicro_m, only: newmicro
    use nr_util, only: assert
    use nuage_m, only: nuage
    USE orbite_m, ONLY: orbite
    USE ozonecm_m, ONLY: ozonecm
    USE phyetat0_m, ONLY: phyetat0
    USE phyredem_m, ONLY: phyredem
    USE phyredem0_m, ONLY: phyredem0
    USE phytrac_m, ONLY: phytrac
    use radlwsw_m, only: radlwsw
    use yoegwd, only: sugwd
    USE suphec_m, ONLY: rcpd, retv, rg, rlvtt, romega, rsigma, rtt, rmo3, md
    use time_phylmdz, only: itap, increment_itap
    use transp_m, only: transp
    use transp_lay_m, only: transp_lay
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
    ! vitesse dans la direction X (de O a E) en m / s

    REAL, intent(in):: v(:, :) ! (klon, llm) vitesse Y (de S a N) en m / s
    REAL, intent(in):: t(:, :) ! (klon, llm) temperature (K)

    REAL, intent(in):: qx(:, :, :) ! (klon, llm, nqmx)
    ! (humidit\'e sp\'ecifique et fractions massiques des autres traceurs)

    REAL, intent(in):: omega(:, :) ! (klon, llm) vitesse verticale en Pa / s
    REAL, intent(out):: d_u(:, :) ! (klon, llm) tendance physique de "u" (m s-2)
    REAL, intent(out):: d_v(:, :) ! (klon, llm) tendance physique de "v" (m s-2)
    REAL, intent(out):: d_t(:, :) ! (klon, llm) tendance physique de "t" (K / s)

    REAL, intent(out):: d_qx(:, :, :) ! (klon, llm, nqmx)
    ! tendance physique de "qx" (s-1)

    ! Local:

    LOGICAL:: firstcal = .true.

    LOGICAL, PARAMETER:: ok_stratus = .FALSE.
    ! Ajouter artificiellement les stratus

    ! pour phystoke avec thermiques
    REAL fm_therm(klon, llm + 1)
    REAL entr_therm(klon, llm)
    real, save:: q2(klon, llm + 1, nbsrf)

    INTEGER, PARAMETER:: ivap = 1 ! indice de traceur pour vapeur d'eau
    INTEGER, PARAMETER:: iliq = 2 ! indice de traceur pour eau liquide

    REAL, save:: t_ancien(klon, llm), q_ancien(klon, llm)
    LOGICAL, save:: ancien_ok

    REAL d_t_dyn(klon, llm) ! tendance dynamique pour "t" (K / s)
    REAL d_q_dyn(klon, llm) ! tendance dynamique pour "q" (kg / kg / s)

    real da(klon, llm), phi(klon, llm, llm), mp(klon, llm)

    REAL, save:: swdn0(klon, llm + 1), swdn(klon, llm + 1)
    REAL, save:: swup0(klon, llm + 1), swup(klon, llm + 1)

    REAL, save:: lwdn0(klon, llm + 1), lwdn(klon, llm + 1)
    REAL, save:: lwup0(klon, llm + 1), lwup(klon, llm + 1)

    ! prw: precipitable water
    real prw(klon)

    ! flwp, fiwp = Liquid Water Path & Ice Water Path (kg / m2)
    ! flwc, fiwc = Liquid Water Content & Ice Water Content (kg / kg)
    REAL flwp(klon), fiwp(klon)
    REAL flwc(klon, llm), fiwc(klon, llm)

    ! Variables propres a la physique

    INTEGER, save:: radpas
    ! Radiative transfer computations are made every "radpas" call to
    ! "physiq".

    REAL, save:: radsol(klon) ! bilan radiatif au sol calcule par code radiatif
    REAL, save:: ftsol(klon, nbsrf) ! skin temperature of surface fraction

    REAL, save:: ftsoil(klon, nsoilmx, nbsrf)
    ! soil temperature of surface fraction

    REAL, save:: fevap(klon, nbsrf) ! evaporation
    REAL fluxlat(klon, nbsrf)

    REAL, save:: fqsurf(klon, nbsrf)
    ! humidite de l'air au contact de la surface

    REAL, save:: qsol(klon) ! column-density of water in soil, in kg m-2
    REAL, save:: fsnow(klon, nbsrf) ! \'epaisseur neigeuse
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
    INTEGER igwd, itest(klon)

    REAL, save:: agesno(klon, nbsrf) ! age de la neige
    REAL, save:: run_off_lic_0(klon)

    ! Variables li\'ees \`a la convection d'Emanuel :
    REAL, save:: Ma(klon, llm) ! undilute upward mass flux
    REAL, save:: qcondc(klon, llm) ! in-cld water content from convect
    REAL, save:: sig1(klon, llm), w01(klon, llm)

    ! Variables pour la couche limite (Alain Lahellec) :
    REAL cdragh(klon) ! drag coefficient pour T and Q
    REAL cdragm(klon) ! drag coefficient pour vent

    REAL ycoefh(klon, llm) ! coef d'echange pour phytrac

    REAL, save:: ffonte(klon, nbsrf)
    ! flux thermique utilise pour fondre la neige

    REAL, save:: fqcalving(klon, nbsrf)
    ! flux d'eau "perdue" par la surface et necessaire pour limiter la
    ! hauteur de neige, en kg / m2 / s

    REAL zxffonte(klon), zxfqcalving(klon)

    REAL, save:: pfrac_impa(klon, llm)! Produits des coefs lessivage impaction
    REAL, save:: pfrac_nucl(klon, llm)! Produits des coefs lessivage nucleation

    REAL, save:: pfrac_1nucl(klon, llm)
    ! Produits des coefs lessi nucl (alpha = 1)

    REAL frac_impa(klon, llm) ! fraction d'a\'erosols lessiv\'es (impaction)
    REAL frac_nucl(klon, llm) ! idem (nucleation)

    REAL, save:: rain_fall(klon)
    ! liquid water mass flux (kg / m2 / s), positive down

    REAL, save:: snow_fall(klon)
    ! solid water mass flux (kg / m2 / s), positive down

    REAL rain_tiedtke(klon), snow_tiedtke(klon)

    REAL evap(klon) ! flux d'\'evaporation au sol
    real devap(klon) ! derivative of the evaporation flux at the surface
    REAL sens(klon) ! flux de chaleur sensible au sol
    real dsens(klon) ! derivee du flux de chaleur sensible au sol
    REAL, save:: dlw(klon) ! derivative of infra-red flux
    REAL bils(klon) ! bilan de chaleur au sol
    REAL fder(klon) ! Derive de flux (sensible et latente)
    REAL ve(klon) ! integr. verticale du transport meri. de l'energie
    REAL vq(klon) ! integr. verticale du transport meri. de l'eau
    REAL ue(klon) ! integr. verticale du transport zonal de l'energie
    REAL uq(klon) ! integr. verticale du transport zonal de l'eau

    REAL, save:: frugs(klon, nbsrf) ! longueur de rugosite
    REAL zxrugs(klon) ! longueur de rugosite

    ! Conditions aux limites

    INTEGER julien
    REAL, save:: pctsrf(klon, nbsrf) ! percentage of surface
    REAL, save:: albsol(klon) ! albedo du sol total, visible, moyen par maille
    REAL, SAVE:: wo(klon, llm) ! column density of ozone in a cell, in kDU
    real, parameter:: dobson_u = 2.1415e-05 ! Dobson unit, in kg m-2

    real, save:: clwcon(klon, llm), rnebcon(klon, llm)
    real, save:: clwcon0(klon, llm), rnebcon0(klon, llm)

    REAL rhcl(klon, llm) ! humiditi relative ciel clair
    REAL dialiq(klon, llm) ! eau liquide nuageuse
    REAL diafra(klon, llm) ! fraction nuageuse
    REAL cldliq(klon, llm) ! eau liquide nuageuse
    REAL cldfra(klon, llm) ! fraction nuageuse
    REAL cldtau(klon, llm) ! epaisseur optique
    REAL cldemi(klon, llm) ! emissivite infrarouge

    REAL flux_q(klon, nbsrf) ! flux turbulent d'humidite à la surface
    REAL flux_t(klon, nbsrf) ! flux turbulent de chaleur à la surface

    REAL flux_u(klon, nbsrf), flux_v(klon, nbsrf)
    ! tension du vent (flux turbulent de vent) à la surface, en Pa

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
    REAL fsollw(klon, nbsrf) ! bilan flux IR pour chaque sous-surface
    REAL fsolsw(klon, nbsrf) ! flux solaire absorb\'e pour chaque sous-surface

    REAL conv_q(klon, llm) ! convergence de l'humidite (kg / kg / s)
    REAL conv_t(klon, llm) ! convergence of temperature (K / s)

    REAL cldl(klon), cldm(klon), cldh(klon) ! nuages bas, moyen et haut
    REAL cldt(klon), cldq(klon) ! nuage total, eau liquide integree

    REAL zxfluxlat(klon)
    REAL dist, mu0(klon), fract(klon)
    real longi
    REAL z_avant(klon), z_apres(klon), z_factor(klon)
    REAL zb
    REAL zx_t, zx_qs, zcor
    real zqsat(klon, llm)
    INTEGER i, k, iq, nsrf
    REAL zphi(klon, llm)

    ! cf. Anne Mathieu, variables pour la couche limite atmosphérique (hbtm)

    REAL, SAVE:: pblh(klon, nbsrf) ! Hauteur de couche limite
    REAL, SAVE:: plcl(klon, nbsrf) ! Niveau de condensation de la CLA
    REAL, SAVE:: capCL(klon, nbsrf) ! CAPE de couche limite
    REAL, SAVE:: oliqCL(klon, nbsrf) ! eau_liqu integree de couche limite
    REAL, SAVE:: cteiCL(klon, nbsrf) ! cloud top instab. crit. couche limite
    REAL, SAVE:: pblt(klon, nbsrf) ! T \`a la hauteur de couche limite
    REAL, SAVE:: therm(klon, nbsrf)
    REAL, SAVE:: trmb1(klon, nbsrf) ! deep_cape
    REAL, SAVE:: trmb2(klon, nbsrf) ! inhibition
    REAL, SAVE:: trmb3(klon, nbsrf) ! Point Omega
    ! Grandeurs de sorties
    REAL s_pblh(klon), s_lcl(klon), s_capCL(klon)
    REAL s_oliqCL(klon), s_cteiCL(klon), s_pblt(klon)
    REAL s_therm(klon), s_trmb1(klon), s_trmb2(klon)
    REAL s_trmb3(klon)

    ! Variables pour la convection de K. Emanuel :

    REAL upwd(klon, llm) ! saturated updraft mass flux
    REAL dnwd(klon, llm) ! saturated downdraft mass flux
    REAL, save:: cape(klon)

    INTEGER iflagctrl(klon) ! flag fonctionnement de convect

    ! Variables du changement

    ! con: convection
    ! lsc: large scale condensation
    ! ajs: ajustement sec
    ! eva: \'evaporation de l'eau liquide nuageuse
    ! vdf: vertical diffusion in boundary layer
    REAL d_t_con(klon, llm), d_q_con(klon, llm)
    REAL, save:: d_u_con(klon, llm), d_v_con(klon, llm)
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
    real ema_pct(klon) ! Emanuel pressure at cloud top, in Pa

    REAL, save:: rain_con(klon)
    real rain_lsc(klon)
    REAL, save:: snow_con(klon) ! neige (mm / s)
    real snow_lsc(klon)
    REAL d_ts(klon, nbsrf) ! variation of ftsol

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

    ! Variables pour effectuer les appels en s\'erie :

    REAL t_seri(klon, llm), q_seri(klon, llm)
    REAL ql_seri(klon, llm)
    REAL u_seri(klon, llm), v_seri(klon, llm)
    REAL tr_seri(klon, llm, nqmx - 2)

    REAL zx_rh(klon, llm)

    REAL zustrdr(klon), zvstrdr(klon)
    REAL zustrli(klon), zvstrli(klon)
    REAL aam, torsfc

    REAL ve_lay(klon, llm) ! transport meri. de l'energie a chaque niveau vert.
    REAL vq_lay(klon, llm) ! transport meri. de l'eau a chaque niveau vert.
    REAL ue_lay(klon, llm) ! transport zonal de l'energie a chaque niveau vert.
    REAL uq_lay(klon, llm) ! transport zonal de l'eau a chaque niveau vert.

    real date0
    REAL tsol(klon)

    REAL d_t_ec(klon, llm)
    ! tendance due \`a la conversion d'\'energie cin\'etique en
    ! énergie thermique

    REAL, save:: t2m(klon, nbsrf), q2m(klon, nbsrf) 
    ! temperature and humidity at 2 m

    REAL, save:: u10m_srf(klon, nbsrf), v10m_srf(klon, nbsrf)
    ! composantes du vent \`a 10 m
    
    REAL zt2m(klon), zq2m(klon) ! température, humidité 2 m moyenne sur 1 maille
    REAL u10m(klon), v10m(klon) ! vent \`a 10 m moyenn\' sur les sous-surfaces

    ! Aerosol effects:

    REAL, save:: topswad(klon), solswad(klon) ! aerosol direct effect
    LOGICAL:: ok_ade = .false. ! apply aerosol direct effect

    REAL:: bl95_b0 = 2., bl95_b1 = 0.2
    ! Parameters in equation (D) of Boucher and Lohmann (1995, Tellus
    ! B). They link cloud droplet number concentration to aerosol mass
    ! concentration.

    real zmasse(klon, llm)
    ! (column-density of mass of air in a cell, in kg m-2)

    integer, save:: ncid_startphy

    namelist /physiq_nml/ fact_cldcon, facttemps, ok_newmicro, iflag_cldcon, &
         ratqsbas, ratqshaut, ok_ade, bl95_b0, bl95_b1, iflag_thermals, &
         nsplit_thermals

    !----------------------------------------------------------------

    IF (nqmx < 2) CALL abort_gcm('physiq', &
         'eaux vapeur et liquide sont indispensables')

    test_firstcal: IF (firstcal) THEN
       ! initialiser
       u10m_srf = 0.
       v10m_srf = 0.
       t2m = 0.
       q2m = 0.
       ffonte = 0.
       fqcalving = 0.
       rain_con = 0.
       snow_con = 0.
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
       pblt =0.
       therm =0.
       trmb1 =0. ! deep_cape
       trmb2 =0. ! inhibition
       trmb3 =0. ! Point Omega

       iflag_thermals = 0
       nsplit_thermals = 1
       print *, "Enter namelist 'physiq_nml'."
       read(unit=*, nml=physiq_nml)
       write(unit_nml, nml=physiq_nml)

       call conf_phys

       ! Initialiser les compteurs:

       frugs = 0.
       CALL phyetat0(pctsrf, ftsol, ftsoil, fqsurf, qsol, fsnow, falbe, &
            fevap, rain_fall, snow_fall, solsw, sollw, dlw, radsol, frugs, &
            agesno, zmea, zstd, zsig, zgam, zthe, zpic, zval, t_ancien, &
            q_ancien, ancien_ok, rnebcon, ratqs, clwcon, run_off_lic_0, sig1, &
            w01, ncid_startphy)

       ! ATTENTION : il faudra a terme relire q2 dans l'etat initial
       q2 = 1e-8

       radpas = lmt_pas / nbapp_rad
       print *, "radpas = ", radpas

       ! Initialisation pour le sch\'ema de convection d'Emanuel :
       IF (conv_emanuel) THEN
          ibas_con = 1
          itop_con = 1
       ENDIF

       IF (ok_orodr) THEN
          rugoro = MAX(1e-5, zstd * zsig / 2)
          CALL SUGWD(paprs, play)
       else
          rugoro = 0.
       ENDIF

       ecrit_ins = NINT(ecrit_ins / dtphys)

       ! Initialisation des sorties

       call ini_histins(dtphys, ok_newmicro)
       CALL ymds2ju(annee_ref, 1, day_ref, 0., date0)
       ! Positionner date0 pour initialisation de ORCHIDEE
       print *, 'physiq date0: ', date0
       CALL phyredem0
    ENDIF test_firstcal

    ! We will modify variables *_seri and we will not touch variables
    ! u, v, t, qx:
    t_seri = t
    u_seri = u
    v_seri = v
    q_seri = qx(:, :, ivap)
    ql_seri = qx(:, :, iliq)
    tr_seri = qx(:, :, 3:nqmx)

    tsol = sum(ftsol * pctsrf, dim = 2)

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

    call increment_itap
    julien = MOD(dayvrai, 360)
    if (julien == 0) julien = 360

    forall (k = 1: llm) zmasse(:, k) = (paprs(:, k) - paprs(:, k + 1)) / rg

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

    frugs = MAX(frugs, 0.000015)
    zxrugs = sum(frugs * pctsrf, dim = 2)

    ! Calculs n\'ecessaires au calcul de l'albedo dans l'interface avec
    ! la surface.

    CALL orbite(REAL(julien), longi, dist)
    CALL zenang(longi, time, dtphys * radpas, mu0, fract)
    albsol = sum(falbe * pctsrf, dim = 2)

    ! R\'epartition sous maille des flux longwave et shortwave
    ! R\'epartition du longwave par sous-surface lin\'earis\'ee

    forall (nsrf = 1: nbsrf)
       fsollw(:, nsrf) = sollw + 4. * RSIGMA * tsol**3 &
            * (tsol - ftsol(:, nsrf))
       fsolsw(:, nsrf) = solsw * (1. - falbe(:, nsrf)) / (1. - albsol)
    END forall

    CALL clmain(dtphys, pctsrf, t_seri, q_seri, u_seri, v_seri, julien, mu0, &
         ftsol, cdmmax, cdhmax, ksta, ksta_ter, ok_kzmin, ftsoil, qsol, &
         paprs, play, fsnow, fqsurf, fevap, falbe, fluxlat, rain_fall, &
         snow_fall, fsolsw, fsollw, frugs, agesno, rugoro, d_t_vdf, d_q_vdf, &
         d_u_vdf, d_v_vdf, d_ts, flux_t, flux_q, flux_u, flux_v, cdragh, &
         cdragm, q2, dsens, devap, ycoefh, t2m, q2m, u10m_srf, v10m_srf, &
         pblh, capCL, oliqCL, cteiCL, pblT, therm, trmb1, trmb2, trmb3, plcl, &
         fqcalving, ffonte, run_off_lic_0)

    ! Incr\'ementation des flux

    sens = - sum(flux_t * pctsrf, dim = 2)
    evap = - sum(flux_q * pctsrf, dim = 2)
    fder = dlw + dsens + devap

    DO k = 1, llm
       DO i = 1, klon
          t_seri(i, k) = t_seri(i, k) + d_t_vdf(i, k)
          q_seri(i, k) = q_seri(i, k) + d_q_vdf(i, k)
          u_seri(i, k) = u_seri(i, k) + d_u_vdf(i, k)
          v_seri(i, k) = v_seri(i, k) + d_v_vdf(i, k)
       ENDDO
    ENDDO

    ! Update surface temperature:

    call assert(abs(sum(pctsrf, dim = 2) - 1.) <= EPSFRA, 'physiq: pctsrf')
    ftsol = ftsol + d_ts
    tsol = sum(ftsol * pctsrf, dim = 2)
    zxfluxlat = sum(fluxlat * pctsrf, dim = 2)
    zt2m = sum(t2m * pctsrf, dim = 2)
    zq2m = sum(q2m * pctsrf, dim = 2)
    u10m = sum(u10m_srf * pctsrf, dim = 2)
    v10m = sum(v10m_srf * pctsrf, dim = 2)
    zxffonte = sum(ffonte * pctsrf, dim = 2)
    zxfqcalving = sum(fqcalving * pctsrf, dim = 2)
    s_pblh = sum(pblh * pctsrf, dim = 2)
    s_lcl = sum(plcl * pctsrf, dim = 2)
    s_capCL = sum(capCL * pctsrf, dim = 2)
    s_oliqCL = sum(oliqCL * pctsrf, dim = 2)
    s_cteiCL = sum(cteiCL * pctsrf, dim = 2)
    s_pblT = sum(pblT * pctsrf, dim = 2)
    s_therm = sum(therm * pctsrf, dim = 2)
    s_trmb1 = sum(trmb1 * pctsrf, dim = 2)
    s_trmb2 = sum(trmb2 * pctsrf, dim = 2)
    s_trmb3 = sum(trmb3 * pctsrf, dim = 2)

    ! Si une sous-fraction n'existe pas, elle prend la valeur moyenne :
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          IF (pctsrf(i, nsrf) < epsfra) then
             ftsol(i, nsrf) = tsol(i)
             t2m(i, nsrf) = zt2m(i)
             q2m(i, nsrf) = zq2m(i)
             u10m_srf(i, nsrf) = u10m(i)
             v10m_srf(i, nsrf) = v10m(i)
             ffonte(i, nsrf) = zxffonte(i)
             fqcalving(i, nsrf) = zxfqcalving(i)
             pblh(i, nsrf) = s_pblh(i)
             plcl(i, nsrf) = s_lcl(i)
             capCL(i, nsrf) = s_capCL(i)
             oliqCL(i, nsrf) = s_oliqCL(i)
             cteiCL(i, nsrf) = s_cteiCL(i)
             pblT(i, nsrf) = s_pblT(i)
             therm(i, nsrf) = s_therm(i)
             trmb1(i, nsrf) = s_trmb1(i)
             trmb2(i, nsrf) = s_trmb2(i)
             trmb3(i, nsrf) = s_trmb3(i)
          end IF
       ENDDO
    ENDDO

    dlw = - 4. * RSIGMA * tsol**3

    ! Appeler la convection

    if (conv_emanuel) then
       CALL concvl(paprs, play, t_seri, q_seri, u_seri, v_seri, sig1, w01, &
            d_t_con, d_q_con, d_u_con, d_v_con, rain_con, ibas_con, itop_con, &
            upwd, dnwd, Ma, cape, iflagctrl, qcondc, pmflxr, da, phi, mp)
       snow_con = 0.
       clwcon0 = qcondc
       mfu = upwd + dnwd

       zqsat = MIN(0.5, r2es * FOEEW(t_seri, rtt >= t_seri) / play)
       zqsat = zqsat / (1. - retv * zqsat)

       ! Properties of convective clouds
       clwcon0 = fact_cldcon * clwcon0
       call clouds_gno(klon, llm, q_seri, zqsat, clwcon0, ptconv, ratqsc, &
            rnebcon0)

       forall (i = 1:klon) ema_pct(i) = paprs(i, itop_con(i) + 1)
       mfd = 0.
       pen_u = 0.
       pen_d = 0.
       pde_d = 0.
       pde_u = 0.
    else
       conv_q = d_q_dyn + d_q_vdf / dtphys
       conv_t = d_t_dyn + d_t_vdf / dtphys
       z_avant = sum((q_seri + ql_seri) * zmasse, dim=2)
       CALL conflx(dtphys, paprs, play, t_seri(:, llm:1:- 1), &
            q_seri(:, llm:1:- 1), conv_t, conv_q, - evap, omega, &
            d_t_con, d_q_con, rain_con, snow_con, mfu(:, llm:1:- 1), &
            mfd(:, llm:1:- 1), pen_u, pde_u, pen_d, pde_d, kcbot, kctop, &
            kdtop, pmflxr, pmflxs)
       WHERE (rain_con < 0.) rain_con = 0.
       WHERE (snow_con < 0.) snow_con = 0.
       ibas_con = llm + 1 - kcbot
       itop_con = llm + 1 - kctop
    END if

    DO k = 1, llm
       DO i = 1, klon
          t_seri(i, k) = t_seri(i, k) + d_t_con(i, k)
          q_seri(i, k) = q_seri(i, k) + d_q_con(i, k)
          u_seri(i, k) = u_seri(i, k) + d_u_con(i, k)
          v_seri(i, k) = v_seri(i, k) + d_v_con(i, k)
       ENDDO
    ENDDO

    IF (.not. conv_emanuel) THEN
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
       call calltherm(dtphys, play, paprs, pphi, u_seri, v_seri, t_seri, &
            q_seri, d_u_ajs, d_v_ajs, d_t_ajs, d_q_ajs, fm_therm, entr_therm)
    endif

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

    ! PRESCRIPTION DES NUAGES POUR LE RAYONNEMENT

    ! 1. NUAGES CONVECTIFS

    IF (iflag_cldcon <= - 1) THEN
       ! seulement pour Tiedtke
       snow_tiedtke = 0.
       if (iflag_cldcon == - 1) then
          rain_tiedtke = rain_con
       else
          rain_tiedtke = 0.
          do k = 1, llm
             do i = 1, klon
                if (d_q_con(i, k) < 0.) then
                   rain_tiedtke(i) = rain_tiedtke(i) - d_q_con(i, k) / dtphys &
                        * zmasse(i, k)
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
       cldliq = cldliq + rnebcon * clwcon
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

    ! Humidit\'e relative pour diagnostic :
    DO k = 1, llm
       DO i = 1, klon
          zx_t = t_seri(i, k)
          zx_qs = r2es * FOEEW(zx_t, rtt >= zx_t) / play(i, k)
          zx_qs = MIN(0.5, zx_qs)
          zcor = 1. / (1. - retv * zx_qs)
          zx_qs = zx_qs * zcor
          zx_rh(i, k) = q_seri(i, k) / zx_qs
          zqsat(i, k) = zx_qs
       ENDDO
    ENDDO

    ! Param\`etres optiques des nuages et quelques param\`etres pour
    ! diagnostics :
    if (ok_newmicro) then
       CALL newmicro(paprs, play, t_seri, cldliq, cldfra, cldtau, cldemi, &
            cldh, cldl, cldm, cldt, cldq, flwp, fiwp, flwc, fiwc)
    else
       CALL nuage(paprs, play, t_seri, cldliq, cldfra, cldtau, cldemi, cldh, &
            cldl, cldm, cldt, cldq)
    endif

    IF (MOD(itap - 1, radpas) == 0) THEN
       wo = ozonecm(REAL(julien), paprs)
       albsol = sum(falbe * pctsrf, dim = 2)
       CALL radlwsw(dist, mu0, fract, paprs, play, tsol, albsol, t_seri, &
            q_seri, wo, cldfra, cldemi, cldtau, heat, heat0, cool, cool0, &
            radsol, albpla, topsw, toplw, solsw, sollw, sollwdown, topsw0, &
            toplw0, solsw0, sollw0, lwdn0, lwdn, lwup0, lwup, swdn0, swdn, &
            swup0, swup, ok_ade, topswad, solswad)
    ENDIF

    ! Ajouter la tendance des rayonnements (tous les pas)
    DO k = 1, llm
       DO i = 1, klon
          t_seri(i, k) = t_seri(i, k) + (heat(i, k) - cool(i, k)) * dtphys &
               / 86400.
       ENDDO
    ENDDO

    ! Calculer le bilan du sol et la d\'erive de temp\'erature (couplage)
    DO i = 1, klon
       bils(i) = radsol(i) - sens(i) + zxfluxlat(i)
    ENDDO

    ! Param\'etrisation de l'orographie \`a l'\'echelle sous-maille :

    IF (ok_orodr) THEN
       ! S\'election des points pour lesquels le sch\'ema est actif :
       igwd = 0
       DO i = 1, klon
          itest(i) = 0
          IF (zpic(i) - zmea(i) > 100. .AND. zstd(i) > 10.) THEN
             itest(i) = 1
             igwd = igwd + 1
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
          IF (zpic(i) - zmea(i) > 100.) THEN
             itest(i) = 1
             igwd = igwd + 1
          ENDIF
       ENDDO

       CALL lift_noro(dtphys, paprs, play, zmea, zstd, zpic, itest, t_seri, &
            u_seri, v_seri, zulow, zvlow, zustrli, zvstrli, d_t_lif, &
            d_u_lif, d_v_lif)

       ! Ajout des tendances :
       DO k = 1, llm
          DO i = 1, klon
             t_seri(i, k) = t_seri(i, k) + d_t_lif(i, k)
             u_seri(i, k) = u_seri(i, k) + d_u_lif(i, k)
             v_seri(i, k) = v_seri(i, k) + d_v_lif(i, k)
          ENDDO
       ENDDO
    ENDIF

    CALL aaam_bud(rg, romega, pphis, zustrdr, zustrli, &
         sum((u_seri - u) / dtphys * zmasse, dim = 2), zvstrdr, &
         zvstrli, sum((v_seri - v) / dtphys * zmasse, dim = 2), paprs, u, v, &
         aam, torsfc)

    ! Calcul des tendances traceurs
    call phytrac(julien, time, firstcal, lafin, dtphys, t, paprs, play, mfu, &
         mfd, pde_u, pen_d, ycoefh, fm_therm, entr_therm, u(:, 1), v(:, 1), &
         ftsol, pctsrf, frac_impa, frac_nucl, da, phi, mp, upwd, dnwd, &
         tr_seri, zmasse, ncid_startphy)

    ! Calculer le transport de l'eau et de l'energie (diagnostique)
    CALL transp(paprs, t_seri, q_seri, u_seri, v_seri, zphi, ve, vq, ue, uq)

    ! diag. bilKP

    CALL transp_lay(paprs, t_seri, q_seri, u_seri, v_seri, zphi, &
         ve_lay, vq_lay, ue_lay, uq_lay)

    ! Accumuler les variables a stocker dans les fichiers histoire:

    ! conversion Ec en énergie thermique
    DO k = 1, llm
       DO i = 1, klon
          d_t_ec(i, k) = 0.5 / (RCPD * (1. + RVTMP2 * q_seri(i, k))) &
               * (u(i, k)**2 + v(i, k)**2 - u_seri(i, k)**2 - v_seri(i, k)**2)
          t_seri(i, k) = t_seri(i, k) + d_t_ec(i, k)
          d_t_ec(i, k) = d_t_ec(i, k) / dtphys
       END DO
    END DO

    ! SORTIES

    ! prw = eau precipitable
    DO i = 1, klon
       prw(i) = 0.
       DO k = 1, llm
          prw(i) = prw(i) + q_seri(i, k) * zmasse(i, k)
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
             d_qx(i, k, iq) = (tr_seri(i, k, iq - 2) - qx(i, k, iq)) / dtphys
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

    CALL histwrite_phy("phis", pphis)
    CALL histwrite_phy("aire", airephy)
    CALL histwrite_phy("psol", paprs(:, 1))
    CALL histwrite_phy("precip", rain_fall + snow_fall)
    CALL histwrite_phy("plul", rain_lsc + snow_lsc)
    CALL histwrite_phy("pluc", rain_con + snow_con)
    CALL histwrite_phy("tsol", tsol)
    CALL histwrite_phy("t2m", zt2m)
    CALL histwrite_phy("q2m", zq2m)
    CALL histwrite_phy("u10m", u10m)
    CALL histwrite_phy("v10m", v10m)
    CALL histwrite_phy("snow", snow_fall)
    CALL histwrite_phy("cdrm", cdragm)
    CALL histwrite_phy("cdrh", cdragh)
    CALL histwrite_phy("topl", toplw)
    CALL histwrite_phy("evap", evap)
    CALL histwrite_phy("sols", solsw)
    CALL histwrite_phy("soll", sollw)
    CALL histwrite_phy("solldown", sollwdown)
    CALL histwrite_phy("bils", bils)
    CALL histwrite_phy("sens", - sens)
    CALL histwrite_phy("fder", fder)
    CALL histwrite_phy("dtsvdfo", d_ts(:, is_oce))
    CALL histwrite_phy("dtsvdft", d_ts(:, is_ter))
    CALL histwrite_phy("dtsvdfg", d_ts(:, is_lic))
    CALL histwrite_phy("dtsvdfi", d_ts(:, is_sic))

    DO nsrf = 1, nbsrf
       CALL histwrite_phy("pourc_"//clnsurf(nsrf), pctsrf(:, nsrf) * 100.)
       CALL histwrite_phy("fract_"//clnsurf(nsrf), pctsrf(:, nsrf))
       CALL histwrite_phy("sens_"//clnsurf(nsrf), flux_t(:, nsrf))
       CALL histwrite_phy("lat_"//clnsurf(nsrf), fluxlat(:, nsrf))
       CALL histwrite_phy("tsol_"//clnsurf(nsrf), ftsol(:, nsrf))
       CALL histwrite_phy("taux_"//clnsurf(nsrf), flux_u(:, nsrf))
       CALL histwrite_phy("tauy_"//clnsurf(nsrf), flux_v(:, nsrf))
       CALL histwrite_phy("rugs_"//clnsurf(nsrf), frugs(:, nsrf))
       CALL histwrite_phy("albe_"//clnsurf(nsrf), falbe(:, nsrf))
       CALL histwrite_phy("u10m_"//clnsurf(nsrf), u10m_srf(:, nsrf))
       CALL histwrite_phy("v10m_"//clnsurf(nsrf), v10m_srf(:, nsrf))
    END DO

    CALL histwrite_phy("albs", albsol)
    CALL histwrite_phy("tro3", wo * dobson_u * 1e3 / zmasse / rmo3 * md)
    CALL histwrite_phy("rugs", zxrugs)
    CALL histwrite_phy("s_pblh", s_pblh)
    CALL histwrite_phy("s_pblt", s_pblt)
    CALL histwrite_phy("s_lcl", s_lcl)
    CALL histwrite_phy("s_capCL", s_capCL)
    CALL histwrite_phy("s_oliqCL", s_oliqCL)
    CALL histwrite_phy("s_cteiCL", s_cteiCL)
    CALL histwrite_phy("s_therm", s_therm)
    CALL histwrite_phy("s_trmb1", s_trmb1)
    CALL histwrite_phy("s_trmb2", s_trmb2)
    CALL histwrite_phy("s_trmb3", s_trmb3)

    if (conv_emanuel) then
       CALL histwrite_phy("ptop", ema_pct)
       CALL histwrite_phy("dnwd0", - mp)
    end if

    CALL histwrite_phy("temp", t_seri)
    CALL histwrite_phy("vitu", u_seri)
    CALL histwrite_phy("vitv", v_seri)
    CALL histwrite_phy("geop", zphi)
    CALL histwrite_phy("pres", play)
    CALL histwrite_phy("dtvdf", d_t_vdf)
    CALL histwrite_phy("dqvdf", d_q_vdf)
    CALL histwrite_phy("rhum", zx_rh)
    CALL histwrite_phy("d_t_ec", d_t_ec)
    CALL histwrite_phy("dtsw0", heat0 / 86400.)
    CALL histwrite_phy("dtlw0", - cool0 / 86400.)
    CALL histwrite_phy("msnow", sum(fsnow * pctsrf, dim = 2))
    call histwrite_phy("qsurf", sum(fqsurf * pctsrf, dim = 2))

    if (ok_instan) call histsync(nid_ins)

    IF (lafin) then
       call NF95_CLOSE(ncid_startphy)
       CALL phyredem(pctsrf, ftsol, ftsoil, fqsurf, qsol, &
            fsnow, falbe, fevap, rain_fall, snow_fall, solsw, sollw, dlw, &
            radsol, frugs, agesno, zmea, zstd, zsig, zgam, zthe, zpic, zval, &
            t_ancien, q_ancien, rnebcon, ratqs, clwcon, run_off_lic_0, sig1, &
            w01)
    end IF

    firstcal = .FALSE.

  END SUBROUTINE physiq

end module physiq_m
