module physiq_m

  IMPLICIT none

contains

  SUBROUTINE physiq(lafin, dayvrai, time, paprs, play, pphi, pphis, u, v, t, &
       qx, omega, d_u, d_v, d_t, d_qx)

    ! From phylmd/physiq.F, version 1.22 2006/02/20 09:38:28
    ! (subversion revision 678)

    ! Author: Z. X. Li (LMD/CNRS) 1993

    ! This is the main procedure for the "physics" part of the program.

    ! Libraries:
    use netcdf95, only: NF95_CLOSE
    use nr_util, only: assert

    use aaam_bud_m, only: aaam_bud
    USE abort_gcm_m, ONLY: abort_gcm
    use comgeom, only:  aire_2d
    use ajsec_m, only: ajsec
    use calltherm_m, only: calltherm
    USE clesphys, ONLY: cdhmax, cdmmax, ok_instan
    USE clesphys2, ONLY: conv_emanuel, nbapp_rad, ok_orodr, ok_orolf
    USE conf_interface_m, ONLY: conf_interface
    USE pbl_surface_m, ONLY: pbl_surface
    use clouds_gno_m, only: clouds_gno
    USE concvl_m, ONLY: concvl
    USE conf_gcm_m, ONLY: nday, lmt_pas, dtphys
    USE conf_phys_m, ONLY: conf_phys
    use conflx_m, only: conflx
    USE ctherm_m, ONLY: iflag_thermals, ctherm
    use cv30_param_m, only: cv30_param
    use diagcld1_m, only: diagcld1
    USE dimensions, ONLY: llm, nqmx
    USE dimphy, ONLY: klon
    USE dimsoil, ONLY: nsoilmx
    use drag_noro_m, only: drag_noro
    USE fcttre, ONLY: foeew
    use fisrtilp_m, only: fisrtilp
    use grid_change, only: dyn_phy
    USE hgardfou_m, ONLY: hgardfou
    USE histsync_m, ONLY: histsync
    USE histwrite_phy_m, ONLY: histwrite_phy
    USE indicesol, ONLY: clnsurf, epsfra, nbsrf
    USE ini_histins_m, ONLY: ini_histins, nid_ins
    use lift_noro_m, only: lift_noro
    use newmicro_m, only: newmicro
    USE orbite_m, ONLY: orbite
    USE ozonecm_m, ONLY: ozonecm
    USE phyetat0_m, ONLY: phyetat0, itau_phy
    USE phyredem_m, ONLY: phyredem
    USE phyredem0_m, ONLY: phyredem0
    USE phytrac_m, ONLY: phytrac
    use radlwsw_m, only: radlwsw
    use yoegwd, only: sugwd
    USE suphec_m, ONLY: rcpd, retv, rg, rlvtt, romega, rtt, rmo3, md
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
    REAL, intent(in):: u(:, :) ! (klon, llm) zonal wind, in m / s
    REAL, intent(in):: v(:, :) ! (klon, llm) meridional wind, in m / s
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

    ! pour phystoke avec thermiques
    REAL fm_therm(klon, llm + 1)
    REAL entr_therm(klon, llm)
    real, save, allocatable:: q2(:, :, :) ! (klon, llm + 1, nbsrf)

    INTEGER, PARAMETER:: ivap = 1 ! indice de traceur pour vapeur d'eau
    INTEGER, PARAMETER:: iliq = 2 ! indice de traceur pour eau liquide

    REAL, save, allocatable:: t_ancien(:, :), q_ancien(:, :) ! (klon, llm)
    LOGICAL, save:: ancien_ok

    REAL d_t_dyn(klon, llm) ! tendance dynamique pour "t" (K / s)
    REAL d_q_dyn(klon, llm) ! tendance dynamique pour "q" (kg / kg / s)

    real da(klon, llm), phi(klon, llm, llm), mp(klon, llm)

    REAL, save, allocatable:: swdn0(:, :), swdn(:, :) ! (klon, llm + 1)
    REAL, save, allocatable:: swup0(:, :), swup(:, :) ! (klon, llm + 1)

    REAL, save, allocatable:: lwdn0(:, :), lwdn(:, :) ! (klon, llm + 1)
    REAL, save, allocatable:: lwup0(:, :), lwup(:, :) ! (klon, llm + 1)

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

    REAL, save, allocatable:: radsol(:) ! (klon)
    ! Bilan radiatif net au sol (W/m2), positif vers le bas. Must be
    ! saved because radlwsw is not called at every time step.

    REAL, save, allocatable:: ftsol(:, :) ! (klon, nbsrf)
    ! skin temperature of surface fraction, in K

    REAL, save, allocatable:: ftsoil(:, :, :) ! (klon, nsoilmx, nbsrf)
    ! temperature of surface fraction inside the ground, in K, layer 1
    ! nearest to the surface

    REAL fluxlat(klon, nbsrf) ! flux de chaleur latente, en W m-2

    REAL, save, allocatable:: fqsurf(:, :) ! (klon, nbsrf)
    ! humidite de l'air au contact de la surface

    REAL, save, allocatable:: qsol(:) ! (klon)
    ! column-density of water in soil, in kg m-2

    REAL, save, allocatable:: fsnow(:, :) ! (klon, nbsrf)
    ! column-density of mass of snow at the surface, in kg m-2

    REAL, save, allocatable:: falbe(:, :) ! (klon, nbsrf)
    ! albedo visible par type de surface

    ! Param\`etres de l'orographie \`a l'\'echelle sous-maille (OESM) :

    REAL, save, allocatable:: zmea(:) ! (klon) ! orographie moyenne
    REAL, save, allocatable:: zstd(:) ! (klon) ! deviation standard de l'OESM
    REAL, save, allocatable:: zsig(:) ! (klon) ! pente de l'OESM
    REAL, save, allocatable:: zgam(:) ! (klon) ! anisotropie de l'OESM
    REAL, save, allocatable:: zthe(:) ! (klon) ! orientation de l'OESM
    REAL, save, allocatable:: zpic(:) ! (klon) ! Maximum de l'OESM
    REAL, save, allocatable:: zval(:) ! (klon) ! Minimum de l'OESM

    REAL, save, allocatable:: rugoro(:) ! (klon)
    ! longueur de rugosite de l'OESM

    REAL zulow(klon), zvlow(klon)
    INTEGER ktest(klon)

    REAL, save, allocatable:: agesno(:, :) ! (klon, nbsrf) ! age de la neige
    REAL, save, allocatable:: run_off_lic_0(:) ! (klon)

    ! Variables li\'ees \`a la convection d'Emanuel :
    REAL, save, allocatable:: Ma(:, :) ! (klon, llm) ! undilute upward mass flux
    REAL, save, allocatable:: sig1(:, :), w01(:, :) ! (klon, llm)

    ! Variables pour la couche limite :
    REAL cdragh(klon) ! drag coefficient for T and Q
    REAL cdragm(klon) ! drag coefficient for wind

    REAL coefh(klon, 2:llm) ! coef d'echange pour phytrac

    REAL, save, allocatable:: ffonte(:, :) ! (klon, nbsrf)
    ! flux thermique utilise pour fondre la neige

    REAL fqcalving(klon, nbsrf)
    ! flux d'eau "perdue" par la surface et n\'ecessaire pour limiter
    ! la hauteur de neige, en kg / m2 / s

    REAL zxffonte(klon)

    REAL, save, allocatable:: pfrac_impa(:, :) ! (klon, llm)
    ! Produits des coefs lessivage impaction
    
    REAL, save, allocatable:: pfrac_nucl(:, :) ! (klon, llm)
    ! Produits des coefs lessivage nucleation

    REAL, save, allocatable:: pfrac_1nucl(:, :) ! (klon, llm)
    ! Produits des coefs lessi nucl (alpha = 1)

    REAL frac_impa(klon, llm) ! fraction d'a\'erosols lessiv\'es (impaction)
    REAL frac_nucl(klon, llm) ! idem (nucleation)

    REAL, save, allocatable:: rain_fall(:) ! (klon)
    ! liquid water mass flux (kg / m2 / s), positive down

    REAL, save, allocatable:: snow_fall(:) ! (klon)
    ! solid water mass flux (kg / m2 / s), positive down

    REAL rain_tiedtke(klon), snow_tiedtke(klon)

    REAL evap(klon) ! flux d'\'evaporation au sol
    REAL sens(klon) ! flux de chaleur sensible au sol
    REAL, save, allocatable:: dlw(:) ! (klon) derivative of infra-red flux
    REAL ve(klon) ! integr. verticale du transport meri. de l'energie
    REAL vq(klon) ! integr. verticale du transport meri. de l'eau
    REAL ue(klon) ! integr. verticale du transport zonal de l'energie
    REAL uq(klon) ! integr. verticale du transport zonal de l'eau

    REAL, save, allocatable:: frugs(:, :) ! (klon, nbsrf) ! longueur de rugosite
    REAL zxrugs(klon) ! longueur de rugosite

    ! Conditions aux limites

    INTEGER julien
    REAL, save, allocatable:: pctsrf(:, :) ! (klon, nbsrf) percentage of surface
    
    REAL, save, allocatable:: albsol(:) ! (klon)
    ! albedo du sol total, visible, moyen par maille
    
    REAL, SAVE, ALLOCATABLE:: wo(:, :) ! (klon, llm)
    ! column density of ozone in a cell, in kDU
    
    real, parameter:: dobson_u = 2.1415e-05 ! Dobson unit, in kg m-2

    real, save, allocatable:: clwcon(:, :), rnebcon(:, :) ! (klon, llm)
    real, save, allocatable:: clwcon0(:, :), rnebcon0(:, :) ! (klon, llm)

    REAL rhcl(klon, llm) ! humidit\'e relative ciel clair
    REAL dialiq(klon, llm) ! eau liquide nuageuse
    REAL diafra(klon, llm) ! fraction nuageuse
    REAL cldliq(klon, llm) ! eau liquide nuageuse
    REAL cldfra(klon, llm) ! fraction nuageuse
    REAL cldtau(klon, llm) ! \'epaisseur optique
    REAL cldemi(klon, llm) ! \'emissivit\'e infrarouge

    REAL flux_q(klon, nbsrf) ! flux turbulent d'humidite à la surface

    REAL flux_t(klon, nbsrf)
    ! flux de chaleur sensible (c_p T) (W / m2) (orientation positive
    ! vers le bas) à la surface

    REAL flux_u(klon, nbsrf), flux_v(klon, nbsrf)
    ! tension du vent (flux turbulent de vent) à la surface, en Pa

    ! Le rayonnement n'est pas calcul\'e tous les pas, il faut donc que
    ! les variables soient r\'emanentes.
    REAL, save, allocatable:: heat(:, :) ! (klon, llm) ! chauffage solaire

    REAL, save, allocatable:: heat0(:, :) ! (klon, llm)
    ! chauffage solaire ciel clair

    REAL, save, allocatable:: cool(:, :) ! (klon, llm)
    ! refroidissement infrarouge

    REAL, save, allocatable:: cool0(:, :) ! (klon, llm)
    ! refroidissement infrarouge ciel clair

    REAL, save, allocatable:: topsw(:), toplw(:), solsw(:) ! (klon)

    REAL, save, allocatable:: sollw(:) ! (klon)
    ! surface net downward longwave flux, in W m-2

    real, save, allocatable:: sollwdown(:) ! (klon)
    ! downwelling longwave flux at surface

    REAL, save, allocatable:: topsw0(:), toplw0(:), solsw0(:), sollw0(:)
    ! (klon)

    REAL conv_q(klon, llm) ! convergence de l'humidite (kg / kg / s)
    REAL conv_t(klon, llm) ! convergence of temperature (K / s)

    REAL cldl(klon), cldm(klon), cldh(klon) ! nuages bas, moyen et haut
    REAL cldt(klon), cldq(klon) ! nuage total, eau liquide integree

    REAL zxfluxlat(klon)
    REAL dist ! distance Terre-Soleil, en ua
    real mu0(klon), fract(klon)
    real longi
    REAL z_avant(klon), z_apres(klon), z_factor(klon)
    REAL zb
    REAL zx_qs, zcor
    real zqsat(klon, llm)
    INTEGER i, k, iq, nsrf
    REAL zphi(klon, llm)

    ! Variables pour la convection de K. Emanuel :

    REAL upwd(klon, llm) ! saturated updraft mass flux
    REAL dnwd(klon, llm) ! saturated downdraft mass flux
    REAL, save, allocatable:: cape(:) ! (klon)

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

    INTEGER, save, allocatable:: ibas_con(:), itop_con(:) ! (klon)
    real ema_pct(klon) ! Emanuel pressure at cloud top, in Pa

    REAL rain_con(klon)
    real rain_lsc(klon)
    REAL snow_con(klon) ! neige (mm / s)
    real snow_lsc(klon)

    REAL d_u_vdf(klon, llm), d_v_vdf(klon, llm)
    REAL d_t_vdf(klon, llm), d_q_vdf(klon, llm)

    REAL d_u_oro(klon, llm), d_v_oro(klon, llm)
    REAL d_t_oro(klon, llm)
    REAL d_u_lif(klon, llm), d_v_lif(klon, llm)
    REAL d_t_lif(klon, llm)

    REAL, save, allocatable:: ratqs(:, :) ! (klon, llm)
    real ratqss(klon, llm), ratqsc(klon, llm)
    real:: ratqsbas = 0.01, ratqshaut = 0.3

    ! Param\`etres li\'es au nouveau sch\'ema de nuages :
    real:: fact_cldcon = 0.375

    real:: facttemps = 1.e-4 ! in s-1
    ! 1 / facttemps est le temps de relaxation des ratqs.

    real facteur

    integer:: iflag_cldcon = 1 ! allowed values: - 2, ..., 3
    logical ptconv(klon, llm)

    ! Variables pour effectuer les appels en s\'erie :
    REAL t_seri(klon, llm)
    real q_seri(klon, llm) ! mass fraction of water vapor
    REAL ql_seri(klon, llm)
    REAL u_seri(klon, llm), v_seri(klon, llm) ! wind, in m s-1
    REAL tr_seri(klon, llm, nqmx - 2)

    REAL zx_rh(klon, llm)

    REAL zustrdr(klon), zvstrdr(klon)
    REAL zustrli(klon), zvstrli(klon)
    REAL aam, torsfc

    REAL ve_lay(klon, llm) ! transport meri. de l'energie a chaque niveau vert.
    REAL vq_lay(klon, llm) ! transport meri. de l'eau a chaque niveau vert.
    REAL ue_lay(klon, llm) ! transport zonal de l'energie a chaque niveau vert.
    REAL uq_lay(klon, llm) ! transport zonal de l'eau a chaque niveau vert.

    REAL tsol(klon)

    REAL d_t_ec(klon, llm)
    ! tendance due \`a la conversion d'\'energie cin\'etique en
    ! énergie thermique

    REAL, save, allocatable:: t2m(:, :), q2m(:, :) ! (klon, nbsrf) 
    ! temperature and humidity at 2 m

    REAL, save, allocatable:: u10m_srf(:, :), v10m_srf(:, :) ! (klon, nbsrf)
    ! composantes du vent \`a 10 m

    REAL zt2m(klon), zq2m(klon) ! température, humidité 2 m moyenne sur 1 maille
    REAL u10m(klon), v10m(klon) ! vent \`a 10 m moyenn\' sur les sous-surfaces

    ! Aerosol effects:

    REAL, save, allocatable:: topswad(:), solswad(:) ! (klon)
    ! aerosol direct effect

    LOGICAL:: ok_ade = .false. ! apply aerosol direct effect

    REAL:: bl95_b0 = 2., bl95_b1 = 0.2
    ! Parameters in equation (D) of Boucher and Lohmann (1995, Tellus
    ! B). They link cloud droplet number concentration to aerosol mass
    ! concentration.

    real zmasse(klon, llm)
    ! (column-density of mass of air in a cell, in kg m-2)

    integer, save:: ncid_startphy
    real, save, allocatable:: airephy(:) ! (klon)

    namelist /physiq_nml/ fact_cldcon, facttemps, iflag_cldcon, ratqsbas, &
         ratqshaut, ok_ade, bl95_b0, bl95_b1

    !----------------------------------------------------------------

    IF (nqmx < 2) CALL abort_gcm('physiq', &
         'eaux vapeur et liquide sont indispensables')

    test_firstcal: IF (firstcal) THEN
       allocate(q2(klon, llm + 1, nbsrf))
       allocate(t_ancien(klon, llm), q_ancien(klon, llm))
       allocate(swdn0(klon, llm + 1), swdn(klon, llm + 1))
       allocate(swup0(klon, llm + 1), swup(klon, llm + 1))
       allocate(lwdn0(klon, llm + 1), lwdn(klon, llm + 1))
       allocate(lwup0(klon, llm + 1), lwup(klon, llm + 1))
       allocate(radsol(klon))
       allocate(ftsol(klon, nbsrf))
       allocate(ftsoil(klon, nsoilmx, nbsrf))
       allocate(fqsurf(klon, nbsrf))
       allocate(qsol(klon))
       allocate(fsnow(klon, nbsrf))
       allocate(falbe(klon, nbsrf))
       allocate(zmea(klon))
       allocate(zstd(klon))
       allocate(zsig(klon))
       allocate(zgam(klon))
       allocate(zthe(klon))
       allocate(zpic(klon))
       allocate(zval(klon))
       allocate(rugoro(klon))
       allocate(agesno(klon, nbsrf))
       allocate(run_off_lic_0(klon))
       allocate(Ma(klon, llm))
       allocate(sig1(klon, llm), w01(klon, llm))
       allocate(ffonte(klon, nbsrf))
       allocate(pfrac_impa(klon, llm))
       allocate(pfrac_nucl(klon, llm))
       allocate(pfrac_1nucl(klon, llm))
       allocate(rain_fall(klon))
       allocate(snow_fall(klon))
       allocate(dlw(klon))
       allocate(frugs(klon, nbsrf))
       allocate(pctsrf(klon, nbsrf))
       allocate(albsol(klon))
       ALLOCATE(wo(klon, llm))
       allocate(clwcon(klon, llm), rnebcon(klon, llm))
       allocate(clwcon0(klon, llm), rnebcon0(klon, llm))
       allocate(cape(klon))
       allocate(ibas_con(klon), itop_con(klon))
       allocate(ratqs(klon, llm))
       allocate(t2m(klon, nbsrf), q2m(klon, nbsrf))
       allocate(u10m_srf(klon, nbsrf), v10m_srf(klon, nbsrf))
       allocate(topswad(klon), solswad(klon))
       allocate(airephy(klon))
       allocate(heat(klon, llm))
       allocate(heat0(klon, llm))
       allocate(cool(klon, llm))
       allocate(cool0(klon, llm))
       allocate(topsw(klon), toplw(klon), solsw(klon))
       allocate(sollw(klon))
       allocate(sollwdown(klon))
       allocate(topsw0(klon), toplw0(klon), solsw0(klon), sollw0(klon))

       ! initialiser
       u10m_srf = 0.
       v10m_srf = 0.
       t2m = 0.
       q2m = 0.
       ffonte = 0.
       rnebcon0 = 0.
       clwcon0 = 0.
       clwcon = 0.

       print *, "Enter namelist 'physiq_nml'."
       read(unit=*, nml=physiq_nml)
       write(unit_nml, nml=physiq_nml)
       call assert(iflag_cldcon >= - 2 .and. iflag_cldcon <= 3, &
            "physiq iflag_cldcon")

       call ctherm
       call conf_phys
       frugs = 0.
       CALL phyetat0(pctsrf, ftsol, ftsoil, fqsurf, qsol, fsnow, falbe, &
            rain_fall, snow_fall, solsw, sollw, dlw, radsol, frugs, agesno, &
            zmea, zstd, zsig, zgam, zthe, zpic, zval, t_ancien, q_ancien, &
            ancien_ok, rnebcon, ratqs, clwcon, run_off_lic_0, sig1, w01, &
            ncid_startphy)

       ! ATTENTION : il faudra a terme relire q2 dans l'etat initial
       q2 = 1e-8

       radpas = lmt_pas / nbapp_rad
       print *, "radpas = ", radpas

       ! Initialisation pour le sch\'ema de convection d'Emanuel :
       IF (conv_emanuel) THEN
          ibas_con = 1
          itop_con = 1
          CALL cv30_param
       ENDIF

       IF (ok_orodr) THEN
          rugoro = MAX(1e-5, zstd * zsig / 2)
          CALL SUGWD(paprs, play)
       else
          rugoro = 0.
       ENDIF

       ! Initialisation des sorties
       call ini_histins
       CALL phyredem0(itau_phy + nday * lmt_pas)
       call conf_interface
       airephy = pack(aire_2d, dyn_phy)
    ENDIF test_firstcal

    ! We will modify variables *_seri and we will not touch variables
    ! u, v, t, qx:
    t_seri = t
    u_seri = u
    v_seri = v
    q_seri = qx(:, :, ivap)
    ql_seri = qx(:, :, iliq)
    tr_seri = qx(:, :, 3:nqmx)

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

    CALL pbl_surface(pctsrf, t_seri, q_seri, u_seri, v_seri, julien, mu0, &
         ftsol, cdmmax, cdhmax, ftsoil, qsol, paprs, play, fsnow, fqsurf, &
         falbe, fluxlat, rain_fall, snow_fall, frugs, agesno, rugoro, d_t_vdf, &
         d_q_vdf, d_u_vdf, d_v_vdf, flux_t, flux_q, flux_u, flux_v, cdragh, &
         cdragm, q2, coefh, t2m, q2m, u10m_srf, v10m_srf, fqcalving, ffonte, &
         run_off_lic_0, albsol, sollw, solsw, tsol, dlw)

    ! Incr\'ementation des flux :
    sens = sum(flux_t * pctsrf, dim = 2)
    evap = - sum(flux_q * pctsrf, dim = 2)

    DO k = 1, llm
       DO i = 1, klon
          t_seri(i, k) = t_seri(i, k) + d_t_vdf(i, k)
          q_seri(i, k) = q_seri(i, k) + d_q_vdf(i, k)
          u_seri(i, k) = u_seri(i, k) + d_u_vdf(i, k)
          v_seri(i, k) = v_seri(i, k) + d_v_vdf(i, k)
       ENDDO
    ENDDO

    call assert(abs(sum(pctsrf, dim = 2) - 1.) <= EPSFRA, 'physiq: pctsrf')
    zxfluxlat = sum(fluxlat * pctsrf, dim = 2)
    zt2m = sum(t2m * pctsrf, dim = 2)
    zq2m = sum(q2m * pctsrf, dim = 2)
    u10m = sum(u10m_srf * pctsrf, dim = 2)
    v10m = sum(v10m_srf * pctsrf, dim = 2)
    zxffonte = sum(ffonte * pctsrf, dim = 2)

    ! Si une sous-fraction n'existe pas, elle prend la valeur moyenne :
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          IF (pctsrf(i, nsrf) < epsfra) then
             t2m(i, nsrf) = zt2m(i)
             q2m(i, nsrf) = zq2m(i)
             u10m_srf(i, nsrf) = u10m(i)
             v10m_srf(i, nsrf) = v10m(i)
             ffonte(i, nsrf) = zxffonte(i)
          end IF
       ENDDO
    ENDDO

    ! Appeler la convection

    if (conv_emanuel) then
       CALL concvl(paprs, play, t_seri, q_seri, u_seri, v_seri, sig1, w01, &
            d_t_con, d_q_con, d_u_con, d_v_con, rain_con, ibas_con, itop_con, &
            upwd, dnwd, Ma, cape, iflagctrl, clwcon0, pmflxr, da, phi, mp)
       snow_con = 0.
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
       u_seri = u_seri + d_u_con
       v_seri = v_seri + d_v_con
    else
       conv_q = d_q_dyn + d_q_vdf / dtphys
       conv_t = d_t_dyn + d_t_vdf / dtphys
       z_avant = sum((q_seri + ql_seri) * zmasse, dim=2)
       CALL conflx(paprs, play, t_seri(:, llm:1:- 1), q_seri(:, llm:1:- 1), &
            conv_t, conv_q, - evap, omega, d_t_con, d_q_con, rain_con, &
            snow_con, mfu(:, llm:1:- 1), mfd(:, llm:1:- 1), pen_u, pde_u, &
            pen_d, pde_d, kcbot, kctop, kdtop, pmflxr, pmflxs)
       WHERE (rain_con < 0.) rain_con = 0.
       WHERE (snow_con < 0.) snow_con = 0.
       ibas_con = llm + 1 - kcbot
       itop_con = llm + 1 - kctop
    END if

    t_seri = t_seri + d_t_con
    q_seri = q_seri + d_q_con

    IF (.not. conv_emanuel) THEN
       z_apres = sum((q_seri + ql_seri) * zmasse, dim=2)
       z_factor = (z_avant - (rain_con + snow_con) * dtphys) / z_apres
       DO k = 1, llm
          DO i = 1, klon
             IF (z_factor(i) /= 1.) THEN
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

    if (iflag_thermals) then
       call calltherm(play, paprs, pphi, u_seri, v_seri, t_seri, q_seri, &
            d_u_ajs, d_v_ajs, d_t_ajs, d_q_ajs, fm_therm, entr_therm)
    else
       CALL ajsec(paprs, play, t_seri, q_seri, d_t_ajs, d_q_ajs)
       t_seri = t_seri + d_t_ajs
       q_seri = q_seri + d_q_ajs
    endif

    ! Caclul des ratqs

    if (iflag_cldcon == 1) then
       ! ratqs convectifs \`a l'ancienne en fonction de (q(z = 0) - q) / q
       ! on \'ecrase le tableau ratqsc calcul\'e par clouds_gno
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
       ratqs = max(ratqs * exp(- dtphys * facttemps), ratqss)
       ratqs = max(ratqs, ratqsc)
    else
       ! on ne prend que le ratqs stable pour fisrtilp
       ratqs = ratqss
    endif

    CALL fisrtilp(paprs, play, t_seri, q_seri, ptconv, ratqs, d_t_lsc, &
         d_q_lsc, d_ql_lsc, rneb, cldliq, rain_lsc, snow_lsc, pfrac_impa, &
         pfrac_nucl, pfrac_1nucl, frac_impa, frac_nucl, prfl, psfl, rhcl)

    WHERE (rain_lsc < 0) rain_lsc = 0.
    WHERE (snow_lsc < 0) snow_lsc = 0.
    DO k = 1, llm
       DO i = 1, klon
          t_seri(i, k) = t_seri(i, k) + d_t_lsc(i, k)
          q_seri(i, k) = q_seri(i, k) + d_q_lsc(i, k)
          ql_seri(i, k) = ql_seri(i, k) + d_ql_lsc(i, k)
          cldfra(i, k) = rneb(i, k)
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

    ! Precipitation totale
    DO i = 1, klon
       rain_fall(i) = rain_con(i) + rain_lsc(i)
       snow_fall(i) = snow_con(i) + snow_lsc(i)
    ENDDO

    ! Humidit\'e relative pour diagnostic :
    DO k = 1, llm
       DO i = 1, klon
          zx_qs = r2es * FOEEW(t_seri(i, k), rtt >= t_seri(i, k)) / play(i, k)
          zx_qs = MIN(0.5, zx_qs)
          zcor = 1. / (1. - retv * zx_qs)
          zx_qs = zx_qs * zcor
          zx_rh(i, k) = q_seri(i, k) / zx_qs
          zqsat(i, k) = zx_qs
       ENDDO
    ENDDO

    ! Param\`etres optiques des nuages et quelques param\`etres pour
    ! diagnostics :
    CALL newmicro(paprs, play, t_seri, cldliq, cldfra, cldtau, cldemi, cldh, &
         cldl, cldm, cldt, cldq, flwp, fiwp, flwc, fiwc)

    IF (MOD(itap - 1, radpas) == 0) THEN
       wo = ozonecm(REAL(julien), paprs)
       albsol = sum(falbe * pctsrf, dim = 2)
       CALL radlwsw(dist, mu0, fract, paprs, play, tsol, albsol, t_seri, &
            q_seri, wo, cldfra, cldemi, cldtau, heat, heat0, cool, cool0, &
            radsol, topsw, toplw, solsw, sollw, sollwdown, topsw0, toplw0, &
            solsw0, sollw0, lwdn0, lwdn, lwup0, lwup, swdn0, swdn, swup0, &
            swup, ok_ade, topswad, solswad)
    ENDIF

    ! Ajouter la tendance des rayonnements (tous les pas)
    DO k = 1, llm
       DO i = 1, klon
          t_seri(i, k) = t_seri(i, k) + (heat(i, k) - cool(i, k)) * dtphys &
               / 86400.
       ENDDO
    ENDDO

    ! Param\'etrisation de l'orographie \`a l'\'echelle sous-maille :

    IF (ok_orodr) THEN
       ! S\'election des points pour lesquels le sch\'ema est actif :
       DO i = 1, klon
          IF (zpic(i) - zmea(i) > 100. .AND. zstd(i) > 10.) THEN
             ktest(i) = 1
          else
             ktest(i) = 0
          ENDIF
       ENDDO

       CALL drag_noro(paprs, play, zmea, zstd, zsig, zgam, zthe, zpic, zval, &
            ktest, t_seri, u_seri, v_seri, zulow, zvlow, zustrdr, zvstrdr, &
            d_t_oro, d_u_oro, d_v_oro)

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
       DO i = 1, klon
          IF (zpic(i) - zmea(i) > 100.) THEN
             ktest(i) = 1
          else
             ktest(i) = 0
          ENDIF
       ENDDO

       CALL lift_noro(paprs, play, zmea, zstd, zpic, ktest, t_seri, u_seri, &
            v_seri, zulow, zvlow, zustrli, zvstrli, d_t_lif, d_u_lif, d_v_lif)

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
    call phytrac(julien, time, firstcal, lafin, t, paprs, play, mfu, mfd, &
         pde_u, pen_d, coefh, cdragh, fm_therm, entr_therm, u(:, 1), v(:, 1), &
         ftsol, pctsrf, frac_impa, frac_nucl, da, phi, mp, upwd, dnwd, &
         tr_seri, zmasse, ncid_startphy)

    ! Calculer le transport de l'eau et de l'energie (diagnostique)
    CALL transp(paprs, t_seri, q_seri, u_seri, v_seri, zphi, ve, vq, ue, uq)

    ! diag. bilKP

    CALL transp_lay(paprs, t_seri, q_seri, u_seri, v_seri, zphi, ve_lay, &
         vq_lay, ue_lay, uq_lay)

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
    CALL histwrite_phy("rls", sollw)
    CALL histwrite_phy("solldown", sollwdown)
    CALL histwrite_phy("bils", radsol + sens + zxfluxlat)
    CALL histwrite_phy("sens", sens)
    CALL histwrite_phy("zxfqcalving", sum(fqcalving * pctsrf, dim = 2))
    CALL histwrite_phy("albs", albsol)
    CALL histwrite_phy("tro3", wo * dobson_u * 1e3 / zmasse / rmo3 * md)
    CALL histwrite_phy("rugs", zxrugs)
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
    call histwrite_phy("pmflxr", pmflxr(:, :llm))
    CALL histwrite_phy("msnow", sum(fsnow * pctsrf, dim = 2))
    call histwrite_phy("qsurf", sum(fqsurf * pctsrf, dim = 2))
    call histwrite_phy("flat", zxfluxlat)
    call histwrite_phy("rld", lwdn)
    call histwrite_phy("rldcs", lwdn0)
    call histwrite_phy("ffonte", zxffonte)

    DO nsrf = 1, nbsrf
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

    if (conv_emanuel) then
       CALL histwrite_phy("ptop", ema_pct)
       CALL histwrite_phy("dnwd0", - mp)
    end if

    if (ok_instan) call histsync(nid_ins)

    IF (lafin) then
       call NF95_CLOSE(ncid_startphy)
       CALL phyredem(pctsrf, ftsol, ftsoil, fqsurf, qsol, fsnow, falbe, &
            rain_fall, snow_fall, solsw, sollw, dlw, radsol, frugs, agesno, &
            zmea, zstd, zsig, zgam, zthe, zpic, zval, t_ancien, q_ancien, &
            rnebcon, ratqs, clwcon, run_off_lic_0, sig1, w01)
    end IF

    firstcal = .FALSE.

  END SUBROUTINE physiq

end module physiq_m
