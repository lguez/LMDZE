module pbl_surface_m

  IMPLICIT NONE

contains

  SUBROUTINE pbl_surface(pctsrf, t_seri, q_seri, u_seri, v_seri, julien, mu0, &
       ftsol, cdmmax, cdhmax, ftsoil, qsol, paprs, play, fsnow, fqsurf, falbe, &
       fluxlat, rain_fall, snow_fall, frugs, agesno, rugoro, d_t, d_q, d_u, &
       d_v, flux_t, flux_q, flux_u, flux_v, cdragh, cdragm, coefh, t2m, q2m, &
       u10m_srf, v10m_srf, fqcalving, ffonte, run_off_lic_0, albsol, sollw, &
       solsw, tsol, dlw)

    ! From phylmd/clmain.F, version 1.6, 2005/11/16 14:47:19
    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: Aug. 18th, 1993
    ! Objet : interface de couche limite (diffusion verticale)

    ! Tout ce qui a trait aux traceurs est dans "phytrac". Le calcul
    ! de la couche limite pour les traceurs se fait avec "cltrac" et
    ! ne tient pas compte de la diff\'erenciation des sous-fractions
    ! de sol.

    ! Libraries:
    use jumble, only: assert

    use cdrag_m, only: cdrag
    use clqh_m, only: clqh
    use clvent_m, only: clvent
    use coef_diff_turb_m, only: coef_diff_turb
    USE conf_gcm_m, ONLY: lmt_pas
    USE conf_phys_m, ONLY: iflag_pbl
    USE dimphy, ONLY: klev, klon
    USE dimsoil, ONLY: nsoilmx
    use hbtm_m, only: hbtm
    USE histwrite_phy_m, ONLY: histwrite_phy
    USE indicesol, ONLY: epsfra, is_lic, is_oce, is_sic, is_ter, nbsrf
    USE interfoce_lim_m, ONLY: interfoce_lim
    use phyetat0_m, only: masque
    use stdlevvar_m, only: stdlevvar
    USE suphec_m, ONLY: rd, rg, rsigma
    use time_phylmdz, only: itap

    REAL, INTENT(inout):: pctsrf(:, :) ! (klon, nbsrf)
    ! pourcentages de surface de chaque maille

    REAL, INTENT(IN):: t_seri(:, :) ! (klon, klev) air temperature, in K
    REAL, INTENT(IN):: q_seri(:, :) ! (klon, klev) mass fraction of water vapor
    REAL, INTENT(IN):: u_seri(:, :), v_seri(:, :) ! (klon, klev) wind, in m s -1
    INTEGER, INTENT(IN):: julien ! jour de l'annee en cours
    REAL, intent(in):: mu0(:) ! (klon) cosinus de l'angle solaire zenithal

    REAL, INTENT(INout):: ftsol(:, :) ! (klon, nbsrf)
    ! skin temperature of surface fraction, in K

    REAL, INTENT(IN):: cdmmax, cdhmax ! seuils cdrm, cdrh

    REAL, INTENT(inout):: ftsoil(:, :, :) ! (klon, nsoilmx, nbsrf)
    ! temperature of surface fraction inside the ground, in K, layer 1
    ! nearest to the surface

    REAL, INTENT(inout):: qsol(:) ! (klon)
    ! column-density of water in soil, in kg m-2

    REAL, INTENT(IN):: paprs(:, :) ! (klon, klev + 1)
    ! pression a intercouche (Pa)

    REAL, INTENT(IN):: play(:, :) ! (klon, klev)
    ! pression au milieu de couche (Pa)

    REAL, INTENT(inout):: fsnow(:, :) ! (klon, nbsrf)
    ! column-density of mass of snow at the surface, in kg m-2

    REAL, INTENT(inout):: fqsurf(:, :) ! (klon, nbsrf)
    REAL, intent(inout):: falbe(:, :) ! (klon, nbsrf)

    REAL, intent(out):: fluxlat(:, :) ! (klon, nbsrf)
    ! flux de chaleur latente, en W m-2

    REAL, intent(in):: rain_fall(:) ! (klon)
    ! liquid water mass flux (kg / m2 / s), positive down

    REAL, intent(in):: snow_fall(:) ! (klon)
    ! solid water mass flux (kg / m2 / s), positive down

    REAL, intent(inout):: frugs(:, :) ! (klon, nbsrf)
    ! longueur de rugosit\'e (en m)

    real, intent(inout):: agesno(:, :) ! (klon, nbsrf)
    REAL, INTENT(IN):: rugoro(:) ! (klon) longueur de rugosit\'e orographique

    REAL, intent(out):: d_t(:, :), d_q(:, :) ! (klon, klev)
    ! changement pour t_seri et q_seri

    REAL, intent(out):: d_u(:, :), d_v(:, :) ! (klon, klev)
    ! changement pour "u_seri" et "v_seri"

    REAL, intent(out):: flux_t(:, :) ! (klon, nbsrf)
    ! flux de chaleur sensible (c_p T) (W / m2) (orientation positive
    ! vers le bas) à la surface

    REAL, intent(out):: flux_q(:, :) ! (klon, nbsrf)
    ! flux de vapeur d'eau (kg / m2 / s) à la surface

    REAL, intent(out):: flux_u(:, :), flux_v(:, :) ! (klon, nbsrf)
    ! tension du vent (flux turbulent de vent) à la surface, en Pa

    REAL, INTENT(out):: cdragh(:) ! (klon)
    ! drag coefficient for latent and sensible heat fluxes

    REAL, INTENT(out):: cdragm(:) ! (klon)

    ! Ocean slab:

    REAL, intent(out):: coefh(:, 2:) ! (klon, 2:klev) Diffusion
    ! coefficient at layer interface, for heat and humidity, in m2
    ! s-1. Pour pouvoir extraire les coefficients d'\'echange, le
    ! champ "coefh" a \'et\'e cr\'e\'e. Nous avons moyenn\'e les
    ! valeurs de ce champ sur les quatre sous-surfaces du mod\`ele.

    REAL, INTENT(inout):: t2m(:, :), q2m(:, :) ! (klon, nbsrf)

    REAL, INTENT(inout):: u10m_srf(:, :), v10m_srf(:, :) ! (klon, nbsrf)
    ! composantes du vent \`a 10m sans spirale d'Ekman

    REAL, intent(out):: fqcalving(:, :) ! (klon, nbsrf)
    ! flux d'eau "perdue" par la surface et n\'ecessaire pour limiter
    ! la hauteur de neige, en kg / m2 / s

    real, INTENT(inout):: ffonte(:, :) ! (klon, nbsrf)
    ! flux thermique utilis\'e pour fondre la neige

    REAL, intent(inout):: run_off_lic_0(:) ! (klon)

    REAL, intent(in):: albsol(:) ! (klon)
    ! alb\'edo du sol total, dans le domaine spectral visible, moyen
    ! par maille

    REAL, intent(in):: sollw(:) ! (klon)
    ! surface net downward longwave flux, in W m-2

    REAL, intent(in):: solsw(:) ! (klon)
    ! surface net downward shortwave flux, in W m-2

    REAL, intent(out):: tsol(:) ! (klon)
    REAL, intent(inout):: dlw(:) ! (klon) derivative of infra-red flux

    ! Local:

    REAL dflux_t(klon) ! d\'eriv\'ee du flux de chaleur sensible au sol
    REAL dflux_q(klon) ! derive du flux latent at the surface
    REAL d_ts(klon, nbsrf) ! variation of ftsol
    REAL fsollw(klon, nbsrf) ! bilan flux IR pour chaque sous-surface
    REAL fsolsw(klon, nbsrf) ! flux solaire absorb\'e pour chaque sous-surface

    ! la nouvelle repartition des surfaces sortie de l'interface
    REAL, save, allocatable:: pctsrf_new_oce(:), pctsrf_new_sic(:) ! (klon)

    REAL y_fqcalving(klon), y_ffonte(klon)
    real y_run_off_lic_0(klon), y_run_off_lic(klon)
    REAL run_off_lic(klon) ! ruissellement total
    REAL rugmer(klon)
    REAL ytsoil(klon, nsoilmx)
    REAL yts(klon), ypctsrf(klon), yz0_new(klon)
    real yrugos(klon) ! longueur de rugosit\'e, en m
    REAL yalbedo(klon)
    REAL snow(klon) ! column-density of mass of snow at the surface, in kg m-2
    real yqsurf(klon), yagesno(klon)
    real yqsol(klon) ! column-density of water in soil, in kg m-2
    REAL yrain_fall(klon) ! liquid water mass flux (kg / m2 / s), positive down
    REAL ysnow_fall(klon) ! solid water mass flux (kg / m2 / s), positive down
    REAL yrugm(klon), yrugoro(klon)
    REAL yfluxlat(klon)
    REAL tsurf_new(klon)
    REAL y_d_t(klon, klev), y_d_q(klon, klev)
    REAL y_d_u(klon, klev), y_d_v(klon, klev)
    REAL y_flux_t(klon), y_flux_q(klon)
    REAL y_flux_u(klon), y_flux_v(klon)
    REAL y_dflux_t(klon), y_dflux_q(klon)
    REAL ycoefh(klon, 2:klev), ycoefm(klon, 2:klev)
    real ycdragh(klon), ycdragm(klon)
    REAL yu(klon, klev), yv(klon, klev)
    REAL yt(klon, klev), yq(klon, klev)
    REAL ypaprs(klon, klev + 1), ypplay(klon, klev), ydelp(klon, klev)
    REAL yq2(klon, klev + 1)
    REAL delp(klon, klev) ! \'epaisseur de couche
    INTEGER i, k, nisrf
    INTEGER ni(klon), knon, j

    REAL pctsrf_pot(klon, nbsrf)
    ! "pourcentage potentiel" pour tenir compte des \'eventuelles
    ! apparitions ou disparitions de la glace de mer

    REAL yt2m(klon), yq2m(klon), wind10m(klon)
    REAL ustar(klon)
    REAL ypblh(klon)
    REAL ylcl(klon)
    REAL ycapcl(klon)
    REAL yoliqcl(klon)
    REAL ycteicl(klon)
    REAL ypblt(klon)
    REAL ytherm(klon)
    REAL u1(klon), v1(klon)
    REAL t1(klon)
    REAL zgeop(klon, klev) ! geopotential at mid-layer, in m2 s-2

    ! Ionela Musat. Cf. Anne Mathieu : planetary boundary layer, hbtm.
    REAL pblh(klon, nbsrf) ! height of planetary boundary layer
    REAL capcl(klon, nbsrf) ! CAPE de couche limite
    REAL oliqcl(klon, nbsrf) ! eau_liqu integree de couche limite
    REAL cteicl(klon, nbsrf) ! cloud top instab. crit. couche limite
    REAL pblt(klon, nbsrf) ! temp\'erature au nveau HCL
    REAL therm(klon, nbsrf)
    REAL lcl(klon, nbsrf) ! Niveau de condensation de la CLA

    real, save, allocatable:: q2(:, :, :) ! (klon, llm + 1, nbsrf)

    !------------------------------------------------------------

    tsol = sum(ftsol * pctsrf, dim = 2)

    ! R\'epartition sous maille des flux longwave et shortwave
    ! R\'epartition du longwave par sous-surface lin\'earis\'ee

    forall (nisrf = 1:nbsrf)
       fsollw(:, nisrf) = sollw + 4. * RSIGMA * tsol**3 &
            * (tsol - ftsol(:, nisrf))
       fsolsw(:, nisrf) = solsw * (1. - falbe(:, nisrf)) / (1. - albsol)
    END forall

    forall (k = 1:klev) delp(:, k) = paprs(:, k) - paprs(:, k + 1)

    ! Initialization:
    rugmer = 0.
    cdragh = 0.
    cdragm = 0.
    dflux_t = 0.
    dflux_q = 0.
    ypaprs = 0.
    ypplay = 0.
    ydelp = 0.
    yrugoro = 0.
    d_ts = 0.
    flux_t = 0.
    flux_q = 0.
    flux_u = 0.
    flux_v = 0.
    fluxlat = 0.
    d_t = 0.
    d_q = 0.
    d_u = 0.
    d_v = 0.
    coefh = 0.
    fqcalving = 0.
    run_off_lic = 0.
    pblh = huge(0.)
    lcl = huge(0.)
    capCL = huge(0.)
    oliqCL = huge(0.)
    cteiCL = huge(0.)
    pblt = huge(0.)
    therm = huge(0.)

    ! Initialisation des "pourcentages potentiels". On consid\`ere ici qu'on
    ! peut avoir potentiellement de la glace sur tout le domaine oc\'eanique
    ! (\`a affiner).

    pctsrf_pot(:, is_ter) = pctsrf(:, is_ter)
    pctsrf_pot(:, is_lic) = pctsrf(:, is_lic)
    pctsrf_pot(:, is_oce) = 1. - masque
    pctsrf_pot(:, is_sic) = 1. - masque

    if (itap == 1) then
       allocate(pctsrf_new_oce(klon), pctsrf_new_sic(klon))
       allocate(q2(klon, klev + 1, nbsrf))
       ! ATTENTION : il faudra a terme relire q2 dans l'etat initial
       q2 = 1e-8
    end if

    ! Tester si c'est le moment de lire le fichier :
    if (mod(itap - 1, lmt_pas) == 0) &
         CALL interfoce_lim(julien, pctsrf_new_oce, pctsrf_new_sic)

    ! Boucler sur toutes les sous-fractions du sol:

    loop_surface: DO nisrf = 1, nbsrf
       ! Define ni and knon:

       ni = 0
       knon = 0

       DO i = 1, klon
          ! Pour d\'eterminer le domaine \`a traiter, on utilise les surfaces
          ! "potentielles"
          IF (pctsrf_pot(i, nisrf) > epsfra) THEN
             knon = knon + 1
             ni(knon) = i
          END IF
       END DO

       if_knon: IF (knon /= 0) then
          ypctsrf(:knon) = pctsrf(ni(:knon), nisrf)
          yts(:knon) = ftsol(ni(:knon), nisrf)
          snow(:knon) = fsnow(ni(:knon), nisrf)
          yqsurf(:knon) = fqsurf(ni(:knon), nisrf)
          yrain_fall(:knon) = rain_fall(ni(:knon))
          ysnow_fall(:knon) = snow_fall(ni(:knon))
          yagesno(:knon) = agesno(ni(:knon), nisrf)
          yrugos(:knon) = frugs(ni(:knon), nisrf)
          yrugoro(:knon) = rugoro(ni(:knon))
          ypaprs(:knon, :) = paprs(ni(:knon), :)
          y_run_off_lic_0(:knon) = run_off_lic_0(ni(:knon))

          ! For continent, copy soil water content
          IF (nisrf == is_ter) yqsol(:knon) = qsol(ni(:knon))

          ytsoil(:knon, :) = ftsoil(ni(:knon), :, nisrf)

          DO k = 1, klev
             DO j = 1, knon
                i = ni(j)
                ypplay(j, k) = play(i, k)
                ydelp(j, k) = delp(i, k)
                yu(j, k) = u_seri(i, k)
                yv(j, k) = v_seri(i, k)
                yt(j, k) = t_seri(i, k)
                yq(j, k) = q_seri(i, k)
             END DO
          END DO

          ! Calculer les géopotentiels de chaque couche:

          zgeop(:knon, 1) = RD * yt(:knon, 1) / (0.5 * (ypaprs(:knon, 1) &
               + ypplay(:knon, 1))) * (ypaprs(:knon, 1) - ypplay(:knon, 1))

          DO k = 2, klev
             zgeop(:knon, k) = zgeop(:knon, k - 1) + RD * 0.5 &
                  * (yt(:knon, k - 1) + yt(:knon, k)) / ypaprs(:knon, k) &
                  * (ypplay(:knon, k - 1) - ypplay(:knon, k))
          ENDDO

          CALL cdrag(nisrf, sqrt(yu(:knon, 1)**2 + yv(:knon, 1)**2), &
               yt(:knon, 1), yq(:knon, 1), zgeop(:knon, 1), ypaprs(:knon, 1), &
               yts(:knon), yqsurf(:knon), yrugos(:knon), ycdragm(:knon), &
               ycdragh(:knon))

          IF (nisrf == is_oce) THEN
             ! On met un seuil pour ycdragm et ycdragh :
             ycdragm(:knon) = min(ycdragm(:knon), cdmmax)
             ycdragh(:knon) = min(ycdragh(:knon), cdhmax)
          END IF

          IF (iflag_pbl >= 6) yq2(:knon, :) = q2(ni(:knon), :, nisrf)
          call coef_diff_turb(nisrf, ni(:knon), ypaprs(:knon, :), &
               ypplay(:knon, :), yu(:knon, :), yv(:knon, :), yq(:knon, :), &
               yt(:knon, :), yts(:knon), ycdragm(:knon), zgeop(:knon, :), &
               ycoefm(:knon, :), ycoefh(:knon, :), yq2(:knon, :))
          CALL clvent(yu(:knon, 1), yv(:knon, 1), ycoefm(:knon, :), &
               ycdragm(:knon), yt(:knon, :), yu(:knon, :), ypaprs(:knon, :), &
               ypplay(:knon, :), ydelp(:knon, :), y_d_u(:knon, :), &
               y_flux_u(:knon))
          CALL clvent(yu(:knon, 1), yv(:knon, 1), ycoefm(:knon, :), &
               ycdragm(:knon), yt(:knon, :), yv(:knon, :), ypaprs(:knon, :), &
               ypplay(:knon, :), ydelp(:knon, :), y_d_v(:knon, :), &
               y_flux_v(:knon))
          CALL clqh(julien, nisrf, ni(:knon), ytsoil(:knon, :), yqsol(:knon), &
               mu0(ni(:knon)), yrugos(:knon), yrugoro(:knon), yu(:knon, 1), &
               yv(:knon, 1), ycoefh(:knon, :), ycdragh(:knon), yt(:knon, :), &
               yq(:knon, :), yts(:knon), ypaprs(:knon, :), ypplay(:knon, :), &
               ydelp(:knon, :), &
               fsolsw(ni(:knon), nisrf) + fsollw(ni(:knon), nisrf), &
               yalbedo(:knon), snow(:knon), yqsurf(:knon), yrain_fall(:knon), &
               ysnow_fall(:knon), yfluxlat(:knon), pctsrf_new_sic(ni(:knon)), &
               yagesno(:knon), y_d_t(:knon, :), y_d_q(:knon, :), &
               tsurf_new(:knon), yz0_new(:knon), y_flux_t(:knon), &
               y_flux_q(:knon), y_dflux_t(:knon), y_dflux_q(:knon), &
               y_fqcalving(:knon), y_ffonte(:knon), y_run_off_lic_0(:knon), &
               y_run_off_lic(:knon))

          ! calculer la longueur de rugosite sur ocean

          yrugm = 0.

          IF (nisrf == is_oce) THEN
             DO j = 1, knon
                yrugm(j) = 0.018 * ycdragm(j) * (yu(j, 1)**2 + yv(j, 1)**2) &
                     / rg + 0.11 * 14E-6 &
                     / sqrt(ycdragm(j) * (yu(j, 1)**2 + yv(j, 1)**2))
                yrugm(j) = max(1.5E-05, yrugm(j))
             END DO
          END IF

          DO k = 1, klev
             DO j = 1, knon
                y_d_t(j, k) = y_d_t(j, k) * ypctsrf(j)
                y_d_q(j, k) = y_d_q(j, k) * ypctsrf(j)
                y_d_u(j, k) = y_d_u(j, k) * ypctsrf(j)
                y_d_v(j, k) = y_d_v(j, k) * ypctsrf(j)
             END DO
          END DO

          flux_t(ni(:knon), nisrf) = y_flux_t(:knon)
          flux_q(ni(:knon), nisrf) = y_flux_q(:knon)
          flux_u(ni(:knon), nisrf) = y_flux_u(:knon)
          flux_v(ni(:knon), nisrf) = y_flux_v(:knon)

          falbe(:, nisrf) = 0.
          fsnow(:, nisrf) = 0.
          fqsurf(:, nisrf) = 0.
          frugs(:, nisrf) = 0.

          DO j = 1, knon
             i = ni(j)
             d_ts(i, nisrf) = tsurf_new(j) - yts(j)
             ftsol(i, nisrf) = tsurf_new(j) ! update surface temperature
             falbe(i, nisrf) = yalbedo(j)
             fsnow(i, nisrf) = snow(j)
             fqsurf(i, nisrf) = yqsurf(j)
             fluxlat(i, nisrf) = yfluxlat(j)

             IF (nisrf == is_oce) THEN
                rugmer(i) = yrugm(j)
                frugs(i, is_oce) = yrugm(j)
             else
                frugs(i, nisrf) = yz0_new(j)
             END IF

             agesno(i, nisrf) = yagesno(j)
             fqcalving(i, nisrf) = y_fqcalving(j)
             ffonte(i, nisrf) = y_ffonte(j)
             cdragh(i) = cdragh(i) + ycdragh(j) * ypctsrf(j)
             cdragm(i) = cdragm(i) + ycdragm(j) * ypctsrf(j)
             dflux_t(i) = dflux_t(i) + y_dflux_t(j) * ypctsrf(j)
             dflux_q(i) = dflux_q(i) + y_dflux_q(j) * ypctsrf(j)
          END DO

          IF (nisrf == is_ter) THEN
             qsol(ni(:knon)) = yqsol(:knon)
          else IF (nisrf == is_lic) THEN
             DO j = 1, knon
                i = ni(j)
                run_off_lic_0(i) = y_run_off_lic_0(j)
                run_off_lic(i) = y_run_off_lic(j)
             END DO
          END IF

          ftsoil(:, :, nisrf) = 0.
          ftsoil(ni(:knon), :, nisrf) = ytsoil(:knon, :)

          DO j = 1, knon
             i = ni(j)
             DO k = 1, klev
                d_t(i, k) = d_t(i, k) + y_d_t(j, k)
                d_q(i, k) = d_q(i, k) + y_d_q(j, k)
                d_u(i, k) = d_u(i, k) + y_d_u(j, k)
                d_v(i, k) = d_v(i, k) + y_d_v(j, k)
             END DO
          END DO

          forall (k = 2:klev) coefh(ni(:knon), k) &
               = coefh(ni(:knon), k) + ycoefh(:knon, k) * ypctsrf(:knon)

          ! Diagnostic temp\'erature, q \`a 2 m et u, v \`a 10 m:

          u1(:knon) = yu(:knon, 1) + y_d_u(:knon, 1)
          v1(:knon) = yv(:knon, 1) + y_d_v(:knon, 1)
          t1(:knon) = yt(:knon, 1) + y_d_t(:knon, 1)

          CALL stdlevvar(nisrf, u1(:knon), v1(:knon), t1(:knon), &
               yq(:knon, 1) + y_d_q(:knon, 1), &
               rd * t1(:knon) / (0.5 * (ypaprs(:knon, 1) + ypplay(:knon, 1))) &
               * (ypaprs(:knon, 1) - ypplay(:knon, 1)), &
               tsurf_new(:knon), yqsurf(:knon), merge(frugs(ni(:knon), &
               is_oce), yrugos(:knon), nisrf == is_oce), ypaprs(:knon, 1), &
               ypplay(:knon, 1), yt2m(:knon), yq2m(:knon), wind10m(:knon), &
               ustar(:knon))

          DO j = 1, knon
             i = ni(j)
             t2m(i, nisrf) = yt2m(j)
             q2m(i, nisrf) = yq2m(j)

             u10m_srf(i, nisrf) = (wind10m(j) * u1(j)) &
                  / sqrt(u1(j)**2 + v1(j)**2)
             v10m_srf(i, nisrf) = (wind10m(j) * v1(j)) &
                  / sqrt(u1(j)**2 + v1(j)**2)
          END DO

          CALL hbtm(ypaprs, ypplay, yt2m, yq2m, ustar(:knon), y_flux_t(:knon), &
               y_flux_q(:knon), yu(:knon, :), yv(:knon, :), yt(:knon, :), &
               yq(:knon, :), ypblh(:knon), ycapcl(:knon), yoliqcl(:knon), &
               ycteicl(:knon), ypblt(:knon), ytherm(:knon), ylcl(:knon))

          DO j = 1, knon
             i = ni(j)
             pblh(i, nisrf) = ypblh(j)
             lcl(i, nisrf) = ylcl(j)
             capcl(i, nisrf) = ycapcl(j)
             oliqcl(i, nisrf) = yoliqcl(j)
             cteicl(i, nisrf) = ycteicl(j)
             pblt(i, nisrf) = ypblt(j)
             therm(i, nisrf) = ytherm(j)
          END DO

          IF (iflag_pbl >= 6) q2(ni(:knon), :, nisrf) = yq2(:knon, :)
       else
          fsnow(:, nisrf) = 0.
       end IF if_knon
    END DO loop_surface

    ! On utilise les nouvelles surfaces
    frugs(:, is_oce) = rugmer
    pctsrf(:, is_oce) = pctsrf_new_oce
    pctsrf(:, is_sic) = pctsrf_new_sic

    CALL histwrite_phy("fder", dlw + dflux_t + dflux_q)
    CALL histwrite_phy("run_off_lic", run_off_lic)
    CALL histwrite_phy("dtsvdfo", d_ts(:, is_oce))
    CALL histwrite_phy("dtsvdft", d_ts(:, is_ter))
    CALL histwrite_phy("dtsvdfg", d_ts(:, is_lic))
    CALL histwrite_phy("dtsvdfi", d_ts(:, is_sic))
    CALL histwrite_phy("s_pblh", sum(pblh * pctsrf, dim = 2))
    CALL histwrite_phy("s_pblt", sum(pblT * pctsrf, dim = 2))
    CALL histwrite_phy("s_lcl", sum(lcl * pctsrf, dim = 2))
    CALL histwrite_phy("s_capCL", sum(capCL * pctsrf, dim = 2))
    CALL histwrite_phy("s_oliqCL", sum(oliqCL * pctsrf, dim = 2))
    CALL histwrite_phy("s_cteiCL", sum(cteiCL * pctsrf, dim = 2))
    CALL histwrite_phy("s_therm", sum(therm * pctsrf, dim = 2))

    tsol = sum(ftsol * pctsrf, dim = 2)
    dlw = - 4. * RSIGMA * tsol**3
    call assert(abs(sum(pctsrf, dim = 2) - 1.) <= EPSFRA, 'pbl_surface: pctsrf')

  END SUBROUTINE pbl_surface

end module pbl_surface_m
