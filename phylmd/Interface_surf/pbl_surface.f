module pbl_surface_m

  IMPLICIT NONE

contains

  SUBROUTINE pbl_surface(pctsrf, t, q, u, v, julien, mu0, ftsol, cdmmax, &
       cdhmax, ftsoil, qsol, paprs, pplay, fsnow, qsurf, evap, falbe, fluxlat, &
       rain_fall, snow_f, fsolsw, fsollw, frugs, agesno, rugoro, d_t, d_q, &
       d_u, d_v, d_ts, flux_t, flux_q, flux_u, flux_v, cdragh, cdragm, q2, &
       dflux_t, dflux_q, coefh, t2m, q2m, u10m_srf, v10m_srf, pblh, capcl, &
       oliqcl, cteicl, pblt, therm, plcl, fqcalving, ffonte, run_off_lic_0)

    ! From phylmd/clmain.F, version 1.6, 2005/11/16 14:47:19
    ! Author: Z. X. Li (LMD/CNRS), date: 1993 Aug. 18th
    ! Objet : interface de couche limite (diffusion verticale)

    ! Tout ce qui a trait aux traceurs est dans "phytrac". Le calcul
    ! de la couche limite pour les traceurs se fait avec "cltrac" et
    ! ne tient pas compte de la diff\'erentiation des sous-fractions
    ! de sol.

    use cdrag_m, only: cdrag
    use clqh_m, only: clqh
    use clvent_m, only: clvent
    use coef_diff_turb_m, only: coef_diff_turb
    USE conf_gcm_m, ONLY: lmt_pas
    USE conf_phys_m, ONLY: iflag_pbl
    USE dimphy, ONLY: klev, klon
    USE dimsoil, ONLY: nsoilmx
    use hbtm_m, only: hbtm
    USE indicesol, ONLY: epsfra, is_lic, is_oce, is_sic, is_ter, nbsrf
    USE interfoce_lim_m, ONLY: interfoce_lim
    use phyetat0_m, only: zmasq
    use stdlevvar_m, only: stdlevvar
    USE suphec_m, ONLY: rd, rg
    use time_phylmdz, only: itap

    REAL, INTENT(inout):: pctsrf(klon, nbsrf)
    ! tableau des pourcentages de surface de chaque maille

    REAL, INTENT(IN):: t(klon, klev) ! temperature (K)
    REAL, INTENT(IN):: q(klon, klev) ! vapeur d'eau (kg / kg)
    REAL, INTENT(IN):: u(klon, klev), v(klon, klev) ! vitesse
    INTEGER, INTENT(IN):: julien ! jour de l'annee en cours
    REAL, intent(in):: mu0(klon) ! cosinus de l'angle solaire zenithal     
    REAL, INTENT(IN):: ftsol(:, :) ! (klon, nbsrf) temp\'erature du sol (en K)
    REAL, INTENT(IN):: cdmmax, cdhmax ! seuils cdrm, cdrh

    REAL, INTENT(inout):: ftsoil(klon, nsoilmx, nbsrf)
    ! soil temperature of surface fraction

    REAL, INTENT(inout):: qsol(:) ! (klon)
    ! column-density of water in soil, in kg m-2

    REAL, INTENT(IN):: paprs(klon, klev + 1) ! pression a intercouche (Pa)
    REAL, INTENT(IN):: pplay(klon, klev) ! pression au milieu de couche (Pa)
    REAL, INTENT(inout):: fsnow(:, :) ! (klon, nbsrf) \'epaisseur neigeuse
    REAL qsurf(klon, nbsrf)
    REAL evap(klon, nbsrf)
    REAL, intent(inout):: falbe(klon, nbsrf)
    REAL, intent(out):: fluxlat(:, :) ! (klon, nbsrf)

    REAL, intent(in):: rain_fall(klon)
    ! liquid water mass flux (kg / m2 / s), positive down

    REAL, intent(in):: snow_f(klon)
    ! solid water mass flux (kg / m2 / s), positive down

    REAL, INTENT(IN):: fsolsw(klon, nbsrf), fsollw(klon, nbsrf)
    REAL, intent(inout):: frugs(klon, nbsrf) ! longueur de rugosit\'e (en m)
    real agesno(klon, nbsrf)
    REAL, INTENT(IN):: rugoro(klon)

    REAL, intent(out):: d_t(:, :), d_q(:, :) ! (klon, klev)
    ! changement pour t et q

    REAL, intent(out):: d_u(klon, klev), d_v(klon, klev)
    ! changement pour "u" et "v"

    REAL, intent(out):: d_ts(:, :) ! (klon, nbsrf) variation of ftsol

    REAL, intent(out):: flux_t(klon, nbsrf)
    ! flux de chaleur sensible (Cp T) (W / m2) (orientation positive vers
    ! le bas) à la surface

    REAL, intent(out):: flux_q(klon, nbsrf) 
    ! flux de vapeur d'eau (kg / m2 / s) à la surface

    REAL, intent(out):: flux_u(klon, nbsrf), flux_v(klon, nbsrf)
    ! tension du vent (flux turbulent de vent) à la surface, en Pa

    REAL, INTENT(out):: cdragh(klon), cdragm(klon)
    real q2(klon, klev + 1, nbsrf)

    REAL, INTENT(out):: dflux_t(klon), dflux_q(klon)
    ! dflux_t derive du flux sensible
    ! dflux_q derive du flux latent
    ! IM "slab" ocean

    REAL, intent(out):: coefh(:, 2:) ! (klon, 2:klev)
    ! Pour pouvoir extraire les coefficients d'\'echange, le champ
    ! "coefh" a \'et\'e cr\'e\'e. Nous avons moyenn\'e les valeurs de
    ! ce champ sur les quatre sous-surfaces du mod\`ele.

    REAL, INTENT(inout):: t2m(klon, nbsrf), q2m(klon, nbsrf)

    REAL, INTENT(inout):: u10m_srf(:, :), v10m_srf(:, :) ! (klon, nbsrf)
    ! composantes du vent \`a 10m sans spirale d'Ekman

    ! Ionela Musat. Cf. Anne Mathieu : planetary boundary layer, hbtm.
    ! Comme les autres diagnostics on cumule dans physiq ce qui permet
    ! de sortir les grandeurs par sous-surface.
    REAL pblh(klon, nbsrf) ! height of planetary boundary layer
    REAL capcl(klon, nbsrf)
    REAL oliqcl(klon, nbsrf)
    REAL cteicl(klon, nbsrf)
    REAL, INTENT(inout):: pblt(klon, nbsrf) ! T au nveau HCL
    REAL therm(klon, nbsrf)
    REAL plcl(klon, nbsrf)

    REAL, intent(out):: fqcalving(klon, nbsrf)
    ! flux d'eau "perdue" par la surface et necessaire pour limiter la
    ! hauteur de neige, en kg / m2 / s

    real ffonte(klon, nbsrf)
    ! ffonte----Flux thermique utilise pour fondre la neige
    REAL run_off_lic_0(klon)

    ! Local:

    LOGICAL:: firstcal = .true.

    ! la nouvelle repartition des surfaces sortie de l'interface
    REAL, save:: pctsrf_new_oce(klon)
    REAL, save:: pctsrf_new_sic(klon)

    REAL y_fqcalving(klon), y_ffonte(klon)
    real y_run_off_lic_0(klon)
    REAL rugmer(klon)
    REAL ytsoil(klon, nsoilmx)
    REAL yts(klon), ypct(klon), yz0_new(klon)
    real yrugos(klon) ! longueur de rugosite (en m)
    REAL yalb(klon)
    REAL snow(klon), yqsurf(klon), yagesno(klon)
    real yqsol(klon) ! column-density of water in soil, in kg m-2
    REAL yrain_f(klon) ! liquid water mass flux (kg / m2 / s), positive down
    REAL ysnow_f(klon) ! solid water mass flux (kg / m2 / s), positive down
    REAL yrugm(klon), yrads(klon), yrugoro(klon)
    REAL yfluxlat(klon)
    REAL y_d_ts(klon)
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
    REAL delp(klon, klev)
    INTEGER i, k, nsrf
    INTEGER ni(klon), knon, j

    REAL pctsrf_pot(klon, nbsrf)
    ! "pourcentage potentiel" pour tenir compte des \'eventuelles
    ! apparitions ou disparitions de la glace de mer

    REAL yt2m(klon), yq2m(klon), wind10m(klon)
    REAL ustar(klon)

    REAL yt10m(klon), yq10m(klon)
    REAL ypblh(klon)
    REAL ylcl(klon)
    REAL ycapcl(klon)
    REAL yoliqcl(klon)
    REAL ycteicl(klon)
    REAL ypblt(klon)
    REAL ytherm(klon)
    REAL u1(klon), v1(klon)
    REAL tair1(klon), qair1(klon), tairsol(klon)
    REAL psfce(klon), patm(klon)

    REAL qairsol(klon), zgeo1(klon)
    REAL rugo1(klon)
    REAL zgeop(klon, klev)

    !------------------------------------------------------------

    ytherm = 0.

    DO k = 1, klev ! epaisseur de couche
       DO i = 1, klon
          delp(i, k) = paprs(i, k) - paprs(i, k + 1)
       END DO
    END DO

    ! Initialization:
    rugmer = 0.
    cdragh = 0.
    cdragm = 0.
    dflux_t = 0.
    dflux_q = 0.
    ypct = 0.
    yqsurf = 0.
    yrain_f = 0.
    ysnow_f = 0.
    yrugos = 0.
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

    ! Initialisation des "pourcentages potentiels". On consid\`ere ici qu'on
    ! peut avoir potentiellement de la glace sur tout le domaine oc\'eanique
    ! (\`a affiner)

    pctsrf_pot(:, is_ter) = pctsrf(:, is_ter)
    pctsrf_pot(:, is_lic) = pctsrf(:, is_lic)
    pctsrf_pot(:, is_oce) = 1. - zmasq
    pctsrf_pot(:, is_sic) = 1. - zmasq

    ! Tester si c'est le moment de lire le fichier:
    if (mod(itap - 1, lmt_pas) == 0) then
       CALL interfoce_lim(julien, pctsrf_new_oce, pctsrf_new_sic)
    endif

    ! Boucler sur toutes les sous-fractions du sol:

    loop_surface: DO nsrf = 1, nbsrf
       ! Chercher les indices :
       ni = 0
       knon = 0
       DO i = 1, klon
          ! Pour d\'eterminer le domaine \`a traiter, on utilise les surfaces
          ! "potentielles"
          IF (pctsrf_pot(i, nsrf) > epsfra) THEN
             knon = knon + 1
             ni(knon) = i
          END IF
       END DO

       if_knon: IF (knon /= 0) then
          DO j = 1, knon
             i = ni(j)
             ypct(j) = pctsrf(i, nsrf)
             yts(j) = ftsol(i, nsrf)
             snow(j) = fsnow(i, nsrf)
             yqsurf(j) = qsurf(i, nsrf)
             yalb(j) = falbe(i, nsrf)
             yrain_f(j) = rain_fall(i)
             ysnow_f(j) = snow_f(i)
             yagesno(j) = agesno(i, nsrf)
             yrugos(j) = frugs(i, nsrf)
             yrugoro(j) = rugoro(i)
             yrads(j) = fsolsw(i, nsrf) + fsollw(i, nsrf)
             ypaprs(j, klev + 1) = paprs(i, klev + 1)
             y_run_off_lic_0(j) = run_off_lic_0(i)
          END DO

          ! For continent, copy soil water content
          IF (nsrf == is_ter) yqsol(:knon) = qsol(ni(:knon))

          ytsoil(:knon, :) = ftsoil(ni(:knon), :, nsrf)

          DO k = 1, klev
             DO j = 1, knon
                i = ni(j)
                ypaprs(j, k) = paprs(i, k)
                ypplay(j, k) = pplay(i, k)
                ydelp(j, k) = delp(i, k)
                yu(j, k) = u(i, k)
                yv(j, k) = v(i, k)
                yt(j, k) = t(i, k)
                yq(j, k) = q(i, k)
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

          CALL cdrag(nsrf, sqrt(yu(:knon, 1)**2 + yv(:knon, 1)**2), &
               yt(:knon, 1), yq(:knon, 1), zgeop(:knon, 1), ypaprs(:knon, 1), &
               yts(:knon), yqsurf(:knon), yrugos(:knon), ycdragm(:knon), &
               ycdragh(:knon)) 

          IF (iflag_pbl == 1) THEN
             ycdragm(:knon) = max(ycdragm(:knon), 0.)
             ycdragh(:knon) = max(ycdragh(:knon), 0.)
          end IF

          ! on met un seuil pour ycdragm et ycdragh
          IF (nsrf == is_oce) THEN
             ycdragm(:knon) = min(ycdragm(:knon), cdmmax)
             ycdragh(:knon) = min(ycdragh(:knon), cdhmax)
          END IF

          IF (iflag_pbl >= 6) then
             DO k = 1, klev + 1
                DO j = 1, knon
                   i = ni(j)
                   yq2(j, k) = q2(i, k, nsrf)
                END DO
             END DO
          end IF

          call coef_diff_turb(nsrf, ni(:knon), ypaprs(:knon, :), &
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

          CALL clqh(julien, firstcal, nsrf, ni(:knon), ytsoil(:knon, :), &
               yqsol(:knon), mu0, yrugos(:knon), yrugoro(:knon), yu(:knon, 1), &
               yv(:knon, 1), ycoefh(:knon, :), ycdragh(:knon), yt(:knon, :), &
               yq(:knon, :), yts(:knon), ypaprs(:knon, :), ypplay(:knon, :), &
               ydelp(:knon, :), yrads(:knon), yalb(:knon), snow(:knon), &
               yqsurf(:knon), yrain_f, ysnow_f, yfluxlat(:knon), &
               pctsrf_new_sic, yagesno(:knon), y_d_t(:knon, :), &
               y_d_q(:knon, :), y_d_ts(:knon), yz0_new(:knon), &
               y_flux_t(:knon), y_flux_q(:knon), y_dflux_t(:knon), &
               y_dflux_q(:knon), y_fqcalving(:knon), y_ffonte, y_run_off_lic_0)

          ! calculer la longueur de rugosite sur ocean

          yrugm = 0.

          IF (nsrf == is_oce) THEN
             DO j = 1, knon
                yrugm(j) = 0.018 * ycdragm(j) * (yu(j, 1)**2 + yv(j, 1)**2) &
                     / rg + 0.11 * 14E-6 &
                     / sqrt(ycdragm(j) * (yu(j, 1)**2 + yv(j, 1)**2))
                yrugm(j) = max(1.5E-05, yrugm(j))
             END DO
          END IF

          DO k = 1, klev
             DO j = 1, knon
                i = ni(j)
                y_d_t(j, k) = y_d_t(j, k) * ypct(j)
                y_d_q(j, k) = y_d_q(j, k) * ypct(j)
                y_d_u(j, k) = y_d_u(j, k) * ypct(j)
                y_d_v(j, k) = y_d_v(j, k) * ypct(j)
             END DO
          END DO

          flux_t(ni(:knon), nsrf) = y_flux_t(:knon)
          flux_q(ni(:knon), nsrf) = y_flux_q(:knon)
          flux_u(ni(:knon), nsrf) = y_flux_u(:knon)
          flux_v(ni(:knon), nsrf) = y_flux_v(:knon)

          evap(:, nsrf) = -flux_q(:, nsrf)

          falbe(:, nsrf) = 0.
          fsnow(:, nsrf) = 0.
          qsurf(:, nsrf) = 0.
          frugs(:, nsrf) = 0.
          DO j = 1, knon
             i = ni(j)
             d_ts(i, nsrf) = y_d_ts(j)
             falbe(i, nsrf) = yalb(j)
             fsnow(i, nsrf) = snow(j)
             qsurf(i, nsrf) = yqsurf(j)
             frugs(i, nsrf) = yz0_new(j)
             fluxlat(i, nsrf) = yfluxlat(j)
             IF (nsrf == is_oce) THEN
                rugmer(i) = yrugm(j)
                frugs(i, nsrf) = yrugm(j)
             END IF
             agesno(i, nsrf) = yagesno(j)
             fqcalving(i, nsrf) = y_fqcalving(j)
             ffonte(i, nsrf) = y_ffonte(j)
             cdragh(i) = cdragh(i) + ycdragh(j) * ypct(j)
             cdragm(i) = cdragm(i) + ycdragm(j) * ypct(j)
             dflux_t(i) = dflux_t(i) + y_dflux_t(j) * ypct(j)
             dflux_q(i) = dflux_q(i) + y_dflux_q(j) * ypct(j)
          END DO
          IF (nsrf == is_ter) THEN
             qsol(ni(:knon)) = yqsol(:knon)
          else IF (nsrf == is_lic) THEN
             DO j = 1, knon
                i = ni(j)
                run_off_lic_0(i) = y_run_off_lic_0(j)
             END DO
          END IF

          ftsoil(:, :, nsrf) = 0.
          ftsoil(ni(:knon), :, nsrf) = ytsoil(:knon, :)

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
               = coefh(ni(:knon), k) + ycoefh(:knon, k) * ypct(:knon)

          ! diagnostic t, q a 2m et u, v a 10m

          DO j = 1, knon
             i = ni(j)
             u1(j) = yu(j, 1) + y_d_u(j, 1)
             v1(j) = yv(j, 1) + y_d_v(j, 1)
             tair1(j) = yt(j, 1) + y_d_t(j, 1)
             qair1(j) = yq(j, 1) + y_d_q(j, 1)
             zgeo1(j) = rd * tair1(j) / (0.5 * (ypaprs(j, 1) + ypplay(j, &
                  1))) * (ypaprs(j, 1)-ypplay(j, 1))
             tairsol(j) = yts(j) + y_d_ts(j)
             rugo1(j) = yrugos(j)
             IF (nsrf == is_oce) THEN
                rugo1(j) = frugs(i, nsrf)
             END IF
             psfce(j) = ypaprs(j, 1)
             patm(j) = ypplay(j, 1)

             qairsol(j) = yqsurf(j)
          END DO

          CALL stdlevvar(nsrf, u1(:knon), v1(:knon), tair1(:knon), qair1, &
               zgeo1, tairsol, qairsol, rugo1, psfce, patm, yt2m, yq2m, yt10m, &
               yq10m, wind10m(:knon), ustar(:knon))

          DO j = 1, knon
             i = ni(j)
             t2m(i, nsrf) = yt2m(j)
             q2m(i, nsrf) = yq2m(j)

             u10m_srf(i, nsrf) = (wind10m(j) * u1(j)) &
                  / sqrt(u1(j)**2 + v1(j)**2)
             v10m_srf(i, nsrf) = (wind10m(j) * v1(j)) &
                  / sqrt(u1(j)**2 + v1(j)**2)
          END DO

          CALL hbtm(ypaprs, ypplay, yt2m, yq2m, ustar(:knon), y_flux_t(:knon), &
               y_flux_q(:knon), yu(:knon, :), yv(:knon, :), yt(:knon, :), &
               yq(:knon, :), ypblh(:knon), ycapcl, yoliqcl, ycteicl, ypblt, &
               ytherm, ylcl)

          DO j = 1, knon
             i = ni(j)
             pblh(i, nsrf) = ypblh(j)
             plcl(i, nsrf) = ylcl(j)
             capcl(i, nsrf) = ycapcl(j)
             oliqcl(i, nsrf) = yoliqcl(j)
             cteicl(i, nsrf) = ycteicl(j)
             pblt(i, nsrf) = ypblt(j)
             therm(i, nsrf) = ytherm(j)
          END DO

          DO j = 1, knon
             DO k = 1, klev + 1
                i = ni(j)
                q2(i, k, nsrf) = yq2(j, k)
             END DO
          END DO
       else
          fsnow(:, nsrf) = 0.
       end IF if_knon
    END DO loop_surface

    ! On utilise les nouvelles surfaces
    frugs(:, is_oce) = rugmer
    pctsrf(:, is_oce) = pctsrf_new_oce
    pctsrf(:, is_sic) = pctsrf_new_sic

    firstcal = .false.

  END SUBROUTINE pbl_surface

end module pbl_surface_m
