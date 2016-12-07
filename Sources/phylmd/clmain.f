module clmain_m

  IMPLICIT NONE

contains

  SUBROUTINE clmain(dtime, pctsrf, t, q, u, v, jour, rmu0, ftsol, cdmmax, &
       cdhmax, ksta, ksta_ter, ok_kzmin, ftsoil, qsol, paprs, pplay, snow, &
       qsurf, evap, falbe, fluxlat, rain_fall, snow_f, solsw, sollw, fder, &
       rugos, agesno, rugoro, d_t, d_q, d_u, d_v, d_ts, flux_t, flux_q, &
       flux_u, flux_v, cdragh, cdragm, q2, dflux_t, dflux_q, ycoefh, zu1, &
       zv1, t2m, q2m, u10m, v10m, pblh, capcl, oliqcl, cteicl, pblt, therm, &
       trmb1, trmb2, trmb3, plcl, fqcalving, ffonte, run_off_lic_0)

    ! From phylmd/clmain.F, version 1.6, 2005/11/16 14:47:19
    ! Author: Z. X. Li (LMD/CNRS), date: 1993/08/18
    ! Objet : interface de couche limite (diffusion verticale)

    ! Tout ce qui a trait aux traceurs est dans "phytrac". Le calcul
    ! de la couche limite pour les traceurs se fait avec "cltrac" et
    ! ne tient pas compte de la diff\'erentiation des sous-fractions
    ! de sol.

    ! Pour pouvoir extraire les coefficients d'\'echanges et le vent
    ! dans la premi\`ere couche, trois champs ont \'et\'e cr\'e\'es : "ycoefh",
    ! "zu1" et "zv1". Nous avons moyenn\'e les valeurs de ces trois
    ! champs sur les quatre sous-surfaces du mod\`ele.

    use clqh_m, only: clqh
    use clvent_m, only: clvent
    use coefkz_m, only: coefkz
    use coefkzmin_m, only: coefkzmin
    USE conf_gcm_m, ONLY: prt_level, lmt_pas
    USE conf_phys_m, ONLY: iflag_pbl
    USE dimphy, ONLY: klev, klon, zmasq
    USE dimsoil, ONLY: nsoilmx
    use hbtm_m, only: hbtm
    USE indicesol, ONLY: epsfra, is_lic, is_oce, is_sic, is_ter, nbsrf
    USE interfoce_lim_m, ONLY: interfoce_lim
    use stdlevvar_m, only: stdlevvar
    USE suphec_m, ONLY: rd, rg, rkappa
    use time_phylmdz, only: itap
    use ustarhb_m, only: ustarhb
    use vdif_kcay_m, only: vdif_kcay
    use yamada4_m, only: yamada4

    REAL, INTENT(IN):: dtime ! interval du temps (secondes)

    REAL, INTENT(inout):: pctsrf(klon, nbsrf)
    ! tableau des pourcentages de surface de chaque maille

    REAL, INTENT(IN):: t(klon, klev) ! temperature (K)
    REAL, INTENT(IN):: q(klon, klev) ! vapeur d'eau (kg/kg)
    REAL, INTENT(IN):: u(klon, klev), v(klon, klev) ! vitesse
    INTEGER, INTENT(IN):: jour ! jour de l'annee en cours
    REAL, intent(in):: rmu0(klon) ! cosinus de l'angle solaire zenithal     
    REAL, INTENT(IN):: ftsol(klon, nbsrf) ! temp\'erature du sol (en K)
    REAL, INTENT(IN):: cdmmax, cdhmax ! seuils cdrm, cdrh
    REAL, INTENT(IN):: ksta, ksta_ter
    LOGICAL, INTENT(IN):: ok_kzmin

    REAL, INTENT(inout):: ftsoil(klon, nsoilmx, nbsrf)
    ! soil temperature of surface fraction

    REAL, INTENT(inout):: qsol(klon)
    ! column-density of water in soil, in kg m-2

    REAL, INTENT(IN):: paprs(klon, klev+1) ! pression a intercouche (Pa)
    REAL, INTENT(IN):: pplay(klon, klev) ! pression au milieu de couche (Pa)
    REAL, INTENT(inout):: snow(klon, nbsrf)
    REAL qsurf(klon, nbsrf)
    REAL evap(klon, nbsrf)
    REAL, intent(inout):: falbe(klon, nbsrf)

    REAL fluxlat(klon, nbsrf)

    REAL, intent(in):: rain_fall(klon)
    ! liquid water mass flux (kg/m2/s), positive down

    REAL, intent(in):: snow_f(klon)
    ! solid water mass flux (kg/m2/s), positive down

    REAL, INTENT(IN):: solsw(klon, nbsrf), sollw(klon, nbsrf)
    REAL, intent(in):: fder(klon)
    REAL, intent(inout):: rugos(klon, nbsrf) ! longueur de rugosit\'e (en m)
    real agesno(klon, nbsrf)
    REAL, INTENT(IN):: rugoro(klon)

    REAL d_t(klon, klev), d_q(klon, klev)
    ! d_t------output-R- le changement pour "t"
    ! d_q------output-R- le changement pour "q"

    REAL, intent(out):: d_u(klon, klev), d_v(klon, klev)
    ! changement pour "u" et "v"

    REAL, intent(out):: d_ts(klon, nbsrf) ! le changement pour ftsol

    REAL, intent(out):: flux_t(klon, nbsrf)
    ! flux de chaleur sensible (Cp T) (W/m2) (orientation positive vers
    ! le bas) à la surface

    REAL, intent(out):: flux_q(klon, nbsrf) 
    ! flux de vapeur d'eau (kg/m2/s) à la surface

    REAL, intent(out):: flux_u(klon, nbsrf), flux_v(klon, nbsrf)
    ! tension du vent à la surface, en Pa

    REAL, INTENT(out):: cdragh(klon), cdragm(klon)
    real q2(klon, klev+1, nbsrf)

    REAL, INTENT(out):: dflux_t(klon), dflux_q(klon)
    ! dflux_t derive du flux sensible
    ! dflux_q derive du flux latent
    ! IM "slab" ocean

    REAL, intent(out):: ycoefh(klon, klev)
    REAL, intent(out):: zu1(klon)
    REAL zv1(klon)
    REAL t2m(klon, nbsrf), q2m(klon, nbsrf)
    REAL u10m(klon, nbsrf), v10m(klon, nbsrf)

    ! Ionela Musat cf. Anne Mathieu : planetary boundary layer, hbtm
    ! (Comme les autres diagnostics on cumule dans physiq ce qui
    ! permet de sortir les grandeurs par sous-surface)
    REAL pblh(klon, nbsrf) ! height of planetary boundary layer
    REAL capcl(klon, nbsrf)
    REAL oliqcl(klon, nbsrf)
    REAL cteicl(klon, nbsrf)
    REAL pblt(klon, nbsrf)
    ! pblT------- T au nveau HCL
    REAL therm(klon, nbsrf)
    REAL trmb1(klon, nbsrf)
    ! trmb1-------deep_cape
    REAL trmb2(klon, nbsrf)
    ! trmb2--------inhibition
    REAL trmb3(klon, nbsrf)
    ! trmb3-------Point Omega
    REAL plcl(klon, nbsrf)
    REAL fqcalving(klon, nbsrf), ffonte(klon, nbsrf)
    ! ffonte----Flux thermique utilise pour fondre la neige
    ! fqcalving-Flux d'eau "perdue" par la surface et necessaire pour limiter la
    !           hauteur de neige, en kg/m2/s
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
    REAL yts(klon), yrugos(klon), ypct(klon), yz0_new(klon)
    REAL yalb(klon)
    REAL yu1(klon), yv1(klon)
    ! on rajoute en output yu1 et yv1 qui sont les vents dans
    ! la premiere couche
    REAL ysnow(klon), yqsurf(klon), yagesno(klon)

    real yqsol(klon)
    ! column-density of water in soil, in kg m-2

    REAL yrain_f(klon)
    ! liquid water mass flux (kg/m2/s), positive down

    REAL ysnow_f(klon)
    ! solid water mass flux (kg/m2/s), positive down

    REAL yfder(klon)
    REAL yrugm(klon), yrads(klon), yrugoro(klon)

    REAL yfluxlat(klon)

    REAL y_d_ts(klon)
    REAL y_d_t(klon, klev), y_d_q(klon, klev)
    REAL y_d_u(klon, klev), y_d_v(klon, klev)
    REAL y_flux_t(klon), y_flux_q(klon)
    REAL y_flux_u(klon), y_flux_v(klon)
    REAL y_dflux_t(klon), y_dflux_q(klon)
    REAL coefh(klon, klev), coefm(klon, klev)
    REAL yu(klon, klev), yv(klon, klev)
    REAL yt(klon, klev), yq(klon, klev)
    REAL ypaprs(klon, klev+1), ypplay(klon, klev), ydelp(klon, klev)

    REAL ycoefm0(klon, klev), ycoefh0(klon, klev)

    REAL yzlay(klon, klev), yzlev(klon, klev+1), yteta(klon, klev)
    REAL ykmm(klon, klev+1), ykmn(klon, klev+1)
    REAL ykmq(klon, klev+1)
    REAL yq2(klon, klev+1)
    REAL q2diag(klon, klev+1)

    REAL u1lay(klon), v1lay(klon)
    REAL delp(klon, klev)
    INTEGER i, k, nsrf

    INTEGER ni(klon), knon, j

    REAL pctsrf_pot(klon, nbsrf)
    ! "pourcentage potentiel" pour tenir compte des \'eventuelles
    ! apparitions ou disparitions de la glace de mer

    REAL zx_alf1, zx_alf2 ! valeur ambiante par extrapolation

    REAL yt2m(klon), yq2m(klon), yu10m(klon)
    REAL yustar(klon)

    REAL yt10m(klon), yq10m(klon)
    REAL ypblh(klon)
    REAL ylcl(klon)
    REAL ycapcl(klon)
    REAL yoliqcl(klon)
    REAL ycteicl(klon)
    REAL ypblt(klon)
    REAL ytherm(klon)
    REAL ytrmb1(klon)
    REAL ytrmb2(klon)
    REAL ytrmb3(klon)
    REAL uzon(klon), vmer(klon)
    REAL tair1(klon), qair1(klon), tairsol(klon)
    REAL psfce(klon), patm(klon)

    REAL qairsol(klon), zgeo1(klon)
    REAL rugo1(klon)

    ! utiliser un jeu de fonctions simples               
    LOGICAL zxli
    PARAMETER (zxli=.FALSE.)

    !------------------------------------------------------------

    ytherm = 0.

    DO k = 1, klev ! epaisseur de couche
       DO i = 1, klon
          delp(i, k) = paprs(i, k) - paprs(i, k+1)
       END DO
    END DO
    DO i = 1, klon ! vent de la premiere couche
       zx_alf1 = 1.0
       zx_alf2 = 1.0 - zx_alf1
       u1lay(i) = u(i, 1)*zx_alf1 + u(i, 2)*zx_alf2
       v1lay(i) = v(i, 1)*zx_alf1 + v(i, 2)*zx_alf2
    END DO

    ! Initialization:
    rugmer = 0.
    cdragh = 0.
    cdragm = 0.
    dflux_t = 0.
    dflux_q = 0.
    zu1 = 0.
    zv1 = 0.
    ypct = 0.
    yts = 0.
    ysnow = 0.
    yqsurf = 0.
    yrain_f = 0.
    ysnow_f = 0.
    yfder = 0.
    yrugos = 0.
    yu1 = 0.
    yv1 = 0.
    yrads = 0.
    ypaprs = 0.
    ypplay = 0.
    ydelp = 0.
    yu = 0.
    yv = 0.
    yt = 0.
    yq = 0.
    y_dflux_t = 0.
    y_dflux_q = 0.
    yrugoro = 0.
    d_ts = 0.
    yfluxlat = 0.
    flux_t = 0.
    flux_q = 0.
    flux_u = 0.
    flux_v = 0.
    d_t = 0.
    d_q = 0.
    d_u = 0.
    d_v = 0.
    ycoefh = 0.

    ! Initialisation des "pourcentages potentiels". On consid\`ere ici qu'on
    ! peut avoir potentiellement de la glace sur tout le domaine oc\'eanique
    ! (\`a affiner)

    pctsrf_pot(:, is_ter) = pctsrf(:, is_ter)
    pctsrf_pot(:, is_lic) = pctsrf(:, is_lic)
    pctsrf_pot(:, is_oce) = 1. - zmasq
    pctsrf_pot(:, is_sic) = 1. - zmasq

    ! Tester si c'est le moment de lire le fichier:
    if (mod(itap - 1, lmt_pas) == 0) then
       CALL interfoce_lim(jour, pctsrf_new_oce, pctsrf_new_sic)
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
             ysnow(j) = snow(i, nsrf)
             yqsurf(j) = qsurf(i, nsrf)
             yalb(j) = falbe(i, nsrf)
             yrain_f(j) = rain_fall(i)
             ysnow_f(j) = snow_f(i)
             yagesno(j) = agesno(i, nsrf)
             yfder(j) = fder(i)
             yrugos(j) = rugos(i, nsrf)
             yrugoro(j) = rugoro(i)
             yu1(j) = u1lay(i)
             yv1(j) = v1lay(i)
             yrads(j) = solsw(i, nsrf) + sollw(i, nsrf)
             ypaprs(j, klev+1) = paprs(i, klev+1)
             y_run_off_lic_0(j) = run_off_lic_0(i)
          END DO

          ! For continent, copy soil water content
          IF (nsrf == is_ter) THEN
             yqsol(:knon) = qsol(ni(:knon))
          ELSE
             yqsol = 0.
          END IF

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

          ! calculer Cdrag et les coefficients d'echange
          CALL coefkz(nsrf, ypaprs, ypplay, ksta, ksta_ter, yts, yrugos, yu, &
               yv, yt, yq, yqsurf, coefm(:knon, :), coefh(:knon, :))
          IF (iflag_pbl == 1) THEN
             CALL coefkz2(nsrf, knon, ypaprs, ypplay, yt, ycoefm0, ycoefh0)
             coefm(:knon, :) = max(coefm(:knon, :), ycoefm0(:knon, :))
             coefh(:knon, :) = max(coefh(:knon, :), ycoefh0(:knon, :))
          END IF

          ! on met un seuil pour coefm et coefh
          IF (nsrf == is_oce) THEN
             coefm(:knon, 1) = min(coefm(:knon, 1), cdmmax)
             coefh(:knon, 1) = min(coefh(:knon, 1), cdhmax)
          END IF

          IF (ok_kzmin) THEN
             ! Calcul d'une diffusion minimale pour les conditions tres stables
             CALL coefkzmin(knon, ypaprs, ypplay, yu, yv, yt, yq, &
                  coefm(:knon, 1), ycoefm0, ycoefh0)
             coefm(:knon, :) = max(coefm(:knon, :), ycoefm0(:knon, :))
             coefh(:knon, :) = max(coefh(:knon, :), ycoefh0(:knon, :))
          END IF

          IF (iflag_pbl >= 3) THEN
             ! Mellor et Yamada adapt\'e \`a Mars, Richard Fournier et
             ! Fr\'ed\'eric Hourdin
             yzlay(:knon, 1) = rd * yt(:knon, 1) / (0.5 * (ypaprs(:knon, 1) &
                  + ypplay(:knon, 1))) &
                  * (ypaprs(:knon, 1) - ypplay(:knon, 1)) / rg
             DO k = 2, klev
                yzlay(1:knon, k) = yzlay(1:knon, k-1) &
                     + rd * 0.5 * (yt(1:knon, k-1) + yt(1:knon, k)) &
                     / ypaprs(1:knon, k) &
                     * (ypplay(1:knon, k-1) - ypplay(1:knon, k)) / rg
             END DO
             DO k = 1, klev
                yteta(1:knon, k) = yt(1:knon, k)*(ypaprs(1:knon, 1) &
                     / ypplay(1:knon, k))**rkappa * (1.+0.61*yq(1:knon, k))
             END DO
             yzlev(1:knon, 1) = 0.
             yzlev(:knon, klev+1) = 2. * yzlay(:knon, klev) &
                  - yzlay(:knon, klev - 1)
             DO k = 2, klev
                yzlev(1:knon, k) = 0.5*(yzlay(1:knon, k)+yzlay(1:knon, k-1))
             END DO
             DO k = 1, klev + 1
                DO j = 1, knon
                   i = ni(j)
                   yq2(j, k) = q2(i, k, nsrf)
                END DO
             END DO

             CALL ustarhb(knon, yu, yv, coefm(:knon, 1), yustar)
             IF (prt_level > 9) PRINT *, 'USTAR = ', yustar

             ! iflag_pbl peut \^etre utilis\'e comme longueur de m\'elange

             IF (iflag_pbl >= 11) THEN
                CALL vdif_kcay(knon, dtime, rg, ypaprs, yzlev, yzlay, yu, yv, &
                     yteta, coefm(:knon, 1), yq2, q2diag, ykmm, ykmn, yustar, &
                     iflag_pbl)
             ELSE
                CALL yamada4(knon, dtime, rg, yzlev, yzlay, yu, yv, yteta, &
                     coefm(:knon, 1), yq2, ykmm, ykmn, ykmq, yustar, iflag_pbl)
             END IF

             coefm(:knon, 2:) = ykmm(:knon, 2:klev)
             coefh(:knon, 2:) = ykmn(:knon, 2:klev)
          END IF

          ! calculer la diffusion des vitesses "u" et "v"
          CALL clvent(knon, dtime, yu1, yv1, coefm(:knon, :), yt, yu, ypaprs, &
               ypplay, ydelp, y_d_u, y_flux_u(:knon))
          CALL clvent(knon, dtime, yu1, yv1, coefm(:knon, :), yt, yv, ypaprs, &
               ypplay, ydelp, y_d_v, y_flux_v(:knon))

          ! calculer la diffusion de "q" et de "h"
          CALL clqh(dtime, jour, firstcal, nsrf, ni(:knon), ytsoil(:knon, :), &
               yqsol, rmu0, yrugos, yrugoro, yu1, yv1, coefh(:knon, :), yt, &
               yq, yts(:knon), ypaprs, ypplay, ydelp, yrads, yalb(:knon), &
               ysnow, yqsurf, yrain_f, ysnow_f, yfder, yfluxlat, &
               pctsrf_new_sic, yagesno(:knon), y_d_t, y_d_q, y_d_ts(:knon), &
               yz0_new, y_flux_t(:knon), y_flux_q(:knon), y_dflux_t, &
               y_dflux_q, y_fqcalving, y_ffonte, y_run_off_lic_0)

          ! calculer la longueur de rugosite sur ocean
          yrugm = 0.
          IF (nsrf == is_oce) THEN
             DO j = 1, knon
                yrugm(j) = 0.018*coefm(j, 1)*(yu1(j)**2+yv1(j)**2)/rg + &
                     0.11*14E-6/sqrt(coefm(j, 1)*(yu1(j)**2+yv1(j)**2))
                yrugm(j) = max(1.5E-05, yrugm(j))
             END DO
          END IF
          DO j = 1, knon
             y_dflux_t(j) = y_dflux_t(j)*ypct(j)
             y_dflux_q(j) = y_dflux_q(j)*ypct(j)
             yu1(j) = yu1(j)*ypct(j)
             yv1(j) = yv1(j)*ypct(j)
          END DO

          DO k = 1, klev
             DO j = 1, knon
                i = ni(j)
                coefh(j, k) = coefh(j, k)*ypct(j)
                coefm(j, k) = coefm(j, k)*ypct(j)
                y_d_t(j, k) = y_d_t(j, k)*ypct(j)
                y_d_q(j, k) = y_d_q(j, k)*ypct(j)
                y_d_u(j, k) = y_d_u(j, k)*ypct(j)
                y_d_v(j, k) = y_d_v(j, k)*ypct(j)
             END DO
          END DO

          DO j = 1, knon
             i = ni(j)
             flux_t(i, nsrf) = y_flux_t(j)
             flux_q(i, nsrf) = y_flux_q(j)
             flux_u(i, nsrf) = y_flux_u(j)
             flux_v(i, nsrf) = y_flux_v(j)
          END DO

          evap(:, nsrf) = -flux_q(:, nsrf)

          falbe(:, nsrf) = 0.
          snow(:, nsrf) = 0.
          qsurf(:, nsrf) = 0.
          rugos(:, nsrf) = 0.
          fluxlat(:, nsrf) = 0.
          DO j = 1, knon
             i = ni(j)
             d_ts(i, nsrf) = y_d_ts(j)
             falbe(i, nsrf) = yalb(j)
             snow(i, nsrf) = ysnow(j)
             qsurf(i, nsrf) = yqsurf(j)
             rugos(i, nsrf) = yz0_new(j)
             fluxlat(i, nsrf) = yfluxlat(j)
             IF (nsrf == is_oce) THEN
                rugmer(i) = yrugm(j)
                rugos(i, nsrf) = yrugm(j)
             END IF
             agesno(i, nsrf) = yagesno(j)
             fqcalving(i, nsrf) = y_fqcalving(j)
             ffonte(i, nsrf) = y_ffonte(j)
             cdragh(i) = cdragh(i) + coefh(j, 1)
             cdragm(i) = cdragm(i) + coefm(j, 1)
             dflux_t(i) = dflux_t(i) + y_dflux_t(j)
             dflux_q(i) = dflux_q(i) + y_dflux_q(j)
             zu1(i) = zu1(i) + yu1(j)
             zv1(i) = zv1(i) + yv1(j)
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
                ycoefh(i, k) = ycoefh(i, k) + coefh(j, k)
             END DO
          END DO

          ! diagnostic t, q a 2m et u, v a 10m

          DO j = 1, knon
             i = ni(j)
             uzon(j) = yu(j, 1) + y_d_u(j, 1)
             vmer(j) = yv(j, 1) + y_d_v(j, 1)
             tair1(j) = yt(j, 1) + y_d_t(j, 1)
             qair1(j) = yq(j, 1) + y_d_q(j, 1)
             zgeo1(j) = rd*tair1(j)/(0.5*(ypaprs(j, 1)+ypplay(j, &
                  1)))*(ypaprs(j, 1)-ypplay(j, 1))
             tairsol(j) = yts(j) + y_d_ts(j)
             rugo1(j) = yrugos(j)
             IF (nsrf == is_oce) THEN
                rugo1(j) = rugos(i, nsrf)
             END IF
             psfce(j) = ypaprs(j, 1)
             patm(j) = ypplay(j, 1)

             qairsol(j) = yqsurf(j)
          END DO

          CALL stdlevvar(klon, knon, nsrf, zxli, uzon, vmer, tair1, qair1, &
               zgeo1, tairsol, qairsol, rugo1, psfce, patm, yt2m, yq2m, &
               yt10m, yq10m, yu10m, yustar)

          DO j = 1, knon
             i = ni(j)
             t2m(i, nsrf) = yt2m(j)
             q2m(i, nsrf) = yq2m(j)

             ! u10m, v10m : composantes du vent a 10m sans spirale de Ekman
             u10m(i, nsrf) = (yu10m(j)*uzon(j))/sqrt(uzon(j)**2+vmer(j)**2)
             v10m(i, nsrf) = (yu10m(j)*vmer(j))/sqrt(uzon(j)**2+vmer(j)**2)
          END DO

          CALL hbtm(ypaprs, ypplay, yt2m, yq2m, yustar, y_flux_t(:knon), &
               y_flux_q(:knon), yu, yv, yt, yq, ypblh(:knon), ycapcl, &
               yoliqcl, ycteicl, ypblt, ytherm, ytrmb1, ytrmb2, ytrmb3, ylcl)

          DO j = 1, knon
             i = ni(j)
             pblh(i, nsrf) = ypblh(j)
             plcl(i, nsrf) = ylcl(j)
             capcl(i, nsrf) = ycapcl(j)
             oliqcl(i, nsrf) = yoliqcl(j)
             cteicl(i, nsrf) = ycteicl(j)
             pblt(i, nsrf) = ypblt(j)
             therm(i, nsrf) = ytherm(j)
             trmb1(i, nsrf) = ytrmb1(j)
             trmb2(i, nsrf) = ytrmb2(j)
             trmb3(i, nsrf) = ytrmb3(j)
          END DO

          DO j = 1, knon
             DO k = 1, klev + 1
                i = ni(j)
                q2(i, k, nsrf) = yq2(j, k)
             END DO
          END DO
       end IF if_knon
    END DO loop_surface

    ! On utilise les nouvelles surfaces
    rugos(:, is_oce) = rugmer
    pctsrf(:, is_oce) = pctsrf_new_oce
    pctsrf(:, is_sic) = pctsrf_new_sic

    firstcal = .false.

  END SUBROUTINE clmain

end module clmain_m
