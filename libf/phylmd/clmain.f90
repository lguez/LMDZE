module clmain_m

  IMPLICIT NONE

contains

  SUBROUTINE clmain(dtime, itap, date0, pctsrf, pctsrf_new, t, q, u, v,&
       jour, rmu0, co2_ppm, ok_veget, ocean, npas, nexca, ts,&
       soil_model, cdmmax, cdhmax, ksta, ksta_ter, ok_kzmin, ftsoil,&
       qsol, paprs, pplay, snow, qsurf, evap, albe, alblw, fluxlat,&
       rain_f, snow_f, solsw, sollw, sollwdown, fder, rlon, rlat, cufi,&
       cvfi, rugos, debut, lafin, agesno, rugoro, d_t, d_q, d_u, d_v,&
       d_ts, flux_t, flux_q, flux_u, flux_v, cdragh, cdragm, q2,&
       dflux_t, dflux_q, zcoefh, zu1, zv1, t2m, q2m, u10m, v10m, pblh,&
       capcl, oliqcl, cteicl, pblt, therm, trmb1, trmb2, trmb3, plcl,&
       fqcalving, ffonte, run_off_lic_0, flux_o, flux_g, tslab, seaice)

    ! From phylmd/clmain.F, version 1.6 2005/11/16 14:47:19
    ! Author: Z.X. Li (LMD/CNRS), date: 1993/08/18
    ! Objet : interface de "couche limite" (diffusion verticale)

    ! Tout ce qui a trait aux traceurs est dans "phytrac" maintenant.
    ! Pour l'instant le calcul de la couche limite pour les traceurs
    ! se fait avec "cltrac" et ne tient pas compte de la différentiation
    ! des sous-fractions de sol.

    ! Pour pouvoir extraire les coefficients d'échanges et le vent
    ! dans la première couche, trois champs supplémentaires ont été
    ! créés : "zcoefh", "zu1" et "zv1". Pour l'instant nous avons
    ! moyenné les valeurs de ces trois champs sur les 4 sous-surfaces
    ! du modèle. Dans l'avenir, si les informations des sous-surfaces
    ! doivent être prises en compte, il faudra sortir ces mêmes champs
    ! en leur ajoutant une dimension, c'est-à-dire "nbsrf" (nombre de
    ! sous-surfaces).

    use calendar, ONLY : ymds2ju
    use clqh_m, only: clqh
    use coefkz_m, only: coefkz
    use coefkzmin_m, only: coefkzmin
    USE conf_phys_m, ONLY : iflag_pbl
    USE dimens_m, ONLY : iim, jjm
    USE dimphy, ONLY : klev, klon, zmasq
    USE dimsoil, ONLY : nsoilmx
    USE dynetat0_m, ONLY : day_ini
    USE gath_cpl, ONLY : gath2cpl
    use hbtm_m, only: hbtm
    USE histcom, ONLY : histbeg_totreg, histdef, histend, histsync
    use histwrite_m, only: histwrite
    USE indicesol, ONLY : epsfra, is_lic, is_oce, is_sic, is_ter, nbsrf
    USE conf_gcm_m, ONLY : prt_level
    USE suphec_m, ONLY : rd, rg, rkappa
    USE temps, ONLY : annee_ref, itau_phy
    use yamada4_m, only: yamada4

    ! Arguments:

    REAL, INTENT (IN) :: dtime ! interval du temps (secondes)
    REAL date0
    ! date0----input-R- jour initial
    INTEGER, INTENT (IN) :: itap
    ! itap-----input-I- numero du pas de temps
    REAL, INTENT(IN):: t(klon, klev), q(klon, klev)
    ! t--------input-R- temperature (K)
    ! q--------input-R- vapeur d'eau (kg/kg)
    REAL, INTENT (IN):: u(klon, klev), v(klon, klev)
    ! u--------input-R- vitesse u
    ! v--------input-R- vitesse v
    REAL, INTENT (IN):: paprs(klon, klev+1)
    ! paprs----input-R- pression a intercouche (Pa)
    REAL, INTENT (IN):: pplay(klon, klev)
    ! pplay----input-R- pression au milieu de couche (Pa)
    REAL, INTENT (IN):: rlon(klon), rlat(klon)
    ! rlat-----input-R- latitude en degree
    REAL cufi(klon), cvfi(klon)
    ! cufi-----input-R- resolution des mailles en x (m)
    ! cvfi-----input-R- resolution des mailles en y (m)
    REAL d_t(klon, klev), d_q(klon, klev)
    ! d_t------output-R- le changement pour "t"
    ! d_q------output-R- le changement pour "q"
    REAL d_u(klon, klev), d_v(klon, klev)
    ! d_u------output-R- le changement pour "u"
    ! d_v------output-R- le changement pour "v"
    REAL flux_t(klon, klev, nbsrf), flux_q(klon, klev, nbsrf)
    ! flux_t---output-R- flux de chaleur sensible (CpT) J/m**2/s (W/m**2)
    !                    (orientation positive vers le bas)
    ! flux_q---output-R- flux de vapeur d'eau (kg/m**2/s)
    REAL dflux_t(klon), dflux_q(klon)
    ! dflux_t derive du flux sensible
    ! dflux_q derive du flux latent
    !IM "slab" ocean
    REAL flux_o(klon), flux_g(klon)
    !IM "slab" ocean
    ! flux_g---output-R-  flux glace (pour OCEAN='slab  ')
    ! flux_o---output-R-  flux ocean (pour OCEAN='slab  ')
    REAL y_flux_o(klon), y_flux_g(klon)
    REAL tslab(klon), ytslab(klon)
    ! tslab-in/output-R temperature du slab ocean (en Kelvin) 
    ! uniqmnt pour slab
    REAL seaice(klon), y_seaice(klon)
    ! seaice---output-R-  glace de mer (kg/m2) (pour OCEAN='slab  ')
    REAL y_fqcalving(klon), y_ffonte(klon)
    REAL fqcalving(klon, nbsrf), ffonte(klon, nbsrf)
    ! ffonte----Flux thermique utilise pour fondre la neige
    ! fqcalving-Flux d'eau "perdue" par la surface et necessaire pour limiter la
    !           hauteur de neige, en kg/m2/s
    REAL run_off_lic_0(klon), y_run_off_lic_0(klon)

    REAL flux_u(klon, klev, nbsrf), flux_v(klon, klev, nbsrf)
    ! flux_u---output-R- tension du vent X: (kg m/s)/(m**2 s) ou Pascal
    ! flux_v---output-R- tension du vent Y: (kg m/s)/(m**2 s) ou Pascal
    REAL rugmer(klon), agesno(klon, nbsrf)
    REAL, INTENT (IN) :: rugoro(klon)
    REAL cdragh(klon), cdragm(klon)
    ! jour de l'annee en cours                
    INTEGER jour
    REAL rmu0(klon) ! cosinus de l'angle solaire zenithal     
    ! taux CO2 atmosphere                     
    REAL co2_ppm
    LOGICAL, INTENT (IN) :: debut
    LOGICAL, INTENT (IN) :: lafin
    LOGICAL ok_veget
    CHARACTER (len=*), INTENT (IN) :: ocean
    INTEGER npas, nexca

    REAL pctsrf(klon, nbsrf)
    REAL ts(klon, nbsrf)
    ! ts-------input-R- temperature du sol (en Kelvin)
    REAL d_ts(klon, nbsrf)
    ! d_ts-----output-R- le changement pour "ts"
    REAL snow(klon, nbsrf)
    REAL qsurf(klon, nbsrf)
    REAL evap(klon, nbsrf)
    REAL albe(klon, nbsrf)
    REAL alblw(klon, nbsrf)

    REAL fluxlat(klon, nbsrf)

    REAL rain_f(klon), snow_f(klon)
    REAL fder(klon)

    REAL sollw(klon, nbsrf), solsw(klon, nbsrf), sollwdown(klon)
    REAL rugos(klon, nbsrf)
    ! rugos----input-R- longeur de rugosite (en m)
    ! la nouvelle repartition des surfaces sortie de l'interface
    REAL pctsrf_new(klon, nbsrf)

    REAL zcoefh(klon, klev)
    REAL zu1(klon)
    REAL zv1(klon)

    !$$$ PB ajout pour soil
    LOGICAL, INTENT (IN) :: soil_model
    !IM ajout seuils cdrm, cdrh
    REAL cdmmax, cdhmax

    REAL ksta, ksta_ter
    LOGICAL ok_kzmin

    REAL ftsoil(klon, nsoilmx, nbsrf)
    REAL ytsoil(klon, nsoilmx)
    REAL qsol(klon)

    EXTERNAL clvent, calbeta, cltrac

    REAL yts(klon), yrugos(klon), ypct(klon), yz0_new(klon)
    REAL yalb(klon)
    REAL yalblw(klon)
    REAL yu1(klon), yv1(klon)
    ! on rajoute en output yu1 et yv1 qui sont les vents dans
    ! la premiere couche
    REAL ysnow(klon), yqsurf(klon), yagesno(klon), yqsol(klon)
    REAL yrain_f(klon), ysnow_f(klon)
    REAL ysollw(klon), ysolsw(klon), ysollwdown(klon)
    REAL yfder(klon), ytaux(klon), ytauy(klon)
    REAL yrugm(klon), yrads(klon), yrugoro(klon)

    REAL yfluxlat(klon)

    REAL y_d_ts(klon)
    REAL y_d_t(klon, klev), y_d_q(klon, klev)
    REAL y_d_u(klon, klev), y_d_v(klon, klev)
    REAL y_flux_t(klon, klev), y_flux_q(klon, klev)
    REAL y_flux_u(klon, klev), y_flux_v(klon, klev)
    REAL y_dflux_t(klon), y_dflux_q(klon)
    REAL ycoefh(klon, klev), ycoefm(klon, klev)
    REAL yu(klon, klev), yv(klon, klev)
    REAL yt(klon, klev), yq(klon, klev)
    REAL ypaprs(klon, klev+1), ypplay(klon, klev), ydelp(klon, klev)

    LOGICAL ok_nonloc
    PARAMETER (ok_nonloc=.FALSE.)
    REAL ycoefm0(klon, klev), ycoefh0(klon, klev)

    REAL yzlay(klon, klev), yzlev(klon, klev+1), yteta(klon, klev)
    REAL ykmm(klon, klev+1), ykmn(klon, klev+1)
    REAL ykmq(klon, klev+1)
    REAL yq2(klon, klev+1), q2(klon, klev+1, nbsrf)
    REAL q2diag(klon, klev+1)

    REAL u1lay(klon), v1lay(klon)
    REAL delp(klon, klev)
    INTEGER i, k, nsrf

    INTEGER ni(klon), knon, j

    REAL pctsrf_pot(klon, nbsrf)
    ! "pourcentage potentiel" pour tenir compte des éventuelles
    ! apparitions ou disparitions de la glace de mer

    REAL zx_alf1, zx_alf2 !valeur ambiante par extrapola.

    ! maf pour sorties IOISPL en cas de debugagage

    CHARACTER (80) cldebug
    SAVE cldebug
    CHARACTER (8) cl_surf(nbsrf)
    SAVE cl_surf
    INTEGER nhoridbg, nidbg
    SAVE nhoridbg, nidbg
    INTEGER ndexbg(iim*(jjm+1))
    REAL zx_lon(iim, jjm+1), zx_lat(iim, jjm+1), zjulian
    REAL tabindx(klon)
    REAL debugtab(iim, jjm+1)
    LOGICAL first_appel
    SAVE first_appel
    DATA first_appel/ .TRUE./
    LOGICAL :: debugindex = .FALSE.
    INTEGER idayref
    REAL t2m(klon, nbsrf), q2m(klon, nbsrf)
    REAL u10m(klon, nbsrf), v10m(klon, nbsrf)

    REAL yt2m(klon), yq2m(klon), yu10m(klon)
    REAL yustar(klon)
    ! -- LOOP
    REAL yu10mx(klon)
    REAL yu10my(klon)
    REAL ywindsp(klon)
    ! -- LOOP

    REAL yt10m(klon), yq10m(klon)
    !IM cf. AM : pbl, hbtm (Comme les autres diagnostics on cumule ds
    ! physiq ce qui permet de sortir les grdeurs par sous surface)
    REAL pblh(klon, nbsrf)
    ! pblh------- HCL
    REAL plcl(klon, nbsrf)
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
    REAL y_cd_h(klon), y_cd_m(klon)
    REAL uzon(klon), vmer(klon)
    REAL tair1(klon), qair1(klon), tairsol(klon)
    REAL psfce(klon), patm(klon)

    REAL qairsol(klon), zgeo1(klon)
    REAL rugo1(klon)

    ! utiliser un jeu de fonctions simples               
    LOGICAL zxli
    PARAMETER (zxli=.FALSE.)

    REAL zt, zqs, zdelta, zcor
    REAL t_coup
    PARAMETER (t_coup=273.15)

    CHARACTER (len=20) :: modname = 'clmain'

    !------------------------------------------------------------

    ytherm = 0.

    IF (debugindex .AND. first_appel) THEN
       first_appel = .FALSE.

       ! initialisation sorties netcdf

       idayref = day_ini
       CALL ymds2ju(annee_ref, 1, idayref, 0., zjulian)
       CALL gr_fi_ecrit(1, klon, iim, jjm+1, rlon, zx_lon)
       DO i = 1, iim
          zx_lon(i, 1) = rlon(i+1)
          zx_lon(i, jjm+1) = rlon(i+1)
       END DO
       CALL gr_fi_ecrit(1, klon, iim, jjm+1, rlat, zx_lat)
       cldebug = 'sous_index'
       CALL histbeg_totreg(cldebug, zx_lon(:, 1), zx_lat(1, :), 1, &
            iim, 1, jjm+1, itau_phy, zjulian, dtime, nhoridbg, nidbg)
       ! no vertical axis
       cl_surf(1) = 'ter'
       cl_surf(2) = 'lic'
       cl_surf(3) = 'oce'
       cl_surf(4) = 'sic'
       DO nsrf = 1, nbsrf
          CALL histdef(nidbg, cl_surf(nsrf), cl_surf(nsrf), '-', iim, jjm+1, &
               nhoridbg, 1, 1, 1, -99, 'inst', dtime, dtime)
       END DO
       CALL histend(nidbg)
       CALL histsync(nidbg)
    END IF

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
    yalb = 0.
    yalblw = 0.
    yrain_f = 0.
    ysnow_f = 0.
    yfder = 0.
    ytaux = 0.
    ytauy = 0.
    ysolsw = 0.
    ysollw = 0.
    ysollwdown = 0.
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
    pctsrf_new = 0.
    y_flux_u = 0.
    y_flux_v = 0.
    !$$ PB
    y_dflux_t = 0.
    y_dflux_q = 0.
    ytsoil = 999999.
    yrugoro = 0.
    ! -- LOOP
    yu10mx = 0.
    yu10my = 0.
    ywindsp = 0.
    ! -- LOOP
    d_ts = 0.
    !§§§ PB
    yfluxlat = 0.
    flux_t = 0.
    flux_q = 0.
    flux_u = 0.
    flux_v = 0.
    d_t = 0.
    d_q = 0.
    d_u = 0.
    d_v = 0.
    zcoefh = 0.

    ! Boucler sur toutes les sous-fractions du sol:

    ! Initialisation des "pourcentages potentiels". On considère ici qu'on
    ! peut avoir potentiellement de la glace sur tout le domaine océanique
    ! (à affiner)

    pctsrf_pot = pctsrf
    pctsrf_pot(:, is_oce) = 1. - zmasq
    pctsrf_pot(:, is_sic) = 1. - zmasq

    loop_surface: DO nsrf = 1, nbsrf
       ! Chercher les indices :
       ni = 0
       knon = 0
       DO i = 1, klon
          ! Pour déterminer le domaine à traiter, on utilise les surfaces
          ! "potentielles"
          IF (pctsrf_pot(i, nsrf) > epsfra) THEN
             knon = knon + 1
             ni(knon) = i
          END IF
       END DO

       ! variables pour avoir une sortie IOIPSL des INDEX
       IF (debugindex) THEN
          tabindx = 0.
          DO i = 1, knon
             tabindx(i) = real(i)
          END DO
          debugtab = 0.
          ndexbg = 0
          CALL gath2cpl(tabindx, debugtab, klon, knon, iim, jjm, ni)
          CALL histwrite(nidbg, cl_surf(nsrf), itap, debugtab)
       END IF

       IF (knon == 0) CYCLE

       DO j = 1, knon
          i = ni(j)
          ypct(j) = pctsrf(i, nsrf)
          yts(j) = ts(i, nsrf)
          ytslab(i) = tslab(i)
          ysnow(j) = snow(i, nsrf)
          yqsurf(j) = qsurf(i, nsrf)
          yalb(j) = albe(i, nsrf)
          yalblw(j) = alblw(i, nsrf)
          yrain_f(j) = rain_f(i)
          ysnow_f(j) = snow_f(i)
          yagesno(j) = agesno(i, nsrf)
          yfder(j) = fder(i)
          ytaux(j) = flux_u(i, 1, nsrf)
          ytauy(j) = flux_v(i, 1, nsrf)
          ysolsw(j) = solsw(i, nsrf)
          ysollw(j) = sollw(i, nsrf)
          ysollwdown(j) = sollwdown(i)
          yrugos(j) = rugos(i, nsrf)
          yrugoro(j) = rugoro(i)
          yu1(j) = u1lay(i)
          yv1(j) = v1lay(i)
          yrads(j) = ysolsw(j) + ysollw(j)
          ypaprs(j, klev+1) = paprs(i, klev+1)
          y_run_off_lic_0(j) = run_off_lic_0(i)
          yu10mx(j) = u10m(i, nsrf)
          yu10my(j) = v10m(i, nsrf)
          ywindsp(j) = sqrt(yu10mx(j)*yu10mx(j)+yu10my(j)*yu10my(j))
       END DO

       ! IF bucket model for continent, copy soil water content
       IF (nsrf == is_ter .AND. .NOT. ok_veget) THEN
          DO j = 1, knon
             i = ni(j)
             yqsol(j) = qsol(i)
          END DO
       ELSE
          yqsol = 0.
       END IF
       !$$$ PB ajour pour soil
       DO k = 1, nsoilmx
          DO j = 1, knon
             i = ni(j)
             ytsoil(j, k) = ftsoil(i, k, nsrf)
          END DO
       END DO
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
       CALL coefkz(nsrf, knon, ypaprs, ypplay, ksta, ksta_ter, yts,&
            yrugos, yu, yv, yt, yq, yqsurf, ycoefm, ycoefh)
       IF (iflag_pbl == 1) THEN
          CALL coefkz2(nsrf, knon, ypaprs, ypplay, yt, ycoefm0, ycoefh0)
          DO k = 1, klev
             DO i = 1, knon
                ycoefm(i, k) = max(ycoefm(i, k), ycoefm0(i, k))
                ycoefh(i, k) = max(ycoefh(i, k), ycoefh0(i, k))
             END DO
          END DO
       END IF

       ! on seuille ycoefm et ycoefh
       IF (nsrf == is_oce) THEN
          DO j = 1, knon
             ycoefm(j, 1) = min(ycoefm(j, 1), cdmmax)
             ycoefh(j, 1) = min(ycoefh(j, 1), cdhmax)
          END DO
       END IF

       IF (ok_kzmin) THEN
          ! Calcul d'une diffusion minimale pour les conditions tres stables
          CALL coefkzmin(knon, ypaprs, ypplay, yu, yv, yt, yq, ycoefm(:, 1), &
               ycoefm0, ycoefh0)

          DO k = 1, klev
             DO i = 1, knon
                ycoefm(i, k) = max(ycoefm(i, k), ycoefm0(i, k))
                ycoefh(i, k) = max(ycoefh(i, k), ycoefh0(i, k))
             END DO
          END DO
       END IF

       IF (iflag_pbl >= 3) THEN
          ! MELLOR ET YAMADA adapté à Mars, Richard Fournier et Frédéric Hourdin
          yzlay(1:knon, 1) = rd*yt(1:knon, 1)/(0.5*(ypaprs(1:knon, &
               1)+ypplay(1:knon, 1)))*(ypaprs(1:knon, 1)-ypplay(1:knon, 1))/rg
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
          yzlev(1:knon, klev+1) = 2.*yzlay(1:knon, klev) - yzlay(1:knon, klev-1)
          DO k = 2, klev
             yzlev(1:knon, k) = 0.5*(yzlay(1:knon, k)+yzlay(1:knon, k-1))
          END DO
          DO k = 1, klev + 1
             DO j = 1, knon
                i = ni(j)
                yq2(j, k) = q2(i, k, nsrf)
             END DO
          END DO

          y_cd_m(1:knon) = ycoefm(1:knon, 1)
          y_cd_h(1:knon) = ycoefh(1:knon, 1)
          CALL ustarhb(knon, yu, yv, y_cd_m, yustar)

          IF (prt_level>9) THEN
             PRINT *, 'USTAR = ', yustar
          END IF

          ! iflag_pbl peut être utilisé comme longueur de mélange

          IF (iflag_pbl >= 11) THEN
             CALL vdif_kcay(knon, dtime, rg, rd, ypaprs, yt, yzlev, yzlay, &
                  yu, yv, yteta, y_cd_m, yq2, q2diag, ykmm, ykmn, yustar, &
                  iflag_pbl)
          ELSE
             CALL yamada4(knon, dtime, rg, yzlev, yzlay, yu, yv, yteta, &
                  y_cd_m, yq2, ykmm, ykmn, ykmq, yustar, iflag_pbl)
          END IF

          ycoefm(1:knon, 1) = y_cd_m(1:knon)
          ycoefh(1:knon, 1) = y_cd_h(1:knon)
          ycoefm(1:knon, 2:klev) = ykmm(1:knon, 2:klev)
          ycoefh(1:knon, 2:klev) = ykmn(1:knon, 2:klev)
       END IF

       ! calculer la diffusion des vitesses "u" et "v"
       CALL clvent(knon, dtime, yu1, yv1, ycoefm, yt, yu, ypaprs, ypplay, &
            ydelp, y_d_u, y_flux_u)
       CALL clvent(knon, dtime, yu1, yv1, ycoefm, yt, yv, ypaprs, ypplay, &
            ydelp, y_d_v, y_flux_v)

       ! pour le couplage
       ytaux = y_flux_u(:, 1)
       ytauy = y_flux_v(:, 1)

       ! calculer la diffusion de "q" et de "h"
       CALL clqh(dtime, itap, date0, jour, debut, lafin, rlon, rlat,&
            cufi, cvfi, knon, nsrf, ni, pctsrf, soil_model, ytsoil,&
            yqsol, ok_veget, ocean, npas, nexca, rmu0, co2_ppm, yrugos,&
            yrugoro, yu1, yv1, ycoefh, yt, yq, yts, ypaprs, ypplay,&
            ydelp, yrads, yalb, yalblw, ysnow, yqsurf, yrain_f, ysnow_f, &
            yfder, ytaux, ytauy, ywindsp, ysollw, ysollwdown, ysolsw,&
            yfluxlat, pctsrf_new, yagesno, y_d_t, y_d_q, y_d_ts,&
            yz0_new, y_flux_t, y_flux_q, y_dflux_t, y_dflux_q,&
            y_fqcalving, y_ffonte, y_run_off_lic_0, y_flux_o, y_flux_g,&
            ytslab, y_seaice)

       ! calculer la longueur de rugosite sur ocean
       yrugm = 0.
       IF (nsrf == is_oce) THEN
          DO j = 1, knon
             yrugm(j) = 0.018*ycoefm(j, 1)*(yu1(j)**2+yv1(j)**2)/rg + &
                  0.11*14E-6/sqrt(ycoefm(j, 1)*(yu1(j)**2+yv1(j)**2))
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
             ycoefh(j, k) = ycoefh(j, k)*ypct(j)
             ycoefm(j, k) = ycoefm(j, k)*ypct(j)
             y_d_t(j, k) = y_d_t(j, k)*ypct(j)
             y_d_q(j, k) = y_d_q(j, k)*ypct(j)
             flux_t(i, k, nsrf) = y_flux_t(j, k)
             flux_q(i, k, nsrf) = y_flux_q(j, k)
             flux_u(i, k, nsrf) = y_flux_u(j, k)
             flux_v(i, k, nsrf) = y_flux_v(j, k)
             y_d_u(j, k) = y_d_u(j, k)*ypct(j)
             y_d_v(j, k) = y_d_v(j, k)*ypct(j)
          END DO
       END DO

       evap(:, nsrf) = -flux_q(:, 1, nsrf)

       albe(:, nsrf) = 0.
       alblw(:, nsrf) = 0.
       snow(:, nsrf) = 0.
       qsurf(:, nsrf) = 0.
       rugos(:, nsrf) = 0.
       fluxlat(:, nsrf) = 0.
       DO j = 1, knon
          i = ni(j)
          d_ts(i, nsrf) = y_d_ts(j)
          albe(i, nsrf) = yalb(j)
          alblw(i, nsrf) = yalblw(j)
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
          cdragh(i) = cdragh(i) + ycoefh(j, 1)
          cdragm(i) = cdragm(i) + ycoefm(j, 1)
          dflux_t(i) = dflux_t(i) + y_dflux_t(j)
          dflux_q(i) = dflux_q(i) + y_dflux_q(j)
          zu1(i) = zu1(i) + yu1(j)
          zv1(i) = zv1(i) + yv1(j)
       END DO
       IF (nsrf == is_ter) THEN
          DO j = 1, knon
             i = ni(j)
             qsol(i) = yqsol(j)
          END DO
       END IF
       IF (nsrf == is_lic) THEN
          DO j = 1, knon
             i = ni(j)
             run_off_lic_0(i) = y_run_off_lic_0(j)
          END DO
       END IF
       !$$$ PB ajout pour soil
       ftsoil(:, :, nsrf) = 0.
       DO k = 1, nsoilmx
          DO j = 1, knon
             i = ni(j)
             ftsoil(i, k, nsrf) = ytsoil(j, k)
          END DO
       END DO

       DO j = 1, knon
          i = ni(j)
          DO k = 1, klev
             d_t(i, k) = d_t(i, k) + y_d_t(j, k)
             d_q(i, k) = d_q(i, k) + y_d_q(j, k)
             d_u(i, k) = d_u(i, k) + y_d_u(j, k)
             d_v(i, k) = d_v(i, k) + y_d_v(j, k)
             zcoefh(i, k) = zcoefh(i, k) + ycoefh(j, k)
          END DO
       END DO

       !cc diagnostic t, q a 2m et u, v a 10m

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

       CALL stdlevvar(klon, knon, nsrf, zxli, uzon, vmer, tair1, qair1, zgeo1, &
            tairsol, qairsol, rugo1, psfce, patm, yt2m, yq2m, yt10m, yq10m, &
            yu10m, yustar)

       DO j = 1, knon
          i = ni(j)
          t2m(i, nsrf) = yt2m(j)
          q2m(i, nsrf) = yq2m(j)

          ! u10m, v10m : composantes du vent a 10m sans spirale de Ekman
          u10m(i, nsrf) = (yu10m(j)*uzon(j))/sqrt(uzon(j)**2+vmer(j)**2)
          v10m(i, nsrf) = (yu10m(j)*vmer(j))/sqrt(uzon(j)**2+vmer(j)**2)

       END DO

       DO i = 1, knon
          y_cd_h(i) = ycoefh(i, 1)
          y_cd_m(i) = ycoefm(i, 1)
       END DO
       CALL hbtm(knon, ypaprs, ypplay, yt2m, yt10m, yq2m, yq10m, yustar, &
            y_flux_t, y_flux_q, yu, yv, yt, yq, ypblh, ycapcl, yoliqcl, &
            ycteicl, ypblt, ytherm, ytrmb1, ytrmb2, ytrmb3, ylcl)

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
       !IM "slab" ocean
       IF (nsrf == is_oce) THEN
          DO j = 1, knon
             ! on projette sur la grille globale
             i = ni(j)
             IF (pctsrf_new(i, is_oce)>epsfra) THEN
                flux_o(i) = y_flux_o(j)
             ELSE
                flux_o(i) = 0.
             END IF
          END DO
       END IF

       IF (nsrf == is_sic) THEN
          DO j = 1, knon
             i = ni(j)
             ! On pondère lorsque l'on fait le bilan au sol :
             IF (pctsrf_new(i, is_sic)>epsfra) THEN
                flux_g(i) = y_flux_g(j)
             ELSE
                flux_g(i) = 0.
             END IF
          END DO

       END IF
       IF (ocean == 'slab  ') THEN
          IF (nsrf == is_oce) THEN
             tslab(1:klon) = ytslab(1:klon)
             seaice(1:klon) = y_seaice(1:klon)
          END IF
       END IF
    END DO loop_surface

    ! On utilise les nouvelles surfaces

    rugos(:, is_oce) = rugmer
    pctsrf = pctsrf_new

  END SUBROUTINE clmain

end module clmain_m
