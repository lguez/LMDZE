module etat0_mod

  use indicesol, only: nbsrf
  use dimphy, only: klon

  IMPLICIT NONE

  REAL pctsrf(klon, nbsrf)

  private nbsrf, klon

contains

  SUBROUTINE etat0

    ! From "etat0_netcdf.F", version 1.3 2005/05/25 13:10:09

    ! This subroutine creates "masque".

    USE ioipsl, only: flinget, flinclo, flinopen_nozoom, flininfo, histclo

    USE start_init_orog_m, only: start_init_orog, masque, phis
    use start_init_phys_m, only: qsol_2d
    use startdyn, only: start_inter_3d, start_init_dyn
    use dimens_m, only: iim, jjm, llm, nqmx
    use paramet_m, only: ip1jm, ip1jmp1
    use comconst, only: dtvr, daysec, cpp, kappa, pi
    use comdissnew, only: lstardis, nitergdiv, nitergrot, niterh, &
         tetagdiv, tetagrot, tetatemp
    use indicesol, only: is_oce, is_sic, is_ter, is_lic, epsfra
    use comvert, only: ap, bp, preff, pa
    use dimphy, only: zmasq
    use conf_gcm_m, only: day_step, iphysiq, dayref, anneeref
    use comgeom, only: rlatu, rlonv, rlonu, rlatv, aire_2d, apoln, apols, &
         cu_2d, cv_2d
    use serre, only: alphax
    use dimsoil, only: nsoilmx
    use temps, only: itau_dyn, itau_phy, annee_ref, day_ref, dt
    use clesphys2, only: ok_orodr, nbapp_rad
    use grid_atob, only: grille_m
    use grid_change, only: init_dyn_phy, dyn_phy
    use q_sat_m, only: q_sat
    use exner_hyb_m, only: exner_hyb
    use advtrac_m, only: iniadvtrac
    use pressure_var, only: pls, p3d
    use dynredem0_m, only: dynredem0
    use regr_lat_time_coefoz_m, only: regr_lat_time_coefoz
    use regr_pr_o3_m, only: regr_pr_o3

    ! Variables local to the procedure:

    REAL latfi(klon), lonfi(klon)
    ! (latitude and longitude of a point of the scalar grid identified
    ! by a simple index, in °)

    REAL, dimension(iim + 1, jjm + 1, llm):: uvent, t3d, tpot
    REAL vvent(iim + 1, jjm, llm)

    REAL q3d(iim + 1, jjm + 1, llm, nqmx)
    ! (mass fractions of trace species
    ! "q3d(i, j, l)" is at longitude "rlonv(i)", latitude "rlatu(j)"
    ! and pressure level "pls(i, j, l)".)

    real qsat(iim + 1, jjm + 1, llm) ! mass fraction of saturating water vapor
    REAL tsol(klon), qsol(klon), sn(klon)
    REAL tsolsrf(klon, nbsrf), qsolsrf(klon, nbsrf), snsrf(klon, nbsrf) 
    REAL albe(klon, nbsrf), evap(klon, nbsrf)
    REAL alblw(klon, nbsrf)
    REAL tsoil(klon, nsoilmx, nbsrf) 
    REAL radsol(klon), rain_fall(klon), snow_fall(klon)
    REAL solsw(klon), sollw(klon), fder(klon)
    !IM "slab" ocean
    REAL tslab(klon)
    real seaice(klon) ! kg m-2
    REAL frugs(klon, nbsrf), agesno(klon, nbsrf)
    REAL rugmer(klon)
    real, dimension(iim + 1, jjm + 1):: relief, zstd_2d, zsig_2d, zgam_2d
    real, dimension(iim + 1, jjm + 1):: zthe_2d, zpic_2d, zval_2d
    real, dimension(iim + 1, jjm + 1):: tsol_2d, psol
    REAL zmea(klon), zstd(klon)
    REAL zsig(klon), zgam(klon)
    REAL zthe(klon)
    REAL zpic(klon), zval(klon)
    REAL rugsrel(klon)
    REAL t_ancien(klon, llm), q_ancien(klon, llm)      !
    REAL run_off_lic_0(klon)
    real clwcon(klon, llm), rnebcon(klon, llm), ratqs(klon, llm)
    ! déclarations pour lecture glace de mer
    INTEGER iml_lic, jml_lic, llm_tmp, ttm_tmp
    INTEGER itaul(1), fid
    REAL lev(1), date
    REAL, ALLOCATABLE:: lon_lic(:, :), lat_lic(:, :)
    REAL, ALLOCATABLE:: dlon_lic(:), dlat_lic(:)
    REAL, ALLOCATABLE:: fraclic(:, :) ! fraction land ice
    REAL flic_tmp(iim + 1, jjm + 1) !fraction land ice temporary

    INTEGER l, ji

    REAL pk(iim + 1, jjm + 1, llm) ! fonction d'Exner aux milieux des couches 
    real pks(iim + 1, jjm + 1)

    REAL masse(iim + 1, jjm + 1, llm)
    REAL phi(ip1jmp1, llm)
    REAL pbaru(ip1jmp1, llm), pbarv(ip1jm, llm)
    REAL w(ip1jmp1, llm)
    REAL phystep
    INTEGER radpas

    !---------------------------------

    print *, "Call sequence information: etat0"

    ! Construct a grid:

    dtvr = daysec / real(day_step)
    print *, 'dtvr = ', dtvr

    pa = 5e4
    CALL iniconst
    CALL inigeom
    CALL inifilr

    latfi(1) = 90.
    latfi(2:klon-1) = pack(spread(rlatu(2:jjm), 1, iim), .true.) * 180. / pi
    ! (with conversion to degrees)
    latfi(klon) = - 90.

    lonfi(1) = 0.
    lonfi(2:klon-1) = pack(spread(rlonv(:iim), 2, jjm - 1), .true.) * 180. / pi
    ! (with conversion to degrees)
    lonfi(klon) = 0.

    call start_init_orog(relief, zstd_2d, zsig_2d, zgam_2d, zthe_2d, zpic_2d, &
         zval_2d) ! also compute "masque" and "phis"
    call init_dyn_phy ! define the mask "dyn_phy" for distinct grid points
    zmasq = pack(masque, dyn_phy)
    PRINT *, 'Masque construit'

    CALL start_init_dyn(tsol_2d, psol) ! also compute "qsol_2d"

    ! Compute pressure on intermediate levels:
    forall(l = 1: llm + 1) p3d(:, :, l) = ap(l) + bp(l) * psol(:, :)
    CALL exner_hyb(psol, p3d, pks, pk)
    IF (MINVAL(pk) == MAXVAL(pk)) stop '"pk" should not be constant'

    pls(:, :, :) = preff * (pk(:, :, :) / cpp)**(1. / kappa)
    PRINT *, "minval(pls(:, :, :)) = ", minval(pls(:, :, :))
    print *, "maxval(pls(:, :, :)) = ", maxval(pls(:, :, :))

    uvent(:, :, :) = start_inter_3d('U', rlonv, rlatv, pls)
    forall (l = 1: llm) uvent(:iim, :, l) = uvent(:iim, :, l) * cu_2d(:iim, :)
    uvent(iim+1, :, :) = uvent(1, :, :)

    vvent(:, :, :) = start_inter_3d('V', rlonu, rlatu(:jjm), pls(:, :jjm, :))
    forall (l = 1: llm) vvent(:iim, :, l) = vvent(:iim, :, l) * cv_2d(:iim, :)
    vvent(iim + 1, :, :) = vvent(1, :, :)

    t3d(:, :, :) = start_inter_3d('TEMP', rlonu, rlatv, pls)
    PRINT *,  'minval(t3d(:, :, :)) = ', minval(t3d(:, :, :))
    print *, "maxval(t3d(:, :, :)) = ", maxval(t3d(:, :, :))

    tpot(:iim, :, :) = t3d(:iim, :, :) * cpp / pk(:iim, :, :)
    tpot(iim + 1, :, :) = tpot(1, :, :)
    DO l=1, llm
       tpot(:, 1, l) = SUM(aire_2d(:, 1) * tpot(:, 1, l)) / apoln
       tpot(:, jjm + 1, l) = SUM(aire_2d(:, jjm + 1) * tpot(:, jjm + 1, l)) &
            / apols
    ENDDO

    ! Calcul de l'humidité à saturation :
    qsat(:, :, :) = q_sat(t3d, pls)
    PRINT *, "minval(qsat(:, :, :)) = ", minval(qsat(:, :, :))
    print *, "maxval(qsat(:, :, :)) = ", maxval(qsat(:, :, :))
    IF (MINVAL(qsat) == MAXVAL(qsat)) stop '"qsat" should not be constant'

    ! Water vapor:
    q3d(:, :, :, 1) = 0.01 * start_inter_3d('R', rlonu, rlatv, pls) * qsat
    WHERE (q3d(:, :, :, 1) < 0.) q3d(:, :, :, 1) = 1E-10
    DO l = 1, llm
       q3d(:, 1, l, 1) = SUM(aire_2d(:, 1) * q3d(:, 1, l, 1)) / apoln
       q3d(:, jjm + 1, l, 1) &
            = SUM(aire_2d(:, jjm + 1) * q3d(:, jjm + 1, l, 1)) / apols
    ENDDO

    q3d(:, :, :, 2:4) = 0. ! liquid water, radon and lead

    if (nqmx >= 5) then
       ! Ozone:
       call regr_lat_time_coefoz
       call regr_pr_o3(q3d(:, :, :, 5))
       ! Convert from mole fraction to mass fraction:
       q3d(:, :, :, 5) = q3d(:, :, :, 5)  * 48. / 29.
    end if

    tsol(:) = pack(tsol_2d, dyn_phy)
    qsol(:) = pack(qsol_2d, dyn_phy)
    sn(:) = 0. ! snow
    radsol(:) = 0.
    tslab(:) = 0. ! IM "slab" ocean
    seaice(:) = 0.
    rugmer(:) = 0.001
    zmea(:) = pack(relief, dyn_phy)
    zstd(:) = pack(zstd_2d, dyn_phy)
    zsig(:) = pack(zsig_2d, dyn_phy)
    zgam(:) = pack(zgam_2d, dyn_phy)
    zthe(:) = pack(zthe_2d, dyn_phy)
    zpic(:) = pack(zpic_2d, dyn_phy)
    zval(:) = pack(zval_2d, dyn_phy)

    rugsrel(:) = 0.
    IF (ok_orodr) rugsrel(:) = MAX(1.e-05, zstd(:) * zsig(:) / 2)

    ! On initialise les sous-surfaces:
    ! Lecture du fichier glace de terre pour fixer la fraction de terre 
    ! et de glace de terre:
    CALL flininfo("landiceref.nc", iml_lic, jml_lic, llm_tmp, &
         ttm_tmp, fid)
    ALLOCATE(lat_lic(iml_lic, jml_lic))
    ALLOCATE(lon_lic(iml_lic, jml_lic))
    ALLOCATE(dlon_lic(iml_lic))
    ALLOCATE(dlat_lic(jml_lic))
    ALLOCATE(fraclic(iml_lic, jml_lic))
    CALL flinopen_nozoom("landiceref.nc", iml_lic, jml_lic, &
         llm_tmp, lon_lic, lat_lic, lev, ttm_tmp, itaul, date, dt,  &
         fid)
    CALL flinget(fid, 'landice', iml_lic, jml_lic, llm_tmp, ttm_tmp &
         , 1, 1, fraclic)
    CALL flinclo(fid)

    ! Interpolation sur la grille T du modèle :
    PRINT *, 'Dimensions de "landice"'
    print *, "iml_lic = ", iml_lic
    print *, "jml_lic = ", jml_lic

    ! Si les coordonnées sont en degrés, on les transforme :
    IF (MAXVAL( lon_lic(:, :) ) > pi)  THEN
       lon_lic(:, :) = lon_lic(:, :) * pi / 180.
    ENDIF
    IF (maxval( lat_lic(:, :) ) > pi) THEN 
       lat_lic(:, :) = lat_lic(:, :) * pi/ 180.
    ENDIF

    dlon_lic(:) = lon_lic(:, 1)
    dlat_lic(:) = lat_lic(1, :) 

    flic_tmp(:iim, :) = grille_m(dlon_lic, dlat_lic, fraclic, rlonv(:iim), &
         rlatu)
    flic_tmp(iim + 1, :) = flic_tmp(1, :)

    ! Passage sur la grille physique
    pctsrf(:, :)=0.
    pctsrf(:, is_lic) = pack(flic_tmp, dyn_phy)
    ! Adéquation avec le maque terre/mer
    WHERE (pctsrf(:, is_lic) < EPSFRA ) pctsrf(:, is_lic) = 0.
    WHERE (zmasq(:) < EPSFRA) pctsrf(:, is_lic) = 0.
    pctsrf(:, is_ter) = zmasq(:)
    where (zmasq(:) > EPSFRA)
       where (pctsrf(:, is_lic) >= zmasq(:))
          pctsrf(:, is_lic) = zmasq(:)
          pctsrf(:, is_ter) = 0.
       elsewhere
          pctsrf(:, is_ter) = zmasq(:) - pctsrf(:, is_lic)
          where (pctsrf(:, is_ter) < EPSFRA)
             pctsrf(:, is_ter) = 0.
             pctsrf(:, is_lic) = zmasq(:)
          end where
       end where
    end where

    ! Sous-surface océan et glace de mer (pour démarrer on met glace
    ! de mer à 0) :
    pctsrf(:, is_oce) = 1. - zmasq(:)
    WHERE (pctsrf(:, is_oce) < EPSFRA) pctsrf(:, is_oce) = 0.

    ! Vérification que somme des sous-surfaces vaut 1:
    ji = count(abs(sum(pctsrf(:, :), dim = 2) - 1. ) > EPSFRA)
    IF (ji /= 0) PRINT *, 'Problème répartition sous maille pour', ji, 'points'

    ! Calcul intermédiaire:
    CALL massdair(p3d, masse)

    print *, 'ALPHAX = ', alphax

    forall  (l = 1:llm)
       masse(:, 1, l) = SUM(aire_2d(:iim, 1) * masse(:iim, 1, l)) / apoln
       masse(:, jjm + 1, l) = &
            SUM(aire_2d(:iim, jjm + 1) * masse(:iim, jjm + 1, l)) / apols
    END forall

    ! Initialisation pour traceurs:
    call iniadvtrac
    ! Ecriture:
    CALL inidissip(lstardis, nitergdiv, nitergrot, niterh, tetagdiv, &
         tetagrot, tetatemp)
    itau_dyn = 0
    itau_phy = 0
    day_ref = dayref
    annee_ref = anneeref

    CALL geopot(ip1jmp1, tpot, pk , pks,  phis  , phi)
    CALL caldyn0(0, uvent, vvent, tpot, psol, masse, pk, phis, phi, w, &
         pbaru, pbarv, 0)
    CALL dynredem0("start.nc", dayref, phis)
    CALL dynredem1("start.nc", 0., vvent, uvent, tpot, q3d, masse, psol)

    ! Ecriture état initial physique:
    print *, 'dtvr = ', dtvr
    print *, "iphysiq = ", iphysiq
    print *, "nbapp_rad = ", nbapp_rad
    phystep   = dtvr * REAL(iphysiq)
    radpas    = NINT (86400./phystep/ nbapp_rad)
    print *, 'phystep = ', phystep
    print *, "radpas = ", radpas

    ! Initialisations :
    tsolsrf(:, is_ter) = tsol
    tsolsrf(:, is_lic) = tsol
    tsolsrf(:, is_oce) = tsol
    tsolsrf(:, is_sic) = tsol
    snsrf(:, is_ter) = sn
    snsrf(:, is_lic) = sn
    snsrf(:, is_oce) = sn
    snsrf(:, is_sic) = sn
    albe(:, is_ter) = 0.08
    albe(:, is_lic) = 0.6
    albe(:, is_oce) = 0.5
    albe(:, is_sic) = 0.6
    alblw = albe
    evap(:, :) = 0.
    qsolsrf(:, is_ter) = 150.
    qsolsrf(:, is_lic) = 150.
    qsolsrf(:, is_oce) = 150.
    qsolsrf(:, is_sic) = 150.
    tsoil = spread(spread(tsol, 2, nsoilmx), 3, nbsrf)
    rain_fall = 0.
    snow_fall = 0.
    solsw = 165.
    sollw = -53.
    t_ancien = 273.15
    q_ancien = 0.
    agesno = 0.
    !IM "slab" ocean
    tslab(:) = tsolsrf(:, is_oce)
    seaice = 0.

    frugs(:, is_oce) = rugmer(:)
    frugs(:, is_ter) = MAX(1.e-05, zstd(:) * zsig(:) / 2)
    frugs(:, is_lic) = MAX(1.e-05, zstd(:) * zsig(:) / 2)
    frugs(:, is_sic) = 0.001
    fder = 0.
    clwcon = 0.
    rnebcon = 0.
    ratqs = 0.
    run_off_lic_0 = 0.

    call phyredem("startphy.nc", radpas, latfi, lonfi, pctsrf, &
         tsolsrf, tsoil, tslab, seaice, qsolsrf, qsol, snsrf, albe, alblw, &
         evap, rain_fall, snow_fall, solsw, sollw, fder, radsol, frugs, &
         agesno, zmea, zstd, zsig, zgam, zthe, zpic, zval, rugsrel, &
         t_ancien, q_ancien, rnebcon, ratqs, clwcon, run_off_lic_0)
    CALL histclo

  END SUBROUTINE etat0

end module etat0_mod
