module etat0_mod

  use indicesol, only: nbsrf
  use dimphy, only: klon

  IMPLICIT NONE

  REAL pctsrf(klon, nbsrf)
  ! ("pctsrf(i, :)" is the composition of the surface at horizontal
  !  position "i")

  private nbsrf, klon

contains

  SUBROUTINE etat0

    ! From "etat0_netcdf.F", version 1.3 2005/05/25 13:10:09

    use caldyn0_m, only: caldyn0
    use comconst, only: dtvr, daysec, cpp, kappa
    use comgeom, only: rlatu, rlonv, rlonu, rlatv, aire_2d, apoln, apols, &
         cu_2d, cv_2d
    use conf_gcm_m, only: day_step, iphysiq, dayref, anneeref
    use dimens_m, only: iim, jjm, llm, nqmx
    use dimphy, only: zmasq
    use dimsoil, only: nsoilmx
    use disvert_m, only: ap, bp, preff, pa
    use dynredem0_m, only: dynredem0
    use dynredem1_m, only: dynredem1
    use exner_hyb_m, only: exner_hyb
    use geopot_m, only: geopot
    use grid_atob, only: grille_m
    use grid_change, only: init_dyn_phy, dyn_phy
    use histclo_m, only: histclo
    use indicesol, only: is_oce, is_sic, is_ter, is_lic, epsfra
    use iniadvtrac_m, only: iniadvtrac
    use inifilr_m, only: inifilr
    use inigeom_m, only: inigeom
    use massdair_m, only: massdair
    use netcdf, only: nf90_nowrite
    use netcdf95, only: nf95_close, nf95_get_var, nf95_gw_var, &
         nf95_inq_varid, nf95_open
    use nr_util, only: pi
    use paramet_m, only: ip1jm, ip1jmp1
    use phyredem_m, only: phyredem
    use pressure_var, only: pls, p3d
    use q_sat_m, only: q_sat
    use regr_lat_time_coefoz_m, only: regr_lat_time_coefoz
    use regr_pr_o3_m, only: regr_pr_o3
    use serre, only: alphax
    use startdyn, only: start_init_dyn
    USE start_init_orog_m, only: start_init_orog, mask, phis
    use start_init_phys_m, only: start_init_phys
    use start_inter_3d_m, only: start_inter_3d
    use temps, only: itau_phy, annee_ref, day_ref

    ! Variables local to the procedure:

    REAL latfi(klon), lonfi(klon)
    ! (latitude and longitude of a point of the scalar grid identified
    ! by a simple index, in °)

    REAL, dimension(iim + 1, jjm + 1, llm):: ucov, t3d, tpot
    REAL vcov(iim + 1, jjm, llm)

    REAL q(iim + 1, jjm + 1, llm, nqmx)
    ! (mass fractions of trace species
    ! "q(i, j, l)" is at longitude "rlonv(i)", latitude "rlatu(j)"
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
    real, dimension(iim + 1, jjm + 1):: tsol_2d, qsol_2d, psol
    REAL zmea(klon), zstd(klon)
    REAL zsig(klon), zgam(klon)
    REAL zthe(klon)
    REAL zpic(klon), zval(klon)
    REAL t_ancien(klon, llm), q_ancien(klon, llm)      !
    REAL run_off_lic_0(klon)
    real clwcon(klon, llm), rnebcon(klon, llm), ratqs(klon, llm)

    ! Déclarations pour lecture glace de mer :
    INTEGER iml_lic, jml_lic
    INTEGER ncid, varid
    REAL, pointer:: dlon_lic(:), dlat_lic(:)
    REAL, ALLOCATABLE:: fraclic(:, :) ! fraction land ice
    REAL flic_tmp(iim + 1, jjm + 1) ! fraction land ice temporary

    INTEGER l, ji

    REAL pk(iim + 1, jjm + 1, llm) ! fonction d'Exner aux milieux des couches 
    real pks(iim + 1, jjm + 1)

    REAL masse(iim + 1, jjm + 1, llm)
    REAL phi(iim + 1, jjm + 1, llm)
    REAL pbaru(ip1jmp1, llm), pbarv(ip1jm, llm)
    REAL w(ip1jmp1, llm)
    REAL phystep

    real sig1(klon, llm) ! section adiabatic updraft
    real w01(klon, llm) ! vertical velocity within adiabatic updraft

    !---------------------------------

    print *, "Call sequence information: etat0"

    dtvr = daysec / real(day_step)
    print *, 'dtvr = ', dtvr

    ! Construct a grid:

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
         zval_2d) ! also compute "mask" and "phis"
    call init_dyn_phy ! define the mask "dyn_phy" for distinct grid points
    zmasq = pack(mask, dyn_phy)
    PRINT *, 'Masque construit'

    call start_init_phys(tsol_2d, qsol_2d)
    CALL start_init_dyn(tsol_2d, psol)

    ! Compute pressure on intermediate levels:
    forall(l = 1: llm + 1) p3d(:, :, l) = ap(l) + bp(l) * psol
    CALL exner_hyb(psol, p3d, pks, pk)
    IF (MINVAL(pk) == MAXVAL(pk)) then
       print *, '"pk" should not be constant'
       stop 1
    end IF

    pls = preff * (pk / cpp)**(1. / kappa)
    PRINT *, "minval(pls) = ", minval(pls)
    print *, "maxval(pls) = ", maxval(pls)

    call start_inter_3d('U', rlonv, rlatv, pls, ucov)
    forall (l = 1: llm) ucov(:iim, :, l) = ucov(:iim, :, l) * cu_2d(:iim, :)
    ucov(iim+1, :, :) = ucov(1, :, :)

    call start_inter_3d('V', rlonu, rlatu(:jjm), pls(:, :jjm, :), vcov)
    forall (l = 1: llm) vcov(:iim, :, l) = vcov(:iim, :, l) * cv_2d(:iim, :)
    vcov(iim + 1, :, :) = vcov(1, :, :)

    call start_inter_3d('TEMP', rlonu, rlatv, pls, t3d)
    PRINT *,  'minval(t3d) = ', minval(t3d)
    print *, "maxval(t3d) = ", maxval(t3d)

    tpot(:iim, :, :) = t3d(:iim, :, :) * cpp / pk(:iim, :, :)
    tpot(iim + 1, :, :) = tpot(1, :, :)
    DO l=1, llm
       tpot(:, 1, l) = SUM(aire_2d(:, 1) * tpot(:, 1, l)) / apoln
       tpot(:, jjm + 1, l) = SUM(aire_2d(:, jjm + 1) * tpot(:, jjm + 1, l)) &
            / apols
    ENDDO

    ! Calcul de l'humidité à saturation :
    qsat = q_sat(t3d, pls)
    PRINT *, "minval(qsat) = ", minval(qsat)
    print *, "maxval(qsat) = ", maxval(qsat)
    IF (MINVAL(qsat) == MAXVAL(qsat)) stop '"qsat" should not be constant'

    ! Water vapor:
    call start_inter_3d('R', rlonu, rlatv, pls, q(:, :, :, 1))
    q(:, :, :, 1) = 0.01 * q(:, :, :, 1) * qsat
    WHERE (q(:, :, :, 1) < 0.) q(:, :, :, 1) = 1E-10
    DO l = 1, llm
       q(:, 1, l, 1) = SUM(aire_2d(:, 1) * q(:, 1, l, 1)) / apoln
       q(:, jjm + 1, l, 1) &
            = SUM(aire_2d(:, jjm + 1) * q(:, jjm + 1, l, 1)) / apols
    ENDDO

    q(:, :, :, 2:4) = 0. ! liquid water, radon and lead

    if (nqmx >= 5) then
       ! Ozone:
       call regr_lat_time_coefoz
       call regr_pr_o3(q(:, :, :, 5))
       ! Convert from mole fraction to mass fraction:
       q(:, :, :, 5) = q(:, :, :, 5)  * 48. / 29.
    end if

    tsol = pack(tsol_2d, dyn_phy)
    qsol = pack(qsol_2d, dyn_phy)
    sn = 0. ! snow
    radsol = 0.
    tslab = 0. ! IM "slab" ocean
    seaice = 0.
    rugmer = 0.001
    zmea = pack(relief, dyn_phy)
    zstd = pack(zstd_2d, dyn_phy)
    zsig = pack(zsig_2d, dyn_phy)
    zgam = pack(zgam_2d, dyn_phy)
    zthe = pack(zthe_2d, dyn_phy)
    zpic = pack(zpic_2d, dyn_phy)
    zval = pack(zval_2d, dyn_phy)

    ! On initialise les sous-surfaces.
    ! Lecture du fichier glace de terre pour fixer la fraction de terre 
    ! et de glace de terre :

    call nf95_open("landiceref.nc", nf90_nowrite, ncid)

    call nf95_inq_varid(ncid, 'longitude', varid)
    call nf95_gw_var(ncid, varid, dlon_lic)
    iml_lic = size(dlon_lic)

    call nf95_inq_varid(ncid, 'latitude', varid)
    call nf95_gw_var(ncid, varid, dlat_lic)
    jml_lic = size(dlat_lic)

    call nf95_inq_varid(ncid, 'landice', varid)
    ALLOCATE(fraclic(iml_lic, jml_lic))
    call nf95_get_var(ncid, varid, fraclic)

    call nf95_close(ncid)

    ! Interpolation sur la grille T du modèle :
    PRINT *, 'Dimensions de "landiceref.nc"'
    print *, "iml_lic = ", iml_lic
    print *, "jml_lic = ", jml_lic

    ! Si les coordonnées sont en degrés, on les transforme :
    IF (MAXVAL( dlon_lic ) > pi)  THEN
       dlon_lic = dlon_lic * pi / 180.
    ENDIF
    IF (maxval( dlat_lic ) > pi) THEN 
       dlat_lic = dlat_lic * pi/ 180.
    ENDIF

    flic_tmp(:iim, :) = grille_m(dlon_lic, dlat_lic, fraclic, rlonv(:iim), &
         rlatu)
    flic_tmp(iim + 1, :) = flic_tmp(1, :)

    deallocate(dlon_lic, dlat_lic) ! pointers

    ! Passage sur la grille physique
    pctsrf = 0.
    pctsrf(:, is_lic) = pack(flic_tmp, dyn_phy)
    ! Adéquation avec le maque terre/mer
    WHERE (pctsrf(:, is_lic) < EPSFRA ) pctsrf(:, is_lic) = 0.
    WHERE (zmasq < EPSFRA) pctsrf(:, is_lic) = 0.
    pctsrf(:, is_ter) = zmasq
    where (zmasq > EPSFRA)
       where (pctsrf(:, is_lic) >= zmasq)
          pctsrf(:, is_lic) = zmasq
          pctsrf(:, is_ter) = 0.
       elsewhere
          pctsrf(:, is_ter) = zmasq - pctsrf(:, is_lic)
          where (pctsrf(:, is_ter) < EPSFRA)
             pctsrf(:, is_ter) = 0.
             pctsrf(:, is_lic) = zmasq
          end where
       end where
    end where

    ! Sous-surface océan et glace de mer (pour démarrer on met glace
    ! de mer à 0) :
    pctsrf(:, is_oce) = 1. - zmasq
    WHERE (pctsrf(:, is_oce) < EPSFRA) pctsrf(:, is_oce) = 0.

    ! Vérification que somme des sous-surfaces vaut 1:
    ji = count(abs(sum(pctsrf, dim = 2) - 1.) > EPSFRA)
    IF (ji /= 0) then
       PRINT *, 'Problème répartition sous maille pour ', ji, 'points'
    end IF

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
    itau_phy = 0
    day_ref = dayref
    annee_ref = anneeref

    CALL geopot(tpot, pk , pks,  phis, phi)
    CALL caldyn0(ucov, vcov, tpot, psol, masse, pk, phis, phi, w, pbaru, &
         pbarv)
    CALL dynredem0("start.nc", dayref, phis)
    CALL dynredem1("start.nc", vcov, ucov, tpot, q, masse, psol, itau=0)

    ! Ecriture état initial physique:
    print *, "iphysiq = ", iphysiq
    phystep   = dtvr * REAL(iphysiq)
    print *, 'phystep = ', phystep

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
    evap = 0.
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
    tslab = tsolsrf(:, is_oce)
    seaice = 0.

    frugs(:, is_oce) = rugmer
    frugs(:, is_ter) = MAX(1.e-05, zstd * zsig / 2)
    frugs(:, is_lic) = MAX(1.e-05, zstd * zsig / 2)
    frugs(:, is_sic) = 0.001
    fder = 0.
    clwcon = 0.
    rnebcon = 0.
    ratqs = 0.
    run_off_lic_0 = 0.
    sig1 = 0.
    w01 = 0.

    call phyredem("startphy.nc", latfi, lonfi, pctsrf, &
         tsolsrf, tsoil, tslab, seaice, qsolsrf, qsol, snsrf, albe, alblw, &
         evap, rain_fall, snow_fall, solsw, sollw, fder, radsol, frugs, &
         agesno, zmea, zstd, zsig, zgam, zthe, zpic, zval, &
         t_ancien, q_ancien, rnebcon, ratqs, clwcon, run_off_lic_0, sig1, w01)
    CALL histclo

  END SUBROUTINE etat0

end module etat0_mod
