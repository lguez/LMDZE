module etat0_mod

  use indicesol, only: nbsrf
  use dimphy, only: klon

  IMPLICIT NONE

  REAL pctsrf(klon, nbsrf)
  ! ("pctsrf(i, :)" is the composition of the surface at horizontal
  ! position "i")

  private nbsrf, klon

contains

  SUBROUTINE etat0(phis)

    ! From "etat0_netcdf.F", version 1.3, 2005/05/25 13:10:09

    use caldyn0_m, only: caldyn0
    use comconst, only: cpp, kappa, iniconst
    use comgeom, only: aire_2d, apoln, apols, cu_2d, cv_2d, inigeom
    use conf_gcm_m, only: nday, day_step, iphysiq
    use dimens_m, only: iim, jjm, llm, nqmx
    use dimphy, only: zmasq
    use dimsoil, only: nsoilmx
    use disvert_m, only: ap, bp, preff, pa, disvert
    use dynetat0_m, only: day_ref, annee_ref, xprimp025, xprimm025, rlatu1, &
         rlatu2, rlatu, rlatv, yprimu1, yprimu2, rlonu, rlonv, xprimu, xprimv
    use dynredem0_m, only: dynredem0
    use dynredem1_m, only: dynredem1
    use exner_hyb_m, only: exner_hyb
    use fxhyp_m, only: fxhyp
    use fyhyp_m, only: fyhyp
    use geopot_m, only: geopot
    use grid_atob, only: grille_m
    use grid_change, only: init_dyn_phy, dyn_phy
    use indicesol, only: is_oce, is_sic, is_ter, is_lic, epsfra
    use iniadvtrac_m, only: iniadvtrac
    use inifilr_m, only: inifilr
    use massdair_m, only: massdair
    use netcdf, only: nf90_nowrite
    use netcdf95, only: nf95_close, nf95_get_var, nf95_gw_var, nf95_put_var, &
         nf95_inq_varid, nf95_open
    use nr_util, only: pi, assert
    use phyetat0_m, only: rlat, rlon
    use phyredem0_m, only: phyredem0, ncid_restartphy
    use phyredem_m, only: phyredem
    use q_sat_m, only: q_sat
    use regr_lat_time_coefoz_m, only: regr_lat_time_coefoz
    use regr_pr_o3_m, only: regr_pr_o3
    use startdyn, only: start_init_dyn
    USE start_init_orog_m, only: start_init_orog, mask
    use start_init_phys_m, only: start_init_phys
    use start_inter_3d_m, only: start_inter_3d
    use test_disvert_m, only: test_disvert
    use unit_nml_m, only: unit_nml

    REAL, intent(out):: phis(:, :) ! (iim + 1, jjm + 1)
    ! surface geopotential, in m2 s-2

    ! Local:

    REAL, dimension(iim + 1, jjm + 1, llm):: ucov, t3d, teta
    REAL vcov(iim + 1, jjm, llm)

    REAL q(iim + 1, jjm + 1, llm, nqmx)
    ! (mass fractions of trace species
    ! "q(i, j, l)" is at longitude "rlonv(i)", latitude "rlatu(j)"
    ! and pressure level "pls(i, j, l)".)

    real qsat(iim + 1, jjm + 1, llm) ! mass fraction of saturating water vapor
    REAL qsolsrf(klon, nbsrf), snsrf(klon, nbsrf) 
    REAL albe(klon, nbsrf), evap(klon, nbsrf)
    REAL tsoil(klon, nsoilmx, nbsrf) 
    REAL null_array(klon)
    REAL solsw(klon), sollw(klon)
    !IM "slab" ocean
    REAL frugs(klon, nbsrf), agesno(klon, nbsrf)
    REAL rugmer(klon)
    real, dimension(iim + 1, jjm + 1):: zmea_2d, zstd_2d, zsig_2d, zgam_2d
    real, dimension(iim + 1, jjm + 1):: zthe_2d, zpic_2d, zval_2d
    real, dimension(iim + 1, jjm + 1):: tsol_2d, qsol_2d, ps
    REAL zmea(klon), zstd(klon)
    REAL zsig(klon), zgam(klon)
    REAL zthe(klon)
    REAL zpic(klon), zval(klon)
    REAL t_ancien(klon, llm), q_ancien(klon, llm)
    real clwcon(klon, llm), rnebcon(klon, llm), ratqs(klon, llm)

    ! D\'eclarations pour lecture glace de mer :
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
    real sig1(klon, llm) ! section adiabatic updraft
    real w01(klon, llm) ! vertical velocity within adiabatic updraft

    real pls(iim + 1, jjm + 1, llm)
    ! (pressure at mid-layer of LMDZ grid, in Pa)
    ! "pls(i, j, l)" is at longitude "rlonv(i)", latitude "rlatu(j)",
    ! for layer "l")

    REAL p3d(iim + 1, jjm + 1, llm+1) ! pressure at layer interfaces, in Pa
    ! ("p3d(i, j, l)" is at longitude "rlonv(i)", latitude "rlatu(j)",
    ! for interface "l")

    namelist /etat0_nml/ day_ref, annee_ref

    !---------------------------------

    print *, "Call sequence information: etat0"

    print *, "Enter namelist 'etat0_nml'."
    read(unit=*, nml=etat0_nml)
    write(unit_nml, nml=etat0_nml)

    CALL iniconst

    ! Construct a grid:

    pa = 5e4
    CALL disvert
    call test_disvert

    CALL fyhyp(rlatu, rlatv, rlatu2, yprimu2, rlatu1, yprimu1)
    CALL fxhyp(xprimm025, rlonv, xprimv, rlonu, xprimu, xprimp025)

    rlatu(1) = pi / 2.
    rlatu(jjm + 1) = -rlatu(1)

    CALL inigeom
    CALL inifilr

    rlat(1) = 90.
    rlat(2:klon-1) = pack(spread(rlatu(2:jjm), 1, iim), .true.) * 180. / pi
    ! (with conversion to degrees)
    rlat(klon) = - 90.

    rlon(1) = 0.
    rlon(2:klon-1) = pack(spread(rlonv(:iim), 2, jjm - 1), .true.) * 180. / pi
    ! (with conversion to degrees)
    rlon(klon) = 0.

    call start_init_orog(phis, zmea_2d, zstd_2d, zsig_2d, zgam_2d, zthe_2d, &
         zpic_2d, zval_2d) ! also compute "mask"
    call init_dyn_phy ! define the mask "dyn_phy" for distinct grid points
    zmasq = pack(mask, dyn_phy)
    PRINT *, 'Masque construit'

    call start_init_phys(tsol_2d, qsol_2d)
    CALL start_init_dyn(tsol_2d, phis, ps)

    ! Compute pressure on intermediate levels:
    forall(l = 1: llm + 1) p3d(:, :, l) = ap(l) + bp(l) * ps
    CALL exner_hyb(ps, p3d, pks, pk)
    call assert(MINVAL(pk) /= MAXVAL(pk), '"pk" should not be constant')

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
    PRINT *, 'minval(t3d) = ', minval(t3d)
    print *, "maxval(t3d) = ", maxval(t3d)

    teta(:iim, :, :) = t3d(:iim, :, :) * cpp / pk(:iim, :, :)
    teta(iim + 1, :, :) = teta(1, :, :)
    DO l = 1, llm
       teta(:, 1, l) = SUM(aire_2d(:, 1) * teta(:, 1, l)) / apoln
       teta(:, jjm + 1, l) = SUM(aire_2d(:, jjm + 1) * teta(:, jjm + 1, l)) &
            / apols
    ENDDO

    ! Calcul de l'humidit\'e \`a saturation :
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
       call regr_pr_o3(p3d, q(:, :, :, 5))
       ! Convert from mole fraction to mass fraction:
       q(:, :, :, 5) = q(:, :, :, 5) * 48. / 29.
    end if

    null_array = 0.
    rugmer = 0.001
    zmea = pack(zmea_2d, dyn_phy)
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

    ! Interpolation sur la grille T du mod\`ele :
    PRINT *, 'Dimensions de "landiceref.nc"'
    print *, "iml_lic = ", iml_lic
    print *, "jml_lic = ", jml_lic

    ! Si les coordonn\'ees sont en degr\'es, on les transforme :
    IF (MAXVAL(dlon_lic) > pi) THEN
       dlon_lic = dlon_lic * pi / 180.
    ENDIF
    IF (maxval(dlat_lic) > pi) THEN 
       dlat_lic = dlat_lic * pi/ 180.
    ENDIF

    flic_tmp(:iim, :) = grille_m(dlon_lic, dlat_lic, fraclic, rlonv(:iim), &
         rlatu)
    flic_tmp(iim + 1, :) = flic_tmp(1, :)

    deallocate(dlon_lic, dlat_lic) ! pointers

    ! Passage sur la grille physique
    pctsrf = 0.
    pctsrf(:, is_lic) = pack(flic_tmp, dyn_phy)
    ! Ad\'equation avec le maque terre/mer
    WHERE (pctsrf(:, is_lic) < EPSFRA) pctsrf(:, is_lic) = 0.
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

    ! Sous-surface oc\'ean et glace de mer (pour d\'emarrer on met glace
    ! de mer \`a 0) :
    pctsrf(:, is_oce) = 1. - zmasq
    WHERE (pctsrf(:, is_oce) < EPSFRA) pctsrf(:, is_oce) = 0.

    ! V\'erification que la somme des sous-surfaces vaut 1 :
    ji = count(abs(sum(pctsrf, dim = 2) - 1.) > EPSFRA)
    IF (ji /= 0) then
       PRINT *, 'Bad surface percentages for ', ji, 'points'
    end IF

    ! Calcul interm\'ediaire :
    CALL massdair(p3d, masse)

    forall (l = 1:llm)
       masse(:, 1, l) = SUM(aire_2d(:iim, 1) * masse(:iim, 1, l)) / apoln
       masse(:, jjm + 1, l) = &
            SUM(aire_2d(:iim, jjm + 1) * masse(:iim, jjm + 1, l)) / apols
    END forall

    call iniadvtrac
    CALL geopot(teta, pk , pks, phis, phi)
    CALL caldyn0(ucov, vcov, teta, ps, pk, phis, phi)
    CALL dynredem0(day_ref, phis)
    CALL dynredem1(vcov, ucov, teta, q, masse, ps, itau = 0)

    ! Initialisations :
    snsrf = 0.
    albe(:, is_ter) = 0.08
    albe(:, is_lic) = 0.6
    albe(:, is_oce) = 0.5
    albe(:, is_sic) = 0.6
    evap = 0.
    qsolsrf = 150.
    tsoil = spread(spread(pack(tsol_2d, dyn_phy), 2, nsoilmx), 3, nbsrf)
    solsw = 165.
    sollw = -53.
    t_ancien = 273.15
    q_ancien = 0.
    agesno = 0.

    frugs(:, is_oce) = rugmer
    frugs(:, is_ter) = MAX(1e-5, zstd * zsig / 2)
    frugs(:, is_lic) = MAX(1e-5, zstd * zsig / 2)
    frugs(:, is_sic) = 0.001
    clwcon = 0.
    rnebcon = 0.
    ratqs = 0.
    sig1 = 0.
    w01 = 0.

    nday = 0
    call phyredem0(lmt_pas = day_step / iphysiq, itau_phy = 0)

    call nf95_inq_varid(ncid_restartphy, "trs", varid)
    call nf95_put_var(ncid_restartphy, varid, null_array)

    call phyredem(pctsrf, tsoil(:, 1, :), tsoil, qsolsrf, &
         pack(qsol_2d, dyn_phy), snsrf, albe, evap, null_array, null_array, &
         solsw, sollw, null_array, null_array, frugs, agesno, zmea, zstd, &
         zsig, zgam, zthe, zpic, zval, t_ancien, q_ancien, rnebcon, ratqs, &
         clwcon, null_array, sig1, w01)

  END SUBROUTINE etat0

end module etat0_mod
