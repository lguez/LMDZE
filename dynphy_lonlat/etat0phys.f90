module etat0phys_m

  IMPLICIT NONE

contains

  SUBROUTINE etat0phys(tsol_2d, phis, pctsrf)

    ! From "etat0_netcdf.F", version 1.3, 2005/05/25 13:10:09

    use netcdf95, only: nf95_put_var, nf95_inq_varid

    use dimensions, only: iim, jjm, llm
    use dimphy, only: klon
    use dimsoil, only: nsoilmx
    use dynetat0_chosen_m, only: day_ref
    use grid_change, only: init_dyn_phy, dyn_phy
    use indicesol, only: is_oce, is_sic, is_ter, is_lic, nbsrf
    use phyetat0_m, only: set_lat, set_lon, set_masque
    use phyredem0_m, only: phyredem0, ncid_restartphy
    use phyredem_m, only: phyredem
    USE start_init_orog_m, only: start_init_orog
    use start_init_phys_m, only: start_init_phys
    use start_init_subsurf_m, only: start_init_subsurf

    real, intent(out):: tsol_2d(:, :) ! (iim + 1, jjm + 1)
    ! both soil temperature and surface temperature, in K

    REAL, intent(out):: phis(:, :) ! (iim + 1, jjm + 1)
    ! surface geopotential, in m2 s-2

    REAL, intent(out):: pctsrf(:, :) ! (klon, nbsrf)
    ! "pctsrf(i, :)" is the composition of the surface at horizontal
    ! position "i".

    ! Local:

    REAL fqsurf(klon, nbsrf), fsnow(klon, nbsrf) 
    REAL falbe(klon, nbsrf)
    REAL ftsoil(klon, nsoilmx, nbsrf) 
    REAL null_array(klon)
    REAL solsw(klon), sollw(klon)
    REAL frugs(klon, nbsrf), agesno(klon, nbsrf)
    REAL rugmer(klon)

    real qsol_2d(iim + 1, jjm + 1) ! column-density of water in soil, in kg m-2
    REAL zmea(klon), zstd(klon)
    REAL zsig(klon), zgam(klon)
    REAL zthe(klon)
    REAL zpic(klon), zval(klon)
    REAL t_ancien(klon, llm), q_ancien(klon, llm)
    real clwcon(klon, llm), rnebcon(klon, llm), ratqs(klon, llm)
    INTEGER varid
    real sig1(klon, llm) ! section adiabatic updraft
    real w01(klon, llm) ! vertical velocity within adiabatic updraft

    !---------------------------------

    print *, "Call sequence information: etat0phys"
    call init_dyn_phy ! define the mask "dyn_phy" for distinct grid points

    call start_init_orog(phis, zmea, zstd, zsig, zgam, zthe, zpic, zval)
    ! (also compute "mask")

    call set_lat
    call set_lon
    call set_masque
    call start_init_phys(tsol_2d, qsol_2d)
    null_array = 0.
    rugmer = 0.001
    call start_init_subsurf(pctsrf)

    ! Initialisations :
    fsnow = 0.
    falbe(:, is_ter) = 0.08
    falbe(:, is_lic) = 0.6
    falbe(:, is_oce) = 0.5
    falbe(:, is_sic) = 0.6
    fqsurf = 0.
    ftsoil = spread(spread(pack(tsol_2d, dyn_phy), 2, nsoilmx), 3, nbsrf)
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
    call phyredem0(0)
    call nf95_inq_varid(ncid_restartphy, "trs", varid)
    call nf95_put_var(ncid_restartphy, varid, null_array)
    call phyredem(pctsrf, ftsoil(:, 1, :), ftsoil, fqsurf, &
         pack(qsol_2d, dyn_phy), fsnow, falbe, null_array, null_array, solsw, &
         sollw, null_array, null_array, frugs, agesno, zmea, zstd, zsig, zgam, &
         zthe, zpic, zval, t_ancien, q_ancien, rnebcon, ratqs, clwcon, &
         null_array, sig1, w01)

  END SUBROUTINE etat0phys

end module etat0phys_m
