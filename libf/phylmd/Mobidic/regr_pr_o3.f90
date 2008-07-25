module regr_pr_o3_m

  implicit none

contains

  subroutine regr_pr_o3(o3_mob_regr)

    ! "regr_pr_o3" stands for "regrid pressure ozone".
    ! This procedure reads Mobidic ozone mole fraction from
    ! "coefoz_LMDZ.nc" at the initial day and regrids it in pressure.
    ! Regridding is by averaging, assuming a step function.
    ! We assume that, in the input file, the pressure levels are in hPa
    ! and strictly increasing.

    use conf_gcm_m, only: dayref
    use dimens_m, only: iim, jjm, llm
    use netcdf95, only: nf95_open, nf95_close, nf95_inq_varid, handle_err, &
         nf95_get_coord
    use netcdf, only:  nf90_nowrite, nf90_get_var
    use regr_pr, only: regr_pr_av
    use numer_rec, only: assert

    real, intent(out):: o3_mob_regr(:, :, :) ! (iim + 1, jjm + 1, llm)
    ! (ozone mole fraction from Mobidic adapted to the LMDZ grid)
    ! ("o3_mob_regr(i, j, l)" is at longitude "rlonv(i)", latitude
    ! "rlatu(j)" and pressure level "pls(i, j, l)")

    ! Variables local to the procedure:

    real, pointer:: plev(:)
    ! (pressure levels of Mobidic data, in Pa, in strictly increasing order)

    real, allocatable:: press_in_edg(:)
    ! (edges of pressure intervals for Mobidic data, in Pa, in strictly
    ! increasing order)

    integer ncid, varid, ncerr ! for NetCDF
    integer n_plev ! number of pressure levels in Mobidic data
    integer j

    real, allocatable:: r_mob(:, :)! (jjm + 1, n_plev)
    ! (ozone mole fraction from Mobidic at day "dayref")
    ! (r_mob(j, k) is at latitude "rlatu(j)" and pressure level "plev(k)".)

    !------------------------------------------------------------

    print *, "Call sequence information: regr_pr_o3"
    call assert(shape(o3_mob_regr) == (/iim + 1, jjm + 1, llm/), "regr_pr_o3")

    call nf95_open("coefoz_LMDZ.nc", nf90_nowrite, ncid)

    call nf95_get_coord(ncid, "plev", plev)
    ! Convert from hPa to Pa because "regr_pr_av" requires so:
    plev = plev * 100.
    n_plev = size(plev)

    ! Compute edges of pressure intervals:
    allocate(press_in_edg(n_plev + 1))
    press_in_edg(1) = 0.
    ! We choose edges halfway in logarithm:
    forall (j = 2:n_plev) press_in_edg(j) = sqrt(plev(j - 1) * plev(j))
    press_in_edg(n_plev + 1) = huge(0.)
    ! (infinity, but any value guaranteed to be greater than the
    ! surface pressure would do)

    deallocate(plev) ! pointer

    call nf95_inq_varid(ncid, "r_Mob", varid)
    allocate(r_mob(jjm + 1, n_plev))

    ! Get data at the right day from the input file:
    ncerr = nf90_get_var(ncid, varid, r_mob, start=(/1, 1, dayref/))
    call handle_err("nf90_get_var r_Mob", ncerr)
    ! Latitudes are in increasing order in the input file while
    ! "rlatu" is in decreasing order so we need to invert order:
    r_mob = r_mob(jjm+1:1:-1, :)

    call nf95_close(ncid)

    o3_mob_regr = regr_pr_av(r_mob, press_in_edg)

  end subroutine regr_pr_o3

end module regr_pr_o3_m
