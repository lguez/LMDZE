module press_coefoz_m

  implicit none

  real, allocatable, save:: plev(:)
  ! (pressure level of Mobidic input data, converted to Pa, in
  ! ascending order)

  real, allocatable, save:: press_in_edg(:)
  ! (edges of pressure intervals for Mobidic input data, in Pa, in
  ! ascending order)

contains

  subroutine press_coefoz

    ! This procedure is called once per run.
    ! It reads the pressure levels from "coefoz_LMDZ.nc".
    ! We assume that, in "coefoz_LMDZ.nc", the pressure levels are in hPa
    ! and strictly increasing.

    use netcdf95, only: nf95_open, nf95_close, nf95_inq_varid, nf95_gw_var
    use netcdf, only: nf90_nowrite

    ! Variables local to the procedure:
    integer ncid, varid ! for NetCDF
    integer n_plev ! number of pressure levels in the input data
    integer k

    !---------------------------------------

    print *, "Call sequence information: press_coefoz"

    call nf95_open("coefoz_LMDZ.nc", nf90_nowrite, ncid)

    call nf95_inq_varid(ncid, "plev", varid)
    call nf95_gw_var(ncid, varid, plev)
    ! Convert from hPa to Pa because "regr_pr_av" and "regr_pr_int"
    ! require so:
    plev = plev * 100.
    n_plev = size(plev)

    call nf95_close(ncid)

    ! Compute edges of pressure intervals:
    allocate(press_in_edg(n_plev + 1))
    press_in_edg(1) = 0.
    ! We choose edges halfway in logarithm:
    forall (k = 2:n_plev) press_in_edg(k) = sqrt(plev(k - 1) * plev(k))
    press_in_edg(n_plev + 1) = huge(0.)
    ! (infinity, but any value guaranteed to be greater than the
    ! surface pressure would do)

  end subroutine press_coefoz

end module press_coefoz_m
