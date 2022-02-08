module regr_pr_o3_m

  implicit none

contains

  subroutine regr_pr_o3(p3d, o3_mob_regr)

    ! "regr_pr_o3" stands for "regrid pressure ozone".
    ! This procedure reads Mobidic ozone mole fraction from
    ! "coefoz_LMDZ.nc" at the initial day and regrids it in pressure.
    ! Ozone mole fraction from "coefoz_LMDZ.nc" is a 2D latitude --
    ! pressure variable.
    ! The target horizontal LMDZ grid is the "scalar" grid: "rlonv", "rlatu".
    ! The target vertical LMDZ grid is the grid of layer boundaries.
    ! We assume that the input variable is already on the LMDZ "rlatu"
    ! latitude grid.
    ! The input variable does not depend on longitude, but the
    ! pressure at LMDZ layers does.
    ! Therefore, the values on the LMDZ grid do depend on longitude.
    ! Regridding is by averaging, assuming a step function.
    ! We assume that, in the input file, the pressure levels are in
    ! hPa and strictly increasing.

    use jumble, only: assert
    use netcdf, only:  nf90_nowrite
    use netcdf95, only: nf95_open, nf95_close, nf95_inq_varid, nf95_get_var, &
         nf95_gw_var
    use numer_rec_95, only: regr1_step_av

    use dimensions, only: iim, jjm, llm
    use dynetat0_chosen_m, only: day_ref
    use grid_change, only: dyn_phy

    REAL, intent(in):: p3d(:, :, :) ! (iim + 1, jjm + 1, llm+1) 
    ! pressure at layer interfaces, in Pa
    ! ("p3d(i, j, l)" is at longitude "rlonv(i)", latitude "rlatu(j)",
    ! for interface "l")

    real, intent(out):: o3_mob_regr(:, :, :) ! (iim + 1, jjm + 1, llm)
    ! (ozone mole fraction from Mobidic adapted to the LMDZ grid)
    ! ("o3_mob_regr(i, j, l)" is at longitude "rlonv(i)", latitude
    ! "rlatu(j)" and pressure level "pls(i, j, l)")

    ! Variables local to the procedure:

    real, allocatable:: plev(:)
    ! (pressure levels of Mobidic data, in Pa, in strictly increasing order)

    real, allocatable:: press_in_edg(:)
    ! (edges of pressure intervals for Mobidic data, in Pa, in strictly
    ! increasing order)

    integer ncid, varid ! for NetCDF
    integer n_plev ! number of pressure levels in Mobidic data
    integer i, j

    real, allocatable:: r_mob(:, :)! (jjm + 1, n_plev)
    ! (ozone mole fraction from Mobidic at day "day_ref")
    ! (r_mob(j, k) is at latitude "rlatu(j)" and pressure level "plev(k)".)

    !------------------------------------------------------------

    print *, "Call sequence information: regr_pr_o3"

    call assert(shape(o3_mob_regr) == (/iim + 1, jjm + 1, llm/), &
         "regr_pr_o3 o3_mob_regr")
    call assert(shape(p3d) == (/iim + 1, jjm + 1, llm + 1/), &
         "regr_pr_o3 p3d")

    call nf95_open("coefoz_LMDZ.nc", nf90_nowrite, ncid)

    call nf95_inq_varid(ncid, "plev", varid)
    call nf95_gw_var(ncid, varid, plev)
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

    call nf95_inq_varid(ncid, "r_Mob", varid)
    allocate(r_mob(jjm + 1, n_plev))

    ! Get data at the right day from the input file:
    call nf95_get_var(ncid, varid, r_mob, start=(/1, 1, day_ref/))
    ! Latitudes are in increasing order in the input file while
    ! "rlatu" is in decreasing order so we need to invert order:
    r_mob = r_mob(jjm+1:1:-1, :)

    call nf95_close(ncid)

    ! Regrid in pressure by averaging a step function of pressure:
    do j = 1, jjm + 1
       do i = 1, iim
          if (dyn_phy(i, j)) then
             o3_mob_regr(i, j, llm:1:-1) &
                  = regr1_step_av(r_mob(j, :), press_in_edg, &
                  p3d(i, j, llm+1:1:-1))
             ! (invert order of indices because "p3d" is decreasing)
          end if
       end do
    end do

    ! Duplicate pole values on all longitudes:
    o3_mob_regr(2:, 1, :) = spread(o3_mob_regr(1, 1, :), dim=1, ncopies=iim)
    o3_mob_regr(2:, jjm + 1, :) &
         = spread(o3_mob_regr(1, jjm + 1, :), dim=1, ncopies=iim)

    ! Duplicate first longitude to last longitude:
    o3_mob_regr(iim + 1, 2:jjm, :) = o3_mob_regr(1, 2:jjm, :)

  end subroutine regr_pr_o3

end module regr_pr_o3_m
