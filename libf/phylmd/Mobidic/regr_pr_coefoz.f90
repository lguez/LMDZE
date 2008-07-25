module regr_pr_coefoz

  implicit none

contains

  subroutine regr_pr_av_coefoz(ncid, name, julien, press_in_edg, v3)

    ! "regr_pr_av_coefoz" stands for "regrid pressure averaging
    ! coefficient ozone".
    ! This procedure reads a single Mobidic ozone coefficient from
    ! "coefoz_LMDZ.nc", at the current day, regrids this parameter in
    ! pressure to the LMDZ vertical grid and packs it to the LMDZ
    ! horizontal "physics" grid.
    ! Regridding in pressure is done by averaging a step function.

    use dimens_m, only: iim, jjm, llm
    use dimphy, only: klon
    use netcdf95, only: nf95_inq_varid, handle_err
    use netcdf, only: nf90_get_var
    use grid_change, only: dyn_phy
    use regr_pr, only: regr_pr_av
    use numer_rec, only: assert

    integer, intent(in):: ncid ! NetCDF ID of the file
    character(len=*), intent(in):: name ! of the NetCDF variable
    integer, intent(in):: julien ! jour julien, 1 <= julien <= 360

    real, intent(in):: press_in_edg(:)
    ! (edges of pressure intervals for Mobidic data, in Pa, in
    ! strictly increasing order)

    real, intent(out):: v3(:, :) ! (klon, llm)
    ! (ozone coefficient from Mobidic on the "physics" grid)
    ! ("v3(i, k)" is at longitude "xlon(i)", latitude
    ! "xlat(i)", middle of layer "k".)

    ! Variables local to the procedure:
    integer varid, ncerr
    integer k

    real  v1(jjm + 1, size(press_in_edg) - 1)
    ! (ozone coefficient from "coefoz_LMDZ.nc" at day "julien")
    ! ("v1(j, k)" is at latitude "rlatu(j)" and for
    ! pressure interval "[press_in_edg(k), press_in_edg(k+1)]".)

    real v2(iim + 1, jjm + 1, llm)
    ! (ozone parameter from Mobidic on the "dynamics" grid)
    ! "v2(i, j, k)" is at longitude "rlonv(i)", latitude
    ! "rlatu(j)", middle of layer "k".)

    !--------------------------------------------

    call assert(shape(v3) == (/klon, llm/), "regr_pr_av_coefoz")

    call nf95_inq_varid(ncid, name, varid)

    ! Get data at the right day from the input file:
    ncerr = nf90_get_var(ncid, varid, v1, start=(/1, 1, julien/))
    call handle_err("regr_pr_av_coefoz nf90_get_var " // name, ncerr, ncid)
    ! Latitudes are in increasing order in the input file while
    ! "rlatu" is in decreasing order so we need to invert order:
    v1 = v1(jjm+1:1:-1, :)

    ! Regrid in pressure at each horizontal position:
    v2 = regr_pr_av(v1, press_in_edg)

    forall (k = 1:llm) v3(:, k) = pack(v2(:, :, k), dyn_phy)

  end subroutine regr_pr_av_coefoz

  !***************************************************************

  subroutine regr_pr_int_coefoz(ncid, name, julien, plev, top_value, v3)

    ! This procedure reads a single Mobidic ozone coefficient from
    !"coefoz_LMDZ.nc", at the current day, regrids this parameter in
    ! pressure to the LMDZ vertical grid and packs it to the LMDZ
    ! horizontal "physics" grid.
    ! Regridding is by linear interpolation.

    use dimens_m, only: iim, jjm, llm
    use dimphy, only: klon
    use netcdf95, only: nf95_inq_varid, handle_err
    use netcdf, only: nf90_get_var
    use grid_change, only: dyn_phy
    use regr_pr, only: regr_pr_int
    use numer_rec, only: assert

    integer, intent(in):: ncid ! NetCDF ID of the file
    character(len=*), intent(in):: name ! of the NetCDF variable
    integer, intent(in):: julien ! jour julien, 1 <= julien <= 360

    real, intent(in):: plev(:)
    ! (pressure levels of Mobidic data, in Pa, in strictly increasing order)

    real, intent(in):: top_value
    ! (extra value of ozone coefficient at 0 pressure)

    real, intent(out):: v3(:, :) ! (klon, llm)
    ! (ozone parameter from Mobidic on the "physics" grid)
    ! ("v3(i, k)" is at longitude "xlon(i)", latitude
    ! "xlat(i)", middle of layer "k".)

    ! Variables local to the procedure:
    integer varid, ncerr
    integer k

    real  v1(jjm + 1, 0:size(plev))
    ! (ozone coefficient from "coefoz_LMDZ.nc" at day "julien")
    ! ("v1(j, k >=1)" is at latitude "rlatu(j)" and pressure "plev(k)".)

    real v2(iim + 1, jjm + 1, llm)
    ! (ozone parameter from Mobidic on the "dynamics" grid)
    ! "v2(i, j, k)" is at longitude "rlonv(i)", latitude
    ! "rlatu(j)", middle of layer "k".)

    !--------------------------------------------

    call assert(shape(v3) == (/klon, llm/), "regr_pr_int_coefoz")

    call nf95_inq_varid(ncid, name, varid)

    ! Get data at the right day from the input file:
    ncerr = nf90_get_var(ncid, varid, v1(:, 1:), start=(/1, 1, julien/))
    call handle_err("regr_pr_int_coefoz nf90_get_var " // name, ncerr, ncid)
    ! Latitudes are in increasing order in the input file while
    ! "rlatu" is in decreasing order so we need to invert order:
    v1(:, 1:) = v1(jjm+1:1:-1, 1:)

    ! Complete "v1" with the value at 0 pressure:
    v1(:, 0) = top_value

    ! Regrid in pressure at each horizontal position:
    v2 = regr_pr_int(v1, (/0., plev/))

    forall (k = 1:llm) v3(:, k) = pack(v2(:, :, k), dyn_phy)

  end subroutine regr_pr_int_coefoz

end module regr_pr_coefoz
