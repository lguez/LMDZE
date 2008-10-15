module regr_pr_coefoz

  ! Both procedures of this module read a single Mobidic ozone
  ! coefficient from "coefoz_LMDZ.nc", at the current day, regrid this
  ! coefficient in pressure to the LMDZ vertical grid and pack it to the LMDZ
  ! horizontal "physics" grid.
  ! The input data is a 2D latitude -- pressure field.
  ! The target horizontal LMDZ grid is the "scalar" grid: "rlonv", "rlatu".
  ! We assume that the input data is already on the LMDZ "rlatu"
  ! latitude grid.

  implicit none

contains

  subroutine regr_pr_av_coefoz(ncid, name, julien, v3)

    ! "regr_pr_av_coefoz" stands for "regrid pressure averaging
    ! coefficient ozone".
    ! The target vertical LMDZ grid is the grid of layer boundaries.
    ! The input data does not depend on longitude, but the pressure
    ! at LMDZ layers does.
    ! Therefore, the values on the LMDZ grid do depend on longitude.
    ! Regridding in pressure is done by averaging a step function of pressure.

    use dimens_m, only: iim, jjm, llm
    use dimphy, only: klon
    use netcdf95, only: nf95_inq_varid, handle_err
    use netcdf, only: nf90_get_var
    use grid_change, only: dyn_phy
    use numer_rec, only: assert
    use press_coefoz_m, only: press_in_edg
    use regr1_step_av_m, only: regr1_step_av
    use pressure_var, only: p3d

    integer, intent(in):: ncid ! NetCDF ID of the file
    character(len=*), intent(in):: name ! of the NetCDF variable
    integer, intent(in):: julien ! jour julien, 1 <= julien <= 360

    real, intent(out):: v3(:, :) ! (klon, llm)
    ! (ozone coefficient from Mobidic on the "physics" grid)
    ! ("v3(i, k)" is at longitude "xlon(i)", latitude
    ! "xlat(i)", in layer "k".)

    ! Variables local to the procedure:

    integer varid, ncerr ! for NetCDF

    real  v1(jjm + 1, size(press_in_edg) - 1)
    ! (ozone coefficient from "coefoz_LMDZ.nc" at day "julien")
    ! ("v1(j, k)" is at latitude "rlatu(j)" and for
    ! pressure interval "[press_in_edg(k), press_in_edg(k+1)]".)

    real v2(iim + 1, jjm + 1, llm)
    ! (ozone coefficient on the "dynamics" grid)
    ! ("v2(i, j, k)" is at longitude "rlonv(i)", latitude
    ! "rlatu(j)" and for pressure interval "[p3d(i, j, k+1), p3d(i, j, k)]".)

    integer i, j, k

    !--------------------------------------------

    call assert(shape(v3) == (/klon, llm/), "regr_pr_av_coefoz")

    call nf95_inq_varid(ncid, name, varid)

    ! Get data at the right day from the input file:
    ncerr = nf90_get_var(ncid, varid, v1, start=(/1, 1, julien/))
    call handle_err("regr_pr_av_coefoz nf90_get_var " // name, ncerr, ncid)
    ! Latitudes are in ascending order in the input file while
    ! "rlatu" is in descending order so we need to invert order:
    v1 = v1(jjm+1:1:-1, :)

    ! Regrid in pressure at each horizontal position:
    do j = 1, jjm + 1
       do i = 1, iim
          if (dyn_phy(i, j)) then
             v2(i, j, llm:1:-1) &
                  = regr1_step_av(v1(j, :), press_in_edg, &
                  p3d(i, j, llm+1:1:-1))
             ! (invert order of indices because "p3d" is in descending order)
          end if
       end do
    end do

    forall (k = 1:llm) v3(:, k) = pack(v2(:, :, k), dyn_phy)

  end subroutine regr_pr_av_coefoz

  !***************************************************************

  subroutine regr_pr_int_coefoz(ncid, name, julien, top_value, v3)

    ! "regr_pr_int_coefoz" stands for "regrid pressure interpolation
    ! coefficient ozone".
    ! The target vertical LMDZ grid is the grid of mid-layers.
    ! The input data does not depend on longitude, but the pressure
    ! at LMDZ mid-layers does.
    ! Therefore, the values on the LMDZ grid do depend on longitude.
    ! Regridding is by linear interpolation.

    use dimens_m, only: iim, jjm, llm
    use dimphy, only: klon
    use netcdf95, only: nf95_inq_varid, handle_err
    use netcdf, only: nf90_get_var
    use grid_change, only: dyn_phy
    use numer_rec, only: assert
    use press_coefoz_m, only: plev
    use regr1_lint_m, only: regr1_lint
    use pressure_var, only: pls

    integer, intent(in):: ncid ! NetCDF ID of the file
    character(len=*), intent(in):: name ! of the NetCDF variable
    integer, intent(in):: julien ! jour julien, 1 <= julien <= 360

    real, intent(in):: top_value
    ! (extra value of ozone coefficient at 0 pressure)

    real, intent(out):: v3(:, :) ! (klon, llm)
    ! (ozone coefficient from Mobidic on the "physics" grid)
    ! ("v3(i, k)" is at longitude "xlon(i)", latitude
    ! "xlat(i)", middle of layer "k".)

    ! Variables local to the procedure:

    integer varid, ncerr ! for NetCDF

    real  v1(jjm + 1, 0:size(plev))
    ! (ozone coefficient from "coefoz_LMDZ.nc" at day "julien")
    ! ("v1(j, k >=1)" is at latitude "rlatu(j)" and pressure "plev(k)".)

    real v2(iim + 1, jjm + 1, llm)
    ! (ozone coefficient on the "dynamics" grid)
    ! "v2(i, j, k)" is at longitude "rlonv(i)", latitude
    ! "rlatu(j)" and pressure "pls(i, j, k)".)

    integer i, j, k

    !--------------------------------------------

    call assert(shape(v3) == (/klon, llm/), "regr_pr_int_coefoz")

    call nf95_inq_varid(ncid, name, varid)

    ! Get data at the right day from the input file:
    ncerr = nf90_get_var(ncid, varid, v1(:, 1:), start=(/1, 1, julien/))
    call handle_err("regr_pr_int_coefoz nf90_get_var " // name, ncerr, ncid)
    ! Latitudes are in ascending order in the input file while
    ! "rlatu" is in descending order so we need to invert order:
    v1(:, 1:) = v1(jjm+1:1:-1, 1:)

    ! Complete "v1" with the value at 0 pressure:
    v1(:, 0) = top_value

    ! Regrid in pressure at each horizontal position:
    do j = 1, jjm + 1
       do i = 1, iim
          if (dyn_phy(i, j)) then
             v2(i, j, llm:1:-1) &
                  = regr1_lint(v1(j, :), (/0., plev/), pls(i, j, llm:1:-1))
             ! (invert order of indices because "pls" is in descending order)
          end if
       end do
    end do

    forall (k = 1:llm) v3(:, k) = pack(v2(:, :, k), dyn_phy)

  end subroutine regr_pr_int_coefoz

end module regr_pr_coefoz
