module regr_pr_int_m

  ! Author: Lionel GUEZ

  implicit none

contains

  subroutine regr_pr_int(ncid, name, julien, pplay, top_value, v3)

    ! "regr_pr_int" stands for "regrid pressure interpolation".

    ! This procedure reads a 2D latitude-pressure field from a NetCDF
    ! file, at a given day, regrids the field in pressure to the LMDZ
    ! vertical grid and packs it to the LMDZ horizontal "physics"
    ! grid. We assume that, in the input file, the field has 3
    ! dimensions: latitude, pressure, julian day. We assume that
    ! latitudes are in ascending order in the input file. We assume
    ! that the input data is already on the LMDZ "rlatu" latitude
    ! grid.

    ! The target vertical LMDZ grid is the grid of mid-layers.
    ! The input data does not depend on longitude, but the pressure
    ! at LMDZ mid-layers does.
    ! Therefore, the values on the LMDZ grid do depend on longitude.
    ! Regridding is by linear interpolation.

    use dimens_m, only: iim, jjm, llm
    use dimphy, only: klon
    use grid_change, only: gr_dyn_phy
    use netcdf95, only: nf95_inq_varid, nf95_get_var
    use nr_util, only: assert
    use numer_rec_95, only: regr1_lint
    use press_coefoz_m, only: plev

    integer, intent(in):: ncid ! NetCDF ID of the file
    character(len=*), intent(in):: name ! of the NetCDF variable
    integer, intent(in):: julien ! jour julien, 1 <= julien <= 360

    real, intent(in):: pplay(:, :) ! (klon, llm)
    ! (pression pour le mileu de chaque couche, en Pa)

    real, intent(in):: top_value
    ! (extra value of field at 0 pressure)

    real, intent(out):: v3(:, :) ! (klon, llm)
    ! Regridded field on the partial "physics" grid. "v3(i, k)" is at
    ! longitude "xlon(i)", latitude "xlat(i)", middle of layer "k".

    ! Variables local to the procedure:

    integer varid ! for NetCDF

    real v1(jjm + 1, 0:size(plev))
    ! Input field at day "julien". "v1(j, k >=1)" is at latitude
    ! "rlatu(j)" and pressure "plev(k)".

    real v2(klon, 0:size(plev))
    ! Field on the "physics" horizontal grid.  "v2(i, k >= 1)" is at
    ! longitude "xlon(i)", latitude "xlat(i)" and pressure "plev(k)".)

    integer i

    !--------------------------------------------

    call assert(shape(v3) == (/klon, llm/), "regr_pr_int")

    call nf95_inq_varid(ncid, name, varid)

    ! Get data at the right day from the input file:
    call nf95_get_var(ncid, varid, v1(:, 1:), start=(/1, 1, julien/))
    ! Latitudes are in ascending order in the input file while
    ! "rlatu" is in descending order so we need to invert order:
    v1(:, 1:) = v1(jjm+1:1:-1, 1:)

    ! Complete "v1" with the value at 0 pressure:
    v1(:, 0) = top_value

    v2 = gr_dyn_phy(spread(v1, dim = 1, ncopies = iim + 1))

    ! Regrid in pressure at each horizontal position:
    do i = 1, klon
       v3(i, llm:1:-1) = regr1_lint(v2(i, :), (/0., plev/), pplay(i, llm:1:-1))
       ! (invert order of indices because "pplay" is in descending order)
    end do

  end subroutine regr_pr_int

end module regr_pr_int_m
