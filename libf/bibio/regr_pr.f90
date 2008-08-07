module regr_pr

  implicit none

contains

  function regr_pr_av(v, press_in)

    ! "regr_pr_av" stands for "regrid pressure averaging".
    ! This function regrids a 2D latitude -- pressure variable to the
    ! LMDZ 3D grid.
    ! The target horizontal LMDZ grid is the "scalar" grid: "rlonv", "rlatu".
    ! The target vertical LMDZ grid is the grid of layer boundaries.
    ! We assume that the variable is already on the LMDZ "rlatu" latitude grid.
    ! There only remains to regrid in pressure at each horizontal
    ! position.
    ! The input variable does not depend on longitude, but the pressure
    ! at LMDZ layers does.
    ! Therefore, the values on the LMDZ grid do depend on longitude.
    ! The variable is regridded by averaging.

    use dimens_m, only: iim, jjm, llm
    use numer_rec, only: assert
    use regr1_step_av_m, only: regr1_step_av
    use pressure_var, only: p3d
    use grid_change, only: dyn_phy

    real, intent(in):: v(:, :)
    ! ("v(j, l)" is at latitude "rlatu(j)" and for pressure interval
    ! "[press_in(l), press_in(l+1)]".)

    real, intent(in):: press_in(:)
    ! (edges of pressure intervals, on input grid, in Pa, in strictly
    ! increasing order)

    real regr_pr_av(iim + 1, jjm + 1, llm)
    ! (variable adapted to the LMDZ grid
    ! "regr_pr_av(i, j, l)" is at longitude "rlonv(i)", latitude
    ! "rlatu(j)" and in pressure interval "[p3d(i, j, l+1), p3d(i, j, l)]")

    ! Variables local to the procedure:
    integer i, j

    !---------------------------------------------

    call assert(size(v, 1) == jjm + 1, "regr_pr_av 1")
    call assert(size(press_in) == size(v, 2) + 1, "regr_pr_av 2")

    ! Regrid in pressure by averaging a step function of pressure.
    do j = 1, jjm + 1
       do i = 1, iim
          if (dyn_phy(i, j)) then
             regr_pr_av(i, j, llm:1:-1) &
                  = regr1_step_av(v(j, :), press_in, p3d(i, j, llm+1:1:-1))
             ! (invert order of indices because "p3d" is decreasing)
          end if
       end do
    end do

    ! Duplicate pole values on all longitudes:
    regr_pr_av(2:, 1, :) = spread(regr_pr_av(1, 1, :), dim=1, ncopies=iim)
    regr_pr_av(2:, jjm + 1, :) &
         = spread(regr_pr_av(1, jjm + 1, :), dim=1, ncopies=iim)

    ! Duplicate first longitude to last longitude:
    regr_pr_av(iim + 1, 2:jjm, :) = regr_pr_av(1, 2:jjm, :)

  end function regr_pr_av

  !************************************************************

  function regr_pr_int(v, press_in)

    ! "regr_pr_int" stands for "regrid pressure interpolation".
    ! This function regrids a 2D latitude -- pressure variable to the
    ! LMDZ 3D grid.
    ! The target horizontal LMDZ grid is the "scalar" grid: "rlonv", "rlatu".
    ! The target vertical LMDZ grid is the grid of mid-layers.
    ! We assume that the variable is already on the LMDZ latitude grid.
    ! There only remains to regrid in pressure at each horizontal
    ! position.
    ! The input variable does not depend on longitude, but the pressure
    ! at LMDZ mid-layers does.
    ! Therefore, the values on the LMDZ grid do depend on longitude.
    ! The variable is regridded by interpolation.

    use dimens_m, only: iim, jjm, llm
    use numer_rec, only: assert
    use regr1_lint_m, only: regr1_lint
    use pressure_var, only: pls
    use grid_change, only: dyn_phy

    real, intent(in):: v(:, :)
    ! ("v(j, l)" is at latitude "rlatu(j)" and pressure level "press_in(l)".)

    real, intent(in):: press_in(:)
    ! (pressure levels on input grid, in Pa, in strictly increasing order)

    real regr_pr_int(iim + 1, jjm + 1, llm)
    ! (variable adapted to the LMDZ grid
    ! "regr_pr_int(i, j, l)" is at longitude "rlonv(i)", latitude
    ! "rlatu(j)" and pressure level "pls(i, j, l)")

    ! Variables local to the procedure:
    integer i, j

    !---------------------------------------------

    call assert(size(v, 1) == jjm + 1, "regr_pr_int 1")
    call assert(size(press_in) == size(v, 2), "regr_pr_int 2")

    ! Regrid in pressure by linear interpolation
    do j = 1, jjm + 1
       do i = 1, iim
          if (dyn_phy(i, j)) then
             regr_pr_int(i, j, llm:1:-1) &
                  = regr1_lint(v(j, :), press_in, pls(i, j, llm:1:-1))
             ! (invert order of indices because "pls" is decreasing)
          end if
       end do
    end do

    ! Duplicate pole values on all longitudes:
    regr_pr_int(2:, 1, :) &
         = spread(regr_pr_int(1, 1, :), dim=1, ncopies=iim)
    regr_pr_int(2:, jjm + 1, :) &
         = spread(regr_pr_int(1, jjm + 1, :), dim=1, ncopies=iim)

    ! Duplicate first longitude to last longitude:
    regr_pr_int(iim + 1, 2:jjm, :) = regr_pr_int(1, 2:jjm, :)

  end function regr_pr_int

end module regr_pr
