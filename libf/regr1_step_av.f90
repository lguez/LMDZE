module regr1_step_av_m

  ! Author: Lionel GUEZ

  implicit none

  interface regr1_step_av

     ! Each procedure regrids a step function by averaging it.
     ! The regridding operation is done on the first dimension of the
     ! input array.
     ! Source grid contains edges of steps.
     ! Target grid contains positions of cell edges.
     ! The target grid should be included in the source grid: no
     ! extrapolation is allowed.
     ! The difference between the procedures is the rank of the first argument.

     module procedure regr11_step_av, regr12_step_av, regr13_step_av, &
          regr14_step_av
  end interface

  private
  public regr1_step_av

contains

  function regr11_step_av(vs, xs, xt) result(vt)

    ! "vs" has rank 1.

    use nr_util, only: assert_eq, assert
    use numer_rec, only: locate

    real, intent(in):: vs(:) ! values of steps on the source grid
    ! (Step "is" is between "xs(is)" and "xs(is + 1)".)

    real, intent(in):: xs(:)
    ! (edges of of steps on the source grid, in strictly increasing order)

    real, intent(in):: xt(:)
    ! (edges of cells of the target grid, in strictly increasing order)

    real vt(size(xt) - 1) ! average values on the target grid
    ! (Cell "it" is between "xt(it)" and "xt(it + 1)".)

    ! Variables local to the procedure:
    integer is, it, ns, nt
    real left_edge

    !---------------------------------------------

    ns = assert_eq(size(vs), size(xs) - 1, "regr11_step_av ns")
    nt = size(xt) - 1
    ! Quick check on sort order:
    call assert(xs(1) < xs(2), "regr11_step_av xs bad order")
    call assert(xt(1) < xt(2), "regr11_step_av xt bad order")

    call assert(xs(1) <= xt(1) .and. xt(nt + 1) <= xs(ns + 1), &
         "regr11_step_av extrapolation")

    is = locate(xs, xt(1)) ! 1 <= is <= ns, because we forbid extrapolation
    do it = 1, nt
       ! 1 <= is <= ns
       ! xs(is) <= xt(it) < xs(is + 1)
       ! Compute "vt(it)":
       left_edge = xt(it)
       vt(it) = 0.
       do while (xs(is + 1) < xt(it + 1))
          ! 1 <= is <= ns - 1
          vt(it) = vt(it) + (xs(is + 1) - left_edge) * vs(is)
          is = is + 1
          left_edge = xs(is)
       end do
       ! 1 <= is <= ns
       vt(it) = (vt(it) + (xt(it + 1) - left_edge) * vs(is)) &
            / (xt(it + 1) - xt(it))
       if (xs(is + 1) == xt(it + 1)) is = is + 1
       ! 1 <= is <= ns .or. it == nt
    end do

  end function regr11_step_av

  !********************************************

  function regr12_step_av(vs, xs, xt) result(vt)

    ! "vs" has rank 2.

    use nr_util, only: assert_eq, assert
    use numer_rec, only: locate

    real, intent(in):: vs(:, :) ! values of steps on the source grid
    ! (Step "is" is between "xs(is)" and "xs(is + 1)".)

    real, intent(in):: xs(:)
    ! (edges of steps on the source grid, in strictly increasing order)

    real, intent(in):: xt(:)
    ! (edges of cells of the target grid, in strictly increasing order)

    real vt(size(xt) - 1, size(vs, 2)) ! average values on the target grid
    ! (Cell "it" is between "xt(it)" and "xt(it + 1)".)

    ! Variables local to the procedure:
    integer is, it, ns, nt
    real left_edge

    !---------------------------------------------

    ns = assert_eq(size(vs, 1), size(xs) - 1, "regr12_step_av ns")
    nt = size(xt) - 1

    ! Quick check on sort order:
    call assert(xs(1) < xs(2), "regr12_step_av xs bad order")
    call assert(xt(1) < xt(2), "regr12_step_av xt bad order")

    call assert(xs(1) <= xt(1) .and. xt(nt + 1) <= xs(ns + 1), &
         "regr12_step_av extrapolation")

    is = locate(xs, xt(1)) ! 1 <= is <= ns, because we forbid extrapolation
    do it = 1, nt
       ! 1 <= is <= ns
       ! xs(is) <= xt(it) < xs(is + 1)
       ! Compute "vt(it, :)":
       left_edge = xt(it)
       vt(it, :) = 0.
       do while (xs(is + 1) < xt(it + 1))
          ! 1 <= is <= ns - 1
          vt(it, :) = vt(it, :) + (xs(is + 1) - left_edge) * vs(is, :)
          is = is + 1
          left_edge = xs(is)
       end do
       ! 1 <= is <= ns
       vt(it, :) = (vt(it, :) + (xt(it + 1) - left_edge) * vs(is, :)) &
            / (xt(it + 1) - xt(it))
       if (xs(is + 1) == xt(it + 1)) is = is + 1
       ! 1 <= is <= ns .or. it == nt
    end do

  end function regr12_step_av

  !********************************************

  function regr13_step_av(vs, xs, xt) result(vt)

    ! "vs" has rank 3.

    use nr_util, only: assert_eq, assert
    use numer_rec, only: locate

    real, intent(in):: vs(:, :, :) ! values of steps on the source grid
    ! (Step "is" is between "xs(is)" and "xs(is + 1)".)

    real, intent(in):: xs(:)
    ! (edges of steps on the source grid, in strictly increasing order)

    real, intent(in):: xt(:)
    ! (edges of cells of the target grid, in strictly increasing order)

    real vt(size(xt) - 1, size(vs, 2), size(vs, 3)) 
    ! (average values on the target grid)
    ! (Cell "it" is between "xt(it)" and "xt(it + 1)".)

    ! Variables local to the procedure:
    integer is, it, ns, nt
    real left_edge

    !---------------------------------------------

    ns = assert_eq(size(vs, 1), size(xs) - 1, "regr13_step_av ns")
    nt = size(xt) - 1

    ! Quick check on sort order:
    call assert(xs(1) < xs(2), "regr13_step_av xs bad order")
    call assert(xt(1) < xt(2), "regr13_step_av xt bad order")

    call assert(xs(1) <= xt(1) .and. xt(nt + 1) <= xs(ns + 1), &
         "regr13_step_av extrapolation")

    is = locate(xs, xt(1)) ! 1 <= is <= ns, because we forbid extrapolation
    do it = 1, nt
       ! 1 <= is <= ns
       ! xs(is) <= xt(it) < xs(is + 1)
       ! Compute "vt(it, :, :)":
       left_edge = xt(it)
       vt(it, :, :) = 0.
       do while (xs(is + 1) < xt(it + 1))
          ! 1 <= is <= ns - 1
          vt(it, :, :) = vt(it, :, :) + (xs(is + 1) - left_edge) * vs(is, :, :)
          is = is + 1
          left_edge = xs(is)
       end do
       ! 1 <= is <= ns
       vt(it, :, :) = (vt(it, :, :) &
            + (xt(it + 1) - left_edge) * vs(is, :, :)) / (xt(it + 1) - xt(it))
       if (xs(is + 1) == xt(it + 1)) is = is + 1
       ! 1 <= is <= ns .or. it == nt
    end do

  end function regr13_step_av

  !********************************************

  function regr14_step_av(vs, xs, xt) result(vt)

    ! "vs" has rank 4.

    use nr_util, only: assert_eq, assert
    use numer_rec, only: locate

    real, intent(in):: vs(:, :, :, :) ! values of steps on the source grid
    ! (Step "is" is between "xs(is)" and "xs(is + 1)".)

    real, intent(in):: xs(:)
    ! (edges of steps on the source grid, in strictly increasing order)

    real, intent(in):: xt(:)
    ! (edges of cells of the target grid, in strictly increasing order)

    real vt(size(xt) - 1, size(vs, 2), size(vs, 3), size(vs, 4))
    ! (average values on the target grid)
    ! (Cell "it" is between "xt(it)" and "xt(it + 1)".)

    ! Variables local to the procedure:
    integer is, it, ns, nt
    real left_edge

    !---------------------------------------------

    ns = assert_eq(size(vs, 1), size(xs) - 1, "regr14_step_av ns")
    nt = size(xt) - 1

    ! Quick check on sort order:
    call assert(xs(1) < xs(2), "regr14_step_av xs bad order")
    call assert(xt(1) < xt(2), "regr14_step_av xt bad order")

    call assert(xs(1) <= xt(1) .and. xt(nt + 1) <= xs(ns + 1), &
         "regr14_step_av extrapolation")

    is = locate(xs, xt(1)) ! 1 <= is <= ns, because we forbid extrapolation
    do it = 1, nt
       ! 1 <= is <= ns
       ! xs(is) <= xt(it) < xs(is + 1)
       ! Compute "vt(it, :, :, :)":
       left_edge = xt(it)
       vt(it, :, :, :) = 0.
       do while (xs(is + 1) < xt(it + 1))
          ! 1 <= is <= ns - 1
          vt(it, :, :, :) = vt(it, :, :, :) + (xs(is + 1) - left_edge) &
               * vs(is, :, :, :)
          is = is + 1
          left_edge = xs(is)
       end do
       ! 1 <= is <= ns
       vt(it, :, :, :) = (vt(it, :, :, :) + (xt(it + 1) - left_edge) &
            * vs(is, :, :, :)) / (xt(it + 1) - xt(it))
       if (xs(is + 1) == xt(it + 1)) is = is + 1
       ! 1 <= is <= ns .or. it == nt
    end do

  end function regr14_step_av

end module regr1_step_av_m
