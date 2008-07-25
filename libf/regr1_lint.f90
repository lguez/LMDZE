module regr1_lint_m

  implicit none

  interface regr1_lint
     ! Each procedure regrids by linear interpolation.
     ! The regridding operation is done on the first dimension of the
     ! input array.
     ! The difference betwwen the procedures is the rank of the first argument.
     module procedure regr11_lint, regr12_lint
  end interface

  private
  public regr1_lint

contains

  function regr11_lint(vs, xs, xt) result(vt)

    ! "vs" has rank 1.

    use numer_rec, only: assert_eq, hunt !!, polint

    real, intent(in):: vs(:)
    ! (values of the function at source points "xs")

    real, intent(in):: xs(:)
    ! (abscissas of points in source grid, in strictly monotonic order)

    real, intent(in):: xt(:)
    ! (abscissas of points in target grid)

    real vt(size(xt)) ! values of the function on the target grid

    ! Variables local to the procedure:
    integer is, it, ns
    integer is_b ! "is" bound between 1 and "ns - 1"

    !--------------------------------------

    ns = assert_eq(size(vs), size(xs), "regr11_lint ns")

    is = -1 ! go immediately to bisection on first call to "hunt"

    do it = 1, size(xt)
       call hunt(xs, xt(it), is)
       is_b = min(max(is, 1), ns - 1)
!!       call polint(xs(is_b:is_b+1), vs(is_b:is_b+1), xt(it), vt(it))
       vt(it) = ((xs(is_b+1) - xt(it)) * vs(is_b) &
            + (xt(it) - xs(is_b)) * vs(is_b+1)) / (xs(is_b+1) - xs(is_b))
    end do

  end function regr11_lint

  !*********************************************************

  function regr12_lint(vs, xs, xt) result(vt)

    ! "vs" has rank 2.

    use numer_rec, only: assert_eq, hunt

    real, intent(in):: vs(:, :)
    ! (values of the function at source points "xs")

    real, intent(in):: xs(:)
    ! (abscissas of points in source grid, in strictly monotonic order)

    real, intent(in):: xt(:)
    ! (abscissas of points in target grid)

    real vt(size(xt), size(vs, 2)) ! values of the function on the target grid

    ! Variables local to the procedure:
    integer is, it, ns
    integer is_b ! "is" bound between 1 and "ns - 1"

    !--------------------------------------

    ns = assert_eq(size(vs, 1), size(xs), "regr12_lint ns")

    is = -1 ! go immediately to bisection on first call to "hunt"

    do it = 1, size(xt)
       call hunt(xs, xt(it), is)
       is_b = min(max(is, 1), ns - 1)
       vt(it, :) = ((xs(is_b+1) - xt(it)) * vs(is_b, :) &
            + (xt(it) - xs(is_b)) * vs(is_b+1, :)) / (xs(is_b+1) - xs(is_b))
    end do

  end function regr12_lint

end module regr1_lint_m
