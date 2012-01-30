subroutine scopy(n, sx, incx, sy, incy)

  ! From cray.F, version 1.1.1.1, 2004/05/19 12:53:05

  ! This subroutine should not exist in a Fortran 95 program. If the
  ! actual arguments 2 and 4 are of rank 1 then replace each call to
  ! this subroutine by the simple statement in this subroutine.

  IMPLICIT NONE

  integer, intent(in):: n, incx, incy
  real, intent(in):: sx((n-1)*incx+1)
  real, intent(inout):: sy((n-1)*incy+1)

  !-------------------------

  sy(::incy) = sx(::incx)

end subroutine scopy
