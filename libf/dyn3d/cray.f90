subroutine scopy(n, sx, incx, sy, incy)

  ! From cray.F, v 1.1.1.1 2004/05/19 12:53:05

  ! This subroutine should not exist in a Fortran 95 program. If the
  ! actual arguments are of rank 1 then replace each call to this
  ! subroutine by the simple statement in this subroutine.

  IMPLICIT NONE

  integer, intent(in):: n, incx, incy
  real, intent(in):: sx((n-1)*incx+1)
  real, intent(inout):: sy((n-1)*incy+1)

  !-------------------------

  sy(::incy) = sx(::incx)

end subroutine scopy

!***********************************

real function ssum(n,sx,incx)

  ! From cray.F, v 1.1.1.1 2004/05/19 12:53:05

  IMPLICIT NONE

  integer, intent(in):: n, incx
  real, intent(in):: sx((n-1)*incx+1)

  !-----------------------

  ssum=sum(sx(::incx))

end function ssum
