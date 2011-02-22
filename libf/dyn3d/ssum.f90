real function ssum(n,sx,incx)

  ! From cray.F, version 1.1.1.1 2004/05/19 12:53:05

  IMPLICIT NONE

  integer, intent(in):: n, incx
  real, intent(in):: sx((n - 1) * incx + 1)

  !-----------------------

  ssum = sum(sx(::incx))

end function ssum
