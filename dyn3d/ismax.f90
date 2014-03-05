
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/ismax.F,v 1.1.1.1 2004/05/19
! 12:53:07 lmdzadmin Exp $

FUNCTION ismax(n, sx, incx)

  IMPLICIT NONE

  INTEGER n, i, incx, ismax, ix
  REAL sx((n-1)*incx+1), sxmax

  ix = 1
  ismax = 1
  sxmax = sx(1)
  DO i = 1, n - 1
    ix = ix + incx
    IF (sx(ix)>sxmax) THEN
      sxmax = sx(ix)
      ismax = i + 1
    END IF
  END DO

  RETURN
END FUNCTION ismax

