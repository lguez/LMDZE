
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/ismin.F,v 1.1.1.1 2004/05/19
! 12:53:07 lmdzadmin Exp $

FUNCTION ismin(n, sx, incx)

  IMPLICIT NONE

  INTEGER n, i, incx, ismin, ix
  REAL sx((n-1)*incx+1), sxmin

  ix = 1
  ismin = 1
  sxmin = sx(1)
  DO i = 1, n - 1
    ix = ix + incx
    IF (sx(ix)<sxmin) THEN
      sxmin = sx(ix)
      ismin = i + 1
    END IF
  END DO

  RETURN
END FUNCTION ismin

