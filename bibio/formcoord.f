
! $Header: /home/cvsroot/LMDZ4/libf/bibio/formcoord.F,v 1.1.1.1 2004/05/19
! 12:53:05 lmdzadmin Exp $

SUBROUTINE formcoord(unit, n, x, a, rev, text)
  IMPLICIT NONE
  INTEGER n, unit, ndec
  LOGICAL rev
  REAL x(n), a
  CHARACTER(len=4) text

  INTEGER i, id, i1, i2, in
  REAL dx, dxmin

  IF (rev) THEN
    id = -1
    i1 = n
    i2 = n - 1
    in = 1
    WRITE (unit, 3000) text(1:1)
  ELSE
    id = 1
    i1 = 1
    i2 = 2
    in = n
  END IF

  IF (n<2) THEN
    ndec = 1
    WRITE (unit, 1000) text, n, x(1)*a
  ELSE
    dxmin = abs(x(2)-x(1))
    DO i = 2, n - 1
      dx = abs(x(i+1)-x(i))
      IF (dx<dxmin) dxmin = dx
    END DO

    ndec = -log10(dxmin) + 2
    IF (mod(n,6)==1) THEN
      WRITE (unit, 1000) text, n, x(i1)*a
      WRITE (unit, 2000)(x(i)*a, i=i2, in, id)
    ELSE
      WRITE (unit, 1000) text, n
      WRITE (unit, 2000)(x(i)*a, i=i1, in, id)
    END IF
  END IF

1000 FORMAT (A4, 2X, I4, ' LEVELS', 43X, F12.2)
2000 FORMAT (6F12.2)
  ! 1000  format(a4,2x,i4,' LEVELS',43x,f12.<ndec>)
  ! 2000  format(6f12.<ndec>)
3000 FORMAT ('FORMAT ', A1, 'REV')
  RETURN

END SUBROUTINE formcoord
