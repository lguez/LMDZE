
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/zilch.F,v 1.1.1.1 2004/05/19
! 12:53:09 lmdzadmin Exp $

SUBROUTINE zilch(x, m)

  ! Zero the real array x dimensioned m.

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: m
  INTEGER i
  REAL x(m)

  DO i = 1, m
    x(i) = 0.0
  END DO
  RETURN
END SUBROUTINE zilch
