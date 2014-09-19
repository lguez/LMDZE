
! $Header: /home/cvsroot/LMDZ4/libf/filtrez/eigen_sort.F,v 1.1.1.1 2004/05/19
! 12:53:09 lmdzadmin Exp $

SUBROUTINE eigen_sort(d, v, n, np)
  INTEGER n, np
  REAL d(np), v(np, np)
  INTEGER i, j, k
  REAL p

  DO i = 1, n - 1
    k = i
    p = d(i)
    DO j = i + 1, n
      IF (d(j)>=p) THEN
        k = j
        p = d(j)
      END IF
    END DO

    IF (k/=i) THEN
      d(k) = d(i)
      d(i) = p
      DO j = 1, n
        p = v(j, i)
        v(j, i) = v(j, k)
        v(j, k) = p
      END DO
    END IF
  END DO

  RETURN
END SUBROUTINE eigen_sort