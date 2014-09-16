
! $Header: /home/cvsroot/LMDZ4/libf/filtrez/acc.F,v 1.1.1.1 2004/05/19
! 12:53:09 lmdzadmin Exp $

SUBROUTINE acc(vec, d, im)
  implicit none
  integer im, i, j
  real vec(im, im), d(im)
  real sum, ssum

  DO j = 1, im
    DO i = 1, im
      d(i) = vec(i, j)*vec(i, j)
    END DO
    sum = ssum(im, d, 1)
    sum = sqrt(sum)
    DO i = 1, im
      vec(i, j) = vec(i, j)/sum
    END DO
  END DO
  RETURN
END SUBROUTINE acc
