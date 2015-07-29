SUBROUTINE a2c(u, v, imr, jmr, j1, j2, crx, cry, dtdx5, dtdy5)
  DIMENSION u(imr, *), v(imr, *), crx(imr, *), cry(imr, *), dtdx5(*)

  DO j = j1, j2
    DO i = 2, imr
      crx(i, j) = dtdx5(j)*(u(i,j)+u(i-1,j))
    END DO
  END DO

  DO j = j1, j2
    crx(1, j) = dtdx5(j)*(u(1,j)+u(imr,j))
  END DO

  DO i = 1, imr*jmr
    cry(i, 2) = dtdy5*(v(i,2)+v(i,1))
  END DO
  RETURN
END SUBROUTINE a2c

