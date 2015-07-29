SUBROUTINE yadv(imr, jnp, j1, j2, p, va, ady, wk, iad)
  REAL p(imr, jnp), ady(imr, jnp), va(imr, jnp)
  REAL wk(imr, -1:jnp+2)

  jmr = jnp - 1
  imh = imr/2
  DO j = 1, jnp
    DO i = 1, imr
      wk(i, j) = p(i, j)
    END DO
  END DO
  ! Poles:
  DO i = 1, imh
    wk(i, -1) = p(i+imh, 3)
    wk(i+imh, -1) = p(i, 3)
    wk(i, 0) = p(i+imh, 2)
    wk(i+imh, 0) = p(i, 2)
    wk(i, jnp+1) = p(i+imh, jmr)
    wk(i+imh, jnp+1) = p(i, jmr)
    wk(i, jnp+2) = p(i+imh, jnp-2)
    wk(i+imh, jnp+2) = p(i, jnp-2)
  END DO

  IF (iad==2) THEN
    DO j = j1 - 1, j2 + 1
      DO i = 1, imr
        jp = nint(va(i,j))
        rv = jp - va(i, j)
        jp = j - jp
        a1 = 0.5*(wk(i,jp+1)+wk(i,jp-1)) - wk(i, jp)
        b1 = 0.5*(wk(i,jp+1)-wk(i,jp-1))
        ady(i, j) = wk(i, jp) + rv*(a1*rv+b1) - wk(i, j)
      END DO
    END DO

  ELSE IF (iad==1) THEN
    DO j = j1 - 1, j2 + 1
      DO i = 1, imr
        jp = float(j) - va(i, j)
        ady(i, j) = va(i, j)*(wk(i,jp)-wk(i,jp+1))
      END DO
    END DO
  END IF

  IF (j1/=2) THEN
    sum1 = 0.
    sum2 = 0.
    DO i = 1, imr
      sum1 = sum1 + ady(i, 2)
      sum2 = sum2 + ady(i, jmr)
    END DO
    sum1 = sum1/imr
    sum2 = sum2/imr

    DO i = 1, imr
      ady(i, 2) = sum1
      ady(i, jmr) = sum2
      ady(i, 1) = sum1
      ady(i, jnp) = sum2
    END DO
  ELSE
    ! Poles:
    sum1 = 0.
    sum2 = 0.
    DO i = 1, imr
      sum1 = sum1 + ady(i, 1)
      sum2 = sum2 + ady(i, jnp)
    END DO
    sum1 = sum1/imr
    sum2 = sum2/imr

    DO i = 1, imr
      ady(i, 1) = sum1
      ady(i, jnp) = sum2
    END DO
  END IF

  RETURN
END SUBROUTINE yadv

