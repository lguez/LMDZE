SUBROUTINE ytp(imr, jnp, j1, j2, acosp, rcap, dq, p, vc, dc2, ymass, fx, a6, &
    ar, al, jord)
  DIMENSION p(imr, jnp), vc(imr, jnp), ymass(imr, jnp), dc2(imr, jnp), &
    dq(imr, jnp), acosp(jnp)
  ! Work array
  DIMENSION fx(imr, jnp), ar(imr, jnp), al(imr, jnp), a6(imr, jnp)

  jmr = jnp - 1
  len = imr*(j2-j1+2)

  IF (jord==1) THEN
    DO i = 1, len
      jt = float(j1) - vc(i, j1)
      fx(i, j1) = p(i, jt)
    END DO
  ELSE

    CALL ymist(imr, jnp, j1, p, dc2, 4)

    IF (jord<=0 .OR. jord>=3) THEN

      CALL fyppm(vc, p, dc2, fx, imr, jnp, j1, j2, a6, ar, al, jord)

    ELSE
      DO i = 1, len
        jt = float(j1) - vc(i, j1)
        fx(i, j1) = p(i, jt) + (sign(1.,vc(i,j1))-vc(i,j1))*dc2(i, jt)
      END DO
    END IF
  END IF

  DO i = 1, len
    fx(i, j1) = fx(i, j1)*ymass(i, j1)
  END DO

  DO j = j1, j2
    DO i = 1, imr
      dq(i, j) = dq(i, j) + (fx(i,j)-fx(i,j+1))*acosp(j)
    END DO
  END DO

  ! Poles
  sum1 = fx(imr, j1)
  sum2 = fx(imr, j2+1)
  DO i = 1, imr - 1
    sum1 = sum1 + fx(i, j1)
    sum2 = sum2 + fx(i, j2+1)
  END DO

  sum1 = dq(1, 1) - sum1*rcap
  sum2 = dq(1, jnp) + sum2*rcap
  DO i = 1, imr
    dq(i, 1) = sum1
    dq(i, jnp) = sum2
  END DO

  IF (j1/=2) THEN
    DO i = 1, imr
      dq(i, 2) = sum1
      dq(i, jmr) = sum2
    END DO
  END IF

  RETURN
END SUBROUTINE ytp

