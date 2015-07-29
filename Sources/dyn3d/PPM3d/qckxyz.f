SUBROUTINE qckxyz(q, qtmp, imr, jnp, nlay, j1, j2, cosp, acosp, cross, ic, &
    nstep)

  PARAMETER (tiny=1.E-60)
  DIMENSION q(imr, jnp, nlay), qtmp(imr, jnp), cosp(*), acosp(*)
  LOGICAL cross

  nlaym1 = nlay - 1
  len = imr*(j2-j1+1)
  ip = 0

  ! Top layer
  l = 1
  icr = 1
  CALL filns(q(1,1,l), imr, jnp, j1, j2, cosp, acosp, ipy, tiny)
  IF (ipy==0) GO TO 50
  CALL filew(q(1,1,l), qtmp, imr, jnp, j1, j2, ipx, tiny)
  IF (ipx==0) GO TO 50

  IF (cross) THEN
    CALL filcr(q(1,1,l), imr, jnp, j1, j2, cosp, acosp, icr, tiny)
  END IF
  IF (icr==0) GO TO 50

  ! Vertical filling...
  DO i = 1, len
    IF (q(i,j1,1)<0.) THEN
      ip = ip + 1
      q(i, j1, 2) = q(i, j1, 2) + q(i, j1, 1)
      q(i, j1, 1) = 0.
    END IF
  END DO

50 CONTINUE
  DO l = 2, nlaym1
    icr = 1

    CALL filns(q(1,1,l), imr, jnp, j1, j2, cosp, acosp, ipy, tiny)
    IF (ipy==0) GO TO 225
    CALL filew(q(1,1,l), qtmp, imr, jnp, j1, j2, ipx, tiny)
    IF (ipx==0) GO TO 225
    IF (cross) THEN
      CALL filcr(q(1,1,l), imr, jnp, j1, j2, cosp, acosp, icr, tiny)
    END IF
    IF (icr==0) GO TO 225

    DO i = 1, len
      IF (q(i,j1,l)<0.) THEN

        ip = ip + 1
        ! From above
        qup = q(i, j1, l-1)
        qly = -q(i, j1, l)
        dup = min(qly, qup)
        q(i, j1, l-1) = qup - dup
        q(i, j1, l) = dup - qly
        ! Below
        q(i, j1, l+1) = q(i, j1, l+1) + q(i, j1, l)
        q(i, j1, l) = 0.
      END IF
    END DO
225 END DO

  ! BOTTOM LAYER
  sum = 0.
  l = nlay

  CALL filns(q(1,1,l), imr, jnp, j1, j2, cosp, acosp, ipy, tiny)
  IF (ipy==0) GO TO 911
  CALL filew(q(1,1,l), qtmp, imr, jnp, j1, j2, ipx, tiny)
  IF (ipx==0) GO TO 911

  CALL filcr(q(1,1,l), imr, jnp, j1, j2, cosp, acosp, icr, tiny)
  IF (icr==0) GO TO 911

  DO i = 1, len
    IF (q(i,j1,l)<0.) THEN
      ip = ip + 1

      ! From above

      qup = q(i, j1, nlaym1)
      qly = -q(i, j1, l)
      dup = min(qly, qup)
      q(i, j1, nlaym1) = qup - dup
      ! From "below" the surface.
      sum = sum + qly - dup
      q(i, j1, l) = 0.
    END IF
  END DO

911 CONTINUE

  IF (ip>imr) THEN
    WRITE (6, *) 'IC=', ic, ' STEP=', nstep, ' Vertical filling pts=', ip
  END IF

  IF (sum>1.E-25) THEN
    WRITE (6, *) ic, nstep, ' Mass source from the ground=', sum
  END IF
  RETURN
END SUBROUTINE qckxyz

