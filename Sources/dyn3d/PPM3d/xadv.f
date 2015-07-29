SUBROUTINE xadv(imr, jnp, j1, j2, p, ua, js, jn, iml, adx, iad)
  REAL p(imr, jnp), adx(imr, jnp), qtmp(-imr:imr+imr), ua(imr, jnp)

  jmr = jnp - 1
  DO j = j1, j2
    IF (j>js .AND. j<jn) GO TO 1309

    DO i = 1, imr
      qtmp(i) = p(i, j)
    END DO

    DO i = -iml, 0
      qtmp(i) = p(imr+i, j)
      qtmp(imr+1-i) = p(1-i, j)
    END DO

    IF (iad==2) THEN
      DO i = 1, imr
        ip = nint(ua(i,j))
        ru = ip - ua(i, j)
        ip = i - ip
        a1 = 0.5*(qtmp(ip+1)+qtmp(ip-1)) - qtmp(ip)
        b1 = 0.5*(qtmp(ip+1)-qtmp(ip-1))
        adx(i, j) = qtmp(ip) + ru*(a1*ru+b1)
      END DO
    ELSE IF (iad==1) THEN
      DO i = 1, imr
        iu = ua(i, j)
        ru = ua(i, j) - iu
        iiu = i - iu
        IF (ua(i,j)>=0.) THEN
          adx(i, j) = qtmp(iiu) + ru*(qtmp(iiu-1)-qtmp(iiu))
        ELSE
          adx(i, j) = qtmp(iiu) + ru*(qtmp(iiu)-qtmp(iiu+1))
        END IF
      END DO
    END IF

    DO i = 1, imr
      adx(i, j) = adx(i, j) - p(i, j)
    END DO
1309 END DO

  ! Eulerian upwind

  DO j = js + 1, jn - 1

    DO i = 1, imr
      qtmp(i) = p(i, j)
    END DO

    qtmp(0) = p(imr, j)
    qtmp(imr+1) = p(1, j)

    IF (iad==2) THEN
      qtmp(-1) = p(imr-1, j)
      qtmp(imr+2) = p(2, j)
      DO i = 1, imr
        ip = nint(ua(i,j))
        ru = ip - ua(i, j)
        ip = i - ip
        a1 = 0.5*(qtmp(ip+1)+qtmp(ip-1)) - qtmp(ip)
        b1 = 0.5*(qtmp(ip+1)-qtmp(ip-1))
        adx(i, j) = qtmp(ip) - p(i, j) + ru*(a1*ru+b1)
      END DO
    ELSE IF (iad==1) THEN
      ! 1st order
      DO i = 1, imr
        ip = i - ua(i, j)
        adx(i, j) = ua(i, j)*(qtmp(ip)-qtmp(ip+1))
      END DO
    END IF
  END DO

  IF (j1/=2) THEN
    DO i = 1, imr
      adx(i, 2) = 0.
      adx(i, jmr) = 0.
    END DO
  END IF
  ! set cross term due to x-adv at the poles to zero.
  DO i = 1, imr
    adx(i, 1) = 0.
    adx(i, jnp) = 0.
  END DO
  RETURN
END SUBROUTINE xadv

