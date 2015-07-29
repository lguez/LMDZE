SUBROUTINE filew(q, qtmp, imr, jnp, j1, j2, ipx, tiny)
  DIMENSION q(imr, *), qtmp(jnp, imr)

  ipx = 0
  ! Copy & swap direction for vectorization.
  DO i = 1, imr
    DO j = j1, j2
      qtmp(j, i) = q(i, j)
    END DO
  END DO

  DO i = 2, imr - 1
    DO j = j1, j2
      IF (qtmp(j,i)<0.) THEN
        ipx = 1
        ! west
        d0 = max(0., qtmp(j,i-1))
        d1 = min(-qtmp(j,i), d0)
        qtmp(j, i-1) = qtmp(j, i-1) - d1
        qtmp(j, i) = qtmp(j, i) + d1
        ! east
        d0 = max(0., qtmp(j,i+1))
        d2 = min(-qtmp(j,i), d0)
        qtmp(j, i+1) = qtmp(j, i+1) - d2
        qtmp(j, i) = qtmp(j, i) + d2 + tiny
      END IF
    END DO
  END DO

  i = 1
  DO j = j1, j2
    IF (qtmp(j,i)<0.) THEN
      ipx = 1
      ! west
      d0 = max(0., qtmp(j,imr))
      d1 = min(-qtmp(j,i), d0)
      qtmp(j, imr) = qtmp(j, imr) - d1
      qtmp(j, i) = qtmp(j, i) + d1
      ! east
      d0 = max(0., qtmp(j,i+1))
      d2 = min(-qtmp(j,i), d0)
      qtmp(j, i+1) = qtmp(j, i+1) - d2

      qtmp(j, i) = qtmp(j, i) + d2 + tiny
    END IF
  END DO
  i = imr
  DO j = j1, j2
    IF (qtmp(j,i)<0.) THEN
      ipx = 1
      ! west
      d0 = max(0., qtmp(j,i-1))
      d1 = min(-qtmp(j,i), d0)
      qtmp(j, i-1) = qtmp(j, i-1) - d1
      qtmp(j, i) = qtmp(j, i) + d1
      ! east
      d0 = max(0., qtmp(j,1))
      d2 = min(-qtmp(j,i), d0)
      qtmp(j, 1) = qtmp(j, 1) - d2

      qtmp(j, i) = qtmp(j, i) + d2 + tiny
    END IF
  END DO

  IF (ipx/=0) THEN
    DO j = j1, j2
      DO i = 1, imr
        q(i, j) = qtmp(j, i)
      END DO
    END DO
  ELSE

    ! Poles.
    IF (q(1,1)<0 .OR. q(1,jnp)<0.) ipx = 1
  END IF
  RETURN
END SUBROUTINE filew
