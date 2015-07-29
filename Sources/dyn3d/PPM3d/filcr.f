SUBROUTINE filcr(q, imr, jnp, j1, j2, cosp, acosp, icr, tiny)
  DIMENSION q(imr, *), cosp(*), acosp(*)

  icr = 0
  DO j = j1 + 1, j2 - 1
    DO i = 1, imr - 1
      IF (q(i,j)<0.) THEN
        icr = 1
        dq = -q(i, j)*cosp(j)
        ! N-E
        dn = q(i+1, j+1)*cosp(j+1)
        d0 = max(0., dn)
        d1 = min(dq, d0)
        q(i+1, j+1) = (dn-d1)*acosp(j+1)
        dq = dq - d1
        ! S-E
        ds = q(i+1, j-1)*cosp(j-1)
        d0 = max(0., ds)
        d2 = min(dq, d0)
        q(i+1, j-1) = (ds-d2)*acosp(j-1)
        q(i, j) = (d2-dq)*acosp(j) + tiny
      END IF
    END DO
    IF (icr==0 .AND. q(imr,j)>=0.) GO TO 65
    DO i = 2, imr
      IF (q(i,j)<0.) THEN
        icr = 1
        dq = -q(i, j)*cosp(j)
        ! N-W
        dn = q(i-1, j+1)*cosp(j+1)
        d0 = max(0., dn)
        d1 = min(dq, d0)
        q(i-1, j+1) = (dn-d1)*acosp(j+1)
        dq = dq - d1
        ! S-W
        ds = q(i-1, j-1)*cosp(j-1)
        d0 = max(0., ds)
        d2 = min(dq, d0)
        q(i-1, j-1) = (ds-d2)*acosp(j-1)
        q(i, j) = (d2-dq)*acosp(j) + tiny
      END IF
    END DO
    ! *****************************************
    ! i=1
    i = 1
    IF (q(i,j)<0.) THEN
      icr = 1
      dq = -q(i, j)*cosp(j)
      ! N-W
      dn = q(imr, j+1)*cosp(j+1)
      d0 = max(0., dn)
      d1 = min(dq, d0)
      q(imr, j+1) = (dn-d1)*acosp(j+1)
      dq = dq - d1
      ! S-W
      ds = q(imr, j-1)*cosp(j-1)
      d0 = max(0., ds)
      d2 = min(dq, d0)
      q(imr, j-1) = (ds-d2)*acosp(j-1)
      q(i, j) = (d2-dq)*acosp(j) + tiny
    END IF
    ! *****************************************
    ! i=IMR
    i = imr
    IF (q(i,j)<0.) THEN
      icr = 1
      dq = -q(i, j)*cosp(j)
      ! N-E
      dn = q(1, j+1)*cosp(j+1)
      d0 = max(0., dn)
      d1 = min(dq, d0)
      q(1, j+1) = (dn-d1)*acosp(j+1)
      dq = dq - d1
      ! S-E
      ds = q(1, j-1)*cosp(j-1)
      d0 = max(0., ds)
      d2 = min(dq, d0)
      q(1, j-1) = (ds-d2)*acosp(j-1)
      q(i, j) = (d2-dq)*acosp(j) + tiny
    END IF
    ! *****************************************
65 END DO

  DO i = 1, imr
    IF (q(i,j1)<0. .OR. q(i,j2)<0.) THEN
      icr = 1
      GO TO 80
    END IF
  END DO

80 CONTINUE

  IF (q(1,1)<0. .OR. q(1,jnp)<0.) THEN
    icr = 1
  END IF

  RETURN
END SUBROUTINE filcr

