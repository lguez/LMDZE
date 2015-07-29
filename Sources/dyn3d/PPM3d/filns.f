SUBROUTINE filns(q, imr, jnp, j1, j2, cosp, acosp, ipy, tiny)
  DIMENSION q(imr, *), cosp(*), acosp(*)
  ! logical first
  ! data first /.true./
  ! save cap1

  ! if(first) then
  dp = 4.*atan(1.)/float(jnp-1)
  cap1 = imr*(1.-cos((j1-1.5)*dp))/dp
  ! first = .false.
  ! endif

  ipy = 0
  DO j = j1 + 1, j2 - 1
    DO i = 1, imr
      IF (q(i,j)<0.) THEN
        ipy = 1
        dq = -q(i, j)*cosp(j)
        ! North
        dn = q(i, j+1)*cosp(j+1)
        d0 = max(0., dn)
        d1 = min(dq, d0)
        q(i, j+1) = (dn-d1)*acosp(j+1)
        dq = dq - d1
        ! South
        ds = q(i, j-1)*cosp(j-1)
        d0 = max(0., ds)
        d2 = min(dq, d0)
        q(i, j-1) = (ds-d2)*acosp(j-1)
        q(i, j) = (d2-dq)*acosp(j) + tiny
      END IF
    END DO
  END DO

  DO i = 1, imr
    IF (q(i,j1)<0.) THEN
      ipy = 1
      dq = -q(i, j1)*cosp(j1)
      ! North
      dn = q(i, j1+1)*cosp(j1+1)
      d0 = max(0., dn)
      d1 = min(dq, d0)
      q(i, j1+1) = (dn-d1)*acosp(j1+1)
      q(i, j1) = (d1-dq)*acosp(j1) + tiny
    END IF
  END DO

  j = j2
  DO i = 1, imr
    IF (q(i,j)<0.) THEN
      ipy = 1
      dq = -q(i, j)*cosp(j)
      ! South
      ds = q(i, j-1)*cosp(j-1)
      d0 = max(0., ds)
      d2 = min(dq, d0)
      q(i, j-1) = (ds-d2)*acosp(j-1)
      q(i, j) = (d2-dq)*acosp(j) + tiny
    END IF
  END DO

  ! Check Poles.
  IF (q(1,1)<0.) THEN
    dq = q(1, 1)*cap1/float(imr)*acosp(j1)
    DO i = 1, imr
      q(i, 1) = 0.
      q(i, j1) = q(i, j1) + dq
      IF (q(i,j1)<0.) ipy = 1
    END DO
  END IF

  IF (q(1,jnp)<0.) THEN
    dq = q(1, jnp)*cap1/float(imr)*acosp(j2)
    DO i = 1, imr
      q(i, jnp) = 0.
      q(i, j2) = q(i, j2) + dq
      IF (q(i,j2)<0.) ipy = 1
    END DO
  END IF

  RETURN
END SUBROUTINE filns

