SUBROUTINE cosa(cosp, cose, jnp, pi, dp)
  DIMENSION cosp(*), cose(*)

  jmr = jnp - 1
  DO j = 2, jnp
    ph5 = -0.5*pi + (float(j-1)-0.5)*dp
    cose(j) = cos(ph5)
  END DO

  jeq = (jnp+1)/2
  IF (jmr==2*(jmr/2)) THEN
    DO j = jnp, jeq + 1, -1
      cose(j) = cose(jnp+2-j)
    END DO
  ELSE
    ! cell edge at equator.
    cose(jeq+1) = 1.
    DO j = jnp, jeq + 2, -1
      cose(j) = cose(jnp+2-j)
    END DO
  END IF

  DO j = 2, jmr
    cosp(j) = 0.5*(cose(j)+cose(j+1))
  END DO
  cosp(1) = 0.
  cosp(jnp) = 0.
  RETURN
END SUBROUTINE cosa

