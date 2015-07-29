SUBROUTINE cosc(cosp, cose, jnp, pi, dp)
  DIMENSION cosp(*), cose(*)

  phi = -0.5*pi
  DO j = 2, jnp - 1
    phi = phi + dp
    cosp(j) = cos(phi)
  END DO
  cosp(1) = 0.
  cosp(jnp) = 0.

  DO j = 2, jnp
    cose(j) = 0.5*(cosp(j)+cosp(j-1))
  END DO

  DO j = 2, jnp - 1
    cosp(j) = 0.5*(cose(j)+cose(j+1))
  END DO
  RETURN
END SUBROUTINE cosc

