SUBROUTINE ymist(imr, jnp, j1, p, dc, id)
  PARAMETER (r24=1./24.)
  DIMENSION p(imr, jnp), dc(imr, jnp)

  imh = imr/2
  jmr = jnp - 1
  ijm3 = imr*(jmr-3)

  IF (id==2) THEN
    DO i = 1, imr*(jmr-1)
      tmp = 0.25*(p(i,3)-p(i,1))
      pmax = max(p(i,1), p(i,2), p(i,3)) - p(i, 2)
      pmin = p(i, 2) - min(p(i,1), p(i,2), p(i,3))
      dc(i, 2) = sign(min(abs(tmp),pmin,pmax), tmp)
    END DO
  ELSE
    DO i = 1, imh
      ! J=2
      tmp = (8.*(p(i,3)-p(i,1))+p(i+imh,2)-p(i,4))*r24
      pmax = max(p(i,1), p(i,2), p(i,3)) - p(i, 2)
      pmin = p(i, 2) - min(p(i,1), p(i,2), p(i,3))
      dc(i, 2) = sign(min(abs(tmp),pmin,pmax), tmp)
      ! J=JMR
      tmp = (8.*(p(i,jnp)-p(i,jmr-1))+p(i,jmr-2)-p(i+imh,jmr))*r24
      pmax = max(p(i,jmr-1), p(i,jmr), p(i,jnp)) - p(i, jmr)
      pmin = p(i, jmr) - min(p(i,jmr-1), p(i,jmr), p(i,jnp))
      dc(i, jmr) = sign(min(abs(tmp),pmin,pmax), tmp)
    END DO
    DO i = imh + 1, imr
      ! J=2
      tmp = (8.*(p(i,3)-p(i,1))+p(i-imh,2)-p(i,4))*r24
      pmax = max(p(i,1), p(i,2), p(i,3)) - p(i, 2)
      pmin = p(i, 2) - min(p(i,1), p(i,2), p(i,3))
      dc(i, 2) = sign(min(abs(tmp),pmin,pmax), tmp)
      ! J=JMR
      tmp = (8.*(p(i,jnp)-p(i,jmr-1))+p(i,jmr-2)-p(i-imh,jmr))*r24
      pmax = max(p(i,jmr-1), p(i,jmr), p(i,jnp)) - p(i, jmr)
      pmin = p(i, jmr) - min(p(i,jmr-1), p(i,jmr), p(i,jnp))
      dc(i, jmr) = sign(min(abs(tmp),pmin,pmax), tmp)
    END DO

    DO i = 1, ijm3
      tmp = (8.*(p(i,4)-p(i,2))+p(i,1)-p(i,5))*r24
      pmax = max(p(i,2), p(i,3), p(i,4)) - p(i, 3)
      pmin = p(i, 3) - min(p(i,2), p(i,3), p(i,4))
      dc(i, 3) = sign(min(abs(tmp),pmin,pmax), tmp)
    END DO
  END IF

  IF (j1/=2) THEN
    DO i = 1, imr
      dc(i, 1) = 0.
      dc(i, jnp) = 0.
    END DO
  ELSE
    ! Determine slopes in polar caps for scalars!

    DO i = 1, imh
      ! South
      tmp = 0.25*(p(i,2)-p(i+imh,2))
      pmax = max(p(i,2), p(i,1), p(i+imh,2)) - p(i, 1)
      pmin = p(i, 1) - min(p(i,2), p(i,1), p(i+imh,2))
      dc(i, 1) = sign(min(abs(tmp),pmax,pmin), tmp)
      ! North.
      tmp = 0.25*(p(i+imh,jmr)-p(i,jmr))
      pmax = max(p(i+imh,jmr), p(i,jnp), p(i,jmr)) - p(i, jnp)
      pmin = p(i, jnp) - min(p(i+imh,jmr), p(i,jnp), p(i,jmr))
      dc(i, jnp) = sign(min(abs(tmp),pmax,pmin), tmp)
    END DO

    DO i = imh + 1, imr
      dc(i, 1) = -dc(i-imh, 1)
      dc(i, jnp) = -dc(i-imh, jnp)
    END DO
  END IF
  RETURN
END SUBROUTINE ymist

