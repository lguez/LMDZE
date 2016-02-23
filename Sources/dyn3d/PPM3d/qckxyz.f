SUBROUTINE qckxyz(q, qtmp, imr, jnp, nlay, j1, j2, cosp, acosp, cross, ic, &
     nstep)

  implicit none

  real, PARAMETER:: my_tiny=tiny(0.)
  integer imr, jnp, nlay, j1, j2, ic, nstep
  real q(imr, jnp, nlay), qtmp(imr, jnp), cosp(*), acosp(*)
  LOGICAL cross
  real dup, qly, qup, sum
  integer i, icr, ip, ipx, ipy, l, len, nlaym1

  nlaym1 = nlay - 1
  len = imr*(j2-j1+1)
  ip = 0

  ! Top layer
  l = 1
  icr = 1
  CALL filns(q(1,1,l), imr, jnp, j1, j2, cosp, acosp, ipy, my_tiny)
  IF (ipy/=0) then
     CALL filew(q(1,1,l), qtmp, imr, jnp, j1, j2, ipx, my_tiny)
     IF (ipx/=0) then

        IF (cross) THEN
           CALL filcr(q(1,1,l), imr, jnp, j1, j2, cosp, acosp, icr, my_tiny)
        END IF
        IF (icr/=0) then

           ! Vertical filling...
           DO i = 1, len
              IF (q(i,j1,1)<0.) THEN
                 ip = ip + 1
                 q(i, j1, 2) = q(i, j1, 2) + q(i, j1, 1)
                 q(i, j1, 1) = 0.
              END IF
           END DO
        end IF
     end IF
  end IF

  DO l = 2, nlaym1
     icr = 1

     CALL filns(q(1,1,l), imr, jnp, j1, j2, cosp, acosp, ipy, my_tiny)
     IF (ipy/=0) then
        CALL filew(q(1,1,l), qtmp, imr, jnp, j1, j2, ipx, my_tiny)
        IF (ipx/=0) then
           IF (cross) THEN
              CALL filcr(q(1,1,l), imr, jnp, j1, j2, cosp, acosp, icr, my_tiny)
           END IF
           IF (icr/=0) then

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
           end IF
        end IF
     end IF
  END DO

  ! BOTTOM LAYER
  sum = 0.
  l = nlay

  CALL filns(q(1,1,l), imr, jnp, j1, j2, cosp, acosp, ipy, my_tiny)
  IF (ipy/=0) then
     CALL filew(q(1,1,l), qtmp, imr, jnp, j1, j2, ipx, my_tiny)
     IF (ipx/=0) then

        CALL filcr(q(1,1,l), imr, jnp, j1, j2, cosp, acosp, icr, my_tiny)
        IF (icr/=0) then

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
        end IF
     end IF
  end IF

  IF (ip>imr) THEN
     WRITE (6, *) 'IC=', ic, ' STEP=', nstep, ' Vertical filling pts=', ip
  END IF

  IF (sum>1.E-25) THEN
     WRITE (6, *) ic, nstep, ' Mass source from the ground=', sum
  END IF

END SUBROUTINE qckxyz
