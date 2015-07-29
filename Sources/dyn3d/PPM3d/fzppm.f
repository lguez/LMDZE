SUBROUTINE fzppm(imr, jnp, nlay, j1, dq, wz, p, dc, dqdt, ar, al, a6, flux, &
    wk1, wk2, wz2, delp, kord)
  PARAMETER (kmax=150)
  PARAMETER (r23=2./3., r3=1./3.)
  REAL wz(imr, jnp, nlay), p(imr, jnp, nlay), dc(imr, jnp, nlay), &
    wk1(imr, *), delp(imr, jnp, nlay), dq(imr, jnp, nlay), &
    dqdt(imr, jnp, nlay)
  ! Assuming JNP >= NLAY
  REAL ar(imr, *), al(imr, *), a6(imr, *), flux(imr, *), wk2(imr, *), &
    wz2(imr, *)

  jmr = jnp - 1
  imjm = imr*jnp
  nlaym1 = nlay - 1

  lmt = kord - 3

  ! ****6***0*********0*********0*********0*********0*********0**********72
  ! Compute DC for PPM
  ! ****6***0*********0*********0*********0*********0*********0**********72

  DO k = 1, nlaym1
    DO i = 1, imjm
      dqdt(i, 1, k) = p(i, 1, k+1) - p(i, 1, k)
    END DO
  END DO

  DO k = 2, nlaym1
    DO i = 1, imjm
      c0 = delp(i, 1, k)/(delp(i,1,k-1)+delp(i,1,k)+delp(i,1,k+1))
      c1 = (delp(i,1,k-1)+0.5*delp(i,1,k))/(delp(i,1,k+1)+delp(i,1,k))
      c2 = (delp(i,1,k+1)+0.5*delp(i,1,k))/(delp(i,1,k-1)+delp(i,1,k))
      tmp = c0*(c1*dqdt(i,1,k)+c2*dqdt(i,1,k-1))
      qmax = max(p(i,1,k-1), p(i,1,k), p(i,1,k+1)) - p(i, 1, k)
      qmin = p(i, 1, k) - min(p(i,1,k-1), p(i,1,k), p(i,1,k+1))
      dc(i, 1, k) = sign(min(abs(tmp),qmax,qmin), tmp)
    END DO
  END DO


  ! ****6***0*********0*********0*********0*********0*********0**********72
  ! Loop over latitudes  (to save memory)
  ! ****6***0*********0*********0*********0*********0*********0**********72

  DO j = 1, jnp
    IF ((j==2 .OR. j==jmr) .AND. j1/=2) GO TO 2000

    DO k = 1, nlay
      DO i = 1, imr
        wz2(i, k) = wz(i, j, k)
        wk1(i, k) = p(i, j, k)
        wk2(i, k) = delp(i, j, k)
        flux(i, k) = dc(i, j, k) !this flux is actually the monotone slope
      END DO
    END DO

    ! ****6***0*********0*********0*********0*********0*********0**********72
    ! Compute first guesses at cell interfaces
    ! First guesses are required to be continuous.
    ! ****6***0*********0*********0*********0*********0*********0**********72

    ! three-cell parabolic subgrid distribution at model top
    ! two-cell parabolic with zero gradient subgrid distribution
    ! at the surface.

    ! First guess top edge value
    DO i = 1, imr
      ! three-cell PPM
      ! Compute a,b, and c of q = aP**2 + bP + c using cell averages and delp
      a = 3.*(dqdt(i,j,2)-dqdt(i,j,1)*(wk2(i,2)+wk2(i,3))/(wk2(i,1)+wk2(i, &
        2)))/((wk2(i,2)+wk2(i,3))*(wk2(i,1)+wk2(i,2)+wk2(i,3)))
      b = 2.*dqdt(i, j, 1)/(wk2(i,1)+wk2(i,2)) - r23*a*(2.*wk2(i,1)+wk2(i,2))
      al(i, 1) = wk1(i, 1) - wk2(i, 1)*(r3*a*wk2(i,1)+0.5*b)
      al(i, 2) = wk2(i, 1)*(a*wk2(i,1)+b) + al(i, 1)

      ! Check if change sign
      IF (wk1(i,1)*al(i,1)<=0.) THEN
        al(i, 1) = 0.
        flux(i, 1) = 0.
      ELSE
        flux(i, 1) = wk1(i, 1) - al(i, 1)
      END IF
    END DO

    ! Bottom
    DO i = 1, imr
      ! 2-cell PPM with zero gradient right at the surface

      fct = dqdt(i, j, nlaym1)*wk2(i, nlay)**2/((wk2(i,nlay)+wk2(i, &
        nlaym1))*(2.*wk2(i,nlay)+wk2(i,nlaym1)))
      ar(i, nlay) = wk1(i, nlay) + fct
      al(i, nlay) = wk1(i, nlay) - (fct+fct)
      IF (wk1(i,nlay)*ar(i,nlay)<=0.) ar(i, nlay) = 0.
      flux(i, nlay) = ar(i, nlay) - wk1(i, nlay)
    END DO


    ! ****6***0*********0*********0*********0*********0*********0**********72
    ! 4th order interpolation in the interior.
    ! ****6***0*********0*********0*********0*********0*********0**********72

    DO k = 3, nlaym1
      DO i = 1, imr
        c1 = dqdt(i, j, k-1)*wk2(i, k-1)/(wk2(i,k-1)+wk2(i,k))
        c2 = 2./(wk2(i,k-2)+wk2(i,k-1)+wk2(i,k)+wk2(i,k+1))
        a1 = (wk2(i,k-2)+wk2(i,k-1))/(2.*wk2(i,k-1)+wk2(i,k))
        a2 = (wk2(i,k)+wk2(i,k+1))/(2.*wk2(i,k)+wk2(i,k-1))
        al(i, k) = wk1(i, k-1) + c1 + c2*(wk2(i,k)*(c1*(a1-a2)+a2*flux(i, &
          k-1))-wk2(i,k-1)*a1*flux(i,k))
      END DO
    END DO

    DO i = 1, imr*nlaym1
      ar(i, 1) = al(i, 2)
    END DO

    DO i = 1, imr*nlay
      a6(i, 1) = 3.*(wk1(i,1)+wk1(i,1)-(al(i,1)+ar(i,1)))
    END DO

    ! ****6***0*********0*********0*********0*********0*********0**********72
    ! Top & Bot always monotonic
    CALL lmtppm(flux(1,1), a6(1,1), ar(1,1), al(1,1), wk1(1,1), imr, 0)
    CALL lmtppm(flux(1,nlay), a6(1,nlay), ar(1,nlay), al(1,nlay), &
      wk1(1,nlay), imr, 0)

    ! Interior depending on KORD
    IF (lmt<=2) CALL lmtppm(flux(1,2), a6(1,2), ar(1,2), al(1,2), wk1(1,2), &
      imr*(nlay-2), lmt)

    ! ****6***0*********0*********0*********0*********0*********0**********72

    DO i = 1, imr*nlaym1
      IF (wz2(i,1)>0.) THEN
        cm = wz2(i, 1)/wk2(i, 1)
        flux(i, 2) = ar(i, 1) + 0.5*cm*(al(i,1)-ar(i,1)+a6(i,1)*(1.-r23*cm))
      ELSE
        cp = wz2(i, 1)/wk2(i, 2)
        flux(i, 2) = al(i, 2) + 0.5*cp*(al(i,2)-ar(i,2)-a6(i,2)*(1.+r23*cp))
      END IF
    END DO

    DO i = 1, imr*nlaym1
      flux(i, 2) = wz2(i, 1)*flux(i, 2)
    END DO

    DO i = 1, imr
      dq(i, j, 1) = dq(i, j, 1) - flux(i, 2)
      dq(i, j, nlay) = dq(i, j, nlay) + flux(i, nlay)
    END DO

    DO k = 2, nlaym1
      DO i = 1, imr
        dq(i, j, k) = dq(i, j, k) + flux(i, k) - flux(i, k+1)
      END DO
    END DO
2000 END DO
  RETURN
END SUBROUTINE fzppm

