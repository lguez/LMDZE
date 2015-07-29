SUBROUTINE fyppm(vc, p, dc, flux, imr, jnp, j1, j2, a6, ar, al, jord)
  PARAMETER (r3=1./3., r23=2./3.)
  REAL vc(imr, *), flux(imr, *), p(imr, *), dc(imr, *)
  ! Local work arrays.
  REAL ar(imr, jnp), al(imr, jnp), a6(imr, jnp)
  INTEGER lmt
  ! logical first
  ! data first /.true./
  ! SAVE LMT

  imh = imr/2
  jmr = jnp - 1
  j11 = j1 - 1
  imjm1 = imr*(j2-j1+2)
  len = imr*(j2-j1+3)
  ! if(first) then
  ! IF(JORD.LE.0) then
  ! if(JMR.GE.90) then
  ! LMT = 0
  ! elseif(JMR.GE.45) then
  ! LMT = 1
  ! else
  ! LMT = 2
  ! endif
  ! else
  ! LMT = JORD - 3
  ! endif

  ! first = .false.
  ! endif

  ! modifs pour pouvoir choisir plusieurs schemas PPM
  lmt = jord - 3

  DO i = 1, imr*jmr
    al(i, 2) = 0.5*(p(i,1)+p(i,2)) + (dc(i,1)-dc(i,2))*r3
    ar(i, 1) = al(i, 2)
  END DO

  ! Poles:

  DO i = 1, imh
    al(i, 1) = al(i+imh, 2)
    al(i+imh, 1) = al(i, 2)

    ar(i, jnp) = ar(i+imh, jmr)
    ar(i+imh, jnp) = ar(i, jmr)
  END DO

  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Rajout pour LMDZ.3.3
  ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ar(imr, 1) = al(1, 1)
  ar(imr, jnp) = al(1, jnp)
  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  DO i = 1, len
    a6(i, j11) = 3.*(p(i,j11)+p(i,j11)-(al(i,j11)+ar(i,j11)))
  END DO

  IF (lmt<=2) CALL lmtppm(dc(1,j11), a6(1,j11), ar(1,j11), al(1,j11), &
    p(1,j11), len, lmt)


  DO i = 1, imjm1
    IF (vc(i,j1)>0.) THEN
      flux(i, j1) = ar(i, j11) + 0.5*vc(i, j1)*(al(i,j11)-ar(i,j11)+a6(i,j11) &
        *(1.-r23*vc(i,j1)))
    ELSE
      flux(i, j1) = al(i, j1) - 0.5*vc(i, j1)*(ar(i,j1)-al(i,j1)+a6(i,j1)*(1. &
        +r23*vc(i,j1)))
    END IF
  END DO
  RETURN
END SUBROUTINE fyppm

