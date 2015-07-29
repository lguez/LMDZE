SUBROUTINE fxppm(imr, iml, ut, p, dc, flux, iord)
  PARAMETER (r3=1./3., r23=2./3.)
  DIMENSION ut(*), flux(*), p(-iml:imr+iml+1), dc(-iml:imr+iml+1)
  DIMENSION ar(0:imr), al(0:imr), a6(0:imr)
  INTEGER lmt
  ! logical first
  ! data first /.true./
  ! SAVE LMT
  ! if(first) then

  ! correction calcul de LMT a chaque passage pour pouvoir choisir
  ! plusieurs schemas PPM pour differents traceurs
  ! IF (IORD.LE.0) then
  ! if(IMR.GE.144) then
  ! LMT = 0
  ! elseif(IMR.GE.72) then
  ! LMT = 1
  ! else
  ! LMT = 2
  ! endif
  ! else
  ! LMT = IORD - 3
  ! endif

  lmt = iord - 3

  DO i = 1, imr
    al(i) = 0.5*(p(i-1)+p(i)) + (dc(i-1)-dc(i))*r3
  END DO

  DO i = 1, imr - 1
    ar(i) = al(i+1)
  END DO
  ar(imr) = al(1)

  DO i = 1, imr
    a6(i) = 3.*(p(i)+p(i)-(al(i)+ar(i)))
  END DO

  IF (lmt<=2) CALL lmtppm(dc(1), a6(1), ar(1), al(1), p(1), imr, lmt)

  al(0) = al(imr)
  ar(0) = ar(imr)
  a6(0) = a6(imr)

  DO i = 1, imr
    IF (ut(i)>0.) THEN
      flux(i) = ar(i-1) + 0.5*ut(i)*(al(i-1)-ar(i-1)+a6(i-1)*(1.-r23*ut(i)))
    ELSE
      flux(i) = al(i) - 0.5*ut(i)*(ar(i)-al(i)+a6(i)*(1.+r23*ut(i)))
    END IF
  END DO
  RETURN
END SUBROUTINE fxppm

