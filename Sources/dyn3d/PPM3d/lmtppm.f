SUBROUTINE lmtppm(dc, a6, ar, al, p, im, lmt)

  ! A6 =  CURVATURE OF THE TEST PARABOLA
  ! AR =  RIGHT EDGE VALUE OF THE TEST PARABOLA
  ! AL =  LEFT  EDGE VALUE OF THE TEST PARABOLA
  ! DC =  0.5 * MISMATCH
  ! P  =  CELL-AVERAGED VALUE
  ! IM =  VECTOR LENGTH

  ! OPTIONS:

  ! LMT = 0: FULL MONOTONICITY
  ! LMT = 1: SEMI-MONOTONIC CONSTRAINT (NO UNDERSHOOTS)
  ! LMT = 2: POSITIVE-DEFINITE CONSTRAINT

  PARAMETER (r12=1./12.)
  DIMENSION a6(im), ar(im), al(im), p(im), dc(im)

  IF (lmt==0) THEN
    ! Full constraint
    DO i = 1, im
      IF (dc(i)==0.) THEN
        ar(i) = p(i)
        al(i) = p(i)
        a6(i) = 0.
      ELSE
        da1 = ar(i) - al(i)
        da2 = da1**2
        a6da = a6(i)*da1
        IF (a6da<-da2) THEN
          a6(i) = 3.*(al(i)-p(i))
          ar(i) = al(i) - a6(i)
        ELSE IF (a6da>da2) THEN
          a6(i) = 3.*(ar(i)-p(i))
          al(i) = ar(i) - a6(i)
        END IF
      END IF
    END DO
  ELSE IF (lmt==1) THEN
    ! Semi-monotonic constraint
    DO i = 1, im
      IF (abs(ar(i)-al(i))>=-a6(i)) GO TO 150
      IF (p(i)<ar(i) .AND. p(i)<al(i)) THEN
        ar(i) = p(i)
        al(i) = p(i)
        a6(i) = 0.
      ELSE IF (ar(i)>al(i)) THEN
        a6(i) = 3.*(al(i)-p(i))
        ar(i) = al(i) - a6(i)
      ELSE
        a6(i) = 3.*(ar(i)-p(i))
        al(i) = ar(i) - a6(i)
      END IF
150 END DO
  ELSE IF (lmt==2) THEN
    DO i = 1, im
      IF (abs(ar(i)-al(i))>=-a6(i)) GO TO 250
      fmin = p(i) + 0.25*(ar(i)-al(i))**2/a6(i) + a6(i)*r12
      IF (fmin>=0.) GO TO 250
      IF (p(i)<ar(i) .AND. p(i)<al(i)) THEN
        ar(i) = p(i)
        al(i) = p(i)
        a6(i) = 0.
      ELSE IF (ar(i)>al(i)) THEN
        a6(i) = 3.*(al(i)-p(i))
        ar(i) = al(i) - a6(i)
      ELSE
        a6(i) = 3.*(ar(i)-p(i))
        al(i) = ar(i) - a6(i)
      END IF
250 END DO
  END IF
  RETURN
END SUBROUTINE lmtppm

