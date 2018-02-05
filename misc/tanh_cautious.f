module tanh_cautious_m

  implicit none

contains

  elemental double precision function tanh_cautious(fa, fb)

    ! Compute tanh(fa / fb) with caution.

    double precision, intent(in):: fa, fb

    !-------------------------------------------------------------

    IF (200d0 * fb < - fa) THEN
       tanh_cautious = - 1d0
    ELSE IF (200d0 * fb < fa) THEN
       tanh_cautious = 1d0
    ELSE
       IF (ABS(fa) < 1d-13 .AND. ABS(fb) < 1d-13) THEN
          IF (200d0 * fb + fa < 1d-10) THEN
             tanh_cautious = - 1d0
          ELSE IF (200d0 * fb - fa < 1d-10) THEN
             tanh_cautious = 1d0
          END IF
       ELSE
          tanh_cautious = TANH(fa / fb)
       END IF
    END IF

  end function tanh_cautious

end module tanh_cautious_m
