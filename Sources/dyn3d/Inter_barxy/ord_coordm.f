module ord_coordm_m

  implicit none

contains


  !***********************************

  function ord_coordm(xi)

    ! From dyn3d/ord_coordm.F, version 1.1.1.1 2004/05/19 12:53:06
    ! Author : P. Le Van

    ! This procedure converts to degrees, if necessary, and inverts the
    ! order.

    use nr_util, only: pi


    REAL, intent(in):: xi(:) ! angle, in rad or degrees
    REAL ord_coordm(size(xi)) ! angle, in degrees

    !-----------------------------

    IF (xi(1) < 6.5) THEN
       ! "xi" is in rad
       ord_coordm(:) = xi(size(xi):1:-1) * 180. / pi
    else
       ! "xi" is in degrees
       ord_coordm(:) = xi(size(xi):1:-1)
    ENDIF

  END function ord_coordm

end module ord_coordm_m
