module principal_cshift_m

  IMPLICIT NONE

contains

  subroutine principal_cshift(is2, xlon, xprimm)

    ! Add or subtract 2 pi so that xlon is near [-pi, pi], then cshift
    ! so that xlon is in ascending order. Make the same cshift on
    ! xprimm.

    ! Libraries:
    use jumble, only: twopi

    use dynetat0_chosen_m, only: clon
    USE dimensions, ONLY: iim

    integer, intent(in):: is2
    real, intent(inout):: xlon(:), xprimm(:) ! (iim + 1)

    !-----------------------------------------------------

    if (is2 /= 0) then
       IF (clon <= 0.) THEN
          IF (is2 /= 1) THEN
             xlon(:is2 - 1) = xlon(:is2 - 1) + twopi
             xlon(:iim) = cshift(xlon(:iim), shift = is2 - 1)
             xprimm(:iim) = cshift(xprimm(:iim), shift = is2 - 1)
          END IF
       else
          xlon(is2 + 1:iim) = xlon(is2 + 1:iim) - twopi
          xlon(:iim) = cshift(xlon(:iim), shift = is2)
          xprimm(:iim) = cshift(xprimm(:iim), shift = is2)
       end IF
    end if

    xlon(iim + 1) = xlon(1) + twopi
    xprimm(iim + 1) = xprimm(1)

  end subroutine principal_cshift

end module principal_cshift_m
