module cv30_tracer_m

  implicit none

contains

  SUBROUTINE cv30_tracer(nloc, ncum, na, ment, sij, da, phi)

    ! Passive tracers.

    integer, intent(in):: ncum, na, nloc
    real, intent(in):: ment(nloc, na, na), sij(nloc, na, na)

    ! Ouputs:
    real, intent(out):: da(nloc, na)
    real phi(nloc, na, na)

    ! Local:
    integer i, j, k

    !------------------------------------------------------------

    da = 0.

    do j = 1, na
       do k = 1, na
          do i = 1, ncum
             da(i, j) = da(i, j) + (1. - sij(i, k, j)) * ment(i, k, j)
             phi(i, j, k) = sij(i, k, j) * ment(i, k, j)
          end do
       end do
    end do

  end SUBROUTINE cv30_tracer

end module cv30_tracer_m
