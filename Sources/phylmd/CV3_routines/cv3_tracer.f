module cv3_tracer_m

  implicit none

contains

  SUBROUTINE cv3_tracer(nloc, ncum, na, ment, sij, da, phi)

    implicit none

    ! Inputs:
    integer, intent(in):: ncum, na, nloc
    real ment(nloc, na, na), sij(nloc, na, na)

    ! Ouputs:
    real da(nloc, na), phi(nloc, na, na)

    ! Local variables:
    integer i, j, k

    !------------------------------------------------------------

    da = 0.

    do j = 1, na
       do k = 1, na
          do i = 1, ncum
             da(i, j) = da(i, j)+(1.-sij(i, k, j))*ment(i, k, j)
             phi(i, j, k) = sij(i, k, j) * ment(i, k, j)
          end do
       end do
    end do

  end SUBROUTINE cv3_tracer

end module cv3_tracer_m