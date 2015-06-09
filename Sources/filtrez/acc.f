module acc_m

  implicit none

contains

  SUBROUTINE acc(vec)

    ! From LMDZ4/libf/filtrez/acc.F, v 1.1.1.1 2004/05/19 12:53:09

    ! Normalize each column of vec.

    real, intent(inout):: vec(:, :)

    ! Local:
    integer j

    !--------------------------------------------------

    forall (j = 1:size(vec, 2)) vec(:, j) = vec(:, j) / sqrt(sum(vec(:, j)**2))

  END SUBROUTINE acc

end module acc_m
