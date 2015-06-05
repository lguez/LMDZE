module filtreg_hemisph_m

  implicit none

contains

  subroutine filtreg_hemisph(champ, sdd, matri)

    USE dimens_m, ONLY: iim

    REAL, intent(inout):: champ(:, :, :) ! (iim + 1, :, :)
    REAL, intent(in):: sdd(:) ! (iim)
    real, intent(in), dimension(:, :, :):: matri ! (iim, iim, :)

    ! Local:
    integer l, j

    !-----------------------------------------------------------------

    DO l = 1, size(champ, 3)
       DO j = 1, size(champ, 2)
          champ(:iim, j, l) = champ(:iim, j, l) * sdd
          champ(:iim, j, l) = (champ(:iim, j, l) &
               + matmul(matri(:, :, j), champ(:iim, j, l))) / sdd
          champ(iim + 1, j, l) = champ(1, j, l)
       END DO
    END DO

  end subroutine filtreg_hemisph

end module filtreg_hemisph_m
