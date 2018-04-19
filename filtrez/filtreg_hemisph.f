module filtreg_hemisph_m

  implicit none

contains

  subroutine filtreg_hemisph(champ, sdd, matri)

    USE dimensions, ONLY: iim

    REAL, intent(inout):: champ(:, :, :) ! (iim + 1, :, :)
    REAL, intent(in):: sdd(:) ! (iim) xprim[uv]^{\pm 1/2}

    real, intent(in):: matri(:, :, :) ! (iim, iim, :) 
    ! filtering matrix, last dimension is latitude

    ! Local:
    integer l, j

    !-----------------------------------------------------------------

    forall (j = 1:size(champ, 2), l = 1:size(champ, 3)) &
         champ(:iim, j, l) = champ(:iim, j, l) &
         + matmul(matri(:, :, j), champ(:iim, j, l) * sdd) / sdd
    champ(iim + 1, :, :) = champ(1, :, :)

  end subroutine filtreg_hemisph

end module filtreg_hemisph_m
