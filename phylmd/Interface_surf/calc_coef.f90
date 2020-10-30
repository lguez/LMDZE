module calc_coef_m

  implicit none

contains

  subroutine calc_coef(c, d, a, b, x, gamma, delp, zx_coef)

    USE dimphy, ONLY: klev
    USE suphec_m, ONLY: rg

    REAL, intent(out):: c(:, 2:), d(:, 2:) ! (knon, 2:klev)
    REAL, intent(out):: a(:), b(:) ! (knon)
    REAL, intent(in):: x(:, :) ! (knon, klev)
    REAL, intent(in):: gamma(:, 2:) ! (knon, 2:klev)

    REAL, intent(in):: delp(:, :) ! (knon, klev)
    ! epaisseur de couche en pression (Pa)

    REAL, intent(in):: zx_coef(:, 2:) ! (knon, 2:klev)

    ! Local:
    
    REAL buf(size(x, 1)) ! (knon)
    INTEGER k

    !-------------------------------------------------------------------

    buf = delp(:, klev) + zx_coef(:, klev)
    c(:, klev) = (x(:, klev) * delp(:, klev) + zx_coef(:, klev) &
         * gamma(:, klev)) / buf
    d(:, klev) = zx_coef(:, klev) / buf

    DO k = klev - 1, 2, - 1
       buf = delp(:, k) + zx_coef(:, k) + zx_coef(:, k + 1) * (1. - d(:, k + 1))
       c(:, k) = (x(:, k) * delp(:, k) + zx_coef(:, k + 1) * c(:, k + 1) &
            - zx_coef(:, k + 1) * gamma(:, k + 1) + zx_coef(:, k) &
            * gamma(:, k)) / buf
       d(:, k) = zx_coef(:, k) / buf
    ENDDO

    buf = delp(:, 1) + zx_coef(:, 2) * (1. - d(:, 2))
    a = (x(:, 1) * delp(:, 1) + zx_coef(:, 2) * (c(:, 2) - gamma(:, 2))) / buf
    b = - RG / buf

  end subroutine calc_coef

end module calc_coef_m
