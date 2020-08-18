module coef_diff_turb_m

  implicit none

contains

  subroutine coef_diff_turb(nsrf, ni, paprs, pplay, u, v, q, t, ts, cdragm, &
       zgeop, coefm, coefh, q2)

    ! Computes coefficients for turbulent diffusion in the atmosphere.

    use nr_util, only: assert

    USE clesphys, ONLY: ok_kzmin
    use coefkz_m, only: coefkz
    use coefkzmin_m, only: coefkzmin
    use coefkz2_m, only: coefkz2
    USE conf_phys_m, ONLY: iflag_pbl
    USE dimphy, ONLY: klev
    use indicesol, only: is_oce
    USE suphec_m, ONLY: rd, rg, rkappa
    use ustarhb_m, only: ustarhb
    use yamada4_m, only: yamada4

    INTEGER, INTENT(IN):: nsrf
    INTEGER, INTENT(IN):: ni(:) ! (knon)
    REAL, INTENT(IN):: paprs(:, :) ! (knon, klev + 1)
    REAL, INTENT(IN):: pplay(:, :) ! (knon, klev)
    REAL, INTENT(IN):: u(:, :), v(:, :) ! (knon, klev) wind, in m s-1
    REAL, INTENT(IN):: q(:, :), t(:, :) ! (knon, klev)
    REAL, INTENT(IN):: ts(:), cdragm(:) ! (knon)
    REAL, INTENT(IN):: zgeop(:, :) ! (knon, klev) geopotential, in m2 s-2
    REAL, intent(out):: coefm(:, 2:) ! (knon, 2:klev) coefficient, vitesse

    real, intent(out):: coefh(:, 2:) ! (knon, 2:klev) 
    ! coefficient, chaleur et humidité

    real, intent(inout):: q2(:, :) ! (knon, klev + 1)

    ! Local:
    REAL coefm0(size(ni), 2:klev), coefh0(size(ni), 2:klev) ! (knon, 2:klev)
    REAL zlay(size(ni), klev), teta(size(ni), klev) ! (knon, klev)
    real zlev(size(ni), klev + 1) ! (knon, klev + 1)
    integer k

    !-------------------------------------------------------------------------

    call assert(size(ni) == [size(paprs, 1), size(pplay, 1), size(u, 1), &
         size(v, 1), size(q, 1), size(t, 1), size(ts), size(cdragm), &
         size(zgeop, 1), size(coefm, 1), size(coefh, 1), size(q2, 1)], &
         "coef_diff_turb knon")

    IF (iflag_pbl >= 6) THEN
       ! Mellor et Yamada adapt\'e \`a Mars, Richard Fournier et
       ! Fr\'ed\'eric Hourdin
       zlay(:, 1) = rd * t(:, 1) / (0.5 * (paprs(:, 1) &
            + pplay(:, 1))) * (paprs(:, 1) - pplay(:, 1)) / rg

       DO k = 2, klev
          zlay(:, k) = zlay(:, k-1) + rd * 0.5 * (t(:, k-1) + t(:, k)) &
               / paprs(:, k) * (pplay(:, k-1) - pplay(:, k)) / rg
       END DO

       forall (k = 1:klev) teta(:, k) = t(:, k) &
            * (paprs(:, 1) / pplay(:, k))**rkappa * (1. + 0.61 * q(:, k))

       zlev(:, 1) = 0.
       forall (k = 2:klev) zlev(:, k) = 0.5 * (zlay(:, k) + zlay(:, k-1))
       zlev(:, klev + 1) = 2. * zlay(:, klev) - zlev(:, klev)

       CALL yamada4(zlev, zlay, u, v, teta, q2, coefm, coefh, &
            ustarhb(u(:, 1), v(:, 1), cdragm))
    else
       CALL coefkz(nsrf, paprs, pplay, ts, u, v, t, q, zgeop, coefm, coefh)

       IF (iflag_pbl == 1 .and. nsrf == is_oce) THEN
          CALL coefkz2(paprs, pplay, t, coefm0, coefh0)
          coefm = max(coefm, coefm0)
          coefh = max(coefh, coefh0)
       END IF

       IF (ok_kzmin) THEN
          ! Calcul d'une diffusion minimale pour les conditions tres stables
          CALL coefkzmin(paprs, pplay, u, v, t, q, cdragm, coefh0)
          coefm = max(coefm, coefh0)
          coefh = max(coefh, coefh0)
       END IF
    END IF

  end subroutine coef_diff_turb

end module coef_diff_turb_m
