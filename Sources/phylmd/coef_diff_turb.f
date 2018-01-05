module coef_diff_turb_m

  implicit none

contains

  subroutine coef_diff_turb(dtime, nsrf, ni, ypaprs, ypplay, yu, yv, yq, yt, &
       yts, ycdragm, zgeop, ycoefm, ycoefh, yq2)

    USE clesphys, ONLY: ok_kzmin
    use coefkz_m, only: coefkz
    use coefkzmin_m, only: coefkzmin
    use coefkz2_m, only: coefkz2
    USE conf_phys_m, ONLY: iflag_pbl
    USE dimphy, ONLY: klev, klon
    USE suphec_m, ONLY: rd, rg, rkappa
    use ustarhb_m, only: ustarhb
    use yamada4_m, only: yamada4

    REAL, INTENT(IN):: dtime ! interval du temps (secondes)
    INTEGER, INTENT(IN):: nsrf
    INTEGER, INTENT(IN):: ni(:) ! (knon)
    REAL, INTENT(IN):: ypaprs(:, :) ! (klon, klev + 1)
    REAL, INTENT(IN):: ypplay(:, :) ! (klon, klev)
    REAL, INTENT(IN), dimension(:, :):: yu, yv, yq, yt ! (klon, klev)
    REAL, INTENT(IN):: yts(:), ycdragm(:) ! (knon)
    REAL, INTENT(IN):: zgeop(:, :) ! (knon, klev)
    REAL, intent(out):: ycoefm(:, 2:) ! (knon, 2:klev) coefficient, vitesse

    real, intent(out):: ycoefh(:, 2:) ! (knon, 2:klev) 
    ! coefficient, chaleur et humidité

    real, intent(inout):: yq2(:, :) ! (klon, klev + 1)

    ! Local:
    REAL ycoefm0(klon, 2:klev), ycoefh0(klon, 2:klev)
    REAL yzlay(klon, klev), zlev(klon, klev + 1), yteta(klon, klev)
    REAL ustar(klon)
    integer k, knon

    !-------------------------------------------------------------------------

    knon = size(ni)
    CALL coefkz(nsrf, ypaprs(:knon, :), ypplay(:knon, :), yts(:knon), &
         yu(:knon, :), yv(:knon, :), yt(:knon, :), yq(:knon, :), zgeop, &
         ycoefm, ycoefh)

    IF (iflag_pbl == 1) THEN
       CALL coefkz2(nsrf, knon, ypaprs, ypplay, yt, ycoefm0(:knon, :), &
            ycoefh0(:knon, :))
       ycoefm = max(ycoefm, ycoefm0(:knon, :))
       ycoefh = max(ycoefh, ycoefh0(:knon, :))
    END IF

    IF (ok_kzmin) THEN
       ! Calcul d'une diffusion minimale pour les conditions tres stables
       CALL coefkzmin(knon, ypaprs, ypplay, yu, yv, yt, yq, ycdragm(:knon), &
            ycoefh0(:knon, :))
       ycoefm0(:knon, :) = ycoefh0(:knon, :)
       ycoefm = max(ycoefm, ycoefm0(:knon, :))
       ycoefh = max(ycoefh, ycoefh0(:knon, :))
    END IF

    IF (iflag_pbl >= 6) THEN
       ! Mellor et Yamada adapt\'e \`a Mars, Richard Fournier et
       ! Fr\'ed\'eric Hourdin
       yzlay(:knon, 1) = rd * yt(:knon, 1) / (0.5 * (ypaprs(:knon, 1) &
            + ypplay(:knon, 1))) &
            * (ypaprs(:knon, 1) - ypplay(:knon, 1)) / rg

       DO k = 2, klev
          yzlay(:knon, k) = yzlay(:knon, k-1) &
               + rd * 0.5 * (yt(1:knon, k-1) + yt(1:knon, k)) &
               / ypaprs(1:knon, k) &
               * (ypplay(1:knon, k-1) - ypplay(1:knon, k)) / rg
       END DO

       DO k = 1, klev
          yteta(1:knon, k) = yt(1:knon, k) * (ypaprs(1:knon, 1) &
               / ypplay(1:knon, k))**rkappa * (1. + 0.61 * yq(1:knon, k))
       END DO

       zlev(:knon, 1) = 0.
       zlev(:knon, klev + 1) = 2. * yzlay(:knon, klev) &
            - yzlay(:knon, klev - 1)

       DO k = 2, klev
          zlev(:knon, k) = 0.5 * (yzlay(:knon, k) + yzlay(:knon, k-1))
       END DO

       ustar(:knon) = ustarhb(yu(:knon, 1), yv(:knon, 1), ycdragm(:knon))
       CALL yamada4(dtime, rg, zlev(:knon, :), yzlay(:knon, :), &
            yu(:knon, :), yv(:knon, :), yteta(:knon, :), yq2(:knon, :), &
            ycoefm, ycoefh, ustar(:knon))
    END IF

  end subroutine coef_diff_turb

end module coef_diff_turb_m
