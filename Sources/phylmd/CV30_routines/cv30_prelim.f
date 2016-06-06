module cv30_prelim_m

  implicit none

contains

  SUBROUTINE cv30_prelim(t1, q1, p1, ph1, lv1, cpn1, tv1, gz1, h1, hm1, th1)

    USE cv30_param_m, ONLY: nl
    USE cv_thermo_m, ONLY: clmcpv, eps
    USE dimphy, ONLY: klev, klon
    use SUPHEC_M, only: rcw, rlvtt, rcpd, rcpv, rd, rv

    ! Calculate arrays of geopotential, heat capacity and static energy

    real, intent(in):: t1(:, :) ! (klon, klev) temperature, in K
    real, intent(in):: q1(:, :) ! (klon, klev) specific humidity
    real, intent(in):: p1(:, :) ! (klon, klev) full level pressure, in hPa
    real, intent(in):: ph1(:, :) ! (klon, klev + 1) half level pressure, in hPa

    ! outputs:

    real, intent(out):: lv1(:, :) ! (klon, nl)
    ! specific latent heat of vaporization of water, in J kg-1
    
    real, intent(out):: cpn1(:, :) ! (klon, nl)
    ! specific heat capacity at constant pressure of humid air, in J K-1 kg-1

    real tv1(:, :) ! (klon, klev)
    real gz1(klon, klev), h1(klon, klev), hm1(klon, klev)
    real, intent(out):: th1(:, :) ! (klon, nl) potential temperature, in K

    ! Local:
    integer k, i
    real kappa
    real tvx, tvy 
    real cpx(klon, klev)

    !--------------------------------------------------------------

    do k = 1, nl 
       do i = 1, klon
          lv1(i, k) =  rlvtt - clmcpv * (t1(i, k) - 273.15)
          cpn1(i, k) = rcpd * (1. - q1(i, k)) + rcpv * q1(i, k)
          cpx(i, k) = rcpd * (1. - q1(i, k)) + rcw * q1(i, k)
          tv1(i, k) = t1(i, k) * (1. + q1(i, k) / eps - q1(i, k))
          kappa = (rd * (1. - q1(i, k)) + q1(i, k) * rv) / cpn1(i, k)
          th1(i, k) = t1(i, k) * (1000. / p1(i, k))**kappa
       end do
    end do

    ! gz1 = phi at the full levels (same as p1).

    do i = 1, klon
       gz1(i, 1) = 0.
    end do

    do k = 2, nl 
       do i = 1, klon
          tvx = t1(i, k) * (1. + q1(i, k) / eps - q1(i, k))
          tvy = t1(i, k - 1) * (1. + q1(i, k - 1) / eps - q1(i, k - 1))
          gz1(i, k) = gz1(i, k - 1) + 0.5 * rd * (tvx + tvy) &
               * (p1(i, k - 1) - p1(i, k)) / ph1(i, k)
       end do
    end do

    ! h1 = phi + cpT (dry static energy).
    ! hm1 = phi + cp(T1 - Tbase) + Lq

    do k = 1, nl 
       do i = 1, klon
          h1(i, k) = gz1(i, k) + cpn1(i, k) * t1(i, k)
          hm1(i, k) = gz1(i, k) + cpx(i, k) * (t1(i, k) - t1(i, 1)) &
               + lv1(i, k) * q1(i, k)
       end do
    end do

  end SUBROUTINE cv30_prelim

end module cv30_prelim_m
