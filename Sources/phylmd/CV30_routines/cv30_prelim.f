module cv30_prelim_m

  implicit none

contains

  SUBROUTINE cv30_prelim(t1, q1, p1, ph1, lv1, cpn1, tv1, gz1, h1, hm1, th1)

    USE cv30_param_m, ONLY: nl
    USE cv_thermo_m, ONLY: cl, clmcpv, cpd, cpv, eps, rrd, rrv
    USE dimphy, ONLY: klev, klon
    use SUPHEC_M, only: rlvtt

    ! Calculate arrays of geopotential, heat capacity and static energy

    real, intent(in):: t1(klon, klev)
    real, intent(in):: q1(klon, klev)
    real, intent(in):: p1(klon, klev), ph1(klon, klev + 1)

    ! outputs:
    real lv1(klon, klev), cpn1(klon, klev), tv1(klon, klev)
    real gz1(klon, klev), h1(klon, klev), hm1(klon, klev)
    real th1(klon, klev) ! potential temperature

    ! Local:
    integer k, i
    real rdcp
    real tvx, tvy 
    real cpx(klon, klev)

    !--------------------------------------------------------------

    do k = 1, nl 
       do i = 1, klon
          lv1(i, k) =  rlvtt - clmcpv * (t1(i, k) - 273.15)
          cpn1(i, k) = cpd * (1. - q1(i, k)) + cpv * q1(i, k)
          cpx(i, k) = cpd * (1. - q1(i, k)) + cl * q1(i, k)
          tv1(i, k) = t1(i, k) * (1. + q1(i, k)/eps - q1(i, k))
          rdcp = (rrd * (1. - q1(i, k)) + q1(i, k) * rrv)/cpn1(i, k)
          th1(i, k) = t1(i, k) * (1000./p1(i, k))**rdcp
       end do
    end do

    ! gz1 = phi at the full levels (same as p1).

    do i = 1, klon
       gz1(i, 1) = 0.
    end do

    do k = 2, nl 
       do i = 1, klon
          tvx = t1(i, k) * (1. + q1(i, k)/eps - q1(i, k))
          tvy = t1(i, k - 1) * (1. + q1(i, k - 1)/eps - q1(i, k - 1))
          gz1(i, k) = gz1(i, k - 1) + 0.5 * rrd * (tvx + tvy) &
               * (p1(i, k - 1) - p1(i, k))/ph1(i, k)
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
