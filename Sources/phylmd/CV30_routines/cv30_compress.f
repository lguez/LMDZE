module cv30_compress_m

  implicit none

contains

  SUBROUTINE cv30_compress(iflag1, nk1, icb1, icbs1, plcl1, tnk1, qnk1, &
       gznk1, pbase1, buoybase1, t1, q1, qs1, u1, v1, gz1, th1, h1, lv1, cpn1, &
       p1, ph1, tv1, tp1, tvp1, clw1, sig1, w01, nk, icb, icbs, plcl, tnk, &
       qnk, gznk, pbase, buoybase, t, q, qs, u, v, gz, th, h, lv, cpn, p, ph, &
       tv, tp, tvp, clw, sig, w0)

    ! Compress the fields (vectorization over convective gridpoints).

    use cv30_param_m, only: nl
    USE dimphy, ONLY: klev, klon
    use nr_util, only: assert

    ! inputs:
    integer, intent(in):: iflag1(:), nk1(:), icb1(:), icbs1(:) ! (klon)
    real, intent(in):: plcl1(klon), tnk1(klon), qnk1(klon), gznk1(klon)
    real pbase1(klon), buoybase1(klon)
    real, intent(in):: t1(klon, klev)
    real, intent(in):: q1(klon, klev), qs1(klon, klev)
    real, intent(in):: u1(klon, klev), v1(klon, klev)
    real gz1(klon, klev), h1(klon, klev), lv1(klon, klev), cpn1(klon, klev)
    real, intent(in):: p1(klon, klev), ph1(klon, klev + 1)
    real, intent(in):: tv1(klon, klev), tp1(klon, klev)
    real tvp1(klon, klev), clw1(klon, klev)
    real th1(klon, klev)
    real sig1(klon, klev), w01(klon, klev)

    ! outputs:
    integer nk(:) ! (klon)
    integer, intent(out):: icb(:) ! (ncum) {2 <= icb <= nl - 3}
    integer icbs(klon)
    real, intent(out):: plcl(:) ! (ncum)
    real tnk(:), qnk(:), gznk(:) ! (klon)
    real pbase(klon), buoybase(klon)
    real t(klon, klev), q(klon, klev), qs(klon, klev)
    real u(klon, klev), v(klon, klev)
    real gz(klon, klev), h(klon, klev), lv(klon, klev), cpn(klon, klev)
    real p(klon, klev)
    real ph(:, :) ! (klon, klev + 1)
    real tv(klon, klev), tp(klon, klev)
    real tvp(klon, klev), clw(klon, klev)
    real th(klon, klev)
    real sig(klon, klev), w0(klon, klev)

    ! Local:
    integer i, k, nn, ncum

    !---------------------------------------------------------------

    ncum = size(icb)

    do k = 1, nl + 1
       nn = 0
       do i = 1, klon
          if (iflag1(i) == 0) then
             nn = nn + 1
             sig(nn, k) = sig1(i, k)
             w0(nn, k) = w01(i, k)
             t(nn, k) = t1(i, k)
             q(nn, k) = q1(i, k)
             qs(nn, k) = qs1(i, k)
             u(nn, k) = u1(i, k)
             v(nn, k) = v1(i, k)
             gz(nn, k) = gz1(i, k)
             h(nn, k) = h1(i, k)
             lv(nn, k) = lv1(i, k)
             cpn(nn, k) = cpn1(i, k)
             p(nn, k) = p1(i, k)
             ph(nn, k) = ph1(i, k)
             tv(nn, k) = tv1(i, k)
             tp(nn, k) = tp1(i, k)
             tvp(nn, k) = tvp1(i, k)
             clw(nn, k) = clw1(i, k)
             th(nn, k) = th1(i, k)
          endif
       end do
    end do

    nn = 0

    do i = 1, klon
       if (iflag1(i) == 0) then
          nn = nn + 1
          pbase(nn) = pbase1(i)
          buoybase(nn) = buoybase1(i)
          plcl(nn) = plcl1(i)
          tnk(nn) = tnk1(i)
          qnk(nn) = qnk1(i)
          gznk(nn) = gznk1(i)
          nk(nn) = nk1(i)
          icb(nn) = icb1(i)
          icbs(nn) = icbs1(i)
       endif
    end do

    do i = 1, ncum
       call assert(2 <= icb(i) .and. icb(i) <= nl - 3 .and. ph(i, icb(i) + 1) &
            < plcl(i) .and. (plcl(i) <= ph(i, icb(i)) .or. icb(i) == 2), &
            "cv30_compress")
    end do

  end SUBROUTINE cv30_compress

end module cv30_compress_m
