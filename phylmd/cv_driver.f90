module cv_driver_m

  implicit none

contains

  SUBROUTINE cv_driver(t1, q1, qs1, u1, v1, p1, ph1, iflag1, ft1, fq1, fu1, &
       fv1, rain, VPrecip1, sig1, w01, icb1, inb1, Ma1, upwd1, dnwd1, &
       qcondc1, cape1, da1, phi1, mp1)

    ! From LMDZ4/libf/phylmd/cv_driver.F, version 1.3, 2005/04/15 12:36:17
    ! Main driver for convection
    ! Author: S. Bony, March 2002

    ! Several modules corresponding to different physical processes

    use conf_gcm_m, only: dtphys
    use cv30_closure_m, only: cv30_closure
    use cv30_compress_m, only: cv30_compress
    use cv30_feed_m, only: cv30_feed
    use cv30_mixing_m, only: cv30_mixing
    use cv30_param_m, only: cv30_param, nl
    use cv30_prelim_m, only: cv30_prelim
    use cv30_tracer_m, only: cv30_tracer
    use cv30_trigger_m, only: cv30_trigger
    use cv30_uncompress_m, only: cv30_uncompress
    use cv30_undilute1_m, only: cv30_undilute1
    use cv30_undilute2_m, only: cv30_undilute2
    use cv30_unsat_m, only: cv30_unsat
    use cv30_yield_m, only: cv30_yield
    USE dimphy, ONLY: klev, klon

    real, intent(in):: t1(klon, klev) ! temperature, in K
    real, intent(in):: q1(klon, klev) ! specific humidity
    real, intent(in):: qs1(klon, klev) ! saturation specific humidity

    real, intent(in):: u1(klon, klev), v1(klon, klev)
    ! zonal wind and meridional velocity (m/s)

    real, intent(in):: p1(klon, klev) ! full level pressure, in hPa

    real, intent(in):: ph1(klon, klev + 1) 
    ! Half level pressure, in hPa. These pressures are defined at levels
    ! intermediate between those of P1, T1, Q1 and QS1. The first
    ! value of PH should be greater than (i.e. at a lower level than)
    ! the first value of the array P1.

    integer, intent(out):: iflag1(:) ! (klon)
    ! Flag for Emanuel conditions.

    ! 0: Moist convection occurs.

    ! 1: Moist convection occurs, but a CFL condition on the
    ! subsidence warming is violated. This does not cause the scheme
    ! to terminate.

    ! 2: Moist convection, but no precipitation because ep(inb) < 1e-4

    ! 3: No moist convection because new cbmf is 0 and old cbmf is 0.

    ! 4: No moist convection; atmosphere is not unstable.

    ! 6: No moist convection because ihmin <= minorig.

    ! 7: No moist convection because unreasonable parcel level
    ! temperature or specific humidity.

    ! 8: No moist convection: lifted condensation level is above the
    ! 200 mbar level.

    ! 9: No moist convection: cloud base is higher than the level NL-1.

    real, intent(out):: ft1(klon, klev) ! temperature tendency (K/s)
    real, intent(out):: fq1(klon, klev) ! specific humidity tendency (s-1)

    real, intent(out):: fu1(klon, klev), fv1(klon, klev)
    ! forcing (tendency) of zonal and meridional velocity (m/s^2)

    real, intent(out):: rain(klon) ! convective precipitation rate (mm/day)

    real, intent(out):: VPrecip1(klon, klev + 1)
    ! vertical profile of convective precipitation (kg/m2/s)

    real, intent(inout):: sig1(klon, klev) ! section of adiabatic updraft

    real, intent(inout):: w01(klon, klev) 
    ! vertical velocity within adiabatic updraft

    integer, intent(out):: icb1(klon)
    integer, intent(inout):: inb1(klon)
    real, intent(out):: Ma1(klon, klev) ! mass flux of adiabatic updraft

    real, intent(out):: upwd1(klon, klev) 
    ! total upward mass flux (adiabatic + mixed)

    real, intent(out):: dnwd1(klon, klev) ! saturated downward mass flux (mixed)

    real, intent(out):: qcondc1(klon, klev)
    ! in-cloud mixing ratio of condensed water

    real, intent(out):: cape1(klon)
    real, intent(out):: da1(:, :) ! (klon, klev)
    real, intent(out):: phi1(:, :, :) ! (klon, klev, klev)

    real, intent(out):: mp1(:, :) ! (klon, klev) Mass flux of the
    ! unsaturated downdraft, defined positive downward, in kg m-2
    ! s-1. M_p in Emanuel (1991 928).

    ! Local:

    real da(klon, klev), phi(klon, klev, klev)

    real, allocatable:: mp(:, :) ! (ncum, nl) Mass flux of the
    ! unsaturated downdraft, defined positive downward, in kg m-2
    ! s-1. M_p in Emanuel (1991 928).

    integer i, k, il
    integer icbs1(klon)
    real plcl1(klon)
    real tnk1(klon)
    real qnk1(klon)
    real gznk1(klon)
    real pbase1(klon)
    real buoybase1(klon)

    real lv1(klon, nl)
    ! specific latent heat of vaporization of water, in J kg-1

    real cpn1(klon, nl)
    ! specific heat capacity at constant pressure of humid air, in J K-1 kg-1

    real tv1(klon, klev)
    real gz1(klon, klev)
    real hm1(klon, klev)
    real h1(klon, klev)
    real tp1(klon, klev)
    real tvp1(klon, klev)
    real clw1(klon, klev)
    real th1(klon, nl) ! potential temperature, in K
    integer ncum

    ! Compressed fields:
    integer, allocatable:: idcum(:), iflag(:) ! (ncum)
    integer, allocatable:: icb(:) ! (ncum)
    integer, allocatable:: nent(:, :) ! (ncum, 2:nl - 1)
    integer icbs(klon)

    integer, allocatable:: inb(:) ! (ncum)
    ! first model level above the level of neutral buoyancy of the
    ! parcel (1 <= inb <= nl - 1)

    real, allocatable:: plcl(:) ! (ncum)
    real tnk(klon), qnk(klon), gznk(klon)
    real t(klon, klev), q(klon, klev), qs(klon, klev)
    real u(klon, klev), v(klon, klev)
    real gz(klon, klev), h(klon, klev)

    real, allocatable:: lv(:, :) ! (ncum, nl)
    ! specific latent heat of vaporization of water, in J kg-1

    real, allocatable:: cpn(:, :) ! (ncum, nl)
    ! specific heat capacity at constant pressure of humid air, in J K-1 kg-1

    real p(klon, klev) ! pressure at full level, in hPa
    real ph(klon, klev + 1), tv(klon, klev), tp(klon, klev)
    real clw(klon, klev)
    real pbase(klon), buoybase(klon)
    real, allocatable:: th(:, :) ! (ncum, nl)
    real tvp(klon, klev)
    real sig(klon, klev), w0(klon, klev)
    real hp(klon, klev), ep(klon, klev)
    real buoy(klon, klev)
    real cape(klon)
    real m(klon, klev), ment(klon, klev, klev), qent(klon, klev, klev)
    real uent(klon, klev, klev), vent(klon, klev, klev)
    real ments(klon, klev, klev), qents(klon, klev, klev)
    real sij(klon, klev, klev), elij(klon, klev, klev)
    real qp(klon, klev), up(klon, klev), vp(klon, klev)
    real wt(klon, klev), water(klon, klev)
    real, allocatable:: evap(:, :) ! (ncum, nl)
    real, allocatable:: b(:, :) ! (ncum, nl - 1)
    real ft(klon, klev), fq(klon, klev)
    real fu(klon, klev), fv(klon, klev)
    real upwd(klon, klev), dnwd(klon, klev)
    real Ma(klon, klev), mike(klon, klev), tls(klon, klev)
    real tps(klon, klev)
    real precip(klon)
    real VPrecip(klon, klev + 1)
    real qcondc(klon, klev) ! cloud

    !-------------------------------------------------------------------

    ! INITIALIZE OUTPUT ARRAYS AND PARAMETERS

    da1 = 0.
    mp1 = 0.
    phi1 = 0.

    do k = 1, klev
       do i = 1, klon
          ft1(i, k) = 0.
          fq1(i, k) = 0.
          fu1(i, k) = 0.
          fv1(i, k) = 0.
          tvp1(i, k) = 0.
          tp1(i, k) = 0.
          clw1(i, k) = 0.
          clw(i, k) = 0.
          gz1(i, k) = 0.
          VPrecip1(i, k) = 0.
          Ma1(i, k) = 0.
          upwd1(i, k) = 0.
          dnwd1(i, k) = 0.
          qcondc1(i, k) = 0.
       end do
    end do

    rain = 0.
    cape1 = 0.
    VPrecip1(:, klev + 1) = 0.

    do il = 1, klon
       sig1(il, klev) = sig1(il, klev) + 1.
       sig1(il, klev) = min(sig1(il, klev), 12.1)
    enddo

    CALL cv30_prelim(t1, q1, p1, ph1, lv1, cpn1, tv1, gz1, h1, hm1, th1)
    CALL cv30_feed(t1, q1, qs1, p1, ph1, gz1, icb1, iflag1, tnk1, qnk1, &
         gznk1, plcl1)
    CALL cv30_undilute1(t1, q1, qs1, gz1, plcl1, p1, icb1, tp1, tvp1, clw1, &
         icbs1)
    CALL cv30_trigger(icb1, plcl1, p1, th1, tv1, tvp1, pbase1, buoybase1, &
         iflag1, sig1, w01)

    ncum = count(iflag1 == 0)

    IF (ncum > 0) THEN
       ! Moist convective adjustment is necessary
       allocate(idcum(ncum), plcl(ncum), inb(ncum))
       allocate(b(ncum, nl - 1), evap(ncum, nl), icb(ncum), iflag(ncum))
       allocate(th(ncum, nl), lv(ncum, nl), cpn(ncum, nl), mp(ncum, nl))
       allocate(nent(ncum, 2:nl - 1))
       idcum = pack((/(i, i = 1, klon)/), iflag1 == 0)
       CALL cv30_compress(idcum, iflag1, icb1, icbs1, plcl1, tnk1, qnk1, &
            gznk1, pbase1, buoybase1, t1, q1, qs1, u1, v1, gz1, th1, h1, lv1, &
            cpn1, p1, ph1, tv1, tp1, tvp1, clw1, sig1, w01, icb, icbs, plcl, &
            tnk, qnk, gznk, pbase, buoybase, t, q, qs, u, v, gz, th, h, lv, &
            cpn, p, ph, tv, tp, tvp, clw, sig, w0)
       CALL cv30_undilute2(icb, icbs(:ncum), tnk, qnk, gznk, t, qs, gz, p, h, &
            tv, lv, pbase(:ncum), buoybase(:ncum), plcl, inb, tp, tvp, &
            clw, hp, ep, buoy)
       CALL cv30_closure(icb, inb, pbase, p, ph(:ncum, :), tv, buoy, sig, w0, &
            cape, m)
       CALL cv30_mixing(icb, inb, t, q, qs, u, v, h, lv, hp, ep, clw, m, sig, &
            ment, qent, uent, vent, nent, sij, elij, ments, qents)
       CALL cv30_unsat(icb, inb, t(:ncum, :nl), q(:ncum, :nl), &
            qs(:ncum, :nl), gz, u(:ncum, :nl), v(:ncum, :nl), p, &
            ph(:ncum, :), th(:ncum, :nl - 1), tv, lv, cpn, ep(:ncum, :), &
            clw(:ncum, :), m(:ncum, :), ment(:ncum, :, :), elij(:ncum, :, :), &
            dtphys, plcl, mp, qp(:ncum, :nl), up(:ncum, :nl), vp(:ncum, :nl), &
            wt(:ncum, :nl), water(:ncum, :nl), evap, b)
       CALL cv30_yield(icb, inb, dtphys, t, q, u, v, gz, p, ph, h, hp, &
            lv, cpn, th, ep, clw, m, tp, mp, qp, up, vp(:ncum, 2:nl), &
            wt(:ncum, :nl - 1), water(:ncum, :nl), evap, b, ment, qent, uent, &
            vent, nent, elij, sig, tv, tvp, iflag, precip, VPrecip, ft, fq, &
            fu, fv, upwd, dnwd, ma, mike, tls, tps, qcondc)
       CALL cv30_tracer(klon, ncum, klev, ment, sij, da, phi)
       CALL cv30_uncompress(idcum, iflag, precip, VPrecip, sig, w0, ft, fq, &
            fu, fv, inb, Ma, upwd, dnwd, qcondc, cape, da, phi, mp, iflag1, &
            rain, VPrecip1, sig1, w01, ft1, fq1, fu1, fv1, inb1, Ma1, &
            upwd1, dnwd1, qcondc1, cape1, da1, phi1, mp1)
    ENDIF

  end SUBROUTINE cv_driver

end module cv_driver_m
