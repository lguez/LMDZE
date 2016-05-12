module cv_driver_m

  implicit none

contains

  SUBROUTINE cv_driver(t1, q1, qs1, u1, v1, p1, ph1, iflag1, ft1, fq1, fu1, &
       fv1, precip1, VPrecip1, sig1, w01, icb1, inb1, delt, Ma1, upwd1, dnwd1, &
       dnwd01, qcondc1, cape1, da1, phi1, mp1)

    ! From LMDZ4/libf/phylmd/cv_driver.F, version 1.3, 2005/04/15 12:36:17
    ! Main driver for convection
    ! Author: S. Bony, March 2002

    ! Several modules corresponding to different physical processes

    use cv30_closure_m, only: cv30_closure
    use cv30_compress_m, only: cv30_compress
    use cv30_feed_m, only: cv30_feed
    use cv30_mixing_m, only: cv30_mixing
    use cv30_param_m, only: cv30_param, nl
    use cv30_prelim_m, only: cv30_prelim
    use cv30_tracer_m, only: cv30_tracer
    use cv30_trigger_m, only: cv30_trigger
    use cv30_uncompress_m, only: cv30_uncompress
    use cv30_undilute2_m, only: cv30_undilute2
    use cv30_unsat_m, only: cv30_unsat
    use cv30_yield_m, only: cv30_yield
    use cv_thermo_m, only: cv_thermo
    USE dimphy, ONLY: klev, klon

    real, intent(in):: t1(klon, klev) ! temperature (K)
    real, intent(in):: q1(klon, klev) ! specific humidity
    real, intent(in):: qs1(klon, klev) ! saturation specific humidity

    real, intent(in):: u1(klon, klev), v1(klon, klev)
    ! zonal wind and meridional velocity (m/s)

    real, intent(in):: p1(klon, klev) ! full level pressure (hPa)

    real, intent(in):: ph1(klon, klev + 1) 
    ! Half level pressure (hPa). These pressures are defined at levels
    ! intermediate between those of P1, T1, Q1 and QS1. The first
    ! value of PH should be greater than (i.e. at a lower level than)
    ! the first value of the array P1.

    integer, intent(out):: iflag1(klon)
    ! Flag for Emanuel conditions.

    ! 0: Moist convection occurs.

    ! 1: Moist convection occurs, but a CFL condition on the
    ! subsidence warming is violated. This does not cause the scheme
    ! to terminate.

    ! 2: Moist convection, but no precipitation because ep(inb) < 1e-4

    ! 3: No moist convection because new cbmf is 0 and old cbmf is 0.

    ! 4: No moist convection; atmosphere is not unstable

    ! 6: No moist convection because ihmin le minorig.

    ! 7: No moist convection because unreasonable parcel level
    ! temperature or specific humidity.

    ! 8: No moist convection: lifted condensation level is above the
    ! 200 mb level.

    ! 9: No moist convection: cloud base is higher then the level NL-1.

    real, intent(out):: ft1(klon, klev) ! temperature tendency (K/s)
    real, intent(out):: fq1(klon, klev) ! specific humidity tendency (s-1)

    real, intent(out):: fu1(klon, klev), fv1(klon, klev)
    ! forcing (tendency) of zonal and meridional velocity (m/s^2)

    real, intent(out):: precip1(klon) ! convective precipitation rate (mm/day)

    real, intent(out):: VPrecip1(klon, klev + 1)
    ! vertical profile of convective precipitation (kg/m2/s)

    real, intent(inout):: sig1(klon, klev) ! section of adiabatic updraft

    real, intent(inout):: w01(klon, klev) 
    ! vertical velocity within adiabatic updraft

    integer, intent(out):: icb1(klon)
    integer, intent(inout):: inb1(klon)
    real, intent(in):: delt ! the model time step (sec) between calls

    real, intent(out):: Ma1(klon, klev) ! mass flux of adiabatic updraft

    real, intent(out):: upwd1(klon, klev) 
    ! total upward mass flux (adiabatic + mixed)

    real, intent(out):: dnwd1(klon, klev) ! saturated downward mass flux (mixed)
    real, intent(out):: dnwd01(klon, klev) ! unsaturated downward mass flux

    real, intent(out):: qcondc1(klon, klev)
    ! in-cloud mixing ratio of condensed water

    real, intent(out):: cape1(klon)
    real, intent(inout):: da1(klon, klev), phi1(klon, klev, klev)
    real, intent(inout):: mp1(klon, klev)

    ! Local:

    real da(klon, klev), phi(klon, klev, klev), mp(klon, klev)
    integer i, k, il
    integer icbmax
    integer nk1(klon)
    integer icbs1(klon)
    real plcl1(klon)
    real tnk1(klon)
    real qnk1(klon)
    real gznk1(klon)
    real pbase1(klon)
    real buoybase1(klon)
    real lv1(klon, klev)
    real cpn1(klon, klev)
    real tv1(klon, klev)
    real gz1(klon, klev)
    real hm1(klon, klev)
    real h1(klon, klev)
    real tp1(klon, klev)
    real tvp1(klon, klev)
    real clw1(klon, klev)
    real th1(klon, klev)
    integer ncum

    ! Compressed fields:
    integer idcum(klon)
    integer iflag(klon), nk(klon), icb(klon)
    integer nent(klon, klev)
    integer icbs(klon)
    integer inb(klon)
    real plcl(klon), tnk(klon), qnk(klon), gznk(klon)
    real t(klon, klev), q(klon, klev), qs(klon, klev)
    real u(klon, klev), v(klon, klev)
    real gz(klon, klev), h(klon, klev), lv(klon, klev), cpn(klon, klev)
    real p(klon, klev), ph(klon, klev + 1), tv(klon, klev), tp(klon, klev)
    real clw(klon, klev)
    real pbase(klon), buoybase(klon), th(klon, klev)
    real tvp(klon, klev)
    real sig(klon, klev), w0(klon, klev)
    real hp(klon, klev), ep(klon, klev), sigp(klon, klev)
    real buoy(klon, klev)
    real cape(klon)
    real m(klon, klev), ment(klon, klev, klev), qent(klon, klev, klev)
    real uent(klon, klev, klev), vent(klon, klev, klev)
    real ments(klon, klev, klev), qents(klon, klev, klev)
    real sij(klon, klev, klev), elij(klon, klev, klev)
    real qp(klon, klev), up(klon, klev), vp(klon, klev)
    real wt(klon, klev), water(klon, klev), evap(klon, klev)
    real, allocatable:: b(:, :) ! (ncum, nl - 1)
    real ft(klon, klev), fq(klon, klev)
    real fu(klon, klev), fv(klon, klev)
    real upwd(klon, klev), dnwd(klon, klev), dnwd0(klon, klev)
    real Ma(klon, klev), mike(klon, klev), tls(klon, klev)
    real tps(klon, klev)
    real precip(klon)
    real VPrecip(klon, klev + 1)
    real qcondc(klon, klev) ! cld

    !-------------------------------------------------------------------

    ! SET CONSTANTS AND PARAMETERS

    ! set thermodynamical constants:
    CALL cv_thermo

    ! set convect parameters
    ! includes microphysical parameters and parameters that
    ! control the rate of approach to quasi-equilibrium)
    ! (common cvparam)
    CALL cv30_param(delt)

    ! INITIALIZE OUTPUT ARRAYS AND PARAMETERS

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
          dnwd01(i, k) = 0.
          qcondc1(i, k) = 0.
       end do
    end do

    do i = 1, klon
       precip1(i) = 0.
       iflag1(i) = 0
       cape1(i) = 0.
       VPrecip1(i, klev + 1) = 0.
    end do

    do il = 1, klon
       sig1(il, klev) = sig1(il, klev) + 1.
       sig1(il, klev) = min(sig1(il, klev), 12.1)
    enddo

    ! CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY & STATIC ENERGY
    CALL cv30_prelim(klon, klev, klev + 1, t1, q1, p1, ph1, lv1, cpn1, tv1, &
         gz1, h1, hm1, th1)

    ! CONVECTIVE FEED
    CALL cv30_feed(klon, klev, t1, q1, qs1, p1, ph1, gz1, nk1, icb1, &
         icbmax, iflag1, tnk1, qnk1, gznk1, plcl1) ! klev->na

    CALL cv30_undilute1(klon, klev, t1, q1, qs1, gz1, plcl1, p1, nk1, icb1, &
         tp1, tvp1, clw1, icbs1) ! klev->na

    ! TRIGGERING
    CALL cv30_trigger(klon, klev, icb1, plcl1, p1, th1, tv1, tvp1, pbase1, &
         buoybase1, iflag1, sig1, w01) ! klev->na

    ! Moist convective adjustment is necessary

    ncum = 0
    do i = 1, klon
       if (iflag1(i) == 0) then
          ncum = ncum + 1
          idcum(ncum) = i
       endif
    end do

    IF (ncum > 0) THEN
       allocate(b(ncum, nl - 1))
       CALL cv30_compress(ncum, iflag1, nk1, icb1, icbs1, plcl1, tnk1, qnk1, &
            gznk1, pbase1, buoybase1, t1, q1, qs1, u1, v1, gz1, th1, h1, lv1, &
            cpn1, p1, ph1, tv1, tp1, tvp1, clw1, sig1, w01, iflag, nk, icb, &
            icbs, plcl, tnk, qnk, gznk, pbase, buoybase, t, q, qs, u, v, gz, &
            th, h, lv, cpn, p, ph, tv, tp, tvp, clw, sig, w0)
       CALL cv30_undilute2(ncum, icb, icbs, nk, tnk, qnk, gznk, t, qs, gz, p, &
            h, tv, lv, pbase, buoybase, plcl, inb(:ncum), tp, tvp, clw, hp, &
            ep, sigp, buoy)

       ! CLOSURE
       CALL cv30_closure(klon, ncum, klev, icb, inb, pbase, p, ph, tv, &
            buoy, sig, w0, cape, m) ! na->klev

       ! MIXING
       CALL cv30_mixing(klon, ncum, klev, klev, icb, nk, inb, t, q, qs, u, &
            v, h, lv, hp, ep, clw, m, sig, ment, qent, uent, vent, nent, &
            sij, elij, ments, qents)

       ! Unsaturated (precipitating) downdrafts
       CALL cv30_unsat(icb(:ncum), inb(:ncum), t, q, qs, gz, u, v, p, ph, th, &
            tv, lv, cpn, ep(:ncum, :), sigp(:ncum, :), clw(:ncum, :), &
            m(:ncum, :), ment(:ncum, :, :), elij(:ncum, :, :), delt, plcl, &
            mp, qp(:ncum, :nl), up(:ncum, :nl), vp(:ncum, :nl), wt, water, &
            evap, b)

       ! Yield (tendencies, precipitation, variables of interface with
       ! other processes, etc)
       CALL cv30_yield(icb(:ncum), inb(:ncum), delt, t, q, u, v, gz, p, ph, &
            h, hp, lv, cpn, th, ep, clw, m, tp, mp, qp, up, vp, wt, &
            water(:ncum, :nl), evap(:ncum, :nl), b, ment, qent, uent, vent, &
            nent, elij, sig, tv, tvp, iflag, precip, VPrecip, ft, fq, fu, fv, &
            upwd, dnwd, dnwd0, ma, mike, tls, tps, qcondc)

       CALL cv30_tracer(klon, ncum, klev, ment, sij, da, phi)

       ! UNCOMPRESS THE FIELDS
       iflag1 = 42 ! for non convective points
       CALL cv30_uncompress(idcum(:ncum), iflag, precip, VPrecip, sig, w0, &
            ft, fq, fu, fv, inb, Ma, upwd, dnwd, dnwd0, qcondc, cape, &
            da, phi, mp, iflag1, precip1, VPrecip1, sig1, w01, ft1, fq1, &
            fu1, fv1, inb1, Ma1, upwd1, dnwd1, dnwd01, qcondc1, cape1, da1, &
            phi1, mp1)
    ENDIF

  end SUBROUTINE cv_driver

end module cv_driver_m
