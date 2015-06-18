module cv_driver_m

  implicit none

contains

  SUBROUTINE cv_driver(t1, q1, qs1, u1, v1, p1, ph1, iflag1, ft1, &
       fq1, fu1, fv1, precip1, VPrecip1, cbmf1, sig1, w01, icb1, inb1, delt, &
       Ma1, upwd1, dnwd1, dnwd01, qcondc1, wd1, cape1, da1, phi1, mp1)

    ! From LMDZ4/libf/phylmd/cv_driver.F, version 1.3, 2005/04/15 12:36:17
    ! Main driver for convection
    ! Author: S. Bony, March 2002

    ! Several modules corresponding to different physical processes

    ! Several versions of convect may be used:
    ! - iflag_con = 3: version lmd
    ! - iflag_con = 4: version 4.3b

    use clesphys2, only: iflag_con
    use cv3_compress_m, only: cv3_compress
    use cv3_feed_m, only: cv3_feed
    use cv3_mixing_m, only: cv3_mixing
    use cv3_param_m, only: cv3_param
    use cv3_prelim_m, only: cv3_prelim
    use cv3_tracer_m, only: cv3_tracer
    use cv3_uncompress_m, only: cv3_uncompress
    use cv3_unsat_m, only: cv3_unsat
    use cv3_yield_m, only: cv3_yield
    use cv_feed_m, only: cv_feed
    use cv_uncompress_m, only: cv_uncompress
    USE dimphy, ONLY: klev, klon

    real, intent(in):: t1(klon, klev) ! temperature
    real, intent(in):: q1(klon, klev) ! specific hum
    real, intent(in):: qs1(klon, klev) ! sat specific hum
    real, intent(in):: u1(klon, klev) ! u-wind
    real, intent(in):: v1(klon, klev) ! v-wind
    real, intent(in):: p1(klon, klev) ! full level pressure
    real, intent(in):: ph1(klon, klev + 1) ! half level pressure
    integer, intent(out):: iflag1(klon) ! flag for Emanuel conditions
    real, intent(out):: ft1(klon, klev) ! temp tend
    real, intent(out):: fq1(klon, klev) ! spec hum tend
    real, intent(out):: fu1(klon, klev) ! u-wind tend
    real, intent(out):: fv1(klon, klev) ! v-wind tend
    real, intent(out):: precip1(klon) ! precipitation

    real, intent(out):: VPrecip1(klon, klev+1)
    ! vertical profile of precipitation

    real, intent(inout):: cbmf1(klon) ! cloud base mass flux
    real, intent(inout):: sig1(klon, klev) ! section adiabatic updraft

    real, intent(inout):: w01(klon, klev) 
    ! vertical velocity within adiabatic updraft

    integer, intent(out):: icb1(klon)
    integer, intent(inout):: inb1(klon)
    real, intent(in):: delt ! time step
    real Ma1(klon, klev)
    ! Ma1 Real Output mass flux adiabatic updraft
    real, intent(out):: upwd1(klon, klev) ! total upward mass flux (adiab+mixed)
    real, intent(out):: dnwd1(klon, klev) ! saturated downward mass flux (mixed)
    real, intent(out):: dnwd01(klon, klev) ! unsaturated downward mass flux

    real qcondc1(klon, klev) ! cld
    ! qcondc1 Real Output in-cld mixing ratio of condensed water
    real wd1(klon) ! gust
    ! wd1 Real Output downdraft velocity scale for sfc fluxes
    real cape1(klon)
    ! cape1 Real Output CAPE

    real, intent(inout):: da1(klon, klev), phi1(klon, klev, klev)
    real, intent(inout):: mp1(klon, klev)

    ! --- ARGUMENTS

    ! --- On input:

    ! t: Array of absolute temperature (K) of dimension KLEV, with first
    ! index corresponding to lowest model level. Note that this array
    ! will be altered by the subroutine if dry convective adjustment
    ! occurs and if IPBL is not equal to 0.

    ! q: Array of specific humidity (gm/gm) of dimension KLEV, with first
    ! index corresponding to lowest model level. Must be defined
    ! at same grid levels as T. Note that this array will be altered
    ! if dry convective adjustment occurs and if IPBL is not equal to 0.

    ! qs: Array of saturation specific humidity of dimension KLEV, with first
    ! index corresponding to lowest model level. Must be defined
    ! at same grid levels as T. Note that this array will be altered
    ! if dry convective adjustment occurs and if IPBL is not equal to 0.

    ! u: Array of zonal wind velocity (m/s) of dimension KLEV, witth first
    ! index corresponding with the lowest model level. Defined at
    ! same levels as T. Note that this array will be altered if
    ! dry convective adjustment occurs and if IPBL is not equal to 0.

    ! v: Same as u but for meridional velocity.

    ! p: Array of pressure (mb) of dimension KLEV, with first
    ! index corresponding to lowest model level. Must be defined
    ! at same grid levels as T.

    ! ph: Array of pressure (mb) of dimension KLEV+1, with first index
    ! corresponding to lowest level. These pressures are defined at
    ! levels intermediate between those of P, T, Q and QS. The first
    ! value of PH should be greater than (i.e. at a lower level than)
    ! the first value of the array P.

    ! nl: The maximum number of levels to which convection can penetrate, plus 1
    ! NL MUST be less than or equal to KLEV-1.

    ! delt: The model time step (sec) between calls to CONVECT

    ! --- On Output:

    ! iflag: An output integer whose value denotes the following:
    ! VALUE INTERPRETATION
    ! ----- --------------
    ! 0 Moist convection occurs.
    ! 1 Moist convection occurs, but a CFL condition
    ! on the subsidence warming is violated. This
    ! does not cause the scheme to terminate.
    ! 2 Moist convection, but no precip because ep(inb) lt 0.0001
    ! 3 No moist convection because new cbmf is 0 and old cbmf is 0.
    ! 4 No moist convection; atmosphere is not
    ! unstable
    ! 6 No moist convection because ihmin le minorig.
    ! 7 No moist convection because unreasonable
    ! parcel level temperature or specific humidity.
    ! 8 No moist convection: lifted condensation
    ! level is above the 200 mb level.
    ! 9 No moist convection: cloud base is higher
    ! then the level NL-1.

    ! ft: Array of temperature tendency (K/s) of dimension KLEV, defined at same
    ! grid levels as T, Q, QS and P.

    ! fq: Array of specific humidity tendencies ((gm/gm)/s) of dimension KLEV,
    ! defined at same grid levels as T, Q, QS and P.

    ! fu: Array of forcing of zonal velocity (m/s^2) of dimension KLEV,
    ! defined at same grid levels as T.

    ! fv: Same as FU, but for forcing of meridional velocity.

    ! precip: Scalar convective precipitation rate (mm/day).

    ! VPrecip: Vertical profile of convective precipitation (kg/m2/s).

    ! wd: A convective downdraft velocity scale. For use in surface
    ! flux parameterizations. See convect.ps file for details.

    ! tprime: A convective downdraft temperature perturbation scale (K).
    ! For use in surface flux parameterizations. See convect.ps
    ! file for details.

    ! qprime: A convective downdraft specific humidity
    ! perturbation scale (gm/gm).
    ! For use in surface flux parameterizations. See convect.ps
    ! file for details.

    ! cbmf: The cloud base mass flux ((kg/m**2)/s). THIS SCALAR VALUE MUST
    ! BE STORED BY THE CALLING PROGRAM AND RETURNED TO CONVECT AT
    ! ITS NEXT CALL. That is, the value of CBMF must be "remembered"
    ! by the calling program between calls to CONVECT.

    ! det: Array of detrainment mass flux of dimension KLEV.

    ! Local arrays

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

    ! (local) compressed fields:

    integer idcum(klon)
    integer iflag(klon), nk(klon), icb(klon)
    integer nent(klon, klev)
    integer icbs(klon)
    integer inb(klon), inbis(klon)

    real cbmf(klon), plcl(klon), tnk(klon), qnk(klon), gznk(klon)
    real t(klon, klev), q(klon, klev), qs(klon, klev)
    real u(klon, klev), v(klon, klev)
    real gz(klon, klev), h(klon, klev), lv(klon, klev), cpn(klon, klev)
    real p(klon, klev), ph(klon, klev+1), tv(klon, klev), tp(klon, klev)
    real clw(klon, klev)
    real dph(klon, klev)
    real pbase(klon), buoybase(klon), th(klon, klev)
    real tvp(klon, klev)
    real sig(klon, klev), w0(klon, klev)
    real hp(klon, klev), ep(klon, klev), sigp(klon, klev)
    real frac(klon), buoy(klon, klev)
    real cape(klon)
    real m(klon, klev), ment(klon, klev, klev), qent(klon, klev, klev)
    real uent(klon, klev, klev), vent(klon, klev, klev)
    real ments(klon, klev, klev), qents(klon, klev, klev)
    real sij(klon, klev, klev), elij(klon, klev, klev)
    real qp(klon, klev), up(klon, klev), vp(klon, klev)
    real wt(klon, klev), water(klon, klev), evap(klon, klev)
    real b(klon, klev), ft(klon, klev), fq(klon, klev)
    real fu(klon, klev), fv(klon, klev)
    real upwd(klon, klev), dnwd(klon, klev), dnwd0(klon, klev)
    real Ma(klon, klev), mike(klon, klev), tls(klon, klev)
    real tps(klon, klev), qprime(klon), tprime(klon)
    real precip(klon)
    real VPrecip(klon, klev+1)
    real qcondc(klon, klev) ! cld
    real wd(klon) ! gust

    !-------------------------------------------------------------------
    ! --- SET CONSTANTS AND PARAMETERS

    ! -- set simulation flags:
    ! (common cvflag)

    CALL cv_flag

    ! -- set thermodynamical constants:
    ! (common cvthermo)

    CALL cv_thermo

    ! -- set convect parameters

    ! includes microphysical parameters and parameters that
    ! control the rate of approach to quasi-equilibrium)
    ! (common cvparam)

    if (iflag_con == 3) CALL cv3_param(klev, delt)

    ! --- INITIALIZE OUTPUT ARRAYS AND PARAMETERS

    do k = 1, klev
       do i = 1, klon
          ft1(i, k) = 0.0
          fq1(i, k) = 0.0
          fu1(i, k) = 0.0
          fv1(i, k) = 0.0
          tvp1(i, k) = 0.0
          tp1(i, k) = 0.0
          clw1(i, k) = 0.0
          !ym
          clw(i, k) = 0.0
          gz1(i, k) = 0.
          VPrecip1(i, k) = 0.
          Ma1(i, k) = 0.0
          upwd1(i, k) = 0.0
          dnwd1(i, k) = 0.0
          dnwd01(i, k) = 0.0
          qcondc1(i, k) = 0.0
       end do
    end do

    do i = 1, klon
       precip1(i) = 0.0
       iflag1(i) = 0
       wd1(i) = 0.0
       cape1(i) = 0.0
       VPrecip1(i, klev+1) = 0.0
    end do

    if (iflag_con == 3) then
       do il = 1, klon
          sig1(il, klev) = sig1(il, klev) + 1.
          sig1(il, klev) = min(sig1(il, klev), 12.1)
       enddo
    endif

    ! --- CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY & STATIC ENERGY

    if (iflag_con == 3) then
       CALL cv3_prelim(klon, klev, klev + 1, t1, q1, p1, ph1, lv1, cpn1, tv1, &
            gz1, h1, hm1, th1)
    else
       ! iflag_con == 4
       CALL cv_prelim(klon, klev, klev + 1, t1, q1, p1, ph1, lv1, cpn1, tv1, &
            gz1, h1, hm1)
    endif

    ! --- CONVECTIVE FEED

    if (iflag_con == 3) then
       CALL cv3_feed(klon, klev, t1, q1, qs1, p1, ph1, gz1, nk1, icb1, &
            icbmax, iflag1, tnk1, qnk1, gznk1, plcl1) ! klev->na
    else
       ! iflag_con == 4
       CALL cv_feed(klon, klev, t1, q1, qs1, p1, hm1, gz1, nk1, icb1, icbmax, &
            iflag1, tnk1, qnk1, gznk1, plcl1)
    endif

    ! --- UNDILUTE (ADIABATIC) UPDRAFT / 1st part
    ! (up through ICB for convect4, up through ICB+1 for convect3)
    ! Calculates the lifted parcel virtual temperature at nk, the
    ! actual temperature, and the adiabatic liquid water content.

    if (iflag_con == 3) then
       CALL cv3_undilute1(klon, klev, t1, q1, qs1, gz1, plcl1, p1, nk1, icb1, &
            tp1, tvp1, clw1, icbs1) ! klev->na
    else
       ! iflag_con == 4
       CALL cv_undilute1(klon, klev, t1, q1, qs1, gz1, p1, nk1, icb1, icbmax, &
            tp1, tvp1, clw1)
    endif

    ! --- TRIGGERING

    if (iflag_con == 3) then
       CALL cv3_trigger(klon, klev, icb1, plcl1, p1, th1, tv1, tvp1, pbase1, &
            buoybase1, iflag1, sig1, w01) ! klev->na
    else
       ! iflag_con == 4
       CALL cv_trigger(klon, klev, icb1, cbmf1, tv1, tvp1, iflag1)
    end if

    ! --- IF THIS POINT IS REACHED, MOIST CONVECTIVE ADJUSTMENT IS NECESSARY

    ncum = 0
    do i = 1, klon
       if(iflag1(i) == 0)then
          ncum = ncum+1
          idcum(ncum) = i
       endif
    end do

    IF (ncum > 0) THEN
       ! --- COMPRESS THE FIELDS
       ! (-> vectorization over convective gridpoints)

       if (iflag_con == 3) then
          CALL cv3_compress(klon, klon, ncum, klev, iflag1, nk1, icb1, icbs1, &
               plcl1, tnk1, qnk1, gznk1, pbase1, buoybase1, t1, q1, qs1, u1, &
               v1, gz1, th1, h1, lv1, cpn1, p1, ph1, tv1, tp1, tvp1, clw1, &
               sig1, w01, iflag, nk, icb, icbs, plcl, tnk, qnk, gznk, pbase, &
               buoybase, t, q, qs, u, v, gz, th, h, lv, cpn, p, ph, tv, tp, &
               tvp, clw, sig, w0)
       else
          ! iflag_con == 4
          CALL cv_compress(klon, klon, ncum, klev, iflag1, nk1, icb1, cbmf1, &
               plcl1, tnk1, qnk1, gznk1, t1, q1, qs1, u1, v1, gz1, h1, lv1, &
               cpn1, p1, ph1, tv1, tp1, tvp1, clw1, iflag, nk, icb, cbmf, &
               plcl, tnk, qnk, gznk, t, q, qs, u, v, gz, h, lv, cpn, p, ph, &
               tv, tp, tvp, clw, dph)
       endif

       ! --- UNDILUTE (ADIABATIC) UPDRAFT / second part :
       ! --- FIND THE REST OF THE LIFTED PARCEL TEMPERATURES
       ! --- &
       ! --- COMPUTE THE PRECIPITATION EFFICIENCIES AND THE
       ! --- FRACTION OF PRECIPITATION FALLING OUTSIDE OF CLOUD
       ! --- &
       ! --- FIND THE LEVEL OF NEUTRAL BUOYANCY

       if (iflag_con == 3) then
          CALL cv3_undilute2(klon, ncum, klev, icb, icbs, nk, tnk, qnk, gznk, &
               t, qs, gz, p, h, tv, lv, pbase, buoybase, plcl, inb, tp, &
               tvp, clw, hp, ep, sigp, buoy) !na->klev
       else
          ! iflag_con == 4
          CALL cv_undilute2(klon, ncum, klev, icb, nk, tnk, qnk, gznk, t, &
               qs, gz, p, dph, h, tv, lv, inb, inbis, tp, tvp, clw, hp, ep, &
               sigp, frac)
       endif

       ! --- CLOSURE

       if (iflag_con == 3) then
          CALL cv3_closure(klon, ncum, klev, icb, inb, pbase, p, ph, tv, &
               buoy, sig, w0, cape, m) ! na->klev
       else
          ! iflag_con == 4
          CALL cv_closure(klon, ncum, klev, nk, icb, tv, tvp, p, ph, dph, &
               plcl, cpn, iflag, cbmf)
       endif

       ! --- MIXING

       if (iflag_con == 3) then
          CALL cv3_mixing(klon, ncum, klev, klev, icb, nk, inb, t, q, qs, u, &
               v, h, lv, hp, ep, clw, m, sig, ment, qent, uent, vent, nent, &
               sij, elij, ments, qents)
       else
          ! iflag_con == 4
          CALL cv_mixing(klon, ncum, klev, icb, nk, inb, inbis, ph, t, q, qs, &
               u, v, h, lv, qnk, hp, tv, tvp, ep, clw, cbmf, m, ment, qent, &
               uent, vent, nent, sij, elij)
       endif

       ! --- UNSATURATED (PRECIPITATING) DOWNDRAFTS

       if (iflag_con == 3) then
          CALL cv3_unsat(klon, ncum, klev, klev, icb, inb, t, q, qs, gz, u, &
               v, p, ph, th, tv, lv, cpn, ep, sigp, clw, m, ment, elij, delt, &
               plcl, mp, qp, up, vp, wt, water, evap, b)! na->klev
       else
          ! iflag_con == 4
          CALL cv_unsat(klon, ncum, klev, inb, t, q, qs, gz, u, v, p, ph, h, &
               lv, ep, sigp, clw, m, ment, elij, iflag, mp, qp, up, vp, wt, &
               water, evap)
       endif

       ! --- YIELD
       ! (tendencies, precipitation, variables of interface with other
       ! processes, etc)

       if (iflag_con == 3) then
          CALL cv3_yield(klon, ncum, klev, klev, icb, inb, delt, t, q, u, v, &
               gz, p, ph, h, hp, lv, cpn, th, ep, clw, m, tp, mp, qp, up, vp, &
               wt, water, evap, b, ment, qent, uent, vent, nent, elij, sig, &
               tv, tvp, iflag, precip, VPrecip, ft, fq, fu, fv, upwd, dnwd, &
               dnwd0, ma, mike, tls, tps, qcondc, wd)! na->klev
       else
          ! iflag_con == 4
          CALL cv_yield(klon, ncum, klev, nk, icb, inb, delt, t, q, u, v, gz, &
               p, ph, h, hp, lv, cpn, ep, clw, frac, m, mp, qp, up, vp, wt, &
               water, evap, ment, qent, uent, vent, nent, elij, tv, tvp, &
               iflag, wd, qprime, tprime, precip, cbmf, ft, fq, fu, fv, Ma, &
               qcondc)
       endif

       ! --- passive tracers

       if (iflag_con == 3) CALL cv3_tracer(klon, ncum, klev, ment, sij, da, phi)

       ! --- UNCOMPRESS THE FIELDS

       ! set iflag1 = 42 for non convective points
       do i = 1, klon
          iflag1(i) = 42
       end do

       if (iflag_con == 3) then
          CALL cv3_uncompress(idcum(:ncum), iflag, precip, VPrecip, sig, w0, &
               ft, fq, fu, fv, inb, Ma, upwd, dnwd, dnwd0, qcondc, wd, cape, &
               da, phi, mp, iflag1, precip1, VPrecip1, sig1, w01, ft1, fq1, &
               fu1, fv1, inb1, Ma1, upwd1, dnwd1, dnwd01, qcondc1, wd1, &
               cape1, da1, phi1, mp1)
       else
          ! iflag_con == 4
          CALL cv_uncompress(idcum(:ncum), iflag, precip, cbmf, ft, fq, fu, &
               fv, Ma, qcondc, iflag1, precip1, cbmf1, ft1, fq1, fu1, fv1, &
               Ma1, qcondc1)
       endif
    ENDIF ! ncum>0

  end SUBROUTINE cv_driver

end module cv_driver_m
