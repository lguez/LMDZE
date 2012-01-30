module cv_driver_m

  implicit none

contains

  SUBROUTINE cv_driver(len, nd, ndp1, ntra, iflag_con, t1, q1, qs1, u1, v1, &
       tra1, p1, ph1, iflag1, ft1, fq1, fu1, fv1, ftra1, precip1, VPrecip1, &
       cbmf1, sig1, w01, icb1, inb1, delt, Ma1, upwd1, dnwd1, dnwd01, &
       qcondc1, wd1, cape1, da1, phi1, mp1)

    ! From LMDZ4/libf/phylmd/cv_driver.F, version 1.3 2005/04/15 12:36:17

    use dimens_m
    use dimphy
    !
    ! PARAMETERS:
    !      Name            Type         Usage            Description
    !   ----------      ----------     -------  ----------------------------
    !
    !      len           Integer        Input        first (i) dimension
    !      nd            Integer        Input        vertical (k) dimension
    !      ndp1          Integer        Input        nd + 1
    !      ntra          Integer        Input        number of tracors
    !      iflag_con     Integer        Input        version of convect (3/4)
    !      t1            Real           Input        temperature
    !      q1            Real           Input        specific hum
    !      qs1           Real           Input        sat specific hum
    !      u1            Real           Input        u-wind
    !      v1            Real           Input        v-wind
    !      tra1          Real           Input        tracors
    !      p1            Real           Input        full level pressure
    !      ph1           Real           Input        half level pressure
    !      iflag1        Integer        Output       flag for Emanuel conditions
    !      ft1           Real           Output       temp tend
    !      fq1           Real           Output       spec hum tend
    !      fu1           Real           Output       u-wind tend
    !      fv1           Real           Output       v-wind tend
    !      ftra1         Real           Output       tracor tend
    !      precip1       Real           Output       precipitation
    !      VPrecip1      Real           Output       vertical profile of precipitations
    !      cbmf1         Real           Output       cloud base mass flux
    !      sig1          Real           In/Out       section adiabatic updraft
    !      w01           Real           In/Out       vertical velocity within adiab updraft
    !      delt          Real           Input        time step
    !      Ma1           Real           Output       mass flux adiabatic updraft
    !      upwd1         Real           Output       total upward mass flux (adiab+mixed)
    !      dnwd1         Real           Output       saturated downward mass flux (mixed)
    !      dnwd01        Real           Output       unsaturated downward mass flux
    !      qcondc1       Real           Output       in-cld mixing ratio of condensed water
    !      wd1           Real           Output       downdraft velocity scale for sfc fluxes
    !      cape1         Real           Output       CAPE
    !
    ! S. Bony, Mar 2002:
    !     * Several modules corresponding to different physical processes
    !     * Several versions of convect may be used:
    !        - iflag_con=3: version lmd  (previously named convect3)
    !        - iflag_con=4: version 4.3b (vect. version, previously convect1/2)
    !   + tard:    - iflag_con=5: version lmd with ice (previously named convectg)
    ! S. Bony, Oct 2002:
    !     * Vectorization of convect3 (ie version lmd)
    !
    !..............................END PROLOGUE.............................
    !
    !

    integer len
    integer nd
    integer ndp1
    integer noff
    integer, intent(in):: iflag_con
    integer ntra
    real, intent(in):: t1(len, nd)
    real q1(len, nd)
    real qs1(len, nd)
    real u1(len, nd)
    real v1(len, nd)
    real p1(len, nd)
    real ph1(len, ndp1)
    integer iflag1(len)
    real ft1(len, nd)
    real fq1(len, nd)
    real fu1(len, nd)
    real fv1(len, nd)
    real precip1(len)
    real cbmf1(len)
    real VPrecip1(len, nd+1)
    real Ma1(len, nd)
    real upwd1(len, nd)
    real dnwd1(len, nd)
    real dnwd01(len, nd)

    real qcondc1(len, nd)     ! cld
    real wd1(len)            ! gust
    real cape1(len)

    real da1(len, nd), phi1(len, nd, nd), mp1(len, nd)
    real da(len, nd), phi(len, nd, nd), mp(len, nd)
    real, intent(in):: tra1(len, nd, ntra)
    real ftra1(len, nd, ntra)

    real, intent(in):: delt

    !-------------------------------------------------------------------
    ! --- ARGUMENTS
    !-------------------------------------------------------------------
    ! --- On input:
    !
    !  t:   Array of absolute temperature (K) of dimension ND, with first
    !       index corresponding to lowest model level. Note that this array
    !       will be altered by the subroutine if dry convective adjustment
    !       occurs and if IPBL is not equal to 0.
    !
    !  q:   Array of specific humidity (gm/gm) of dimension ND, with first
    !       index corresponding to lowest model level. Must be defined
    !       at same grid levels as T. Note that this array will be altered
    !       if dry convective adjustment occurs and if IPBL is not equal to 0.
    !
    !  qs:  Array of saturation specific humidity of dimension ND, with first
    !       index corresponding to lowest model level. Must be defined
    !       at same grid levels as T. Note that this array will be altered
    !       if dry convective adjustment occurs and if IPBL is not equal to 0.
    !
    !  u:   Array of zonal wind velocity (m/s) of dimension ND, witth first
    !       index corresponding with the lowest model level. Defined at
    !       same levels as T. Note that this array will be altered if
    !       dry convective adjustment occurs and if IPBL is not equal to 0.
    !
    !  v:   Same as u but for meridional velocity.
    !
    !  tra: Array of passive tracer mixing ratio, of dimensions (ND, NTRA),
    !       where NTRA is the number of different tracers. If no
    !       convective tracer transport is needed, define a dummy
    !       input array of dimension (ND, 1). Tracers are defined at
    !       same vertical levels as T. Note that this array will be altered
    !       if dry convective adjustment occurs and if IPBL is not equal to 0.
    !
    !  p:   Array of pressure (mb) of dimension ND, with first
    !       index corresponding to lowest model level. Must be defined
    !       at same grid levels as T.
    !
    !  ph:  Array of pressure (mb) of dimension ND+1, with first index
    !       corresponding to lowest level. These pressures are defined at
    !       levels intermediate between those of P, T, Q and QS. The first
    !       value of PH should be greater than (i.e. at a lower level than)
    !       the first value of the array P.
    !
    !  nl:  The maximum number of levels to which convection can penetrate, plus 1.
    !       NL MUST be less than or equal to ND-1.
    !
    !  delt: The model time step (sec) between calls to CONVECT
    !
    !----------------------------------------------------------------------------
    ! ---   On Output:
    !
    !  iflag: An output integer whose value denotes the following:
    !       VALUE   INTERPRETATION
    !       -----   --------------
    !         0     Moist convection occurs.
    !         1     Moist convection occurs, but a CFL condition
    !               on the subsidence warming is violated. This
    !               does not cause the scheme to terminate.
    !         2     Moist convection, but no precip because ep(inb) lt 0.0001
    !         3     No moist convection because new cbmf is 0 and old cbmf is 0.
    !         4     No moist convection; atmosphere is not
    !               unstable
    !         6     No moist convection because ihmin le minorig.
    !         7     No moist convection because unreasonable
    !               parcel level temperature or specific humidity.
    !         8     No moist convection: lifted condensation
    !               level is above the 200 mb level.
    !         9     No moist convection: cloud base is higher
    !               then the level NL-1.
    !
    !  ft:   Array of temperature tendency (K/s) of dimension ND, defined at same
    !        grid levels as T, Q, QS and P.
    !
    !  fq:   Array of specific humidity tendencies ((gm/gm)/s) of dimension ND,
    !        defined at same grid levels as T, Q, QS and P.
    !
    !  fu:   Array of forcing of zonal velocity (m/s^2) of dimension ND,
    !        defined at same grid levels as T.
    !
    !  fv:   Same as FU, but for forcing of meridional velocity.
    !
    !  ftra: Array of forcing of tracer content, in tracer mixing ratio per
    !        second, defined at same levels as T. Dimensioned (ND, NTRA).
    !
    !  precip: Scalar convective precipitation rate (mm/day).
    !
    !  VPrecip: Vertical profile of convective precipitation (kg/m2/s).
    !
    !  wd:   A convective downdraft velocity scale. For use in surface
    !        flux parameterizations. See convect.ps file for details.
    !
    !  tprime: A convective downdraft temperature perturbation scale (K).
    !          For use in surface flux parameterizations. See convect.ps
    !          file for details.
    !
    !  qprime: A convective downdraft specific humidity
    !          perturbation scale (gm/gm).
    !          For use in surface flux parameterizations. See convect.ps
    !          file for details.
    !
    !  cbmf: The cloud base mass flux ((kg/m**2)/s). THIS SCALAR VALUE MUST
    !        BE STORED BY THE CALLING PROGRAM AND RETURNED TO CONVECT AT
    !        ITS NEXT CALL. That is, the value of CBMF must be "remembered"
    !        by the calling program between calls to CONVECT.
    !
    !  det:   Array of detrainment mass flux of dimension ND.
    !
    !-------------------------------------------------------------------
    !
    !  Local arrays
    !

    integer i, k, n, il, j
    integer icbmax
    integer nk1(klon)
    integer icb1(klon)
    integer inb1(klon)
    integer icbs1(klon)

    real plcl1(klon)
    real tnk1(klon)
    real qnk1(klon)
    real gznk1(klon)
    real pnk1(klon)
    real qsnk1(klon)
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
    real sig1(klon, klev)
    real w01(klon, klev)
    real th1(klon, klev)
    !
    integer ncum
    !
    ! (local) compressed fields:
    !
    integer nloc
    parameter (nloc=klon) ! pour l'instant

    integer idcum(nloc)
    integer iflag(nloc), nk(nloc), icb(nloc)
    integer nent(nloc, klev)
    integer icbs(nloc)
    integer inb(nloc), inbis(nloc)

    real cbmf(nloc), plcl(nloc), tnk(nloc), qnk(nloc), gznk(nloc)
    real t(nloc, klev), q(nloc, klev), qs(nloc, klev)
    real u(nloc, klev), v(nloc, klev)
    real gz(nloc, klev), h(nloc, klev), lv(nloc, klev), cpn(nloc, klev)
    real p(nloc, klev), ph(nloc, klev+1), tv(nloc, klev), tp(nloc, klev)
    real clw(nloc, klev)
    real dph(nloc, klev)
    real pbase(nloc), buoybase(nloc), th(nloc, klev)
    real tvp(nloc, klev)
    real sig(nloc, klev), w0(nloc, klev)
    real hp(nloc, klev), ep(nloc, klev), sigp(nloc, klev)
    real frac(nloc), buoy(nloc, klev)
    real cape(nloc)
    real m(nloc, klev), ment(nloc, klev, klev), qent(nloc, klev, klev)
    real uent(nloc, klev, klev), vent(nloc, klev, klev)
    real ments(nloc, klev, klev), qents(nloc, klev, klev)
    real sij(nloc, klev, klev), elij(nloc, klev, klev)
    real qp(nloc, klev), up(nloc, klev), vp(nloc, klev)
    real wt(nloc, klev), water(nloc, klev), evap(nloc, klev)
    real b(nloc, klev), ft(nloc, klev), fq(nloc, klev)
    real fu(nloc, klev), fv(nloc, klev)
    real upwd(nloc, klev), dnwd(nloc, klev), dnwd0(nloc, klev)
    real Ma(nloc, klev), mike(nloc, klev), tls(nloc, klev)
    real tps(nloc, klev), qprime(nloc), tprime(nloc)
    real precip(nloc)
    real VPrecip(nloc, klev+1)
    real tra(nloc, klev, ntra), trap(nloc, klev, ntra)
    real ftra(nloc, klev, ntra), traent(nloc, klev, klev, ntra)
    real qcondc(nloc, klev)  ! cld
    real wd(nloc)           ! gust

    !-------------------------------------------------------------------
    ! --- SET CONSTANTS AND PARAMETERS
    !-------------------------------------------------------------------

    ! -- set simulation flags:
    !   (common cvflag)

    CALL cv_flag

    ! -- set thermodynamical constants:
    !     (common cvthermo)

    CALL cv_thermo(iflag_con)

    ! -- set convect parameters
    !
    !     includes microphysical parameters and parameters that
    !     control the rate of approach to quasi-equilibrium)
    !     (common cvparam)

    if (iflag_con.eq.3) then
       CALL cv3_param(nd, delt)
    endif

    if (iflag_con.eq.4) then
       CALL cv_param(nd)
    endif

    !---------------------------------------------------------------------
    ! --- INITIALIZE OUTPUT ARRAYS AND PARAMETERS
    !---------------------------------------------------------------------

    do k=1, nd
       do  i=1, len
          ft1(i, k)=0.0
          fq1(i, k)=0.0
          fu1(i, k)=0.0
          fv1(i, k)=0.0
          tvp1(i, k)=0.0
          tp1(i, k)=0.0
          clw1(i, k)=0.0
          !ym
          clw(i, k)=0.0
          gz1(i, k) = 0.
          VPrecip1(i, k) = 0.
          Ma1(i, k)=0.0
          upwd1(i, k)=0.0
          dnwd1(i, k)=0.0
          dnwd01(i, k)=0.0
          qcondc1(i, k)=0.0
       end do
    end do

    do  j=1, ntra
       do  k=1, nd
          do  i=1, len
             ftra1(i, k, j)=0.0
          end do
       end do
    end do

    do  i=1, len
       precip1(i)=0.0
       iflag1(i)=0
       wd1(i)=0.0
       cape1(i)=0.0
       VPrecip1(i, nd+1)=0.0
    end do

    if (iflag_con.eq.3) then
       do il=1, len
          sig1(il, nd)=sig1(il, nd)+1.
          sig1(il, nd)=amin1(sig1(il, nd), 12.1)
       enddo
    endif

    !--------------------------------------------------------------------
    ! --- CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY & STATIC ENERGY
    !--------------------------------------------------------------------

    if (iflag_con.eq.3) then
       CALL cv3_prelim(len, nd, ndp1, t1, q1, p1, ph1, lv1, cpn1, tv1, gz1, &
            h1, hm1, th1)! nd->na
    endif

    if (iflag_con.eq.4) then
       CALL cv_prelim(len, nd, ndp1, t1, q1, p1, ph1 &
            , lv1, cpn1, tv1, gz1, h1, hm1)
    endif

    !--------------------------------------------------------------------
    ! --- CONVECTIVE FEED
    !--------------------------------------------------------------------

    if (iflag_con.eq.3) then
       CALL cv3_feed(len, nd, t1, q1, qs1, p1, ph1, hm1, gz1            &
            , nk1, icb1, icbmax, iflag1, tnk1, qnk1, gznk1, plcl1) ! nd->na
    endif

    if (iflag_con.eq.4) then
       CALL cv_feed(len, nd, t1, q1, qs1, p1, hm1, gz1 &
            , nk1, icb1, icbmax, iflag1, tnk1, qnk1, gznk1, plcl1)
    endif

    !--------------------------------------------------------------------
    ! --- UNDILUTE (ADIABATIC) UPDRAFT / 1st part
    ! (up through ICB for convect4, up through ICB+1 for convect3)
    !     Calculates the lifted parcel virtual temperature at nk, the
    !     actual temperature, and the adiabatic liquid water content.
    !--------------------------------------------------------------------

    if (iflag_con.eq.3) then
       CALL cv3_undilute1(len, nd, t1, q1, qs1, gz1, plcl1, p1, nk1, icb1   &
            , tp1, tvp1, clw1, icbs1) ! nd->na
    endif

    if (iflag_con.eq.4) then
       CALL cv_undilute1(len, nd, t1, q1, qs1, gz1, p1, nk1, icb1, icbmax &
            , tp1, tvp1, clw1)
    endif

    !-------------------------------------------------------------------
    ! --- TRIGGERING
    !-------------------------------------------------------------------

    if (iflag_con.eq.3) then
       CALL cv3_trigger(len, nd, icb1, plcl1, p1, th1, tv1, tvp1       &
            , pbase1, buoybase1, iflag1, sig1, w01) ! nd->na
    endif

    if (iflag_con.eq.4) then
       CALL cv_trigger(len, nd, icb1, cbmf1, tv1, tvp1, iflag1)
    endif

    !=====================================================================
    ! --- IF THIS POINT IS REACHED, MOIST CONVECTIVE ADJUSTMENT IS NECESSARY
    !=====================================================================

    ncum=0
    do  i=1, len
       if(iflag1(i).eq.0)then
          ncum=ncum+1
          idcum(ncum)=i
       endif
    end do

    !       print*, 'klon, ncum = ', len, ncum

    IF (ncum.gt.0) THEN

       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       ! --- COMPRESS THE FIELDS
       !        (-> vectorization over convective gridpoints)
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       if (iflag_con.eq.3) then
          CALL cv3_compress( len, nloc, ncum, nd, ntra &
               , iflag1, nk1, icb1, icbs1 &
               , plcl1, tnk1, qnk1, gznk1, pbase1, buoybase1 &
               , t1, q1, qs1, u1, v1, gz1, th1 &
               , tra1 &
               , h1, lv1, cpn1, p1, ph1, tv1, tp1, tvp1, clw1  &
               , sig1, w01 &
               , iflag, nk, icb, icbs &
               , plcl, tnk, qnk, gznk, pbase, buoybase &
               , t, q, qs, u, v, gz, th &
               , tra &
               , h, lv, cpn, p, ph, tv, tp, tvp, clw  &
               , sig, w0  )
       endif

       if (iflag_con.eq.4) then
          CALL cv_compress( len, nloc, ncum, nd &
               , iflag1, nk1, icb1 &
               , cbmf1, plcl1, tnk1, qnk1, gznk1 &
               , t1, q1, qs1, u1, v1, gz1 &
               , h1, lv1, cpn1, p1, ph1, tv1, tp1, tvp1, clw1 &
               , iflag, nk, icb &
               , cbmf, plcl, tnk, qnk, gznk &
               , t, q, qs, u, v, gz, h, lv, cpn, p, ph, tv, tp, tvp, clw  &
               , dph )
       endif

       !-------------------------------------------------------------------
       ! --- UNDILUTE (ADIABATIC) UPDRAFT / second part :
       ! ---   FIND THE REST OF THE LIFTED PARCEL TEMPERATURES
       ! ---   &
       ! ---   COMPUTE THE PRECIPITATION EFFICIENCIES AND THE
       ! ---   FRACTION OF PRECIPITATION FALLING OUTSIDE OF CLOUD
       ! ---   &
       ! ---   FIND THE LEVEL OF NEUTRAL BUOYANCY
       !-------------------------------------------------------------------

       if (iflag_con.eq.3) then
          CALL cv3_undilute2(nloc, ncum, nd, icb, icbs, nk         &
               , tnk, qnk, gznk, t, q, qs, gz &
               , p, h, tv, lv, pbase, buoybase, plcl &
               , inb, tp, tvp, clw, hp, ep, sigp, buoy) !na->nd
       endif

       if (iflag_con.eq.4) then
          CALL cv_undilute2(nloc, ncum, nd, icb, nk &
               , tnk, qnk, gznk, t, q, qs, gz &
               , p, dph, h, tv, lv &
               , inb, inbis, tp, tvp, clw, hp, ep, sigp, frac)
       endif

       !-------------------------------------------------------------------
       ! --- CLOSURE
       !-------------------------------------------------------------------

       if (iflag_con.eq.3) then
          CALL cv3_closure(nloc, ncum, nd, icb, inb               &
               , pbase, p, ph, tv, buoy &
               , sig, w0, cape, m) ! na->nd
       endif

       if (iflag_con.eq.4) then
          CALL cv_closure(nloc, ncum, nd, nk, icb &
               , tv, tvp, p, ph, dph, plcl, cpn &
               , iflag, cbmf)
       endif

       !-------------------------------------------------------------------
       ! --- MIXING
       !-------------------------------------------------------------------

       if (iflag_con.eq.3) then
          CALL cv3_mixing(nloc, ncum, nd, nd, ntra, icb, nk, inb     &
               , ph, t, q, qs, u, v, tra, h, lv, qnk &
               , hp, tv, tvp, ep, clw, m, sig &
               , ment, qent, uent, vent, nent, sij, elij, ments, qents, traent)! na->nd
       endif

       if (iflag_con.eq.4) then
          CALL cv_mixing(nloc, ncum, nd, icb, nk, inb, inbis &
               , ph, t, q, qs, u, v, h, lv, qnk &
               , hp, tv, tvp, ep, clw, cbmf &
               , m, ment, qent, uent, vent, nent, sij, elij)
       endif

       !-------------------------------------------------------------------
       ! --- UNSATURATED (PRECIPITATING) DOWNDRAFTS
       !-------------------------------------------------------------------

       if (iflag_con.eq.3) then
          CALL cv3_unsat(nloc, ncum, nd, nd, ntra, icb, inb     &
               , t, q, qs, gz, u, v, tra, p, ph &
               , th, tv, lv, cpn, ep, sigp, clw &
               , m, ment, elij, delt, plcl &
               , mp, qp, up, vp, trap, wt, water, evap, b)! na->nd
       endif

       if (iflag_con.eq.4) then
          CALL cv_unsat(nloc, ncum, nd, inb, t, q, qs, gz, u, v, p, ph &
               , h, lv, ep, sigp, clw, m, ment, elij &
               , iflag, mp, qp, up, vp, wt, water, evap)
       endif

       !-------------------------------------------------------------------
       ! --- YIELD
       !     (tendencies, precipitation, variables of interface with other
       !      processes, etc)
       !-------------------------------------------------------------------

       if (iflag_con.eq.3) then
          CALL cv3_yield(nloc, ncum, nd, nd, ntra             &
               , icb, inb, delt &
               , t, q, u, v, tra, gz, p, ph, h, hp, lv, cpn, th &
               , ep, clw, m, tp, mp, qp, up, vp, trap &
               , wt, water, evap, b &
               , ment, qent, uent, vent, nent, elij, traent, sig &
               , tv, tvp &
               , iflag, precip, VPrecip, ft, fq, fu, fv, ftra &
               , upwd, dnwd, dnwd0, ma, mike, tls, tps, qcondc, wd)! na->nd
       endif

       if (iflag_con.eq.4) then
          CALL cv_yield(nloc, ncum, nd, nk, icb, inb, delt &
               , t, q, u, v, gz, p, ph, h, hp, lv, cpn &
               , ep, clw, frac, m, mp, qp, up, vp &
               , wt, water, evap &
               , ment, qent, uent, vent, nent, elij &
               , tv, tvp &
               , iflag, wd, qprime, tprime &
               , precip, cbmf, ft, fq, fu, fv, Ma, qcondc)
       endif

       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       ! --- passive tracers
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       if (iflag_con.eq.3) then
          CALL cv3_tracer(nloc, len, ncum, nd, nd, &
               ment, sij, da, phi)
       endif

       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       ! --- UNCOMPRESS THE FIELDS
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       ! set iflag1 =42 for non convective points
       do  i=1, len
          iflag1(i)=42
       end do
       !
       if (iflag_con.eq.3) then
          CALL cv3_uncompress(nloc, len, ncum, nd, ntra, idcum &
               , iflag &
               , precip, VPrecip, sig, w0 &
               , ft, fq, fu, fv, ftra &
               , inb  &
               , Ma, upwd, dnwd, dnwd0, qcondc, wd, cape &
               , da, phi, mp &
               , iflag1 &
               , precip1, VPrecip1, sig1, w01 &
               , ft1, fq1, fu1, fv1, ftra1 &
               , inb1 &
               , Ma1, upwd1, dnwd1, dnwd01, qcondc1, wd1, cape1  &
               , da1, phi1, mp1)
       endif

       if (iflag_con.eq.4) then
          CALL cv_uncompress(nloc, len, ncum, nd, idcum &
               , iflag &
               , precip, cbmf &
               , ft, fq, fu, fv &
               , Ma, qcondc             &
               , iflag1 &
               , precip1, cbmf1 &
               , ft1, fq1, fu1, fv1 &
               , Ma1, qcondc1 )
       endif
    ENDIF ! ncum>0

  end SUBROUTINE cv_driver

end module cv_driver_m