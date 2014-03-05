
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/ppm3d.F,v 1.1.1.1 2004/05/19
! 12:53:07 lmdzadmin Exp $


! From lin@explorer.gsfc.nasa.gov Wed Apr 15 17:44:44 1998
! Date: Wed, 15 Apr 1998 11:37:03 -0400
! From: lin@explorer.gsfc.nasa.gov
! To: Frederic.Hourdin@lmd.jussieu.fr
! Subject: 3D transport module of the GSFC CTM and GEOS GCM


! This code is sent to you by S-J Lin, DAO, NASA-GSFC

! Note: this version is intended for machines like CRAY
! -90. No multitasking directives implemented.


! ********************************************************************

! TransPort Core for Goddard Chemistry Transport Model (G-CTM), Goddard
! Earth Observing System General Circulation Model (GEOS-GCM), and Data
! Assimilation System (GEOS-DAS).

! ********************************************************************

! Purpose: given horizontal winds on  a hybrid sigma-p surfaces,
! one call to tpcore updates the 3-D mixing ratio
! fields one time step (NDT). [vertical mass flux is computed
! internally consistent with the discretized hydrostatic mass
! continuity equation of the C-Grid GEOS-GCM (for IGD=1)].

! Schemes: Multi-dimensional Flux Form Semi-Lagrangian (FFSL) scheme based
! on the van Leer or PPM.
! (see Lin and Rood 1996).
! Version 4.5
! Last modified: Dec. 5, 1996
! Major changes from version 4.0: a more general vertical hybrid sigma-
! pressure coordinate.
! Subroutines modified: xtp, ytp, fzppm, qckxyz
! Subroutines deleted: vanz

! Author: Shian-Jiann Lin
! mail address:
! Shian-Jiann Lin*
! Code 910.3, NASA/GSFC, Greenbelt, MD 20771
! Phone: 301-286-9540
! E-mail: lin@dao.gsfc.nasa.gov

! *affiliation:
! Joint Center for Earth Systems Technology
! The University of Maryland Baltimore County
! NASA - Goddard Space Flight Center
! References:

! 1. Lin, S.-J., and R. B. Rood, 1996: Multidimensional flux form semi-
! Lagrangian transport schemes. Mon. Wea. Rev., 124, 2046-2070.

! 2. Lin, S.-J., W. C. Chao, Y. C. Sud, and G. K. Walker, 1994: A class of
! the van Leer-type transport schemes and its applications to the moist-
! ure transport in a General Circulation Model. Mon. Wea. Rev., 122,
! 1575-1593.

! ****6***0*********0*********0*********0*********0*********0**********72

SUBROUTINE ppm3d(igd, q, ps1, ps2, u, v, w, ndt, iord, jord, kord, nc, imr, &
    jnp, j1, nlay, ap, bp, pt, ae, fill, dum, umax)

  ! implicit none

  ! rajout de déclarations
  ! integer Jmax,kmax,ndt0,nstep,k,j,i,ic,l,js,jn,imh,iad,jad,krd
  ! integer iu,iiu,j2,jmr,js0,jt
  ! real dtdy,dtdy5,rcap,iml,jn0,imjm,pi,dl,dp
  ! real dt,cr1,maxdt,ztc,d5,sum1,sum2,ru

  ! ********************************************************************

  ! =============
  ! INPUT:
  ! =============

  ! Q(IMR,JNP,NLAY,NC): mixing ratios at current time (t)
  ! NC: total number of constituents
  ! IMR: first dimension (E-W); number of Grid intervals in E-W is IMR
  ! JNP: 2nd dimension (N-S); number of Grid intervals in N-S is JNP-1
  ! NLAY: 3rd dimension (number of layers); vertical index increases from 1
  ! at
  ! the model top to NLAY near the surface (see fig. below).
  ! It is assumed that 6 <= NLAY <= JNP (for dynamic memory allocation)

  ! PS1(IMR,JNP): surface pressure at current time (t)
  ! PS2(IMR,JNP): surface pressure at mid-time-level (t+NDT/2)
  ! PS2 is replaced by the predicted PS (at t+NDT) on output.
  ! Note: surface pressure can have any unit or can be multiplied by any
  ! const.

  ! The pressure at layer edges are defined as follows:

  ! p(i,j,k) = AP(k)*PT  +  BP(k)*PS(i,j)          (1)

  ! Where PT is a constant having the same unit as PS.
  ! AP and BP are unitless constants given at layer edges
  ! defining the vertical coordinate.
  ! BP(1) = 0., BP(NLAY+1) = 1.
  ! The pressure at the model top is PTOP = AP(1)*PT

  ! For pure sigma system set AP(k) = 1 for all k, PT = PTOP,
  ! BP(k) = sige(k) (sigma at edges), PS = Psfc - PTOP.

  ! Note: the sigma-P coordinate is a subset of Eq. 1, which in turn
  ! is a subset of the following even more general sigma-P-thelta coord.
  ! currently under development.
  ! p(i,j,k) = (AP(k)*PT + BP(k)*PS(i,j))/(D(k)-C(k)*TE**(-1/kapa))

  ! /////////////////////////////////
  ! / \ ------------- PTOP --------------  AP(1), BP(1)
  ! |
  ! delp(1)    |  ........... Q(i,j,1) ............
  ! |
  ! W(1)    \ / ---------------------------------  AP(2), BP(2)



  ! W(k-1)   / \ ---------------------------------  AP(k), BP(k)
  ! |
  ! delp(K)    |  ........... Q(i,j,k) ............
  ! |
  ! W(k)    \ / ---------------------------------  AP(k+1), BP(k+1)



  ! / \ ---------------------------------  AP(NLAY), BP(NLAY)
  ! |
  ! delp(NLAY)   |  ........... Q(i,j,NLAY) .........
  ! |
  ! W(NLAY)=0  \ / ------------- surface ----------- AP(NLAY+1), BP(NLAY+1)
  ! //////////////////////////////////

  ! U(IMR,JNP,NLAY) & V(IMR,JNP,NLAY):winds (m/s) at mid-time-level (t+NDT/2)
  ! U and V may need to be polar filtered in advance in some cases.

  ! IGD:      grid type on which winds are defined.
  ! IGD = 0:  A-Grid  [all variables defined at the same point from south
  ! pole (j=1) to north pole (j=JNP) ]

  ! IGD = 1  GEOS-GCM C-Grid
  ! [North]

  ! V(i,j)
  ! |
  ! |
  ! |
  ! U(i-1,j)---Q(i,j)---U(i,j) [EAST]
  ! |
  ! |
  ! |
  ! V(i,j-1)

  ! U(i,  1) is defined at South Pole.
  ! V(i,  1) is half grid north of the South Pole.
  ! V(i,JMR) is half grid south of the North Pole.

  ! V must be defined at j=1 and j=JMR if IGD=1
  ! V at JNP need not be given.

  ! NDT: time step in seconds (need not be constant during the course of
  ! the integration). Suggested value: 30 min. for 4x5, 15 min. for 2x2.5
  ! (Lat-Lon) resolution. Smaller values are recommanded if the model
  ! has a well-resolved stratosphere.

  ! J1 defines the size of the polar cap:
  ! South polar cap edge is located at -90 + (j1-1.5)*180/(JNP-1) deg.
  ! North polar cap edge is located at  90 - (j1-1.5)*180/(JNP-1) deg.
  ! There are currently only two choices (j1=2 or 3).
  ! IMR must be an even integer if j1 = 2. Recommended value: J1=3.

  ! IORD, JORD, and KORD are integers controlling various options in E-W,
  ! N-S,
  ! and vertical transport, respectively. Recommended values for positive
  ! definite scalars: IORD=JORD=3, KORD=5. Use KORD=3 for non-
  ! positive definite scalars or when linear correlation between constituents
  ! is to be maintained.

  ! _ORD=
  ! 1: 1st order upstream scheme (too diffusive, not a useful option; it
  ! can be used for debugging purposes; this is THE only known "linear"
  ! monotonic advection scheme.).
  ! 2: 2nd order van Leer (full monotonicity constraint;
  ! see Lin et al 1994, MWR)
  ! 3: monotonic PPM* (slightly improved PPM of Collela & Woodward 1984)
  ! 4: semi-monotonic PPM (same as 3, but overshoots are allowed)
  ! 5: positive-definite PPM (constraint on the subgrid distribution is
  ! only strong enough to prevent generation of negative values;
  ! both overshoots & undershoots are possible).
  ! 6: un-constrained PPM (nearly diffusion free; slightly faster but
  ! positivity not quaranteed. Use this option only when the fields
  ! and winds are very smooth).

  ! *PPM: Piece-wise Parabolic Method

  ! Note that KORD <=2 options are no longer supported. DO not use option 4
  ! or 5.
  ! for non-positive definite scalars (such as Ertel Potential Vorticity).

  ! The implicit numerical diffusion decreases as _ORD increases.
  ! The last two options (ORDER=5, 6) should only be used when there is
  ! significant explicit diffusion (such as a turbulence parameterization).
  ! You
  ! might get dispersive results otherwise.
  ! No filter of any kind is applied to the constituent fields here.

  ! AE: Radius of the sphere (meters).
  ! Recommended value for the planet earth: 6.371E6

  ! fill(logical):   flag to do filling for negatives (see note below).

  ! Umax: Estimate (upper limit) of the maximum U-wind speed (m/s).
  ! (220 m/s is a good value for troposphere model; 280 m/s otherwise)

  ! =============
  ! Output
  ! =============

  ! Q: mixing ratios at future time (t+NDT) (original values are
  ! over-written)
  ! W(NLAY): large-scale vertical mass flux as diagnosed from the hydrostatic
  ! relationship. W will have the same unit as PS1 and PS2 (eg, mb).
  ! W must be divided by NDT to get the correct mass-flux unit.
  ! The vertical Courant number C = W/delp_UPWIND, where delp_UPWIND
  ! is the pressure thickness in the "upwind" direction. For example,
  ! C(k) = W(k)/delp(k)   if W(k) > 0;
  ! C(k) = W(k)/delp(k+1) if W(k) < 0.
  ! ( W > 0 is downward, ie, toward surface)
  ! PS2: predicted PS at t+NDT (original values are over-written)

  ! ********************************************************************
  ! NOTES:
  ! This forward-in-time upstream-biased transport scheme reduces to
  ! the 2nd order center-in-time center-in-space mass continuity eqn.
  ! if Q = 1 (constant fields will remain constant). This also ensures
  ! that the computed vertical velocity to be identical to GEOS-1 GCM
  ! for on-line transport.

  ! A larger polar cap is used if j1=3 (recommended for C-Grid winds or when
  ! winds are noisy near poles).

  ! Flux-Form Semi-Lagrangian transport in the East-West direction is used
  ! when and where Courant number is greater than one.

  ! The user needs to change the parameter Jmax or Kmax if the resolution
  ! is greater than 0.5 deg in N-S or 150 layers in the vertical direction.
  ! (this TransPort Core is otherwise resolution independent and can be used
  ! as a library routine).

  ! PPM is 4th order accurate when grid spacing is uniform (x & y); 3rd
  ! order accurate for non-uniform grid (vertical sigma coord.).

  ! Time step is limitted only by transport in the meridional direction.
  ! (the FFSL scheme is not implemented in the meridional direction).

  ! Since only 1-D limiters are applied, negative values could
  ! potentially be generated when large time step is used and when the
  ! initial fields contain discontinuities.
  ! This does not necessarily imply the integration is unstable.
  ! These negatives are typically very small. A filling algorithm is
  ! activated if the user set "fill" to be true.

  ! The van Leer scheme used here is nearly as accurate as the original PPM
  ! due to the use of a 4th order accurate reference slope. The PPM imple-
  ! mented here is an improvement over the original and is also based on
  ! the 4th order reference slope.

  ! ****6***0*********0*********0*********0*********0*********0**********72

  ! User modifiable parameters

  PARAMETER (jmax=361, kmax=150)

  ! ****6***0*********0*********0*********0*********0*********0**********72

  ! Input-Output arrays


  REAL q(imr, jnp, nlay, nc), ps1(imr, jnp), ps2(imr, jnp), &
    u(imr, jnp, nlay), v(imr, jnp, nlay), ap(nlay+1), bp(nlay+1), &
    w(imr, jnp, nlay), ndt, val(nlay), umax
  INTEGER igd, iord, jord, kord, nc, imr, jnp, j1, nlay, ae
  INTEGER imrd2
  REAL pt
  LOGICAL cross, fill, dum

  ! Local dynamic arrays

  REAL crx(imr, jnp), cry(imr, jnp), xmass(imr, jnp), ymass(imr, jnp), &
    fx1(imr+1), dpi(imr, jnp, nlay), delp1(imr, jnp, nlay), &
    wk1(imr, jnp, nlay), pu(imr, jnp), pv(imr, jnp), dc2(imr, jnp), &
    delp2(imr, jnp, nlay), dq(imr, jnp, nlay, nc), va(imr, jnp), &
    ua(imr, jnp), qtmp(-imr:2*imr)

  ! Local static  arrays

  REAL dtdx(jmax), dtdx5(jmax), acosp(jmax), cosp(jmax), cose(jmax), &
    dap(kmax), dbk(kmax)
  DATA ndt0, nstep/0, 0/
  DATA cross/.TRUE./
  SAVE dtdy, dtdy5, rcap, js0, jn0, iml, dtdx, dtdx5, acosp, cosp, cose, dap, &
    dbk


  jmr = jnp - 1
  imjm = imr*jnp
  j2 = jnp - j1 + 1
  nstep = nstep + 1

  ! *********** Initialization **********************
  IF (nstep==1) THEN

    WRITE (6, *) '------------------------------------ '
    WRITE (6, *) 'NASA/GSFC Transport Core Version 4.5'
    WRITE (6, *) '------------------------------------ '

    WRITE (6, *) 'IMR=', imr, ' JNP=', jnp, ' NLAY=', nlay, ' j1=', j1
    WRITE (6, *) 'NC=', nc, iord, jord, kord, ndt

    ! controles sur les parametres
    IF (nlay<6) THEN
      WRITE (6, *) 'NLAY must be >= 6'
      STOP
    END IF
    IF (jnp<nlay) THEN
      WRITE (6, *) 'JNP must be >= NLAY'
      STOP
    END IF
    imrd2 = mod(imr, 2)
    IF (j1==2 .AND. imrd2/=0) THEN
      WRITE (6, *) 'if j1=2 IMR must be an even integer'
      STOP
    END IF


    IF (jmax<jnp .OR. kmax<nlay) THEN
      WRITE (6, *) 'Jmax or Kmax is too small'
      STOP
    END IF

    DO k = 1, nlay
      dap(k) = (ap(k+1)-ap(k))*pt
      dbk(k) = bp(k+1) - bp(k)
    END DO

    pi = 4.*atan(1.)
    dl = 2.*pi/float(imr)
    dp = pi/float(jmr)

    IF (igd==0) THEN
      ! Compute analytic cosine at cell edges
      CALL cosa(cosp, cose, jnp, pi, dp)
    ELSE
      ! Define cosine consistent with GEOS-GCM (using dycore2.0 or later)
      CALL cosc(cosp, cose, jnp, pi, dp)
    END IF

    DO j = 2, jmr
      acosp(j) = 1./cosp(j)
    END DO

    ! Inverse of the Scaled polar cap area.

    rcap = dp/(imr*(1.-cos((j1-1.5)*dp)))
    acosp(1) = rcap
    acosp(jnp) = rcap
  END IF

  IF (ndt0/=ndt) THEN
    dt = ndt
    ndt0 = ndt

    IF (umax<180.) THEN
      WRITE (6, *) 'Umax may be too small!'
    END IF
    cr1 = abs(umax*dt)/(dl*ae)
    maxdt = dp*ae/abs(umax) + 0.5
    WRITE (6, *) 'Largest time step for max(V)=', umax, ' is ', maxdt
    IF (maxdt<abs(ndt)) THEN
      WRITE (6, *) 'Warning!!! NDT maybe too large!'
    END IF

    IF (cr1>=0.95) THEN
      js0 = 0
      jn0 = 0
      iml = imr - 2
      ztc = 0.
    ELSE
      ztc = acos(cr1)*(180./pi)

      js0 = float(jmr)*(90.-ztc)/180. + 2
      js0 = max(js0, j1+1)
      iml = min(6*js0/(j1-1)+2, 4*imr/5)
      jn0 = jnp - js0 + 1
    END IF


    DO j = 2, jmr
      dtdx(j) = dt/(dl*ae*cosp(j))

      dtdx5(j) = 0.5*dtdx(j)
    END DO


    dtdy = dt/(ae*dp)
    dtdy5 = 0.5*dtdy

  END IF

  ! *********** End Initialization **********************

  ! delp = pressure thickness: the psudo-density in a hydrostatic system.
  DO k = 1, nlay
    DO j = 1, jnp
      DO i = 1, imr
        delp1(i, j, k) = dap(k) + dbk(k)*ps1(i, j)
        delp2(i, j, k) = dap(k) + dbk(k)*ps2(i, j)
      END DO
    END DO
  END DO


  IF (j1/=2) THEN
    DO ic = 1, nc
      DO l = 1, nlay
        DO i = 1, imr
          q(i, 2, l, ic) = q(i, 1, l, ic)
          q(i, jmr, l, ic) = q(i, jnp, l, ic)
        END DO
      END DO
    END DO
  END IF

  ! Compute "tracer density"
  DO ic = 1, nc
    DO k = 1, nlay
      DO j = 1, jnp
        DO i = 1, imr
          dq(i, j, k, ic) = q(i, j, k, ic)*delp1(i, j, k)
        END DO
      END DO
    END DO
  END DO

  DO k = 1, nlay

    IF (igd==0) THEN
      ! Convert winds on A-Grid to Courant number on C-Grid.
      CALL a2c(u(1,1,k), v(1,1,k), imr, jmr, j1, j2, crx, cry, dtdx5, dtdy5)
    ELSE
      ! Convert winds on C-grid to Courant number
      DO j = j1, j2
        DO i = 2, imr
          crx(i, j) = dtdx(j)*u(i-1, j, k)
        END DO
      END DO


      DO j = j1, j2
        crx(1, j) = dtdx(j)*u(imr, j, k)
      END DO

      DO i = 1, imr*jmr
        cry(i, 2) = dtdy*v(i, 1, k)
      END DO
    END IF

    ! Determine JS and JN
    js = j1
    jn = j2

    DO j = js0, j1 + 1, -1
      DO i = 1, imr
        IF (abs(crx(i,j))>1.) THEN
          js = j
          GO TO 2222
        END IF
      END DO
    END DO

2222 CONTINUE
    DO j = jn0, j2 - 1
      DO i = 1, imr
        IF (abs(crx(i,j))>1.) THEN
          jn = j
          GO TO 2233
        END IF
      END DO
    END DO
2233 CONTINUE

    IF (j1/=2) THEN ! Enlarged polar cap.
      DO i = 1, imr
        dpi(i, 2, k) = 0.
        dpi(i, jmr, k) = 0.
      END DO
    END IF

    ! ******* Compute horizontal mass fluxes ************

    ! N-S component
    DO j = j1, j2 + 1
      d5 = 0.5*cose(j)
      DO i = 1, imr
        ymass(i, j) = cry(i, j)*d5*(delp2(i,j,k)+delp2(i,j-1,k))
      END DO
    END DO

    DO j = j1, j2
      DO i = 1, imr
        dpi(i, j, k) = (ymass(i,j)-ymass(i,j+1))*acosp(j)
      END DO
    END DO

    ! Poles
    sum1 = ymass(imr, j1)
    sum2 = ymass(imr, j2+1)
    DO i = 1, imr - 1
      sum1 = sum1 + ymass(i, j1)
      sum2 = sum2 + ymass(i, j2+1)
    END DO

    sum1 = -sum1*rcap
    sum2 = sum2*rcap
    DO i = 1, imr
      dpi(i, 1, k) = sum1
      dpi(i, jnp, k) = sum2
    END DO

    ! E-W component

    DO j = j1, j2
      DO i = 2, imr
        pu(i, j) = 0.5*(delp2(i,j,k)+delp2(i-1,j,k))
      END DO
    END DO

    DO j = j1, j2
      pu(1, j) = 0.5*(delp2(1,j,k)+delp2(imr,j,k))
    END DO

    DO j = j1, j2
      DO i = 1, imr
        xmass(i, j) = pu(i, j)*crx(i, j)
      END DO
    END DO

    DO j = j1, j2
      DO i = 1, imr - 1
        dpi(i, j, k) = dpi(i, j, k) + xmass(i, j) - xmass(i+1, j)
      END DO
    END DO

    DO j = j1, j2
      dpi(imr, j, k) = dpi(imr, j, k) + xmass(imr, j) - xmass(1, j)
    END DO

    DO j = j1, j2
      DO i = 1, imr - 1
        ua(i, j) = 0.5*(crx(i,j)+crx(i+1,j))
      END DO
    END DO

    DO j = j1, j2
      ua(imr, j) = 0.5*(crx(imr,j)+crx(1,j))
    END DO
    ! cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! Rajouts pour LMDZ.3.3
    ! cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    DO i = 1, imr
      DO j = 1, jnp
        va(i, j) = 0.
      END DO
    END DO

    DO i = 1, imr*(jmr-1)
      va(i, 2) = 0.5*(cry(i,2)+cry(i,3))
    END DO

    IF (j1==2) THEN
      imh = imr/2
      DO i = 1, imh
        va(i, 1) = 0.5*(cry(i,2)-cry(i+imh,2))
        va(i+imh, 1) = -va(i, 1)
        va(i, jnp) = 0.5*(cry(i,jnp)-cry(i+imh,jmr))
        va(i+imh, jnp) = -va(i, jnp)
      END DO
      va(imr, 1) = va(1, 1)
      va(imr, jnp) = va(1, jnp)
    END IF

    ! ****6***0*********0*********0*********0*********0*********0**********72
    DO ic = 1, nc

      DO i = 1, imjm
        wk1(i, 1, 1) = 0.
        wk1(i, 1, 2) = 0.
      END DO

      ! E-W advective cross term
      DO j = j1, j2
        IF (j>js .AND. j<jn) GO TO 250

        DO i = 1, imr
          qtmp(i) = q(i, j, k, ic)
        END DO

        DO i = -iml, 0
          qtmp(i) = q(imr+i, j, k, ic)
          qtmp(imr+1-i) = q(1-i, j, k, ic)
        END DO

        DO i = 1, imr
          iu = ua(i, j)
          ru = ua(i, j) - iu
          iiu = i - iu
          IF (ua(i,j)>=0.) THEN
            wk1(i, j, 1) = qtmp(iiu) + ru*(qtmp(iiu-1)-qtmp(iiu))
          ELSE
            wk1(i, j, 1) = qtmp(iiu) + ru*(qtmp(iiu)-qtmp(iiu+1))
          END IF
          wk1(i, j, 1) = wk1(i, j, 1) - qtmp(i)
        END DO
250   END DO

      IF (jn/=0) THEN
        DO j = js + 1, jn - 1

          DO i = 1, imr
            qtmp(i) = q(i, j, k, ic)
          END DO

          qtmp(0) = q(imr, j, k, ic)
          qtmp(imr+1) = q(1, j, k, ic)

          DO i = 1, imr
            iu = i - ua(i, j)
            wk1(i, j, 1) = ua(i, j)*(qtmp(iu)-qtmp(iu+1))
          END DO
        END DO
      END IF
      ! ****6***0*********0*********0*********0*********0*********0**********72
      ! Contribution from the N-S advection
      DO i = 1, imr*(j2-j1+1)
        jt = float(j1) - va(i, j1)
        wk1(i, j1, 2) = va(i, j1)*(q(i,jt,k,ic)-q(i,jt+1,k,ic))
      END DO

      DO i = 1, imjm
        wk1(i, 1, 1) = q(i, 1, k, ic) + 0.5*wk1(i, 1, 1)
        wk1(i, 1, 2) = q(i, 1, k, ic) + 0.5*wk1(i, 1, 2)
      END DO

      IF (cross) THEN
        ! Add cross terms in the vertical direction.
        IF (iord>=2) THEN
          iad = 2
        ELSE
          iad = 1
        END IF

        IF (jord>=2) THEN
          jad = 2
        ELSE
          jad = 1
        END IF
        CALL xadv(imr, jnp, j1, j2, wk1(1,1,2), ua, js, jn, iml, dc2, iad)
        CALL yadv(imr, jnp, j1, j2, wk1(1,1,1), va, pv, w, jad)
        DO j = 1, jnp
          DO i = 1, imr
            q(i, j, k, ic) = q(i, j, k, ic) + dc2(i, j) + pv(i, j)
          END DO
        END DO
      END IF

      CALL xtp(imr, jnp, iml, j1, j2, jn, js, pu, dq(1,1,k,ic), wk1(1,1,2), &
        crx, fx1, xmass, iord)

      CALL ytp(imr, jnp, j1, j2, acosp, rcap, dq(1,1,k,ic), wk1(1,1,1), cry, &
        dc2, ymass, wk1(1,1,3), wk1(1,1,4), wk1(1,1,5), wk1(1,1,6), jord)

    END DO
  END DO

  ! ******* Compute vertical mass flux (same unit as PS) ***********

  ! 1st step: compute total column mass CONVERGENCE.

  DO j = 1, jnp
    DO i = 1, imr
      cry(i, j) = dpi(i, j, 1)
    END DO
  END DO

  DO k = 2, nlay
    DO j = 1, jnp
      DO i = 1, imr
        cry(i, j) = cry(i, j) + dpi(i, j, k)
      END DO
    END DO
  END DO

  DO j = 1, jnp
    DO i = 1, imr

      ! 2nd step: compute PS2 (PS at n+1) using the hydrostatic assumption.
      ! Changes (increases) to surface pressure = total column mass
      ! convergence

      ps2(i, j) = ps1(i, j) + cry(i, j)

      ! 3rd step: compute vertical mass flux from mass conservation
      ! principle.

      w(i, j, 1) = dpi(i, j, 1) - dbk(1)*cry(i, j)
      w(i, j, nlay) = 0.
    END DO
  END DO

  DO k = 2, nlay - 1
    DO j = 1, jnp
      DO i = 1, imr
        w(i, j, k) = w(i, j, k-1) + dpi(i, j, k) - dbk(k)*cry(i, j)
      END DO
    END DO
  END DO

  DO k = 1, nlay
    DO j = 1, jnp
      DO i = 1, imr
        delp2(i, j, k) = dap(k) + dbk(k)*ps2(i, j)
      END DO
    END DO
  END DO

  krd = max(3, kord)
  DO ic = 1, nc

    ! ****6***0*********0*********0*********0*********0*********0**********72

    CALL fzppm(imr, jnp, nlay, j1, dq(1,1,1,ic), w, q(1,1,1,ic), wk1, dpi, &
      dc2, crx, cry, pu, pv, xmass, ymass, delp1, krd)


    IF (fill) CALL qckxyz(dq(1,1,1,ic), dc2, imr, jnp, nlay, j1, j2, cosp, &
      acosp, .FALSE., ic, nstep)

    ! Recover tracer mixing ratio from "density" using predicted
    ! "air density" (pressure thickness) at time-level n+1

    DO k = 1, nlay
      DO j = 1, jnp
        DO i = 1, imr
          q(i, j, k, ic) = dq(i, j, k, ic)/delp2(i, j, k)
        END DO
      END DO
    END DO

    IF (j1/=2) THEN
      DO k = 1, nlay
        DO i = 1, imr
          ! j=1 c'est le pôle Sud, j=JNP c'est le pôle Nord
          q(i, 2, k, ic) = q(i, 1, k, ic)
          q(i, jmr, k, ic) = q(i, jmp, k, ic)
        END DO
      END DO
    END IF
  END DO

  IF (j1/=2) THEN
    DO k = 1, nlay
      DO i = 1, imr
        w(i, 2, k) = w(i, 1, k)
        w(i, jmr, k) = w(i, jnp, k)
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE ppm3d

! ****6***0*********0*********0*********0*********0*********0**********72
SUBROUTINE fzppm(imr, jnp, nlay, j1, dq, wz, p, dc, dqdt, ar, al, a6, flux, &
    wk1, wk2, wz2, delp, kord)
  PARAMETER (kmax=150)
  PARAMETER (r23=2./3., r3=1./3.)
  REAL wz(imr, jnp, nlay), p(imr, jnp, nlay), dc(imr, jnp, nlay), &
    wk1(imr, *), delp(imr, jnp, nlay), dq(imr, jnp, nlay), &
    dqdt(imr, jnp, nlay)
  ! Assuming JNP >= NLAY
  REAL ar(imr, *), al(imr, *), a6(imr, *), flux(imr, *), wk2(imr, *), &
    wz2(imr, *)

  jmr = jnp - 1
  imjm = imr*jnp
  nlaym1 = nlay - 1

  lmt = kord - 3

  ! ****6***0*********0*********0*********0*********0*********0**********72
  ! Compute DC for PPM
  ! ****6***0*********0*********0*********0*********0*********0**********72

  DO k = 1, nlaym1
    DO i = 1, imjm
      dqdt(i, 1, k) = p(i, 1, k+1) - p(i, 1, k)
    END DO
  END DO

  DO k = 2, nlaym1
    DO i = 1, imjm
      c0 = delp(i, 1, k)/(delp(i,1,k-1)+delp(i,1,k)+delp(i,1,k+1))
      c1 = (delp(i,1,k-1)+0.5*delp(i,1,k))/(delp(i,1,k+1)+delp(i,1,k))
      c2 = (delp(i,1,k+1)+0.5*delp(i,1,k))/(delp(i,1,k-1)+delp(i,1,k))
      tmp = c0*(c1*dqdt(i,1,k)+c2*dqdt(i,1,k-1))
      qmax = max(p(i,1,k-1), p(i,1,k), p(i,1,k+1)) - p(i, 1, k)
      qmin = p(i, 1, k) - min(p(i,1,k-1), p(i,1,k), p(i,1,k+1))
      dc(i, 1, k) = sign(min(abs(tmp),qmax,qmin), tmp)
    END DO
  END DO


  ! ****6***0*********0*********0*********0*********0*********0**********72
  ! Loop over latitudes  (to save memory)
  ! ****6***0*********0*********0*********0*********0*********0**********72

  DO j = 1, jnp
    IF ((j==2 .OR. j==jmr) .AND. j1/=2) GO TO 2000

    DO k = 1, nlay
      DO i = 1, imr
        wz2(i, k) = wz(i, j, k)
        wk1(i, k) = p(i, j, k)
        wk2(i, k) = delp(i, j, k)
        flux(i, k) = dc(i, j, k) !this flux is actually the monotone slope
      END DO
    END DO

    ! ****6***0*********0*********0*********0*********0*********0**********72
    ! Compute first guesses at cell interfaces
    ! First guesses are required to be continuous.
    ! ****6***0*********0*********0*********0*********0*********0**********72

    ! three-cell parabolic subgrid distribution at model top
    ! two-cell parabolic with zero gradient subgrid distribution
    ! at the surface.

    ! First guess top edge value
    DO i = 1, imr
      ! three-cell PPM
      ! Compute a,b, and c of q = aP**2 + bP + c using cell averages and delp
      a = 3.*(dqdt(i,j,2)-dqdt(i,j,1)*(wk2(i,2)+wk2(i,3))/(wk2(i,1)+wk2(i, &
        2)))/((wk2(i,2)+wk2(i,3))*(wk2(i,1)+wk2(i,2)+wk2(i,3)))
      b = 2.*dqdt(i, j, 1)/(wk2(i,1)+wk2(i,2)) - r23*a*(2.*wk2(i,1)+wk2(i,2))
      al(i, 1) = wk1(i, 1) - wk2(i, 1)*(r3*a*wk2(i,1)+0.5*b)
      al(i, 2) = wk2(i, 1)*(a*wk2(i,1)+b) + al(i, 1)

      ! Check if change sign
      IF (wk1(i,1)*al(i,1)<=0.) THEN
        al(i, 1) = 0.
        flux(i, 1) = 0.
      ELSE
        flux(i, 1) = wk1(i, 1) - al(i, 1)
      END IF
    END DO

    ! Bottom
    DO i = 1, imr
      ! 2-cell PPM with zero gradient right at the surface

      fct = dqdt(i, j, nlaym1)*wk2(i, nlay)**2/((wk2(i,nlay)+wk2(i, &
        nlaym1))*(2.*wk2(i,nlay)+wk2(i,nlaym1)))
      ar(i, nlay) = wk1(i, nlay) + fct
      al(i, nlay) = wk1(i, nlay) - (fct+fct)
      IF (wk1(i,nlay)*ar(i,nlay)<=0.) ar(i, nlay) = 0.
      flux(i, nlay) = ar(i, nlay) - wk1(i, nlay)
    END DO


    ! ****6***0*********0*********0*********0*********0*********0**********72
    ! 4th order interpolation in the interior.
    ! ****6***0*********0*********0*********0*********0*********0**********72

    DO k = 3, nlaym1
      DO i = 1, imr
        c1 = dqdt(i, j, k-1)*wk2(i, k-1)/(wk2(i,k-1)+wk2(i,k))
        c2 = 2./(wk2(i,k-2)+wk2(i,k-1)+wk2(i,k)+wk2(i,k+1))
        a1 = (wk2(i,k-2)+wk2(i,k-1))/(2.*wk2(i,k-1)+wk2(i,k))
        a2 = (wk2(i,k)+wk2(i,k+1))/(2.*wk2(i,k)+wk2(i,k-1))
        al(i, k) = wk1(i, k-1) + c1 + c2*(wk2(i,k)*(c1*(a1-a2)+a2*flux(i, &
          k-1))-wk2(i,k-1)*a1*flux(i,k))
      END DO
    END DO

    DO i = 1, imr*nlaym1
      ar(i, 1) = al(i, 2)
    END DO

    DO i = 1, imr*nlay
      a6(i, 1) = 3.*(wk1(i,1)+wk1(i,1)-(al(i,1)+ar(i,1)))
    END DO

    ! ****6***0*********0*********0*********0*********0*********0**********72
    ! Top & Bot always monotonic
    CALL lmtppm(flux(1,1), a6(1,1), ar(1,1), al(1,1), wk1(1,1), imr, 0)
    CALL lmtppm(flux(1,nlay), a6(1,nlay), ar(1,nlay), al(1,nlay), &
      wk1(1,nlay), imr, 0)

    ! Interior depending on KORD
    IF (lmt<=2) CALL lmtppm(flux(1,2), a6(1,2), ar(1,2), al(1,2), wk1(1,2), &
      imr*(nlay-2), lmt)

    ! ****6***0*********0*********0*********0*********0*********0**********72

    DO i = 1, imr*nlaym1
      IF (wz2(i,1)>0.) THEN
        cm = wz2(i, 1)/wk2(i, 1)
        flux(i, 2) = ar(i, 1) + 0.5*cm*(al(i,1)-ar(i,1)+a6(i,1)*(1.-r23*cm))
      ELSE
        cp = wz2(i, 1)/wk2(i, 2)
        flux(i, 2) = al(i, 2) + 0.5*cp*(al(i,2)-ar(i,2)-a6(i,2)*(1.+r23*cp))
      END IF
    END DO

    DO i = 1, imr*nlaym1
      flux(i, 2) = wz2(i, 1)*flux(i, 2)
    END DO

    DO i = 1, imr
      dq(i, j, 1) = dq(i, j, 1) - flux(i, 2)
      dq(i, j, nlay) = dq(i, j, nlay) + flux(i, nlay)
    END DO

    DO k = 2, nlaym1
      DO i = 1, imr
        dq(i, j, k) = dq(i, j, k) + flux(i, k) - flux(i, k+1)
      END DO
    END DO
2000 END DO
  RETURN
END SUBROUTINE fzppm

SUBROUTINE xtp(imr, jnp, iml, j1, j2, jn, js, pu, dq, q, uc, fx1, xmass, &
    iord)
  DIMENSION uc(imr, *), dc(-iml:imr+iml+1), xmass(imr, jnp), fx1(imr+1), &
    dq(imr, jnp), qtmp(-iml:imr+1+iml)
  DIMENSION pu(imr, jnp), q(imr, jnp), isave(imr)

  imp = imr + 1

  ! van Leer at high latitudes
  jvan = max(1, jnp/18)
  j1vl = j1 + jvan
  j2vl = j2 - jvan

  DO j = j1, j2

    DO i = 1, imr
      qtmp(i) = q(i, j)
    END DO

    IF (j>=jn .OR. j<=js) GO TO 2222
    ! ************* Eulerian **********

    qtmp(0) = q(imr, j)
    qtmp(-1) = q(imr-1, j)
    qtmp(imp) = q(1, j)
    qtmp(imp+1) = q(2, j)

    IF (iord==1 .OR. j==j1 .OR. j==j2) THEN
      DO i = 1, imr
        iu = float(i) - uc(i, j)
        fx1(i) = qtmp(iu)
      END DO
    ELSE
      CALL xmist(imr, iml, qtmp, dc)
      dc(0) = dc(imr)

      IF (iord==2 .OR. j<=j1vl .OR. j>=j2vl) THEN
        DO i = 1, imr
          iu = float(i) - uc(i, j)
          fx1(i) = qtmp(iu) + dc(iu)*(sign(1.,uc(i,j))-uc(i,j))
        END DO
      ELSE
        CALL fxppm(imr, iml, uc(1,j), qtmp, dc, fx1, iord)
      END IF

    END IF

    DO i = 1, imr
      fx1(i) = fx1(i)*xmass(i, j)
    END DO

    GO TO 1309

    ! ***** Conservative (flux-form) Semi-Lagrangian transport *****

2222 CONTINUE

    DO i = -iml, 0
      qtmp(i) = q(imr+i, j)
      qtmp(imp-i) = q(1-i, j)
    END DO

    IF (iord==1 .OR. j==j1 .OR. j==j2) THEN
      DO i = 1, imr
        itmp = int(uc(i,j))
        isave(i) = i - itmp
        iu = i - uc(i, j)
        fx1(i) = (uc(i,j)-itmp)*qtmp(iu)
      END DO
    ELSE
      CALL xmist(imr, iml, qtmp, dc)

      DO i = -iml, 0
        dc(i) = dc(imr+i)
        dc(imp-i) = dc(1-i)
      END DO

      DO i = 1, imr
        itmp = int(uc(i,j))
        rut = uc(i, j) - itmp
        isave(i) = i - itmp
        iu = i - uc(i, j)
        fx1(i) = rut*(qtmp(iu)+dc(iu)*(sign(1.,rut)-rut))
      END DO
    END IF

    DO i = 1, imr
      IF (uc(i,j)>1.) THEN
        ! DIR$ NOVECTOR
        DO ist = isave(i), i - 1
          fx1(i) = fx1(i) + qtmp(ist)
        END DO
      ELSE IF (uc(i,j)<-1.) THEN
        DO ist = i, isave(i) - 1
          fx1(i) = fx1(i) - qtmp(ist)
        END DO
        ! DIR$ VECTOR
      END IF
    END DO
    DO i = 1, imr
      fx1(i) = pu(i, j)*fx1(i)
    END DO

    ! ***************************************

1309 fx1(imp) = fx1(1)
    DO i = 1, imr
      dq(i, j) = dq(i, j) + fx1(i) - fx1(i+1)
    END DO

    ! ***************************************

  END DO
  RETURN
END SUBROUTINE xtp

SUBROUTINE fxppm(imr, iml, ut, p, dc, flux, iord)
  PARAMETER (r3=1./3., r23=2./3.)
  DIMENSION ut(*), flux(*), p(-iml:imr+iml+1), dc(-iml:imr+iml+1)
  DIMENSION ar(0:imr), al(0:imr), a6(0:imr)
  INTEGER lmt
  ! logical first
  ! data first /.true./
  ! SAVE LMT
  ! if(first) then

  ! correction calcul de LMT a chaque passage pour pouvoir choisir
  ! plusieurs schemas PPM pour differents traceurs
  ! IF (IORD.LE.0) then
  ! if(IMR.GE.144) then
  ! LMT = 0
  ! elseif(IMR.GE.72) then
  ! LMT = 1
  ! else
  ! LMT = 2
  ! endif
  ! else
  ! LMT = IORD - 3
  ! endif

  lmt = iord - 3

  DO i = 1, imr
    al(i) = 0.5*(p(i-1)+p(i)) + (dc(i-1)-dc(i))*r3
  END DO

  DO i = 1, imr - 1
    ar(i) = al(i+1)
  END DO
  ar(imr) = al(1)

  DO i = 1, imr
    a6(i) = 3.*(p(i)+p(i)-(al(i)+ar(i)))
  END DO

  IF (lmt<=2) CALL lmtppm(dc(1), a6(1), ar(1), al(1), p(1), imr, lmt)

  al(0) = al(imr)
  ar(0) = ar(imr)
  a6(0) = a6(imr)

  DO i = 1, imr
    IF (ut(i)>0.) THEN
      flux(i) = ar(i-1) + 0.5*ut(i)*(al(i-1)-ar(i-1)+a6(i-1)*(1.-r23*ut(i)))
    ELSE
      flux(i) = al(i) - 0.5*ut(i)*(ar(i)-al(i)+a6(i)*(1.+r23*ut(i)))
    END IF
  END DO
  RETURN
END SUBROUTINE fxppm

SUBROUTINE xmist(imr, iml, p, dc)
  PARAMETER (r24=1./24.)
  DIMENSION p(-iml:imr+1+iml), dc(-iml:imr+1+iml)

  DO i = 1, imr
    tmp = r24*(8.*(p(i+1)-p(i-1))+p(i-2)-p(i+2))
    pmax = max(p(i-1), p(i), p(i+1)) - p(i)
    pmin = p(i) - min(p(i-1), p(i), p(i+1))
    dc(i) = sign(min(abs(tmp),pmax,pmin), tmp)
  END DO
  RETURN
END SUBROUTINE xmist

SUBROUTINE ytp(imr, jnp, j1, j2, acosp, rcap, dq, p, vc, dc2, ymass, fx, a6, &
    ar, al, jord)
  DIMENSION p(imr, jnp), vc(imr, jnp), ymass(imr, jnp), dc2(imr, jnp), &
    dq(imr, jnp), acosp(jnp)
  ! Work array
  DIMENSION fx(imr, jnp), ar(imr, jnp), al(imr, jnp), a6(imr, jnp)

  jmr = jnp - 1
  len = imr*(j2-j1+2)

  IF (jord==1) THEN
    DO i = 1, len
      jt = float(j1) - vc(i, j1)
      fx(i, j1) = p(i, jt)
    END DO
  ELSE

    CALL ymist(imr, jnp, j1, p, dc2, 4)

    IF (jord<=0 .OR. jord>=3) THEN

      CALL fyppm(vc, p, dc2, fx, imr, jnp, j1, j2, a6, ar, al, jord)

    ELSE
      DO i = 1, len
        jt = float(j1) - vc(i, j1)
        fx(i, j1) = p(i, jt) + (sign(1.,vc(i,j1))-vc(i,j1))*dc2(i, jt)
      END DO
    END IF
  END IF

  DO i = 1, len
    fx(i, j1) = fx(i, j1)*ymass(i, j1)
  END DO

  DO j = j1, j2
    DO i = 1, imr
      dq(i, j) = dq(i, j) + (fx(i,j)-fx(i,j+1))*acosp(j)
    END DO
  END DO

  ! Poles
  sum1 = fx(imr, j1)
  sum2 = fx(imr, j2+1)
  DO i = 1, imr - 1
    sum1 = sum1 + fx(i, j1)
    sum2 = sum2 + fx(i, j2+1)
  END DO

  sum1 = dq(1, 1) - sum1*rcap
  sum2 = dq(1, jnp) + sum2*rcap
  DO i = 1, imr
    dq(i, 1) = sum1
    dq(i, jnp) = sum2
  END DO

  IF (j1/=2) THEN
    DO i = 1, imr
      dq(i, 2) = sum1
      dq(i, jmr) = sum2
    END DO
  END IF

  RETURN
END SUBROUTINE ytp

SUBROUTINE ymist(imr, jnp, j1, p, dc, id)
  PARAMETER (r24=1./24.)
  DIMENSION p(imr, jnp), dc(imr, jnp)

  imh = imr/2
  jmr = jnp - 1
  ijm3 = imr*(jmr-3)

  IF (id==2) THEN
    DO i = 1, imr*(jmr-1)
      tmp = 0.25*(p(i,3)-p(i,1))
      pmax = max(p(i,1), p(i,2), p(i,3)) - p(i, 2)
      pmin = p(i, 2) - min(p(i,1), p(i,2), p(i,3))
      dc(i, 2) = sign(min(abs(tmp),pmin,pmax), tmp)
    END DO
  ELSE
    DO i = 1, imh
      ! J=2
      tmp = (8.*(p(i,3)-p(i,1))+p(i+imh,2)-p(i,4))*r24
      pmax = max(p(i,1), p(i,2), p(i,3)) - p(i, 2)
      pmin = p(i, 2) - min(p(i,1), p(i,2), p(i,3))
      dc(i, 2) = sign(min(abs(tmp),pmin,pmax), tmp)
      ! J=JMR
      tmp = (8.*(p(i,jnp)-p(i,jmr-1))+p(i,jmr-2)-p(i+imh,jmr))*r24
      pmax = max(p(i,jmr-1), p(i,jmr), p(i,jnp)) - p(i, jmr)
      pmin = p(i, jmr) - min(p(i,jmr-1), p(i,jmr), p(i,jnp))
      dc(i, jmr) = sign(min(abs(tmp),pmin,pmax), tmp)
    END DO
    DO i = imh + 1, imr
      ! J=2
      tmp = (8.*(p(i,3)-p(i,1))+p(i-imh,2)-p(i,4))*r24
      pmax = max(p(i,1), p(i,2), p(i,3)) - p(i, 2)
      pmin = p(i, 2) - min(p(i,1), p(i,2), p(i,3))
      dc(i, 2) = sign(min(abs(tmp),pmin,pmax), tmp)
      ! J=JMR
      tmp = (8.*(p(i,jnp)-p(i,jmr-1))+p(i,jmr-2)-p(i-imh,jmr))*r24
      pmax = max(p(i,jmr-1), p(i,jmr), p(i,jnp)) - p(i, jmr)
      pmin = p(i, jmr) - min(p(i,jmr-1), p(i,jmr), p(i,jnp))
      dc(i, jmr) = sign(min(abs(tmp),pmin,pmax), tmp)
    END DO

    DO i = 1, ijm3
      tmp = (8.*(p(i,4)-p(i,2))+p(i,1)-p(i,5))*r24
      pmax = max(p(i,2), p(i,3), p(i,4)) - p(i, 3)
      pmin = p(i, 3) - min(p(i,2), p(i,3), p(i,4))
      dc(i, 3) = sign(min(abs(tmp),pmin,pmax), tmp)
    END DO
  END IF

  IF (j1/=2) THEN
    DO i = 1, imr
      dc(i, 1) = 0.
      dc(i, jnp) = 0.
    END DO
  ELSE
    ! Determine slopes in polar caps for scalars!

    DO i = 1, imh
      ! South
      tmp = 0.25*(p(i,2)-p(i+imh,2))
      pmax = max(p(i,2), p(i,1), p(i+imh,2)) - p(i, 1)
      pmin = p(i, 1) - min(p(i,2), p(i,1), p(i+imh,2))
      dc(i, 1) = sign(min(abs(tmp),pmax,pmin), tmp)
      ! North.
      tmp = 0.25*(p(i+imh,jmr)-p(i,jmr))
      pmax = max(p(i+imh,jmr), p(i,jnp), p(i,jmr)) - p(i, jnp)
      pmin = p(i, jnp) - min(p(i+imh,jmr), p(i,jnp), p(i,jmr))
      dc(i, jnp) = sign(min(abs(tmp),pmax,pmin), tmp)
    END DO

    DO i = imh + 1, imr
      dc(i, 1) = -dc(i-imh, 1)
      dc(i, jnp) = -dc(i-imh, jnp)
    END DO
  END IF
  RETURN
END SUBROUTINE ymist

SUBROUTINE fyppm(vc, p, dc, flux, imr, jnp, j1, j2, a6, ar, al, jord)
  PARAMETER (r3=1./3., r23=2./3.)
  REAL vc(imr, *), flux(imr, *), p(imr, *), dc(imr, *)
  ! Local work arrays.
  REAL ar(imr, jnp), al(imr, jnp), a6(imr, jnp)
  INTEGER lmt
  ! logical first
  ! data first /.true./
  ! SAVE LMT

  imh = imr/2
  jmr = jnp - 1
  j11 = j1 - 1
  imjm1 = imr*(j2-j1+2)
  len = imr*(j2-j1+3)
  ! if(first) then
  ! IF(JORD.LE.0) then
  ! if(JMR.GE.90) then
  ! LMT = 0
  ! elseif(JMR.GE.45) then
  ! LMT = 1
  ! else
  ! LMT = 2
  ! endif
  ! else
  ! LMT = JORD - 3
  ! endif

  ! first = .false.
  ! endif

  ! modifs pour pouvoir choisir plusieurs schemas PPM
  lmt = jord - 3

  DO i = 1, imr*jmr
    al(i, 2) = 0.5*(p(i,1)+p(i,2)) + (dc(i,1)-dc(i,2))*r3
    ar(i, 1) = al(i, 2)
  END DO

  ! Poles:

  DO i = 1, imh
    al(i, 1) = al(i+imh, 2)
    al(i+imh, 1) = al(i, 2)

    ar(i, jnp) = ar(i+imh, jmr)
    ar(i+imh, jnp) = ar(i, jmr)
  END DO

  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Rajout pour LMDZ.3.3
  ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ar(imr, 1) = al(1, 1)
  ar(imr, jnp) = al(1, jnp)
  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  DO i = 1, len
    a6(i, j11) = 3.*(p(i,j11)+p(i,j11)-(al(i,j11)+ar(i,j11)))
  END DO

  IF (lmt<=2) CALL lmtppm(dc(1,j11), a6(1,j11), ar(1,j11), al(1,j11), &
    p(1,j11), len, lmt)


  DO i = 1, imjm1
    IF (vc(i,j1)>0.) THEN
      flux(i, j1) = ar(i, j11) + 0.5*vc(i, j1)*(al(i,j11)-ar(i,j11)+a6(i,j11) &
        *(1.-r23*vc(i,j1)))
    ELSE
      flux(i, j1) = al(i, j1) - 0.5*vc(i, j1)*(ar(i,j1)-al(i,j1)+a6(i,j1)*(1. &
        +r23*vc(i,j1)))
    END IF
  END DO
  RETURN
END SUBROUTINE fyppm

SUBROUTINE yadv(imr, jnp, j1, j2, p, va, ady, wk, iad)
  REAL p(imr, jnp), ady(imr, jnp), va(imr, jnp)
  REAL wk(imr, -1:jnp+2)

  jmr = jnp - 1
  imh = imr/2
  DO j = 1, jnp
    DO i = 1, imr
      wk(i, j) = p(i, j)
    END DO
  END DO
  ! Poles:
  DO i = 1, imh
    wk(i, -1) = p(i+imh, 3)
    wk(i+imh, -1) = p(i, 3)
    wk(i, 0) = p(i+imh, 2)
    wk(i+imh, 0) = p(i, 2)
    wk(i, jnp+1) = p(i+imh, jmr)
    wk(i+imh, jnp+1) = p(i, jmr)
    wk(i, jnp+2) = p(i+imh, jnp-2)
    wk(i+imh, jnp+2) = p(i, jnp-2)
  END DO

  IF (iad==2) THEN
    DO j = j1 - 1, j2 + 1
      DO i = 1, imr
        jp = nint(va(i,j))
        rv = jp - va(i, j)
        jp = j - jp
        a1 = 0.5*(wk(i,jp+1)+wk(i,jp-1)) - wk(i, jp)
        b1 = 0.5*(wk(i,jp+1)-wk(i,jp-1))
        ady(i, j) = wk(i, jp) + rv*(a1*rv+b1) - wk(i, j)
      END DO
    END DO

  ELSE IF (iad==1) THEN
    DO j = j1 - 1, j2 + 1
      DO i = 1, imr
        jp = float(j) - va(i, j)
        ady(i, j) = va(i, j)*(wk(i,jp)-wk(i,jp+1))
      END DO
    END DO
  END IF

  IF (j1/=2) THEN
    sum1 = 0.
    sum2 = 0.
    DO i = 1, imr
      sum1 = sum1 + ady(i, 2)
      sum2 = sum2 + ady(i, jmr)
    END DO
    sum1 = sum1/imr
    sum2 = sum2/imr

    DO i = 1, imr
      ady(i, 2) = sum1
      ady(i, jmr) = sum2
      ady(i, 1) = sum1
      ady(i, jnp) = sum2
    END DO
  ELSE
    ! Poles:
    sum1 = 0.
    sum2 = 0.
    DO i = 1, imr
      sum1 = sum1 + ady(i, 1)
      sum2 = sum2 + ady(i, jnp)
    END DO
    sum1 = sum1/imr
    sum2 = sum2/imr

    DO i = 1, imr
      ady(i, 1) = sum1
      ady(i, jnp) = sum2
    END DO
  END IF

  RETURN
END SUBROUTINE yadv

SUBROUTINE xadv(imr, jnp, j1, j2, p, ua, js, jn, iml, adx, iad)
  REAL p(imr, jnp), adx(imr, jnp), qtmp(-imr:imr+imr), ua(imr, jnp)

  jmr = jnp - 1
  DO j = j1, j2
    IF (j>js .AND. j<jn) GO TO 1309

    DO i = 1, imr
      qtmp(i) = p(i, j)
    END DO

    DO i = -iml, 0
      qtmp(i) = p(imr+i, j)
      qtmp(imr+1-i) = p(1-i, j)
    END DO

    IF (iad==2) THEN
      DO i = 1, imr
        ip = nint(ua(i,j))
        ru = ip - ua(i, j)
        ip = i - ip
        a1 = 0.5*(qtmp(ip+1)+qtmp(ip-1)) - qtmp(ip)
        b1 = 0.5*(qtmp(ip+1)-qtmp(ip-1))
        adx(i, j) = qtmp(ip) + ru*(a1*ru+b1)
      END DO
    ELSE IF (iad==1) THEN
      DO i = 1, imr
        iu = ua(i, j)
        ru = ua(i, j) - iu
        iiu = i - iu
        IF (ua(i,j)>=0.) THEN
          adx(i, j) = qtmp(iiu) + ru*(qtmp(iiu-1)-qtmp(iiu))
        ELSE
          adx(i, j) = qtmp(iiu) + ru*(qtmp(iiu)-qtmp(iiu+1))
        END IF
      END DO
    END IF

    DO i = 1, imr
      adx(i, j) = adx(i, j) - p(i, j)
    END DO
1309 END DO

  ! Eulerian upwind

  DO j = js + 1, jn - 1

    DO i = 1, imr
      qtmp(i) = p(i, j)
    END DO

    qtmp(0) = p(imr, j)
    qtmp(imr+1) = p(1, j)

    IF (iad==2) THEN
      qtmp(-1) = p(imr-1, j)
      qtmp(imr+2) = p(2, j)
      DO i = 1, imr
        ip = nint(ua(i,j))
        ru = ip - ua(i, j)
        ip = i - ip
        a1 = 0.5*(qtmp(ip+1)+qtmp(ip-1)) - qtmp(ip)
        b1 = 0.5*(qtmp(ip+1)-qtmp(ip-1))
        adx(i, j) = qtmp(ip) - p(i, j) + ru*(a1*ru+b1)
      END DO
    ELSE IF (iad==1) THEN
      ! 1st order
      DO i = 1, imr
        ip = i - ua(i, j)
        adx(i, j) = ua(i, j)*(qtmp(ip)-qtmp(ip+1))
      END DO
    END IF
  END DO

  IF (j1/=2) THEN
    DO i = 1, imr
      adx(i, 2) = 0.
      adx(i, jmr) = 0.
    END DO
  END IF
  ! set cross term due to x-adv at the poles to zero.
  DO i = 1, imr
    adx(i, 1) = 0.
    adx(i, jnp) = 0.
  END DO
  RETURN
END SUBROUTINE xadv

SUBROUTINE lmtppm(dc, a6, ar, al, p, im, lmt)

  ! A6 =  CURVATURE OF THE TEST PARABOLA
  ! AR =  RIGHT EDGE VALUE OF THE TEST PARABOLA
  ! AL =  LEFT  EDGE VALUE OF THE TEST PARABOLA
  ! DC =  0.5 * MISMATCH
  ! P  =  CELL-AVERAGED VALUE
  ! IM =  VECTOR LENGTH

  ! OPTIONS:

  ! LMT = 0: FULL MONOTONICITY
  ! LMT = 1: SEMI-MONOTONIC CONSTRAINT (NO UNDERSHOOTS)
  ! LMT = 2: POSITIVE-DEFINITE CONSTRAINT

  PARAMETER (r12=1./12.)
  DIMENSION a6(im), ar(im), al(im), p(im), dc(im)

  IF (lmt==0) THEN
    ! Full constraint
    DO i = 1, im
      IF (dc(i)==0.) THEN
        ar(i) = p(i)
        al(i) = p(i)
        a6(i) = 0.
      ELSE
        da1 = ar(i) - al(i)
        da2 = da1**2
        a6da = a6(i)*da1
        IF (a6da<-da2) THEN
          a6(i) = 3.*(al(i)-p(i))
          ar(i) = al(i) - a6(i)
        ELSE IF (a6da>da2) THEN
          a6(i) = 3.*(ar(i)-p(i))
          al(i) = ar(i) - a6(i)
        END IF
      END IF
    END DO
  ELSE IF (lmt==1) THEN
    ! Semi-monotonic constraint
    DO i = 1, im
      IF (abs(ar(i)-al(i))>=-a6(i)) GO TO 150
      IF (p(i)<ar(i) .AND. p(i)<al(i)) THEN
        ar(i) = p(i)
        al(i) = p(i)
        a6(i) = 0.
      ELSE IF (ar(i)>al(i)) THEN
        a6(i) = 3.*(al(i)-p(i))
        ar(i) = al(i) - a6(i)
      ELSE
        a6(i) = 3.*(ar(i)-p(i))
        al(i) = ar(i) - a6(i)
      END IF
150 END DO
  ELSE IF (lmt==2) THEN
    DO i = 1, im
      IF (abs(ar(i)-al(i))>=-a6(i)) GO TO 250
      fmin = p(i) + 0.25*(ar(i)-al(i))**2/a6(i) + a6(i)*r12
      IF (fmin>=0.) GO TO 250
      IF (p(i)<ar(i) .AND. p(i)<al(i)) THEN
        ar(i) = p(i)
        al(i) = p(i)
        a6(i) = 0.
      ELSE IF (ar(i)>al(i)) THEN
        a6(i) = 3.*(al(i)-p(i))
        ar(i) = al(i) - a6(i)
      ELSE
        a6(i) = 3.*(ar(i)-p(i))
        al(i) = ar(i) - a6(i)
      END IF
250 END DO
  END IF
  RETURN
END SUBROUTINE lmtppm

SUBROUTINE a2c(u, v, imr, jmr, j1, j2, crx, cry, dtdx5, dtdy5)
  DIMENSION u(imr, *), v(imr, *), crx(imr, *), cry(imr, *), dtdx5(*)

  DO j = j1, j2
    DO i = 2, imr
      crx(i, j) = dtdx5(j)*(u(i,j)+u(i-1,j))
    END DO
  END DO

  DO j = j1, j2
    crx(1, j) = dtdx5(j)*(u(1,j)+u(imr,j))
  END DO

  DO i = 1, imr*jmr
    cry(i, 2) = dtdy5*(v(i,2)+v(i,1))
  END DO
  RETURN
END SUBROUTINE a2c

SUBROUTINE cosa(cosp, cose, jnp, pi, dp)
  DIMENSION cosp(*), cose(*)

  jmr = jnp - 1
  DO j = 2, jnp
    ph5 = -0.5*pi + (float(j-1)-0.5)*dp
    cose(j) = cos(ph5)
  END DO

  jeq = (jnp+1)/2
  IF (jmr==2*(jmr/2)) THEN
    DO j = jnp, jeq + 1, -1
      cose(j) = cose(jnp+2-j)
    END DO
  ELSE
    ! cell edge at equator.
    cose(jeq+1) = 1.
    DO j = jnp, jeq + 2, -1
      cose(j) = cose(jnp+2-j)
    END DO
  END IF

  DO j = 2, jmr
    cosp(j) = 0.5*(cose(j)+cose(j+1))
  END DO
  cosp(1) = 0.
  cosp(jnp) = 0.
  RETURN
END SUBROUTINE cosa

SUBROUTINE cosc(cosp, cose, jnp, pi, dp)
  DIMENSION cosp(*), cose(*)

  phi = -0.5*pi
  DO j = 2, jnp - 1
    phi = phi + dp
    cosp(j) = cos(phi)
  END DO
  cosp(1) = 0.
  cosp(jnp) = 0.

  DO j = 2, jnp
    cose(j) = 0.5*(cosp(j)+cosp(j-1))
  END DO

  DO j = 2, jnp - 1
    cosp(j) = 0.5*(cose(j)+cose(j+1))
  END DO
  RETURN
END SUBROUTINE cosc

SUBROUTINE qckxyz(q, qtmp, imr, jnp, nlay, j1, j2, cosp, acosp, cross, ic, &
    nstep)

  PARAMETER (tiny=1.E-60)
  DIMENSION q(imr, jnp, nlay), qtmp(imr, jnp), cosp(*), acosp(*)
  LOGICAL cross

  nlaym1 = nlay - 1
  len = imr*(j2-j1+1)
  ip = 0

  ! Top layer
  l = 1
  icr = 1
  CALL filns(q(1,1,l), imr, jnp, j1, j2, cosp, acosp, ipy, tiny)
  IF (ipy==0) GO TO 50
  CALL filew(q(1,1,l), qtmp, imr, jnp, j1, j2, ipx, tiny)
  IF (ipx==0) GO TO 50

  IF (cross) THEN
    CALL filcr(q(1,1,l), imr, jnp, j1, j2, cosp, acosp, icr, tiny)
  END IF
  IF (icr==0) GO TO 50

  ! Vertical filling...
  DO i = 1, len
    IF (q(i,j1,1)<0.) THEN
      ip = ip + 1
      q(i, j1, 2) = q(i, j1, 2) + q(i, j1, 1)
      q(i, j1, 1) = 0.
    END IF
  END DO

50 CONTINUE
  DO l = 2, nlaym1
    icr = 1

    CALL filns(q(1,1,l), imr, jnp, j1, j2, cosp, acosp, ipy, tiny)
    IF (ipy==0) GO TO 225
    CALL filew(q(1,1,l), qtmp, imr, jnp, j1, j2, ipx, tiny)
    IF (ipx==0) GO TO 225
    IF (cross) THEN
      CALL filcr(q(1,1,l), imr, jnp, j1, j2, cosp, acosp, icr, tiny)
    END IF
    IF (icr==0) GO TO 225

    DO i = 1, len
      IF (q(i,j1,l)<0.) THEN

        ip = ip + 1
        ! From above
        qup = q(i, j1, l-1)
        qly = -q(i, j1, l)
        dup = min(qly, qup)
        q(i, j1, l-1) = qup - dup
        q(i, j1, l) = dup - qly
        ! Below
        q(i, j1, l+1) = q(i, j1, l+1) + q(i, j1, l)
        q(i, j1, l) = 0.
      END IF
    END DO
225 END DO

  ! BOTTOM LAYER
  sum = 0.
  l = nlay

  CALL filns(q(1,1,l), imr, jnp, j1, j2, cosp, acosp, ipy, tiny)
  IF (ipy==0) GO TO 911
  CALL filew(q(1,1,l), qtmp, imr, jnp, j1, j2, ipx, tiny)
  IF (ipx==0) GO TO 911

  CALL filcr(q(1,1,l), imr, jnp, j1, j2, cosp, acosp, icr, tiny)
  IF (icr==0) GO TO 911

  DO i = 1, len
    IF (q(i,j1,l)<0.) THEN
      ip = ip + 1

      ! From above

      qup = q(i, j1, nlaym1)
      qly = -q(i, j1, l)
      dup = min(qly, qup)
      q(i, j1, nlaym1) = qup - dup
      ! From "below" the surface.
      sum = sum + qly - dup
      q(i, j1, l) = 0.
    END IF
  END DO

911 CONTINUE

  IF (ip>imr) THEN
    WRITE (6, *) 'IC=', ic, ' STEP=', nstep, ' Vertical filling pts=', ip
  END IF

  IF (sum>1.E-25) THEN
    WRITE (6, *) ic, nstep, ' Mass source from the ground=', sum
  END IF
  RETURN
END SUBROUTINE qckxyz

SUBROUTINE filcr(q, imr, jnp, j1, j2, cosp, acosp, icr, tiny)
  DIMENSION q(imr, *), cosp(*), acosp(*)

  icr = 0
  DO j = j1 + 1, j2 - 1
    DO i = 1, imr - 1
      IF (q(i,j)<0.) THEN
        icr = 1
        dq = -q(i, j)*cosp(j)
        ! N-E
        dn = q(i+1, j+1)*cosp(j+1)
        d0 = max(0., dn)
        d1 = min(dq, d0)
        q(i+1, j+1) = (dn-d1)*acosp(j+1)
        dq = dq - d1
        ! S-E
        ds = q(i+1, j-1)*cosp(j-1)
        d0 = max(0., ds)
        d2 = min(dq, d0)
        q(i+1, j-1) = (ds-d2)*acosp(j-1)
        q(i, j) = (d2-dq)*acosp(j) + tiny
      END IF
    END DO
    IF (icr==0 .AND. q(imr,j)>=0.) GO TO 65
    DO i = 2, imr
      IF (q(i,j)<0.) THEN
        icr = 1
        dq = -q(i, j)*cosp(j)
        ! N-W
        dn = q(i-1, j+1)*cosp(j+1)
        d0 = max(0., dn)
        d1 = min(dq, d0)
        q(i-1, j+1) = (dn-d1)*acosp(j+1)
        dq = dq - d1
        ! S-W
        ds = q(i-1, j-1)*cosp(j-1)
        d0 = max(0., ds)
        d2 = min(dq, d0)
        q(i-1, j-1) = (ds-d2)*acosp(j-1)
        q(i, j) = (d2-dq)*acosp(j) + tiny
      END IF
    END DO
    ! *****************************************
    ! i=1
    i = 1
    IF (q(i,j)<0.) THEN
      icr = 1
      dq = -q(i, j)*cosp(j)
      ! N-W
      dn = q(imr, j+1)*cosp(j+1)
      d0 = max(0., dn)
      d1 = min(dq, d0)
      q(imr, j+1) = (dn-d1)*acosp(j+1)
      dq = dq - d1
      ! S-W
      ds = q(imr, j-1)*cosp(j-1)
      d0 = max(0., ds)
      d2 = min(dq, d0)
      q(imr, j-1) = (ds-d2)*acosp(j-1)
      q(i, j) = (d2-dq)*acosp(j) + tiny
    END IF
    ! *****************************************
    ! i=IMR
    i = imr
    IF (q(i,j)<0.) THEN
      icr = 1
      dq = -q(i, j)*cosp(j)
      ! N-E
      dn = q(1, j+1)*cosp(j+1)
      d0 = max(0., dn)
      d1 = min(dq, d0)
      q(1, j+1) = (dn-d1)*acosp(j+1)
      dq = dq - d1
      ! S-E
      ds = q(1, j-1)*cosp(j-1)
      d0 = max(0., ds)
      d2 = min(dq, d0)
      q(1, j-1) = (ds-d2)*acosp(j-1)
      q(i, j) = (d2-dq)*acosp(j) + tiny
    END IF
    ! *****************************************
65 END DO

  DO i = 1, imr
    IF (q(i,j1)<0. .OR. q(i,j2)<0.) THEN
      icr = 1
      GO TO 80
    END IF
  END DO

80 CONTINUE

  IF (q(1,1)<0. .OR. q(1,jnp)<0.) THEN
    icr = 1
  END IF

  RETURN
END SUBROUTINE filcr

SUBROUTINE filns(q, imr, jnp, j1, j2, cosp, acosp, ipy, tiny)
  DIMENSION q(imr, *), cosp(*), acosp(*)
  ! logical first
  ! data first /.true./
  ! save cap1

  ! if(first) then
  dp = 4.*atan(1.)/float(jnp-1)
  cap1 = imr*(1.-cos((j1-1.5)*dp))/dp
  ! first = .false.
  ! endif

  ipy = 0
  DO j = j1 + 1, j2 - 1
    DO i = 1, imr
      IF (q(i,j)<0.) THEN
        ipy = 1
        dq = -q(i, j)*cosp(j)
        ! North
        dn = q(i, j+1)*cosp(j+1)
        d0 = max(0., dn)
        d1 = min(dq, d0)
        q(i, j+1) = (dn-d1)*acosp(j+1)
        dq = dq - d1
        ! South
        ds = q(i, j-1)*cosp(j-1)
        d0 = max(0., ds)
        d2 = min(dq, d0)
        q(i, j-1) = (ds-d2)*acosp(j-1)
        q(i, j) = (d2-dq)*acosp(j) + tiny
      END IF
    END DO
  END DO

  DO i = 1, imr
    IF (q(i,j1)<0.) THEN
      ipy = 1
      dq = -q(i, j1)*cosp(j1)
      ! North
      dn = q(i, j1+1)*cosp(j1+1)
      d0 = max(0., dn)
      d1 = min(dq, d0)
      q(i, j1+1) = (dn-d1)*acosp(j1+1)
      q(i, j1) = (d1-dq)*acosp(j1) + tiny
    END IF
  END DO

  j = j2
  DO i = 1, imr
    IF (q(i,j)<0.) THEN
      ipy = 1
      dq = -q(i, j)*cosp(j)
      ! South
      ds = q(i, j-1)*cosp(j-1)
      d0 = max(0., ds)
      d2 = min(dq, d0)
      q(i, j-1) = (ds-d2)*acosp(j-1)
      q(i, j) = (d2-dq)*acosp(j) + tiny
    END IF
  END DO

  ! Check Poles.
  IF (q(1,1)<0.) THEN
    dq = q(1, 1)*cap1/float(imr)*acosp(j1)
    DO i = 1, imr
      q(i, 1) = 0.
      q(i, j1) = q(i, j1) + dq
      IF (q(i,j1)<0.) ipy = 1
    END DO
  END IF

  IF (q(1,jnp)<0.) THEN
    dq = q(1, jnp)*cap1/float(imr)*acosp(j2)
    DO i = 1, imr
      q(i, jnp) = 0.
      q(i, j2) = q(i, j2) + dq
      IF (q(i,j2)<0.) ipy = 1
    END DO
  END IF

  RETURN
END SUBROUTINE filns

SUBROUTINE filew(q, qtmp, imr, jnp, j1, j2, ipx, tiny)
  DIMENSION q(imr, *), qtmp(jnp, imr)

  ipx = 0
  ! Copy & swap direction for vectorization.
  DO i = 1, imr
    DO j = j1, j2
      qtmp(j, i) = q(i, j)
    END DO
  END DO

  DO i = 2, imr - 1
    DO j = j1, j2
      IF (qtmp(j,i)<0.) THEN
        ipx = 1
        ! west
        d0 = max(0., qtmp(j,i-1))
        d1 = min(-qtmp(j,i), d0)
        qtmp(j, i-1) = qtmp(j, i-1) - d1
        qtmp(j, i) = qtmp(j, i) + d1
        ! east
        d0 = max(0., qtmp(j,i+1))
        d2 = min(-qtmp(j,i), d0)
        qtmp(j, i+1) = qtmp(j, i+1) - d2
        qtmp(j, i) = qtmp(j, i) + d2 + tiny
      END IF
    END DO
  END DO

  i = 1
  DO j = j1, j2
    IF (qtmp(j,i)<0.) THEN
      ipx = 1
      ! west
      d0 = max(0., qtmp(j,imr))
      d1 = min(-qtmp(j,i), d0)
      qtmp(j, imr) = qtmp(j, imr) - d1
      qtmp(j, i) = qtmp(j, i) + d1
      ! east
      d0 = max(0., qtmp(j,i+1))
      d2 = min(-qtmp(j,i), d0)
      qtmp(j, i+1) = qtmp(j, i+1) - d2

      qtmp(j, i) = qtmp(j, i) + d2 + tiny
    END IF
  END DO
  i = imr
  DO j = j1, j2
    IF (qtmp(j,i)<0.) THEN
      ipx = 1
      ! west
      d0 = max(0., qtmp(j,i-1))
      d1 = min(-qtmp(j,i), d0)
      qtmp(j, i-1) = qtmp(j, i-1) - d1
      qtmp(j, i) = qtmp(j, i) + d1
      ! east
      d0 = max(0., qtmp(j,1))
      d2 = min(-qtmp(j,i), d0)
      qtmp(j, 1) = qtmp(j, 1) - d2

      qtmp(j, i) = qtmp(j, i) + d2 + tiny
    END IF
  END DO

  IF (ipx/=0) THEN
    DO j = j1, j2
      DO i = 1, imr
        q(i, j) = qtmp(j, i)
      END DO
    END DO
  ELSE

    ! Poles.
    IF (q(1,1)<0 .OR. q(1,jnp)<0.) ipx = 1
  END IF
  RETURN
END SUBROUTINE filew
