module ppm3d_m

  implicit none

contains

  SUBROUTINE ppm3d(igd, q, ps1, ps2, u, v, w, ndt, iord, jord, kord, nc, imr, &
       jnp, j1, nlay, ap, bp, pt, ae, fill, umax)

    ! rajout de d\'eclarations
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

    integer, PARAMETER:: jmax=361, kmax=150

    ! ****6***0*********0*********0*********0*********0*********0**********72

    ! Input-Output arrays

    integer imr
    INTEGER igd, iord, jord, kord, nc, jnp, j1, nlay, ae
    REAL q(imr, jnp, nlay, nc), ps1(imr, jnp), ps2(imr, jnp), &
         u(imr, jnp, nlay), v(imr, jnp, nlay), ap(nlay+1), bp(nlay+1), &
         w(imr, jnp, nlay), ndt, umax
    INTEGER imrd2
    REAL pt
    LOGICAL cross, fill

    ! Local dynamic arrays

    REAL crx(imr, jnp), cry(imr, jnp), xmass(imr, jnp), ymass(imr, jnp), &
         fx1(imr+1), dpi(imr, jnp, nlay), delp1(imr, jnp, nlay), &
         wk1(imr, jnp, nlay), pu(imr, jnp), pv(imr, jnp), dc2(imr, jnp), &
         delp2(imr, jnp, nlay), dq(imr, jnp, nlay, nc), va(imr, jnp), &
         ua(imr, jnp), qtmp(-imr:2*imr)

    ! Local static  arrays

    REAL dtdx(jmax), dtdx5(jmax), acosp(jmax), cosp(jmax), cose(jmax), &
         dap(kmax), dbk(kmax)
    integer ndt0, nstep
    DATA ndt0, nstep/0, 0/
    DATA cross/.TRUE./
    SAVE dtdy, dtdy5, rcap, js0, jn0, iml, dtdx, dtdx5, acosp, cosp, cose, dap, dbk
    real cr1, d5, dl, dp, dt, dtdy, dtdy5, pi, rcap, ru, sum1, sum2, ztc
    integer i, iad, ic, iiu, imh, imjm, iml, iu, j, j2, jad, jmp, jmr, jn, jn0
    integer js, js0, jt, k, krd, l, maxdt

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

2222   CONTINUE
       DO j = jn0, j2 - 1
          DO i = 1, imr
             IF (abs(crx(i,j))>1.) THEN
                jn = j
                GO TO 2233
             END IF
          END DO
       END DO
2233   CONTINUE

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
             IF (j>js .AND. j<jn) cycle

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
          END DO

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
                ! j=1 c'est le p\^ole Sud, j=JNP c'est le p\^ole Nord
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

  END SUBROUTINE ppm3d

end module ppm3d_m
