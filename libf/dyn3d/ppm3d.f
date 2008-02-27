!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/ppm3d.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!

cFrom lin@explorer.gsfc.nasa.gov Wed Apr 15 17:44:44 1998
cDate: Wed, 15 Apr 1998 11:37:03 -0400
cFrom: lin@explorer.gsfc.nasa.gov
cTo: Frederic.Hourdin@lmd.jussieu.fr
cSubject: 3D transport module of the GSFC CTM and GEOS GCM


cThis code is sent to you by S-J Lin, DAO, NASA-GSFC

cNote: this version is intended for machines like CRAY
C-90. No multitasking directives implemented.

      
C ********************************************************************
C
C TransPort Core for Goddard Chemistry Transport Model (G-CTM), Goddard
C Earth Observing System General Circulation Model (GEOS-GCM), and Data
C Assimilation System (GEOS-DAS).
C
C ********************************************************************
C
C Purpose: given horizontal winds on  a hybrid sigma-p surfaces,
C          one call to tpcore updates the 3-D mixing ratio
C          fields one time step (NDT). [vertical mass flux is computed
C          internally consistent with the discretized hydrostatic mass
C          continuity equation of the C-Grid GEOS-GCM (for IGD=1)].
C
C Schemes: Multi-dimensional Flux Form Semi-Lagrangian (FFSL) scheme based
C          on the van Leer or PPM.
C          (see Lin and Rood 1996).
C Version 4.5
C Last modified: Dec. 5, 1996
C Major changes from version 4.0: a more general vertical hybrid sigma-
C pressure coordinate.
C Subroutines modified: xtp, ytp, fzppm, qckxyz
C Subroutines deleted: vanz
C
C Author: Shian-Jiann Lin
C mail address:
C                 Shian-Jiann Lin*
C                 Code 910.3, NASA/GSFC, Greenbelt, MD 20771
C                 Phone: 301-286-9540
C                 E-mail: lin@dao.gsfc.nasa.gov
C
C *affiliation:
C                 Joint Center for Earth Systems Technology
C                 The University of Maryland Baltimore County
C                 NASA - Goddard Space Flight Center
C References:
C
C 1. Lin, S.-J., and R. B. Rood, 1996: Multidimensional flux form semi-
C    Lagrangian transport schemes. Mon. Wea. Rev., 124, 2046-2070.
C
C 2. Lin, S.-J., W. C. Chao, Y. C. Sud, and G. K. Walker, 1994: A class of
C    the van Leer-type transport schemes and its applications to the moist-
C    ure transport in a General Circulation Model. Mon. Wea. Rev., 122,
C    1575-1593.
C
C ****6***0*********0*********0*********0*********0*********0**********72
C
      subroutine ppm3d(IGD,Q,PS1,PS2,U,V,W,NDT,IORD,JORD,KORD,NC,IMR,
     &                  JNP,j1,NLAY,AP,BP,PT,AE,fill,dum,Umax)

c      implicit none

c     rajout de déclarations
c      integer Jmax,kmax,ndt0,nstep,k,j,i,ic,l,js,jn,imh,iad,jad,krd
c      integer iu,iiu,j2,jmr,js0,jt
c      real dtdy,dtdy5,rcap,iml,jn0,imjm,pi,dl,dp
c      real dt,cr1,maxdt,ztc,d5,sum1,sum2,ru
C
C ********************************************************************
C
C =============
C INPUT:
C =============
C
C Q(IMR,JNP,NLAY,NC): mixing ratios at current time (t)
C NC: total # of constituents
C IMR: first dimension (E-W); # of Grid intervals in E-W is IMR
C JNP: 2nd dimension (N-S); # of Grid intervals in N-S is JNP-1
C NLAY: 3rd dimension (# of layers); vertical index increases from 1 at
C       the model top to NLAY near the surface (see fig. below).
C       It is assumed that 6 <= NLAY <= JNP (for dynamic memory allocation)
C
C PS1(IMR,JNP): surface pressure at current time (t)
C PS2(IMR,JNP): surface pressure at mid-time-level (t+NDT/2)
C PS2 is replaced by the predicted PS (at t+NDT) on output.
C Note: surface pressure can have any unit or can be multiplied by any
C       const.
C
C The pressure at layer edges are defined as follows:
C
C        p(i,j,k) = AP(k)*PT  +  BP(k)*PS(i,j)          (1)
C
C Where PT is a constant having the same unit as PS.
C AP and BP are unitless constants given at layer edges
C defining the vertical coordinate. 
C BP(1) = 0., BP(NLAY+1) = 1.
C The pressure at the model top is PTOP = AP(1)*PT
C
C For pure sigma system set AP(k) = 1 for all k, PT = PTOP,
C BP(k) = sige(k) (sigma at edges), PS = Psfc - PTOP.
C
C Note: the sigma-P coordinate is a subset of Eq. 1, which in turn
C is a subset of the following even more general sigma-P-thelta coord.
C currently under development.
C  p(i,j,k) = (AP(k)*PT + BP(k)*PS(i,j))/(D(k)-C(k)*TE**(-1/kapa))
C
C                  /////////////////////////////////
C              / \ ------------- PTOP --------------  AP(1), BP(1)
C               |
C    delp(1)    |  ........... Q(i,j,1) ............  
C               |
C      W(1)    \ / ---------------------------------  AP(2), BP(2)
C
C
C
C     W(k-1)   / \ ---------------------------------  AP(k), BP(k)
C               |
C    delp(K)    |  ........... Q(i,j,k) ............ 
C               |
C      W(k)    \ / ---------------------------------  AP(k+1), BP(k+1)
C
C
C
C              / \ ---------------------------------  AP(NLAY), BP(NLAY)
C               |
C  delp(NLAY)   |  ........... Q(i,j,NLAY) .........  
C               |
C   W(NLAY)=0  \ / ------------- surface ----------- AP(NLAY+1), BP(NLAY+1)
C                 //////////////////////////////////
C
C U(IMR,JNP,NLAY) & V(IMR,JNP,NLAY):winds (m/s) at mid-time-level (t+NDT/2)
C U and V may need to be polar filtered in advance in some cases.
C 
C IGD:      grid type on which winds are defined.
C IGD = 0:  A-Grid  [all variables defined at the same point from south
C                   pole (j=1) to north pole (j=JNP) ]
C
C IGD = 1  GEOS-GCM C-Grid
C                                      [North]
C
C                                       V(i,j)
C                                          |
C                                          |
C                                          |
C                             U(i-1,j)---Q(i,j)---U(i,j) [EAST]
C                                          |
C                                          |
C                                          |
C                                       V(i,j-1)
C
C         U(i,  1) is defined at South Pole.
C         V(i,  1) is half grid north of the South Pole.
C         V(i,JMR) is half grid south of the North Pole.
C
C         V must be defined at j=1 and j=JMR if IGD=1
C         V at JNP need not be given.
C
C NDT: time step in seconds (need not be constant during the course of
C      the integration). Suggested value: 30 min. for 4x5, 15 min. for 2x2.5
C      (Lat-Lon) resolution. Smaller values are recommanded if the model
C      has a well-resolved stratosphere.
C
C J1 defines the size of the polar cap:
C South polar cap edge is located at -90 + (j1-1.5)*180/(JNP-1) deg.
C North polar cap edge is located at  90 - (j1-1.5)*180/(JNP-1) deg.
C There are currently only two choices (j1=2 or 3).
C IMR must be an even integer if j1 = 2. Recommended value: J1=3.
C
C IORD, JORD, and KORD are integers controlling various options in E-W, N-S,
C and vertical transport, respectively. Recommended values for positive
C definite scalars: IORD=JORD=3, KORD=5. Use KORD=3 for non-
C positive definite scalars or when linear correlation between constituents
C is to be maintained.
C
C  _ORD= 
C        1: 1st order upstream scheme (too diffusive, not a useful option; it
C           can be used for debugging purposes; this is THE only known "linear"
C           monotonic advection scheme.).
C        2: 2nd order van Leer (full monotonicity constraint;
C           see Lin et al 1994, MWR)
C        3: monotonic PPM* (slightly improved PPM of Collela & Woodward 1984)
C        4: semi-monotonic PPM (same as 3, but overshoots are allowed)
C        5: positive-definite PPM (constraint on the subgrid distribution is
C           only strong enough to prevent generation of negative values;
C           both overshoots & undershoots are possible).
C        6: un-constrained PPM (nearly diffusion free; slightly faster but
C           positivity not quaranteed. Use this option only when the fields
C           and winds are very smooth).
C
C *PPM: Piece-wise Parabolic Method
C
C Note that KORD <=2 options are no longer supported. DO not use option 4 or 5.
C for non-positive definite scalars (such as Ertel Potential Vorticity).
C
C The implicit numerical diffusion decreases as _ORD increases.
C The last two options (ORDER=5, 6) should only be used when there is
C significant explicit diffusion (such as a turbulence parameterization). You
C might get dispersive results otherwise.
C No filter of any kind is applied to the constituent fields here.
C
C AE: Radius of the sphere (meters).
C     Recommended value for the planet earth: 6.371E6
C
C fill(logical):   flag to do filling for negatives (see note below).
C
C Umax: Estimate (upper limit) of the maximum U-wind speed (m/s).
C (220 m/s is a good value for troposphere model; 280 m/s otherwise)
C
C =============
C Output
C =============
C
C Q: mixing ratios at future time (t+NDT) (original values are over-written)
C W(NLAY): large-scale vertical mass flux as diagnosed from the hydrostatic
C          relationship. W will have the same unit as PS1 and PS2 (eg, mb).
C          W must be divided by NDT to get the correct mass-flux unit.
C          The vertical Courant number C = W/delp_UPWIND, where delp_UPWIND
C          is the pressure thickness in the "upwind" direction. For example,
C          C(k) = W(k)/delp(k)   if W(k) > 0;
C          C(k) = W(k)/delp(k+1) if W(k) < 0.
C              ( W > 0 is downward, ie, toward surface)
C PS2: predicted PS at t+NDT (original values are over-written)
C
C ********************************************************************
C NOTES:
C This forward-in-time upstream-biased transport scheme reduces to
C the 2nd order center-in-time center-in-space mass continuity eqn.
C if Q = 1 (constant fields will remain constant). This also ensures
C that the computed vertical velocity to be identical to GEOS-1 GCM
C for on-line transport.
C
C A larger polar cap is used if j1=3 (recommended for C-Grid winds or when
C winds are noisy near poles).
C
C Flux-Form Semi-Lagrangian transport in the East-West direction is used
C when and where Courant # is greater than one.
C
C The user needs to change the parameter Jmax or Kmax if the resolution
C is greater than 0.5 deg in N-S or 150 layers in the vertical direction.
C (this TransPort Core is otherwise resolution independent and can be used
C as a library routine).
C
C PPM is 4th order accurate when grid spacing is uniform (x & y); 3rd
C order accurate for non-uniform grid (vertical sigma coord.).
C
C Time step is limitted only by transport in the meridional direction.
C (the FFSL scheme is not implemented in the meridional direction).
C
C Since only 1-D limiters are applied, negative values could
C potentially be generated when large time step is used and when the
C initial fields contain discontinuities.
C This does not necessarily imply the integration is unstable.
C These negatives are typically very small. A filling algorithm is
C activated if the user set "fill" to be true.
C
C The van Leer scheme used here is nearly as accurate as the original PPM
C due to the use of a 4th order accurate reference slope. The PPM imple-
C mented here is an improvement over the original and is also based on
C the 4th order reference slope.
C
C ****6***0*********0*********0*********0*********0*********0**********72
C
C     User modifiable parameters
C
      parameter (Jmax = 361, kmax = 150)
C
C ****6***0*********0*********0*********0*********0*********0**********72
C
C Input-Output arrays
C
      
      real Q(IMR,JNP,NLAY,NC),PS1(IMR,JNP),PS2(IMR,JNP),
     &     U(IMR,JNP,NLAY),V(IMR,JNP,NLAY),AP(NLAY+1),
     &     BP(NLAY+1),W(IMR,JNP,NLAY),NDT,val(NLAY),Umax
      integer IGD,IORD,JORD,KORD,NC,IMR,JNP,j1,NLAY,AE
      integer IMRD2
      real    PT       
      logical  cross, fill, dum
C
C Local dynamic arrays
C
      real CRX(IMR,JNP),CRY(IMR,JNP),xmass(IMR,JNP),ymass(IMR,JNP),
     &     fx1(IMR+1),DPI(IMR,JNP,NLAY),delp1(IMR,JNP,NLAY),
     &     WK1(IMR,JNP,NLAY),PU(IMR,JNP),PV(IMR,JNP),DC2(IMR,JNP),
     &     delp2(IMR,JNP,NLAY),DQ(IMR,JNP,NLAY,NC),VA(IMR,JNP),
     &     UA(IMR,JNP),qtmp(-IMR:2*IMR)
C
C Local static  arrays
C
      real DTDX(Jmax), DTDX5(Jmax), acosp(Jmax),
     &     cosp(Jmax), cose(Jmax), DAP(kmax),DBK(Kmax)
      data NDT0, NSTEP /0, 0/
      data cross /.true./
      SAVE DTDY, DTDY5, RCAP, JS0, JN0, IML,
     &     DTDX, DTDX5, ACOSP, COSP, COSE, DAP,DBK
C
            
      JMR = JNP -1
      IMJM  = IMR*JNP
      j2 = JNP - j1 + 1
      NSTEP = NSTEP + 1
C
C *********** Initialization **********************
      if(NSTEP.eq.1) then
c
      write(6,*) '------------------------------------ '
      write(6,*) 'NASA/GSFC Transport Core Version 4.5'
      write(6,*) '------------------------------------ '
c
      WRITE(6,*) 'IMR=',IMR,' JNP=',JNP,' NLAY=',NLAY,' j1=',j1
      WRITE(6,*) 'NC=',NC,IORD,JORD,KORD,NDT
C
C controles sur les parametres
      if(NLAY.LT.6) then
        write(6,*) 'NLAY must be >= 6'
        stop
      endif
      if (JNP.LT.NLAY) then
         write(6,*) 'JNP must be >= NLAY'
        stop
      endif
      IMRD2=mod(IMR,2)
      if (j1.eq.2.and.IMRD2.NE.0) then
         write(6,*) 'if j1=2 IMR must be an even integer'
        stop
      endif

C
      if(Jmax.lt.JNP .or. Kmax.lt.NLAY) then
        write(6,*) 'Jmax or Kmax is too small'
        stop
      endif
C
      DO k=1,NLAY
      DAP(k) = (AP(k+1) - AP(k))*PT
      DBK(k) =  BP(k+1) - BP(k)
      ENDDO     
C
      PI = 4. * ATAN(1.)
      DL = 2.*PI / float(IMR)
      DP =    PI / float(JMR)
C
      if(IGD.eq.0) then
C Compute analytic cosine at cell edges
            call cosa(cosp,cose,JNP,PI,DP)
      else
C Define cosine consistent with GEOS-GCM (using dycore2.0 or later)
            call cosc(cosp,cose,JNP,PI,DP)
      endif
C
      do 15 J=2,JMR
15    acosp(j) = 1. / cosp(j)
C
C Inverse of the Scaled polar cap area.
C
      RCAP  = DP / (IMR*(1.-COS((j1-1.5)*DP)))
      acosp(1)   = RCAP
      acosp(JNP) = RCAP
      endif
C
      if(NDT0 .ne. NDT) then
      DT   = NDT
      NDT0 = NDT

	if(Umax .lt. 180.) then
         write(6,*) 'Umax may be too small!'
	endif
      CR1  = abs(Umax*DT)/(DL*AE)
      MaxDT = DP*AE / abs(Umax) + 0.5
      write(6,*)'Largest time step for max(V)=',Umax,' is ',MaxDT
      if(MaxDT .lt. abs(NDT)) then
            write(6,*) 'Warning!!! NDT maybe too large!'
      endif
C
      if(CR1.ge.0.95) then
      JS0 = 0
      JN0 = 0
      IML = IMR-2
      ZTC = 0.
      else
      ZTC  = acos(CR1) * (180./PI)
C
      JS0 = float(JMR)*(90.-ZTC)/180. + 2
      JS0 = max(JS0, J1+1)
      IML = min(6*JS0/(J1-1)+2, 4*IMR/5)
      JN0 = JNP-JS0+1
      endif
C     
C
      do J=2,JMR
      DTDX(j)  = DT / ( DL*AE*COSP(J) )

c     print*,'dtdx=',dtdx(j)
      DTDX5(j) = 0.5*DTDX(j)
      enddo
C
      
      DTDY  = DT /(AE*DP)
c      print*,'dtdy=',dtdy
      DTDY5 = 0.5*DTDY
C
c      write(6,*) 'J1=',J1,' J2=', J2
      endif
C
C *********** End Initialization **********************
C
C delp = pressure thickness: the psudo-density in a hydrostatic system.
      do  k=1,NLAY
         do  j=1,JNP
            do  i=1,IMR
               delp1(i,j,k)=DAP(k)+DBK(k)*PS1(i,j)
               delp2(i,j,k)=DAP(k)+DBK(k)*PS2(i,j)       
            enddo
         enddo
      enddo
          
C
      if(j1.ne.2) then
      DO 40 IC=1,NC
      DO 40 L=1,NLAY
      DO 40 I=1,IMR
      Q(I,  2,L,IC) = Q(I,  1,L,IC)
40    Q(I,JMR,L,IC) = Q(I,JNP,L,IC)
      endif
C
C Compute "tracer density"
      DO 550 IC=1,NC
      DO 44 k=1,NLAY
      DO 44 j=1,JNP
      DO 44 i=1,IMR
44    DQ(i,j,k,IC) = Q(i,j,k,IC)*delp1(i,j,k)
550	continue
C
      do 1500 k=1,NLAY
C
      if(IGD.eq.0) then
C Convert winds on A-Grid to Courant # on C-Grid.
      call A2C(U(1,1,k),V(1,1,k),IMR,JMR,j1,j2,CRX,CRY,dtdx5,DTDY5)
      else
C Convert winds on C-grid to Courant #
      do 45 j=j1,j2
      do 45 i=2,IMR
45    CRX(i,J) = dtdx(j)*U(i-1,j,k)
   
C
      do 50 j=j1,j2
50    CRX(1,J) = dtdx(j)*U(IMR,j,k)
C
      do 55 i=1,IMR*JMR
55    CRY(i,2) = DTDY*V(i,1,k)
      endif
C     
C Determine JS and JN
      JS = j1
      JN = j2
C
      do j=JS0,j1+1,-1
      do i=1,IMR
      if(abs(CRX(i,j)).GT.1.) then
            JS = j
            go to 2222
      endif
      enddo
      enddo
C
2222  continue
      do j=JN0,j2-1
      do i=1,IMR
      if(abs(CRX(i,j)).GT.1.) then
            JN = j
            go to 2233
      endif
      enddo
      enddo
2233  continue
C
      if(j1.ne.2) then           ! Enlarged polar cap.
      do i=1,IMR
      DPI(i,  2,k) = 0.
      DPI(i,JMR,k) = 0.
      enddo
      endif
C
C ******* Compute horizontal mass fluxes ************
C
C N-S component
      do j=j1,j2+1
      D5 = 0.5 * COSE(j)
      do i=1,IMR
      ymass(i,j) = CRY(i,j)*D5*(delp2(i,j,k) + delp2(i,j-1,k))
      enddo
      enddo
C
      do 95 j=j1,j2
      DO 95 i=1,IMR
95    DPI(i,j,k) = (ymass(i,j) - ymass(i,j+1)) * acosp(j)
C
C Poles
      sum1 = ymass(IMR,j1  )
      sum2 = ymass(IMR,J2+1)
      do i=1,IMR-1
      sum1 = sum1 + ymass(i,j1  )
      sum2 = sum2 + ymass(i,J2+1)
      enddo
C
      sum1 = - sum1 * RCAP
      sum2 =   sum2 * RCAP
      do i=1,IMR
      DPI(i,  1,k) = sum1
      DPI(i,JNP,k) = sum2
      enddo
C
C E-W component
C
      do j=j1,j2
      do i=2,IMR
      PU(i,j) = 0.5 * (delp2(i,j,k) + delp2(i-1,j,k))
      enddo
      enddo
C
      do j=j1,j2
      PU(1,j) = 0.5 * (delp2(1,j,k) + delp2(IMR,j,k))
      enddo
C
      do 110 j=j1,j2
      DO 110 i=1,IMR
110   xmass(i,j) = PU(i,j)*CRX(i,j)
C
      DO 120 j=j1,j2
      DO 120 i=1,IMR-1
120   DPI(i,j,k) = DPI(i,j,k) + xmass(i,j) - xmass(i+1,j)
C
      DO 130 j=j1,j2
130   DPI(IMR,j,k) = DPI(IMR,j,k) + xmass(IMR,j) - xmass(1,j)
C
      DO j=j1,j2
      do i=1,IMR-1
      UA(i,j) = 0.5 * (CRX(i,j)+CRX(i+1,j))
      enddo
      enddo
C
      DO j=j1,j2
      UA(imr,j) = 0.5 * (CRX(imr,j)+CRX(1,j))
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Rajouts pour LMDZ.3.3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,IMR
         do j=1,JNP
             VA(i,j)=0.
         enddo
      enddo

      do i=1,imr*(JMR-1)
      VA(i,2) = 0.5*(CRY(i,2)+CRY(i,3))
      enddo
C
      if(j1.eq.2) then
	IMH = IMR/2
      do i=1,IMH
      VA(i,      1) = 0.5*(CRY(i,2)-CRY(i+IMH,2))
      VA(i+IMH,  1) = -VA(i,1)
      VA(i,    JNP) = 0.5*(CRY(i,JNP)-CRY(i+IMH,JMR))
      VA(i+IMH,JNP) = -VA(i,JNP)
      enddo
      VA(IMR,1)=VA(1,1)
      VA(IMR,JNP)=VA(1,JNP)
      endif
C
C ****6***0*********0*********0*********0*********0*********0**********72
      do 1000 IC=1,NC
C
      do i=1,IMJM
      wk1(i,1,1) = 0.
      wk1(i,1,2) = 0.
      enddo
C
C E-W advective cross term
      do 250 j=J1,J2
      if(J.GT.JS  .and. J.LT.JN) GO TO 250
C
      do i=1,IMR
      qtmp(i) = q(i,j,k,IC)
      enddo
C
      do i=-IML,0
      qtmp(i)       = q(IMR+i,j,k,IC)
      qtmp(IMR+1-i) = q(1-i,j,k,IC)
      enddo
C
      DO 230 i=1,IMR
      iu = UA(i,j)
      ru = UA(i,j) - iu
      iiu = i-iu
      if(UA(i,j).GE.0.) then
      wk1(i,j,1) = qtmp(iiu)+ru*(qtmp(iiu-1)-qtmp(iiu))
      else
      wk1(i,j,1) = qtmp(iiu)+ru*(qtmp(iiu)-qtmp(iiu+1))
      endif
      wk1(i,j,1) = wk1(i,j,1) - qtmp(i)
230   continue
250   continue
C
      if(JN.ne.0) then
      do j=JS+1,JN-1
C
      do i=1,IMR
      qtmp(i) = q(i,j,k,IC)
      enddo
C
      qtmp(0)     = q(IMR,J,k,IC)
      qtmp(IMR+1) = q(  1,J,k,IC)
C
      do i=1,imr
      iu = i - UA(i,j)
      wk1(i,j,1) = UA(i,j)*(qtmp(iu) - qtmp(iu+1))
      enddo
      enddo
      endif
C ****6***0*********0*********0*********0*********0*********0**********72
C Contribution from the N-S advection
      do i=1,imr*(j2-j1+1)
      JT = float(J1) - VA(i,j1)
      wk1(i,j1,2) = VA(i,j1) * (q(i,jt,k,IC) - q(i,jt+1,k,IC))
      enddo
C
      do i=1,IMJM
      wk1(i,1,1) = q(i,1,k,IC) + 0.5*wk1(i,1,1)
      wk1(i,1,2) = q(i,1,k,IC) + 0.5*wk1(i,1,2)
      enddo
C
	if(cross) then
C Add cross terms in the vertical direction.
	if(IORD .GE. 2) then
		iad = 2
	else
		iad = 1
	endif
C
	if(JORD .GE. 2) then
		jad = 2
	else
		jad = 1
	endif
      call xadv(IMR,JNP,j1,j2,wk1(1,1,2),UA,JS,JN,IML,DC2,iad)
      call yadv(IMR,JNP,j1,j2,wk1(1,1,1),VA,PV,W,jad)
      do j=1,JNP
      do i=1,IMR
      q(i,j,k,IC) = q(i,j,k,IC) + DC2(i,j) + PV(i,j)
      enddo
      enddo
      endif
C
      call xtp(IMR,JNP,IML,j1,j2,JN,JS,PU,DQ(1,1,k,IC),wk1(1,1,2)
     &        ,CRX,fx1,xmass,IORD)

      call ytp(IMR,JNP,j1,j2,acosp,RCAP,DQ(1,1,k,IC),wk1(1,1,1),CRY,
     &  DC2,ymass,WK1(1,1,3),wk1(1,1,4),WK1(1,1,5),WK1(1,1,6),JORD)
C
1000  continue
1500  continue
C
C ******* Compute vertical mass flux (same unit as PS) ***********
C
C 1st step: compute total column mass CONVERGENCE.
C
      do 320 j=1,JNP
      do 320 i=1,IMR
320   CRY(i,j) = DPI(i,j,1)
C
      do 330 k=2,NLAY
      do 330 j=1,JNP
      do 330 i=1,IMR
      CRY(i,j)  = CRY(i,j) + DPI(i,j,k)
330   continue
C
      do 360 j=1,JNP
      do 360 i=1,IMR
C
C 2nd step: compute PS2 (PS at n+1) using the hydrostatic assumption.
C Changes (increases) to surface pressure = total column mass convergence
C
      PS2(i,j)  = PS1(i,j) + CRY(i,j)
C
C 3rd step: compute vertical mass flux from mass conservation principle.
C
      W(i,j,1) = DPI(i,j,1) - DBK(1)*CRY(i,j)
      W(i,j,NLAY) = 0.
360   continue
C
      do 370 k=2,NLAY-1
      do 370 j=1,JNP
      do 370 i=1,IMR
      W(i,j,k) = W(i,j,k-1) + DPI(i,j,k) - DBK(k)*CRY(i,j)
370   continue
C
      DO 380 k=1,NLAY
      DO 380 j=1,JNP
      DO 380 i=1,IMR
      delp2(i,j,k) = DAP(k) + DBK(k)*PS2(i,j)
380   continue
C
	KRD = max(3, KORD)
      do 4000 IC=1,NC
C
C****6***0*********0*********0*********0*********0*********0**********72
   
      call FZPPM(IMR,JNP,NLAY,j1,DQ(1,1,1,IC),W,Q(1,1,1,IC),WK1,DPI,
     &           DC2,CRX,CRY,PU,PV,xmass,ymass,delp1,KRD)
C
    
      if(fill) call qckxyz(DQ(1,1,1,IC),DC2,IMR,JNP,NLAY,j1,j2,
     &                     cosp,acosp,.false.,IC,NSTEP)
C
C Recover tracer mixing ratio from "density" using predicted
C "air density" (pressure thickness) at time-level n+1
C
      DO k=1,NLAY
      DO j=1,JNP
      DO i=1,IMR
            Q(i,j,k,IC) = DQ(i,j,k,IC) / delp2(i,j,k)
c            print*,'i=',i,'j=',j,'k=',k,'Q(i,j,k,IC)=',Q(i,j,k,IC)
      enddo
      enddo
      enddo
C     
      if(j1.ne.2) then
      DO 400 k=1,NLAY
      DO 400 I=1,IMR
c     j=1 c'est le pôle Sud, j=JNP c'est le pôle Nord
      Q(I,  2,k,IC) = Q(I,  1,k,IC)
      Q(I,JMR,k,IC) = Q(I,JMP,k,IC)
400   CONTINUE
      endif
4000  continue
C
      if(j1.ne.2) then
      DO 5000 k=1,NLAY
      DO 5000 i=1,IMR
      W(i,  2,k) = W(i,  1,k)
      W(i,JMR,k) = W(i,JNP,k)
5000  continue
      endif
C
      RETURN
      END
C
C****6***0*********0*********0*********0*********0*********0**********72
      subroutine FZPPM(IMR,JNP,NLAY,j1,DQ,WZ,P,DC,DQDT,AR,AL,A6,
     &                 flux,wk1,wk2,wz2,delp,KORD)
      parameter ( kmax = 150 )
      parameter ( R23 = 2./3., R3 = 1./3.)
      real WZ(IMR,JNP,NLAY),P(IMR,JNP,NLAY),DC(IMR,JNP,NLAY),
     &     wk1(IMR,*),delp(IMR,JNP,NLAY),DQ(IMR,JNP,NLAY),
     &     DQDT(IMR,JNP,NLAY)
C Assuming JNP >= NLAY
      real AR(IMR,*),AL(IMR,*),A6(IMR,*),flux(IMR,*),wk2(IMR,*),
     &     wz2(IMR,*)
C
      JMR = JNP - 1
      IMJM = IMR*JNP
      NLAYM1 = NLAY - 1
C
      LMT = KORD - 3
C
C ****6***0*********0*********0*********0*********0*********0**********72
C Compute DC for PPM
C ****6***0*********0*********0*********0*********0*********0**********72
C
      do 1000 k=1,NLAYM1
      do 1000 i=1,IMJM
      DQDT(i,1,k) = P(i,1,k+1) - P(i,1,k)
1000  continue
C
      DO 1220 k=2,NLAYM1
      DO 1220 I=1,IMJM    
       c0 =  delp(i,1,k) / (delp(i,1,k-1)+delp(i,1,k)+delp(i,1,k+1))
       c1 = (delp(i,1,k-1)+0.5*delp(i,1,k))/(delp(i,1,k+1)+delp(i,1,k))    
       c2 = (delp(i,1,k+1)+0.5*delp(i,1,k))/(delp(i,1,k-1)+delp(i,1,k))
      tmp = c0*(c1*DQDT(i,1,k) + c2*DQDT(i,1,k-1))
      Qmax = max(P(i,1,k-1),P(i,1,k),P(i,1,k+1)) - P(i,1,k)
      Qmin = P(i,1,k) - min(P(i,1,k-1),P(i,1,k),P(i,1,k+1))
      DC(i,1,k) = sign(min(abs(tmp),Qmax,Qmin), tmp)   
1220  CONTINUE
     
C     
C ****6***0*********0*********0*********0*********0*********0**********72
C Loop over latitudes  (to save memory)
C ****6***0*********0*********0*********0*********0*********0**********72
C
      DO 2000 j=1,JNP
      if((j.eq.2 .or. j.eq.JMR) .and. j1.ne.2) goto 2000
C
      DO k=1,NLAY
      DO i=1,IMR
      wz2(i,k) =   WZ(i,j,k)
      wk1(i,k) =    P(i,j,k)
      wk2(i,k) = delp(i,j,k)
      flux(i,k) = DC(i,j,k)  !this flux is actually the monotone slope
      enddo
      enddo
C
C****6***0*********0*********0*********0*********0*********0**********72
C Compute first guesses at cell interfaces
C First guesses are required to be continuous.
C ****6***0*********0*********0*********0*********0*********0**********72
C
C three-cell parabolic subgrid distribution at model top
C two-cell parabolic with zero gradient subgrid distribution 
C at the surface.
C
C First guess top edge value
      DO 10 i=1,IMR
C three-cell PPM
C Compute a,b, and c of q = aP**2 + bP + c using cell averages and delp
      a = 3.*( DQDT(i,j,2) - DQDT(i,j,1)*(wk2(i,2)+wk2(i,3))/
     &         (wk2(i,1)+wk2(i,2)) ) /
     &       ( (wk2(i,2)+wk2(i,3))*(wk2(i,1)+wk2(i,2)+wk2(i,3)) )
      b = 2.*DQDT(i,j,1)/(wk2(i,1)+wk2(i,2)) - 
     &    R23*a*(2.*wk2(i,1)+wk2(i,2))
      AL(i,1) =  wk1(i,1) - wk2(i,1)*(R3*a*wk2(i,1) + 0.5*b)
      AL(i,2) =  wk2(i,1)*(a*wk2(i,1) + b) + AL(i,1)
C
C Check if change sign
      if(wk1(i,1)*AL(i,1).le.0.) then
		 AL(i,1) = 0.
             flux(i,1) = 0.
	else
             flux(i,1) =  wk1(i,1) - AL(i,1)
	endif
10    continue
C
C Bottom
      DO 15 i=1,IMR
C 2-cell PPM with zero gradient right at the surface
C
      fct = DQDT(i,j,NLAYM1)*wk2(i,NLAY)**2 /
     & ( (wk2(i,NLAY)+wk2(i,NLAYM1))*(2.*wk2(i,NLAY)+wk2(i,NLAYM1)))
      AR(i,NLAY) = wk1(i,NLAY) + fct
      AL(i,NLAY) = wk1(i,NLAY) - (fct+fct)
      if(wk1(i,NLAY)*AR(i,NLAY).le.0.) AR(i,NLAY) = 0.
      flux(i,NLAY) = AR(i,NLAY) -  wk1(i,NLAY)
15    continue
     
C
C****6***0*********0*********0*********0*********0*********0**********72
C 4th order interpolation in the interior.
C****6***0*********0*********0*********0*********0*********0**********72
C
      DO 14 k=3,NLAYM1
      DO 12 i=1,IMR
      c1 =  DQDT(i,j,k-1)*wk2(i,k-1) / (wk2(i,k-1)+wk2(i,k))
      c2 =  2. / (wk2(i,k-2)+wk2(i,k-1)+wk2(i,k)+wk2(i,k+1))
      A1   =  (wk2(i,k-2)+wk2(i,k-1)) / (2.*wk2(i,k-1)+wk2(i,k))
      A2   =  (wk2(i,k  )+wk2(i,k+1)) / (2.*wk2(i,k)+wk2(i,k-1))
      AL(i,k) = wk1(i,k-1) + c1 + c2 *
     &        ( wk2(i,k  )*(c1*(A1 - A2)+A2*flux(i,k-1)) -
     &          wk2(i,k-1)*A1*flux(i,k)  )
C      print *,'AL1',i,k, AL(i,k)
12    CONTINUE
14    continue
C
      do 20 i=1,IMR*NLAYM1
      AR(i,1) = AL(i,2)
C      print *,'AR1',i,AR(i,1)
20    continue
C
      do 30 i=1,IMR*NLAY
      A6(i,1) = 3.*(wk1(i,1)+wk1(i,1) - (AL(i,1)+AR(i,1)))
C      print *,'A61',i,A6(i,1)
30    continue
C
C****6***0*********0*********0*********0*********0*********0**********72
C Top & Bot always monotonic
      call lmtppm(flux(1,1),A6(1,1),AR(1,1),AL(1,1),wk1(1,1),IMR,0)
      call lmtppm(flux(1,NLAY),A6(1,NLAY),AR(1,NLAY),AL(1,NLAY),
     &            wk1(1,NLAY),IMR,0)
C
C Interior depending on KORD
      if(LMT.LE.2)
     &  call lmtppm(flux(1,2),A6(1,2),AR(1,2),AL(1,2),wk1(1,2),
     &              IMR*(NLAY-2),LMT)
C
C****6***0*********0*********0*********0*********0*********0**********72
C
      DO 140 i=1,IMR*NLAYM1
      IF(wz2(i,1).GT.0.) then
        CM = wz2(i,1) / wk2(i,1)
        flux(i,2) = AR(i,1)+0.5*CM*(AL(i,1)-AR(i,1)+A6(i,1)*(1.-R23*CM))
      else
C        print *,'test2-0',i,j,wz2(i,1),wk2(i,2)
        CP= wz2(i,1) / wk2(i,2)        
C        print *,'testCP',CP
        flux(i,2) = AL(i,2)+0.5*CP*(AL(i,2)-AR(i,2)-A6(i,2)*(1.+R23*CP))
C        print *,'test2',i, AL(i,2),AR(i,2),A6(i,2),R23
      endif
140   continue
C
      DO 250 i=1,IMR*NLAYM1
      flux(i,2) = wz2(i,1) * flux(i,2)
250   continue
C
      do 350 i=1,IMR
      DQ(i,j,   1) = DQ(i,j,   1) - flux(i,   2)
      DQ(i,j,NLAY) = DQ(i,j,NLAY) + flux(i,NLAY)
350   continue
C
      do 360 k=2,NLAYM1
      do 360 i=1,IMR
360   DQ(i,j,k) = DQ(i,j,k) + flux(i,k) - flux(i,k+1)
2000  continue
      return
      end
C
      subroutine xtp(IMR,JNP,IML,j1,j2,JN,JS,PU,DQ,Q,UC,
     &               fx1,xmass,IORD)
      dimension UC(IMR,*),DC(-IML:IMR+IML+1),xmass(IMR,JNP)
     &    ,fx1(IMR+1),DQ(IMR,JNP),qtmp(-IML:IMR+1+IML)
      dimension PU(IMR,JNP),Q(IMR,JNP),ISAVE(IMR)
C
      IMP = IMR + 1
C
C van Leer at high latitudes
      jvan = max(1,JNP/18)
      j1vl = j1+jvan
      j2vl = j2-jvan
C
      do 1310 j=j1,j2
C
      do i=1,IMR
      qtmp(i) = q(i,j)
      enddo
C
      if(j.ge.JN .or. j.le.JS) goto 2222
C ************* Eulerian **********
C
      qtmp(0)     = q(IMR,J)
      qtmp(-1)    = q(IMR-1,J)
      qtmp(IMP)   = q(1,J)
      qtmp(IMP+1) = q(2,J)
C
      IF(IORD.eq.1 .or. j.eq.j1. or. j.eq.j2) THEN
      DO 1406 i=1,IMR
      iu = float(i) - uc(i,j)
1406  fx1(i) = qtmp(iu)
      ELSE
      call xmist(IMR,IML,Qtmp,DC)
      DC(0) = DC(IMR)
C
      if(IORD.eq.2 .or. j.le.j1vl .or. j.ge.j2vl) then
      DO 1408 i=1,IMR
      iu = float(i) - uc(i,j)
1408  fx1(i) = qtmp(iu) + DC(iu)*(sign(1.,uc(i,j))-uc(i,j))
      else
      call fxppm(IMR,IML,UC(1,j),Qtmp,DC,fx1,IORD)
      endif
C
      ENDIF
C
      DO 1506 i=1,IMR
1506  fx1(i) = fx1(i)*xmass(i,j)
C
      goto 1309
C
C ***** Conservative (flux-form) Semi-Lagrangian transport *****
C
2222  continue
C
      do i=-IML,0
      qtmp(i)     = q(IMR+i,j)
      qtmp(IMP-i) = q(1-i,j)
      enddo
C
      IF(IORD.eq.1 .or. j.eq.j1. or. j.eq.j2) THEN
      DO 1306 i=1,IMR
      itmp = INT(uc(i,j))
      ISAVE(i) = i - itmp
      iu = i - uc(i,j)
1306  fx1(i) = (uc(i,j) - itmp)*qtmp(iu)
      ELSE
      call xmist(IMR,IML,Qtmp,DC)
C
      do i=-IML,0
      DC(i)     = DC(IMR+i)
      DC(IMP-i) = DC(1-i)
      enddo
C
      DO 1307 i=1,IMR
      itmp = INT(uc(i,j))
      rut  = uc(i,j) - itmp
      ISAVE(i) = i - itmp
      iu = i - uc(i,j)
1307  fx1(i) = rut*(qtmp(iu) + DC(iu)*(sign(1.,rut) - rut))
      ENDIF
C
      do 1308 i=1,IMR
      IF(uc(i,j).GT.1.) then
CDIR$ NOVECTOR
        do ist = ISAVE(i),i-1
        fx1(i) = fx1(i) + qtmp(ist)
        enddo
      elseIF(uc(i,j).LT.-1.) then
        do ist = i,ISAVE(i)-1
        fx1(i) = fx1(i) - qtmp(ist)
        enddo
CDIR$ VECTOR
      endif
1308  continue
      do i=1,IMR
      fx1(i) = PU(i,j)*fx1(i)
      enddo
C
C ***************************************
C
1309  fx1(IMP) = fx1(1)
      DO 1215 i=1,IMR
1215  DQ(i,j) =  DQ(i,j) + fx1(i)-fx1(i+1)
C
C ***************************************
C
1310  continue
      return
      end
C
      subroutine fxppm(IMR,IML,UT,P,DC,flux,IORD)
      parameter ( R3 = 1./3., R23 = 2./3. )
      DIMENSION UT(*),flux(*),P(-IML:IMR+IML+1),DC(-IML:IMR+IML+1)
      DIMENSION AR(0:IMR),AL(0:IMR),A6(0:IMR)
      integer LMT 
c      logical first
c      data first /.true./
c      SAVE LMT
c      if(first) then
C
C correction calcul de LMT a chaque passage pour pouvoir choisir
c plusieurs schemas PPM pour differents traceurs
c      IF (IORD.LE.0) then
c            if(IMR.GE.144) then
c                  LMT = 0
c            elseif(IMR.GE.72) then
c                  LMT = 1
c            else
c                  LMT = 2
c            endif
c      else
c            LMT = IORD - 3
c      endif
C
      LMT = IORD - 3
c      write(6,*) 'PPM option in E-W direction = ', LMT
c      first = .false.
C      endif
C
      DO 10 i=1,IMR
10    AL(i) = 0.5*(p(i-1)+p(i)) + (DC(i-1) - DC(i))*R3
C
      do 20 i=1,IMR-1
20    AR(i) = AL(i+1)
      AR(IMR) = AL(1)
C
      do 30 i=1,IMR
30    A6(i) = 3.*(p(i)+p(i)  - (AL(i)+AR(i)))
C
      if(LMT.LE.2) call lmtppm(DC(1),A6(1),AR(1),AL(1),P(1),IMR,LMT)
C
      AL(0) = AL(IMR)
      AR(0) = AR(IMR)
      A6(0) = A6(IMR)
C
      DO i=1,IMR
      IF(UT(i).GT.0.) then
      flux(i) = AR(i-1) + 0.5*UT(i)*(AL(i-1) - AR(i-1) +
     &                 A6(i-1)*(1.-R23*UT(i)) )
      else
      flux(i) = AL(i) - 0.5*UT(i)*(AR(i) - AL(i) +
     &                        A6(i)*(1.+R23*UT(i)))
      endif
      enddo
      return
      end
C
      subroutine xmist(IMR,IML,P,DC)
      parameter( R24 = 1./24.)
      dimension P(-IML:IMR+1+IML),DC(-IML:IMR+1+IML)
C
      do 10  i=1,IMR
      tmp = R24*(8.*(p(i+1) - p(i-1)) + p(i-2) - p(i+2))
      Pmax = max(P(i-1), p(i), p(i+1)) - p(i)
      Pmin = p(i) - min(P(i-1), p(i), p(i+1))
10    DC(i) = sign(min(abs(tmp),Pmax,Pmin), tmp)
      return
      end
C
      subroutine ytp(IMR,JNP,j1,j2,acosp,RCAP,DQ,P,VC,DC2
     &              ,ymass,fx,A6,AR,AL,JORD)
      dimension P(IMR,JNP),VC(IMR,JNP),ymass(IMR,JNP)
     &       ,DC2(IMR,JNP),DQ(IMR,JNP),acosp(JNP)
C Work array
      DIMENSION fx(IMR,JNP),AR(IMR,JNP),AL(IMR,JNP),A6(IMR,JNP)
C
      JMR = JNP - 1
      len = IMR*(J2-J1+2)
C
      if(JORD.eq.1) then
      DO 1000 i=1,len
      JT = float(J1) - VC(i,J1)
1000  fx(i,j1) = p(i,JT)
      else
   
      call ymist(IMR,JNP,j1,P,DC2,4)
C
      if(JORD.LE.0 .or. JORD.GE.3) then
   
      call fyppm(VC,P,DC2,fx,IMR,JNP,j1,j2,A6,AR,AL,JORD)
    
      else
      DO 1200 i=1,len
      JT = float(J1) - VC(i,J1)
1200  fx(i,j1) = p(i,JT) + (sign(1.,VC(i,j1))-VC(i,j1))*DC2(i,JT)
      endif
      endif
C
      DO 1300 i=1,len
1300  fx(i,j1) = fx(i,j1)*ymass(i,j1)
C
      DO 1400 j=j1,j2
      DO 1400 i=1,IMR
1400  DQ(i,j) = DQ(i,j) + (fx(i,j) - fx(i,j+1)) * acosp(j)
C
C Poles
      sum1 = fx(IMR,j1  )
      sum2 = fx(IMR,J2+1)
      do i=1,IMR-1
      sum1 = sum1 + fx(i,j1  )
      sum2 = sum2 + fx(i,J2+1)
      enddo
C
      sum1 = DQ(1,  1) - sum1 * RCAP
      sum2 = DQ(1,JNP) + sum2 * RCAP
      do i=1,IMR
      DQ(i,  1) = sum1
      DQ(i,JNP) = sum2
      enddo
C
      if(j1.ne.2) then
      do i=1,IMR
      DQ(i,  2) = sum1
      DQ(i,JMR) = sum2
      enddo
      endif
C
      return
      end
C
      subroutine  ymist(IMR,JNP,j1,P,DC,ID)
      parameter ( R24 = 1./24. )
      dimension P(IMR,JNP),DC(IMR,JNP)
C
      IMH = IMR / 2
      JMR = JNP - 1
      IJM3 = IMR*(JMR-3)
C
      IF(ID.EQ.2) THEN
      do 10 i=1,IMR*(JMR-1)
      tmp = 0.25*(p(i,3) - p(i,1))
      Pmax = max(p(i,1),p(i,2),p(i,3)) - p(i,2)
      Pmin = p(i,2) - min(p(i,1),p(i,2),p(i,3))
      DC(i,2) = sign(min(abs(tmp),Pmin,Pmax),tmp)
10    CONTINUE
      ELSE
      do 12 i=1,IMH
C J=2
      tmp = (8.*(p(i,3) - p(i,1)) + p(i+IMH,2) - p(i,4))*R24
      Pmax = max(p(i,1),p(i,2),p(i,3)) - p(i,2)
      Pmin = p(i,2) - min(p(i,1),p(i,2),p(i,3))
      DC(i,2) = sign(min(abs(tmp),Pmin,Pmax),tmp)
C J=JMR
      tmp=(8.*(p(i,JNP)-p(i,JMR-1))+p(i,JMR-2)-p(i+IMH,JMR))*R24
      Pmax = max(p(i,JMR-1),p(i,JMR),p(i,JNP)) - p(i,JMR)
      Pmin = p(i,JMR) - min(p(i,JMR-1),p(i,JMR),p(i,JNP))
      DC(i,JMR) = sign(min(abs(tmp),Pmin,Pmax),tmp)
12    CONTINUE
      do 14 i=IMH+1,IMR
C J=2
      tmp = (8.*(p(i,3) - p(i,1)) + p(i-IMH,2) - p(i,4))*R24
      Pmax = max(p(i,1),p(i,2),p(i,3)) - p(i,2)
      Pmin = p(i,2) - min(p(i,1),p(i,2),p(i,3))
      DC(i,2) = sign(min(abs(tmp),Pmin,Pmax),tmp)
C J=JMR
      tmp=(8.*(p(i,JNP)-p(i,JMR-1))+p(i,JMR-2)-p(i-IMH,JMR))*R24
      Pmax = max(p(i,JMR-1),p(i,JMR),p(i,JNP)) - p(i,JMR)
      Pmin = p(i,JMR) - min(p(i,JMR-1),p(i,JMR),p(i,JNP))
      DC(i,JMR) = sign(min(abs(tmp),Pmin,Pmax),tmp)
14    CONTINUE
C
      do 15 i=1,IJM3
      tmp = (8.*(p(i,4) - p(i,2)) + p(i,1) - p(i,5))*R24
      Pmax = max(p(i,2),p(i,3),p(i,4)) - p(i,3)
      Pmin = p(i,3) - min(p(i,2),p(i,3),p(i,4))
      DC(i,3) = sign(min(abs(tmp),Pmin,Pmax),tmp)
15    CONTINUE
      ENDIF
C
      if(j1.ne.2) then
      do i=1,IMR
      DC(i,1) = 0.
      DC(i,JNP) = 0.
      enddo
      else
C Determine slopes in polar caps for scalars!
C
      do 13 i=1,IMH
C South
      tmp = 0.25*(p(i,2) - p(i+imh,2))
      Pmax = max(p(i,2),p(i,1), p(i+imh,2)) - p(i,1)
      Pmin = p(i,1) - min(p(i,2),p(i,1), p(i+imh,2))
      DC(i,1)=sign(min(abs(tmp),Pmax,Pmin),tmp)
C North.
      tmp = 0.25*(p(i+imh,JMR) - p(i,JMR))
      Pmax = max(p(i+imh,JMR),p(i,jnp), p(i,JMR)) - p(i,JNP)
      Pmin = p(i,JNP) - min(p(i+imh,JMR),p(i,jnp), p(i,JMR))
      DC(i,JNP) = sign(min(abs(tmp),Pmax,pmin),tmp)
13    continue
C
      do 25 i=imh+1,IMR
      DC(i,  1) =  - DC(i-imh,  1)
      DC(i,JNP) =  - DC(i-imh,JNP)
25    continue
      endif
      return
      end
C
      subroutine fyppm(VC,P,DC,flux,IMR,JNP,j1,j2,A6,AR,AL,JORD)
      parameter ( R3 = 1./3., R23 = 2./3. )
      real VC(IMR,*),flux(IMR,*),P(IMR,*),DC(IMR,*)
C Local work arrays.
      real AR(IMR,JNP),AL(IMR,JNP),A6(IMR,JNP)
      integer LMT
c      logical first
C      data first /.true./
C      SAVE LMT
C
      IMH = IMR / 2
      JMR = JNP - 1
      j11 = j1-1
      IMJM1 = IMR*(J2-J1+2)
      len   = IMR*(J2-J1+3)
C      if(first) then
C      IF(JORD.LE.0) then
C            if(JMR.GE.90) then
C                  LMT = 0
C            elseif(JMR.GE.45) then
C                  LMT = 1
C            else
C                  LMT = 2
C            endif
C      else
C            LMT = JORD - 3
C      endif
C
C      first = .false.
C      endif
C     
c modifs pour pouvoir choisir plusieurs schemas PPM
      LMT = JORD - 3      
C
      DO 10 i=1,IMR*JMR        
      AL(i,2) = 0.5*(p(i,1)+p(i,2)) + (DC(i,1) - DC(i,2))*R3
      AR(i,1) = AL(i,2)
10    CONTINUE
C
CPoles:
C
      DO i=1,IMH
      AL(i,1) = AL(i+IMH,2)
      AL(i+IMH,1) = AL(i,2)
C
      AR(i,JNP) = AR(i+IMH,JMR)
      AR(i+IMH,JNP) = AR(i,JMR)
      ENDDO

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Rajout pour LMDZ.3.3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      AR(IMR,1)=AL(1,1)
      AR(IMR,JNP)=AL(1,JNP)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
           
      do 30 i=1,len
30    A6(i,j11) = 3.*(p(i,j11)+p(i,j11)  - (AL(i,j11)+AR(i,j11)))
C
      if(LMT.le.2) call lmtppm(DC(1,j11),A6(1,j11),AR(1,j11)
     &                       ,AL(1,j11),P(1,j11),len,LMT)
C
     
      DO 140 i=1,IMJM1
      IF(VC(i,j1).GT.0.) then
      flux(i,j1) = AR(i,j11) + 0.5*VC(i,j1)*(AL(i,j11) - AR(i,j11) +
     &                         A6(i,j11)*(1.-R23*VC(i,j1)) )
      else
      flux(i,j1) = AL(i,j1) - 0.5*VC(i,j1)*(AR(i,j1) - AL(i,j1) +
     &                        A6(i,j1)*(1.+R23*VC(i,j1)))
      endif
140   continue
      return
      end
C
	subroutine yadv(IMR,JNP,j1,j2,p,VA,ady,wk,IAD)
	REAL p(IMR,JNP),ady(IMR,JNP),VA(IMR,JNP)
        REAL WK(IMR,-1:JNP+2)
C
	JMR = JNP-1
	IMH = IMR/2
	do j=1,JNP
	do i=1,IMR
	wk(i,j) = p(i,j)
	enddo
	enddo
C Poles:
	do i=1,IMH
	wk(i,   -1) = p(i+IMH,3)
	wk(i+IMH,-1) = p(i,3)
	wk(i,    0) = p(i+IMH,2)
	wk(i+IMH,0) = p(i,2)
	wk(i,JNP+1) = p(i+IMH,JMR)
	wk(i+IMH,JNP+1) = p(i,JMR)
	wk(i,JNP+2) = p(i+IMH,JNP-2)
	wk(i+IMH,JNP+2) = p(i,JNP-2)
	enddo
c        write(*,*) 'toto 1' 
C --------------------------------
      IF(IAD.eq.2) then
      do j=j1-1,j2+1
      do i=1,IMR
c      write(*,*) 'avt NINT','i=',i,'j=',j
      JP = NINT(VA(i,j))      
      rv = JP - VA(i,j)
c      write(*,*) 'VA=',VA(i,j), 'JP1=',JP,'rv=',rv
      JP = j - JP
c      write(*,*) 'JP2=',JP
      a1 = 0.5*(wk(i,jp+1)+wk(i,jp-1)) - wk(i,jp)
      b1 = 0.5*(wk(i,jp+1)-wk(i,jp-1))
c      write(*,*) 'a1=',a1,'b1=',b1
      ady(i,j) = wk(i,jp) + rv*(a1*rv + b1) - wk(i,j)
      enddo
      enddo
c      write(*,*) 'toto 2'
C
      ELSEIF(IAD.eq.1) then
	do j=j1-1,j2+1
      do i=1,imr
      JP = float(j)-VA(i,j)
      ady(i,j) = VA(i,j)*(wk(i,jp)-wk(i,jp+1))
      enddo
      enddo
      ENDIF
C
	if(j1.ne.2) then
	sum1 = 0.
	sum2 = 0.
      do i=1,imr
      sum1 = sum1 + ady(i,2)
      sum2 = sum2 + ady(i,JMR)
      enddo
	sum1 = sum1 / IMR
	sum2 = sum2 / IMR
C
      do i=1,imr
      ady(i,  2) =  sum1
      ady(i,JMR) =  sum2
      ady(i,  1) =  sum1
      ady(i,JNP) =  sum2
      enddo
	else
C Poles:
	sum1 = 0.
	sum2 = 0.
      do i=1,imr
      sum1 = sum1 + ady(i,1)
      sum2 = sum2 + ady(i,JNP)
      enddo
	sum1 = sum1 / IMR
	sum2 = sum2 / IMR
C
      do i=1,imr
      ady(i,  1) =  sum1
      ady(i,JNP) =  sum2
      enddo
	endif
C
	return
	end
C
	subroutine xadv(IMR,JNP,j1,j2,p,UA,JS,JN,IML,adx,IAD)
	REAL p(IMR,JNP),adx(IMR,JNP),qtmp(-IMR:IMR+IMR),UA(IMR,JNP)
C
	JMR = JNP-1
      do 1309 j=j1,j2
      if(J.GT.JS  .and. J.LT.JN) GO TO 1309
C
      do i=1,IMR
      qtmp(i) = p(i,j)
      enddo
C
      do i=-IML,0
      qtmp(i)       = p(IMR+i,j)
      qtmp(IMR+1-i) = p(1-i,j)
      enddo
C
      IF(IAD.eq.2) THEN
      DO i=1,IMR
      IP = NINT(UA(i,j))
      ru = IP - UA(i,j)
      IP = i - IP
      a1 = 0.5*(qtmp(ip+1)+qtmp(ip-1)) - qtmp(ip)
      b1 = 0.5*(qtmp(ip+1)-qtmp(ip-1))
      adx(i,j) = qtmp(ip) + ru*(a1*ru + b1)
      enddo
      ELSEIF(IAD.eq.1) then
      DO i=1,IMR
      iu = UA(i,j)
      ru = UA(i,j) - iu
      iiu = i-iu
      if(UA(i,j).GE.0.) then
      adx(i,j) = qtmp(iiu)+ru*(qtmp(iiu-1)-qtmp(iiu))
      else
      adx(i,j) = qtmp(iiu)+ru*(qtmp(iiu)-qtmp(iiu+1))
      endif
      enddo
      ENDIF
C
      do i=1,IMR
      adx(i,j) = adx(i,j) - p(i,j)
      enddo
1309  continue
C
C Eulerian upwind
C
      do j=JS+1,JN-1
C
      do i=1,IMR
      qtmp(i) = p(i,j)
      enddo
C
      qtmp(0)     = p(IMR,J)
      qtmp(IMR+1) = p(1,J)
C
      IF(IAD.eq.2) THEN
      qtmp(-1)     = p(IMR-1,J)
      qtmp(IMR+2) = p(2,J)
      do i=1,imr
      IP = NINT(UA(i,j))
      ru = IP - UA(i,j)
      IP = i - IP
      a1 = 0.5*(qtmp(ip+1)+qtmp(ip-1)) - qtmp(ip)
      b1 = 0.5*(qtmp(ip+1)-qtmp(ip-1))
      adx(i,j) = qtmp(ip)- p(i,j) + ru*(a1*ru + b1)
      enddo
      ELSEIF(IAD.eq.1) then
C 1st order
      DO i=1,IMR
      IP = i - UA(i,j)
      adx(i,j) = UA(i,j)*(qtmp(ip)-qtmp(ip+1))
      enddo
      ENDIF
      enddo
C
	if(j1.ne.2) then
      do i=1,IMR
      adx(i,  2) = 0.
      adx(i,JMR) = 0.
      enddo
	endif
C set cross term due to x-adv at the poles to zero.
      do i=1,IMR
      adx(i,  1) = 0.
      adx(i,JNP) = 0.
      enddo
	return
	end
C
      subroutine lmtppm(DC,A6,AR,AL,P,IM,LMT)
C
C A6 =  CURVATURE OF THE TEST PARABOLA
C AR =  RIGHT EDGE VALUE OF THE TEST PARABOLA
C AL =  LEFT  EDGE VALUE OF THE TEST PARABOLA
C DC =  0.5 * MISMATCH
C P  =  CELL-AVERAGED VALUE
C IM =  VECTOR LENGTH
C
C OPTIONS:
C
C LMT = 0: FULL MONOTONICITY
C LMT = 1: SEMI-MONOTONIC CONSTRAINT (NO UNDERSHOOTS)
C LMT = 2: POSITIVE-DEFINITE CONSTRAINT
C
      parameter ( R12 = 1./12. )
      dimension A6(IM),AR(IM),AL(IM),P(IM),DC(IM)
C
      if(LMT.eq.0) then
C Full constraint
      do 100 i=1,IM
      if(DC(i).eq.0.) then
            AR(i) = p(i)
            AL(i) = p(i)
            A6(i) = 0.
      else
      da1  = AR(i) - AL(i)
      da2  = da1**2
      A6DA = A6(i)*da1
      if(A6DA .lt. -da2) then
            A6(i) = 3.*(AL(i)-p(i))
            AR(i) = AL(i) - A6(i)
      elseif(A6DA .gt. da2) then
            A6(i) = 3.*(AR(i)-p(i))
            AL(i) = AR(i) - A6(i)
      endif
      endif
100   continue
      elseif(LMT.eq.1) then
C Semi-monotonic constraint
      do 150 i=1,IM
      if(abs(AR(i)-AL(i)) .GE. -A6(i)) go to 150
      if(p(i).lt.AR(i) .and. p(i).lt.AL(i)) then
            AR(i) = p(i)
            AL(i) = p(i)
            A6(i) = 0.
      elseif(AR(i) .gt. AL(i)) then
            A6(i) = 3.*(AL(i)-p(i))
            AR(i) = AL(i) - A6(i)
      else
            A6(i) = 3.*(AR(i)-p(i))
            AL(i) = AR(i) - A6(i)
      endif
150   continue
      elseif(LMT.eq.2) then
      do 250 i=1,IM
      if(abs(AR(i)-AL(i)) .GE. -A6(i)) go to 250
      fmin = p(i) + 0.25*(AR(i)-AL(i))**2/A6(i) + A6(i)*R12
      if(fmin.ge.0.) go to 250
      if(p(i).lt.AR(i) .and. p(i).lt.AL(i)) then
            AR(i) = p(i)
            AL(i) = p(i)
            A6(i) = 0.
      elseif(AR(i) .gt. AL(i)) then
            A6(i) = 3.*(AL(i)-p(i))
            AR(i) = AL(i) - A6(i)
      else
            A6(i) = 3.*(AR(i)-p(i))
            AL(i) = AR(i) - A6(i)
      endif
250   continue
      endif
      return
      end
C
      subroutine A2C(U,V,IMR,JMR,j1,j2,CRX,CRY,dtdx5,DTDY5)
      dimension U(IMR,*),V(IMR,*),CRX(IMR,*),CRY(IMR,*),DTDX5(*)
C
      do 35 j=j1,j2
      do 35 i=2,IMR
35    CRX(i,J) = dtdx5(j)*(U(i,j)+U(i-1,j))
C
      do 45 j=j1,j2
45    CRX(1,J) = dtdx5(j)*(U(1,j)+U(IMR,j))
C
      do 55 i=1,IMR*JMR
55    CRY(i,2) = DTDY5*(V(i,2)+V(i,1))
      return
      end
C
      subroutine cosa(cosp,cose,JNP,PI,DP)
      dimension cosp(*),cose(*)
      JMR = JNP-1
      do 55 j=2,JNP
        ph5  =  -0.5*PI + (FLOAT(J-1)-0.5)*DP
55      cose(j) = cos(ph5)
C
      JEQ = (JNP+1) / 2
      if(JMR .eq. 2*(JMR/2) ) then
      do j=JNP, JEQ+1, -1
       cose(j) =  cose(JNP+2-j)
      enddo
      else
C cell edge at equator.
       cose(JEQ+1) =  1.
      do j=JNP, JEQ+2, -1
       cose(j) =  cose(JNP+2-j)
       enddo
      endif
C
      do 66 j=2,JMR
66    cosp(j) = 0.5*(cose(j)+cose(j+1))
      cosp(1) = 0.
      cosp(JNP) = 0.
      return
      end
C
      subroutine cosc(cosp,cose,JNP,PI,DP)
      dimension cosp(*),cose(*)
C
      phi = -0.5*PI
      do 55 j=2,JNP-1
      phi  =  phi + DP
55    cosp(j) = cos(phi)
        cosp(  1) = 0.
        cosp(JNP) = 0.
C
      do 66 j=2,JNP
        cose(j) = 0.5*(cosp(j)+cosp(j-1))
66    CONTINUE
C
      do 77 j=2,JNP-1
       cosp(j) = 0.5*(cose(j)+cose(j+1))
77    CONTINUE
      return
      end
C
      SUBROUTINE qckxyz (Q,qtmp,IMR,JNP,NLAY,j1,j2,cosp,acosp,
     &                   cross,IC,NSTEP)
C
      parameter( tiny = 1.E-60 )
      DIMENSION Q(IMR,JNP,NLAY),qtmp(IMR,JNP),cosp(*),acosp(*)
      logical cross
C
      NLAYM1 = NLAY-1
      len = IMR*(j2-j1+1)
      ip = 0
C
C Top layer
      L = 1
	icr = 1
      call filns(q(1,1,L),IMR,JNP,j1,j2,cosp,acosp,ipy,tiny)
      if(ipy.eq.0) goto 50
      call filew(q(1,1,L),qtmp,IMR,JNP,j1,j2,ipx,tiny)
      if(ipx.eq.0) goto 50
C
      if(cross) then
      call filcr(q(1,1,L),IMR,JNP,j1,j2,cosp,acosp,icr,tiny)
      endif
      if(icr.eq.0) goto 50
C
C Vertical filling...
      do i=1,len
      IF( Q(i,j1,1).LT.0.) THEN
      ip = ip + 1
          Q(i,j1,2) = Q(i,j1,2) + Q(i,j1,1)
          Q(i,j1,1) = 0.
      endif
      enddo
C
50    continue
      DO 225 L = 2,NLAYM1
      icr = 1
C
      call filns(q(1,1,L),IMR,JNP,j1,j2,cosp,acosp,ipy,tiny)
      if(ipy.eq.0) goto 225
      call filew(q(1,1,L),qtmp,IMR,JNP,j1,j2,ipx,tiny)
      if(ipx.eq.0) go to 225
      if(cross) then
      call filcr(q(1,1,L),IMR,JNP,j1,j2,cosp,acosp,icr,tiny)
      endif
      if(icr.eq.0) goto 225
C
      do i=1,len
      IF( Q(I,j1,L).LT.0.) THEN
C
      ip = ip + 1
C From above
          qup =  Q(I,j1,L-1)
          qly = -Q(I,j1,L)
          dup  = min(qly,qup)
          Q(I,j1,L-1) = qup - dup
          Q(I,j1,L  ) = dup-qly
C Below
          Q(I,j1,L+1) = Q(I,j1,L+1) + Q(I,j1,L)
          Q(I,j1,L)   = 0.
      ENDIF
      ENDDO
225   CONTINUE
C
C BOTTOM LAYER
      sum = 0.
      L = NLAY
C
      call filns(q(1,1,L),IMR,JNP,j1,j2,cosp,acosp,ipy,tiny)
      if(ipy.eq.0) goto 911
      call filew(q(1,1,L),qtmp,IMR,JNP,j1,j2,ipx,tiny)
      if(ipx.eq.0) goto 911
C
      call filcr(q(1,1,L),IMR,JNP,j1,j2,cosp,acosp,icr,tiny)
      if(icr.eq.0) goto 911
C
      DO  I=1,len
      IF( Q(I,j1,L).LT.0.) THEN
      ip = ip + 1
c
C From above
C
          qup = Q(I,j1,NLAYM1)
          qly = -Q(I,j1,L)
          dup = min(qly,qup)
          Q(I,j1,NLAYM1) = qup - dup
C From "below" the surface.
          sum = sum + qly-dup
          Q(I,j1,L) = 0.
       ENDIF
      ENDDO
C
911   continue
C
      if(ip.gt.IMR) then
      write(6,*) 'IC=',IC,' STEP=',NSTEP,
     &           ' Vertical filling pts=',ip
      endif
C
      if(sum.gt.1.e-25) then
      write(6,*) IC,NSTEP,' Mass source from the ground=',sum
      endif
      RETURN
      END
C
      subroutine filcr(q,IMR,JNP,j1,j2,cosp,acosp,icr,tiny)
      dimension q(IMR,*),cosp(*),acosp(*)
      icr = 0
      do 65 j=j1+1,j2-1
      DO 50 i=1,IMR-1
      IF(q(i,j).LT.0.) THEN
      icr =  1
      dq  = - q(i,j)*cosp(j)
C N-E
      dn = q(i+1,j+1)*cosp(j+1)
      d0 = max(0.,dn)
      d1 = min(dq,d0)
      q(i+1,j+1) = (dn - d1)*acosp(j+1)
      dq = dq - d1
C S-E
      ds = q(i+1,j-1)*cosp(j-1)
      d0 = max(0.,ds)
      d2 = min(dq,d0)
      q(i+1,j-1) = (ds - d2)*acosp(j-1)
      q(i,j) = (d2 - dq)*acosp(j) + tiny
      endif
50    continue
      if(icr.eq.0 .and. q(IMR,j).ge.0.) goto 65
      DO 55 i=2,IMR
      IF(q(i,j).LT.0.) THEN
      icr =  1
      dq  = - q(i,j)*cosp(j)
C N-W
      dn = q(i-1,j+1)*cosp(j+1)
      d0 = max(0.,dn)
      d1 = min(dq,d0)
      q(i-1,j+1) = (dn - d1)*acosp(j+1)
      dq = dq - d1
C S-W
      ds = q(i-1,j-1)*cosp(j-1)
      d0 = max(0.,ds)
      d2 = min(dq,d0)
      q(i-1,j-1) = (ds - d2)*acosp(j-1)
      q(i,j) = (d2 - dq)*acosp(j) + tiny
      endif
55    continue
C *****************************************
C i=1
      i=1
      IF(q(i,j).LT.0.) THEN
      icr =  1
      dq  = - q(i,j)*cosp(j)
C N-W
      dn = q(IMR,j+1)*cosp(j+1)
      d0 = max(0.,dn)
      d1 = min(dq,d0)
      q(IMR,j+1) = (dn - d1)*acosp(j+1)
      dq = dq - d1
C S-W
      ds = q(IMR,j-1)*cosp(j-1)
      d0 = max(0.,ds)
      d2 = min(dq,d0)
      q(IMR,j-1) = (ds - d2)*acosp(j-1)
      q(i,j) = (d2 - dq)*acosp(j) + tiny
      endif
C *****************************************
C i=IMR
      i=IMR
      IF(q(i,j).LT.0.) THEN
      icr =  1
      dq  = - q(i,j)*cosp(j)
C N-E
      dn = q(1,j+1)*cosp(j+1)
      d0 = max(0.,dn)
      d1 = min(dq,d0)
      q(1,j+1) = (dn - d1)*acosp(j+1)
      dq = dq - d1
C S-E
      ds = q(1,j-1)*cosp(j-1)
      d0 = max(0.,ds)
      d2 = min(dq,d0)
      q(1,j-1) = (ds - d2)*acosp(j-1)
      q(i,j) = (d2 - dq)*acosp(j) + tiny
      endif
C *****************************************
65    continue
C
      do i=1,IMR
      if(q(i,j1).lt.0. .or. q(i,j2).lt.0.) then
      icr = 1
      goto 80
      endif
      enddo
C
80    continue
C
      if(q(1,1).lt.0. .or. q(1,jnp).lt.0.) then
      icr = 1
      endif
C
      return
      end
C
      subroutine filns(q,IMR,JNP,j1,j2,cosp,acosp,ipy,tiny)
      dimension q(IMR,*),cosp(*),acosp(*)
c      logical first
c      data first /.true./
c      save cap1
C
c      if(first) then
      DP = 4.*ATAN(1.)/float(JNP-1)
      CAP1 = IMR*(1.-COS((j1-1.5)*DP))/DP
c      first = .false.
c      endif
C
      ipy = 0
      do 55 j=j1+1,j2-1
      DO 55 i=1,IMR
      IF(q(i,j).LT.0.) THEN
      ipy =  1
      dq  = - q(i,j)*cosp(j)
C North
      dn = q(i,j+1)*cosp(j+1)
      d0 = max(0.,dn)
      d1 = min(dq,d0)
      q(i,j+1) = (dn - d1)*acosp(j+1)
      dq = dq - d1
C South
      ds = q(i,j-1)*cosp(j-1)
      d0 = max(0.,ds)
      d2 = min(dq,d0)
      q(i,j-1) = (ds - d2)*acosp(j-1)
      q(i,j) = (d2 - dq)*acosp(j) + tiny
      endif
55    continue
C
      do i=1,imr
      IF(q(i,j1).LT.0.) THEN
      ipy =  1
      dq  = - q(i,j1)*cosp(j1)
C North
      dn = q(i,j1+1)*cosp(j1+1)
      d0 = max(0.,dn)
      d1 = min(dq,d0)
      q(i,j1+1) = (dn - d1)*acosp(j1+1)
      q(i,j1) = (d1 - dq)*acosp(j1) + tiny
      endif
      enddo
C
      j = j2
      do i=1,imr
      IF(q(i,j).LT.0.) THEN
      ipy =  1
      dq  = - q(i,j)*cosp(j)
C South
      ds = q(i,j-1)*cosp(j-1)
      d0 = max(0.,ds)
      d2 = min(dq,d0)
      q(i,j-1) = (ds - d2)*acosp(j-1)
      q(i,j) = (d2 - dq)*acosp(j) + tiny
      endif
      enddo
C
C Check Poles.
      if(q(1,1).lt.0.) then
      dq = q(1,1)*cap1/float(IMR)*acosp(j1)
      do i=1,imr
      q(i,1) = 0.
      q(i,j1) = q(i,j1) + dq
      if(q(i,j1).lt.0.) ipy = 1
      enddo
      endif
C
      if(q(1,JNP).lt.0.) then
      dq = q(1,JNP)*cap1/float(IMR)*acosp(j2)
      do i=1,imr
      q(i,JNP) = 0.
      q(i,j2) = q(i,j2) + dq
      if(q(i,j2).lt.0.) ipy = 1
      enddo
      endif
C
      return
      end
C
      subroutine filew(q,qtmp,IMR,JNP,j1,j2,ipx,tiny)
      dimension q(IMR,*),qtmp(JNP,IMR)
C
      ipx = 0
C Copy & swap direction for vectorization.
      do 25 i=1,imr
      do 25 j=j1,j2
25    qtmp(j,i) = q(i,j)
C
      do 55 i=2,imr-1
      do 55 j=j1,j2
      if(qtmp(j,i).lt.0.) then
      ipx =  1
c west
      d0 = max(0.,qtmp(j,i-1))
      d1 = min(-qtmp(j,i),d0)
      qtmp(j,i-1) = qtmp(j,i-1) - d1
      qtmp(j,i) = qtmp(j,i) + d1
c east
      d0 = max(0.,qtmp(j,i+1))
      d2 = min(-qtmp(j,i),d0)
      qtmp(j,i+1) = qtmp(j,i+1) - d2
      qtmp(j,i) = qtmp(j,i) + d2 + tiny
      endif
55    continue
c
      i=1
      do 65 j=j1,j2
      if(qtmp(j,i).lt.0.) then
      ipx =  1
c west
      d0 = max(0.,qtmp(j,imr))
      d1 = min(-qtmp(j,i),d0)
      qtmp(j,imr) = qtmp(j,imr) - d1
      qtmp(j,i) = qtmp(j,i) + d1
c east
      d0 = max(0.,qtmp(j,i+1))
      d2 = min(-qtmp(j,i),d0)
      qtmp(j,i+1) = qtmp(j,i+1) - d2
c
      qtmp(j,i) = qtmp(j,i) + d2 + tiny
      endif
65    continue
      i=IMR
      do 75 j=j1,j2
      if(qtmp(j,i).lt.0.) then
      ipx =  1
c west
      d0 = max(0.,qtmp(j,i-1))
      d1 = min(-qtmp(j,i),d0)
      qtmp(j,i-1) = qtmp(j,i-1) - d1
      qtmp(j,i) = qtmp(j,i) + d1
c east
      d0 = max(0.,qtmp(j,1))
      d2 = min(-qtmp(j,i),d0)
      qtmp(j,1) = qtmp(j,1) - d2
c
      qtmp(j,i) = qtmp(j,i) + d2 + tiny
      endif
75    continue
C
      if(ipx.ne.0) then
      do 85 j=j1,j2
      do 85 i=1,imr
85    q(i,j) = qtmp(j,i)
      else
C
C Poles.
      if(q(1,1).lt.0. or. q(1,JNP).lt.0.) ipx = 1
      endif
      return
      end
C
      subroutine zflip(q,im,km,nc)
C This routine flip the array q (in the vertical).
      real q(im,km,nc)
C local dynamic array
      real qtmp(im,km)
C
      do 4000 IC = 1, nc
C
      do 1000 k=1,km
      do 1000 i=1,im
      qtmp(i,k) = q(i,km+1-k,IC)
1000  continue
C
      do 2000 i=1,im*km
2000  q(i,1,IC) = qtmp(i,1)
4000  continue
      return
      end
