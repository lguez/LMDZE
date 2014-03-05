
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/convect3.F,v 1.1.1.1 2004/05/19
! 12:53:09 lmdzadmin Exp $

SUBROUTINE convect3(dtime, epmax, ok_adj, t1, r1, rs, u, v, tra, p, ph, nd, &
    ndp1, nl, ntra, delt, iflag, ft, fr, fu, fv, ftra, precip, icb, inb, &
    upwd, dnwd, dnwd0, sig, w0, mike, mke, ma, ments, qents, tps, tls, sigij, &
    cape, tvp, pbase, buoybase, &  ! ccc     *
                                   ! DTVPDT1,DTVPDQ1,DPLCLDT,DPLCLDR)
    dtvpdt1, dtvpdq1, dplcldt, dplcldr, & ! sbl
    ft2, fr2, fu2, fv2, wd, qcond, qcondc) ! sbl

  ! ***  THE PARAMETER NA SHOULD IN GENERAL EQUAL ND   ***


  ! Fleur       Introduction des traceurs dans convect3 le 6 juin 200

  USE dimens_m
  USE dimphy
  USE suphec_m

  REAL, INTENT (IN) :: dtime, delt
  PARAMETER (na=60)

  INTEGER ntrac
  PARAMETER (ntrac=nqmx-2)
  REAL deltac ! cld
  PARAMETER (deltac=0.01) ! cld

  INTEGER nent(na)
  REAL t1(nd), r1(nd), rs(nd), u(nd), v(nd), tra(nd, ntra)
  REAL p(nd), ph(ndp1)
  REAL ft(nd), fr(nd), fu(nd), fv(nd), ftra(nd, ntra)
  REAL sig(nd), w0(nd)
  REAL uent(na, na), vent(na, na), traent(na, na, ntrac), tratm(na)
  REAL up(na), vp(na), trap(na, ntrac)
  REAL m(na), mp(na), ment(na, na), qent(na, na), elij(na, na)
  REAL sij(na, na), tvp(na), tv(na), water(na)
  REAL rp(na), ep(na), th(na), wt(na), evap(na), clw(na)
  REAL sigp(na), b(na), tp(na), cpn(na)
  REAL lv(na), lvcp(na), h(na), hp(na), gz(na)
  REAL t(na), rr(na)

  REAL ft2(nd), fr2(nd), fu2(nd), fv2(nd) ! added sbl
  REAL u1(nd), v1(nd) ! added sbl

  REAL buoy(na) !  Lifted parcel buoyancy
  REAL dtvpdt1(nd), dtvpdq1(nd) ! Derivatives of parcel virtual
  ! temperature wrt T1 and Q1
  REAL clw_new(na), qi(na)

  REAL wd, betad ! for gust factor (sb)
  REAL qcondc(nd) ! interface cld param (sb)
  REAL qcond(nd), nqcond(na), wa(na), maa(na), siga(na), axc(na) ! cld

  LOGICAL ice_conv, ok_adj
  PARAMETER (ice_conv=.TRUE.)

  ! ccccccccccccccccccccccccccccccccccccccccccccc
  ! declaration des variables a sortir
  ! cccccccccccccccccccccccccccccccccccccccccccc
  REAL mke(nd)
  REAL mike(nd)
  REAL ma(nd)
  REAL tps(nd) !temperature dans les ascendances non diluees
  REAL tls(nd) !temperature potentielle
  REAL ments(nd, nd)
  REAL qents(nd, nd)
  REAL sigij(klev, klev)
  REAL pbase ! pressure at the cloud base level
  REAL buoybase ! buoyancy at the cloud base level
  ! ccccccccccccccccccccccccccccccccccccccccccccc




  REAL dnwd0(nd) !  precipitation driven unsaturated downdraft flux
  REAL dnwd(nd), dn1 ! in-cloud saturated downdraft mass flux
  REAL upwd(nd), up1 ! in-cloud saturated updraft mass flux

  ! ***         ASSIGN VALUES OF THERMODYNAMIC CONSTANTS        ***
  ! ***             THESE SHOULD BE CONSISTENT WITH             ***
  ! ***              THOSE USED IN CALLING PROGRAM              ***
  ! ***     NOTE: THESE ARE ALSO SPECIFIED IN SUBROUTINE TLIFT  ***

  ! sb      CPD=1005.7
  ! sb      CPV=1870.0
  ! sb      CL=4190.0
  ! sb      CPVMCL=CL-CPV
  ! sb      RV=461.5
  ! sb      RD=287.04
  ! sb      EPS=RD/RV
  ! sb      ALV0=2.501E6
  ! cccccccccccccccccccccc
  ! constantes coherentes avec le modele du Centre Europeen
  ! sb      RD = 1000.0 * 1.380658E-23 * 6.0221367E+23 / 28.9644
  ! sb      RV = 1000.0 * 1.380658E-23 * 6.0221367E+23 / 18.0153
  ! sb      CPD = 3.5 * RD
  ! sb      CPV = 4.0 * RV
  ! sb      CL = 4218.0
  ! sb      CPVMCL=CL-CPV
  ! sb      EPS=RD/RV
  ! sb      ALV0=2.5008E+06
  ! ccccccccccccccccccccc
  ! on utilise les constantes thermo du Centre Europeen: (SB)


  cpd = rcpd
  cpv = rcpv
  cl = rcw
  cpvmcl = cl - cpv
  eps = rd/rv
  alv0 = rlvtt

  nk = 1 ! origin level of the lifted parcel

  ! ccccccccccccccccccccc

  ! ***  INITIALIZE OUTPUT ARRAYS AND PARAMETERS  ***

  DO i = 1, nd
    ft(i) = 0.0
    fr(i) = 0.0
    fu(i) = 0.0
    fv(i) = 0.0

    ft2(i) = 0.0
    fr2(i) = 0.0
    fu2(i) = 0.0
    fv2(i) = 0.0

    DO j = 1, ntra
      ftra(i, j) = 0.0
    END DO

    qcondc(i) = 0.0 ! cld
    qcond(i) = 0.0 ! cld
    nqcond(i) = 0.0 ! cld

    t(i) = t1(i)
    rr(i) = r1(i)
    u1(i) = u(i) ! added sbl
    v1(i) = v(i) ! added sbl
  END DO
  DO i = 1, nl
    rdcp = (rd*(1.-rr(i))+rr(i)*rv)/(cpd*(1.-rr(i))+rr(i)*cpv)
    th(i) = t(i)*(1000.0/p(i))**rdcp
  END DO

  ! ************************************************************
  ! *    CALCUL DES TEMPERATURES POTENTIELLES A SORTIR
  ! ************************************************************
  DO i = 1, nd
    rdcp = (rd*(1.-rr(i))+rr(i)*rv)/(cpd*(1.-rr(i))+rr(i)*cpv)

    tls(i) = t(i)*(1000.0/p(i))**rdcp
  END DO




  ! ***********************************************************


  precip = 0.0
  wd = 0.0 ! sb
  iflag = 1

  ! ***                    SPECIFY PARAMETERS                        ***
  ! ***  PBCRIT IS THE CRITICAL CLOUD DEPTH (MB) BENEATH WHICH THE   ***
  ! ***       PRECIPITATION EFFICIENCY IS ASSUMED TO BE ZERO         ***
  ! ***  PTCRIT IS THE CLOUD DEPTH (MB) ABOVE WHICH THE PRECIP.      ***
  ! ***            EFFICIENCY IS ASSUMED TO BE UNITY                 ***
  ! ***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
  ! ***  SPFAC IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE      ***
  ! ***                        OF CLOUD                              ***
  ! ***    ALPHA AND BETA ARE PARAMETERS THAT CONTROL THE RATE OF    ***
  ! ***                 APPROACH TO QUASI-EQUILIBRIUM                ***
  ! ***    (THEIR STANDARD VALUES ARE 1.0 AND 0.96, RESPECTIVELY)    ***
  ! ***           (BETA MUST BE LESS THAN OR EQUAL TO 1)             ***
  ! ***    DTCRIT IS THE CRITICAL BUOYANCY (K) USED TO ADJUST THE    ***
  ! ***                 APPROACH TO QUASI-EQUILIBRIUM                ***
  ! ***                     IT MUST BE LESS THAN 0                   ***

  pbcrit = 150.0
  ptcrit = 500.0
  sigd = 0.01
  spfac = 0.15
  ! sb:
  ! EPMAX=0.993 ! precip efficiency less than unity
  ! EPMAX=1. ! precip efficiency less than unity

  ! jyg
  ! CC      BETA=0.96
  ! Beta is now expressed as a function of the characteristic time
  ! of the convective process.
  ! CC        Old value : TAU = 15000.   !(for dtime = 600.s)
  ! CC        Other value (inducing little change) :TAU = 8000.
  tau = 8000.
  beta = 1. - dtime/tau
  ! jyg
  ! CC      ALPHA=1.0
  alpha = 1.5E-3*dtime/tau
  ! Increase alpha in order to compensate W decrease
  alpha = alpha*1.5

  ! jyg (voir CONVECT 3)
  ! CC      DTCRIT=-0.2
  dtcrit = -2.
  ! gf&jyg
  ! CC     DT pour l'overshoot.
  dtovsh = -0.2


  ! ***        INCREMENT THE COUNTER       ***

  sig(nd) = sig(nd) + 1.0
  sig(nd) = amin1(sig(nd), 12.1)

  ! ***    IF NOPT IS AN INTEGER OTHER THAN 0, CONVECT     ***
  ! ***     RETURNS ARRAYS T AND R THAT MAY HAVE BEEN      ***
  ! ***  ALTERED BY DRY ADIABATIC ADJUSTMENT; OTHERWISE    ***
  ! ***        THE RETURNED ARRAYS ARE UNALTERED.          ***

  nopt = 0
  ! !      NOPT=1 ! sbl

  ! ***            PERFORM DRY ADIABATIC ADJUSTMENT            ***

  ! ***  DO NOT BYPASS THIS EVEN IF THE CALLING PROGRAM HAS A  ***
  ! ***                BOUNDARY LAYER SCHEME !!!               ***

  IF (ok_adj) THEN ! added sbl

    DO i = nl - 1, 1, -1
      jn = 0
      DO j = i + 1, nl
        IF (th(j)<th(i)) jn = j
      END DO
      IF (jn==0) GO TO 30
      ahm = 0.0
      rm = 0.0
      um = 0.0
      vm = 0.0
      DO k = 1, ntra
        tratm(k) = 0.0
      END DO
      DO j = i, jn
        ahm = ahm + (cpd*(1.-rr(j))+rr(j)*cpv)*t(j)*(ph(j)-ph(j+1))
        rm = rm + rr(j)*(ph(j)-ph(j+1))
        um = um + u(j)*(ph(j)-ph(j+1))
        vm = vm + v(j)*(ph(j)-ph(j+1))
        DO k = 1, ntra
          tratm(k) = tratm(k) + tra(j, k)*(ph(j)-ph(j+1))
        END DO
      END DO
      dphinv = 1./(ph(i)-ph(jn+1))
      rm = rm*dphinv
      um = um*dphinv
      vm = vm*dphinv
      DO k = 1, ntra
        tratm(k) = tratm(k)*dphinv
      END DO
      a2 = 0.0
      DO j = i, jn
        rr(j) = rm
        u(j) = um
        v(j) = vm
        DO k = 1, ntra
          tra(j, k) = tratm(k)
        END DO
        rdcp = (rd*(1.-rr(j))+rr(j)*rv)/(cpd*(1.-rr(j))+rr(j)*cpv)
        x = (0.001*p(j))**rdcp
        t(j) = x
        a2 = a2 + (cpd*(1.-rr(j))+rr(j)*cpv)*x*(ph(j)-ph(j+1))
      END DO
      DO j = i, jn
        th(j) = ahm/a2
        t(j) = t(j)*th(j)
      END DO
30  END DO

  END IF ! added sbl

  ! ***   RESET INPUT ARRAYS IF ok_adj 0   ***

  IF (ok_adj) THEN
    DO i = 1, nd

      ft2(i) = (t(i)-t1(i))/delt ! sbl
      fr2(i) = (rr(i)-r1(i))/delt ! sbl
      fu2(i) = (u(i)-u1(i))/delt ! sbl
      fv2(i) = (v(i)-v1(i))/delt ! sbl

      ! !            T1(I)=T(I)      ! commente sbl
      ! !            R1(I)=RR(I)     ! commente sbl
    END DO
  END IF

  ! *** CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY AND STATIC ENERGY

  gz(1) = 0.0
  cpn(1) = cpd*(1.-rr(1)) + rr(1)*cpv
  h(1) = t(1)*cpn(1)
  DO i = 2, nl
    tvx = t(i)*(1.+rr(i)/eps-rr(i))
    tvy = t(i-1)*(1.+rr(i-1)/eps-rr(i-1))
    gz(i) = gz(i-1) + 0.5*rd*(tvx+tvy)*(p(i-1)-p(i))/ph(i)
    cpn(i) = cpd*(1.-rr(i)) + cpv*rr(i)
    h(i) = t(i)*cpn(i) + gz(i)
  END DO

  ! ***  CALCULATE LIFTED CONDENSATION LEVEL OF AIR AT LOWEST MODEL LEVEL ***
  ! ***       (WITHIN 0.2% OF FORMULA OF BOLTON, MON. WEA. REV.,1980)     ***

  IF (t(1)<250.0 .OR. rr(1)<=0.0) THEN
    iflag = 0
    ! sb3d         print*,'je suis passe par 366'
    RETURN
  END IF

  ! jyg1 Utilisation de la subroutine CLIFT
  ! C      RH=RR(1)/RS(1)
  ! C      CHI=T(1)/(1669.0-122.0*RH-T(1))
  ! C      PLCL=P(1)*(RH**CHI)
  CALL clift(p(1), t(1), rr(1), rs(1), plcl, dplcldt, dplcldr)
  ! jyg2
  ! sb3d      PRINT *,' em_plcl,p1,t1,r1,rs1,rh '
  ! sb3d     $        ,PLCL,P(1),T(1),RR(1),RS(1),RH

  IF (plcl<200.0 .OR. plcl>=2000.0) THEN
    iflag = 2
    RETURN
  END IF
  ! jyg1
  ! Essais de modification de ICB

  ! ***  CALCULATE FIRST LEVEL ABOVE LCL (=ICB)  ***

  ! C      ICB=NL-1
  ! C      DO 50 I=2,NL-1
  ! C         IF(P(I).LT.PLCL)THEN
  ! C            ICB=MIN(ICB,I)   ! ICB sup ou egal a 2
  ! C         END IF
  ! C   50 CONTINUE
  ! C      IF(ICB.EQ.(NL-1))THEN
  ! C         IFLAG=3
  ! C         RETURN
  ! C      END IF

  ! *** CALCULATE LAYER CONTAINING LCL (=ICB)   ***

  icb = nl - 1
  ! sb      DO 50 I=2,NL-1
  DO i = 3, nl - 1 ! modif sb pour que ICB soit sup/egal a 2
    ! la modification consiste a comparer PLCL a PH et non a P:
    ! ICB est defini par :  PH(ICB)<PLCL<PH(ICB-!)
    IF (ph(i)<plcl) THEN
      icb = min(icb, i)
    END IF
  END DO
  IF (icb==(nl-1)) THEN
    iflag = 3
    RETURN
  END IF
  icb = icb - 1 ! ICB sup ou egal a 2
  ! jyg2



  ! *** SUBROUTINE TLIFT CALCULATES PART OF THE LIFTED PARCEL VIRTUAL
  ! ***
  ! ***  TEMPERATURE, THE ACTUAL TEMPERATURE AND THE ADIABATIC
  ! ***
  ! ***                   LIQUID WATER CONTENT
  ! ***


  ! jyg1
  ! make sure that "Cloud base" seen by TLIFT is actually the
  ! fisrt level where adiabatic ascent is saturated
  IF (plcl>p(icb)) THEN
    ! sb        CALL TLIFT(P,T,RR,RS,GZ,PLCL,ICB,TVP,TP,CLW,ND,NL)
    CALL tlift(p, t, rr, rs, gz, plcl, icb, nk, tvp, tp, clw, nd, nl, &
      dtvpdt1, dtvpdq1)
  ELSE
    ! sb        CALL TLIFT(P,T,RR,RS,GZ,PLCL,ICB+1,TVP,TP,CLW,ND,NL)
    CALL tlift(p, t, rr, rs, gz, plcl, icb+1, nk, tvp, tp, clw, nd, nl, &
      dtvpdt1, dtvpdq1)
  END IF
  ! jyg2

  ! *****************************************************************************
  ! ***     SORTIE DE LA TEMPERATURE DE L ASCENDANCE NON DILUE
  ! *****************************************************************************
  DO i = 1, nd
    tps(i) = tp(i)
  END DO


  ! *****************************************************************************


  ! ***  SET THE PRECIPITATION EFFICIENCIES AND THE FRACTION OF   ***
  ! ***          PRECIPITATION FALLING OUTSIDE OF CLOUD           ***
  ! ***      THESE MAY BE FUNCTIONS OF TP(I), P(I) AND CLW(I)     ***

  DO i = 1, nl
    pden = ptcrit - pbcrit

    ! jyg
    ! cc         EP(I)=(P(ICB)-P(I)-PBCRIT)/PDEN
    ! sb         EP(I)=(PLCL-P(I)-PBCRIT)/PDEN
    ep(i) = (plcl-p(i)-pbcrit)/pden*epmax ! sb

    ep(i) = amax1(ep(i), 0.0)
    ! sb         EP(I)=AMIN1(EP(I),1.0)
    ep(i) = amin1(ep(i), epmax) ! sb
    sigp(i) = spfac

    ! ***       CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL     ***
    ! ***                    VIRTUAL TEMPERATURE                    ***

    tv(i) = t(i)*(1.+rr(i)/eps-rr(i))
    ! cd1
    ! . Keep all liquid water in lifted parcel (-> adiabatic CAPE)

    ! cc    TVP(I)=TVP(I)-TP(I)*(RR(1)-EP(I)*CLW(I))
    ! !!! sb         TVP(I)=TVP(I)-TP(I)*RR(1) ! calcule dans tlift
    ! cd2

    ! ***       Calculate first estimate of buoyancy

    buoy(i) = tvp(i) - tv(i)
  END DO

  ! ***   Set Cloud Base Buoyancy at (Plcl+DPbase) level buoyancy

  dpbase = -40. !That is 400m above LCL
  pbase = plcl + dpbase
  tvpbase = tvp(icb)*(pbase-p(icb+1))/(p(icb)-p(icb+1)) + &
    tvp(icb+1)*(p(icb)-pbase)/(p(icb)-p(icb+1))
  tvbase = tv(icb)*(pbase-p(icb+1))/(p(icb)-p(icb+1)) + &
    tv(icb+1)*(p(icb)-pbase)/(p(icb)-p(icb+1))

  ! test sb:
  ! @      write(*,*) '++++++++++++++++++++++++++++++++++++++++'
  ! @      write(*,*) 'plcl,dpbas,tvpbas,tvbas,tvp(icb),tvp(icb1)
  ! @     :             ,tv(icb),tv(icb1)'
  ! @      write(*,*) plcl,dpbase,tvpbase,tvbase,tvp(icb)
  ! @     L          ,tvp(icb+1),tv(icb),tv(icb+1)
  ! @      write(*,*) '++++++++++++++++++++++++++++++++++++++++'
  ! fin test sb
  buoybase = tvpbase - tvbase

  ! C       Set buoyancy = BUOYBASE for all levels below BASE.
  ! C       For safety, set : BUOY(ICB) = BUOYBASE
  DO i = icb, nl
    IF (p(i)>=pbase) THEN
      buoy(i) = buoybase
    END IF
  END DO
  buoy(icb) = buoybase

  ! sb3d      print *,'buoybase,tvp_tv,tvpbase,tvbase,pbase,plcl'
  ! sb3d     $,        buoybase,tvp(icb)-tv(icb),tvpbase,tvbase,pbase,plcl
  ! sb3d      print *,'TVP ',(tvp(i),i=1,nl)
  ! sb3d      print *,'TV ',(tv(i),i=1,nl)
  ! sb3d      print *, 'P ',(p(i),i=1,nl)
  ! sb3d      print *,'ICB ',icb
  ! test sb:
  ! @      write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
  ! @      write(*,*) 'icb,icbs,inb,buoybase:'
  ! @      write(*,*) icb,icb+1,inb,buoybase
  ! @      write(*,*) 'k,tvp,tv,tp,buoy,ep: '
  ! @      do k=1,nl
  ! @      write(*,*) k,tvp(k),tv(k),tp(k),buoy(k),ep(k)
  ! @      enddo
  ! @      write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
  ! fin test sb



  ! ***   MAKE SURE THAT COLUMN IS DRY ADIABATIC BETWEEN THE SURFACE  ***
  ! ***    AND CLOUD BASE, AND THAT LIFTED AIR IS POSITIVELY BUOYANT  ***
  ! ***                         AT CLOUD BASE                         ***
  ! ***       IF NOT, RETURN TO CALLING PROGRAM AFTER RESETTING       ***
  ! ***                        SIG(I) AND W0(I)                       ***

  ! jyg
  ! CC      TDIF=TVP(ICB)-TV(ICB)
  tdif = buoy(icb)
  ath1 = th(1)
  ! jyg
  ! CC      ATH=TH(ICB-1)-1.0
  ath = th(icb-1) - 5.0
  ! ATH=0.                          ! ajout sb
  ! IF (ICB.GT.1) ATH=TH(ICB-1)-5.0 ! modif sb
  IF (tdif<dtcrit .OR. ath>ath1) THEN
    DO i = 1, nl
      sig(i) = beta*sig(i) - 2.*alpha*tdif*tdif
      sig(i) = amax1(sig(i), 0.0)
      w0(i) = beta*w0(i)
    END DO
    iflag = 0
    RETURN
  END IF



  ! ***  IF THIS POINT IS REACHED, MOIST CONVECTIVE ADJUSTMENT IS NECESSARY
  ! ***
  ! ***        NOW INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS
  ! ***

  DO i = 1, nl
    hp(i) = h(i)
    wt(i) = 0.001
    rp(i) = rr(i)
    up(i) = u(i)
    vp(i) = v(i)
    DO j = 1, ntra
      trap(i, j) = tra(i, j)
    END DO
    nent(i) = 0
    water(i) = 0.0
    evap(i) = 0.0
    b(i) = 0.0
    mp(i) = 0.0
    m(i) = 0.0
    lv(i) = alv0 - cpvmcl*(t(i)-273.15)
    lvcp(i) = lv(i)/cpn(i)
    DO j = 1, nl
      qent(i, j) = rr(j)
      elij(i, j) = 0.0
      ment(i, j) = 0.0
      sij(i, j) = 0.0
      uent(i, j) = u(j)
      vent(i, j) = v(j)
      DO k = 1, ntra
        traent(i, j, k) = tra(j, k)
      END DO
    END DO
  END DO

  delti = 1.0/delt

  ! ***  FIND THE FIRST MODEL LEVEL (INB) ABOVE THE PARCEL'S       ***
  ! ***                LEVEL OF NEUTRAL BUOYANCY                   ***

  inb = nl - 1
  DO i = icb, nl - 1
    ! jyg
    ! CC         IF((TVP(I)-TV(I)).LT.DTCRIT)THEN
    IF (buoy(i)<dtovsh) THEN
      inb = min(inb, i)
    END IF
  END DO





  ! ***          RESET SIG(I) AND W0(I) FOR I>INB AND I<ICB       ***

  IF (inb<(nl-1)) THEN
    DO i = inb + 1, nl - 1
      ! jyg
      ! CC            SIG(I)=BETA*SIG(I)-2.0E-4*ALPHA*(TV(INB)-TVP(INB))*
      ! CC     1              ABS(TV(INB)-TVP(INB))
      sig(i) = beta*sig(i) + 2.*alpha*buoy(inb)*abs(buoy(inb))
      sig(i) = amax1(sig(i), 0.0)
      w0(i) = beta*w0(i)
    END DO
  END IF
  DO i = 1, icb
    ! jyg
    ! CC         SIG(I)=BETA*SIG(I)-2.0E-4*ALPHA*(TV(ICB)-TVP(ICB))*
    ! CC     1           (TV(ICB)-TVP(ICB))
    sig(i) = beta*sig(i) - 2.*alpha*buoy(icb)*buoy(icb)
    sig(i) = amax1(sig(i), 0.0)
    w0(i) = beta*w0(i)
  END DO

  ! ***    RESET FRACTIONAL AREAS OF UPDRAFTS AND W0 AT INITIAL TIME    ***
  ! ***           AND AFTER 10 TIME STEPS OF NO CONVECTION              ***


  IF (sig(nd)<1.5 .OR. sig(nd)>12.0) THEN
    DO i = 1, nl - 1
      sig(i) = 0.0
      w0(i) = 0.0
    END DO
  END IF

  ! ***   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL   ***

  DO i = icb, inb
    hp(i) = h(1) + (lv(i)+(cpd-cpv)*t(i))*ep(i)*clw(i)
  END DO

  ! ***  CALCULATE CONVECTIVE AVAILABLE POTENTIAL ENERGY (CAPE),  ***
  ! ***     VERTICAL VELOCITY (W), FRACTIONAL AREA COVERED BY     ***
  ! ***     UNDILUTE UPDRAFT (SIG),  AND UPDRAFT MASS FLUX (M)    ***

  cape = 0.0

  DO i = icb + 1, inb
    ! jyg1
    ! CC         CAPE=CAPE+RD*(TVP(I-1)-TV(I-1))*(PH(I-1)-PH(I))/P(I-1)
    ! CC         DCAPE=RD*BUOY(I-1)*(PH(I-1)-PH(I))/P(I-1)
    ! CC         DLNP=(PH(I-1)-PH(I))/P(I-1)
    ! The interval on which CAPE is computed starts at PBASE :
    deltap = min(pbase, ph(i-1)) - min(pbase, ph(i))
    cape = cape + rd*buoy(i-1)*deltap/p(i-1)
    dcape = rd*buoy(i-1)*deltap/p(i-1)
    dlnp = deltap/p(i-1)

    cape = amax1(0.0, cape)

    sigold = sig(i)
    dtmin = 100.0
    DO j = icb, i - 1
      ! jyg
      ! CC            DTMIN=AMIN1(DTMIN,(TVP(J)-TV(J)))
      dtmin = amin1(dtmin, buoy(j))
    END DO
    ! sb3d     print *, 'DTMIN, BETA, ALPHA, SIG = ',DTMIN,BETA,ALPHA,SIG(I)
    sig(i) = beta*sig(i) + alpha*dtmin*abs(dtmin)
    sig(i) = amax1(sig(i), 0.0)
    sig(i) = amin1(sig(i), 0.01)
    fac = amin1(((dtcrit-dtmin)/dtcrit), 1.0)
    ! jyg
    ! C    Essais de reduction de la vitesse
    ! C         FAC = FAC*.5

    w = (1.-beta)*fac*sqrt(cape) + beta*w0(i)
    amu = 0.5*(sig(i)+sigold)*w
    m(i) = amu*0.007*p(i)*(ph(i)-ph(i+1))/tv(i)

    w0(i) = w
  END DO
  w0(icb) = 0.5*w0(icb+1)
  m(icb) = 0.5*m(icb+1)*(ph(icb)-ph(icb+1))/(ph(icb+1)-ph(icb+2))
  sig(icb) = sig(icb+1)
  sig(icb-1) = sig(icb)
  ! jyg1
  ! sb3d      print *, 'Cloud base, c. top, CAPE',ICB,INB,cape
  ! sb3d      print *, 'SIG ',(sig(i),i=1,inb)
  ! sb3d      print *, 'W ',(w0(i),i=1,inb)
  ! sb3d      print *, 'M ',(m(i), i=1,inb)
  ! sb3d      print *, 'Dt1 ',(tvp(i)-tv(i),i=1,inb)
  ! sb3d      print *, 'Dt_vrai ',(buoy(i),i=1,inb)
  ! jyg2

  ! ***  CALCULATE ENTRAINED AIR MASS FLUX (MENT), TOTAL WATER MIXING  ***
  ! ***     RATIO (QENT), TOTAL CONDENSED WATER (ELIJ), AND MIXING     ***
  ! ***                        FRACTION (SIJ)                          ***



  DO i = icb, inb
    rti = rr(1) - ep(i)*clw(i)
    DO j = icb - 1, inb
      bf2 = 1. + lv(j)*lv(j)*rs(j)/(rv*t(j)*t(j)*cpd)
      anum = h(j) - hp(i) + (cpv-cpd)*t(j)*(rti-rr(j))
      denom = h(i) - hp(i) + (cpd-cpv)*(rr(i)-rti)*t(j)
      dei = denom
      IF (abs(dei)<0.01) dei = 0.01
      sij(i, j) = anum/dei
      sij(i, i) = 1.0
      altem = sij(i, j)*rr(i) + (1.-sij(i,j))*rti - rs(j)
      altem = altem/bf2
      cwat = clw(j)*(1.-ep(j))
      stemp = sij(i, j)
      IF ((stemp<0.0 .OR. stemp>1.0 .OR. altem>cwat) .AND. j>i) THEN
        anum = anum - lv(j)*(rti-rs(j)-cwat*bf2)
        denom = denom + lv(j)*(rr(i)-rti)
        IF (abs(denom)<0.01) denom = 0.01
        sij(i, j) = anum/denom
        altem = sij(i, j)*rr(i) + (1.-sij(i,j))*rti - rs(j)
        altem = altem - (bf2-1.)*cwat
      END IF


      IF (sij(i,j)>0.0 .AND. sij(i,j)<0.95) THEN
        qent(i, j) = sij(i, j)*rr(i) + (1.-sij(i,j))*rti
        uent(i, j) = sij(i, j)*u(i) + (1.-sij(i,j))*u(nk)
        vent(i, j) = sij(i, j)*v(i) + (1.-sij(i,j))*v(nk)
        DO k = 1, ntra
          traent(i, j, k) = sij(i, j)*tra(i, k) + (1.-sij(i,j))*tra(nk, k)
        END DO
        elij(i, j) = altem
        elij(i, j) = amax1(0.0, elij(i,j))
        ment(i, j) = m(i)/(1.-sij(i,j))
        nent(i) = nent(i) + 1
      END IF
      sij(i, j) = amax1(0.0, sij(i,j))
      sij(i, j) = amin1(1.0, sij(i,j))
    END DO

    ! ***   IF NO AIR CAN ENTRAIN AT LEVEL I ASSUME THAT UPDRAFT DETRAINS
    ! ***
    ! ***   AT THAT LEVEL AND CALCULATE DETRAINED AIR FLUX AND PROPERTIES
    ! ***

    IF (nent(i)==0) THEN
      ment(i, i) = m(i)
      qent(i, i) = rr(nk) - ep(i)*clw(i)
      uent(i, i) = u(nk)
      vent(i, i) = v(nk)
      DO j = 1, ntra
        traent(i, i, j) = tra(nk, j)
      END DO
      elij(i, i) = clw(i)
      sij(i, i) = 1.0
    END IF

    DO j = icb - 1, inb
      sigij(i, j) = sij(i, j)
    END DO

  END DO

  ! ***  NORMALIZE ENTRAINED AIR MASS FLUXES TO REPRESENT EQUAL  ***
  ! ***              PROBABILITIES OF MIXING                     ***


  DO i = icb, inb
    IF (nent(i)/=0) THEN
      qp = rr(1) - ep(i)*clw(i)
      anum = h(i) - hp(i) - lv(i)*(qp-rs(i)) + (cpv-cpd)*t(i)*(qp-rr(i))
      denom = h(i) - hp(i) + lv(i)*(rr(i)-qp) + (cpd-cpv)*t(i)*(rr(i)-qp)
      IF (abs(denom)<0.01) denom = 0.01
      scrit = anum/denom
      alt = qp - rs(i) + scrit*(rr(i)-qp)
      IF (scrit<=0.0 .OR. alt<=0.0) scrit = 1.0
      smax = 0.0
      asij = 0.0
      DO j = inb, icb - 1, -1
        IF (sij(i,j)>1.0E-16 .AND. sij(i,j)<0.95) THEN
          wgh = 1.0
          IF (j>i) THEN
            sjmax = amax1(sij(i,j+1), smax)
            sjmax = amin1(sjmax, scrit)
            smax = amax1(sij(i,j), smax)
            sjmin = amax1(sij(i,j-1), smax)
            sjmin = amin1(sjmin, scrit)
            IF (sij(i,j)<(smax-1.0E-16)) wgh = 0.0
            smid = amin1(sij(i,j), scrit)
          ELSE
            sjmax = amax1(sij(i,j+1), scrit)
            smid = amax1(sij(i,j), scrit)
            sjmin = 0.0
            IF (j>1) sjmin = sij(i, j-1)
            sjmin = amax1(sjmin, scrit)
          END IF
          delp = abs(sjmax-smid)
          delm = abs(sjmin-smid)
          asij = asij + wgh*(delp+delm)
          ment(i, j) = ment(i, j)*(delp+delm)*wgh
        END IF
      END DO
      asij = amax1(1.0E-16, asij)
      asij = 1.0/asij
      DO j = icb - 1, inb
        ment(i, j) = ment(i, j)*asij
      END DO
      asum = 0.0
      bsum = 0.0
      DO j = icb - 1, inb
        asum = asum + ment(i, j)
        ment(i, j) = ment(i, j)*sig(j)
        bsum = bsum + ment(i, j)
      END DO
      bsum = amax1(bsum, 1.0E-16)
      bsum = 1.0/bsum
      DO j = icb - 1, inb
        ment(i, j) = ment(i, j)*asum*bsum
      END DO
      csum = 0.0
      DO j = icb - 1, inb
        csum = csum + ment(i, j)
      END DO

      IF (csum<m(i)) THEN
        nent(i) = 0
        ment(i, i) = m(i)
        qent(i, i) = rr(1) - ep(i)*clw(i)
        uent(i, i) = u(nk)
        vent(i, i) = v(nk)
        DO j = 1, ntra
          traent(i, i, j) = tra(nk, j)
        END DO
        elij(i, i) = clw(i)
        sij(i, i) = 1.0
      END IF
    END IF
  END DO



  ! **************************************************************
  ! *       CALCUL DES MENTS(I,J) ET DES QENTS(I,J)
  ! *************************************************************

  DO im = 1, nd
    DO jm = 1, nd

      qents(im, jm) = qent(im, jm)
      ments(im, jm) = ment(im, jm)
    END DO
  END DO

  ! **********************************************************
  ! --- test sb:
  ! @       write(*,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
  ! @       write(*,*) 'inb,m(inb),ment(inb,inb),sigij(inb,inb):'
  ! @       write(*,*) inb,m(inb),ment(inb,inb),sigij(inb,inb)
  ! @       write(*,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
  ! ---





  ! ***  CHECK WHETHER EP(INB)=0, IF SO, SKIP PRECIPITATING    ***
  ! ***             DOWNDRAFT CALCULATION                      ***

  IF (ep(inb)<0.0001) GO TO 405

  ! ***  INTEGRATE LIQUID WATER EQUATION TO FIND CONDENSED WATER   ***
  ! ***                AND CONDENSED WATER FLUX                    ***

  wflux = 0.0
  tinv = 1./3.

  ! ***                    BEGIN DOWNDRAFT LOOP                    ***

  DO i = inb, 1, -1

    ! ***              CALCULATE DETRAINED PRECIPITATION             ***



    wdtrain = 10.0*ep(i)*m(i)*clw(i)
    IF (i>1) THEN
      DO j = 1, i - 1
        awat = elij(j, i) - (1.-ep(i))*clw(i)
        awat = amax1(awat, 0.0)
        wdtrain = wdtrain + 10.0*awat*ment(j, i)
      END DO
    END IF

    ! ***    FIND RAIN WATER AND EVAPORATION USING PROVISIONAL   ***
    ! ***              ESTIMATES OF RP(I)AND RP(I-1)             ***



    wt(i) = 45.0
    IF (i<inb) THEN
      rp(i) = rp(i+1) + (cpd*(t(i+1)-t(i))+gz(i+1)-gz(i))/lv(i)
      rp(i) = 0.5*(rp(i)+rr(i))
    END IF
    rp(i) = amax1(rp(i), 0.0)
    rp(i) = amin1(rp(i), rs(i))
    rp(inb) = rr(inb)
    IF (i==1) THEN
      afac = p(1)*(rs(1)-rp(1))/(1.0E4+2000.0*p(1)*rs(1))
    ELSE
      rp(i-1) = rp(i) + (cpd*(t(i)-t(i-1))+gz(i)-gz(i-1))/lv(i)
      rp(i-1) = 0.5*(rp(i-1)+rr(i-1))
      rp(i-1) = amin1(rp(i-1), rs(i-1))
      rp(i-1) = amax1(rp(i-1), 0.0)
      afac1 = p(i)*(rs(i)-rp(i))/(1.0E4+2000.0*p(i)*rs(i))
      afac2 = p(i-1)*(rs(i-1)-rp(i-1))/(1.0E4+2000.0*p(i-1)*rs(i-1))
      afac = 0.5*(afac1+afac2)
    END IF
    IF (i==inb) afac = 0.0
    afac = amax1(afac, 0.0)
    bfac = 1./(sigd*wt(i))

    ! jyg1
    ! CC        SIGT=1.0
    ! CC        IF(I.GE.ICB)SIGT=SIGP(I)
    ! Prise en compte de la variation progressive de SIGT dans
    ! les couches ICB et ICB-1:
    ! Pour PLCL<PH(I+1), PR1=0 & PR2=1
    ! Pour PLCL>PH(I),   PR1=1 & PR2=0
    ! Pour PH(I+1)<PLCL<PH(I), PR1 est la proportion a cheval
    ! sur le nuage, et PR2 est la proportion sous la base du
    ! nuage.
    pr1 = (plcl-ph(i+1))/(ph(i)-ph(i+1))
    pr1 = max(0., min(1.,pr1))
    pr2 = (ph(i)-plcl)/(ph(i)-ph(i+1))
    pr2 = max(0., min(1.,pr2))
    sigt = sigp(i)*pr1 + pr2
    ! sb3d         print *,'i,sigt,pr1,pr2', i,sigt,pr1,pr2
    ! jyg2

    b6 = bfac*50.*sigd*(ph(i)-ph(i+1))*sigt*afac
    c6 = water(i+1) + bfac*wdtrain - 50.*sigd*bfac*(ph(i)-ph(i+1))*evap(i+1)
    IF (c6>0.0) THEN
      revap = 0.5*(-b6+sqrt(b6*b6+4.*c6))
      evap(i) = sigt*afac*revap
      water(i) = revap*revap
    ELSE
      evap(i) = -evap(i+1) + 0.02*(wdtrain+sigd*wt(i)*water(i+1))/(sigd*(ph(i &
        )-ph(i+1)))
    END IF



    ! ***  CALCULATE PRECIPITATING DOWNDRAFT MASS FLUX UNDER     ***
    ! ***              HYDROSTATIC APPROXIMATION                 ***

    IF (i==1) GO TO 360
    tevap = amax1(0.0, evap(i))
    delth = amax1(0.001, (th(i)-th(i-1)))
    mp(i) = 10.*lvcp(i)*sigd*tevap*(p(i-1)-p(i))/delth

    ! ***           IF HYDROSTATIC ASSUMPTION FAILS,             ***
    ! ***   SOLVE CUBIC DIFFERENCE EQUATION FOR DOWNDRAFT THETA  ***
    ! ***  AND MASS FLUX FROM TWO SIMULTANEOUS DIFFERENTIAL EQNS ***

    amfac = sigd*sigd*70.0*ph(i)*(p(i-1)-p(i))*(th(i)-th(i-1))/(tv(i)*th(i))
    amp2 = abs(mp(i+1)*mp(i+1)-mp(i)*mp(i))
    IF (amp2>(0.1*amfac)) THEN
      xf = 100.0*sigd*sigd*sigd*(ph(i)-ph(i+1))
      tf = b(i) - 5.0*(th(i)-th(i-1))*t(i)/(lvcp(i)*sigd*th(i))
      af = xf*tf + mp(i+1)*mp(i+1)*tinv
      bf = 2.*(tinv*mp(i+1))**3 + tinv*mp(i+1)*xf*tf + &
        50.*(p(i-1)-p(i))*xf*tevap
      fac2 = 1.0
      IF (bf<0.0) fac2 = -1.0
      bf = abs(bf)
      ur = 0.25*bf*bf - af*af*af*tinv*tinv*tinv
      IF (ur>=0.0) THEN
        sru = sqrt(ur)
        fac = 1.0
        IF ((0.5*bf-sru)<0.0) fac = -1.0
        mp(i) = mp(i+1)*tinv + (0.5*bf+sru)**tinv + &
          fac*(abs(0.5*bf-sru))**tinv
      ELSE
        d = atan(2.*sqrt(-ur)/(bf+1.0E-28))
        IF (fac2<0.0) d = 3.14159 - d
        mp(i) = mp(i+1)*tinv + 2.*sqrt(af*tinv)*cos(d*tinv)
      END IF
      mp(i) = amax1(0.0, mp(i))
      b(i-1) = b(i) + 100.0*(p(i-1)-p(i))*tevap/(mp(i)+sigd*0.1) - &
        10.0*(th(i)-th(i-1))*t(i)/(lvcp(i)*sigd*th(i))
      b(i-1) = amax1(b(i-1), 0.0)
    END IF



    ! ***         LIMIT MAGNITUDE OF MP(I) TO MEET CFL CONDITION      ***

    ampmax = 2.0*(ph(i)-ph(i+1))*delti
    amp2 = 2.0*(ph(i-1)-ph(i))*delti
    ampmax = amin1(ampmax, amp2)
    mp(i) = amin1(mp(i), ampmax)

    ! ***      FORCE MP TO DECREASE LINEARLY TO ZERO                 ***
    ! ***       BETWEEN CLOUD BASE AND THE SURFACE                   ***

    IF (p(i)>p(icb)) THEN
      mp(i) = mp(icb)*(p(1)-p(i))/(p(1)-p(icb))
    END IF
360 CONTINUE

    ! ***       FIND MIXING RATIO OF PRECIPITATING DOWNDRAFT     ***

    IF (i==inb) GO TO 400
    rp(i) = rr(i)
    IF (mp(i)>mp(i+1)) THEN
      rp(i) = rp(i+1)*mp(i+1) + rr(i)*(mp(i)-mp(i+1)) + &
        5.*sigd*(ph(i)-ph(i+1))*(evap(i+1)+evap(i))
      rp(i) = rp(i)/mp(i)
      up(i) = up(i+1)*mp(i+1) + u(i)*(mp(i)-mp(i+1))
      up(i) = up(i)/mp(i)
      vp(i) = vp(i+1)*mp(i+1) + v(i)*(mp(i)-mp(i+1))
      vp(i) = vp(i)/mp(i)
      DO j = 1, ntra
        trap(i, j) = trap(i+1, j)*mp(i+1) + trap(i, j)*(mp(i)-mp(i+1))
        trap(i, j) = trap(i, j)/mp(i)
      END DO
    ELSE
      IF (mp(i+1)>1.0E-16) THEN
        rp(i) = rp(i+1) + 5.0*sigd*(ph(i)-ph(i+1))*(evap(i+1)+evap(i))/mp(i+1 &
          )
        up(i) = up(i+1)
        vp(i) = vp(i+1)
        DO j = 1, ntra
          trap(i, j) = trap(i+1, j)
        END DO
      END IF
    END IF
    rp(i) = amin1(rp(i), rs(i))
    rp(i) = amax1(rp(i), 0.0)
400 END DO

  ! ***  CALCULATE SURFACE PRECIPITATION IN MM/DAY     ***

  precip = wt(1)*sigd*water(1)*8640.0

  ! sb  ***  Calculate downdraft velocity scale and surface temperature and
  ! ***
  ! sb  ***                    water vapor fluctuations
  ! ***
  ! sb		(inspire de convect 4.3)

  ! BETAD=10.0
  betad = 5.0
  wd = betad*abs(mp(icb))*0.01*rd*t(icb)/(sigd*p(icb))

405 CONTINUE

  ! ***  CALCULATE TENDENCIES OF LOWEST LEVEL POTENTIAL TEMPERATURE  ***
  ! ***                      AND MIXING RATIO                        ***

  dpinv = 1.0/(ph(1)-ph(2))
  am = 0.0
  DO k = 2, inb
    am = am + m(k)
  END DO
  IF ((0.1*dpinv*am)>=delti) iflag = 4
  ft(1) = 0.1*dpinv*am*(t(2)-t(1)+(gz(2)-gz(1))/cpn(1))
  ft(1) = ft(1) - 0.5*lvcp(1)*sigd*(evap(1)+evap(2))
  ft(1) = ft(1) - 0.09*sigd*mp(2)*t(1)*b(1)*dpinv
  ft(1) = ft(1) + 0.01*sigd*wt(1)*(cl-cpd)*water(2)*(t(2)-t(1))*dpinv/cpn(1)
  fr(1) = 0.1*mp(2)*(rp(2)-rr(1))* & ! correction bug conservation eau
  ! 1    DPINV+SIGD*0.5*(EVAP(1)+EVAP(2))
    dpinv + sigd*0.5*(evap(1)+evap(2))
  ! IM cf. SBL
  ! 1    DPINV+SIGD*EVAP(1)
  fr(1) = fr(1) + 0.1*am*(rr(2)-rr(1))*dpinv
  fu(1) = fu(1) + 0.1*dpinv*(mp(2)*(up(2)-u(1))+am*(u(2)-u(1)))
  fv(1) = fv(1) + 0.1*dpinv*(mp(2)*(vp(2)-v(1))+am*(v(2)-v(1)))
  DO j = 1, ntra
    ftra(1, j) = ftra(1, j) + 0.1*dpinv*(mp(2)*(trap(2,j)-tra(1, &
      j))+am*(tra(2,j)-tra(1,j)))
  END DO
  amde = 0.0
  DO j = 2, inb
    fr(1) = fr(1) + 0.1*dpinv*ment(j, 1)*(qent(j,1)-rr(1))
    fu(1) = fu(1) + 0.1*dpinv*ment(j, 1)*(uent(j,1)-u(1))
    fv(1) = fv(1) + 0.1*dpinv*ment(j, 1)*(vent(j,1)-v(1))
    DO k = 1, ntra
      ftra(1, k) = ftra(1, k) + 0.1*dpinv*ment(j, 1)*(traent(j,1,k)-tra(1,k))
    END DO
  END DO

  ! ***  CALCULATE TENDENCIES OF POTENTIAL TEMPERATURE AND MIXING RATIO  ***
  ! ***               AT LEVELS ABOVE THE LOWEST LEVEL                   ***

  ! ***  FIRST FIND THE NET SATURATED UPDRAFT AND DOWNDRAFT MASS FLUXES  ***
  ! ***                      THROUGH EACH LEVEL                          ***



  DO i = 2, inb
    dpinv = 1.0/(ph(i)-ph(i+1))
    cpinv = 1.0/cpn(i)
    amp1 = 0.0
    DO k = i + 1, inb + 1
      amp1 = amp1 + m(k)
    END DO
    DO k = 1, i
      DO j = i + 1, inb + 1
        amp1 = amp1 + ment(k, j)
      END DO
    END DO
    IF ((0.1*dpinv*amp1)>=delti) iflag = 4
    ad = 0.0
    DO k = 1, i - 1
      DO j = i, inb
        ad = ad + ment(j, k)
      END DO
    END DO
    ft(i) = 0.1*dpinv*(amp1*(t(i+1)-t(i)+(gz(i+1)-gz(i))*cpinv)-ad*(t(i)-t(i- &
      1)+(gz(i)-gz(i-1))*cpinv)) - 0.5*sigd*lvcp(i)*(evap(i)+evap(i+1))
    rat = cpn(i-1)*cpinv
    ft(i) = ft(i) - 0.09*sigd*(mp(i+1)*t(i)*b(i)-mp(i)*t(i-1)*rat*b(i-1))* &
      dpinv
    ft(i) = ft(i) + 0.1*dpinv*ment(i, i)*(hp(i)-h(i)+t(i)*(cpv-cpd)*(rr(i)- &
      qent(i,i)))*cpinv
    ft(i) = ft(i) + 0.01*sigd*wt(i)*(cl-cpd)*water(i+1)*(t(i+1)-t(i))*dpinv* &
      cpinv
    fr(i) = 0.1*dpinv*(amp1*(rr(i+1)-rr(i))-ad*(rr(i)-rr(i-1)))
    fu(i) = fu(i) + 0.1*dpinv*(amp1*(u(i+1)-u(i))-ad*(u(i)-u(i-1)))
    fv(i) = fv(i) + 0.1*dpinv*(amp1*(v(i+1)-v(i))-ad*(v(i)-v(i-1)))
    DO k = 1, ntra
      ftra(i, k) = ftra(i, k) + 0.1*dpinv*(amp1*(tra(i+1,k)-tra(i, &
        k))-ad*(tra(i,k)-tra(i-1,k)))
    END DO
    DO k = 1, i - 1
      awat = elij(k, i) - (1.-ep(i))*clw(i)
      awat = amax1(awat, 0.0)
      fr(i) = fr(i) + 0.1*dpinv*ment(k, i)*(qent(k,i)-awat-rr(i))
      fu(i) = fu(i) + 0.1*dpinv*ment(k, i)*(uent(k,i)-u(i))
      fv(i) = fv(i) + 0.1*dpinv*ment(k, i)*(vent(k,i)-v(i))
      ! (saturated updrafts resulting from mixing)      ! cld
      qcond(i) = qcond(i) + (elij(k,i)-awat) ! cld
      nqcond(i) = nqcond(i) + 1. ! cld
      DO j = 1, ntra
        ftra(i, j) = ftra(i, j) + 0.1*dpinv*ment(k, i)*(traent(k,i,j)-tra(i,j &
          ))
      END DO
    END DO
    DO k = i, inb
      fr(i) = fr(i) + 0.1*dpinv*ment(k, i)*(qent(k,i)-rr(i))
      fu(i) = fu(i) + 0.1*dpinv*ment(k, i)*(uent(k,i)-u(i))
      fv(i) = fv(i) + 0.1*dpinv*ment(k, i)*(vent(k,i)-v(i))
      DO j = 1, ntra
        ftra(i, j) = ftra(i, j) + 0.1*dpinv*ment(k, i)*(traent(k,i,j)-tra(i,j &
          ))
      END DO
    END DO
    fr(i) = fr(i) + 0.5*sigd*(evap(i)+evap(i+1)) + 0.1*(mp(i+1)*(rp(i+ &
      1)-rr(i))-mp(i)*(rp(i)-rr(i-1)))*dpinv
    fu(i) = fu(i) + 0.1*(mp(i+1)*(up(i+1)-u(i))-mp(i)*(up(i)-u(i-1)))*dpinv
    fv(i) = fv(i) + 0.1*(mp(i+1)*(vp(i+1)-v(i))-mp(i)*(vp(i)-v(i-1)))*dpinv
    DO j = 1, ntra
      ftra(i, j) = ftra(i, j) + 0.1*dpinv*(mp(i+1)*(trap(i+1,j)-tra(i, &
        j))-mp(i)*(trap(i,j)-trap(i-1,j)))
    END DO
    ! (saturated downdrafts resulting from mixing)    ! cld
    DO k = i + 1, inb ! cld
      qcond(i) = qcond(i) + elij(k, i) ! cld
      nqcond(i) = nqcond(i) + 1. ! cld
    END DO ! cld
    ! (particular case: no detraining level is found) ! cld
    IF (nent(i)==0) THEN ! cld
      qcond(i) = qcond(i) + (1-ep(i))*clw(i) ! cld
      nqcond(i) = nqcond(i) + 1. ! cld
    END IF ! cld
    IF (nqcond(i)/=0.) THEN ! cld
      qcond(i) = qcond(i)/nqcond(i) ! cld
    END IF ! cld
  END DO




  ! ***   MOVE THE DETRAINMENT AT LEVEL INB DOWN TO LEVEL INB-1   ***
  ! ***        IN SUCH A WAY AS TO PRESERVE THE VERTICALLY        ***
  ! ***          INTEGRATED ENTHALPY AND WATER TENDENCIES         ***

  ! test sb:
  ! @      write(*,*) '--------------------------------------------'
  ! @      write(*,*) 'inb,ft,hp,h,t,rr,qent,ment,water,waterp,wt,mp,b'
  ! @      write(*,*) inb,ft(inb),hp(inb),h(inb)
  ! @     :   ,t(inb),rr(inb),qent(inb,inb)
  ! @     :   ,ment(inb,inb),water(inb)
  ! @     :   ,water(inb+1),wt(inb),mp(inb),b(inb)
  ! @      write(*,*) '--------------------------------------------'
  ! fin test sb:

  ax = 0.1*ment(inb, inb)*(hp(inb)-h(inb)+t(inb)*(cpv-cpd)*(rr(inb)-qent(inb, &
    inb)))/(cpn(inb)*(ph(inb)-ph(inb+1)))
  ft(inb) = ft(inb) - ax
  ft(inb-1) = ft(inb-1) + ax*cpn(inb)*(ph(inb)-ph(inb+1))/(cpn(inb-1)*(ph(inb &
    -1)-ph(inb)))
  bx = 0.1*ment(inb, inb)*(qent(inb,inb)-rr(inb))/(ph(inb)-ph(inb+1))
  fr(inb) = fr(inb) - bx
  fr(inb-1) = fr(inb-1) + bx*(ph(inb)-ph(inb+1))/(ph(inb-1)-ph(inb))
  cx = 0.1*ment(inb, inb)*(uent(inb,inb)-u(inb))/(ph(inb)-ph(inb+1))
  fu(inb) = fu(inb) - cx
  fu(inb-1) = fu(inb-1) + cx*(ph(inb)-ph(inb+1))/(ph(inb-1)-ph(inb))
  dx = 0.1*ment(inb, inb)*(vent(inb,inb)-v(inb))/(ph(inb)-ph(inb+1))
  fv(inb) = fv(inb) - dx
  fv(inb-1) = fv(inb-1) + dx*(ph(inb)-ph(inb+1))/(ph(inb-1)-ph(inb))
  DO j = 1, ntra
    ex = 0.1*ment(inb, inb)*(traent(inb,inb,j)-tra(inb,j))/ &
      (ph(inb)-ph(inb+1))
    ftra(inb, j) = ftra(inb, j) - ex
    ftra(inb-1, j) = ftra(inb-1, j) + ex*(ph(inb)-ph(inb+1))/(ph(inb-1)-ph( &
      inb))
  END DO

  ! ***    HOMOGINIZE TENDENCIES BELOW CLOUD BASE    ***

  asum = 0.0
  bsum = 0.0
  csum = 0.0
  dsum = 0.0
  DO i = 1, icb - 1
    asum = asum + ft(i)*(ph(i)-ph(i+1))
    bsum = bsum + fr(i)*(lv(i)+(cl-cpd)*(t(i)-t(1)))*(ph(i)-ph(i+1))
    csum = csum + (lv(i)+(cl-cpd)*(t(i)-t(1)))*(ph(i)-ph(i+1))
    dsum = dsum + t(i)*(ph(i)-ph(i+1))/th(i)
  END DO
  DO i = 1, icb - 1
    ft(i) = asum*t(i)/(th(i)*dsum)
    fr(i) = bsum/csum
  END DO

  ! ***           RESET COUNTER AND RETURN           ***

  sig(nd) = 2.0


  DO i = 1, nd
    upwd(i) = 0.0
    dnwd(i) = 0.0
    ! sb       dnwd0(i) = - mp(i)
  END DO

  DO i = 1, nl
    dnwd0(i) = -mp(i)
  END DO
  DO i = nl + 1, nd
    dnwd0(i) = 0.
  END DO

  DO i = icb, inb
    upwd(i) = 0.0
    dnwd(i) = 0.0

    DO k = i, inb
      up1 = 0.0
      dn1 = 0.0
      DO n = 1, i - 1
        up1 = up1 + ment(n, k)
        dn1 = dn1 - ment(k, n)
      END DO
      upwd(i) = upwd(i) + m(k) + up1
      dnwd(i) = dnwd(i) + dn1
    END DO
  END DO

  ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! DETERMINATION DE LA VARIATION DE FLUX ASCENDANT ENTRE
  ! DEUX NIVEAU NON DILUE Mike
  ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  ! sb      do i=1,ND
  ! sb      Mike(i)=M(i)
  ! sb      enddo

  DO i = 1, nl
    mike(i) = m(i)
  END DO
  DO i = nl + 1, nd
    mike(i) = 0.
  END DO

  DO i = 1, nd
    ma(i) = 0
  END DO

  ! sb      do i=1,nd
  ! sb      do j=i,nd
  ! sb      Ma(i)=Ma(i)+M(j)
  ! sb      enddo
  ! sb      enddo

  DO i = 1, nl
    DO j = i, nl
      ma(i) = ma(i) + m(j)
    END DO
  END DO

  DO i = nl + 1, nd
    ma(i) = 0.
  END DO

  DO i = 1, icb - 1
    ma(i) = 0
  END DO



  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! ICB REPRESENTE DE NIVEAU OU SE TROUVE LA
  ! BASE DU NUAGE , ET INB LE TOP DU NUAGE
  ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  DO i = 1, nd
    mke(i) = upwd(i) + dnwd(i)
  END DO


  ! *** Diagnose the in-cloud mixing ratio   ***              ! cld
  ! ***           of condensed water         ***              ! cld
  ! ! cld
  DO i = 1, nd ! cld
    maa(i) = 0.0 ! cld
    wa(i) = 0.0 ! cld
    siga(i) = 0.0 ! cld
  END DO ! cld
  DO i = nk, inb ! cld
    DO k = i + 1, inb + 1 ! cld
      maa(i) = maa(i) + m(k) ! cld
    END DO ! cld
  END DO ! cld
  DO i = icb, inb - 1 ! cld
    axc(i) = 0. ! cld
    DO j = icb, i ! cld
      axc(i) = axc(i) + rd*(tvp(j)-tv(j))*(ph(j)-ph(j+1))/p(j) ! cld
    END DO ! cld
    IF (axc(i)>0.0) THEN ! cld
      wa(i) = sqrt(2.*axc(i)) ! cld
    END IF ! cld
  END DO ! cld
  DO i = 1, nl ! cld
    IF (wa(i)>0.0) &               ! cld
      siga(i) = maa(i)/wa(i)*rd*tvp(i)/p(i)/100./deltac ! cld
    siga(i) = min(siga(i), 1.0) ! cld
    qcondc(i) = siga(i)*clw(i)*(1.-ep(i)) & ! cld
      +(1.-siga(i))*qcond(i) ! cld
  END DO ! cld


  ! $$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! $$$         call writeg1d(1,klev,ma,'ma  ','ma  ')
  ! $$$          call writeg1d(1,klev,upwd,'upwd  ','upwd  ')
  ! $$$          call writeg1d(1,klev,dnwd,'dnwd  ','dnwd  ')
  ! $$$          call writeg1d(1,klev,dnwd0,'dnwd0  ','dnwd0  ')
  ! $$$          call writeg1d(1,klev,tvp,'tvp  ','tvp  ')
  ! $$$          call writeg1d(1,klev,tra(1:klev,3),'tra3  ','tra3  ')
  ! $$$          call writeg1d(1,klev,tra(1:klev,4),'tra4  ','tra4  ')
  ! $$$          call writeg1d(1,klev,tra(1:klev,5),'tra5  ','tra5  ')
  ! $$$          call writeg1d(1,klev,tra(1:klev,6),'tra6  ','tra6  ')
  ! $$$          call writeg1d(1,klev,tra(1:klev,7),'tra7  ','tra7  ')
  ! $$$          call writeg1d(1,klev,tra(1:klev,8),'tra8  ','tra8  ')
  ! $$$          call writeg1d(1,klev,tra(1:klev,9),'tra9  ','tra9  ')
  ! $$$          call writeg1d(1,klev,tra(1:klev,10),'tra10','tra10')
  ! $$$          call writeg1d(1,klev,tra(1:klev,11),'tra11','tra11')
  ! $$$          call writeg1d(1,klev,tra(1:klev,12),'tra12','tra12')
  ! $$$          call writeg1d(1,klev,tra(1:klev,13),'tra13','tra13')
  ! $$$          call writeg1d(1,klev,tra(1:klev,14),'tra14','tra14')
  ! $$$          call writeg1d(1,klev,tra(1:klev,15),'tra15','tra15')
  ! $$$          call writeg1d(1,klev,tra(1:klev,16),'tra16','tra16')
  ! $$$          call writeg1d(1,klev,tra(1:klev,17),'tra17','tra17')
  ! $$$          call writeg1d(1,klev,tra(1:klev,18),'tra18','tra18')
  ! $$$          call writeg1d(1,klev,tra(1:klev,19),'tra19','tra19')
  ! $$$          call writeg1d(1,klev,tra(1:klev,20),'tra20','tra20 ')
  ! $$$          call writeg1d(1,klev,trap(1:klev,1),'trp1','trp1')
  ! $$$          call writeg1d(1,klev,trap(1:klev,2),'trp2','trp2')
  ! $$$          call writeg1d(1,klev,trap(1:klev,3),'trp3','trp3')
  ! $$$          call writeg1d(1,klev,trap(1:klev,4),'trp4','trp4')
  ! $$$          call writeg1d(1,klev,trap(1:klev,5),'trp5','trp5')
  ! $$$          call writeg1d(1,klev,trap(1:klev,10),'trp10','trp10')
  ! $$$          call writeg1d(1,klev,trap(1:klev,12),'trp12','trp12')
  ! $$$          call writeg1d(1,klev,trap(1:klev,15),'trp15','trp15')
  ! $$$          call writeg1d(1,klev,trap(1:klev,20),'trp20','trp20')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,1),'ftr1  ','ftr1  ')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,2),'ftr2  ','ftr2  ')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,3),'ftr3  ','ftr3  ')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,4),'ftr4  ','ftr4  ')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,5),'ftr5  ','ftr5  ')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,6),'ftr6  ','ftr6  ')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,7),'ftr7  ','ftr7  ')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,8),'ftr8  ','ftr8  ')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,9),'ftr9  ','ftr9  ')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,10),'ftr10','ftr10')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,11),'ftr11','ftr11')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,12),'ftr12','ftr12')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,13),'ftr13','ftr13')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,14),'ftr14','ftr14')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,15),'ftr15','ftr15')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,16),'ftr16','ftr16')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,17),'ftr17','ftr17')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,18),'ftr18','ftr18')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,19),'ftr19','ftr19')
  ! $$$          call writeg1d(1,klev,ftra(1:klev,20),'ftr20','ftr20 ')
  ! $$$          call writeg1d(1,klev,mp,'mp  ','mp ')
  ! $$$          call writeg1d(1,klev,Mke,'Mke  ','Mke ')



  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  RETURN
END SUBROUTINE convect3
! ---------------------------------------------------------------------------
