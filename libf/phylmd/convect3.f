!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/convect3.F,v 1.1.1.1 2004/05/19 12:53:09 lmdzadmin Exp $
!
      SUBROUTINE CONVECT3	
     *    (DTIME,EPMAX,ok_adj,
     *     T1,   R1,   RS,    U,  V,  TRA,   P,     PH,
     *     ND,       NDP1,     NL, NTRA,  DELT,  IFLAG,
     *     FT, FR, FU,  FV,  FTRA,  PRECIP,
     *     icb,inb,   upwd,dnwd,dnwd0,SIG, W0,Mike,Mke,
     *     Ma,MENTS,QENTS,TPS,TLS,SIGIJ,CAPE,TVP,PBASE,BUOYBASE,
cccc     *     DTVPDT1,DTVPDQ1,DPLCLDT,DPLCLDR)
     *     DTVPDT1,DTVPDQ1,DPLCLDT,DPLCLDR,   ! sbl
     *     FT2,FR2,FU2,FV2,WD,QCOND,QCONDC)   ! sbl
C
C    ***  THE PARAMETER NA SHOULD IN GENERAL EQUAL ND   ***
C

cFleur       Introduction des traceurs dans convect3 le 6 juin 200

      use dimens_m
      use dimphy
      use SUPHEC_M

      real, intent(in):: dtime, DELT
      PARAMETER (NA=60)

      integer NTRAC
      PARAMETER (NTRAC=nqmx-2)
      REAL DELTAC              ! cld
      PARAMETER (DELTAC=0.01)  ! cld

      INTEGER NENT(NA)
      REAL T1(ND),R1(ND),RS(ND),U(ND),V(ND),TRA(ND,NTRA)
      REAL P(ND),PH(NDP1)
      REAL FT(ND),FR(ND),FU(ND),FV(ND),FTRA(ND,NTRA)
      REAL SIG(ND),W0(ND)
      REAL UENT(NA,NA),VENT(NA,NA),TRAENT(NA,NA,NTRAC),TRATM(NA)
      REAL UP(NA),VP(NA),TRAP(NA,NTRAC)
      REAL M(NA),MP(NA),MENT(NA,NA),QENT(NA,NA),ELIJ(NA,NA)
      REAL SIJ(NA,NA),TVP(NA),TV(NA),WATER(NA)
      REAL RP(NA),EP(NA),TH(NA),WT(NA),EVAP(NA),CLW(NA)
      REAL SIGP(NA),B(NA),TP(NA),CPN(NA)
      REAL LV(NA),LVCP(NA),H(NA),HP(NA),GZ(NA)
      REAL T(NA),RR(NA)
C
      REAL FT2(ND),FR2(ND),FU2(ND),FV2(ND) ! added sbl
      REAL U1(ND),V1(ND) ! added sbl
C
      REAL BUOY(NA)     !  Lifted parcel buoyancy
      REAL DTVPDT1(ND),DTVPDQ1(ND)   ! Derivatives of parcel virtual
C                                      temperature wrt T1 and Q1
      REAL CLW_NEW(NA),QI(NA)
C
      REAL WD, BETAD ! for gust factor (sb)
      REAL QCONDC(ND)  ! interface cld param (sb)
      REAL QCOND(ND),NQCOND(NA),WA(NA),MAA(NA),SIGA(NA),AXC(NA) ! cld
C
      LOGICAL ICE_CONV,ok_adj
      PARAMETER (ICE_CONV=.TRUE.)
 
cccccccccccccccccccccccccccccccccccccccccccccc
c     declaration des variables a sortir
ccccccccccccccccccccccccccccccccccccccccccccc
      real Mke(nd)
      real Mike(nd)
      real Ma(nd)
      real TPS(ND) !temperature dans les ascendances non diluees
      real TLS(ND) !temperature potentielle
      real MENTS(nd,nd)
      real QENTS(nd,nd)
      REAL SIGIJ(KLEV,KLEV)
      REAL PBASE ! pressure at the cloud base level
      REAL BUOYBASE ! buoyancy at the cloud base level
cccccccccccccccccccccccccccccccccccccccccccccc
 
 
 
c
      real dnwd0(nd)  !  precipitation driven unsaturated downdraft flux
      real dnwd(nd), dn1  ! in-cloud saturated downdraft mass flux
      real upwd(nd), up1  ! in-cloud saturated updraft mass flux
C
C   ***         ASSIGN VALUES OF THERMODYNAMIC CONSTANTS        ***
C   ***             THESE SHOULD BE CONSISTENT WITH             ***
C   ***              THOSE USED IN CALLING PROGRAM              ***
C   ***     NOTE: THESE ARE ALSO SPECIFIED IN SUBROUTINE TLIFT  ***
C
c sb      CPD=1005.7
c sb      CPV=1870.0
c sb      CL=4190.0
c sb      CPVMCL=CL-CPV
c sb      RV=461.5
c sb      RD=287.04
c sb      EPS=RD/RV
c sb      ALV0=2.501E6
ccccccccccccccccccccccc
c constantes coherentes avec le modele du Centre Europeen
c sb      RD = 1000.0 * 1.380658E-23 * 6.0221367E+23 / 28.9644
c sb      RV = 1000.0 * 1.380658E-23 * 6.0221367E+23 / 18.0153
c sb      CPD = 3.5 * RD
c sb      CPV = 4.0 * RV
c sb      CL = 4218.0
c sb      CPVMCL=CL-CPV
c sb      EPS=RD/RV
c sb      ALV0=2.5008E+06
cccccccccccccccccccccc
c on utilise les constantes thermo du Centre Europeen: (SB)
c
c
       CPD = RCPD
       CPV = RCPV
       CL = RCW
       CPVMCL = CL-CPV
       EPS = RD/RV
       ALV0 = RLVTT
c
       NK = 1 ! origin level of the lifted parcel
c
cccccccccccccccccccccc
C
C           ***  INITIALIZE OUTPUT ARRAYS AND PARAMETERS  ***
C
      DO 5 I=1,ND
         FT(I)=0.0
         FR(I)=0.0
         FU(I)=0.0
         FV(I)=0.0

         FT2(I)=0.0
         FR2(I)=0.0
         FU2(I)=0.0
         FV2(I)=0.0

         DO 4 J=1,NTRA
          FTRA(I,J)=0.0
    4    CONTINUE

         QCONDC(I)=0.0  ! cld
         QCOND(I)=0.0   ! cld
         NQCOND(I)=0.0  ! cld

         T(I)=T1(I)
         RR(I)=R1(I)
         U1(I)=U(I) ! added sbl
         V1(I)=V(I) ! added sbl
    5 CONTINUE
      DO 7 I=1,NL
         RDCP=(RD*(1.-RR(I))+RR(I)*RV)/ (CPD*(1.-RR(I))+RR(I)*CPV)
         TH(I)=T(I)*(1000.0/P(I))**RDCP
    7 CONTINUE
C
*************************************************************
**    CALCUL DES TEMPERATURES POTENTIELLES A SORTIR
*************************************************************
      do i=1,ND
      RDCP=(RD*(1.-RR(I))+RR(I)*RV)/ (CPD*(1.-RR(I))+RR(I)*CPV)
 
      TLS(i)=T(I)*(1000.0/P(I))**RDCP
      enddo
 
 
 
 
************************************************************
 
 
      PRECIP=0.0
      WD=0.0 ! sb
      IFLAG=1
C
C   ***                    SPECIFY PARAMETERS                        ***
C   ***  PBCRIT IS THE CRITICAL CLOUD DEPTH (MB) BENEATH WHICH THE   ***
C   ***       PRECIPITATION EFFICIENCY IS ASSUMED TO BE ZERO         ***
C   ***  PTCRIT IS THE CLOUD DEPTH (MB) ABOVE WHICH THE PRECIP.      ***
C   ***            EFFICIENCY IS ASSUMED TO BE UNITY                 ***
C   ***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
C   ***  SPFAC IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE      ***
C   ***                        OF CLOUD                              ***
C   ***    ALPHA AND BETA ARE PARAMETERS THAT CONTROL THE RATE OF    ***
C   ***                 APPROACH TO QUASI-EQUILIBRIUM                ***
C   ***    (THEIR STANDARD VALUES ARE 1.0 AND 0.96, RESPECTIVELY)    ***
C   ***           (BETA MUST BE LESS THAN OR EQUAL TO 1)             ***
C   ***    DTCRIT IS THE CRITICAL BUOYANCY (K) USED TO ADJUST THE    ***
C   ***                 APPROACH TO QUASI-EQUILIBRIUM                ***
C   ***                     IT MUST BE LESS THAN 0                   ***
C
      PBCRIT=150.0
      PTCRIT=500.0
      SIGD=0.01
      SPFAC=0.15
c sb:
c     EPMAX=0.993 ! precip efficiency less than unity
c      EPMAX=1. ! precip efficiency less than unity
C
Cjyg
CCC      BETA=0.96
C  Beta is now expressed as a function of the characteristic time
C  of the convective process.
CCC        Old value : TAU = 15000.   !(for dtime = 600.s)
CCC        Other value (inducing little change) :TAU = 8000.
      TAU = 8000.
      BETA = 1.-DTIME/TAU
Cjyg
CCC      ALPHA=1.0
      ALPHA=1.5E-3*DTIME/TAU
C        Increase alpha in order to compensate W decrease
      ALPHA = ALPHA*1.5
C
Cjyg (voir CONVECT 3)
CCC      DTCRIT=-0.2
      DTCRIT=-2.
Cgf&jyg
CCC     DT pour l'overshoot.
      DTOVSH = -0.2
 
C
C           ***        INCREMENT THE COUNTER       ***
C
      SIG(ND)=SIG(ND)+1.0
      SIG(ND)=AMIN1(SIG(ND),12.1)
C
C           ***    IF NOPT IS AN INTEGER OTHER THAN 0, CONVECT     ***
C           ***     RETURNS ARRAYS T AND R THAT MAY HAVE BEEN      ***
C           ***  ALTERED BY DRY ADIABATIC ADJUSTMENT; OTHERWISE    ***
C           ***        THE RETURNED ARRAYS ARE UNALTERED.          ***
C
      NOPT=0
c!      NOPT=1 ! sbl
C
C     ***            PERFORM DRY ADIABATIC ADJUSTMENT            ***
C
C     ***  DO NOT BYPASS THIS EVEN IF THE CALLING PROGRAM HAS A  ***
C     ***                BOUNDARY LAYER SCHEME !!!               ***
C
      IF (ok_adj) THEN ! added sbl

      DO 30 I=NL-1,1,-1
         JN=0
         DO 10 J=I+1,NL
   10    IF(TH(J).LT.TH(I))JN=J
         IF(JN.EQ.0)GOTO 30
         AHM=0.0
         RM=0.0
         UM=0.0
         VM=0.0
         DO K=1,NTRA
          TRATM(K)=0.0
         END DO
         DO 15 J=I,JN
          AHM=AHM+(CPD*(1.-RR(J))+RR(J)*CPV)*T(J)*(PH(J)-PH(J+1))
          RM=RM+RR(J)*(PH(J)-PH(J+1))
          UM=UM+U(J)*(PH(J)-PH(J+1))
          VM=VM+V(J)*(PH(J)-PH(J+1))
          DO K=1,NTRA
           TRATM(K)=TRATM(K)+TRA(J,K)*(PH(J)-PH(J+1))
          END DO
   15    CONTINUE
         DPHINV=1./(PH(I)-PH(JN+1))
         RM=RM*DPHINV
         UM=UM*DPHINV
         VM=VM*DPHINV
         DO K=1,NTRA
          TRATM(K)=TRATM(K)*DPHINV
         END DO
         A2=0.0
         DO 20 J=I,JN
            RR(J)=RM
          U(J)=UM
          V(J)=VM
          DO K=1,NTRA
           TRA(J,K)=TRATM(K)
          END DO
            RDCP=(RD*(1.-RR(J))+RR(J)*RV)/ (CPD*(1.-RR(J))+RR(J)*CPV)
            X=(0.001*P(J))**RDCP
            T(J)=X
            A2=A2+(CPD*(1.-RR(J))+RR(J)*CPV)*X*(PH(J)-PH(J+1))
   20    CONTINUE
         DO 25 J=I,JN
            TH(J)=AHM/A2
            T(J)=T(J)*TH(J)
   25    CONTINUE
   30 CONTINUE

      ENDIF ! added sbl
C
C   ***   RESET INPUT ARRAYS IF ok_adj 0   ***
C
      IF (ok_adj)THEN
         DO 35 I=1,ND

           FT2(I)=(T(I)-T1(I))/DELT  ! sbl
           FR2(I)=(RR(I)-R1(I))/DELT  ! sbl
           FU2(I)=(U(I)-U1(I))/DELT  ! sbl
           FV2(I)=(V(I)-V1(I))/DELT  ! sbl

c!            T1(I)=T(I)      ! commente sbl
c!            R1(I)=RR(I)     ! commente sbl
   35    CONTINUE
      END IF
C
C  *** CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY AND STATIC ENERGY
C
      GZ(1)=0.0
      CPN(1)=CPD*(1.-RR(1))+RR(1)*CPV
      H(1)=T(1)*CPN(1)
      DO 40 I=2,NL
        TVX=T(I)*(1.+RR(I)/EPS-RR(I))
        TVY=T(I-1)*(1.+RR(I-1)/EPS-RR(I-1))
        GZ(I)=GZ(I-1)+0.5*RD*(TVX+TVY)*(P(I-1)-P(I))/PH(I)
        CPN(I)=CPD*(1.-RR(I))+CPV*RR(I)
        H(I)=T(I)*CPN(I)+GZ(I)
   40 CONTINUE
C
C   ***  CALCULATE LIFTED CONDENSATION LEVEL OF AIR AT LOWEST MODEL LEVEL ***
C   ***       (WITHIN 0.2% OF FORMULA OF BOLTON, MON. WEA. REV.,1980)     ***
C
      IF (T(1).LT.250.0.OR.RR(1).LE.0.0)THEN
         IFLAG=0
c sb3d         print*,'je suis passe par 366'
         RETURN
      END IF

cjyg1 Utilisation de la subroutine CLIFT
CC      RH=RR(1)/RS(1)
CC      CHI=T(1)/(1669.0-122.0*RH-T(1))
CC      PLCL=P(1)*(RH**CHI)
      CALL CLIFT(P(1),T(1),RR(1),RS(1),PLCL,DPLCLDT,DPLCLDR)
cjyg2
c sb3d      PRINT *,' em_plcl,p1,t1,r1,rs1,rh '
c sb3d     $        ,PLCL,P(1),T(1),RR(1),RS(1),RH
c
      IF (PLCL.LT.200.0.OR.PLCL.GE.2000.0)THEN
         IFLAG=2
         RETURN
      END IF
Cjyg1
C     Essais de modification de ICB
C
C   ***  CALCULATE FIRST LEVEL ABOVE LCL (=ICB)  ***
C
CC      ICB=NL-1
CC      DO 50 I=2,NL-1
CC         IF(P(I).LT.PLCL)THEN
CC            ICB=MIN(ICB,I)   ! ICB sup ou egal a 2
CC         END IF
CC   50 CONTINUE
CC      IF(ICB.EQ.(NL-1))THEN
CC         IFLAG=3
CC         RETURN
CC      END IF
C
C   *** CALCULATE LAYER CONTAINING LCL (=ICB)   ***
C
      ICB=NL-1
c sb      DO 50 I=2,NL-1
      DO 50 I=3,NL-1 ! modif sb pour que ICB soit sup/egal a 2
C   la modification consiste a comparer PLCL a PH et non a P:
C   ICB est defini par :  PH(ICB)<PLCL<PH(ICB-!)
         IF(PH(I).LT.PLCL)THEN
            ICB=MIN(ICB,I)
         END IF
   50 CONTINUE
      IF(ICB.EQ.(NL-1))THEN
         IFLAG=3
         RETURN
      END IF
      ICB = ICB-1 ! ICB sup ou egal a 2 
Cjyg2
C
C
 
C   *** SUBROUTINE TLIFT CALCULATES PART OF THE LIFTED PARCEL VIRTUAL      ***
C   ***  TEMPERATURE, THE ACTUAL TEMPERATURE AND THE ADIABATIC             ***
C   ***                   LIQUID WATER CONTENT                             ***
C
 
cjyg1
c make sure that "Cloud base" seen by TLIFT is actually the 
c fisrt level where adiabatic ascent is saturated 
       IF (PLCL .GT. P(ICB)) THEN
c sb        CALL TLIFT(P,T,RR,RS,GZ,PLCL,ICB,TVP,TP,CLW,ND,NL)
        CALL TLIFT(P,T,RR,RS,GZ,PLCL,ICB,NK,TVP,TP,CLW,ND,NL
     :            ,DTVPDT1,DTVPDQ1)
       ELSE
c sb        CALL TLIFT(P,T,RR,RS,GZ,PLCL,ICB+1,TVP,TP,CLW,ND,NL)
        CALL TLIFT(P,T,RR,RS,GZ,PLCL,ICB+1,NK,TVP,TP,CLW,ND,NL
     :            ,DTVPDT1,DTVPDQ1)
       ENDIF
cjyg2
 
******************************************************************************
****     SORTIE DE LA TEMPERATURE DE L ASCENDANCE NON DILUE
******************************************************************************
        do i=1,ND
        TPS(i)=TP(i)
        enddo
 
 
******************************************************************************
 
C
C   ***  SET THE PRECIPITATION EFFICIENCIES AND THE FRACTION OF   ***
C   ***          PRECIPITATION FALLING OUTSIDE OF CLOUD           ***
C   ***      THESE MAY BE FUNCTIONS OF TP(I), P(I) AND CLW(I)     ***
C
      DO 55 I=1,NL
         PDEN=PTCRIT-PBCRIT
c
cjyg
ccc         EP(I)=(P(ICB)-P(I)-PBCRIT)/PDEN
c sb         EP(I)=(PLCL-P(I)-PBCRIT)/PDEN
         EP(I)=(PLCL-P(I)-PBCRIT)/PDEN * EPMAX ! sb
c
         EP(I)=AMAX1(EP(I),0.0)
c sb         EP(I)=AMIN1(EP(I),1.0)
         EP(I)=AMIN1(EP(I),EPMAX) ! sb
         SIGP(I)=SPFAC
C
C   ***       CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL     ***
C   ***                    VIRTUAL TEMPERATURE                    ***
C
         TV(I)=T(I)*(1.+RR(I)/EPS-RR(I))
Ccd1
C    . Keep all liquid water in lifted parcel (-> adiabatic CAPE)
C
ccc    TVP(I)=TVP(I)-TP(I)*(RR(1)-EP(I)*CLW(I))
c!!! sb         TVP(I)=TVP(I)-TP(I)*RR(1) ! calcule dans tlift
Ccd2
C
C   ***       Calculate first estimate of buoyancy
C
         BUOY(I) = TVP(I) - TV(I)
   55 CONTINUE
C
C   ***   Set Cloud Base Buoyancy at (Plcl+DPbase) level buoyancy
C
      DPBASE = -40.   !That is 400m above LCL
      PBASE = PLCL + DPBASE
      TVPBASE = TVP(ICB  )*(PBASE -P(ICB+1))/(P(ICB)-P(ICB+1))
     $         +TVP(ICB+1)*(P(ICB)-PBASE   )/(P(ICB)-P(ICB+1))
      TVBASE = TV(ICB  )*(PBASE -P(ICB+1))/(P(ICB)-P(ICB+1))
     $        +TV(ICB+1)*(P(ICB)-PBASE   )/(P(ICB)-P(ICB+1))
C
c test sb:
c@      write(*,*) '++++++++++++++++++++++++++++++++++++++++'
c@      write(*,*) 'plcl,dpbas,tvpbas,tvbas,tvp(icb),tvp(icb1)
c@     :             ,tv(icb),tv(icb1)'
c@      write(*,*) plcl,dpbase,tvpbase,tvbase,tvp(icb)
c@     L          ,tvp(icb+1),tv(icb),tv(icb+1)
c@      write(*,*) '++++++++++++++++++++++++++++++++++++++++'
c fin test sb
      BUOYBASE = TVPBASE - TVBASE
C
CC       Set buoyancy = BUOYBASE for all levels below BASE.
CC       For safety, set : BUOY(ICB) = BUOYBASE
      DO I = ICB,NL
        IF (P(I) .GE. PBASE) THEN
          BUOY(I) = BUOYBASE
        ENDIF
      ENDDO
      BUOY(ICB) = BUOYBASE
C
c sb3d      print *,'buoybase,tvp_tv,tvpbase,tvbase,pbase,plcl'
c sb3d     $,        buoybase,tvp(icb)-tv(icb),tvpbase,tvbase,pbase,plcl
c sb3d      print *,'TVP ',(tvp(i),i=1,nl)
c sb3d      print *,'TV ',(tv(i),i=1,nl)
c sb3d      print *, 'P ',(p(i),i=1,nl)
c sb3d      print *,'ICB ',icb
c test sb:
c@      write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
c@      write(*,*) 'icb,icbs,inb,buoybase:'
c@      write(*,*) icb,icb+1,inb,buoybase
c@      write(*,*) 'k,tvp,tv,tp,buoy,ep: '
c@      do k=1,nl
c@      write(*,*) k,tvp(k),tv(k),tp(k),buoy(k),ep(k)
c@      enddo
c@      write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
c fin test sb


C
C   ***   MAKE SURE THAT COLUMN IS DRY ADIABATIC BETWEEN THE SURFACE  ***
C   ***    AND CLOUD BASE, AND THAT LIFTED AIR IS POSITIVELY BUOYANT  ***
C   ***                         AT CLOUD BASE                         ***
C   ***       IF NOT, RETURN TO CALLING PROGRAM AFTER RESETTING       ***
C   ***                        SIG(I) AND W0(I)                       ***
C
Cjyg
CCC      TDIF=TVP(ICB)-TV(ICB)
      TDIF = BUOY(ICB)
      ATH1=TH(1)
Cjyg
CCC      ATH=TH(ICB-1)-1.0
      ATH=TH(ICB-1)-5.0
c      ATH=0.                          ! ajout sb
c      IF (ICB.GT.1) ATH=TH(ICB-1)-5.0 ! modif sb
      IF(TDIF.LT.DTCRIT.OR.ATH.GT.ATH1)THEN
         DO 60 I=1,NL
            SIG(I)=BETA*SIG(I)-2.*ALPHA*TDIF*TDIF
            SIG(I)=AMAX1(SIG(I),0.0)
            W0(I)=BETA*W0(I)
   60    CONTINUE
         IFLAG=0
         RETURN
      END IF
C
 
 
C   ***  IF THIS POINT IS REACHED, MOIST CONVECTIVE ADJUSTMENT IS NECESSARY ***
C   ***        NOW INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS       ***
C
      DO 70 I=1,NL
         HP(I)=H(I)
         WT(I)=0.001
         RP(I)=RR(I)
         UP(I)=U(I)
         VP(I)=V(I)
         DO 71 J=1,NTRA
          TRAP(I,J)=TRA(I,J)
   71    CONTINUE
         NENT(I)=0
         WATER(I)=0.0
         EVAP(I)=0.0
         B(I)=0.0
         MP(I)=0.0
         M(I)=0.0
         LV(I)=ALV0-CPVMCL*(T(I)-273.15)
         LVCP(I)=LV(I)/CPN(I)
         DO 70 J=1,NL
            QENT(I,J)=RR(J)
            ELIJ(I,J)=0.0
            MENT(I,J)=0.0
            SIJ(I,J)=0.0
          UENT(I,J)=U(J)
          VENT(I,J)=V(J)
          DO 70 K=1,NTRA
           TRAENT(I,J,K)=TRA(J,K)
   70 CONTINUE
 
      DELTI=1.0/DELT
C
C  ***  FIND THE FIRST MODEL LEVEL (INB) ABOVE THE PARCEL'S       ***
C  ***                LEVEL OF NEUTRAL BUOYANCY                   ***
C
      INB=NL-1
      DO 80 I=ICB,NL-1
Cjyg
CCC         IF((TVP(I)-TV(I)).LT.DTCRIT)THEN
         IF(BUOY(I).LT.DTOVSH)THEN
            INB=MIN(INB,I)
         END IF
   80 CONTINUE
 
 
 
 
C
C   ***          RESET SIG(I) AND W0(I) FOR I>INB AND I<ICB       ***
C
      IF(INB.LT.(NL-1))THEN
         DO 85 I=INB+1,NL-1
Cjyg
CCC            SIG(I)=BETA*SIG(I)-2.0E-4*ALPHA*(TV(INB)-TVP(INB))*
CCC     1              ABS(TV(INB)-TVP(INB))
            SIG(I)=BETA*SIG(I)+2.*ALPHA*BUOY(INB)*
     1              ABS(BUOY(INB))
            SIG(I)=AMAX1(SIG(I),0.0)
            W0(I)=BETA*W0(I)
   85    CONTINUE
      END IF
      DO 87 I=1,ICB
Cjyg
CCC         SIG(I)=BETA*SIG(I)-2.0E-4*ALPHA*(TV(ICB)-TVP(ICB))*
CCC     1           (TV(ICB)-TVP(ICB))
         SIG(I)=BETA*SIG(I)-2.*ALPHA*BUOY(ICB)*BUOY(ICB)
         SIG(I)=AMAX1(SIG(I),0.0)
         W0(I)=BETA*W0(I)
   87 CONTINUE
C
C   ***    RESET FRACTIONAL AREAS OF UPDRAFTS AND W0 AT INITIAL TIME    ***
C   ***           AND AFTER 10 TIME STEPS OF NO CONVECTION              ***
C
 
      IF(SIG(ND).LT.1.5.OR.SIG(ND).GT.12.0)THEN
         DO 90 I=1,NL-1
            SIG(I)=0.0
            W0(I)=0.0
   90    CONTINUE
      END IF
C
C   ***   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL   ***
C
      DO 95 I=ICB,INB
         HP(I)=H(1)+(LV(I)+(CPD-CPV)*T(I))*EP(I)*CLW(I)
   95 CONTINUE
C
C   ***  CALCULATE CONVECTIVE AVAILABLE POTENTIAL ENERGY (CAPE),  ***
C   ***     VERTICAL VELOCITY (W), FRACTIONAL AREA COVERED BY     ***
C   ***     UNDILUTE UPDRAFT (SIG),  AND UPDRAFT MASS FLUX (M)    ***
C
      CAPE=0.0
C
      DO 98 I=ICB+1,INB
Cjyg1
CCC         CAPE=CAPE+RD*(TVP(I-1)-TV(I-1))*(PH(I-1)-PH(I))/P(I-1)
CCC         DCAPE=RD*BUOY(I-1)*(PH(I-1)-PH(I))/P(I-1)
CCC         DLNP=(PH(I-1)-PH(I))/P(I-1)
C          The interval on which CAPE is computed starts at PBASE :
         DELTAP = MIN(PBASE,PH(I-1))-MIN(PBASE,PH(I))
         CAPE=CAPE+RD*BUOY(I-1)*DELTAP/P(I-1)
         DCAPE=RD*BUOY(I-1)*DELTAP/P(I-1)
         DLNP=DELTAP/P(I-1)

         CAPE=AMAX1(0.0,CAPE)
C
         SIGOLD=SIG(I)
         DTMIN=100.0
         DO 97 J=ICB,I-1
Cjyg
CCC            DTMIN=AMIN1(DTMIN,(TVP(J)-TV(J)))
            DTMIN=AMIN1(DTMIN,BUOY(J))
   97    CONTINUE
c sb3d     print *, 'DTMIN, BETA, ALPHA, SIG = ',DTMIN,BETA,ALPHA,SIG(I)
         SIG(I)=BETA*SIG(I)+ALPHA*DTMIN*ABS(DTMIN)
         SIG(I)=AMAX1(SIG(I),0.0)
         SIG(I)=AMIN1(SIG(I),0.01)
         FAC=AMIN1(((DTCRIT-DTMIN)/DTCRIT),1.0)
Cjyg
CC    Essais de reduction de la vitesse
CC         FAC = FAC*.5
C
         W=(1.-BETA)*FAC*SQRT(CAPE)+BETA*W0(I)
         AMU=0.5*(SIG(I)+SIGOLD)*W
         M(I)=AMU*0.007*P(I)*(PH(I)-PH(I+1))/TV(I)

         W0(I)=W
   98 CONTINUE
      W0(ICB)=0.5*W0(ICB+1)
      M(ICB)=0.5*M(ICB+1)*(PH(ICB)-PH(ICB+1))/(PH(ICB+1)-PH(ICB+2))
      SIG(ICB)=SIG(ICB+1)
      SIG(ICB-1)=SIG(ICB)
cjyg1
c sb3d      print *, 'Cloud base, c. top, CAPE',ICB,INB,cape
c sb3d      print *, 'SIG ',(sig(i),i=1,inb)
c sb3d      print *, 'W ',(w0(i),i=1,inb)
c sb3d      print *, 'M ',(m(i), i=1,inb)
c sb3d      print *, 'Dt1 ',(tvp(i)-tv(i),i=1,inb)
c sb3d      print *, 'Dt_vrai ',(buoy(i),i=1,inb)
Cjyg2
C
C   ***  CALCULATE ENTRAINED AIR MASS FLUX (MENT), TOTAL WATER MIXING  ***
C   ***     RATIO (QENT), TOTAL CONDENSED WATER (ELIJ), AND MIXING     ***
C   ***                        FRACTION (SIJ)                          ***
C
 
 
      DO 170 I=ICB,INB
         RTI=RR(1)-EP(I)*CLW(I)
         DO 160 J=ICB-1,INB
            BF2=1.+LV(J)*LV(J)*RS(J)/(RV*T(J)*T(J)*CPD)
            ANUM=H(J)-HP(I)+(CPV-CPD)*T(J)*(RTI-RR(J))
            DENOM=H(I)-HP(I)+(CPD-CPV)*(RR(I)-RTI)*T(J)
            DEI=DENOM
            IF(ABS(DEI).LT.0.01)DEI=0.01
            SIJ(I,J)=ANUM/DEI
            SIJ(I,I)=1.0
            ALTEM=SIJ(I,J)*RR(I)+(1.-SIJ(I,J))*RTI-RS(J)
            ALTEM=ALTEM/BF2
            CWAT=CLW(J)*(1.-EP(J))
            STEMP=SIJ(I,J)
            IF((STEMP.LT.0.0.OR.STEMP.GT.1.0.OR.
     1      ALTEM.GT.CWAT).AND.J.GT.I)THEN
            ANUM=ANUM-LV(J)*(RTI-RS(J)-CWAT*BF2)
            DENOM=DENOM+LV(J)*(RR(I)-RTI)
            IF(ABS(DENOM).LT.0.01)DENOM=0.01
            SIJ(I,J)=ANUM/DENOM
            ALTEM=SIJ(I,J)*RR(I)+(1.-SIJ(I,J))*RTI-RS(J)
            ALTEM=ALTEM-(BF2-1.)*CWAT
            END IF
 
 
            IF(SIJ(I,J).GT.0.0.AND.SIJ(I,J).LT.0.95)THEN
               QENT(I,J)=SIJ(I,J)*RR(I)+(1.-SIJ(I,J))*RTI
               UENT(I,J)=SIJ(I,J)*U(I)+(1.-SIJ(I,J))*U(NK)
               VENT(I,J)=SIJ(I,J)*V(I)+(1.-SIJ(I,J))*V(NK)
               DO K=1,NTRA
               TRAENT(I,J,K)=SIJ(I,J)*TRA(I,K)+(1.-SIJ(I,J))*
     1         TRA(NK,K)
               END DO
               ELIJ(I,J)=ALTEM
               ELIJ(I,J)=AMAX1(0.0,ELIJ(I,J))
               MENT(I,J)=M(I)/(1.-SIJ(I,J))
               NENT(I)=NENT(I)+1
            END IF
            SIJ(I,J)=AMAX1(0.0,SIJ(I,J))
            SIJ(I,J)=AMIN1(1.0,SIJ(I,J))
  160    CONTINUE
C
C   ***   IF NO AIR CAN ENTRAIN AT LEVEL I ASSUME THAT UPDRAFT DETRAINS  ***
C   ***   AT THAT LEVEL AND CALCULATE DETRAINED AIR FLUX AND PROPERTIES  ***
C
         IF(NENT(I).EQ.0)THEN
            MENT(I,I)=M(I)
            QENT(I,I)=RR(NK)-EP(I)*CLW(I)
           UENT(I,I)=U(NK)
           VENT(I,I)=V(NK)
           DO J=1,NTRA
            TRAENT(I,I,J)=TRA(NK,J)
           END DO
            ELIJ(I,I)=CLW(I)
            SIJ(I,I)=1.0
         END IF
C
         DO J = ICB-1,INB
           SIGIJ(I,J) = SIJ(I,J)
         ENDDO
C	
  170 CONTINUE
C
C   ***  NORMALIZE ENTRAINED AIR MASS FLUXES TO REPRESENT EQUAL  ***
C   ***              PROBABILITIES OF MIXING                     ***
C
 
      DO 200 I=ICB,INB
      IF(NENT(I).NE.0)THEN
       QP=RR(1)-EP(I)*CLW(I)
       ANUM=H(I)-HP(I)-LV(I)*(QP-RS(I))+(CPV-CPD)*T(I)*
     1    (QP-RR(I))
       DENOM=H(I)-HP(I)+LV(I)*(RR(I)-QP)+
     1    (CPD-CPV)*T(I)*(RR(I)-QP)
       IF(ABS(DENOM).LT.0.01)DENOM=0.01
       SCRIT=ANUM/DENOM
       ALT=QP-RS(I)+SCRIT*(RR(I)-QP)
       IF(SCRIT.LE.0.0.OR.ALT.LE.0.0)SCRIT=1.0
       SMAX=0.0
       ASIJ=0.0
        DO 175 J=INB,ICB-1,-1
        IF(SIJ(I,J).GT.1.0E-16.AND.SIJ(I,J).LT.0.95)THEN
         WGH=1.0
         IF(J.GT.I)THEN
          SJMAX=AMAX1(SIJ(I,J+1),SMAX)
          SJMAX=AMIN1(SJMAX,SCRIT)
          SMAX=AMAX1(SIJ(I,J),SMAX)
          SJMIN=AMAX1(SIJ(I,J-1),SMAX)
          SJMIN=AMIN1(SJMIN,SCRIT)
          IF(SIJ(I,J).LT.(SMAX-1.0E-16))WGH=0.0
          SMID=AMIN1(SIJ(I,J),SCRIT)
         ELSE
          SJMAX=AMAX1(SIJ(I,J+1),SCRIT)
          SMID=AMAX1(SIJ(I,J),SCRIT)
          SJMIN=0.0
          IF(J.GT.1)SJMIN=SIJ(I,J-1)
          SJMIN=AMAX1(SJMIN,SCRIT)
         END IF
         DELP=ABS(SJMAX-SMID)
         DELM=ABS(SJMIN-SMID)
         ASIJ=ASIJ+WGH*(DELP+DELM)
         MENT(I,J)=MENT(I,J)*(DELP+DELM)*WGH
        END IF
  175       CONTINUE
       ASIJ=AMAX1(1.0E-16,ASIJ)
       ASIJ=1.0/ASIJ
       DO 180 J=ICB-1,INB
        MENT(I,J)=MENT(I,J)*ASIJ
  180    CONTINUE
       ASUM=0.0
       BSUM=0.0
       DO 190 J=ICB-1,INB
        ASUM=ASUM+MENT(I,J)
        MENT(I,J)=MENT(I,J)*SIG(J)
        BSUM=BSUM+MENT(I,J)
  190       CONTINUE
       BSUM=AMAX1(BSUM,1.0E-16)
       BSUM=1.0/BSUM
       DO 195 J=ICB-1,INB
        MENT(I,J)=MENT(I,J)*ASUM*BSUM	
  195       CONTINUE
       CSUM=0.0
       DO 197 J=ICB-1,INB
        CSUM=CSUM+MENT(I,J)
  197       CONTINUE
 
       IF(CSUM.LT.M(I))THEN
        NENT(I)=0
        MENT(I,I)=M(I)
        QENT(I,I)=RR(1)-EP(I)*CLW(I)
          UENT(I,I)=U(NK)
          VENT(I,I)=V(NK)
          DO J=1,NTRA
           TRAENT(I,I,J)=TRA(NK,J)
          END DO
        ELIJ(I,I)=CLW(I)
        SIJ(I,I)=1.0
       END IF
      END IF
  200      CONTINUE
 
 
 
***************************************************************
**       CALCUL DES MENTS(I,J) ET DES QENTS(I,J)
**************************************************************
 
         DO im=1,nd
         do jm=1,nd
 
         QENTS(im,jm)=QENT(im,jm)
         MENTS(im,jm)=MENT(im,jm)
         enddo
         enddo
 
***********************************************************
c--- test sb:
c@       write(*,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
c@       write(*,*) 'inb,m(inb),ment(inb,inb),sigij(inb,inb):'
c@       write(*,*) inb,m(inb),ment(inb,inb),sigij(inb,inb)
c@       write(*,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
c---

 
 
 
C
C   ***  CHECK WHETHER EP(INB)=0, IF SO, SKIP PRECIPITATING    ***
C   ***             DOWNDRAFT CALCULATION                      ***
C
        IF(EP(INB).LT.0.0001)GOTO 405
C
C   ***  INTEGRATE LIQUID WATER EQUATION TO FIND CONDENSED WATER   ***
C   ***                AND CONDENSED WATER FLUX                    ***
C
        WFLUX=0.0
        TINV=1./3.
C
C    ***                    BEGIN DOWNDRAFT LOOP                    ***
C
        DO 400 I=INB,1,-1
C
C    ***              CALCULATE DETRAINED PRECIPITATION             ***
C
 
 
        WDTRAIN=10.0*EP(I)*M(I)*CLW(I)
        IF(I.GT.1)THEN
         DO 320 J=1,I-1
       AWAT=ELIJ(J,I)-(1.-EP(I))*CLW(I)
       AWAT=AMAX1(AWAT,0.0)
  320    WDTRAIN=WDTRAIN+10.0*AWAT*MENT(J,I)
        END IF
C
C    ***    FIND RAIN WATER AND EVAPORATION USING PROVISIONAL   ***
C    ***              ESTIMATES OF RP(I)AND RP(I-1)             ***
C
 
 
        WT(I)=45.0
      IF(I.LT.INB)THEN
       RP(I)=RP(I+1)+(CPD*(T(I+1)-T(I))+GZ(I+1)-GZ(I))/LV(I)
       RP(I)=0.5*(RP(I)+RR(I))
      END IF
      RP(I)=AMAX1(RP(I),0.0)
      RP(I)=AMIN1(RP(I),RS(I))
      RP(INB)=RR(INB)
      IF(I.EQ.1)THEN
       AFAC=P(1)*(RS(1)-RP(1))/(1.0E4+2000.0*P(1)*RS(1))
      ELSE
       RP(I-1)=RP(I)+(CPD*(T(I)-T(I-1))+GZ(I)-GZ(I-1))/LV(I)
       RP(I-1)=0.5*(RP(I-1)+RR(I-1))
       RP(I-1)=AMIN1(RP(I-1),RS(I-1))
       RP(I-1)=AMAX1(RP(I-1),0.0)
       AFAC1=P(I)*(RS(I)-RP(I))/(1.0E4+2000.0*P(I)*RS(I))
       AFAC2=P(I-1)*(RS(I-1)-RP(I-1))/(1.0E4+
     1    2000.0*P(I-1)*RS(I-1))
       AFAC=0.5*(AFAC1+AFAC2)
      END IF
      IF(I.EQ.INB)AFAC=0.0
        AFAC=AMAX1(AFAC,0.0)
        BFAC=1./(SIGD*WT(I))
C
Cjyg1
CCC        SIGT=1.0
CCC        IF(I.GE.ICB)SIGT=SIGP(I)
C Prise en compte de la variation progressive de SIGT dans
C les couches ICB et ICB-1:
C 	Pour PLCL<PH(I+1), PR1=0 & PR2=1
C 	Pour PLCL>PH(I),   PR1=1 & PR2=0
C 	Pour PH(I+1)<PLCL<PH(I), PR1 est la proportion a cheval
C    sur le nuage, et PR2 est la proportion sous la base du
C    nuage.
         PR1 =(PLCL-PH(I+1))/(PH(I)-PH(I+1))
         PR1 = MAX(0.,MIN(1.,PR1))
         PR2 = (PH(I)-PLCL)/(PH(I)-PH(I+1))
         PR2 = MAX(0.,MIN(1.,PR2))
         SIGT = SIGP(I)*PR1 + PR2
c sb3d         print *,'i,sigt,pr1,pr2', i,sigt,pr1,pr2
Cjyg2
C
        B6=BFAC*50.*SIGD*(PH(I)-PH(I+1))*SIGT*AFAC
        C6=WATER(I+1)+BFAC*WDTRAIN-50.*SIGD*BFAC*
     1   (PH(I)-PH(I+1))*EVAP(I+1)
      IF(C6.GT.0.0)THEN
         REVAP=0.5*(-B6+SQRT(B6*B6+4.*C6))
         EVAP(I)=SIGT*AFAC*REVAP
         WATER(I)=REVAP*REVAP
      ELSE
       EVAP(I)=-EVAP(I+1)+0.02*(WDTRAIN+SIGD*WT(I)*
     1    WATER(I+1))/(SIGD*(PH(I)-PH(I+1)))
      END IF
 
 
C
C    ***  CALCULATE PRECIPITATING DOWNDRAFT MASS FLUX UNDER     ***
C    ***              HYDROSTATIC APPROXIMATION                 ***
C
        IF(I.EQ.1)GOTO 360
      TEVAP=AMAX1(0.0,EVAP(I))
      DELTH=AMAX1(0.001,(TH(I)-TH(I-1)))
      MP(I)=10.*LVCP(I)*SIGD*TEVAP*(P(I-1)-P(I))/DELTH
C
C    ***           IF HYDROSTATIC ASSUMPTION FAILS,             ***
C    ***   SOLVE CUBIC DIFFERENCE EQUATION FOR DOWNDRAFT THETA  ***
C    ***  AND MASS FLUX FROM TWO SIMULTANEOUS DIFFERENTIAL EQNS ***
C
      AMFAC=SIGD*SIGD*70.0*PH(I)*(P(I-1)-P(I))*
     1   (TH(I)-TH(I-1))/(TV(I)*TH(I))
      AMP2=ABS(MP(I+1)*MP(I+1)-MP(I)*MP(I))
      IF(AMP2.GT.(0.1*AMFAC))THEN
         XF=100.0*SIGD*SIGD*SIGD*(PH(I)-PH(I+1))
         TF=B(I)-5.0*(TH(I)-TH(I-1))*T(I)/(LVCP(I)*SIGD*TH(I))
         AF=XF*TF+MP(I+1)*MP(I+1)*TINV
         BF=2.*(TINV*MP(I+1))**3+TINV*MP(I+1)*XF*TF+50.*
     1    (P(I-1)-P(I))*XF*TEVAP
         FAC2=1.0
         IF(BF.LT.0.0)FAC2=-1.0
         BF=ABS(BF)
         UR=0.25*BF*BF-AF*AF*AF*TINV*TINV*TINV
         IF(UR.GE.0.0)THEN
          SRU=SQRT(UR)
          FAC=1.0
          IF((0.5*BF-SRU).LT.0.0)FAC=-1.0
          MP(I)=MP(I+1)*TINV+(0.5*BF+SRU)**TINV+
     1     FAC*(ABS(0.5*BF-SRU))**TINV
         ELSE
          D=ATAN(2.*SQRT(-UR)/(BF+1.0E-28))
          IF(FAC2.LT.0.0)D=3.14159-D
          MP(I)=MP(I+1)*TINV+2.*SQRT(AF*TINV)*COS(D*TINV)
         END IF
         MP(I)=AMAX1(0.0,MP(I))
         B(I-1)=B(I)+100.0*(P(I-1)-P(I))*TEVAP/(MP(I)+SIGD*0.1)-
     1    10.0*(TH(I)-TH(I-1))*T(I)/(LVCP(I)*SIGD*TH(I))
         B(I-1)=AMAX1(B(I-1),0.0)
      END IF
 
 
C
C   ***         LIMIT MAGNITUDE OF MP(I) TO MEET CFL CONDITION      ***
C
      AMPMAX=2.0*(PH(I)-PH(I+1))*DELTI
      AMP2=2.0*(PH(I-1)-PH(I))*DELTI
      AMPMAX=AMIN1(AMPMAX,AMP2)
      MP(I)=AMIN1(MP(I),AMPMAX)
C
C    ***      FORCE MP TO DECREASE LINEARLY TO ZERO                 ***
C    ***       BETWEEN CLOUD BASE AND THE SURFACE                   ***
C
          IF(P(I).GT.P(ICB))THEN
           MP(I)=MP(ICB)*(P(1)-P(I))/(P(1)-P(ICB))
          END IF
  360   CONTINUE
C
C    ***       FIND MIXING RATIO OF PRECIPITATING DOWNDRAFT     ***
C
        IF(I.EQ.INB)GOTO 400
      RP(I)=RR(I)
        IF(MP(I).GT.MP(I+1))THEN
        RP(I)=RP(I+1)*MP(I+1)+RR(I)*(MP(I)-MP(I+1))+
     1       5.*SIGD*(PH(I)-PH(I+1))*(EVAP(I+1)+EVAP(I))
        RP(I)=RP(I)/MP(I)
          UP(I)=UP(I+1)*MP(I+1)+U(I)*(MP(I)-MP(I+1))
         UP(I)=UP(I)/MP(I)
          VP(I)=VP(I+1)*MP(I+1)+V(I)*(MP(I)-MP(I+1))
         VP(I)=VP(I)/MP(I)
          DO J=1,NTRA
           TRAP(I,J)=TRAP(I+1,J)*MP(I+1)+
     s     TRAP(I,J)*(MP(I)-MP(I+1))
           TRAP(I,J)=TRAP(I,J)/MP(I)
          END DO
        ELSE
        IF(MP(I+1).GT.1.0E-16)THEN
           RP(I)=RP(I+1)+5.0*SIGD*(PH(I)-PH(I+1))*(EVAP(I+1)+
     1      EVAP(I))/MP(I+1)
            UP(I)=UP(I+1)
            VP(I)=VP(I+1)
            DO J=1,NTRA
             TRAP(I,J)=TRAP(I+1,J)
            END DO
        END IF
        END IF
      RP(I)=AMIN1(RP(I),RS(I))
      RP(I)=AMAX1(RP(I),0.0)
  400   CONTINUE
C
C   ***  CALCULATE SURFACE PRECIPITATION IN MM/DAY     ***
C
        PRECIP=WT(1)*SIGD*WATER(1)*8640.0

c sb  ***  Calculate downdraft velocity scale and surface temperature and  ***
c sb  ***                    water vapor fluctuations                      ***
c sb		(inspire de convect 4.3)

c       BETAD=10.0         
       BETAD=5.0         
       WD=BETAD*ABS(MP(ICB))*0.01*RD*T(ICB)/(SIGD*P(ICB))

  405   CONTINUE
C
C   ***  CALCULATE TENDENCIES OF LOWEST LEVEL POTENTIAL TEMPERATURE  ***
C   ***                      AND MIXING RATIO                        ***
C
      DPINV=1.0/(PH(1)-PH(2))
        AM=0.0
        DO 410 K=2,INB
  410   AM=AM+M(K)
      IF((0.1*DPINV*AM).GE.DELTI)IFLAG=4
      FT(1)=0.1*DPINV*AM*(T(2)-T(1)+(GZ(2)-GZ(1))/CPN(1))
        FT(1)=FT(1)-0.5*LVCP(1)*SIGD*(EVAP(1)+EVAP(2))
        FT(1)=FT(1)-0.09*SIGD*MP(2)*T(1)*B(1)*DPINV
      FT(1)=FT(1)+0.01*SIGD*WT(1)*(CL-CPD)*WATER(2)*(T(2)-
     1   T(1))*DPINV/CPN(1)
        FR(1)=0.1*MP(2)*(RP(2)-RR(1))*
Ccorrection bug conservation eau
C    1    DPINV+SIGD*0.5*(EVAP(1)+EVAP(2))
     1    DPINV+SIGD*0.5*(EVAP(1)+EVAP(2))
cIM cf. SBL
C    1    DPINV+SIGD*EVAP(1)
        FR(1)=FR(1)+0.1*AM*(RR(2)-RR(1))*DPINV
        FU(1)=FU(1)+0.1*DPINV*(MP(2)*(UP(2)-U(1))+AM*(U(2)-U(1)))
        FV(1)=FV(1)+0.1*DPINV*(MP(2)*(VP(2)-V(1))+AM*(V(2)-V(1)))
        DO J=1,NTRA
         FTRA(1,J)=FTRA(1,J)+0.1*DPINV*(MP(2)*(TRAP(2,J)-TRA(1,J))+
     1    AM*(TRA(2,J)-TRA(1,J)))
        END DO
        AMDE=0.0
        DO 415 J=2,INB
         FR(1)=FR(1)+0.1*DPINV*MENT(J,1)*(QENT(J,1)-RR(1))
         FU(1)=FU(1)+0.1*DPINV*MENT(J,1)*(UENT(J,1)-U(1))
         FV(1)=FV(1)+0.1*DPINV*MENT(J,1)*(VENT(J,1)-V(1))
         DO K=1,NTRA
          FTRA(1,K)=FTRA(1,K)+0.1*DPINV*MENT(J,1)*(TRAENT(J,1,K)-
     1     TRA(1,K))
         END DO
  415      CONTINUE
C
C   ***  CALCULATE TENDENCIES OF POTENTIAL TEMPERATURE AND MIXING RATIO  ***
C   ***               AT LEVELS ABOVE THE LOWEST LEVEL                   ***
C
C   ***  FIRST FIND THE NET SATURATED UPDRAFT AND DOWNDRAFT MASS FLUXES  ***
C   ***                      THROUGH EACH LEVEL                          ***
C
 
 
        DO 500 I=2,INB
        DPINV=1.0/(PH(I)-PH(I+1))
      CPINV=1.0/CPN(I)
        AMP1=0.0
        DO 440 K=I+1,INB+1
  440   AMP1=AMP1+M(K)
        DO 450 K=1,I
        DO 450 J=I+1,INB+1
         AMP1=AMP1+MENT(K,J)
  450   CONTINUE
      IF((0.1*DPINV*AMP1).GE.DELTI)IFLAG=4
        AD=0.0
        DO 470 K=1,I-1
        DO 470 J=I,INB
  470   AD=AD+MENT(J,K)
      FT(I)=0.1*DPINV*(AMP1*(T(I+1)-T(I)+(GZ(I+1)-GZ(I))*
     1   CPINV)-AD*(T(I)-T(I-1)+(GZ(I)-GZ(I-1))*CPINV))
     2   -0.5*SIGD*LVCP(I)*(EVAP(I)+EVAP(I+1))
      RAT=CPN(I-1)*CPINV
        FT(I)=FT(I)-0.09*SIGD*(MP(I+1)*T(I)*
     1    B(I)-MP(I)*T(I-1)*RAT*B(I-1))*DPINV
      FT(I)=FT(I)+0.1*DPINV*MENT(I,I)*(HP(I)-H(I)+
     1    T(I)*(CPV-CPD)*(RR(I)-QENT(I,I)))*CPINV
      FT(I)=FT(I)+0.01*SIGD*WT(I)*(CL-CPD)*WATER(I+1)*
     1    (T(I+1)-T(I))*DPINV*CPINV
        FR(I)=0.1*DPINV*(AMP1*(RR(I+1)-RR(I))-
     1    AD*(RR(I)-RR(I-1)))
        FU(I)=FU(I)+0.1*DPINV*(AMP1*(U(I+1)-U(I))-
     1    AD*(U(I)-U(I-1)))
        FV(I)=FV(I)+0.1*DPINV*(AMP1*(V(I+1)-V(I))-
     1    AD*(V(I)-V(I-1)))
        DO K=1,NTRA
         FTRA(I,K)=FTRA(I,K)+0.1*DPINV*(AMP1*(TRA(I+1,K)-
     1    TRA(I,K))-AD*(TRA(I,K)-TRA(I-1,K)))
        END DO
        DO 480 K=1,I-1
       AWAT=ELIJ(K,I)-(1.-EP(I))*CLW(I)
       AWAT=AMAX1(AWAT,0.0)
         FR(I)=FR(I)+0.1*DPINV*MENT(K,I)*(QENT(K,I)-AWAT
     1    -RR(I))
         FU(I)=FU(I)+0.1*DPINV*MENT(K,I)*(UENT(K,I)-U(I))
         FV(I)=FV(I)+0.1*DPINV*MENT(K,I)*(VENT(K,I)-V(I))
C (saturated updrafts resulting from mixing)      ! cld   
         QCOND(I)=QCOND(I)+(ELIJ(K,I)-AWAT)       ! cld
         NQCOND(I)=NQCOND(I)+1.                   ! cld
         DO J=1,NTRA
          FTRA(I,J)=FTRA(I,J)+0.1*DPINV*MENT(K,I)*(TRAENT(K,I,J)-
     1     TRA(I,J))
         END DO
  480   CONTINUE
      DO 490 K=I,INB
       FR(I)=FR(I)+0.1*DPINV*MENT(K,I)*(QENT(K,I)-RR(I))
         FU(I)=FU(I)+0.1*DPINV*MENT(K,I)*(UENT(K,I)-U(I))
         FV(I)=FV(I)+0.1*DPINV*MENT(K,I)*(VENT(K,I)-V(I))
         DO J=1,NTRA
          FTRA(I,J)=FTRA(I,J)+0.1*DPINV*MENT(K,I)*(TRAENT(K,I,J)-
     1     TRA(I,J))
         END DO
  490      CONTINUE
        FR(I)=FR(I)+0.5*SIGD*(EVAP(I)+EVAP(I+1))+0.1*(MP(I+1)*
     1    (RP(I+1)-RR(I))-MP(I)*(RP(I)-RR(I-1)))*DPINV
        FU(I)=FU(I)+0.1*(MP(I+1)*(UP(I+1)-U(I))-MP(I)*
     1    (UP(I)-U(I-1)))*DPINV
        FV(I)=FV(I)+0.1*(MP(I+1)*(VP(I+1)-V(I))-MP(I)*
     1    (VP(I)-V(I-1)))*DPINV
        DO J=1,NTRA
         FTRA(I,J)=FTRA(I,J)+0.1*DPINV*(MP(I+1)*(TRAP(I+1,J)-TRA(I,J))-
     1    MP(I)*(TRAP(I,J)-TRAP(I-1,J)))
        END DO
C (saturated downdrafts resulting from mixing)    ! cld
        DO K=I+1,INB                              ! cld
         QCOND(I)=QCOND(I)+ELIJ(K,I)              ! cld
         NQCOND(I)=NQCOND(I)+1.                   ! cld
        ENDDO                                     ! cld
C (particular case: no detraining level is found) ! cld
        IF (NENT(I).EQ.0) THEN                    ! cld
         QCOND(I)=QCOND(I)+(1-EP(I))*CLW(I)       ! cld
         NQCOND(I)=NQCOND(I)+1.                   ! cld
        ENDIF                                     ! cld
        IF (NQCOND(I).NE.0.) THEN                 ! cld
         QCOND(I)=QCOND(I)/NQCOND(I)              ! cld
        ENDIF                                     ! cld
  500   CONTINUE
 
 
 
C
C   ***   MOVE THE DETRAINMENT AT LEVEL INB DOWN TO LEVEL INB-1   ***
C   ***        IN SUCH A WAY AS TO PRESERVE THE VERTICALLY        ***
C   ***          INTEGRATED ENTHALPY AND WATER TENDENCIES         ***
C
c test sb:
c@      write(*,*) '--------------------------------------------'
c@      write(*,*) 'inb,ft,hp,h,t,rr,qent,ment,water,waterp,wt,mp,b'
c@      write(*,*) inb,ft(inb),hp(inb),h(inb)
c@     :   ,t(inb),rr(inb),qent(inb,inb)
c@     :   ,ment(inb,inb),water(inb)
c@     :   ,water(inb+1),wt(inb),mp(inb),b(inb)
c@      write(*,*) '--------------------------------------------'
c fin test sb:

      AX=0.1*MENT(INB,INB)*(HP(INB)-H(INB)+T(INB)*
     1    (CPV-CPD)*(RR(INB)-QENT(INB,INB)))/(CPN(INB)*
     2    (PH(INB)-PH(INB+1)))
      FT(INB)=FT(INB)-AX
      FT(INB-1)=FT(INB-1)+AX*CPN(INB)*(PH(INB)-PH(INB+1))/
     1    (CPN(INB-1)*(PH(INB-1)-PH(INB)))
      BX=0.1*MENT(INB,INB)*(QENT(INB,INB)-RR(INB))/
     1    (PH(INB)-PH(INB+1))
      FR(INB)=FR(INB)-BX
      FR(INB-1)=FR(INB-1)+BX*(PH(INB)-PH(INB+1))/
     1    (PH(INB-1)-PH(INB))
      CX=0.1*MENT(INB,INB)*(UENT(INB,INB)-U(INB))/
     1    (PH(INB)-PH(INB+1))
      FU(INB)=FU(INB)-CX
      FU(INB-1)=FU(INB-1)+CX*(PH(INB)-PH(INB+1))/
     1    (PH(INB-1)-PH(INB))
      DX=0.1*MENT(INB,INB)*(VENT(INB,INB)-V(INB))/
     1    (PH(INB)-PH(INB+1))
      FV(INB)=FV(INB)-DX
      FV(INB-1)=FV(INB-1)+DX*(PH(INB)-PH(INB+1))/
     1    (PH(INB-1)-PH(INB))
      DO J=1,NTRA
      EX=0.1*MENT(INB,INB)*(TRAENT(INB,INB,J)
     1    -TRA(INB,J))/(PH(INB)-PH(INB+1))
      FTRA(INB,J)=FTRA(INB,J)-EX
      FTRA(INB-1,J)=FTRA(INB-1,J)+EX*
     1     (PH(INB)-PH(INB+1))/(PH(INB-1)-PH(INB))
      ENDDO   
C
C   ***    HOMOGINIZE TENDENCIES BELOW CLOUD BASE    ***
C
      ASUM=0.0
      BSUM=0.0
      CSUM=0.0
        DSUM=0.0
      DO 650 I=1,ICB-1
       ASUM=ASUM+FT(I)*(PH(I)-PH(I+1))
         BSUM=BSUM+FR(I)*(LV(I)+(CL-CPD)*(T(I)-T(1)))*
     1    (PH(I)-PH(I+1))
       CSUM=CSUM+(LV(I)+(CL-CPD)*(T(I)-T(1)))*(PH(I)-PH(I+1))
       DSUM=DSUM+T(I)*(PH(I)-PH(I+1))/TH(I)
  650      CONTINUE
      DO 700 I=1,ICB-1
       FT(I)=ASUM*T(I)/(TH(I)*DSUM)
       FR(I)=BSUM/CSUM
  700      CONTINUE
C
C   ***           RESET COUNTER AND RETURN           ***
C
      SIG(ND)=2.0
c
c
      do i = 1, nd
         upwd(i) = 0.0
         dnwd(i) = 0.0
c sb       dnwd0(i) = - mp(i)
      enddo
c
      do i = 1, nl
       dnwd0(i) = - mp(i)
      enddo
      do i = nl+1, nd
       dnwd0(i) = 0.
      enddo
c
      do i = icb, inb
         upwd(i) = 0.0
         dnwd(i) = 0.0

         do k =i, inb
            up1=0.0
            dn1=0.0
            do n = 1, i-1
               up1 = up1 + ment(n,k)
               dn1 = dn1 - ment(k,n)
            enddo
            upwd(i) = upwd(i)+ m(k) + up1
            dnwd(i) = dnwd(i) + dn1
         enddo
        enddo
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        DETERMINATION DE LA VARIATION DE FLUX ASCENDANT ENTRE
C        DEUX NIVEAU NON DILUE Mike
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
c sb      do i=1,ND
c sb      Mike(i)=M(i)
c sb      enddo
 
      do i = 1, NL
       Mike(i) = M(i)
      enddo
      do i = NL+1, ND
       Mike(i) = 0.
      enddo
 
      do i=1,nd
      Ma(i)=0
      enddo
 
c sb      do i=1,nd
c sb      do j=i,nd
c sb      Ma(i)=Ma(i)+M(j)
c sb      enddo
c sb      enddo

      do i = 1, NL
      do j = i, NL
       Ma(i) = Ma(i) + M(j)
      enddo
      enddo
c
      do i = NL+1, ND
       Ma(i) = 0.
      enddo
c 
      do i=1,ICB-1
      Ma(i)=0
      enddo
 
 
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C        ICB REPRESENTE DE NIVEAU OU SE TROUVE LA
c        BASE DU NUAGE , ET INB LE TOP DU NUAGE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
       do i=1,ND
       Mke(i)=upwd(i)+dnwd(i)
       enddo

C
C   *** Diagnose the in-cloud mixing ratio   ***              ! cld
C   ***           of condensed water         ***              ! cld
C                                                             ! cld
       DO I=1,ND                                              ! cld
        MAA(I)=0.0                                            ! cld
        WA(I)=0.0                                             ! cld
        SIGA(I)=0.0                                           ! cld
       ENDDO                                                  ! cld
       DO I=NK,INB                                            ! cld
       DO K=I+1,INB+1                                         ! cld
        MAA(I)=MAA(I)+M(K)                                    ! cld
       ENDDO                                                  ! cld
       ENDDO                                                  ! cld
       DO I=ICB,INB-1                                         ! cld
        AXC(I)=0.                                             ! cld
        DO J=ICB,I                                            ! cld
         AXC(I)=AXC(I)+RD*(TVP(J)-TV(J))*(PH(J)-PH(J+1))/P(J) ! cld
        ENDDO                                                 ! cld
        IF (AXC(I).GT.0.0) THEN                               ! cld
         WA(I)=SQRT(2.*AXC(I))                                ! cld
        ENDIF                                                 ! cld
       ENDDO                                                  ! cld
       DO I=1,NL                                              ! cld
        IF (WA(I).GT.0.0)                                     ! cld
     1    SIGA(I)=MAA(I)/WA(I)*RD*TVP(I)/P(I)/100./DELTAC     ! cld
        SIGA(I) = MIN(SIGA(I),1.0)                            ! cld
        QCONDC(I)=SIGA(I)*CLW(I)*(1.-EP(I))                   ! cld
     1          + (1.-SIGA(I))*QCOND(I)                       ! cld
       ENDDO                                                  ! cld


c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$         call writeg1d(1,klev,ma,'ma  ','ma  ')
c$$$          call writeg1d(1,klev,upwd,'upwd  ','upwd  ')
c$$$          call writeg1d(1,klev,dnwd,'dnwd  ','dnwd  ')
c$$$          call writeg1d(1,klev,dnwd0,'dnwd0  ','dnwd0  ')
c$$$          call writeg1d(1,klev,tvp,'tvp  ','tvp  ')
c$$$          call writeg1d(1,klev,tra(1:klev,3),'tra3  ','tra3  ')
c$$$          call writeg1d(1,klev,tra(1:klev,4),'tra4  ','tra4  ')
c$$$          call writeg1d(1,klev,tra(1:klev,5),'tra5  ','tra5  ')
c$$$          call writeg1d(1,klev,tra(1:klev,6),'tra6  ','tra6  ')
c$$$          call writeg1d(1,klev,tra(1:klev,7),'tra7  ','tra7  ')
c$$$          call writeg1d(1,klev,tra(1:klev,8),'tra8  ','tra8  ')
c$$$          call writeg1d(1,klev,tra(1:klev,9),'tra9  ','tra9  ')
c$$$          call writeg1d(1,klev,tra(1:klev,10),'tra10','tra10')
c$$$          call writeg1d(1,klev,tra(1:klev,11),'tra11','tra11')
c$$$          call writeg1d(1,klev,tra(1:klev,12),'tra12','tra12')
c$$$          call writeg1d(1,klev,tra(1:klev,13),'tra13','tra13')
c$$$          call writeg1d(1,klev,tra(1:klev,14),'tra14','tra14')
c$$$          call writeg1d(1,klev,tra(1:klev,15),'tra15','tra15')
c$$$          call writeg1d(1,klev,tra(1:klev,16),'tra16','tra16')
c$$$          call writeg1d(1,klev,tra(1:klev,17),'tra17','tra17')
c$$$          call writeg1d(1,klev,tra(1:klev,18),'tra18','tra18')
c$$$          call writeg1d(1,klev,tra(1:klev,19),'tra19','tra19')
c$$$          call writeg1d(1,klev,tra(1:klev,20),'tra20','tra20 ')
c$$$          call writeg1d(1,klev,trap(1:klev,1),'trp1','trp1')
c$$$          call writeg1d(1,klev,trap(1:klev,2),'trp2','trp2')
c$$$          call writeg1d(1,klev,trap(1:klev,3),'trp3','trp3')
c$$$          call writeg1d(1,klev,trap(1:klev,4),'trp4','trp4')
c$$$          call writeg1d(1,klev,trap(1:klev,5),'trp5','trp5')
c$$$          call writeg1d(1,klev,trap(1:klev,10),'trp10','trp10')
c$$$          call writeg1d(1,klev,trap(1:klev,12),'trp12','trp12')
c$$$          call writeg1d(1,klev,trap(1:klev,15),'trp15','trp15')
c$$$          call writeg1d(1,klev,trap(1:klev,20),'trp20','trp20')
c$$$          call writeg1d(1,klev,ftra(1:klev,1),'ftr1  ','ftr1  ')
c$$$          call writeg1d(1,klev,ftra(1:klev,2),'ftr2  ','ftr2  ')
c$$$          call writeg1d(1,klev,ftra(1:klev,3),'ftr3  ','ftr3  ')
c$$$          call writeg1d(1,klev,ftra(1:klev,4),'ftr4  ','ftr4  ')
c$$$          call writeg1d(1,klev,ftra(1:klev,5),'ftr5  ','ftr5  ')
c$$$          call writeg1d(1,klev,ftra(1:klev,6),'ftr6  ','ftr6  ')
c$$$          call writeg1d(1,klev,ftra(1:klev,7),'ftr7  ','ftr7  ')
c$$$          call writeg1d(1,klev,ftra(1:klev,8),'ftr8  ','ftr8  ')
c$$$          call writeg1d(1,klev,ftra(1:klev,9),'ftr9  ','ftr9  ')
c$$$          call writeg1d(1,klev,ftra(1:klev,10),'ftr10','ftr10')
c$$$          call writeg1d(1,klev,ftra(1:klev,11),'ftr11','ftr11')
c$$$          call writeg1d(1,klev,ftra(1:klev,12),'ftr12','ftr12')
c$$$          call writeg1d(1,klev,ftra(1:klev,13),'ftr13','ftr13')
c$$$          call writeg1d(1,klev,ftra(1:klev,14),'ftr14','ftr14')
c$$$          call writeg1d(1,klev,ftra(1:klev,15),'ftr15','ftr15')
c$$$          call writeg1d(1,klev,ftra(1:klev,16),'ftr16','ftr16')
c$$$          call writeg1d(1,klev,ftra(1:klev,17),'ftr17','ftr17')
c$$$          call writeg1d(1,klev,ftra(1:klev,18),'ftr18','ftr18')
c$$$          call writeg1d(1,klev,ftra(1:klev,19),'ftr19','ftr19')
c$$$          call writeg1d(1,klev,ftra(1:klev,20),'ftr20','ftr20 ')
c$$$          call writeg1d(1,klev,mp,'mp  ','mp ')
c$$$          call writeg1d(1,klev,Mke,'Mke  ','Mke ')

 
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        RETURN
        END
C ---------------------------------------------------------------------------
