      SUBROUTINE LWTT(PGA,PGB,PUU, PTT)
      use dimens_m
      use dimphy
      use raddim
            use raddimlw
      IMPLICIT none
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
C     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN ALL SIX SPECTRAL
C     INTERVALS.
C
C     METHOD.
C     -------
C
C          1. TRANSMISSION FUNCTION BY H2O AND UNIFORMLY MIXED GASES ARE
C     COMPUTED USING PADE APPROXIMANTS AND HORNER'S ALGORITHM.
C          2. TRANSMISSION BY O3 IS EVALUATED WITH MALKMUS'S BAND MODEL.
C          3. TRANSMISSION BY H2O CONTINUUM AND AEROSOLS FOLLOW AN
C     A SIMPLE EXPONENTIAL DECREASE WITH ABSORBER AMOUNT.
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 88-12-15
C
C-----------------------------------------------------------------------
      REAL*8 O1H, O2H
      PARAMETER (O1H=2230.)
      PARAMETER (O2H=100.)
      REAL*8 RPIALF0
      PARAMETER (RPIALF0=2.0)
C
C* ARGUMENTS:
C
      REAL*8 PUU(KDLON,NUA)
      REAL*8 PTT(KDLON,NTRA)
      REAL*8 PGA(KDLON,8,2)
      REAL*8 PGB(KDLON,8,2)
C
C* LOCAL VARIABLES:
C
      REAL*8 zz, zxd, zxn
      REAL*8 zpu, zpu10, zpu11, zpu12, zpu13
      REAL*8 zeu, zeu10, zeu11, zeu12, zeu13
      REAL*8 zx, zy, zsq1, zsq2, zvxy, zuxy
      REAL*8 zaercn, zto1, zto2, zxch4, zych4, zxn2o, zyn2o
      REAL*8 zsqn21, zodn21, zsqh42, zodh42
      REAL*8 zsqh41, zodh41, zsqn22, zodn22, zttf11, zttf12
      REAL*8 zuu11, zuu12, za11, za12
      INTEGER jl, ja
C     ------------------------------------------------------------------
C
C*         1.     HORNER'S ALGORITHM FOR H2O AND CO2 TRANSMISSION
C                 -----------------------------------------------
C
 100  CONTINUE
C
C
      DO 130 JA = 1 , 8
      DO 120 JL = 1, KDLON
      ZZ      =SQRT(PUU(JL,JA))
c     ZXD(JL,1)=PGB( JL, 1,1) + ZZ(JL, 1)*(PGB( JL, 1,2) + ZZ(JL, 1))
c     ZXN(JL,1)=PGA( JL, 1,1) + ZZ(JL, 1)*(PGA( JL, 1,2) )
c     PTT(JL,1)=ZXN(JL,1)/ZXD(JL,1)
      ZXD      =PGB( JL,JA,1) + ZZ       *(PGB( JL,JA,2) + ZZ       )
      ZXN      =PGA( JL,JA,1) + ZZ       *(PGA( JL,JA,2) )
      PTT(JL,JA)=ZXN      /ZXD
  120 CONTINUE
  130 CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.     CONTINUUM, OZONE AND AEROSOL TRANSMISSION FUNCTIONS
C                 ---------------------------------------------------
C
 200  CONTINUE
C
      DO 201 JL = 1, KDLON
      PTT(JL, 9) = PTT(JL, 8)
C
C-  CONTINUUM ABSORPTION: E- AND P-TYPE
C
      ZPU   = 0.002 * PUU(JL,10)
      ZPU10 = 112. * ZPU
      ZPU11 = 6.25 * ZPU
      ZPU12 = 5.00 * ZPU
      ZPU13 = 80.0 * ZPU
      ZEU   =  PUU(JL,11)
      ZEU10 =  12. * ZEU
      ZEU11 = 6.25 * ZEU
      ZEU12 = 5.00 * ZEU
      ZEU13 = 80.0 * ZEU
C
C-  OZONE ABSORPTION
C
      ZX = PUU(JL,12)
      ZY = PUU(JL,13)
      ZUXY = 4. * ZX * ZX / (RPIALF0 * ZY)
      ZSQ1 = SQRT(1. + O1H * ZUXY ) - 1.
      ZSQ2 = SQRT(1. + O2H * ZUXY ) - 1.
      ZVXY = RPIALF0 * ZY / (2. * ZX)
      ZAERCN = PUU(JL,17) + ZEU12 + ZPU12
      ZTO1 = EXP( - ZVXY * ZSQ1 - ZAERCN )
      ZTO2 = EXP( - ZVXY * ZSQ2 - ZAERCN )
C
C-- TRACE GASES (CH4, N2O, CFC-11, CFC-12)
C
C* CH4 IN INTERVAL 800-970 + 1110-1250 CM-1
C
c     NEXOTIC=1
c     IF (NEXOTIC.EQ.1) THEN
      ZXCH4 = PUU(JL,19)
      ZYCH4 = PUU(JL,20)
      ZUXY = 4. * ZXCH4*ZXCH4/(0.103*ZYCH4)
      ZSQH41 = SQRT(1. + 33.7 * ZUXY) - 1.
      ZVXY = 0.103 * ZYCH4 / (2. * ZXCH4)
      ZODH41 = ZVXY * ZSQH41
C
C* N2O IN INTERVAL 800-970 + 1110-1250 CM-1
C
      ZXN2O = PUU(JL,21)
      ZYN2O = PUU(JL,22)
      ZUXY = 4. * ZXN2O*ZXN2O/(0.416*ZYN2O)
      ZSQN21 = SQRT(1. + 21.3 * ZUXY) - 1.
      ZVXY = 0.416 * ZYN2O / (2. * ZXN2O)
      ZODN21 = ZVXY * ZSQN21
C
C* CH4 IN INTERVAL 1250-1450 + 1880-2820 CM-1
C
      ZUXY = 4. * ZXCH4*ZXCH4/(0.113*ZYCH4)
      ZSQH42 = SQRT(1. + 400. * ZUXY) - 1.
      ZVXY = 0.113 * ZYCH4 / (2. * ZXCH4)
      ZODH42 = ZVXY * ZSQH42
C
C* N2O IN INTERVAL 1250-1450 + 1880-2820 CM-1
C
      ZUXY = 4. * ZXN2O*ZXN2O/(0.197*ZYN2O)
      ZSQN22 = SQRT(1. + 2000. * ZUXY) - 1.
      ZVXY = 0.197 * ZYN2O / (2. * ZXN2O)
      ZODN22 = ZVXY * ZSQN22
C
C* CFC-11 IN INTERVAL 800-970 + 1110-1250 CM-1
C
      ZA11 = 2. * PUU(JL,23) * 4.404E+05
      ZTTF11 = 1. - ZA11 * 0.003225
C
C* CFC-12 IN INTERVAL 800-970 + 1110-1250 CM-1
C
      ZA12 = 2. * PUU(JL,24) * 6.7435E+05
      ZTTF12 = 1. - ZA12 * 0.003225
C
      ZUU11 = - PUU(JL,15) - ZEU10 - ZPU10
      ZUU12 = - PUU(JL,16) - ZEU11 - ZPU11 - ZODH41 - ZODN21
      PTT(JL,10) = EXP( - PUU(JL,14) )
      PTT(JL,11) = EXP( ZUU11 )
      PTT(JL,12) = EXP( ZUU12 ) * ZTTF11 * ZTTF12
      PTT(JL,13) = 0.7554 * ZTO1 + 0.2446 * ZTO2
      PTT(JL,14) = PTT(JL,10) * EXP( - ZEU13 - ZPU13 )
      PTT(JL,15) = EXP ( - PUU(JL,14) - ZODH42 - ZODN22 )
 201  CONTINUE
C
      RETURN
      END
