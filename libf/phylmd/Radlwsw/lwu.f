cIM ctes ds clesphys.h   SUBROUTINE LWU(RCO2, RCH4, RN2O, RCFC11, RCFC12,
      SUBROUTINE LWU(
     S               PAER,PDP,PPMB,PPSOL,POZ,PTAVE,PVIEW,PWV,
     S               PABCU)
      use dimens_m
      use dimphy
      use clesphys
      use YOMCST
      use raddim
      use radepsi
      use radopt
            use raddimlw
      IMPLICIT none
C
C     PURPOSE.
C     --------
C           COMPUTES ABSORBER AMOUNTS INCLUDING PRESSURE AND
C           TEMPERATURE EFFECTS
C
C     METHOD.
C     -------
C
C          1. COMPUTES THE PRESSURE AND TEMPERATURE WEIGHTED AMOUNTS OF
C     ABSORBERS.
C
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
C        ORIGINAL : 89-07-14
C        Voigt lines (loop 404 modified) - JJM & PhD - 01/96
C-----------------------------------------------------------------------
C* ARGUMENTS:
cIM ctes ds clesphys.h
c     REAL*8 RCO2
c     REAL*8 RCH4, RN2O, RCFC11, RCFC12
      REAL*8 PAER(KDLON,KFLEV,5)
      REAL*8 PDP(KDLON,KFLEV)
      REAL*8 PPMB(KDLON,KFLEV+1)
      REAL*8 PPSOL(KDLON)
      REAL*8 POZ(KDLON,KFLEV)
      REAL*8 PTAVE(KDLON,KFLEV)
      REAL*8 PVIEW(KDLON)
      REAL*8 PWV(KDLON,KFLEV)
C
      REAL*8 PABCU(KDLON,NUA,3*KFLEV+1) ! EFFECTIVE ABSORBER AMOUNTS
C
C-----------------------------------------------------------------------
C* LOCAL VARIABLES:
      REAL*8 ZABLY(KDLON,NUA,3*KFLEV+1)
      REAL*8 ZDUC(KDLON,3*KFLEV+1)
      REAL*8 ZPHIO(KDLON)
      REAL*8 ZPSC2(KDLON)
      REAL*8 ZPSC3(KDLON)
      REAL*8 ZPSH1(KDLON)
      REAL*8 ZPSH2(KDLON)
      REAL*8 ZPSH3(KDLON)
      REAL*8 ZPSH4(KDLON)
      REAL*8 ZPSH5(KDLON)
      REAL*8 ZPSH6(KDLON)
      REAL*8 ZPSIO(KDLON)
      REAL*8 ZTCON(KDLON)
      REAL*8 ZPHM6(KDLON)
      REAL*8 ZPSM6(KDLON)
      REAL*8 ZPHN6(KDLON)
      REAL*8 ZPSN6(KDLON)
      REAL*8 ZSSIG(KDLON,3*KFLEV+1)
      REAL*8 ZTAVI(KDLON)
      REAL*8 ZUAER(KDLON,Ninter)
      REAL*8 ZXOZ(KDLON)
      REAL*8 ZXWV(KDLON)
C
      INTEGER jl, jk, jkj, jkjr, jkjp, ig1
      INTEGER jki, jkip1, ja, jj
      INTEGER jkl, jkp1, jkk, jkjpn
      INTEGER jae1, jae2, jae3, jae, jjpn
      INTEGER ir, jc, jcp1
      REAL*8 zdpm, zupm, zupmh2o, zupmco2, zupmo3, zu6, zup
      REAL*8 zfppw, ztx, ztx2, zzably
      REAL*8 zcah1, zcbh1, zcah2, zcbh2, zcah3, zcbh3
      REAL*8 zcah4, zcbh4, zcah5, zcbh5, zcah6, zcbh6
      REAL*8 zcac8, zcbc8
      REAL*8 zalup, zdiff
c
      REAL*8 PVGCO2, PVGH2O, PVGO3
C
      REAL*8 R10E  ! DECIMAL/NATURAL LOG.FACTOR
      PARAMETER (R10E=0.4342945)
c
c Used Data Block:
c
      REAL*8 TREF
      SAVE TREF
      REAL*8 RT1(2)
      SAVE RT1
      REAL*8 RAER(5,5)
      SAVE RAER
      REAL*8 AT(8,3), BT(8,3)
      SAVE AT, BT
      REAL*8 OCT(4)
      SAVE OCT
      DATA TREF /250.0/
      DATA (RT1(IG1),IG1=1,2) / -0.577350269, +0.577350269 /
      DATA RAER / .038520, .037196, .040532, .054934, .038520
     1          , .12613 , .18313 , .10357 , .064106, .126130
     2          , .012579, .013649, .018652, .025181, .012579
     3          , .011890, .016142, .021105, .028908, .011890
     4          , .013792, .026810, .052203, .066338, .013792 /
      DATA (AT(1,IR),IR=1,3) /
     S 0.298199E-02,-.394023E-03,0.319566E-04 /
      DATA (BT(1,IR),IR=1,3) /
     S-0.106432E-04,0.660324E-06,0.174356E-06 /
      DATA (AT(2,IR),IR=1,3) /
     S 0.143676E-01,0.366501E-02,-.160822E-02 /
      DATA (BT(2,IR),IR=1,3) /
     S-0.553979E-04,-.101701E-04,0.920868E-05 /
      DATA (AT(3,IR),IR=1,3) /
     S 0.197861E-01,0.315541E-02,-.174547E-02 /
      DATA (BT(3,IR),IR=1,3) /
     S-0.877012E-04,0.513302E-04,0.523138E-06 /
      DATA (AT(4,IR),IR=1,3) /
     S 0.289560E-01,-.208807E-02,-.121943E-02 /
      DATA (BT(4,IR),IR=1,3) /
     S-0.165960E-03,0.157704E-03,-.146427E-04 /
      DATA (AT(5,IR),IR=1,3) /
     S 0.103800E-01,0.436296E-02,-.161431E-02 /
      DATA (BT(5,IR),IR=1,3) /
     S -.276744E-04,-.327381E-04,0.127646E-04 /
      DATA (AT(6,IR),IR=1,3) /
     S 0.868859E-02,-.972752E-03,0.000000E-00 /
      DATA (BT(6,IR),IR=1,3) /
     S -.278412E-04,-.713940E-06,0.117469E-05 /
      DATA (AT(7,IR),IR=1,3) /
     S 0.250073E-03,0.455875E-03,0.109242E-03 /
      DATA (BT(7,IR),IR=1,3) /
     S 0.199846E-05,-.216313E-05,0.175991E-06 /
      DATA (AT(8,IR),IR=1,3) /
     S 0.307423E-01,0.110879E-02,-.322172E-03 /
      DATA (BT(8,IR),IR=1,3) /
     S-0.108482E-03,0.258096E-05,-.814575E-06 /
c
      DATA OCT /-.326E-03, -.102E-05, .137E-02, -.535E-05/
C-----------------------------------------------------------------------
c
      IF (LEVOIGT) THEN
         PVGCO2= 60.
         PVGH2O= 30.
         PVGO3 =400.
      ELSE
         PVGCO2= 0.
         PVGH2O= 0.
         PVGO3 = 0.
      ENDIF
C
C
C*         2.    PRESSURE OVER GAUSS SUB-LEVELS
C                ------------------------------
C
 200  CONTINUE
C
      DO 201 JL = 1, KDLON
      ZSSIG(JL, 1 ) = PPMB(JL,1) * 100.
 201  CONTINUE
C
      DO 206 JK = 1 , KFLEV
      JKJ=(JK-1)*NG1P1+1
      JKJR = JKJ
      JKJP = JKJ + NG1P1
      DO 203 JL = 1, KDLON
      ZSSIG(JL,JKJP)=PPMB(JL,JK+1)* 100.
 203  CONTINUE
      DO 205 IG1=1,NG1
      JKJ=JKJ+1
      DO 204 JL = 1, KDLON
      ZSSIG(JL,JKJ)= (ZSSIG(JL,JKJR)+ZSSIG(JL,JKJP))*0.5
     S  + RT1(IG1) * (ZSSIG(JL,JKJP) - ZSSIG(JL,JKJR)) * 0.5
 204  CONTINUE
 205  CONTINUE
 206  CONTINUE
C
C-----------------------------------------------------------------------
C
C
C*         4.    PRESSURE THICKNESS AND MEAN PRESSURE OF SUB-LAYERS
C                --------------------------------------------------
C
 400  CONTINUE
C
      DO 402 JKI=1,3*KFLEV
      JKIP1=JKI+1
      DO 401 JL = 1, KDLON
      ZABLY(JL,5,JKI)=(ZSSIG(JL,JKI)+ZSSIG(JL,JKIP1))*0.5
      ZABLY(JL,3,JKI)=(ZSSIG(JL,JKI)-ZSSIG(JL,JKIP1))
     S                                 /(10.*RG)
 401  CONTINUE
 402  CONTINUE
C
      DO 406 JK = 1 , KFLEV
      JKP1=JK+1
      JKL = KFLEV+1 - JK
      DO 403 JL = 1, KDLON
      ZXWV(JL) = MAX (PWV(JL,JK) , ZEPSCQ )
      ZXOZ(JL) = MAX (POZ(JL,JK) / PDP(JL,JK) , ZEPSCO )
 403  CONTINUE
      JKJ=(JK-1)*NG1P1+1
      JKJPN=JKJ+NG1
      DO 405 JKK=JKJ,JKJPN
      DO 404 JL = 1, KDLON
      ZDPM = ZABLY(JL,3,JKK)
      ZUPM = ZABLY(JL,5,JKK)             * ZDPM / 101325.
      ZUPMCO2 = ( ZABLY(JL,5,JKK) + PVGCO2 ) * ZDPM / 101325.
      ZUPMH2O = ( ZABLY(JL,5,JKK) + PVGH2O ) * ZDPM / 101325.
      ZUPMO3  = ( ZABLY(JL,5,JKK) + PVGO3  ) * ZDPM / 101325.
      ZDUC(JL,JKK) = ZDPM
      ZABLY(JL,12,JKK) = ZXOZ(JL) * ZDPM
      ZABLY(JL,13,JKK) = ZXOZ(JL) * ZUPMO3
      ZU6 = ZXWV(JL) * ZUPM
      ZFPPW = 1.6078 * ZXWV(JL) / (1.+0.608*ZXWV(JL))
      ZABLY(JL,6,JKK) = ZXWV(JL) * ZUPMH2O
      ZABLY(JL,11,JKK) = ZU6 * ZFPPW
      ZABLY(JL,10,JKK) = ZU6 * (1.-ZFPPW)
      ZABLY(JL,9,JKK) = RCO2 * ZUPMCO2
      ZABLY(JL,8,JKK) = RCO2 * ZDPM
 404  CONTINUE
 405  CONTINUE
 406  CONTINUE
C
C-----------------------------------------------------------------------
C
C
C*         5.    CUMULATIVE ABSORBER AMOUNTS FROM TOP OF ATMOSPHERE
C                --------------------------------------------------
C
 500  CONTINUE
C
      DO 502 JA = 1, NUA
      DO 501 JL = 1, KDLON
      PABCU(JL,JA,3*KFLEV+1) = 0.
  501 CONTINUE
  502 CONTINUE
C
      DO 529 JK = 1 , KFLEV
      JJ=(JK-1)*NG1P1+1
      JJPN=JJ+NG1
      JKL=KFLEV+1-JK
C
C
C*         5.1  CUMULATIVE AEROSOL AMOUNTS FROM TOP OF ATMOSPHERE
C               --------------------------------------------------
C
 510  CONTINUE
C
      JAE1=3*KFLEV+1-JJ
      JAE2=3*KFLEV+1-(JJ+1)
      JAE3=3*KFLEV+1-JJPN
      DO 512 JAE=1,5
      DO 511 JL = 1, KDLON
      ZUAER(JL,JAE) = (RAER(JAE,1)*PAER(JL,JKL,1)
     S      +RAER(JAE,2)*PAER(JL,JKL,2)+RAER(JAE,3)*PAER(JL,JKL,3)
     S      +RAER(JAE,4)*PAER(JL,JKL,4)+RAER(JAE,5)*PAER(JL,JKL,5))
     S      /(ZDUC(JL,JAE1)+ZDUC(JL,JAE2)+ZDUC(JL,JAE3))
 511  CONTINUE
 512  CONTINUE
C
C
C
C*         5.2  INTRODUCES TEMPERATURE EFFECTS ON ABSORBER AMOUNTS
C               --------------------------------------------------
C
 520  CONTINUE
C
      DO 521 JL = 1, KDLON
      ZTAVI(JL)=PTAVE(JL,JKL)
      ZTCON(JL)=EXP(6.08*(296./ZTAVI(JL)-1.))
      ZTX=ZTAVI(JL)-TREF
      ZTX2=ZTX*ZTX
      ZZABLY = ZABLY(JL,6,JAE1)+ZABLY(JL,6,JAE2)+ZABLY(JL,6,JAE3)
CMAF      ZUP=MIN( MAX( 0.5*R10E*LOG( ZZABLY ) + 5., 0.), 6.0)
      ZUP=MIN( MAX( 0.5*R10E*LOG( ZZABLY ) + 5., 0.d+0), 6.d+0)
      ZCAH1=AT(1,1)+ZUP*(AT(1,2)+ZUP*(AT(1,3)))
      ZCBH1=BT(1,1)+ZUP*(BT(1,2)+ZUP*(BT(1,3)))
      ZPSH1(JL)=EXP( ZCAH1 * ZTX + ZCBH1 * ZTX2 )
      ZCAH2=AT(2,1)+ZUP*(AT(2,2)+ZUP*(AT(2,3)))
      ZCBH2=BT(2,1)+ZUP*(BT(2,2)+ZUP*(BT(2,3)))
      ZPSH2(JL)=EXP( ZCAH2 * ZTX + ZCBH2 * ZTX2 )
      ZCAH3=AT(3,1)+ZUP*(AT(3,2)+ZUP*(AT(3,3)))
      ZCBH3=BT(3,1)+ZUP*(BT(3,2)+ZUP*(BT(3,3)))
      ZPSH3(JL)=EXP( ZCAH3 * ZTX + ZCBH3 * ZTX2 )
      ZCAH4=AT(4,1)+ZUP*(AT(4,2)+ZUP*(AT(4,3)))
      ZCBH4=BT(4,1)+ZUP*(BT(4,2)+ZUP*(BT(4,3)))
      ZPSH4(JL)=EXP( ZCAH4 * ZTX + ZCBH4 * ZTX2 )
      ZCAH5=AT(5,1)+ZUP*(AT(5,2)+ZUP*(AT(5,3)))
      ZCBH5=BT(5,1)+ZUP*(BT(5,2)+ZUP*(BT(5,3)))
      ZPSH5(JL)=EXP( ZCAH5 * ZTX + ZCBH5 * ZTX2 )
      ZCAH6=AT(6,1)+ZUP*(AT(6,2)+ZUP*(AT(6,3)))
      ZCBH6=BT(6,1)+ZUP*(BT(6,2)+ZUP*(BT(6,3)))
      ZPSH6(JL)=EXP( ZCAH6 * ZTX + ZCBH6 * ZTX2 )
      ZPHM6(JL)=EXP(-5.81E-4 * ZTX - 1.13E-6 * ZTX2 )
      ZPSM6(JL)=EXP(-5.57E-4 * ZTX - 3.30E-6 * ZTX2 )
      ZPHN6(JL)=EXP(-3.46E-5 * ZTX + 2.05E-7 * ZTX2 )
      ZPSN6(JL)=EXP( 3.70E-3 * ZTX - 2.30E-6 * ZTX2 )
 521  CONTINUE
C
      DO 522 JL = 1, KDLON
      ZTAVI(JL)=PTAVE(JL,JKL)
      ZTX=ZTAVI(JL)-TREF
      ZTX2=ZTX*ZTX
      ZZABLY = ZABLY(JL,9,JAE1)+ZABLY(JL,9,JAE2)+ZABLY(JL,9,JAE3)
      ZALUP = R10E * LOG ( ZZABLY )
CMAF      ZUP   = MAX( 0.0 , 5.0 + 0.5 * ZALUP )
      ZUP   = MAX( 0.d+0 , 5.0 + 0.5 * ZALUP )
      ZPSC2(JL) = (ZTAVI(JL)/TREF) ** ZUP
      ZCAC8=AT(8,1)+ZUP*(AT(8,2)+ZUP*(AT(8,3)))
      ZCBC8=BT(8,1)+ZUP*(BT(8,2)+ZUP*(BT(8,3)))
      ZPSC3(JL)=EXP( ZCAC8 * ZTX + ZCBC8 * ZTX2 )
      ZPHIO(JL) = EXP( OCT(1) * ZTX + OCT(2) * ZTX2)
      ZPSIO(JL) = EXP( 2.* (OCT(3)*ZTX+OCT(4)*ZTX2))
 522  CONTINUE
C
      DO 524 JKK=JJ,JJPN
      JC=3*KFLEV+1-JKK
      JCP1=JC+1
      DO 523 JL = 1, KDLON
      ZDIFF = PVIEW(JL)
      PABCU(JL,10,JC)=PABCU(JL,10,JCP1)
     S                +ZABLY(JL,10,JC)           *ZDIFF
      PABCU(JL,11,JC)=PABCU(JL,11,JCP1)
     S                +ZABLY(JL,11,JC)*ZTCON(JL)*ZDIFF
C
      PABCU(JL,12,JC)=PABCU(JL,12,JCP1)
     S                +ZABLY(JL,12,JC)*ZPHIO(JL)*ZDIFF
      PABCU(JL,13,JC)=PABCU(JL,13,JCP1)
     S                +ZABLY(JL,13,JC)*ZPSIO(JL)*ZDIFF
C
      PABCU(JL,7,JC)=PABCU(JL,7,JCP1)
     S               +ZABLY(JL,9,JC)*ZPSC2(JL)*ZDIFF
      PABCU(JL,8,JC)=PABCU(JL,8,JCP1)
     S               +ZABLY(JL,9,JC)*ZPSC3(JL)*ZDIFF
      PABCU(JL,9,JC)=PABCU(JL,9,JCP1)
     S               +ZABLY(JL,9,JC)*ZPSC3(JL)*ZDIFF
C
      PABCU(JL,1,JC)=PABCU(JL,1,JCP1)
     S               +ZABLY(JL,6,JC)*ZPSH1(JL)*ZDIFF
      PABCU(JL,2,JC)=PABCU(JL,2,JCP1)
     S               +ZABLY(JL,6,JC)*ZPSH2(JL)*ZDIFF
      PABCU(JL,3,JC)=PABCU(JL,3,JCP1)
     S               +ZABLY(JL,6,JC)*ZPSH5(JL)*ZDIFF
      PABCU(JL,4,JC)=PABCU(JL,4,JCP1)
     S               +ZABLY(JL,6,JC)*ZPSH3(JL)*ZDIFF
      PABCU(JL,5,JC)=PABCU(JL,5,JCP1)
     S               +ZABLY(JL,6,JC)*ZPSH4(JL)*ZDIFF
      PABCU(JL,6,JC)=PABCU(JL,6,JCP1)
     S               +ZABLY(JL,6,JC)*ZPSH6(JL)*ZDIFF
C
      PABCU(JL,14,JC)=PABCU(JL,14,JCP1)
     S                +ZUAER(JL,1)    *ZDUC(JL,JC)*ZDIFF
      PABCU(JL,15,JC)=PABCU(JL,15,JCP1)
     S                +ZUAER(JL,2)    *ZDUC(JL,JC)*ZDIFF
      PABCU(JL,16,JC)=PABCU(JL,16,JCP1)
     S                +ZUAER(JL,3)    *ZDUC(JL,JC)*ZDIFF
      PABCU(JL,17,JC)=PABCU(JL,17,JCP1)
     S                +ZUAER(JL,4)    *ZDUC(JL,JC)*ZDIFF
      PABCU(JL,18,JC)=PABCU(JL,18,JCP1)
     S                +ZUAER(JL,5)    *ZDUC(JL,JC)*ZDIFF
C
      PABCU(JL,19,JC)=PABCU(JL,19,JCP1)
     S               +ZABLY(JL,8,JC)*RCH4/RCO2*ZPHM6(JL)*ZDIFF
      PABCU(JL,20,JC)=PABCU(JL,20,JCP1)
     S               +ZABLY(JL,9,JC)*RCH4/RCO2*ZPSM6(JL)*ZDIFF
      PABCU(JL,21,JC)=PABCU(JL,21,JCP1)
     S               +ZABLY(JL,8,JC)*RN2O/RCO2*ZPHN6(JL)*ZDIFF
      PABCU(JL,22,JC)=PABCU(JL,22,JCP1)
     S               +ZABLY(JL,9,JC)*RN2O/RCO2*ZPSN6(JL)*ZDIFF
C
      PABCU(JL,23,JC)=PABCU(JL,23,JCP1)
     S               +ZABLY(JL,8,JC)*RCFC11/RCO2         *ZDIFF
      PABCU(JL,24,JC)=PABCU(JL,24,JCP1)
     S               +ZABLY(JL,8,JC)*RCFC12/RCO2         *ZDIFF
 523  CONTINUE
 524  CONTINUE
C
 529  CONTINUE
C
C
      RETURN
      END
