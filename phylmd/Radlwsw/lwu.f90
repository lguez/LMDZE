module LWU_m

  IMPLICIT none

contains

  SUBROUTINE LWU(PDP, PPMB, POZ, PTAVE, PVIEW, PWV, PABCU)

    ! Purpose. Computes absorber amounts including pressure and
    ! temperature effects.

    ! Method. Computes the pressure and temperature weighted amounts
    ! of absorbers.

    ! Reference. See radiation's part of the model's documentation and
    ! ECMWF research department documentation of the IFS.

    ! Author. Jean-Jacques Morcrette, ECMWF.

    ! Modifications.
    ! Original : 89-07-14
    ! Voigt lines (loop 404 modified) - JJM & PhD - 01/96

    USE clesphys, ONLY: rcfc11, rcfc12, rch4, rco2, rn2o
    use dimphy, only: klon
    USE suphec_m, ONLY: rg
    use dimensions, only: llm
    USE radepsi, ONLY: zepsco, zepscq
    USE radopt, ONLY: levoigt
    USE raddimlw, ONLY: ng1, ng1p1, ninter, nua

    ! ARGUMENTS:

    DOUBLE PRECISION PDP(KLON, LLM)
    DOUBLE PRECISION PPMB(KLON, LLM + 1)
    DOUBLE PRECISION POZ(KLON, LLM)
    DOUBLE PRECISION PTAVE(KLON, LLM)
    DOUBLE PRECISION PVIEW(KLON)
    DOUBLE PRECISION PWV(KLON, LLM)

    DOUBLE PRECISION PABCU(KLON, NUA, 3 * LLM + 1)
    ! effective absorber amounts

    ! LOCAL VARIABLES:

    DOUBLE PRECISION PAER ! OPTICAL THICKNESS OF THE AEROSOLS
    DOUBLE PRECISION ZABLY(KLON, NUA, 3 * LLM + 1)
    DOUBLE PRECISION ZDUC(KLON, 3 * LLM + 1)
    DOUBLE PRECISION ZPHIO(KLON)
    DOUBLE PRECISION ZPSC2(KLON)
    DOUBLE PRECISION ZPSC3(KLON)
    DOUBLE PRECISION ZPSH1(KLON)
    DOUBLE PRECISION ZPSH2(KLON)
    DOUBLE PRECISION ZPSH3(KLON)
    DOUBLE PRECISION ZPSH4(KLON)
    DOUBLE PRECISION ZPSH5(KLON)
    DOUBLE PRECISION ZPSH6(KLON)
    DOUBLE PRECISION ZPSIO(KLON)
    DOUBLE PRECISION ZTCON(KLON)
    DOUBLE PRECISION ZPHM6(KLON)
    DOUBLE PRECISION ZPSM6(KLON)
    DOUBLE PRECISION ZPHN6(KLON)
    DOUBLE PRECISION ZPSN6(KLON)
    DOUBLE PRECISION ZSSIG(KLON, 3 * LLM + 1)
    DOUBLE PRECISION ZTAVI(KLON)
    DOUBLE PRECISION ZUAER(KLON, Ninter)
    DOUBLE PRECISION ZXOZ(KLON)
    DOUBLE PRECISION ZXWV(KLON)

    INTEGER jl, jk, jkj, jkjr, jkjp, ig1
    INTEGER jki, jkip1, ja, jj
    INTEGER jkl, jkk, jkjpn
    INTEGER jae1, jae2, jae3, jae, jjpn
    INTEGER ir, jc, jcp1
    DOUBLE PRECISION zdpm, zupm, zupmh2o, zupmco2, zupmo3, zu6, zup
    DOUBLE PRECISION zfppw, ztx, ztx2, zzably
    DOUBLE PRECISION zcah1, zcbh1, zcah2, zcbh2, zcah3, zcbh3
    DOUBLE PRECISION zcah4, zcbh4, zcah5, zcbh5, zcah6, zcbh6
    DOUBLE PRECISION zcac8, zcbc8
    DOUBLE PRECISION zalup, zdiff

    DOUBLE PRECISION PVGCO2, PVGH2O, PVGO3

    DOUBLE PRECISION, PARAMETER:: R10E = 0.4342945
    ! decimal / natural logarithm factor

    ! Used Data Block:

    DOUBLE PRECISION:: TREF = 250d0
    DOUBLE PRECISION:: RT1(2) = (/ - 0.577350269d0, 0.577350269d0/)
    DOUBLE PRECISION RAER(5, 5)
    DOUBLE PRECISION AT(8, 3), BT(8, 3)
    DOUBLE PRECISION:: OCT(4) = (/- 0.326D-3, - 0.102D-5, 0.137D-2, - 0.535D-5/)

    DATA RAER / .038520d0, .037196d0, .040532d0, .054934d0, .038520d0, &
         .12613d0, .18313d0, .10357d0, .064106d0, .126130d0, &
         .012579d0, .013649d0, .018652d0, .025181d0, .012579d0, &
         .011890d0, .016142d0, .021105d0, .028908d0, .011890d0, &
         .013792d0, .026810d0, .052203d0, .066338d0, .013792d0 /

    DATA (AT(1, IR), IR = 1, 3) / 0.298199D-02, - .394023D-03, 0.319566D-04 /
    DATA (BT(1, IR), IR = 1, 3) / - 0.106432D-04, 0.660324D-06, 0.174356D-06 /
    DATA (AT(2, IR), IR = 1, 3) / 0.143676D-01, 0.366501D-02, -.160822D-02 /
    DATA (BT(2, IR), IR = 1, 3) / -0.553979D-04, - .101701D-04, 0.920868D-05 /
    DATA (AT(3, IR), IR = 1, 3) / 0.197861D-01, 0.315541D-02, - .174547D-02 /
    DATA (BT(3, IR), IR = 1, 3) / - 0.877012D-04, 0.513302D-04, 0.523138D-06 /
    DATA (AT(4, IR), IR = 1, 3) / 0.289560D-01, - .208807D-02, - .121943D-02 /
    DATA (BT(4, IR), IR = 1, 3) / - 0.165960D-03, 0.157704D-03, - .146427D-04 /
    DATA (AT(5, IR), IR = 1, 3) / 0.103800D-01, 0.436296D-02, - .161431D-02 /
    DATA (BT(5, IR), IR = 1, 3) / - .276744D-04, - .327381D-04, 0.127646D-04 /
    DATA (AT(6, IR), IR = 1, 3) / 0.868859D-02, - .972752D-03, 0.000000D-00 /
    DATA (BT(6, IR), IR = 1, 3) / - .278412D-04, - .713940D-06, 0.117469D-05 /
    DATA (AT(7, IR), IR = 1, 3) / 0.250073D-03, 0.455875D-03, 0.109242D-03 /
    DATA (BT(7, IR), IR = 1, 3) / 0.199846D-05, - .216313D-05, 0.175991D-06 /
    DATA (AT(8, IR), IR = 1, 3) / 0.307423D-01, 0.110879D-02, - .322172D-03 /
    DATA (BT(8, IR), IR = 1, 3) / - 0.108482D-03, 0.258096D-05, - .814575D-06 /

    !-----------------------------------------------------------------------

    paer = 1d-15

    IF (LEVOIGT) THEN
       PVGCO2 = 60.
       PVGH2O = 30.
       PVGO3 = 400.
    ELSE
       PVGCO2 = 0.
       PVGH2O = 0.
       PVGO3 = 0.
    ENDIF

    ! 2. PRESSURE OVER GAUSS SUB-LEVELS

    DO JL = 1, KLON
       ZSSIG(JL, 1) = PPMB(JL, 1) * 100.
    end DO

    DO JK = 1, LLM
       JKJ = (JK - 1) * NG1P1 + 1
       JKJR = JKJ
       JKJP = JKJ + NG1P1
       DO JL = 1, KLON
          ZSSIG(JL, JKJP) = PPMB(JL, JK + 1) * 100.
       end DO
       DO IG1 = 1, NG1
          JKJ = JKJ + 1
          DO JL = 1, KLON
             ZSSIG(JL, JKJ) = (ZSSIG(JL, JKJR) + ZSSIG(JL, JKJP)) * 0.5 &
                  + RT1(IG1) * (ZSSIG(JL, JKJP) - ZSSIG(JL, JKJR)) * 0.5
          end DO
       end DO
    end DO

    ! 4. PRESSURE THICKNESS AND MEAN PRESSURE OF SUB-LAYERS

    DO JKI = 1, 3 * LLM
       JKIP1 = JKI + 1
       DO JL = 1, KLON
          ZABLY(JL, 5, JKI) = (ZSSIG(JL, JKI) + ZSSIG(JL, JKIP1)) * 0.5
          ZABLY(JL, 3, JKI) = (ZSSIG(JL, JKI) - ZSSIG(JL, JKIP1)) / (10. * RG)
       end DO
    end DO

    DO JK = 1, LLM
       JKL = LLM + 1 - JK
       DO JL = 1, KLON
          ZXWV(JL) = MAX(PWV(JL, JK), ZEPSCQ)
          ZXOZ(JL) = MAX(POZ(JL, JK) / PDP(JL, JK), ZEPSCO)
       end DO
       JKJ = (JK - 1) * NG1P1 + 1
       JKJPN = JKJ + NG1
       DO JKK = JKJ, JKJPN
          DO JL = 1, KLON
             ZDPM = ZABLY(JL, 3, JKK)
             ZUPM = ZABLY(JL, 5, JKK) * ZDPM / 101325.
             ZUPMCO2 = (ZABLY(JL, 5, JKK) + PVGCO2) * ZDPM / 101325.
             ZUPMH2O = (ZABLY(JL, 5, JKK) + PVGH2O) * ZDPM / 101325.
             ZUPMO3 = (ZABLY(JL, 5, JKK) + PVGO3) * ZDPM / 101325.
             ZDUC(JL, JKK) = ZDPM
             ZABLY(JL, 12, JKK) = ZXOZ(JL) * ZDPM
             ZABLY(JL, 13, JKK) = ZXOZ(JL) * ZUPMO3
             ZU6 = ZXWV(JL) * ZUPM
             ZFPPW = 1.6078 * ZXWV(JL) / (1. + 0.608 * ZXWV(JL))
             ZABLY(JL, 6, JKK) = ZXWV(JL) * ZUPMH2O
             ZABLY(JL, 11, JKK) = ZU6 * ZFPPW
             ZABLY(JL, 10, JKK) = ZU6 * (1. - ZFPPW)
             ZABLY(JL, 9, JKK) = RCO2 * ZUPMCO2
             ZABLY(JL, 8, JKK) = RCO2 * ZDPM
          end DO
       end DO
    end DO

    ! 5. CUMULATIVE ABSORBER AMOUNTS FROM TOP OF ATMOSPHERE

    DO JA = 1, NUA
       DO JL = 1, KLON
          PABCU(JL, JA, 3 * LLM + 1) = 0.
       end DO
    end DO

    DO JK = 1, LLM
       JJ = (JK - 1) * NG1P1 + 1
       JJPN = JJ + NG1
       JKL = LLM + 1 - JK

       ! 5.1 CUMULATIVE AEROSOL AMOUNTS FROM TOP OF ATMOSPHERE

       JAE1 = 3 * LLM + 1 - JJ
       JAE2 = 3 * LLM + 1 - (JJ + 1)
       JAE3 = 3 * LLM + 1 - JJPN
       DO JAE = 1, 5
          DO JL = 1, KLON
             ZUAER(JL, JAE) = (RAER(JAE, 1) * PAER &
                  + RAER(JAE, 2) * PAER &
                  + RAER(JAE, 3) * PAER &
                  + RAER(JAE, 4) * PAER &
                  + RAER(JAE, 5) * PAER) &
                  / (ZDUC(JL, JAE1) + ZDUC(JL, JAE2) + ZDUC(JL, JAE3))
          end DO
       end DO

       ! 5.2 INTRODUCES TEMPERATURE EFFECTS ON ABSORBER AMOUNTS

       DO JL = 1, KLON
          ZTAVI(JL) = PTAVE(JL, JKL)
          ZTCON(JL) = EXP(6.08 * (296. / ZTAVI(JL) - 1.))
          ZTX = ZTAVI(JL) - TREF
          ZTX2 = ZTX * ZTX
          ZZABLY = ZABLY(JL, 6, JAE1) + ZABLY(JL, 6, JAE2) + ZABLY(JL, 6, JAE3)
          ZUP = MIN(MAX(0.5 * R10E * LOG(ZZABLY) + 5., 0d0), 6d0)
          ZCAH1 = AT(1, 1) + ZUP * (AT(1, 2) + ZUP * (AT(1, 3)))
          ZCBH1 = BT(1, 1) + ZUP * (BT(1, 2) + ZUP * (BT(1, 3)))
          ZPSH1(JL) = EXP(ZCAH1 * ZTX + ZCBH1 * ZTX2)
          ZCAH2 = AT(2, 1) + ZUP * (AT(2, 2) + ZUP * (AT(2, 3)))
          ZCBH2 = BT(2, 1) + ZUP * (BT(2, 2) + ZUP * (BT(2, 3)))
          ZPSH2(JL) = EXP(ZCAH2 * ZTX + ZCBH2 * ZTX2)
          ZCAH3 = AT(3, 1) + ZUP * (AT(3, 2) + ZUP * (AT(3, 3)))
          ZCBH3 = BT(3, 1) + ZUP * (BT(3, 2) + ZUP * (BT(3, 3)))
          ZPSH3(JL) = EXP(ZCAH3 * ZTX + ZCBH3 * ZTX2)
          ZCAH4 = AT(4, 1) + ZUP * (AT(4, 2) + ZUP * (AT(4, 3)))
          ZCBH4 = BT(4, 1) + ZUP * (BT(4, 2) + ZUP * (BT(4, 3)))
          ZPSH4(JL) = EXP(ZCAH4 * ZTX + ZCBH4 * ZTX2)
          ZCAH5 = AT(5, 1) + ZUP * (AT(5, 2) + ZUP * (AT(5, 3)))
          ZCBH5 = BT(5, 1) + ZUP * (BT(5, 2) + ZUP * (BT(5, 3)))
          ZPSH5(JL) = EXP(ZCAH5 * ZTX + ZCBH5 * ZTX2)
          ZCAH6 = AT(6, 1) + ZUP * (AT(6, 2) + ZUP * (AT(6, 3)))
          ZCBH6 = BT(6, 1) + ZUP * (BT(6, 2) + ZUP * (BT(6, 3)))
          ZPSH6(JL) = EXP(ZCAH6 * ZTX + ZCBH6 * ZTX2)
          ZPHM6(JL) = EXP(- 5.81E-4 * ZTX - 1.13E-6 * ZTX2)
          ZPSM6(JL) = EXP(- 5.57E-4 * ZTX - 3.30E-6 * ZTX2)
          ZPHN6(JL) = EXP(- 3.46E-5 * ZTX + 2.05E-7 * ZTX2)
          ZPSN6(JL) = EXP(3.70E-3 * ZTX - 2.30E-6 * ZTX2)
       end DO

       DO JL = 1, KLON
          ZTAVI(JL) = PTAVE(JL, JKL)
          ZTX = ZTAVI(JL) - TREF
          ZTX2 = ZTX * ZTX
          ZZABLY = ZABLY(JL, 9, JAE1) + ZABLY(JL, 9, JAE2) + ZABLY(JL, 9, JAE3)
          ZALUP = R10E * LOG(ZZABLY)
          ZUP = MAX(0d0, 5.0 + 0.5 * ZALUP)
          ZPSC2(JL) = (ZTAVI(JL) / TREF) ** ZUP
          ZCAC8 = AT(8, 1) + ZUP * (AT(8, 2) + ZUP * (AT(8, 3)))
          ZCBC8 = BT(8, 1) + ZUP * (BT(8, 2) + ZUP * (BT(8, 3)))
          ZPSC3(JL) = EXP(ZCAC8 * ZTX + ZCBC8 * ZTX2)
          ZPHIO(JL) = EXP(OCT(1) * ZTX + OCT(2) * ZTX2)
          ZPSIO(JL) = EXP(2. * (OCT(3) * ZTX + OCT(4) * ZTX2))
       end DO

       DO JKK = JJ, JJPN
          JC = 3 * LLM + 1 - JKK
          JCP1 = JC + 1
          DO JL = 1, KLON
             ZDIFF = PVIEW(JL)
             PABCU(JL, 10, JC) = PABCU(JL, 10, JCP1) &
                  + ZABLY(JL, 10, JC) * ZDIFF
             PABCU(JL, 11, JC) = PABCU(JL, 11, JCP1) &
                  + ZABLY(JL, 11, JC) * ZTCON(JL) * ZDIFF

             PABCU(JL, 12, JC) = PABCU(JL, 12, JCP1) &
                  + ZABLY(JL, 12, JC) * ZPHIO(JL) * ZDIFF
             PABCU(JL, 13, JC) = PABCU(JL, 13, JCP1) &
                  + ZABLY(JL, 13, JC) * ZPSIO(JL) * ZDIFF

             PABCU(JL, 7, JC) = PABCU(JL, 7, JCP1) &
                  + ZABLY(JL, 9, JC) * ZPSC2(JL) * ZDIFF
             PABCU(JL, 8, JC) = PABCU(JL, 8, JCP1) &
                  + ZABLY(JL, 9, JC) * ZPSC3(JL) * ZDIFF
             PABCU(JL, 9, JC) = PABCU(JL, 9, JCP1) &
                  + ZABLY(JL, 9, JC) * ZPSC3(JL) * ZDIFF

             PABCU(JL, 1, JC) = PABCU(JL, 1, JCP1) &
                  + ZABLY(JL, 6, JC) * ZPSH1(JL) * ZDIFF
             PABCU(JL, 2, JC) = PABCU(JL, 2, JCP1) &
                  + ZABLY(JL, 6, JC) * ZPSH2(JL) * ZDIFF
             PABCU(JL, 3, JC) = PABCU(JL, 3, JCP1) &
                  + ZABLY(JL, 6, JC) * ZPSH5(JL) * ZDIFF
             PABCU(JL, 4, JC) = PABCU(JL, 4, JCP1) &
                  + ZABLY(JL, 6, JC) * ZPSH3(JL) * ZDIFF
             PABCU(JL, 5, JC) = PABCU(JL, 5, JCP1) &
                  + ZABLY(JL, 6, JC) * ZPSH4(JL) * ZDIFF
             PABCU(JL, 6, JC) = PABCU(JL, 6, JCP1) &
                  + ZABLY(JL, 6, JC) * ZPSH6(JL) * ZDIFF

             PABCU(JL, 14, JC) = PABCU(JL, 14, JCP1) &
                  + ZUAER(JL, 1) * ZDUC(JL, JC) * ZDIFF
             PABCU(JL, 15, JC) = PABCU(JL, 15, JCP1) &
                  + ZUAER(JL, 2) * ZDUC(JL, JC) * ZDIFF
             PABCU(JL, 16, JC) = PABCU(JL, 16, JCP1) &
                  + ZUAER(JL, 3) * ZDUC(JL, JC) * ZDIFF
             PABCU(JL, 17, JC) = PABCU(JL, 17, JCP1) &
                  + ZUAER(JL, 4) * ZDUC(JL, JC) * ZDIFF
             PABCU(JL, 18, JC) = PABCU(JL, 18, JCP1) &
                  + ZUAER(JL, 5) * ZDUC(JL, JC) * ZDIFF

             PABCU(JL, 19, JC) = PABCU(JL, 19, JCP1) &
                  + ZABLY(JL, 8, JC) * RCH4 / RCO2 * ZPHM6(JL) * ZDIFF
             PABCU(JL, 20, JC) = PABCU(JL, 20, JCP1) &
                  + ZABLY(JL, 9, JC) * RCH4 / RCO2 * ZPSM6(JL) * ZDIFF
             PABCU(JL, 21, JC) = PABCU(JL, 21, JCP1) &
                  + ZABLY(JL, 8, JC) * RN2O / RCO2 * ZPHN6(JL) * ZDIFF
             PABCU(JL, 22, JC) = PABCU(JL, 22, JCP1) &
                  + ZABLY(JL, 9, JC) * RN2O / RCO2 * ZPSN6(JL) * ZDIFF

             PABCU(JL, 23, JC) = PABCU(JL, 23, JCP1) &
                  + ZABLY(JL, 8, JC) * RCFC11 / RCO2 * ZDIFF
             PABCU(JL, 24, JC) = PABCU(JL, 24, JCP1) &
                  + ZABLY(JL, 8, JC) * RCFC12 / RCO2 * ZDIFF
          end DO
       end DO
    end DO

  END SUBROUTINE LWU

end module LWU_m
