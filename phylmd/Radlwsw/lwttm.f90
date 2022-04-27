SUBROUTINE lwttm(pga, pgb, puu1, puu2, ptt)
  USE dimensions
  USE dimphy
  USE raddimlw
  IMPLICIT NONE

  ! ------------------------------------------------------------------
  ! PURPOSE.
  ! --------
  ! THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
  ! ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN ALL SIX SPECTRAL
  ! INTERVALS.

  ! METHOD.
  ! -------

  ! 1. TRANSMISSION FUNCTION BY H2O AND UNIFORMLY MIXED GASES ARE
  ! COMPUTED USING PADE APPROXIMANTS AND HORNER'S ALGORITHM.
  ! 2. TRANSMISSION BY O3 IS EVALUATED WITH MALKMUS'S BAND MODEL.
  ! 3. TRANSMISSION BY H2O CONTINUUM AND AEROSOLS FOLLOW AN
  ! A SIMPLE EXPONENTIAL DECREASE WITH ABSORBER AMOUNT.

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
  ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 88-12-15

  ! -----------------------------------------------------------------------
  DOUBLE PRECISION o1h, o2h
  PARAMETER (o1h=2230.)
  PARAMETER (o2h=100.)
  DOUBLE PRECISION rpialf0
  PARAMETER (rpialf0=2.0)

  ! * ARGUMENTS:

  DOUBLE PRECISION pga(klon, 8, 2) ! PADE APPROXIMANTS
  DOUBLE PRECISION pgb(klon, 8, 2) ! PADE APPROXIMANTS
  DOUBLE PRECISION puu1(klon, nua) ! ABSORBER AMOUNTS FROM TOP TO LEVEL 1
  DOUBLE PRECISION puu2(klon, nua) ! ABSORBER AMOUNTS FROM TOP TO LEVEL 2
  DOUBLE PRECISION ptt(klon, ntra) ! TRANSMISSION FUNCTIONS

  ! * LOCAL VARIABLES:

  INTEGER ja, jl
  DOUBLE PRECISION zz, zxd, zxn
  DOUBLE PRECISION zpu, zpu10, zpu11, zpu12, zpu13
  DOUBLE PRECISION zeu, zeu10, zeu11, zeu12, zeu13
  DOUBLE PRECISION zx, zy, zuxy, zsq1, zsq2, zvxy, zaercn, zto1
  DOUBLE PRECISION zto2
  DOUBLE PRECISION zxch4, zych4, zsqh41, zodh41
  DOUBLE PRECISION zxn2o, zyn2o, zsqn21, zodn21, zsqh42, zodh42
  DOUBLE PRECISION zsqn22, zodn22, za11, zttf11, za12, zttf12
  DOUBLE PRECISION zuu11, zuu12
  ! ------------------------------------------------------------------

  ! *         1.     HORNER'S ALGORITHM FOR H2O AND CO2 TRANSMISSION
  ! -----------------------------------------------



  DO ja = 1, 8
    DO jl = 1, klon
      zz = sqrt(puu1(jl,ja)-puu2(jl,ja))
      zxd = pgb(jl, ja, 1) + zz*(pgb(jl,ja,2)+zz)
      zxn = pga(jl, ja, 1) + zz*(pga(jl,ja,2))
      ptt(jl, ja) = zxn/zxd
    END DO
  END DO

  ! ------------------------------------------------------------------

  ! *         2.     CONTINUUM, OZONE AND AEROSOL TRANSMISSION FUNCTIONS
  ! ---------------------------------------------------


  DO jl = 1, klon
    ptt(jl, 9) = ptt(jl, 8)

    ! -  CONTINUUM ABSORPTION: E- AND P-TYPE

    zpu = 0.002*(puu1(jl,10)-puu2(jl,10))
    zpu10 = 112.*zpu
    zpu11 = 6.25*zpu
    zpu12 = 5.00*zpu
    zpu13 = 80.0*zpu
    zeu = (puu1(jl,11)-puu2(jl,11))
    zeu10 = 12.*zeu
    zeu11 = 6.25*zeu
    zeu12 = 5.00*zeu
    zeu13 = 80.0*zeu

    ! -  OZONE ABSORPTION

    zx = (puu1(jl,12)-puu2(jl,12))
    zy = (puu1(jl,13)-puu2(jl,13))
    zuxy = 4.*zx*zx/(rpialf0*zy)
    zsq1 = sqrt(1.+o1h*zuxy) - 1.
    zsq2 = sqrt(1.+o2h*zuxy) - 1.
    zvxy = rpialf0*zy/(2.*zx)
    zaercn = (puu1(jl,17)-puu2(jl,17)) + zeu12 + zpu12
    zto1 = exp(-zvxy*zsq1-zaercn)
    zto2 = exp(-zvxy*zsq2-zaercn)

    ! -- TRACE GASES (CH4, N2O, CFC-11, CFC-12)

    ! * CH4 IN INTERVAL 800-970 + 1110-1250 CM-1

    zxch4 = (puu1(jl,19)-puu2(jl,19))
    zych4 = (puu1(jl,20)-puu2(jl,20))
    zuxy = 4.*zxch4*zxch4/(0.103*zych4)
    zsqh41 = sqrt(1.+33.7*zuxy) - 1.
    zvxy = 0.103*zych4/(2.*zxch4)
    zodh41 = zvxy*zsqh41

    ! * N2O IN INTERVAL 800-970 + 1110-1250 CM-1

    zxn2o = (puu1(jl,21)-puu2(jl,21))
    zyn2o = (puu1(jl,22)-puu2(jl,22))
    zuxy = 4.*zxn2o*zxn2o/(0.416*zyn2o)
    zsqn21 = sqrt(1.+21.3*zuxy) - 1.
    zvxy = 0.416*zyn2o/(2.*zxn2o)
    zodn21 = zvxy*zsqn21

    ! * CH4 IN INTERVAL 1250-1450 + 1880-2820 CM-1

    zuxy = 4.*zxch4*zxch4/(0.113*zych4)
    zsqh42 = sqrt(1.+400.*zuxy) - 1.
    zvxy = 0.113*zych4/(2.*zxch4)
    zodh42 = zvxy*zsqh42

    ! * N2O IN INTERVAL 1250-1450 + 1880-2820 CM-1

    zuxy = 4.*zxn2o*zxn2o/(0.197*zyn2o)
    zsqn22 = sqrt(1.+2000.*zuxy) - 1.
    zvxy = 0.197*zyn2o/(2.*zxn2o)
    zodn22 = zvxy*zsqn22

    ! * CFC-11 IN INTERVAL 800-970 + 1110-1250 CM-1

    za11 = (puu1(jl,23)-puu2(jl,23))*4.404E+05
    zttf11 = 1. - za11*0.003225

    ! * CFC-12 IN INTERVAL 800-970 + 1110-1250 CM-1

    za12 = (puu1(jl,24)-puu2(jl,24))*6.7435E+05
    zttf12 = 1. - za12*0.003225

    zuu11 = -(puu1(jl,15)-puu2(jl,15)) - zeu10 - zpu10
    zuu12 = -(puu1(jl,16)-puu2(jl,16)) - zeu11 - zpu11 - zodh41 - zodn21
    ptt(jl, 10) = exp(-(puu1(jl,14)-puu2(jl,14)))
    ptt(jl, 11) = exp(zuu11)
    ptt(jl, 12) = exp(zuu12)*zttf11*zttf12
    ptt(jl, 13) = 0.7554*zto1 + 0.2446*zto2
    ptt(jl, 14) = ptt(jl, 10)*exp(-zeu13-zpu13)
    ptt(jl, 15) = exp(-(puu1(jl,14)-puu2(jl,14))-zodh42-zodn22)
  END DO

  RETURN
END SUBROUTINE lwttm
