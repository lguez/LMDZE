module sw_m

  IMPLICIT none

contains

  SUBROUTINE SW(PSCT, PRMU0, FRACT, PPMB, PDP, PPSOL, PALBD, PALBP, PTAVE, &
       PWV, PQS, POZON, PCLDSW, PTAU, POMEGA, PCG, PHEAT, PHEAT0, PTOPSW, &
       PSOLSW, PTOPSW0, PSOLSW0, ZFSUP, ZFSDN, ZFSUP0, ZFSDN0, PTOPSWAD, &
       PSOLSWAD)

    ! Purpose.
    ! This routine computes the shortwave radiation fluxes in two
    ! spectral intervals following Fouquart and Bonnel (1980).

    ! Method.
    ! 1. Computes absorber amounts (swu)
    ! 2. Computes fluxes in 1st spectral interval (SW1S)
    ! 3. Computes fluxes in 2nd spectral interval (SW2S)

    ! Reference.
    ! See radiation part of the ECMWF research department
    ! documentation, and Fouquart and Bonnel (1980)

    ! Author.
    ! Jean-Jacques Morcrette *ecmwf*

    ! Modifications.
    ! Original: 89-07-14
    ! 95-01-01 J.-J. Morcrette direct/diffuse albedo
    ! 03-11-27 J. Quaas Introduce aerosol forcings (based on Boucher)

    use comconst, only: daysec
    use dimphy, only: klon
    use dimensions, only: llm
    USE suphec_m, ONLY: rcpd, rg
    use sw1s_m, only: sw1s
    use sw2s_m, only: sw2s
    use swu_m, only: swu

    ! ARGUMENTS:

    DOUBLE PRECISION PSCT ! constante solaire
    DOUBLE PRECISION PRMU0(KLON) ! COSINE OF ZENITHAL ANGLE
    DOUBLE PRECISION, intent(in):: FRACT(KLON) ! fraction de la journee
    DOUBLE PRECISION PPMB(KLON, LLM+1) ! HALF-LEVEL PRESSURE (MB)
    DOUBLE PRECISION PDP(KLON, LLM) ! LAYER THICKNESS (PA)
    DOUBLE PRECISION PPSOL(KLON) ! SURFACE PRESSURE (PA)
    DOUBLE PRECISION PALBD(KLON, 2) ! albedo du sol (lumiere diffuse)
    DOUBLE PRECISION PALBP(KLON, 2) ! albedo du sol (lumiere parallele)
    DOUBLE PRECISION PTAVE(KLON, LLM) ! LAYER TEMPERATURE (K)
    DOUBLE PRECISION PWV(KLON, LLM) ! SPECIFIC HUMIDITY (KG/KG)
    DOUBLE PRECISION PQS(KLON, LLM) ! SATURATED WATER VAPOUR (KG/KG)
    DOUBLE PRECISION POZON(KLON, LLM) ! OZONE CONCENTRATION (KG/KG)
    DOUBLE PRECISION PCLDSW(KLON, LLM) ! CLOUD FRACTION
    DOUBLE PRECISION PTAU(KLON, 2, LLM) ! CLOUD OPTICAL THICKNESS
    DOUBLE PRECISION POMEGA(KLON, 2, LLM) ! SINGLE SCATTERING ALBEDO
    DOUBLE PRECISION PCG(KLON, 2, LLM) ! ASYMETRY FACTOR
    DOUBLE PRECISION PHEAT(KLON, LLM) ! SHORTWAVE HEATING (K/DAY)
    DOUBLE PRECISION PHEAT0(KLON, LLM)! SHORTWAVE HEATING (K/DAY) clear-sky
    DOUBLE PRECISION PTOPSW(KLON) ! SHORTWAVE FLUX AT T.O.A.
    DOUBLE PRECISION PSOLSW(KLON) ! SHORTWAVE FLUX AT SURFACE
    DOUBLE PRECISION PTOPSW0(KLON) ! SHORTWAVE FLUX AT T.O.A. (CLEAR-SKY)
    DOUBLE PRECISION PSOLSW0(KLON) ! SHORTWAVE FLUX AT SURFACE (CLEAR-SKY)
    DOUBLE PRECISION ZFSUP(KLON, LLM+1)
    DOUBLE PRECISION ZFSDN(KLON, LLM+1)
    DOUBLE PRECISION ZFSUP0(KLON, LLM+1)
    DOUBLE PRECISION ZFSDN0(KLON, LLM+1)

    DOUBLE PRECISION, intent(out):: PTOPSWAD(KLON)
    ! (diagnosed aerosol forcing)SHORTWAVE FLUX AT T.O.A.(+AEROSOL DIR)

    DOUBLE PRECISION, intent(out):: PSOLSWAD(KLON)
    ! (diagnosed aerosol forcing)SHORTWAVE FLUX AT SURFACE(+AEROSOL DIR)

    ! Local:

    DOUBLE PRECISION ZOZ(KLON, LLM)
    DOUBLE PRECISION ZAKI(KLON, 2)
    DOUBLE PRECISION ZCLD(KLON, LLM)
    DOUBLE PRECISION ZCLEAR(KLON)
    DOUBLE PRECISION ZDSIG(KLON, LLM)
    DOUBLE PRECISION ZFACT(KLON)
    DOUBLE PRECISION ZFD(KLON, LLM+1)
    DOUBLE PRECISION ZFDOWN(KLON, LLM+1)
    DOUBLE PRECISION ZFU(KLON, LLM+1)
    DOUBLE PRECISION ZFUP(KLON, LLM+1)
    DOUBLE PRECISION ZRMU(KLON)
    DOUBLE PRECISION ZSEC(KLON)
    DOUBLE PRECISION ZUD(KLON, 5, LLM+1)
    DOUBLE PRECISION ZCLDSW0(KLON, LLM)
    INTEGER inu, jl, jk, i, k, kpl1
    INTEGER, PARAMETER:: swpas = 1 ! Every swpas steps, sw is calculated

    INTEGER:: itapsw = 0
    LOGICAL:: appel1er = .TRUE.
    !jq-Introduced for aerosol forcings

    !jq - Fluxes including aerosol effects
    DOUBLE PRECISION, save, allocatable:: ZFSUPAD(:, :) ! (KLON, LLM+1)
    DOUBLE PRECISION, save, allocatable:: ZFSDNAD(:, :) ! (KLON, LLM+1)

    logical:: initialized = .false.
    REAL, PARAMETER :: dobson_u = 2.1415E-05 ! Dobson unit, in kg m-2

    !-------------------------------------------------------------------

    if(.not.initialized) then
       initialized=.TRUE.
       allocate(ZFSUPAD(KLON, LLM+1), ZFSDNAD(KLON, LLM+1))
       ZFSUPAD = 0.
       ZFSDNAD = 0.
    endif
    !rv

    IF (appel1er) THEN
       PRINT*, 'SW calling frequency: ', swpas
       PRINT*, " In general, it should be 1"
       appel1er = .FALSE.
    ENDIF

    IF (MOD(itapsw, swpas) == 0) THEN
       DO JK = 1, LLM
          DO JL = 1, KLON
             ZCLDSW0(JL, JK) = 0.0
             ZOZ(JL, JK) = POZON(JL, JK) / (dobson_u * 1E3 * rg) * PDP(JL, JK)
          ENDDO
       ENDDO

       ! clear-sky:
       CALL SWU(PSCT, ZCLDSW0, PPMB, PPSOL, PRMU0, FRACT, PTAVE, PWV, ZAKI, &
            ZCLD, ZCLEAR, ZDSIG, ZFACT, ZRMU, ZSEC, ZUD)
       INU = 1
       CALL SW1S(INU, PALBD, PALBP, PCG, ZCLD, ZCLEAR, ZDSIG, POMEGA, ZOZ, &
            ZRMU, ZSEC, PTAU, ZUD, ZFD, ZFU)
       INU = 2
       CALL SW2S(INU, ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, ZDSIG, POMEGA, &
            ZOZ, ZRMU, ZSEC, PTAU, ZUD, PWV, PQS, ZFDOWN, ZFUP)
       DO JK = 1, LLM+1
          DO JL = 1, KLON
             ZFSUP0(JL, JK) = (ZFUP(JL, JK) + ZFU(JL, JK)) * ZFACT(JL)
             ZFSDN0(JL, JK) = (ZFDOWN(JL, JK) + ZFD(JL, JK)) * ZFACT(JL)
          ENDDO
       ENDDO

       CALL SWU(PSCT, PCLDSW, PPMB, PPSOL, PRMU0, FRACT, PTAVE, PWV, ZAKI, &
            ZCLD, ZCLEAR, ZDSIG, ZFACT, ZRMU, ZSEC, ZUD)
       INU = 1
       CALL SW1S(INU, PALBD, PALBP, PCG, ZCLD, ZCLEAR, ZDSIG, POMEGA, ZOZ, &
            ZRMU, ZSEC, PTAU, ZUD, ZFD, ZFU)
       INU = 2
       CALL SW2S(INU, ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, ZDSIG, POMEGA, &
            ZOZ, ZRMU, ZSEC, PTAU, ZUD, PWV, PQS, ZFDOWN, ZFUP)

       ! cloudy-sky:

       DO JK = 1, LLM+1
          DO JL = 1, KLON
             ZFSUP(JL, JK) = (ZFUP(JL, JK) + ZFU(JL, JK)) * ZFACT(JL)
             ZFSDN(JL, JK) = (ZFDOWN(JL, JK) + ZFD(JL, JK)) * ZFACT(JL)
          ENDDO
       ENDDO

       itapsw = 0
    ENDIF
    itapsw = itapsw + 1

    DO k = 1, LLM
       kpl1 = k+1
       DO i = 1, KLON
          PHEAT(i, k) = -(ZFSUP(i, kpl1)-ZFSUP(i, k)) &
               -(ZFSDN(i, k)-ZFSDN(i, kpl1))
          PHEAT(i, k) = PHEAT(i, k) * DAYSEC*RG/RCPD / PDP(i, k)
          PHEAT0(i, k) = -(ZFSUP0(i, kpl1)-ZFSUP0(i, k)) &
               -(ZFSDN0(i, k)-ZFSDN0(i, kpl1))
          PHEAT0(i, k) = PHEAT0(i, k) * DAYSEC*RG/RCPD / PDP(i, k)
       ENDDO
    ENDDO
    DO i = 1, KLON
       PSOLSW(i) = ZFSDN(i, 1) - ZFSUP(i, 1)
       PTOPSW(i) = ZFSDN(i, LLM+1) - ZFSUP(i, LLM+1)

       PSOLSW0(i) = ZFSDN0(i, 1) - ZFSUP0(i, 1)
       PTOPSW0(i) = ZFSDN0(i, LLM+1) - ZFSUP0(i, LLM+1)

       PSOLSWAD(i) = ZFSDNAD(i, 1) - ZFSUPAD(i, 1)
       PTOPSWAD(i) = ZFSDNAD(i, LLM+1) - ZFSUPAD(i, LLM+1)
    ENDDO

  END SUBROUTINE SW

end module sw_m
