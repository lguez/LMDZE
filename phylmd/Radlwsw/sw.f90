module sw_m

  IMPLICIT none

contains

  SUBROUTINE SW(PSCT, PRMU0, FRACT, PPMB, PDP, PPSOL, PALBD, PALBP, PTAVE, &
       PWV, PQS, POZON, PCLDSW, PTAU, POMEGA, PCG, PHEAT, PHEAT0, TOPSW, &
       PSOLSW, TOPSW0, PSOLSW0, ZFSUP, ZFSDN, ZFSUP0, ZFSDN0)

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

    DOUBLE PRECISION, intent(in):: PSCT ! constante solaire
    DOUBLE PRECISION, intent(in):: PRMU0(KLON) ! COSINE OF ZENITHAL ANGLE
    DOUBLE PRECISION, intent(in):: FRACT(KLON) ! fraction de la journee
    DOUBLE PRECISION, intent(in):: PPMB(KLON, LLM+1) ! HALF-LEVEL PRESSURE (MB)
    DOUBLE PRECISION, intent(in):: PDP(KLON, LLM) ! LAYER THICKNESS (PA)
    DOUBLE PRECISION, intent(in):: PPSOL(KLON) ! SURFACE PRESSURE (PA)
    DOUBLE PRECISION, intent(in):: PALBD(KLON, 2) ! albedo du sol (lumiere diffuse)
    DOUBLE PRECISION, intent(in):: PALBP(KLON, 2) ! albedo du sol (lumiere parallele)
    DOUBLE PRECISION, intent(in):: PTAVE(KLON, LLM) ! LAYER TEMPERATURE (K)
    DOUBLE PRECISION, intent(in):: PWV(KLON, LLM) ! SPECIFIC HUMIDITY (KG/KG)
    DOUBLE PRECISION, intent(in):: PQS(KLON, LLM) ! SATURATED WATER VAPOUR (KG/KG)
    DOUBLE PRECISION, intent(in):: POZON(KLON, LLM) ! OZONE CONCENTRATION (KG/KG)
    DOUBLE PRECISION, intent(in):: PCLDSW(KLON, LLM) ! CLOUD FRACTION
    DOUBLE PRECISION, intent(in):: PTAU(KLON, 2, LLM) ! CLOUD OPTICAL THICKNESS
    DOUBLE PRECISION, intent(in):: POMEGA(KLON, 2, LLM) ! SINGLE SCATTERING ALBEDO
    DOUBLE PRECISION, intent(in):: PCG(KLON, 2, LLM) ! ASYMETRY FACTOR
    DOUBLE PRECISION PHEAT(KLON, LLM) ! SHORTWAVE HEATING (K/DAY)
    DOUBLE PRECISION PHEAT0(KLON, LLM)! SHORTWAVE HEATING (K/DAY) clear-sky
    DOUBLE PRECISION, intent(out):: TOPSW(KLON) ! SHORTWAVE FLUX AT T.O.A.
    DOUBLE PRECISION PSOLSW(KLON) ! SHORTWAVE FLUX AT SURFACE
    DOUBLE PRECISION TOPSW0(KLON) ! SHORTWAVE FLUX AT T.O.A. (CLEAR-SKY)
    DOUBLE PRECISION PSOLSW0(KLON) ! SHORTWAVE FLUX AT SURFACE (CLEAR-SKY)
    DOUBLE PRECISION ZFSUP(KLON, LLM+1)
    DOUBLE PRECISION ZFSDN(KLON, LLM+1)
    DOUBLE PRECISION ZFSUP0(KLON, LLM+1)
    DOUBLE PRECISION ZFSDN0(KLON, LLM+1)

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
    REAL, PARAMETER :: dobson_u = 2.1415E-05 ! Dobson unit, in kg m-2

    !-------------------------------------------------------------------

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
       TOPSW(i) = ZFSDN(i, LLM+1) - ZFSUP(i, LLM+1)
       PSOLSW0(i) = ZFSDN0(i, 1) - ZFSUP0(i, 1)
       TOPSW0(i) = ZFSDN0(i, LLM+1) - ZFSUP0(i, LLM+1)
    ENDDO

  END SUBROUTINE SW

end module sw_m
