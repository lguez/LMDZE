module sw_m

  IMPLICIT none

contains

  SUBROUTINE SW(PSCT, PRMU0, PFRAC, PPMB, PDP, PPSOL, PALBD, PALBP, PTAVE, &
       PWV, PQS, POZON, PAER, PCLDSW, PTAU, POMEGA, PCG, PHEAT, PHEAT0, &
       PALBPLA, PTOPSW, PSOLSW, PTOPSW0, PSOLSW0, ZFSUP, ZFSDN, ZFSUP0, &
       ZFSDN0, tauae, pizae, cgae, PTAUA, POMEGAA, PTOPSWAD, PSOLSWAD, &
       PTOPSWAI, PSOLSWAI, ok_ade, ok_aie)

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

    USE raddim, ONLY: kdlon, kflev
    USE suphec_m, ONLY: rcpd, rday, rg

    ! ARGUMENTS:

    DOUBLE PRECISION PSCT ! constante solaire (valeur conseillee: 1370)

    DOUBLE PRECISION PPSOL(KDLON) ! SURFACE PRESSURE (PA)
    DOUBLE PRECISION PDP(KDLON, KFLEV) ! LAYER THICKNESS (PA)
    DOUBLE PRECISION PPMB(KDLON, KFLEV+1) ! HALF-LEVEL PRESSURE (MB)

    DOUBLE PRECISION PRMU0(KDLON) ! COSINE OF ZENITHAL ANGLE
    DOUBLE PRECISION PFRAC(KDLON) ! fraction de la journee

    DOUBLE PRECISION PTAVE(KDLON, KFLEV) ! LAYER TEMPERATURE (K)
    DOUBLE PRECISION PWV(KDLON, KFLEV) ! SPECIFIC HUMIDITY (KG/KG)
    DOUBLE PRECISION PQS(KDLON, KFLEV) ! SATURATED WATER VAPOUR (KG/KG)
    DOUBLE PRECISION POZON(KDLON, KFLEV) ! OZONE CONCENTRATION (KG/KG)
    DOUBLE PRECISION PAER(KDLON, KFLEV, 5) ! AEROSOLS' OPTICAL THICKNESS

    DOUBLE PRECISION PALBD(KDLON, 2) ! albedo du sol (lumiere diffuse)
    DOUBLE PRECISION PALBP(KDLON, 2) ! albedo du sol (lumiere parallele)

    DOUBLE PRECISION PCLDSW(KDLON, KFLEV) ! CLOUD FRACTION
    DOUBLE PRECISION PTAU(KDLON, 2, KFLEV) ! CLOUD OPTICAL THICKNESS
    DOUBLE PRECISION PCG(KDLON, 2, KFLEV) ! ASYMETRY FACTOR
    DOUBLE PRECISION POMEGA(KDLON, 2, KFLEV) ! SINGLE SCATTERING ALBEDO

    DOUBLE PRECISION PHEAT(KDLON, KFLEV) ! SHORTWAVE HEATING (K/DAY)
    DOUBLE PRECISION PHEAT0(KDLON, KFLEV)! SHORTWAVE HEATING (K/DAY) clear-sky
    DOUBLE PRECISION PALBPLA(KDLON) ! PLANETARY ALBEDO
    DOUBLE PRECISION PTOPSW(KDLON) ! SHORTWAVE FLUX AT T.O.A.
    DOUBLE PRECISION PSOLSW(KDLON) ! SHORTWAVE FLUX AT SURFACE
    DOUBLE PRECISION PTOPSW0(KDLON) ! SHORTWAVE FLUX AT T.O.A. (CLEAR-SKY)
    DOUBLE PRECISION PSOLSW0(KDLON) ! SHORTWAVE FLUX AT SURFACE (CLEAR-SKY)

    ! LOCAL VARIABLES:

    DOUBLE PRECISION ZOZ(KDLON, KFLEV)
    DOUBLE PRECISION ZAKI(KDLON, 2) 
    DOUBLE PRECISION ZCLD(KDLON, KFLEV)
    DOUBLE PRECISION ZCLEAR(KDLON) 
    DOUBLE PRECISION ZDSIG(KDLON, KFLEV)
    DOUBLE PRECISION ZFACT(KDLON)
    DOUBLE PRECISION ZFD(KDLON, KFLEV+1)
    DOUBLE PRECISION ZFDOWN(KDLON, KFLEV+1)
    DOUBLE PRECISION ZFU(KDLON, KFLEV+1)
    DOUBLE PRECISION ZFUP(KDLON, KFLEV+1)
    DOUBLE PRECISION ZRMU(KDLON)
    DOUBLE PRECISION ZSEC(KDLON)
    DOUBLE PRECISION ZUD(KDLON, 5, KFLEV+1)
    DOUBLE PRECISION ZCLDSW0(KDLON, KFLEV)

    DOUBLE PRECISION ZFSUP(KDLON, KFLEV+1)
    DOUBLE PRECISION ZFSDN(KDLON, KFLEV+1)
    DOUBLE PRECISION ZFSUP0(KDLON, KFLEV+1)
    DOUBLE PRECISION ZFSDN0(KDLON, KFLEV+1)

    INTEGER inu, jl, jk, i, k, kpl1

    INTEGER, PARAMETER:: swpas = 1 ! Every swpas steps, sw is calculated

    INTEGER:: itapsw = 0
    LOGICAL:: appel1er = .TRUE.
    !jq-Introduced for aerosol forcings
    double precision, save:: flag_aer
    logical, intent(in):: ok_ade, ok_aie ! use aerosol forcings or not?
    double precision tauae(kdlon, kflev, 2) ! aerosol optical properties
    double precision pizae(kdlon, kflev, 2) 
    ! aerosol optical properties(see aeropt.F)
    
    double precision cgae(kdlon, kflev, 2) !aerosol optical properties -"-
    DOUBLE PRECISION PTAUA(KDLON, 2, KFLEV) 
    ! CLOUD OPTICAL THICKNESS (pre-industrial value)

    DOUBLE PRECISION POMEGAA(KDLON, 2, KFLEV) ! SINGLE SCATTERING ALBEDO
    DOUBLE PRECISION PTOPSWAD(KDLON) 
    ! (diagnosed aerosol forcing)SHORTWAVE FLUX AT T.O.A.(+AEROSOL DIR)

    DOUBLE PRECISION PSOLSWAD(KDLON) 
    ! (diagnosed aerosol forcing)SHORTWAVE FLUX AT SURFACE(+AEROSOL DIR)

    DOUBLE PRECISION PTOPSWAI(KDLON) 
    ! (diagnosed aerosol forcing)SHORTWAVE FLUX AT T.O.A.(+AEROSOL IND)

    DOUBLE PRECISION PSOLSWAI(KDLON) 
    ! (diagnosed aerosol forcing)SHORTWAVE FLUX AT SURFACE(+AEROSOL IND)

    !jq - Fluxes including aerosol effects
    DOUBLE PRECISION, save:: ZFSUPAD(KDLON, KFLEV+1)
    DOUBLE PRECISION, save:: ZFSDNAD(KDLON, KFLEV+1)
    DOUBLE PRECISION, save:: ZFSUPAI(KDLON, KFLEV+1)
    DOUBLE PRECISION, save:: ZFSDNAI(KDLON, KFLEV+1)

    logical:: initialized = .false.
    REAL, PARAMETER :: dobson_u = 2.1415E-05 ! Dobson unit, in kg m-2

    !-------------------------------------------------------------------

    if(.not.initialized) then
       flag_aer=0.
       initialized=.TRUE.
       ZFSUPAD = 0.
       ZFSDNAD = 0.
       ZFSUPAI = 0.
       ZFSDNAI = 0.
    endif
    !rv

    IF (appel1er) THEN
       PRINT*, 'SW calling frequency: ', swpas
       PRINT*, " In general, it should be 1"
       appel1er = .FALSE.
    ENDIF

    IF (MOD(itapsw, swpas).EQ.0) THEN
       DO JK = 1 , KFLEV
          DO JL = 1, KDLON
             ZCLDSW0(JL, JK) = 0.0
             ZOZ(JL, JK) = POZON(JL, JK) / (dobson_u * 1E3 * rg) * PDP(JL, JK)
          ENDDO
       ENDDO

       ! clear-sky:
       CALL SWU(PSCT, ZCLDSW0, PPMB, PPSOL, &
            PRMU0, PFRAC, PTAVE, PWV, &
            ZAKI, ZCLD, ZCLEAR, ZDSIG, ZFACT, ZRMU, ZSEC, ZUD)
       INU = 1
       CALL SW1S(INU, &
            PAER, flag_aer, tauae, pizae, cgae, &
            PALBD, PALBP, PCG, ZCLD, ZCLEAR, ZCLDSW0, &
            ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD, &
            ZFD, ZFU)
       INU = 2
       CALL SW2S(INU, &
            PAER, flag_aer, tauae, pizae, cgae, &
            ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, ZCLDSW0, &
            ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD, &
            PWV, PQS, &
            ZFDOWN, ZFUP)
       DO JK = 1 , KFLEV+1
          DO JL = 1, KDLON
             ZFSUP0(JL, JK) = (ZFUP(JL, JK) + ZFU(JL, JK)) * ZFACT(JL)
             ZFSDN0(JL, JK) = (ZFDOWN(JL, JK) + ZFD(JL, JK)) * ZFACT(JL)
          ENDDO
       ENDDO

       flag_aer=0.
       CALL SWU(PSCT, PCLDSW, PPMB, PPSOL, &
            PRMU0, PFRAC, PTAVE, PWV, &
            ZAKI, ZCLD, ZCLEAR, ZDSIG, ZFACT, ZRMU, ZSEC, ZUD)
       INU = 1
       CALL SW1S(INU, &
            PAER, flag_aer, tauae, pizae, cgae, &
            PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW, &
            ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD, &
            ZFD, ZFU)
       INU = 2
       CALL SW2S(INU, &
            PAER, flag_aer, tauae, pizae, cgae, &
            ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW, &
            ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD, &
            PWV, PQS, &
            ZFDOWN, ZFUP)

       ! cloudy-sky:

       DO JK = 1 , KFLEV+1
          DO JL = 1, KDLON
             ZFSUP(JL, JK) = (ZFUP(JL, JK) + ZFU(JL, JK)) * ZFACT(JL)
             ZFSDN(JL, JK) = (ZFDOWN(JL, JK) + ZFD(JL, JK)) * ZFACT(JL)
          ENDDO
       ENDDO

       IF (ok_ade) THEN
          ! cloudy-sky + aerosol dir OB
          flag_aer=1.
          CALL SWU(PSCT, PCLDSW, PPMB, PPSOL, &
               PRMU0, PFRAC, PTAVE, PWV, &
               ZAKI, ZCLD, ZCLEAR, ZDSIG, ZFACT, ZRMU, ZSEC, ZUD)
          INU = 1
          CALL SW1S(INU, &
               PAER, flag_aer, tauae, pizae, cgae, &
               PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW, &
               ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD, &
               ZFD, ZFU)
          INU = 2
          CALL SW2S(INU, &
               PAER, flag_aer, tauae, pizae, cgae, &
               ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW, &
               ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD, &
               PWV, PQS, &
               ZFDOWN, ZFUP)
          DO JK = 1 , KFLEV+1
             DO JL = 1, KDLON
                ZFSUPAD(JL, JK) = ZFSUP(JL, JK) 
                ZFSDNAD(JL, JK) = ZFSDN(JL, JK) 
                ZFSUP(JL, JK) = (ZFUP(JL, JK) + ZFU(JL, JK)) * ZFACT(JL)
                ZFSDN(JL, JK) = (ZFDOWN(JL, JK) + ZFD(JL, JK)) * ZFACT(JL)
             ENDDO
          ENDDO
       ENDIF

       IF (ok_aie) THEN
          !jq cloudy-sky + aerosol direct + aerosol indirect
          flag_aer=1.0
          CALL SWU(PSCT, PCLDSW, PPMB, PPSOL, &
               PRMU0, PFRAC, PTAVE, PWV, &
               ZAKI, ZCLD, ZCLEAR, ZDSIG, ZFACT, ZRMU, ZSEC, ZUD)
          INU = 1
          CALL SW1S(INU, &
               PAER, flag_aer, tauae, pizae, cgae, &
               PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW, &
               ZDSIG, POMEGAA, ZOZ, ZRMU, ZSEC, PTAUA, ZUD, &
               ZFD, ZFU)
          INU = 2
          CALL SW2S(INU, &
               PAER, flag_aer, tauae, pizae, cgae, &
               ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW, &
               ZDSIG, POMEGAA, ZOZ, ZRMU, ZSEC, PTAUA, ZUD, &
               PWV, PQS, &
               ZFDOWN, ZFUP)
          DO JK = 1 , KFLEV+1
             DO JL = 1, KDLON
                ZFSUPAI(JL, JK) = ZFSUP(JL, JK) 
                ZFSDNAI(JL, JK) = ZFSDN(JL, JK) 
                ZFSUP(JL, JK) = (ZFUP(JL, JK) + ZFU(JL, JK)) * ZFACT(JL)
                ZFSDN(JL, JK) = (ZFDOWN(JL, JK) + ZFD(JL, JK)) * ZFACT(JL)
             ENDDO
          ENDDO
       ENDIF

       itapsw = 0
    ENDIF
    itapsw = itapsw + 1

    DO k = 1, KFLEV
       kpl1 = k+1
       DO i = 1, KDLON
          PHEAT(i, k) = -(ZFSUP(i, kpl1)-ZFSUP(i, k)) &
               -(ZFSDN(i, k)-ZFSDN(i, kpl1))
          PHEAT(i, k) = PHEAT(i, k) * RDAY*RG/RCPD / PDP(i, k)
          PHEAT0(i, k) = -(ZFSUP0(i, kpl1)-ZFSUP0(i, k)) &
               -(ZFSDN0(i, k)-ZFSDN0(i, kpl1))
          PHEAT0(i, k) = PHEAT0(i, k) * RDAY*RG/RCPD / PDP(i, k)
       ENDDO
    ENDDO
    DO i = 1, KDLON
       PALBPLA(i) = ZFSUP(i, KFLEV+1)/(ZFSDN(i, KFLEV+1)+1.0e-20)

       PSOLSW(i) = ZFSDN(i, 1) - ZFSUP(i, 1)
       PTOPSW(i) = ZFSDN(i, KFLEV+1) - ZFSUP(i, KFLEV+1)

       PSOLSW0(i) = ZFSDN0(i, 1) - ZFSUP0(i, 1)
       PTOPSW0(i) = ZFSDN0(i, KFLEV+1) - ZFSUP0(i, KFLEV+1)

       PSOLSWAD(i) = ZFSDNAD(i, 1) - ZFSUPAD(i, 1)
       PTOPSWAD(i) = ZFSDNAD(i, KFLEV+1) - ZFSUPAD(i, KFLEV+1)

       PSOLSWAI(i) = ZFSDNAI(i, 1) - ZFSUPAI(i, 1)
       PTOPSWAI(i) = ZFSDNAI(i, KFLEV+1) - ZFSUPAI(i, KFLEV+1)
    ENDDO

  END SUBROUTINE SW

end module sw_m
