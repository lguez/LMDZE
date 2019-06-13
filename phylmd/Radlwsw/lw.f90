module lw_m

  IMPLICIT none

contains

  SUBROUTINE LW(PPMB, PDP, PDT0, PEMIS, PTL, PTAVE, PWV, POZON, PAER, PCLDLD, &
       PCLDLU, PVIEW, PCOLR, PCOLR0, PTOPLW, PSOLLW, PTOPLW0, PSOLLW0, &
       psollwdown, plwup, plwdn, plwup0, plwdn0)

    use lwbv_m, only: lwbv
    use LWU_m, only: LWU
    USE raddim, ONLY: kdlon, kflev
    USE raddimlw, ONLY: nua
    USE suphec_m, ONLY: md, rcpd, rday, rg, rmo3

    ! Method.

    ! 1. Computes the pressure and temperature weighted amounts of
    ! absorbers.

    ! 2. Computes the Planck functions on the interfaces and the
    ! gradient of Planck functions in the layers.

    ! 3. Performs the vertical integration distinguishing the
    ! contributions of the adjacent and distant layers and those from
    ! the boundaries.

    ! 4. Computes the clear-sky downward and upward emissivities.

    ! 5. Introduces the effects of the clouds on the fluxes.

    ! Reference: see radiation part of ECMWF documentation of the IFS.

    ! Author:
    ! Jean-Jacques Morcrette ECMWF

    ! Original : July 14th, 1989

    DOUBLE PRECISION PCLDLD(KDLON, KFLEV) ! DOWNWARD EFFECTIVE CLOUD COVER
    DOUBLE PRECISION PCLDLU(KDLON, KFLEV) ! UPWARD EFFECTIVE CLOUD COVER
    DOUBLE PRECISION PDP(KDLON, KFLEV) ! LAYER PRESSURE THICKNESS (Pa)
    DOUBLE PRECISION PDT0(KDLON) ! SURFACE TEMPERATURE DISCONTINUITY (K)
    DOUBLE PRECISION PEMIS(KDLON) ! SURFACE EMISSIVITY
    DOUBLE PRECISION PPMB(KDLON, KFLEV+1) ! HALF LEVEL PRESSURE (mb)
    DOUBLE PRECISION POZON(KDLON, KFLEV) ! O3 CONCENTRATION (kg/kg)
    DOUBLE PRECISION PTL(KDLON, KFLEV+1) ! HALF LEVEL TEMPERATURE (K)
    DOUBLE PRECISION PAER(KDLON, KFLEV, 5) ! OPTICAL THICKNESS OF THE AEROSOLS
    DOUBLE PRECISION PTAVE(KDLON, KFLEV) ! LAYER TEMPERATURE (K)
    DOUBLE PRECISION PVIEW(KDLON) ! COSECANT OF VIEWING ANGLE
    DOUBLE PRECISION PWV(KDLON, KFLEV) ! SPECIFIC HUMIDITY (kg/kg)

    DOUBLE PRECISION PCOLR(KDLON, KFLEV) ! LONG-WAVE TENDENCY (K/day)
    DOUBLE PRECISION PCOLR0(KDLON, KFLEV) ! LONG-WAVE TENDENCY (K/day) clear-sky
    DOUBLE PRECISION PTOPLW(KDLON) ! LONGWAVE FLUX AT T.O.A.
    DOUBLE PRECISION PSOLLW(KDLON) ! LONGWAVE FLUX AT SURFACE
    DOUBLE PRECISION PTOPLW0(KDLON) ! LONGWAVE FLUX AT T.O.A. (CLEAR-SKY)
    DOUBLE PRECISION PSOLLW0(KDLON) ! LONGWAVE FLUX AT SURFACE (CLEAR-SKY)
    ! Rajout LF
    double precision psollwdown(kdlon) ! LONGWAVE downwards flux at surface
    !IM
    DOUBLE PRECISION plwup(KDLON, KFLEV+1) ! LW up total sky
    DOUBLE PRECISION plwup0(KDLON, KFLEV+1) ! LW up clear sky
    DOUBLE PRECISION plwdn(KDLON, KFLEV+1) ! LW down total sky
    DOUBLE PRECISION plwdn0(KDLON, KFLEV+1) ! LW down clear sky

    DOUBLE PRECISION ZABCU(KDLON, NUA, 3*KFLEV+1)
    DOUBLE PRECISION ZOZ(KDLON, KFLEV)

    DOUBLE PRECISION, save:: ZFLUX(KDLON, 2, KFLEV+1)
    ! RADIATIVE FLUXES (1:up; 2:down)

    DOUBLE PRECISION, save:: ZFLUC(KDLON, 2, KFLEV+1)
    ! CLEAR-SKY RADIATIVE FLUXES

    ! Intermediate variables:
    DOUBLE PRECISION, save:: ZBINT(KDLON, KFLEV+1)
    DOUBLE PRECISION, save:: ZBSUI(KDLON)
    DOUBLE PRECISION, save:: ZCTS(KDLON, KFLEV)
    DOUBLE PRECISION, save:: ZCNTRB(KDLON, KFLEV+1, KFLEV+1)

    INTEGER ilim, i, k, kpl1

    INTEGER, PARAMETER:: lw0pas = 1 ! Every lw0pas steps, clear-sky is done
    INTEGER, PARAMETER:: lwpas = 1 ! Every lwpas steps, cloudy-sky is done
    ! In general, lw0pas and lwpas should be 1

    INTEGER:: itaplw0 = 0, itaplw = 0

    ! ------------------------------------------------------------------

    IF (MOD(itaplw0, lw0pas) == 0) THEN
       DO k = 1, KFLEV
          DO i = 1, KDLON
             ! convertir ozone de kg/kg en pa (modif MPL 100505)
             ZOZ(i, k) = POZON(i, k)*PDP(i, k) * MD/RMO3
          ENDDO
       ENDDO
       CALL LWU(PAER, PDP, PPMB, ZOZ, PTAVE, PVIEW, PWV, ZABCU)
       CALL LWBV(ILIM, PDT0, PEMIS, PPMB, PTL, PTAVE, ZABCU, &
            ZFLUC, ZBINT, ZBSUI, ZCTS, ZCNTRB)
       itaplw0 = 0
    ENDIF
    itaplw0 = itaplw0 + 1

    IF (MOD(itaplw, lwpas) == 0) THEN
       CALL LWC(ILIM, PCLDLD, PCLDLU, PEMIS, &
            ZFLUC, ZBINT, ZBSUI, ZCTS, ZCNTRB, &
            ZFLUX)
       itaplw = 0
    ENDIF
    itaplw = itaplw + 1

    DO k = 1, KFLEV
       kpl1 = k+1
       DO i = 1, KDLON
          PCOLR(i, k) = ZFLUX(i, 1, kpl1)+ZFLUX(i, 2, kpl1) &
               - ZFLUX(i, 1, k)- ZFLUX(i, 2, k)
          PCOLR(i, k) = PCOLR(i, k) * RDAY*RG/RCPD / PDP(i, k)
          PCOLR0(i, k) = ZFLUC(i, 1, kpl1)+ZFLUC(i, 2, kpl1) &
               - ZFLUC(i, 1, k)- ZFLUC(i, 2, k)
          PCOLR0(i, k) = PCOLR0(i, k) * RDAY*RG/RCPD / PDP(i, k)
       ENDDO
    ENDDO
    DO i = 1, KDLON
       PSOLLW(i) = -ZFLUX(i, 1, 1)-ZFLUX(i, 2, 1)
       PTOPLW(i) = ZFLUX(i, 1, KFLEV+1) + ZFLUX(i, 2, KFLEV+1)

       PSOLLW0(i) = -ZFLUC(i, 1, 1)-ZFLUC(i, 2, 1)
       PTOPLW0(i) = ZFLUC(i, 1, KFLEV+1) + ZFLUC(i, 2, KFLEV+1)
       psollwdown(i) = -ZFLUX(i, 2, 1)

       !IM attention aux signes !; LWtop >0, LWdn < 0
       DO k = 1, KFLEV+1
          plwup(i, k) = ZFLUX(i, 1, k)
          plwup0(i, k) = ZFLUC(i, 1, k)
          plwdn(i, k) = ZFLUX(i, 2, k)
          plwdn0(i, k) = ZFLUC(i, 2, k)
       ENDDO
    ENDDO

  END SUBROUTINE LW

end module lw_m
