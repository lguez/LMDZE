module lw_m

  IMPLICIT none

contains

  SUBROUTINE LW(PMB, DP, DT0, EMIS, TL, TAVE, WV, OZON, CLDLD, &
       CLDLU, VIEW, COLR, COLR0, TOPLW, SOLLW, TOPLW0, SOLLW0, &
       sollwdown, lwup, lwdn, lwup0, lwdn0)

    use comconst, only: daysec
    use dimphy, only: klon
    use lwbv_m, only: lwbv
    use LWU_m, only: LWU
    use dimensions, only: llm
    USE raddimlw, ONLY: nua
    USE suphec_m, ONLY: md, rcpd, rg, rmo3

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

    DOUBLE PRECISION CLDLD(KLON, LLM) ! DOWNWARD EFFECTIVE CLOUD COVER
    DOUBLE PRECISION CLDLU(KLON, LLM) ! UPWARD EFFECTIVE CLOUD COVER
    DOUBLE PRECISION DP(KLON, LLM) ! LAYER PRESSURE THICKNESS (Pa)
    DOUBLE PRECISION DT0(KLON) ! SURFACE TEMPERATURE DISCONTINUITY (K)
    DOUBLE PRECISION EMIS(KLON) ! SURFACE EMISSIVITY
    DOUBLE PRECISION PMB(KLON, LLM+1) ! HALF LEVEL PRESSURE (mb)
    DOUBLE PRECISION OZON(KLON, LLM) ! O3 CONCENTRATION (kg/kg)
    DOUBLE PRECISION TL(KLON, LLM+1) ! HALF LEVEL TEMPERATURE (K)

    DOUBLE PRECISION TAVE(KLON, LLM) ! LAYER TEMPERATURE (K)
    DOUBLE PRECISION VIEW(KLON) ! COSECANT OF VIEWING ANGLE
    DOUBLE PRECISION WV(KLON, LLM) ! SPECIFIC HUMIDITY (kg/kg)

    DOUBLE PRECISION COLR(KLON, LLM) ! LONG-WAVE TENDENCY (K/day)
    DOUBLE PRECISION COLR0(KLON, LLM) ! LONG-WAVE TENDENCY (K/day) clear-sky
    DOUBLE PRECISION TOPLW(KLON) ! LONGWAVE FLUX AT T.O.A.
    DOUBLE PRECISION SOLLW(KLON) ! LONGWAVE FLUX AT SURFACE
    DOUBLE PRECISION TOPLW0(KLON) ! LONGWAVE FLUX AT T.O.A. (CLEAR-SKY)
    DOUBLE PRECISION SOLLW0(KLON) ! LONGWAVE FLUX AT SURFACE (CLEAR-SKY)
    ! Rajout LF
    double precision sollwdown(klon) ! LONGWAVE downwards flux at surface
    !IM
    DOUBLE PRECISION lwup(KLON, LLM+1) ! LW up total sky
    DOUBLE PRECISION lwup0(KLON, LLM+1) ! LW up clear sky
    DOUBLE PRECISION lwdn(KLON, LLM+1) ! LW down total sky
    DOUBLE PRECISION lwdn0(KLON, LLM+1) ! LW down clear sky

    DOUBLE PRECISION ZABCU(KLON, NUA, 3*LLM+1)
    DOUBLE PRECISION ZOZ(KLON, LLM)

    DOUBLE PRECISION, save, allocatable:: ZFLUX(:, :, :) ! (KLON, 2, LLM+1)
    ! RADIATIVE FLUXES (1:up; 2:down)

    DOUBLE PRECISION, save, allocatable:: ZFLUC(:, :, :) ! (KLON, 2, LLM+1)
    ! CLEAR-SKY RADIATIVE FLUXES

    ! Intermediate variables:
    DOUBLE PRECISION, save, allocatable:: ZBINT(:, :) ! (KLON, LLM+1)
    DOUBLE PRECISION, save, allocatable:: ZBSUI(:) ! (KLON)
    DOUBLE PRECISION, save, allocatable:: ZCTS(:, :) ! (KLON, LLM)

    DOUBLE PRECISION, save, allocatable:: ZCNTRB(:, :, :)
    ! (KLON, LLM+1, LLM+1)

    INTEGER ilim, i, k, kpl1

    INTEGER, PARAMETER:: lw0pas = 1 ! Every lw0pas steps, clear-sky is done
    INTEGER, PARAMETER:: lwpas = 1 ! Every lwpas steps, cloudy-sky is done
    ! In general, lw0pas and lwpas should be 1

    INTEGER:: itaplw0 = 0, itaplw = 0
    logical:: first_call = .true.

    ! ------------------------------------------------------------------

    if (first_call) then
       allocate(ZFLUX(KLON, 2, LLM+1), ZFLUC(KLON, 2, LLM+1), &
            ZBINT(KLON, LLM+1), ZBSUI(KLON), ZCTS(KLON, LLM), &
            ZCNTRB(KLON, LLM+1, LLM+1))
       first_call = .false.
    end if

    IF (MOD(itaplw0, lw0pas) == 0) THEN
       DO k = 1, LLM
          DO i = 1, KLON
             ! convertir ozone de kg/kg en pa (modif MPL 100505)
             ZOZ(i, k) = OZON(i, k)*DP(i, k) * MD/RMO3
          ENDDO
       ENDDO

       CALL LWU(DP, PMB, ZOZ, TAVE, VIEW, WV, ZABCU)
       CALL LWBV(ILIM, DT0, EMIS, PMB, TL, TAVE, ZABCU, &
            ZFLUC, ZBINT, ZBSUI, ZCTS, ZCNTRB)
       itaplw0 = 0
    ENDIF
    itaplw0 = itaplw0 + 1

    IF (MOD(itaplw, lwpas) == 0) THEN
       CALL LWC(ILIM, CLDLD, CLDLU, EMIS, &
            ZFLUC, ZBINT, ZBSUI, ZCTS, ZCNTRB, &
            ZFLUX)
       itaplw = 0
    ENDIF
    itaplw = itaplw + 1

    DO k = 1, LLM
       kpl1 = k+1
       DO i = 1, KLON
          COLR(i, k) = ZFLUX(i, 1, kpl1)+ZFLUX(i, 2, kpl1) &
               - ZFLUX(i, 1, k)- ZFLUX(i, 2, k)
          COLR(i, k) = COLR(i, k) * DAYSEC*RG/RCPD / DP(i, k)
          COLR0(i, k) = ZFLUC(i, 1, kpl1)+ZFLUC(i, 2, kpl1) &
               - ZFLUC(i, 1, k)- ZFLUC(i, 2, k)
          COLR0(i, k) = COLR0(i, k) * DAYSEC*RG/RCPD / DP(i, k)
       ENDDO
    ENDDO
    DO i = 1, KLON
       SOLLW(i) = -ZFLUX(i, 1, 1)-ZFLUX(i, 2, 1)
       TOPLW(i) = ZFLUX(i, 1, LLM+1) + ZFLUX(i, 2, LLM+1)

       SOLLW0(i) = -ZFLUC(i, 1, 1)-ZFLUC(i, 2, 1)
       TOPLW0(i) = ZFLUC(i, 1, LLM+1) + ZFLUC(i, 2, LLM+1)
       sollwdown(i) = -ZFLUX(i, 2, 1)

       ! Attention aux signes : LWtop > 0, LWdn < 0
       DO k = 1, LLM+1
          lwup(i, k) = ZFLUX(i, 1, k)
          lwup0(i, k) = ZFLUC(i, 1, k)
          lwdn(i, k) = ZFLUX(i, 2, k)
          lwdn0(i, k) = ZFLUC(i, 2, k)
       ENDDO
    ENDDO

  END SUBROUTINE LW

end module lw_m
