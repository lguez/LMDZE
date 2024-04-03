module lw_m

  IMPLICIT none

contains

  SUBROUTINE LW(PMB, DP, DT0, EMIS, TL, TAVE, WV, OZON, CLDLD, CLDLU, VIEW, &
       COLR, COLR0, TOPLW, SOLLW, TOPLW0, SOLLW0, sollwdown, lwup, lwdn, &
       lwup0, lwdn0)

    use comconst, only: daysec
    use dimensions, only: llm
    use dimphy, only: klon
    use lwbv_m, only: lwbv
    use lwc_m, only: lwc
    use LWU_m, only: LWU
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
    ! Author: Jean-Jacques Morcrette ECMWF
    ! Original: July 14th, 1989

    DOUBLE PRECISION PMB(KLON, LLM+1) ! HALF LEVEL PRESSURE (mb)
    DOUBLE PRECISION DP(KLON, LLM) ! LAYER PRESSURE THICKNESS (Pa)

    DOUBLE PRECISION, intent(in):: DT0(KLON)
    ! SURFACE TEMPERATURE DISCONTINUITY (K)

    DOUBLE PRECISION EMIS(KLON) ! SURFACE EMISSIVITY
    DOUBLE PRECISION TL(KLON, LLM+1) ! HALF LEVEL TEMPERATURE (K)
    DOUBLE PRECISION TAVE(KLON, LLM) ! LAYER TEMPERATURE (K)
    DOUBLE PRECISION WV(KLON, LLM) ! SPECIFIC HUMIDITY (kg/kg)
    DOUBLE PRECISION OZON(KLON, LLM) ! O3 CONCENTRATION (kg/kg)

    DOUBLE PRECISION, intent(in):: CLDLD(KLON, LLM)
    ! DOWNWARD EFFECTIVE CLOUD COVER

    DOUBLE PRECISION, intent(in):: CLDLU(KLON, LLM)
    ! UPWARD EFFECTIVE CLOUD COVER

    DOUBLE PRECISION VIEW(KLON) ! COSECANT OF VIEWING ANGLE
    DOUBLE PRECISION COLR(KLON, LLM) ! LONG-WAVE TENDENCY (K/day)
    DOUBLE PRECISION COLR0(KLON, LLM) ! LONG-WAVE TENDENCY (K/day) clear-sky
    DOUBLE PRECISION TOPLW(KLON) ! LONGWAVE FLUX AT T.O.A.
    DOUBLE PRECISION SOLLW(KLON) ! LONGWAVE FLUX AT SURFACE
    DOUBLE PRECISION TOPLW0(KLON) ! LONGWAVE FLUX AT T.O.A. (CLEAR-SKY)
    DOUBLE PRECISION SOLLW0(KLON) ! LONGWAVE FLUX AT SURFACE (CLEAR-SKY)
    double precision sollwdown(klon) ! LONGWAVE downwards flux at surface
    DOUBLE PRECISION lwup(KLON, LLM+1) ! LW up total sky
    DOUBLE PRECISION lwdn(KLON, LLM+1) ! LW down total sky
    DOUBLE PRECISION lwup0(KLON, LLM+1) ! LW up clear sky
    DOUBLE PRECISION lwdn0(KLON, LLM+1) ! LW down clear sky

    ! Local:

    DOUBLE PRECISION ZABCU(KLON, NUA, 3*LLM+1)
    DOUBLE PRECISION ZOZ(KLON, LLM)
    DOUBLE PRECISION ZFLUX(KLON, 2, LLM+1) ! RADIATIVE FLUXES (1:up; 2:down)
    DOUBLE PRECISION ZFLUC(KLON, 2, LLM+1) ! CLEAR-SKY RADIATIVE FLUXES

    ! Intermediate variables:
    DOUBLE PRECISION ZBINT(KLON, LLM+1)
    DOUBLE PRECISION ZBSUI(KLON)
    DOUBLE PRECISION ZCTS(KLON, LLM)

    DOUBLE PRECISION ZCNTRB(KLON, LLM+1, LLM+1)
    INTEGER ilim, i, k, kpl1

    !------------------------------------------------------------------

    DO k = 1, LLM
       DO i = 1, KLON
          ! convertir ozone de kg/kg en pa
          ZOZ(i, k) = OZON(i, k)*DP(i, k) * MD/RMO3
       ENDDO
    ENDDO

    CALL LWU(DP, PMB, ZOZ, TAVE, VIEW, WV, ZABCU)
    CALL LWBV(ILIM, DT0, EMIS, PMB, TL, TAVE, ZABCU, ZFLUC, ZBINT, ZBSUI, &
         ZCTS, ZCNTRB)
    CALL LWC(ILIM, CLDLD, CLDLU, EMIS, ZFLUC, ZBINT, ZBSUI, ZCTS, ZCNTRB, ZFLUX)

    DO k = 1, LLM
       kpl1 = k+1

       DO i = 1, KLON
          COLR(i, k) = ZFLUX(i, 1, kpl1)+ZFLUX(i, 2, kpl1) &
               - ZFLUX(i, 1, k) - ZFLUX(i, 2, k)
          COLR(i, k) = COLR(i, k) * DAYSEC*RG/RCPD / DP(i, k)
          COLR0(i, k) = ZFLUC(i, 1, kpl1)+ZFLUC(i, 2, kpl1) &
               - ZFLUC(i, 1, k) - ZFLUC(i, 2, k)
          COLR0(i, k) = COLR0(i, k) * DAYSEC*RG/RCPD / DP(i, k)
       ENDDO
    ENDDO

    DO i = 1, KLON
       SOLLW(i) = - ZFLUX(i, 1, 1) - ZFLUX(i, 2, 1)
       TOPLW(i) = ZFLUX(i, 1, LLM+1) + ZFLUX(i, 2, LLM+1)
       SOLLW0(i) = - ZFLUC(i, 1, 1) - ZFLUC(i, 2, 1)
       TOPLW0(i) = ZFLUC(i, 1, LLM+1) + ZFLUC(i, 2, LLM+1)
       sollwdown(i) = - ZFLUX(i, 2, 1)

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
