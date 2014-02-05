module flxddraf_m

  IMPLICIT none

contains

  SUBROUTINE flxddraf(ptenh, pqenh, pgeoh, paph, prfl, ptd, pqd, pmfd, pmfds, &
       pmfdq, pdmfdp, lddraf, pen_d, pde_d)

    ! This routine calculates cumulus downdraft descent

    ! To produce the vertical profiles for cumulus downdrafts
    ! (i.e. T, q, u and v and fluxes)

    ! Input is T, q, p, Phi, u, v at half levels.
    ! It returns fluxes of s, q and evaporation rate
    ! and u, v at levels where downdraft occurs

    ! Calculate moist descent for entraining/detraining plume by
    ! A) moving air dry-adiabatically to next level below and
    ! B) correcting for evaporation to obtain saturated state.

    USE dimphy, ONLY: klev, klon
    USE flxadjtq_m, ONLY: flxadjtq
    USE suphec_m, ONLY: rcpd, rd, retv, rg
    USE yoecumf, ONLY: cmfcmin, entrdd

    REAL ptenh(klon, klev), pqenh(klon, klev)
    REAL, intent(in):: pgeoh(klon, klev), paph(klon, klev + 1)
    REAL prfl(klon)
    REAL ptd(klon, klev), pqd(klon, klev)
    REAL pmfd(klon, klev), pmfds(klon, klev), pmfdq(klon, klev)
    REAL pdmfdp(klon, klev)
    LOGICAL lddraf(klon)
    REAL pen_d(klon, klev), pde_d(klon, klev)

    ! Local:
    real zcond(klon)
    LOGICAL llo2(klon), llo1
    INTEGER i, k, is, icall, itopde
    REAL zentr, zseen, zqeen, zsdde, zqdde, zmfdsk, zmfdqk, zdmfdp
    REAL zbuo

    !----------------------------------------------------------------------

    ! Calculate moist descent for cumulus downdraft by

    ! (A) calculating entrainment rates, assuming linear decrease of
    ! massflux in PBL

    ! (B) doing moist descent - evaporative cooling and moistening is
    ! calculated in flxadjtq

    ! (C) checking for negative buoyancy and specifying final T, q, u,
    ! v and downward fluxes

    DO k = 3, klev
       is = 0
       DO i = 1, klon
          llo2(i)=lddraf(i).AND.pmfd(i, k-1).LT.0.
          IF (llo2(i)) is = is + 1
       ENDDO
       IF (is.EQ.0) cycle

       DO i = 1, klon
          IF (llo2(i)) THEN
             zentr = ENTRDD*pmfd(i, k-1)*RD*ptenh(i, k-1)/ &
                  (RG*paph(i, k-1))*(paph(i, k)-paph(i, k-1))
             pen_d(i, k) = zentr
             pde_d(i, k) = zentr
          ENDIF
       ENDDO

       itopde = klev-2
       IF (k.GT.itopde) THEN
          DO i = 1, klon
             IF (llo2(i)) THEN
                pen_d(i, k)=0.
                pde_d(i, k) = pmfd(i, itopde) * (paph(i, k) - paph(i, k - 1)) &
                     / (paph(i, klev + 1) - paph(i, itopde))
             ENDIF
          ENDDO
       ENDIF

       DO i = 1, klon
          IF (llo2(i)) THEN
             pmfd(i, k) = pmfd(i, k-1) + pen_d(i, k)-pde_d(i, k)
             zseen = (RCPD*ptenh(i, k-1) + pgeoh(i, k-1))*pen_d(i, k)
             zqeen = pqenh(i, k-1)*pen_d(i, k)
             zsdde = (RCPD*ptd(i, k-1) + pgeoh(i, k-1))*pde_d(i, k)
             zqdde = pqd(i, k-1)*pde_d(i, k)
             zmfdsk = pmfds(i, k-1) + zseen-zsdde
             zmfdqk = pmfdq(i, k-1) + zqeen-zqdde
             pqd(i, k) = zmfdqk*(1./MIN(-CMFCMIN, pmfd(i, k)))
             ptd(i, k) = (zmfdsk*(1./MIN(-CMFCMIN, pmfd(i, k)))- &
                  pgeoh(i, k))/RCPD
             ptd(i, k) = MIN(400., ptd(i, k))
             ptd(i, k) = MAX(100., ptd(i, k))
             zcond(i) = pqd(i, k)
          ENDIF
       ENDDO

       icall = 2
       CALL flxadjtq(paph(:, k), ptd(1, k), pqd(1, k), llo2, icall)

       DO i = 1, klon
          IF (llo2(i)) THEN
             zcond(i) = zcond(i)-pqd(i, k)
             zbuo = ptd(i, k)*(1. + RETV *pqd(i, k))- &
                  ptenh(i, k)*(1. + RETV *pqenh(i, k))
             llo1 = zbuo.LT.0..AND.(prfl(i)-pmfd(i, k)*zcond(i).GT.0.)
             IF (.not.llo1) pmfd(i, k) = 0.0
             pmfds(i, k) = (RCPD*ptd(i, k) + pgeoh(i, k))*pmfd(i, k)
             pmfdq(i, k) = pqd(i, k)*pmfd(i, k)
             zdmfdp = -pmfd(i, k)*zcond(i)
             pdmfdp(i, k-1) = zdmfdp
             prfl(i) = prfl(i) + zdmfdp
          ENDIF
       ENDDO
    end DO

  END SUBROUTINE flxddraf

end module flxddraf_m
