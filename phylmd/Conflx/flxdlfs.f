module flxdlfs_m

  IMPLICIT none

contains

  SUBROUTINE flxdlfs(ptenh, pqenh, pgeoh, paph, ptu, pqu, ldcum, kcbot, &
       kctop, pmfub, prfl, ptd, pqd, pmfd, pmfds, pmfdq, pdmfdp, kdtop, lddraf)

    ! This routine calculates level of free sinking for cumulus
    ! downdrafts and specifies T, q, u and v values

    ! To produce LFS-values for cumulus downdrafts for massflux
    ! cumulus parameterization

    ! Input are environmental values of T, q, u, v, p, Phi and updraft
    ! values T, q, u and v and also cloud base massflux and
    ! cu-precipitation rate.  it returns T, q, u and v values and
    ! massflux at LFS.

    ! Check for negative buoyancy of air of equal parts of moist
    ! environmental air and cloud air.

    USE dimphy, ONLY: klev, klon
    USE flxadjtq_m, ONLY: flxadjtq
    USE suphec_m, ONLY: rcpd, retv
    USE yoecumf, ONLY: cmfdeps

    REAL ptenh(klon, klev)
    REAL pqenh(klon, klev)
    REAL, intent(in):: pgeoh(klon, klev), paph(klon, klev + 1)
    REAL ptu(klon, klev), pqu(klon, klev)
    LOGICAL  ldcum(klon)
    INTEGER  kcbot(klon), kctop(klon)
    REAL pmfub(klon)
    REAL prfl(klon)
    REAL ptd(klon, klev), pqd(klon, klev)
    REAL pmfd(klon, klev), pmfds(klon, klev), pmfdq(klon, klev)
    REAL pdmfdp(klon, klev)
    INTEGER kdtop(klon)
    LOGICAL lddraf(klon)

    ! Local:
    REAL ztenwb(klon, klev), zqenwb(klon, klev), zcond(klon)
    REAL zttest, zqtest, zbuo, zmftop
    LOGICAL  llo2(klon)
    INTEGER i, k, is, icall

    !----------------------------------------------------------------------

    DO i= 1, klon
       lddraf(i)=.FALSE.
       kdtop(i)=klev + 1
    ENDDO

    ! DETERMINE LEVEL OF FREE SINKING BY
    ! DOING A SCAN FROM TOP TO BASE OF CUMULUS CLOUDS

    ! FOR EVERY POINT AND PROCEED AS FOLLOWS:
    !     (1) DETEMINE WET BULB ENVIRONMENTAL T AND Q
    !     (2) DO MIXING WITH CUMULUS CLOUD AIR
    !     (3) CHECK FOR NEGATIVE BUOYANCY

    ! THE ASSUMPTION IS THAT AIR OF DOWNDRAFTS IS MIXTURE
    ! OF 50% CLOUD AIR + 50% ENVIRONMENTAL AIR AT WET BULB
    ! TEMPERATURE (I.E. WHICH BECAME SATURATED DUE TO
    ! EVAPORATION OF RAIN AND CLOUD WATER)

    DO k = 3, klev-3
       is=0
       DO i= 1, klon
          ztenwb(i, k)=ptenh(i, k)
          zqenwb(i, k)=pqenh(i, k)
          llo2(i) = ldcum(i).AND.prfl(i).GT.0. &
               .AND..NOT.lddraf(i) &
               .AND.(k.LT.kcbot(i).AND.k.GT.kctop(i))
          IF ( llo2(i) ) is = is + 1
       end DO
       IF(is.EQ.0) cycle

       icall=2
       CALL flxadjtq(paph(:, k), ztenwb(1, k), zqenwb(1, k), llo2, icall)

       ! DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
       ! AND CHECK FOR NEGATIVE BUOYANCY.
       ! THEN SET VALUES FOR DOWNDRAFT AT LFS.

       DO i= 1, klon
          IF (llo2(i)) THEN
             zttest=0.5*(ptu(i, k) + ztenwb(i, k))
             zqtest=0.5*(pqu(i, k) + zqenwb(i, k))
             zbuo=zttest*(1. + RETV*zqtest)- &
                  ptenh(i, k)*(1. + RETV  *pqenh(i, k))
             zcond(i)=pqenh(i, k)-zqenwb(i, k)
             zmftop=-CMFDEPS*pmfub(i)
             IF (zbuo.LT.0..AND.prfl(i).GT.10.*zmftop*zcond(i)) THEN
                kdtop(i)=k
                lddraf(i)=.TRUE.
                ptd(i, k)=zttest
                pqd(i, k)=zqtest
                pmfd(i, k)=zmftop
                pmfds(i, k)=pmfd(i, k)*(RCPD*ptd(i, k) + pgeoh(i, k))
                pmfdq(i, k)=pmfd(i, k)*pqd(i, k)
                pdmfdp(i, k-1)=-0.5*pmfd(i, k)*zcond(i)
                prfl(i)=prfl(i) + pdmfdp(i, k-1)
             ENDIF
          ENDIF
       end DO
    end DO

  END SUBROUTINE flxdlfs

end module flxdlfs_m
