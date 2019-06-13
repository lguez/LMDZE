module flxbase_m

  IMPLICIT none

contains

  SUBROUTINE flxbase(ptenh, pqenh, pgeoh, paph, ptu, pqu, plu, ldcum, kcbot, &
       klab)

    ! This routine calculates cloud base values (T and Q).
    ! Input are environmental values of T, q, p, Phi at half levels.
    ! It returns cloud base values and flags as follows:
    ! klab = 1 for subcloud levels
    ! klab = 2 for condensation level

    ! Lift surface air dry-adiabatically to cloud base
    ! (non entraining plume, i. e. constant massflux)

    USE dimphy, ONLY: klev, klon
    USE flxadjtq_m, ONLY: flxadjtq
    USE suphec_m, ONLY: rcpd, retv

    REAL ptenh(klon, klev), pqenh(klon, klev)
    REAL, intent(in):: pgeoh(klon, klev), paph(klon, klev+1)
    REAL ptu(klon, klev), pqu(klon, klev), plu(klon, klev)
    LOGICAL ldcum(klon)
    INTEGER kcbot(klon), klab(klon, klev)

    ! Local:
    LOGICAL llflag(klon)
    INTEGER i, k, icall, is
    REAL zbuo, zqold(klon)

    !----------------------------------------------------------------------

    ! INITIALIZE VALUES AT LIFTING LEVEL
    DO i = 1, klon
       klab(i, klev) = 1
       kcbot(i) = klev-1
       ldcum(i) = .FALSE.
    ENDDO

    ! DO ASCENT IN SUBCLOUD LAYER,
    ! CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
    ! ADJUST T, Q AND L ACCORDINGLY
    ! CHECK FOR BUOYANCY AND SET FLAGS

    DO k = klev-1, 2, -1
       is = 0
       DO i = 1, klon
          IF (klab(i, k+1).EQ.1) is = is + 1
          llflag(i) = .FALSE.
          IF (klab(i, k+1).EQ.1) llflag(i) = .TRUE.
       ENDDO

       IF (is.EQ.0) cycle

       DO i = 1, klon
          IF(llflag(i)) THEN
             pqu(i, k) = pqu(i, k+1)
             ptu(i, k) = ptu(i, k+1)+(pgeoh(i, k+1)-pgeoh(i, k))/RCPD
             zbuo = ptu(i, k)*(1.+RETV*pqu(i, k))- &
                  ptenh(i, k)*(1.+RETV*pqenh(i, k))+0.5
             IF (zbuo.GT.0.) klab(i, k) = 1
             zqold(i) = pqu(i, k)
          ENDIF
       ENDDO

       icall = 1
       CALL flxadjtq(paph(:, k), ptu(1, k), pqu(1, k), llflag, icall)

       DO i = 1, klon
          IF (llflag(i).AND.pqu(i, k).NE.zqold(i)) THEN
             klab(i, k) = 2
             plu(i, k) = plu(i, k) + zqold(i)-pqu(i, k)
             zbuo = ptu(i, k)*(1.+RETV*pqu(i, k))- &
                  ptenh(i, k)*(1.+RETV*pqenh(i, k))+0.5
             IF (zbuo.GT.0.) kcbot(i) = k
             IF (zbuo.GT.0.) ldcum(i) = .TRUE.
          ENDIF
       ENDDO
    end DO

  END SUBROUTINE flxbase

end module flxbase_m
