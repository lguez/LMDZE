      SUBROUTINE flxbase(ptenh, pqenh, pgeoh, paph, &
           ptu, pqu, plu, ldcum, kcbot, klab)
      use dimens_m
      use dimphy
      use SUPHEC_M
      use yoethf_m
      IMPLICIT none
!----------------------------------------------------------------------
! THIS ROUTINE CALCULATES CLOUD BASE VALUES (T AND Q)
!
! INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
! IT RETURNS CLOUD BASE VALUES AND FLAGS AS FOLLOWS;
!   klab=1 FOR SUBCLOUD LEVELS
!   klab=2 FOR CONDENSATION LEVEL
!
! LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
! (NON ENTRAINING PLUME,I.E.CONSTANT MASSFLUX)
!----------------------------------------------------------------------
!       ----------------------------------------------------------------
      REAL ptenh(klon,klev), pqenh(klon,klev)
      REAL pgeoh(klon,klev), paph(klon,klev+1)
!
      REAL ptu(klon,klev), pqu(klon,klev), plu(klon,klev)
      INTEGER  klab(klon,klev), kcbot(klon)
!
      LOGICAL llflag(klon), ldcum(klon)
      INTEGER i, k, icall, is
      REAL zbuo, zqold(klon)
!----------------------------------------------------------------------
! INITIALIZE VALUES AT LIFTING LEVEL
!----------------------------------------------------------------------
      DO i = 1, klon
         klab(i,klev)=1
         kcbot(i)=klev-1
         ldcum(i)=.FALSE.
      ENDDO
!----------------------------------------------------------------------
! DO ASCENT IN SUBCLOUD LAYER,
! CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
! ADJUST T,Q AND L ACCORDINGLY
! CHECK FOR BUOYANCY AND SET FLAGS
!----------------------------------------------------------------------
      DO 290 k = klev-1, 2, -1
!
      is = 0
      DO i = 1, klon
         IF (klab(i,k+1).EQ.1) is = is + 1
         llflag(i) = .FALSE.
         IF (klab(i,k+1).EQ.1) llflag(i) = .TRUE.
      ENDDO
      IF (is.EQ.0) GOTO 290
!
      DO i = 1, klon
      IF(llflag(i)) THEN
         pqu(i,k) = pqu(i,k+1)
         ptu(i,k) = ptu(i,k+1)+(pgeoh(i,k+1)-pgeoh(i,k))/RCPD
         zbuo = ptu(i,k)*(1.+RETV*pqu(i,k))- &
                ptenh(i,k)*(1.+RETV*pqenh(i,k))+0.5
         IF (zbuo.GT.0.) klab(i,k) = 1
         zqold(i) = pqu(i,k)
      ENDIF
      ENDDO
!
      icall=1
      CALL flxadjtq(paph(1,k), ptu(1,k), pqu(1,k), llflag, icall)
!
      DO i = 1, klon
      IF (llflag(i).AND.pqu(i,k).NE.zqold(i)) THEN
         klab(i,k) = 2
         plu(i,k) = plu(i,k) + zqold(i)-pqu(i,k)
         zbuo = ptu(i,k)*(1.+RETV*pqu(i,k))- &
                ptenh(i,k)*(1.+RETV*pqenh(i,k))+0.5
         IF (zbuo.GT.0.) kcbot(i) = k
         IF (zbuo.GT.0.) ldcum(i) = .TRUE.
      ENDIF
      ENDDO
!
  290 CONTINUE
!
      RETURN
      END
