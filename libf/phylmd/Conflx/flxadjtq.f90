      SUBROUTINE flxadjtq(pp, pt, pq, ldflag, kcall)
      use dimens_m
      use dimphy
      use SUPHEC_M
      use yoethf_m
      use fcttre
      IMPLICIT none
!======================================================================
! Objet: ajustement entre T et Q
!======================================================================
! NOTE: INPUT PARAMETER kcall DEFINES CALCULATION AS
!        kcall=0    ENV. T AND QS IN*CUINI*
!        kcall=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
!        kcall=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)
!
!
      REAL pt(klon), pq(klon), pp(klon)
      LOGICAL ldflag(klon)
      INTEGER kcall
!
      REAL zcond(klon), zcond1
      REAL Z5alvcp, z5alscp, zalvdcp, zalsdcp
      REAL zdelta, zcvm5, zldcp, zqsat, zcor
      INTEGER is, i
!
      z5alvcp = r5les*RLVTT/RCPD
      z5alscp = r5ies*RLSTT/RCPD
      zalvdcp = rlvtt/RCPD
      zalsdcp = rlstt/RCPD
!

      DO i = 1, klon
         zcond(i) = 0.0
      ENDDO

      DO 210 i =1, klon
      IF (ldflag(i)) THEN
         zdelta = MAX(0.,SIGN(1.,RTT-pt(i)))
         zcvm5 = z5alvcp*(1.-zdelta) + zdelta*z5alscp
         zldcp = zalvdcp*(1.-zdelta) + zdelta*zalsdcp
         zqsat = R2ES*FOEEW(pt(i),zdelta) / pp(i)
         zqsat = MIN(0.5,zqsat)
         zcor = 1./(1.-RETV*zqsat)
         zqsat = zqsat*zcor
         zcond(i) = (pq(i)-zqsat) &
           / (1. + FOEDE(pt(i), zdelta, zcvm5, zqsat, zcor))
         IF (kcall.EQ.1) zcond(i) = MAX(zcond(i),0.)
         IF (kcall.EQ.2) zcond(i) = MIN(zcond(i),0.)
         pt(i) = pt(i) + zldcp*zcond(i)
         pq(i) = pq(i) - zcond(i)
      ENDIF
  210 CONTINUE
!
      is = 0
      DO i =1, klon
         IF (zcond(i).NE.0.) is = is + 1
      ENDDO
      IF (is.EQ.0) GOTO 230
!
      DO 220 i = 1, klon
      IF(ldflag(i).AND.zcond(i).NE.0.) THEN
         zdelta = MAX(0.,SIGN(1.,RTT-pt(i)))
         zcvm5 = z5alvcp*(1.-zdelta) + zdelta*z5alscp
         zldcp = zalvdcp*(1.-zdelta) + zdelta*zalsdcp
         zqsat = R2ES* FOEEW(pt(i),zdelta) / pp(i)
         zqsat = MIN(0.5,zqsat)
         zcor = 1./(1.-RETV*zqsat)
         zqsat = zqsat*zcor
         zcond1 = (pq(i)-zqsat) &
           / (1. + FOEDE(pt(i),zdelta,zcvm5,zqsat,zcor))
         pt(i) = pt(i) + zldcp*zcond1
         pq(i) = pq(i) - zcond1
      ENDIF
  220 CONTINUE
!
  230 CONTINUE
      RETURN
      END
