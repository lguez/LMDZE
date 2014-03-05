module flxadjtq_m

  IMPLICIT none

contains

  SUBROUTINE flxadjtq(pp, pt, pq, ldflag, kcall)

    ! Objet : ajustement entre T et Q

    USE dimphy, ONLY: klon
    USE fcttre, ONLY: foede, foeew
    USE suphec_m, ONLY: rcpd, retv, rlstt, rlvtt, rtt
    USE yoethf_m, ONLY: r2es, r5ies, r5les

    REAL, intent(in):: pp(klon)
    real pt(klon), pq(klon)
    LOGICAL, intent(in):: ldflag(klon)
    INTEGER, intent(in):: kcall
    ! Defines calculation as:
    ! kcall = 0 env. T AND QS IN*CUINI*
    ! kcall = 1 condensation in updrafts (e.g. cubase, cuasc)
    ! kcall = 2 evaporation in downdrafts (e.g. cudlfs, cuddraf)

    ! Local:
    REAL zcond(klon), zcond1
    REAL Z5alvcp, z5alscp, zalvdcp, zalsdcp
    REAL zdelta, zcvm5, zldcp, zqsat, zcor
    INTEGER is, i

    !---------------------------------------------------------------------

    z5alvcp = r5les*RLVTT/RCPD
    z5alscp = r5ies*RLSTT/RCPD
    zalvdcp = rlvtt/RCPD
    zalsdcp = rlstt/RCPD

    DO i = 1, klon
       zcond(i) = 0.0
    ENDDO

    DO i = 1, klon
       IF (ldflag(i)) THEN
          zdelta = MAX(0., SIGN(1., RTT-pt(i)))
          zcvm5 = z5alvcp*(1.-zdelta) + zdelta*z5alscp
          zldcp = zalvdcp*(1.-zdelta) + zdelta*zalsdcp
          zqsat = R2ES * FOEEW(pt(i), zdelta) / pp(i)
          zqsat = MIN(0.5, zqsat)
          zcor = 1./(1.-RETV*zqsat)
          zqsat = zqsat*zcor
          zcond(i) = (pq(i)-zqsat) &
               / (1. + FOEDE(pt(i), zdelta, zcvm5, zqsat, zcor))
          IF (kcall.EQ.1) zcond(i) = MAX(zcond(i), 0.)
          IF (kcall.EQ.2) zcond(i) = MIN(zcond(i), 0.)
          pt(i) = pt(i) + zldcp*zcond(i)
          pq(i) = pq(i) - zcond(i)
       ENDIF
    end DO

    is = 0
    DO i = 1, klon
       IF (zcond(i).NE.0.) is = is + 1
    ENDDO
    IF (is /= 0) then
       DO i = 1, klon
          IF(ldflag(i).AND.zcond(i).NE.0.) THEN
             zdelta = MAX(0., SIGN(1., RTT-pt(i)))
             zcvm5 = z5alvcp*(1.-zdelta) + zdelta*z5alscp
             zldcp = zalvdcp*(1.-zdelta) + zdelta*zalsdcp
             zqsat = R2ES* FOEEW(pt(i), zdelta) / pp(i)
             zqsat = MIN(0.5, zqsat)
             zcor = 1./(1.-RETV*zqsat)
             zqsat = zqsat*zcor
             zcond1 = (pq(i)-zqsat) &
                  / (1. + FOEDE(pt(i), zdelta, zcvm5, zqsat, zcor))
             pt(i) = pt(i) + zldcp*zcond1
             pq(i) = pq(i) - zcond1
          ENDIF
       end DO
    end IF

  END SUBROUTINE flxadjtq

end module flxadjtq_m
