module diagcld2_m

  IMPLICIT none

contains

  SUBROUTINE diagcld2(paprs, pplay, t, q, diafra, dialiq)

    USE dimphy, ONLY : klev, klon
    USE suphec_m, ONLY : rcpd, rd, retv, rtt
    USE yoethf_m, ONLY : r2es
    USE fcttre, ONLY : foeew, qsatl, qsats

    ! Arguments d'entree:
    REAL, intent(in):: paprs(klon, klev+1) ! pression (Pa) a inter-couche
    REAL, intent(in):: pplay(klon, klev) ! pression (Pa) au milieu de couche
    REAL, intent(in):: t(klon, klev) ! temperature (K)
    REAL q(klon, klev) ! humidite specifique (Kg/Kg)

    ! Arguments de sortie:
    REAL diafra(klon, klev) ! fraction nuageuse diagnostiquee
    REAL dialiq(klon, klev) ! eau liquide nuageuse

    REAL, PARAMETER:: CETAMB = 0.8
    REAL CLOIA, CLOIB, CLOIC, CLOID
    PARAMETER (CLOIA=1.0E+02, CLOIB=-10.00, CLOIC=-0.6, CLOID=5.0)
    REAL RGAMMAS
    PARAMETER (RGAMMAS=0.05)
    REAL CRHL
    PARAMETER (CRHL=0.15)

    ! Variables locales:
    INTEGER i, k, kb, invb(klon)
    REAL zqs, zrhb, zcll, zdthmin(klon), zdthdp
    REAL zcor

    !-----------------------------------------------------------

    ! Initialisation:

    DO k = 1, klev
       DO i = 1, klon
          diafra(i, k) = 0.0
          dialiq(i, k) = 0.0
       ENDDO
    ENDDO

    DO i = 1, klon
       invb(i) = klev
       zdthmin(i)=0.0
    ENDDO

    DO k = 2, klev / 2 - 1
       DO i = 1, klon
          zdthdp = (t(i, k) - t(i, k+1)) / (pplay(i, k) - pplay(i, k+1)) &
               - RD * 0.5 * (t(i, k) + t(i, k+1)) / RCPD / paprs(i, k+1)
          zdthdp = zdthdp * CLOIA
          IF (pplay(i, k) > CETAMB * paprs(i, 1) .AND. zdthdp < zdthmin(i)) THEN
             zdthmin(i) = zdthdp
             invb(i) = k
          ENDIF
       ENDDO
    ENDDO

    DO i = 1, klon
       kb=invb(i)
       zqs= R2ES*FOEEW(t(i, kb), RTT >= t(i, kb))/pplay(i, kb)
       zqs=MIN(0.5, zqs)
       zcor=1./(1.-RETV*zqs)
       zqs=zqs*zcor
       zcll = CLOIB * zdthmin(i) + CLOIC
       zcll = MIN(1.0, MAX(0.0, zcll))
       zrhb= q(i, kb)/zqs
       IF (zcll > 0.0.AND.zrhb < CRHL) &
            zcll=zcll*(1.-(CRHL-zrhb)*CLOID)
       zcll=MIN(1.0, MAX(0.0, zcll))
       diafra(i, kb) = MAX(diafra(i, kb), zcll)
       dialiq(i, kb)= diafra(i, kb) * RGAMMAS*zqs
    ENDDO

  END SUBROUTINE diagcld2

end module diagcld2_m
