      SUBROUTINE flxddraf(ptenh, pqenh, pgeoh, paph, prfl, &
                 ptd, pqd, pmfd, pmfds, pmfdq, pdmfdp, &
                 lddraf, pen_d, pde_d)
      use dimens_m
      use dimphy
      use flxadjtq_m, only: flxadjtq
      use SUPHEC_M
      use yoethf_m
            use yoecumf
      IMPLICIT none
!
!----------------------------------------------------------------------
!          THIS ROUTINE CALCULATES CUMULUS DOWNDRAFT DESCENT
!
!          TO PRODUCE THE VERTICAL PROFILES FOR CUMULUS DOWNDRAFTS
!          (I.E. T,Q,U AND V AND FLUXES)
!
!          INPUT IS T,Q,P,PHI,U,V AT HALF LEVELS.
!          IT RETURNS FLUXES OF S,Q AND EVAPORATION RATE
!          AND U,V AT LEVELS WHERE DOWNDRAFT OCCURS
!
!          CALCULATE MOIST DESCENT FOR ENTRAINING/DETRAINING PLUME BY
!          A) MOVING AIR DRY-ADIABATICALLY TO NEXT LEVEL BELOW AND
!          B) CORRECTING FOR EVAPORATION TO OBTAIN SATURATED STATE.
!
!----------------------------------------------------------------------
!
      REAL ptenh(klon,klev), pqenh(klon,klev)
      REAL pgeoh(klon,klev), paph(klon,klev+1)
!
      REAL ptd(klon,klev), pqd(klon,klev)
      REAL pmfd(klon,klev), pmfds(klon,klev), pmfdq(klon,klev)
      REAL pdmfdp(klon,klev)
      REAL prfl(klon)
      LOGICAL lddraf(klon)
!
      REAL pen_d(klon,klev), pde_d(klon,klev), zcond(klon)
      LOGICAL llo2(klon), llo1
      INTEGER i, k, is, icall, itopde
      REAL zentr, zseen, zqeen, zsdde, zqdde, zmfdsk, zmfdqk, zdmfdp
      REAL zbuo
!----------------------------------------------------------------------
! CALCULATE MOIST DESCENT FOR CUMULUS DOWNDRAFT BY
!       (A) CALCULATING ENTRAINMENT RATES, ASSUMING
!           LINEAR DECREASE OF MASSFLUX IN PBL
!       (B) DOING MOIST DESCENT - EVAPORATIVE COOLING
!           AND MOISTENING IS CALCULATED IN *flxadjtq*
!       (C) CHECKING FOR NEGATIVE BUOYANCY AND
!           SPECIFYING FINAL T,Q,U,V AND DOWNWARD FLUXES
!
      DO 180 k = 3, klev
!
      is = 0
      DO i = 1, klon
         llo2(i)=lddraf(i).AND.pmfd(i,k-1).LT.0.
         IF (llo2(i)) is = is + 1
      ENDDO
      IF (is.EQ.0) GOTO 180
!
      DO i = 1, klon
      IF (llo2(i)) THEN
         zentr = ENTRDD*pmfd(i,k-1)*RD*ptenh(i,k-1)/ &
                 (RG*paph(i,k-1))*(paph(i,k)-paph(i,k-1))
         pen_d(i,k) = zentr
         pde_d(i,k) = zentr
      ENDIF
      ENDDO
!
      itopde = klev-2
      IF (k.GT.itopde) THEN
         DO i = 1, klon
         IF (llo2(i)) THEN
            pen_d(i,k)=0.
            pde_d(i,k)=pmfd(i,itopde)* &
            (paph(i,k)-paph(i,k-1))/(paph(i,klev+1)-paph(i,itopde))
         ENDIF
         ENDDO
      ENDIF
!
      DO i = 1, klon
      IF (llo2(i)) THEN
         pmfd(i,k) = pmfd(i,k-1)+pen_d(i,k)-pde_d(i,k)
         zseen = (RCPD*ptenh(i,k-1)+pgeoh(i,k-1))*pen_d(i,k)
         zqeen = pqenh(i,k-1)*pen_d(i,k)
         zsdde = (RCPD*ptd(i,k-1)+pgeoh(i,k-1))*pde_d(i,k)
         zqdde = pqd(i,k-1)*pde_d(i,k)
         zmfdsk = pmfds(i,k-1)+zseen-zsdde
         zmfdqk = pmfdq(i,k-1)+zqeen-zqdde
         pqd(i,k) = zmfdqk*(1./MIN(-CMFCMIN,pmfd(i,k)))
         ptd(i,k) = (zmfdsk*(1./MIN(-CMFCMIN,pmfd(i,k)))- &
                     pgeoh(i,k))/RCPD
         ptd(i,k) = MIN(400.,ptd(i,k))
         ptd(i,k) = MAX(100.,ptd(i,k))
         zcond(i) = pqd(i,k)
      ENDIF
      ENDDO
!
      icall = 2
      CALL flxadjtq(paph(1,k), ptd(1,k), pqd(1,k), llo2, icall)
!
      DO i = 1, klon
      IF (llo2(i)) THEN
         zcond(i) = zcond(i)-pqd(i,k)
         zbuo = ptd(i,k)*(1.+RETV  *pqd(i,k))- &
                ptenh(i,k)*(1.+RETV  *pqenh(i,k))
         llo1 = zbuo.LT.0..AND.(prfl(i)-pmfd(i,k)*zcond(i).GT.0.)
         IF (.not.llo1) pmfd(i,k) = 0.0
         pmfds(i,k) = (RCPD*ptd(i,k)+pgeoh(i,k))*pmfd(i,k)
         pmfdq(i,k) = pqd(i,k)*pmfd(i,k)
         zdmfdp = -pmfd(i,k)*zcond(i)
         pdmfdp(i,k-1) = zdmfdp
         prfl(i) = prfl(i)+zdmfdp
      ENDIF
      ENDDO
!
  180 CONTINUE
      RETURN
      END
