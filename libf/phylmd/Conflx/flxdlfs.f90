      SUBROUTINE flxdlfs(ptenh, pqenh, pgeoh, paph, ptu, pqu, &
           ldcum, kcbot, kctop, pmfub, prfl, ptd, pqd, &
           pmfd, pmfds, pmfdq, pdmfdp, kdtop, lddraf)
      use dimens_m
      use dimphy
      use SUPHEC_M
      use yoethf_m
            use yoecumf
      IMPLICIT none
!
!----------------------------------------------------------------------
! THIS ROUTINE CALCULATES LEVEL OF FREE SINKING FOR
! CUMULUS DOWNDRAFTS AND SPECIFIES T,Q,U AND V VALUES
!
! TO PRODUCE LFS-VALUES FOR CUMULUS DOWNDRAFTS
! FOR MASSFLUX CUMULUS PARAMETERIZATION
!
! INPUT ARE ENVIRONMENTAL VALUES OF T,Q,U,V,P,PHI
! AND UPDRAFT VALUES T,Q,U AND V AND ALSO
! CLOUD BASE MASSFLUX AND CU-PRECIPITATION RATE.
! IT RETURNS T,Q,U AND V VALUES AND MASSFLUX AT LFS.
!
! CHECK FOR NEGATIVE BUOYANCY OF AIR OF EQUAL PARTS OF
! MOIST ENVIRONMENTAL AIR AND CLOUD AIR.
!----------------------------------------------------------------------
!
      REAL ptenh(klon,klev)
      REAL pqenh(klon,klev)
      REAL pgeoh(klon,klev), paph(klon,klev+1)
      REAL ptu(klon,klev), pqu(klon,klev)
      REAL pmfub(klon)
      REAL prfl(klon)
!
      REAL ptd(klon,klev), pqd(klon,klev)
      REAL pmfd(klon,klev), pmfds(klon,klev), pmfdq(klon,klev)
      REAL pdmfdp(klon,klev)
      INTEGER  kcbot(klon), kctop(klon), kdtop(klon)
      LOGICAL  ldcum(klon), lddraf(klon)
!
      REAL ztenwb(klon,klev), zqenwb(klon,klev), zcond(klon)
      REAL zttest, zqtest, zbuo, zmftop
      LOGICAL  llo2(klon)
      INTEGER i, k, is, icall
!----------------------------------------------------------------------
      DO i= 1, klon
         lddraf(i)=.FALSE.
         kdtop(i)=klev+1
      ENDDO
!
!----------------------------------------------------------------------
! DETERMINE LEVEL OF FREE SINKING BY
! DOING A SCAN FROM TOP TO BASE OF CUMULUS CLOUDS
!
! FOR EVERY POINT AND PROCEED AS FOLLOWS:
!     (1) DETEMINE WET BULB ENVIRONMENTAL T AND Q
!     (2) DO MIXING WITH CUMULUS CLOUD AIR
!     (3) CHECK FOR NEGATIVE BUOYANCY
!
! THE ASSUMPTION IS THAT AIR OF DOWNDRAFTS IS MIXTURE
! OF 50% CLOUD AIR + 50% ENVIRONMENTAL AIR AT WET BULB
! TEMPERATURE (I.E. WHICH BECAME SATURATED DUE TO
! EVAPORATION OF RAIN AND CLOUD WATER)
!----------------------------------------------------------------------
!
      DO 290 k = 3, klev-3
!
      is=0
      DO 212 i= 1, klon
         ztenwb(i,k)=ptenh(i,k)
         zqenwb(i,k)=pqenh(i,k)
         llo2(i) = ldcum(i).AND.prfl(i).GT.0. &
                   .AND..NOT.lddraf(i) &
                   .AND.(k.LT.kcbot(i).AND.k.GT.kctop(i))
         IF ( llo2(i) ) is = is + 1
  212 CONTINUE
      IF(is.EQ.0) GO TO 290
!
      icall=2
      CALL flxadjtq(paph(1,k), ztenwb(1,k), zqenwb(1,k), llo2, icall)
!
!----------------------------------------------------------------------
! DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
! AND CHECK FOR NEGATIVE BUOYANCY.
! THEN SET VALUES FOR DOWNDRAFT AT LFS.
!----------------------------------------------------------------------
      DO 222 i= 1, klon
      IF (llo2(i)) THEN
         zttest=0.5*(ptu(i,k)+ztenwb(i,k))
         zqtest=0.5*(pqu(i,k)+zqenwb(i,k))
         zbuo=zttest*(1.+RETV*zqtest)- &
              ptenh(i,k)*(1.+RETV  *pqenh(i,k))
         zcond(i)=pqenh(i,k)-zqenwb(i,k)
         zmftop=-CMFDEPS*pmfub(i)
         IF (zbuo.LT.0..AND.prfl(i).GT.10.*zmftop*zcond(i)) THEN
            kdtop(i)=k
            lddraf(i)=.TRUE.
            ptd(i,k)=zttest
            pqd(i,k)=zqtest
            pmfd(i,k)=zmftop
            pmfds(i,k)=pmfd(i,k)*(RCPD*ptd(i,k)+pgeoh(i,k))
            pmfdq(i,k)=pmfd(i,k)*pqd(i,k)
            pdmfdp(i,k-1)=-0.5*pmfd(i,k)*zcond(i)
            prfl(i)=prfl(i)+pdmfdp(i,k-1)
         ENDIF
      ENDIF
  222 CONTINUE
!
  290 CONTINUE
!
      RETURN
      END
