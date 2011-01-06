    SUBROUTINE orolift(nlon,nlev,ktest,ptsphy,paphm1,pgeom1,ptm1,pum1,pvm1, &
        plat,pmea,pvaror,ppic & 
        ,pulow,pvlow,pvom,pvol,pte)


!**** *OROLIFT: SIMULATE THE GEOSTROPHIC LIFT.

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!          CALLED FROM *lift_noro
!     ----------

!     AUTHOR.
!     -------
!     F.LOTT  LMD 22/11/95

      USE dimens_m
      USE dimphy
      USE suphec_m
      USE yoegwd
      IMPLICIT NONE


!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------


      INTEGER nlon, nlev
      REAL pte(nlon,nlev), pvol(nlon,nlev), pvom(nlon,nlev), pulow(nlon), &
        pvlow(nlon)
      REAL pum1(nlon,nlev), pvm1(nlon,nlev), ptm1(nlon,nlev)
      REAL, INTENT (IN) :: plat(nlon)
      REAL pmea(nlon)
      REAL, INTENT (IN) :: pvaror(nlon)
      REAL ppic(nlon), pgeom1(nlon,nlev), paphm1(nlon,nlev+1)

      INTEGER ktest(nlon)
      REAL, INTENT (IN) :: ptsphy
!-----------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------
      LOGICAL lifthigh
      INTEGER klevm1, jl, ilevh, jk
      REAL zcons1, ztmst, zrtmst, zpi, zhgeo
      REAL zdelp, zslow, zsqua, zscav, zbet
      INTEGER iknub(klon), iknul(klon)
      LOGICAL ll1(klon,klev+1)

      REAL ztau(klon,klev+1), ztav(klon,klev+1), zrho(klon,klev+1)
      REAL zdudt(klon), zdvdt(klon)
      REAL zhcrit(klon,klev)
!-----------------------------------------------------------------------

!*         1.1  INITIALIZATIONS
!               ---------------

      lifthigh = .FALSE.

      IF (nlon/=klon .OR. nlev/=klev) STOP
      zcons1 = 1./rd
      klevm1 = klev - 1
      ztmst = ptsphy
      zrtmst = 1./ztmst
      zpi = acos(-1.)

      DO 1001 jl = 1, klon
        zrho(jl,klev+1) = 0.0
        pulow(jl) = 0.0
        pvlow(jl) = 0.0
        iknub(jl) = klev
        iknul(jl) = klev
        ilevh = klev/3
        ll1(jl,klev+1) = .FALSE.
        DO 1000 jk = 1, klev
          pvom(jl,jk) = 0.0
          pvol(jl,jk) = 0.0
          pte(jl,jk) = 0.0
1000    CONTINUE
1001  CONTINUE


!*         2.1     DEFINE LOW LEVEL WIND, PROJECT WINDS IN PLANE OF
!*                 LOW LEVEL WIND, DETERMINE SECTOR IN WHICH TO TAKE
!*                 THE VARIANCE AND SET INDICATOR FOR CRITICAL LEVELS.



      DO 2006 jk = klev, 1, -1
        DO 2007 jl = 1, klon
          IF (ktest(jl)==1) THEN
            zhcrit(jl,jk) = amax1(ppic(jl)-pmea(jl),100.)
            zhgeo = pgeom1(jl,jk)/rg
            ll1(jl,jk) = (zhgeo>zhcrit(jl,jk))
            IF (ll1(jl,jk) .NEQV. ll1(jl,jk+1)) THEN
              iknub(jl) = jk
            END IF
          END IF
2007    CONTINUE
2006  CONTINUE

      DO 2010 jl = 1, klon
        IF (ktest(jl)==1) THEN
          iknub(jl) = max(iknub(jl),klev/2)
          iknul(jl) = max(iknul(jl),2*klev/3)
          IF (iknub(jl)>nktopg) iknub(jl) = nktopg
          IF (iknub(jl)==nktopg) iknul(jl) = klev
          IF (iknub(jl)==iknul(jl)) iknub(jl) = iknul(jl) - 1
        END IF
2010  CONTINUE

!     do 2011 jl=1,klon
!     IF(KTEST(JL).EQ.1) THEN
!       print *,' iknul= ',iknul(jl),'  iknub=',iknub(jl)
!     ENDIF
!2011 continue

!     PRINT *,'  DANS OROLIFT: 2010'

      DO 223 jk = klev, 2, -1
        DO 222 jl = 1, klon
          zrho(jl,jk) = 2.*paphm1(jl,jk)*zcons1/(ptm1(jl,jk)+ptm1(jl,jk-1))
222     CONTINUE
223   CONTINUE
!     PRINT *,'  DANS OROLIFT: 223'

!********************************************************************

!*     DEFINE LOW LEVEL FLOW
!      -------------------
      DO 2115 jk = klev, 1, -1
        DO 2116 jl = 1, klon
          IF (ktest(jl)==1) THEN
            IF (jk>=iknub(jl) .AND. jk<=iknul(jl)) THEN
              pulow(jl) = pulow(jl) + pum1(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl, &
                jk))
              pvlow(jl) = pvlow(jl) + pvm1(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl, &
                jk))
              zrho(jl,klev+1) = zrho(jl,klev+1) + zrho(jl,jk)*(paphm1(jl,jk+1) &
                -paphm1(jl,jk))
            END IF
          END IF
2116    CONTINUE
2115  CONTINUE
      DO 2110 jl = 1, klon
        IF (ktest(jl)==1) THEN
          pulow(jl) = pulow(jl)/(paphm1(jl,iknul(jl)+1)-paphm1(jl,iknub(jl)))
          pvlow(jl) = pvlow(jl)/(paphm1(jl,iknul(jl)+1)-paphm1(jl,iknub(jl)))
          zrho(jl,klev+1) = zrho(jl,klev+1)/(paphm1(jl,iknul(jl)+1)-paphm1(jl, &
            iknub(jl)))
        END IF
2110  CONTINUE


200   CONTINUE

!***********************************************************

!*         3.      COMPUTE MOUNTAIN LIFT

300   CONTINUE

      DO 301 jl = 1, klon
        IF (ktest(jl)==1) THEN
          ztau(jl,klev+1) = -gklift*zrho(jl,klev+1)*2.*romega*2*pvaror(jl)*sin &
            (zpi/180.*plat(jl))*pvlow(jl)
          ztav(jl,klev+1) = gklift*zrho(jl,klev+1)*2.*romega*2*pvaror(jl)* &
            sin(zpi/180.*plat(jl))*pulow(jl)
        ELSE
          ztau(jl,klev+1) = 0.0
          ztav(jl,klev+1) = 0.0
        END IF
301   CONTINUE


!*         4.      COMPUTE LIFT PROFILE
!*                 --------------------


400   CONTINUE

      DO 401 jk = 1, klev
        DO 401 jl = 1, klon
          IF (ktest(jl)==1) THEN
            ztau(jl,jk) = ztau(jl,klev+1)*paphm1(jl,jk)/paphm1(jl,klev+1)
            ztav(jl,jk) = ztav(jl,klev+1)*paphm1(jl,jk)/paphm1(jl,klev+1)
          ELSE
            ztau(jl,jk) = 0.0
            ztav(jl,jk) = 0.0
          END IF
401   CONTINUE


!*         5.      COMPUTE TENDENCIES.
!*                 -------------------
      IF (lifthigh) THEN

500     CONTINUE
!     PRINT *,'  DANS OROLIFT: 500'

!  EXPLICIT SOLUTION AT ALL LEVELS

        DO 524 jk = 1, klev
          DO 523 jl = 1, klon
            IF (ktest(jl)==1) THEN
              zdelp = paphm1(jl,jk+1) - paphm1(jl,jk)
              zdudt(jl) = -rg*(ztau(jl,jk+1)-ztau(jl,jk))/zdelp
              zdvdt(jl) = -rg*(ztav(jl,jk+1)-ztav(jl,jk))/zdelp
            END IF
523       CONTINUE
524     CONTINUE

!  PROJECT PERPENDICULARLY TO U NOT TO DESTROY ENERGY

        DO 530 jk = 1, klev
          DO 530 jl = 1, klon
            IF (ktest(jl)==1) THEN

              zslow = sqrt(pulow(jl)**2+pvlow(jl)**2)
              zsqua = amax1(sqrt(pum1(jl,jk)**2+pvm1(jl,jk)**2),gvsec)
              zscav = -zdudt(jl)*pvm1(jl,jk) + zdvdt(jl)*pum1(jl,jk)
              IF (zsqua>gvsec) THEN
                pvom(jl,jk) = -zscav*pvm1(jl,jk)/zsqua**2
                pvol(jl,jk) = zscav*pum1(jl,jk)/zsqua**2
              ELSE
                pvom(jl,jk) = 0.0
                pvol(jl,jk) = 0.0
              END IF
              zsqua = sqrt(pum1(jl,jk)**2+pum1(jl,jk)**2)
              IF (zsqua<zslow) THEN
                pvom(jl,jk) = zsqua/zslow*pvom(jl,jk)
                pvol(jl,jk) = zsqua/zslow*pvol(jl,jk)
              END IF

            END IF
530     CONTINUE

!  6.  LOW LEVEL LIFT, SEMI IMPLICIT:
!  ----------------------------------

      ELSE

        DO 601 jl = 1, klon
          IF (ktest(jl)==1) THEN
            DO jk = klev, iknub(jl), -1
              zbet = gklift*2.*romega*sin(zpi/180.*plat(jl))*ztmst* &
                (pgeom1(jl,iknub(jl)-1)-pgeom1(jl,jk))/ &
                (pgeom1(jl,iknub(jl)-1)-pgeom1(jl,klev))
              zdudt(jl) = -pum1(jl,jk)/ztmst/(1+zbet**2)
              zdvdt(jl) = -pvm1(jl,jk)/ztmst/(1+zbet**2)
              pvom(jl,jk) = zbet**2*zdudt(jl) - zbet*zdvdt(jl)
              pvol(jl,jk) = zbet*zdudt(jl) + zbet**2*zdvdt(jl)
            END DO
          END IF
601     CONTINUE

      END IF

      RETURN
    END
