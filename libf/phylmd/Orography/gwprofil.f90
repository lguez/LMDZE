    SUBROUTINE gwprofil(nlon,nlev,kgwd,kdx,ktest,kkcrith,kcrit,paphm1,prho, &
        pstab,pvph,pri,ptau,pdmod,psig,pvar)

!**** *GWPROFIL*

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!          FROM *GWDRAG*

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
!     ==== OUTPUTS ===

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD:
!     -------
!     THE STRESS PROFILE FOR GRAVITY WAVES IS COMPUTED AS FOLLOWS:
!     IT IS CONSTANT (NO GWD) AT THE LEVELS BETWEEN THE GROUND
!     AND THE TOP OF THE BLOCKED LAYER (KKENVH).
!     IT DECREASES LINEARLY WITH HEIGHTS FROM THE TOP OF THE
!     BLOCKED LAYER TO 3*VAROR (kKNU), TO SIMULATES LEE WAVES OR
!     NONLINEAR GRAVITY WAVE BREAKING.
!     ABOVE IT IS CONSTANT, EXCEPT WHEN THE WAVE ENCOUNTERS A CRITICAL
!     LEVEL (KCRIT) OR WHEN IT BREAKS.



!     EXTERNALS.
!     ----------


!     REFERENCE.
!     ----------

!        SEE ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S."

!     AUTHOR.
!     -------

!     MODIFICATIONS.
!     --------------
!     PASSAGE OF THE NEW GWDRAG TO I.F.S. (F. LOTT, 22/11/93)
!-----------------------------------------------------------------------
      USE dimens_m
      USE dimphy
      USE suphec_m
      USE yoegwd
      IMPLICIT NONE





!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

      INTEGER nlon, nlev
      INTEGER kkcrith(nlon), kcrit(nlon), kdx(nlon), ktest(nlon)


      REAL paphm1(nlon,nlev+1), pstab(nlon,nlev+1), prho(nlon,nlev+1), &
        pvph(nlon,nlev+1), pri(nlon,nlev+1), ptau(nlon,nlev+1)

      REAL pdmod(nlon)
      REAL, INTENT (IN) :: psig(nlon)
      REAL, INTENT (IN) :: pvar(nlon)

!-----------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

      INTEGER ilevh, ji, kgwd, jl, jk
      REAL zsqr, zalfa, zriw, zdel, zb, zalpha, zdz2n
      REAL zdelp, zdelpt
      REAL zdz2(klon,klev), znorm(klon), zoro(klon)
      REAL ztau(klon,klev+1)

!-----------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------

!      print *,' entree gwprofil'
100   CONTINUE


!*    COMPUTATIONAL CONSTANTS.
!     ------------- ----------

      ilevh = klev/3

!     DO 400 ji=1,kgwd
!     jl=kdx(ji)
!  Modif vectorisation 02/04/2004
      DO 400 jl = 1, klon
        IF (ktest(jl)==1) THEN
          zoro(jl) = psig(jl)*pdmod(jl)/4./max(pvar(jl),1.0)
          ztau(jl,klev+1) = ptau(jl,klev+1)
        END IF
400   CONTINUE


      DO 430 jk = klev, 2, -1


!*         4.1    CONSTANT WAVE STRESS UNTIL TOP OF THE
!                 BLOCKING LAYER.
410     CONTINUE

!     DO 411 ji=1,kgwd
!     jl=kdx(ji)
!  Modif vectorisation 02/04/2004
        DO 411 jl = 1, klon
          IF (ktest(jl)==1) THEN
            IF (jk>kkcrith(jl)) THEN
              ptau(jl,jk) = ztau(jl,klev+1)
!          ENDIF
!          IF(JK.EQ.KKCRITH(JL)) THEN
            ELSE
              ptau(jl,jk) = grahilo*ztau(jl,klev+1)
            END IF
          END IF
411     CONTINUE

!*         4.15   CONSTANT SHEAR STRESS UNTIL THE TOP OF THE
!                 LOW LEVEL FLOW LAYER.
415     CONTINUE


!*         4.2    WAVE DISPLACEMENT AT NEXT LEVEL.

420     CONTINUE

!     DO 421 ji=1,kgwd
!     jl=kdx(ji)
!  Modif vectorisation 02/04/2004
        DO 421 jl = 1, klon
          IF (ktest(jl)==1) THEN
            IF (jk<kkcrith(jl)) THEN
              znorm(jl) = gkdrag*prho(jl,jk)*sqrt(pstab(jl,jk))*pvph(jl,jk)* &
                zoro(jl)
              zdz2(jl,jk) = ptau(jl,jk+1)/max(znorm(jl),gssec)
            END IF
          END IF
421     CONTINUE

!*         4.3    WAVE RICHARDSON NUMBER, NEW WAVE DISPLACEMENT
!*                AND STRESS:  BREAKING EVALUATION AND CRITICAL
!                 LEVEL


!     DO 431 ji=1,kgwd
!     jl=Kdx(ji)
!  Modif vectorisation 02/04/2004
        DO 431 jl = 1, klon
          IF (ktest(jl)==1) THEN

            IF (jk<kkcrith(jl)) THEN
              IF ((ptau(jl,jk+1)<gtsec) .OR. (jk<=kcrit(jl))) THEN
                ptau(jl,jk) = 0.0
              ELSE
                zsqr = sqrt(pri(jl,jk))
                zalfa = sqrt(pstab(jl,jk)*zdz2(jl,jk))/pvph(jl,jk)
                zriw = pri(jl,jk)*(1.-zalfa)/(1+zalfa*zsqr)**2
                IF (zriw<grcrit) THEN
                  zdel = 4./zsqr/grcrit + 1./grcrit**2 + 4./grcrit
                  zb = 1./grcrit + 2./zsqr
                  zalpha = 0.5*(-zb+sqrt(zdel))
                  zdz2n = (pvph(jl,jk)*zalpha)**2/pstab(jl,jk)
                  ptau(jl,jk) = znorm(jl)*zdz2n
                ELSE
                  ptau(jl,jk) = znorm(jl)*zdz2(jl,jk)
                END IF
                ptau(jl,jk) = min(ptau(jl,jk),ptau(jl,jk+1))
              END IF
            END IF
          END IF
431     CONTINUE

430   CONTINUE
440   CONTINUE

!  REORGANISATION OF THE STRESS PROFILE AT LOW LEVEL

!     DO 530 ji=1,kgwd
!     jl=kdx(ji)
!  Modif vectorisation 02/04/2004
      DO 530 jl = 1, klon
        IF (ktest(jl)==1) THEN
          ztau(jl,kkcrith(jl)) = ptau(jl,kkcrith(jl))
          ztau(jl,nstra) = ptau(jl,nstra)
        END IF
530   CONTINUE

      DO 531 jk = 1, klev

!     DO 532 ji=1,kgwd
!     jl=kdx(ji)
!  Modif vectorisation 02/04/2004
        DO 532 jl = 1, klon
          IF (ktest(jl)==1) THEN


            IF (jk>kkcrith(jl)) THEN

              zdelp = paphm1(jl,jk) - paphm1(jl,klev+1)
              zdelpt = paphm1(jl,kkcrith(jl)) - paphm1(jl,klev+1)
              ptau(jl,jk) = ztau(jl,klev+1) + (ztau(jl,kkcrith(jl))-ztau(jl, &
                klev+1))*zdelp/zdelpt

            END IF

          END IF
532     CONTINUE

!  REORGANISATION IN THE STRATOSPHERE

!     DO 533 ji=1,kgwd
!     jl=kdx(ji)
!  Modif vectorisation 02/04/2004
        DO 533 jl = 1, klon
          IF (ktest(jl)==1) THEN


            IF (jk<nstra) THEN

              zdelp = paphm1(jl,nstra)
              zdelpt = paphm1(jl,jk)
              ptau(jl,jk) = ztau(jl,nstra)*zdelpt/zdelp

            END IF

          END IF
533     CONTINUE

! REORGANISATION IN THE TROPOSPHERE

!      DO 534 ji=1,kgwd
!      jl=kdx(ji)
!  Modif vectorisation 02/04/2004
        DO 534 jl = 1, klon
          IF (ktest(jl)==1) THEN


            IF (jk<kkcrith(jl) .AND. jk>nstra) THEN

              zdelp = paphm1(jl,jk) - paphm1(jl,kkcrith(jl))
              zdelpt = paphm1(jl,nstra) - paphm1(jl,kkcrith(jl))
              ptau(jl,jk) = ztau(jl,kkcrith(jl)) + (ztau(jl,nstra)-ztau(jl, &
                kkcrith(jl)))*zdelp/zdelpt

            END IF
          END IF
534     CONTINUE


531   CONTINUE


      RETURN
    END
