module gwprofil_m

  IMPLICIT NONE

contains

  SUBROUTINE gwprofil(nlon, nlev, ktest, kkcrith, kcrit, paphm1, prho, pstab, &
       pvph, pri, ptau, pdmod, psig, pvar)

    ! Method. The stress profile for gravity waves is computed as
    ! follows: it is constant (no gwd) at the levels between the
    ! ground and the top of the blocked layer (kkenvh).  It decreases
    ! linearly with height from the top of the blocked layer to 3 *
    ! varor (kknu), to simulate lee waves or nonlinear gravity wave
    ! breaking. Above it is constant, except when the wave encounters
    ! a critical level (kcrit) or when it breaks.

    ! Reference. See ECMWF research department documentation of the
    ! "I.F.S."

    ! Modifications. Passage of the new gwdrag TO I.F.S. (F. LOTT,
    ! 22/11/93)

    USE dimphy, ONLY : klev, klon
    USE yoegwd, ONLY : gkdrag, grahilo, grcrit, gssec, gtsec, nstra

    INTEGER, intent(in):: nlon, nlev
    INTEGER, intent(in):: ktest(nlon), kkcrith(nlon), kcrit(nlon)
    REAL, intent(in):: paphm1(nlon, nlev+1), prho(nlon, nlev+1)
    REAL, intent(in):: pstab(nlon, nlev+1)
    real, intent(in):: pvph(nlon, nlev+1), pri(nlon, nlev+1)
    real ptau(nlon, nlev+1)
    REAL, intent(in):: pdmod(nlon)
    REAL, INTENT (IN) :: psig(nlon)
    REAL, INTENT (IN) :: pvar(nlon)

    ! Local:
    INTEGER jl, jk
    REAL zsqr, zalfa, zriw, zdel, zb, zalpha, zdz2n
    REAL zdelp, zdelpt
    REAL zdz2(klon, klev), znorm(klon), zoro(klon)
    REAL ztau(klon, klev+1)

    !-----------------------------------------------------------------------

    ! 1. INITIALIZATION

    ! COMPUTATIONAL CONSTANTS.

    DO jl = 1, klon
       IF (ktest(jl)==1) THEN
          zoro(jl) = psig(jl)*pdmod(jl)/4./max(pvar(jl), 1.0)
          ztau(jl, klev+1) = ptau(jl, klev+1)
       END IF
    end DO

    DO jk = klev, 2, -1
       ! 4.1 CONSTANT WAVE STRESS UNTIL TOP OF THE
       ! BLOCKING LAYER.
       DO jl = 1, klon
          IF (ktest(jl)==1) THEN
             IF (jk>kkcrith(jl)) THEN
                ptau(jl, jk) = ztau(jl, klev+1)
             ELSE
                ptau(jl, jk) = grahilo*ztau(jl, klev+1)
             END IF
          END IF
       end DO

       ! 4.2 WAVE DISPLACEMENT AT NEXT LEVEL.
       DO jl = 1, klon
          IF (ktest(jl)==1) THEN
             IF (jk<kkcrith(jl)) THEN
                znorm(jl) = gkdrag * prho(jl, jk) * sqrt(pstab(jl, jk)) &
                     * pvph(jl, jk)* zoro(jl)
                zdz2(jl, jk) = ptau(jl, jk+1)/max(znorm(jl), gssec)
             END IF
          END IF
       end DO

       ! 4.3 WAVE RICHARDSON NUMBER, NEW WAVE DISPLACEMENT
       ! AND STRESS: BREAKING EVALUATION AND CRITICAL
       ! LEVEL
       DO jl = 1, klon
          IF (ktest(jl)==1) THEN
             IF (jk<kkcrith(jl)) THEN
                IF ((ptau(jl, jk+1)<gtsec) .OR. (jk<=kcrit(jl))) THEN
                   ptau(jl, jk) = 0.0
                ELSE
                   zsqr = sqrt(pri(jl, jk))
                   zalfa = sqrt(pstab(jl, jk)*zdz2(jl, jk))/pvph(jl, jk)
                   zriw = pri(jl, jk)*(1.-zalfa)/(1+zalfa*zsqr)**2
                   IF (zriw<grcrit) THEN
                      zdel = 4./zsqr/grcrit + 1./grcrit**2 + 4./grcrit
                      zb = 1./grcrit + 2./zsqr
                      zalpha = 0.5*(-zb+sqrt(zdel))
                      zdz2n = (pvph(jl, jk)*zalpha)**2/pstab(jl, jk)
                      ptau(jl, jk) = znorm(jl)*zdz2n
                   ELSE
                      ptau(jl, jk) = znorm(jl)*zdz2(jl, jk)
                   END IF
                   ptau(jl, jk) = min(ptau(jl, jk), ptau(jl, jk+1))
                END IF
             END IF
          END IF
       end DO
    end DO

    ! REORGANISATION OF THE STRESS PROFILE AT LOW LEVEL

    DO jl = 1, klon
       IF (ktest(jl)==1) THEN
          ztau(jl, kkcrith(jl)) = ptau(jl, kkcrith(jl))
          ztau(jl, nstra) = ptau(jl, nstra)
       END IF
    end DO

    DO jk = 1, klev
       DO jl = 1, klon
          IF (ktest(jl)==1) THEN
             IF (jk>kkcrith(jl)) THEN
                zdelp = paphm1(jl, jk) - paphm1(jl, klev+1)
                zdelpt = paphm1(jl, kkcrith(jl)) - paphm1(jl, klev+1)
                ptau(jl, jk) = ztau(jl, klev+1) &
                     + (ztau(jl, kkcrith(jl)) - ztau(jl, klev+1))*zdelp/zdelpt
             END IF
          END IF
       end DO

       ! REORGANISATION IN THE STRATOSPHERE
       DO jl = 1, klon
          IF (ktest(jl)==1) THEN
             IF (jk < nstra) THEN
                zdelp = paphm1(jl, nstra)
                zdelpt = paphm1(jl, jk)
                ptau(jl, jk) = ztau(jl, nstra) * zdelpt / zdelp
             END IF
          END IF
       end DO

       ! REORGANISATION IN THE TROPOSPHERE
       DO jl = 1, klon
          IF (ktest(jl)==1) THEN
             IF (jk<kkcrith(jl) .AND. jk > nstra) THEN
                zdelp = paphm1(jl, jk) - paphm1(jl, kkcrith(jl))
                zdelpt = paphm1(jl, nstra) - paphm1(jl, kkcrith(jl))
                ptau(jl, jk) = ztau(jl, kkcrith(jl)) &
                     + (ztau(jl, nstra) - ztau(jl, kkcrith(jl)))*zdelp/zdelpt
             END IF
          END IF
       end DO
    end DO

  END SUBROUTINE gwprofil

end module gwprofil_m
