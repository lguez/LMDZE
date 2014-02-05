module flxasc_m

  IMPLICIT none

contains

  SUBROUTINE flxasc(pdtime, ptenh, pqenh, pten, pqen, pqsen, pgeo, pgeoh, &
       pap, paph, pqte, pvervel, ldland, ldcum, ktype, klab, ptu, pqu, plu, &
       pmfu, pmfub, pentr, pmfus, pmfuq, pmful, plude, pdmfup, kcbot, kctop, &
       kctop0, kcum, pen_u, pde_u)

    ! This routine does the calculations for cloud ascents for cumulus
    ! parameterization.

    USE dimphy, ONLY: klev, klon
    use flxadjtq_m, only: flxadjtq
    USE suphec_m, ONLY: rcpd, rd, retv, rg, rtt
    USE yoecumf, ONLY: cmfcmin, cmfctop, cprcon, entrmid, lmfmid

    REAL, intent(in):: pdtime
    REAL, intent(in):: ptenh(klon, klev)
    REAL, intent(in):: pqenh(klon, klev)
    REAL, intent(in):: pten(klon, klev)
    REAL, intent(in):: pqen(klon, klev)
    REAL, intent(in):: pqsen(klon, klev)
    REAL, intent(in):: pgeo(klon, klev), pgeoh(klon, klev)
    REAL, intent(in):: pap(klon, klev), paph(klon, klev+1)
    REAL, intent(in):: pqte(klon, klev)
    REAL, intent(in):: pvervel(klon, klev) ! vitesse verticale en Pa/s
    LOGICAL, intent(in):: ldland(klon)
    LOGICAL, intent(inout):: ldcum(klon)
    INTEGER, intent(inout):: ktype(klon)
    integer klab(klon, klev)
    REAL ptu(klon, klev), pqu(klon, klev), plu(klon, klev)
    REAL pmfu(klon, klev)
    REAL, intent(inout):: pmfub(klon)
    real pentr(klon)
    real pmfus(klon, klev)
    REAL pmfuq(klon, klev), pmful(klon, klev)
    REAL plude(klon, klev)
    REAL pdmfup(klon, klev)
    integer kcbot(klon), kctop(klon)
    INTEGER kctop0(klon)
    integer, intent(out):: kcum
    REAL pen_u(klon, klev), pde_u(klon, klev)

    ! Local:

    REAL zqold(klon)
    REAL zdland(klon)
    LOGICAL llflag(klon)
    INTEGER k, i, is, icall
    REAL ztglace, zdphi, zqeen, zseen, zscde, zqude
    REAL zmfusk, zmfuqk, zmfulk, zbuo, zdnoprc, zprcon, zlnew

    REAL zpbot(klon), zptop(klon), zrho(klon)
    REAL zdprho, zentr, zpmid, zmftest, zmfmax
    LOGICAL llo1, llo2

    REAL zwmax(klon), zzzmb
    INTEGER klwmin(klon) ! level of maximum vertical velocity
    real fact

    !----------------------------------------------------------------------

    ztglace = RTT - 13.

    ! Chercher le niveau où la vitesse verticale est maximale :

    DO i = 1, klon
       klwmin(i) = klev
       zwmax(i) = 0.0
    ENDDO

    DO k = klev, 3, -1
       DO i = 1, klon
          IF (pvervel(i, k) < zwmax(i)) THEN
             zwmax(i) = pvervel(i, k)
             klwmin(i) = k
          ENDIF
       ENDDO
    ENDDO

    ! Set default values:

    DO i = 1, klon
       IF (.NOT. ldcum(i)) ktype(i)=0
    ENDDO

    DO k=1, klev
       DO i = 1, klon
          plu(i, k)=0.
          pmfu(i, k)=0.
          pmfus(i, k)=0.
          pmfuq(i, k)=0.
          pmful(i, k)=0.
          plude(i, k)=0.
          pdmfup(i, k)=0.
          IF (.NOT. ldcum(i) .OR. ktype(i) == 3) klab(i, k)=0
          IF (.NOT. ldcum(i) .AND. paph(i, k) < 4e4) kctop0(i) = k
       ENDDO
    ENDDO

    DO i = 1, klon
       IF (ldland(i)) THEN
          zdland(i)=3.0E4
          zdphi=pgeoh(i, kctop0(i))-pgeoh(i, kcbot(i))
          IF (ptu(i, kctop0(i)) >= ztglace) zdland(i)=zdphi
          zdland(i)=MAX(3.0E4, zdland(i))
          zdland(i)=MIN(5.0E4, zdland(i))
       ENDIF
    ENDDO

    ! Initialiser les valeurs au niveau d'ascendance

    DO i = 1, klon
       kctop(i) = klev-1
       IF (.NOT. ldcum(i)) THEN
          kcbot(i) = klev-1
          pmfub(i) = 0.
          pqu(i, klev) = 0.
       ENDIF
       pmfu(i, klev) = pmfub(i)
       pmfus(i, klev) = pmfub(i) * (RCPD * ptu(i, klev)+pgeoh(i, klev))
       pmfuq(i, klev) = pmfub(i) * pqu(i, klev)
    ENDDO

    DO i = 1, klon
       ldcum(i) = .FALSE.
    ENDDO

    ! Do ascent: subcloud layer (klab=1), clouds (klab=2) by doing
    ! first dry-adiabatic ascent and then by adjusting t, q and l
    ! accordingly in flxadjtq, then check for buoyancy and set flags
    ! accordingly.

    DO k = klev - 1, 3, -1
       IF (LMFMID .AND. k < klev - 1 .AND. k > klev / 2) THEN
          DO i = 1, klon
             IF (.NOT. ldcum(i) .AND. klab(i, k + 1) == 0 .AND. &
                  pqen(i, k) > 0.9 * pqsen(i, k)) THEN
                ptu(i, k+1) = pten(i, k) +(pgeo(i, k)-pgeoh(i, k+1))/RCPD
                pqu(i, k+1) = pqen(i, k)
                plu(i, k+1) = 0.0
                zzzmb = MAX(CMFCMIN, -pvervel(i, k)/RG)
                zmfmax = (paph(i, k) - paph(i, k-1)) / (RG * pdtime)
                pmfub(i) = MIN(zzzmb, zmfmax)
                pmfu(i, k+1) = pmfub(i)
                pmfus(i, k+1) = pmfub(i) * (RCPD * ptu(i, k+1)+pgeoh(i, k+1))
                pmfuq(i, k+1) = pmfub(i) * pqu(i, k+1)
                pmful(i, k+1) = 0.0
                pdmfup(i, k+1) = 0.0
                kcbot(i) = k
                klab(i, k+1) = 1
                ktype(i) = 3
                pentr(i) = ENTRMID
             ENDIF
          ENDDO
       ENDIF

       is = 0
       DO i = 1, klon
          is = is + klab(i, k+1)
          IF (klab(i, k+1) == 0) klab(i, k) = 0
          llflag(i) = .FALSE.
          IF (klab(i, k+1) > 0) llflag(i) = .TRUE.
       ENDDO
       IF (is == 0) cycle

       ! Calculer le taux d'entraînement et de détraînement :

       DO i = 1, klon
          pen_u(i, k) = 0.0
          pde_u(i, k) = 0.0
          zrho(i) = paph(i, k + 1) / (RD * ptenh(i, k + 1))
          zpbot(i) = paph(i, kcbot(i))
          zptop(i) = paph(i, kctop0(i))
       ENDDO

       DO i = 1, klon
          IF (ldcum(i)) THEN
             zdprho = (paph(i, k + 1) - paph(i, k)) / (RG * zrho(i))
             zentr=pentr(i) * pmfu(i, k+1) * zdprho
             llo1=k < kcbot(i)
             IF (llo1) pde_u(i, k)=zentr
             zpmid=0.5 * (zpbot(i)+zptop(i))
             llo2 = llo1 .AND. ktype(i) == 2 &
                  .AND. (zpbot(i) - paph(i, k) < 0.2E5 .OR. paph(i, k) > zpmid)
             IF (llo2) pen_u(i, k)=zentr
             llo2 = llo1 .AND. (ktype(i) == 1 .OR. ktype(i) == 3) .AND. &
                  (k >= MAX(klwmin(i), kctop0(i) + 2) .OR. pap(i, k) > zpmid)
             IF (llo2) pen_u(i, k)=zentr
             llo1=pen_u(i, k) > 0. .AND. (ktype(i) == 1 .OR. ktype(i) == 2)
             IF (llo1) THEN
                fact = 1. + 3. * (1. - MIN(1., (zpbot(i) - pap(i, k)) / 1.5E4))
                zentr = zentr * fact
                pen_u(i, k)=pen_u(i, k) * fact
                pde_u(i, k)=pde_u(i, k) * fact
             ENDIF
             IF (llo2 .AND. pqenh(i, k+1) > 1e-5) &
                  pen_u(i, k)=zentr+MAX(pqte(i, k), 0.)/pqenh(i, k+1) * &
                  zrho(i) * zdprho
          ENDIF
       end DO

       ! Do adiabatic ascent for entraining/detraining plume

       DO i = 1, klon
          IF (llflag(i)) THEN
             IF (k < kcbot(i)) THEN
                zmftest = pmfu(i, k+1)+pen_u(i, k)-pde_u(i, k)
                zmfmax = MIN(zmftest, &
                     (paph(i, k) - paph(i, k - 1)) / (RG * pdtime))
                pen_u(i, k)=MAX(pen_u(i, k)-MAX(0.0, zmftest-zmfmax), 0.0)
             ENDIF
             pde_u(i, k)=MIN(pde_u(i, k), 0.75 * pmfu(i, k+1))
             ! calculer le flux de masse du niveau k a partir de celui du k+1
             pmfu(i, k)=pmfu(i, k+1)+pen_u(i, k)-pde_u(i, k)
             ! calculer les valeurs Su, Qu et l du niveau k dans le
             ! panache montant
             zqeen=pqenh(i, k+1) * pen_u(i, k)
             zseen=(RCPD * ptenh(i, k+1)+pgeoh(i, k+1)) * pen_u(i, k)
             zscde=(RCPD * ptu(i, k+1)+pgeoh(i, k+1)) * pde_u(i, k)
             zqude=pqu(i, k+1) * pde_u(i, k)
             plude(i, k)=plu(i, k+1) * pde_u(i, k)
             zmfusk=pmfus(i, k+1)+zseen-zscde
             zmfuqk=pmfuq(i, k+1)+zqeen-zqude
             zmfulk=pmful(i, k+1) -plude(i, k)
             plu(i, k)=zmfulk * (1./MAX(CMFCMIN, pmfu(i, k)))
             pqu(i, k)=zmfuqk * (1./MAX(CMFCMIN, pmfu(i, k)))
             ptu(i, k)=(zmfusk * (1./MAX(CMFCMIN, pmfu(i, k)))- &
                  pgeoh(i, k))/RCPD
             ptu(i, k)=MAX(100., ptu(i, k))
             ptu(i, k)=MIN(400., ptu(i, k))
             zqold(i)=pqu(i, k)
          ELSE
             zqold(i)=0.0
          ENDIF
       end DO

       ! Do corrections for moist ascent by adjusting t, q and l

       icall = 1
       CALL flxadjtq(paph(1, k), ptu(1, k), pqu(1, k), llflag, icall)

       DO i = 1, klon
          IF (llflag(i) .AND. pqu(i, k).NE.zqold(i)) THEN
             klab(i, k) = 2
             plu(i, k) = plu(i, k)+zqold(i)-pqu(i, k)
             zbuo = ptu(i, k) * (1.+RETV * pqu(i, k))- &
                  ptenh(i, k) * (1.+RETV * pqenh(i, k))
             IF (klab(i, k+1) == 1) zbuo=zbuo+0.5
             IF (zbuo > 0. .AND. pmfu(i, k) >= 0.1 * pmfub(i)) THEN
                kctop(i) = k
                ldcum(i) = .TRUE.
                zdnoprc = 1.5E4
                IF (ldland(i)) zdnoprc = zdland(i)
                zprcon = CPRCON
                IF ((zpbot(i) - paph(i, k)) < zdnoprc) zprcon = 0.
                zlnew=plu(i, k)/(1.+zprcon * (pgeoh(i, k)-pgeoh(i, k+1)))
                pdmfup(i, k)=MAX(0., (plu(i, k)-zlnew) * pmfu(i, k))
                plu(i, k)=zlnew
             ELSE
                klab(i, k)=0
                pmfu(i, k)=0.
             ENDIF
          ENDIF
       end DO
       DO i = 1, klon
          IF (llflag(i)) THEN
             pmful(i, k)=plu(i, k) * pmfu(i, k)
             pmfus(i, k)=(RCPD * ptu(i, k)+pgeoh(i, k)) * pmfu(i, k)
             pmfuq(i, k)=pqu(i, k) * pmfu(i, k)
          ENDIF
       end DO
    end DO

    ! Determine convective fluxes above non-buoyancy level (note:
    ! cloud variables like t, q and l are not affected by detrainment
    ! and are already known from previous calculations above).

    DO i = 1, klon
       IF (kctop(i) == klev-1) ldcum(i) = .FALSE.
       kcbot(i) = MAX(kcbot(i), kctop(i))
    ENDDO

    ldcum(1)=ldcum(1)

    is = 0
    DO i = 1, klon
       if (ldcum(i)) is = is + 1
    ENDDO
    kcum = is
    IF (is /= 0) then
       DO i = 1, klon
          IF (ldcum(i)) THEN
             k=kctop(i)-1
             pde_u(i, k)=(1.-CMFCTOP) * pmfu(i, k+1)
             plude(i, k)=pde_u(i, k) * plu(i, k+1)
             pmfu(i, k)=pmfu(i, k+1)-pde_u(i, k)
             zlnew=plu(i, k)
             pdmfup(i, k)=MAX(0., (plu(i, k)-zlnew) * pmfu(i, k))
             plu(i, k)=zlnew
             pmfus(i, k)=(RCPD * ptu(i, k)+pgeoh(i, k)) * pmfu(i, k)
             pmfuq(i, k)=pqu(i, k) * pmfu(i, k)
             pmful(i, k)=plu(i, k) * pmfu(i, k)
             plude(i, k-1)=pmful(i, k)
          ENDIF
       end DO
    end IF

  END SUBROUTINE flxasc

end module flxasc_m
