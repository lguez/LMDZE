module flxflux_m

  IMPLICIT none

contains

  SUBROUTINE flxflux(pdtime, pqen, pqsen, ptenh, pqenh, pap, paph, ldland, &
       pgeoh, kcbot, kctop, lddraf, kdtop, ktype, ldcum, pmfu, pmfd, pmfus, &
       pmfds, pmfuq, pmfdq, pmful, plude, pdmfup, pdmfdp, pten, prfl, psfl, &
       pdpmel, ktopm2, pmflxr, pmflxs)

    ! This routine does the final calculation of convective fluxes in
    ! the cloud layer and in the subcloud layer.

    USE dimphy, ONLY: klev, klon
    USE suphec_m, ONLY: rcpd, retv, rg, rlmlt, rtt
    USE yoethf_m, ONLY: r2es
    USE fcttre, ONLY: foeew

    REAL, intent(in):: pdtime
    REAL, intent(in):: pqen(klon, klev)
    real, intent(inout):: pqsen(klon, klev)
    REAL, intent(in):: ptenh(klon, klev)
    REAL, intent(in):: pqenh(klon, klev)
    REAL pap(klon, klev)
    REAL paph(klon, klev+1)
    LOGICAL ldland(klon)
    real pgeoh(klon, klev)
    INTEGER, intent(in):: kcbot(klon), kctop(klon)
    LOGICAL lddraf(klon)
    INTEGER kdtop(klon)
    INTEGER, intent(inout):: ktype(klon)
    LOGICAL, intent(in):: ldcum(klon)
    REAL pmfu(klon, klev)
    REAL pmfd(klon, klev)
    real pmfus(klon, klev)
    real pmfds(klon, klev)
    REAL pmfuq(klon, klev)
    REAL pmfdq(klon, klev)
    real pmful(klon, klev)
    REAL plude(klon, klev)
    REAL pdmfup(klon, klev)
    REAL pdmfdp(klon, klev)
    REAL, intent(in):: pten(klon, klev)
    REAL prfl(klon), psfl(klon)
    real pdpmel(klon, klev)
    INTEGER, intent(out):: ktopm2
    REAL pmflxr(klon, klev+1), pmflxs(klon, klev+1)

    ! Local:
    REAL cevapcu(klev)
    REAL ztmsmlt, zdelta, zqsat

    !jq The variable maxpdmfdp(klon) has been introduced by Olivier Boucher
    !jq 14/11/00 to fix the problem with the negative precipitation.
    real maxpdmfdp(klon, klev)

    INTEGER k, kp, i
    REAL zcons1, zcons2, zcucov, ztmelp2
    real zdp, zzp, zfac, zsnmlt, zrfl, zrnew
    REAL zrmin, zrfln, zdrfl
    REAL zpds, zpdr, zdenom
    INTEGER ikb

    !---------------------------------------------------------------

    DO k=1, klev
       CEVAPCU(k)=1.93E-6*261.*SQRT(1.E3/(38.3*0.293) &
            *SQRT(0.5*(paph(1, k)+paph(1, k+1))/paph(1, klev+1))) * 0.5/RG
    end DO

    ! SPECIFY CONSTANTS

    zcons1 = RCPD/(RLMLT*RG*pdtime)
    zcons2 = 1./(RG*pdtime)
    zcucov = 0.05
    ztmelp2 = RTT + 2.

    ! DETERMINE FINAL CONVECTIVE FLUXES

    DO i = 1, klon
       IF (.NOT. ldcum(i) .OR. kdtop(i) < kctop(i)) lddraf(i)=.FALSE.
       IF (.NOT. ldcum(i)) ktype(i) = 0
    end DO

    ktopm2 = min(klev, minval(kctop)) - 2

    DO k = ktopm2, klev
       DO i = 1, klon
          IF (ldcum(i) .AND. k >= kctop(i) - 1) THEN
             pmfus(i, k) = pmfus(i, k) &
                  - pmfu(i, k) * (RCPD * ptenh(i, k) + pgeoh(i, k))
             pmfuq(i, k)=pmfuq(i, k)-pmfu(i, k)*pqenh(i, k)
             zdp = 1.5E4
             IF (ldland(i)) zdp = 3.E4

             ! L'eau liquide détrainée est precipitée quand certaines
             ! conditions sont réunies (sinon, elle est considérée
             ! évaporée dans l'environnement)

             IF (paph(i, kcbot(i)) - paph(i, kctop(i)) >= zdp &
                  .AND. pqen(i, k - 1) > 0.8 * pqsen(i, k - 1)) &
                  pdmfup(i, k - 1) = pdmfup(i, k - 1) + plude(i, k - 1)

             IF (lddraf(i).AND.k >= kdtop(i)) THEN
                pmfds(i, k)=pmfds(i, k)-pmfd(i, k)* &
                     (RCPD*ptenh(i, k)+pgeoh(i, k))
                pmfdq(i, k)=pmfdq(i, k)-pmfd(i, k)*pqenh(i, k)
             ELSE
                pmfd(i, k)=0.
                pmfds(i, k)=0.
                pmfdq(i, k)=0.
                pdmfdp(i, k-1)=0.
             END IF
          ELSE
             pmfu(i, k)=0.
             pmfus(i, k)=0.
             pmfuq(i, k)=0.
             pmful(i, k)=0.
             pdmfup(i, k-1)=0.
             plude(i, k-1)=0.
             pmfd(i, k)=0.
             pmfds(i, k)=0.
             pmfdq(i, k)=0.
             pdmfdp(i, k-1)=0.
          ENDIF
       end DO
    end DO

    DO k=ktopm2, klev
       DO i = 1, klon
          IF (ldcum(i) .AND. k > kcbot(i)) THEN
             ikb=kcbot(i)
             zzp = ((paph(i, klev + 1) - paph(i, k)) &
                  / (paph(i, klev + 1) - paph(i, ikb)))
             IF (ktype(i) == 3) zzp = zzp**2
             pmfu(i, k)=pmfu(i, ikb)*zzp
             pmfus(i, k)=pmfus(i, ikb)*zzp
             pmfuq(i, k)=pmfuq(i, ikb)*zzp
             pmful(i, k)=pmful(i, ikb)*zzp
          ENDIF
       end DO
    end DO

    ! CALCULATE RAIN/SNOW FALL RATES
    ! CALCULATE MELTING OF SNOW
    ! CALCULATE EVAPORATION OF PRECIP

    DO k = 1, klev+1
       DO i = 1, klon
          pmflxr(i, k) = 0.0
          pmflxs(i, k) = 0.0
       ENDDO
    ENDDO
    DO k = ktopm2, klev
       DO i = 1, klon
          IF (ldcum(i)) THEN
             IF (pmflxs(i, k) > 0.0 .AND. pten(i, k) > ztmelp2) THEN
                zfac=zcons1*(paph(i, k+1)-paph(i, k))
                zsnmlt=MIN(pmflxs(i, k), zfac*(pten(i, k)-ztmelp2))
                pdpmel(i, k)=zsnmlt
                ztmsmlt=pten(i, k)-zsnmlt/zfac
                zdelta=MAX(0., SIGN(1., RTT-ztmsmlt))
                zqsat = R2ES * FOEEW(ztmsmlt, zdelta) / pap(i, k)
                zqsat = MIN(0.5, zqsat)
                zqsat = zqsat / (1. - RETV * zqsat)
                pqsen(i, k) = zqsat
             ENDIF
             IF (pten(i, k) > RTT) THEN
                pmflxr(i, k + 1) = pmflxr(i, k) + pdmfup(i, k) + pdmfdp(i, k) &
                     + pdpmel(i, k)
                pmflxs(i, k+1)=pmflxs(i, k)-pdpmel(i, k)
             ELSE
                pmflxs(i, k+1)=pmflxs(i, k)+pdmfup(i, k)+pdmfdp(i, k)
                pmflxr(i, k+1)=pmflxr(i, k)
             ENDIF
             ! si la precipitation est negative, on ajuste le plux du
             ! panache descendant pour eliminer la negativite
             IF ((pmflxr(i, k+1)+pmflxs(i, k+1)) < 0.0) THEN
                pdmfdp(i, k) = -pmflxr(i, k)-pmflxs(i, k)-pdmfup(i, k)
                pmflxr(i, k+1) = 0.0
                pmflxs(i, k+1) = 0.0
                pdpmel(i, k) = 0.0
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    ! The new variable is initialized here.  It contains the
    ! humidity which is fed to the downdraft by evaporation of
    ! precipitation in the column below the base of convection.

    ! In the former version, this term has been subtracted from precip
    ! as well as the evaporation.

    DO k = 1, klev
       DO i = 1, klon
          maxpdmfdp(i, k)=0.0
       ENDDO
    ENDDO
    DO k = 1, klev
       DO kp = k, klev
          DO i = 1, klon
             maxpdmfdp(i, k)=maxpdmfdp(i, k)+pdmfdp(i, kp)
          ENDDO
       ENDDO
    ENDDO
    ! End of initialization

    DO k = ktopm2, klev
       DO i = 1, klon
          IF (ldcum(i) .AND. k >= kcbot(i)) THEN
             zrfl = pmflxr(i, k) + pmflxs(i, k)
             IF (zrfl > 1.0E-20) THEN
                zrnew = (MAX(0., SQRT(zrfl / zcucov) - CEVAPCU(k) &
                     * (paph(i, k + 1) - paph(i, k)) &
                     * MAX(0., pqsen(i, k) - pqen(i, k))))**2 * zcucov
                zrmin = zrfl - zcucov &
                     * MAX(0., 0.8 * pqsen(i, k) - pqen(i, k)) &
                     * zcons2 * (paph(i, k + 1) - paph(i, k))
                zrnew=MAX(zrnew, zrmin)
                zrfln=MAX(zrnew, 0.)
                zdrfl=MIN(0., zrfln-zrfl)
                ! At least the amount of precipiation needed to feed
                ! the downdraft with humidity below the base of
                ! convection has to be left and can't be evaporated
                ! (surely the evaporation can't be positive):
                zdrfl=MAX(zdrfl, &
                     MIN(-pmflxr(i, k)-pmflxs(i, k)-maxpdmfdp(i, k), 0.0))

                zdenom=1.0/MAX(1.0E-20, pmflxr(i, k)+pmflxs(i, k))
                IF (pten(i, k) > RTT) THEN
                   zpdr = pdmfdp(i, k)
                   zpds = 0.0
                ELSE
                   zpdr = 0.0
                   zpds = pdmfdp(i, k)
                ENDIF
                pmflxr(i, k+1) = pmflxr(i, k) + zpdr + pdpmel(i, k) &
                     + zdrfl*pmflxr(i, k)*zdenom
                pmflxs(i, k+1) = pmflxs(i, k) + zpds - pdpmel(i, k) &
                     + zdrfl*pmflxs(i, k)*zdenom
                pdmfup(i, k) = pdmfup(i, k) + zdrfl
             ELSE
                pmflxr(i, k+1) = 0.0
                pmflxs(i, k+1) = 0.0
                pdmfdp(i, k) = 0.0
                pdpmel(i, k) = 0.0
             ENDIF
             if (pmflxr(i, k) + pmflxs(i, k) < -1.e-26) &
                  write(*, *) 'precip. < 1e-16 ', pmflxr(i, k) + pmflxs(i, k)
          ENDIF
       ENDDO
    ENDDO

    DO i = 1, klon
       prfl(i) = pmflxr(i, klev+1)
       psfl(i) = pmflxs(i, klev+1)
    end DO

    RETURN
  END SUBROUTINE flxflux

end module flxflux_m
