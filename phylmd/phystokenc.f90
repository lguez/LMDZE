module phystokenc_m

  IMPLICIT NONE

contains

  SUBROUTINE phystokenc(pdtphys, rlon, rlat, pt, pmfu, pmfd, pen_u, pde_u, &
       pen_d, pde_d, pfm_therm, pentr_therm, pcoefh, yu1, yv1, ftsol, pctsrf, &
       frac_impa, frac_nucl, pphis, paire, dtime, itap)

    ! From phylmd/phystokenc.F, version 1.2 2004/06/22 11:45:35
    ! Author: Frédéric Hourdin
    ! Objet : écriture des variables pour transport offline

    USE histwrite_m, ONLY: histwrite
    USE histsync_m, ONLY: histsync
    USE dimens_m, ONLY: iim, jjm, nqmx
    USE indicesol, ONLY: nbsrf
    USE dimphy, ONLY: klev, klon
    USE tracstoke, ONLY: istphy

    REAL, INTENT (IN):: pdtphys ! pas d'integration pour la physique (seconde)
    REAL, INTENT (IN):: rlon(klon), rlat(klon)
    REAL, intent(in):: pt(klon, klev)

    ! convection: 

    REAL, INTENT (IN):: pmfu(klon, klev) ! flux de masse dans le panache montant

    REAL, intent(in):: pmfd(klon, klev)
    ! flux de masse dans le panache descendant

    REAL, intent(in):: pen_u(klon, klev) ! flux entraine dans le panache montant
    REAL, intent(in):: pde_u(klon, klev) ! flux detraine dans le panache montant

    REAL, intent(in):: pen_d(klon, klev) 
    ! flux entraine dans le panache descendant

    REAL, intent(in):: pde_d(klon, klev) 
    ! flux detraine dans le panache descendant

    ! Les Thermiques
    REAL pfm_therm(klon, klev+1)
    REAL pentr_therm(klon, klev)

    ! Couche limite: 

    REAL pcoefh(klon, klev) ! coeff melange Couche limite
    REAL yu1(klon)
    REAL yv1(klon)

    ! Arguments necessaires pour les sources et puits de traceur 

    REAL ftsol(klon, nbsrf) ! Temperature du sol (surf)(Kelvin)
    REAL pctsrf(klon, nbsrf) ! Pourcentage de sol f(nature du sol)

    ! Lessivage: 

    REAL frac_impa(klon, klev)
    REAL frac_nucl(klon, klev)

    REAL, INTENT(IN):: pphis(klon)
    real paire(klon)
    REAL, INTENT (IN):: dtime
    INTEGER, INTENT (IN):: itap

    ! Variables local to the procedure:

    real t(klon, klev)
    INTEGER, SAVE:: physid
    REAL zx_tmp_3d(iim, jjm+1, klev), zx_tmp_2d(iim, jjm+1)

    ! Les Thermiques

    REAL fm_therm1(klon, klev)
    REAL entr_therm(klon, klev)
    REAL fm_therm(klon, klev)

    INTEGER i, k

    REAL, save:: mfu(klon, klev) ! flux de masse dans le panache montant
    REAL mfd(klon, klev) ! flux de masse dans le panache descendant
    REAL en_u(klon, klev) ! flux entraine dans le panache montant
    REAL de_u(klon, klev) ! flux detraine dans le panache montant
    REAL en_d(klon, klev) ! flux entraine dans le panache descendant
    REAL de_d(klon, klev) ! flux detraine dans le panache descendant
    REAL coefh(klon, klev) ! flux detraine dans le panache descendant

    REAL pyu1(klon), pyv1(klon)
    REAL pftsol(klon, nbsrf), ppsrf(klon, nbsrf)
    REAL pftsol1(klon), pftsol2(klon), pftsol3(klon), pftsol4(klon)
    REAL ppsrf1(klon), ppsrf2(klon), ppsrf3(klon), ppsrf4(klon)

    REAL dtcum

    INTEGER:: iadvtr = 0, irec = 1
    REAL zmin, zmax
    LOGICAL ok_sync

    SAVE t, mfd, en_u, de_u, en_d, de_d, coefh, dtcum
    SAVE fm_therm, entr_therm
    SAVE pyu1, pyv1, pftsol, ppsrf

    !------------------------------------------------------

    ! Couche limite: 

    ok_sync = .TRUE.

    IF (iadvtr==0) THEN
       CALL initphysto('phystoke', rlon, rlat, dtime, dtime*istphy, &
            dtime*istphy, nqmx, physid)
    END IF

    i = itap
    CALL gr_fi_ecrit(1, klon, iim, jjm+1, pphis, zx_tmp_2d)
    CALL histwrite(physid, 'phis', i, zx_tmp_2d)

    i = itap
    CALL gr_fi_ecrit(1, klon, iim, jjm+1, paire, zx_tmp_2d)
    CALL histwrite(physid, 'aire', i, zx_tmp_2d)

    iadvtr = iadvtr + 1

    IF (mod(iadvtr, istphy) == 1 .OR. istphy == 1) THEN
       PRINT *, 'reinitialisation des champs cumules a iadvtr =', iadvtr
       DO k = 1, klev
          DO i = 1, klon
             mfu(i, k) = 0.
             mfd(i, k) = 0.
             en_u(i, k) = 0.
             de_u(i, k) = 0.
             en_d(i, k) = 0.
             de_d(i, k) = 0.
             coefh(i, k) = 0.
             t(i, k) = 0.
             fm_therm(i, k) = 0.
             entr_therm(i, k) = 0.
          END DO
       END DO
       DO i = 1, klon
          pyv1(i) = 0.
          pyu1(i) = 0.
       END DO
       DO k = 1, nbsrf
          DO i = 1, klon
             pftsol(i, k) = 0.
             ppsrf(i, k) = 0.
          END DO
       END DO

       dtcum = 0.
    END IF

    DO k = 1, klev
       DO i = 1, klon
          mfu(i, k) = mfu(i, k) + pmfu(i, k)*pdtphys
          mfd(i, k) = mfd(i, k) + pmfd(i, k)*pdtphys
          en_u(i, k) = en_u(i, k) + pen_u(i, k)*pdtphys
          de_u(i, k) = de_u(i, k) + pde_u(i, k)*pdtphys
          en_d(i, k) = en_d(i, k) + pen_d(i, k)*pdtphys
          de_d(i, k) = de_d(i, k) + pde_d(i, k)*pdtphys
          coefh(i, k) = coefh(i, k) + pcoefh(i, k)*pdtphys
          t(i, k) = t(i, k) + pt(i, k)*pdtphys
          fm_therm(i, k) = fm_therm(i, k) + pfm_therm(i, k)*pdtphys
          entr_therm(i, k) = entr_therm(i, k) + pentr_therm(i, k)*pdtphys
       END DO
    END DO
    DO i = 1, klon
       pyv1(i) = pyv1(i) + yv1(i)*pdtphys
       pyu1(i) = pyu1(i) + yu1(i)*pdtphys
    END DO
    DO k = 1, nbsrf
       DO i = 1, klon
          pftsol(i, k) = pftsol(i, k) + ftsol(i, k)*pdtphys
          ppsrf(i, k) = ppsrf(i, k) + pctsrf(i, k)*pdtphys
       END DO
    END DO

    dtcum = dtcum + pdtphys

    IF (mod(iadvtr, istphy) == 0) THEN
       ! normalisation par le temps cumule 
       DO k = 1, klev
          DO i = 1, klon
             mfu(i, k) = mfu(i, k)/dtcum
             mfd(i, k) = mfd(i, k)/dtcum
             en_u(i, k) = en_u(i, k)/dtcum
             de_u(i, k) = de_u(i, k)/dtcum
             en_d(i, k) = en_d(i, k)/dtcum
             de_d(i, k) = de_d(i, k)/dtcum
             coefh(i, k) = coefh(i, k)/dtcum
             ! Unitel a enlever
             t(i, k) = t(i, k)/dtcum
             fm_therm(i, k) = fm_therm(i, k)/dtcum
             entr_therm(i, k) = entr_therm(i, k)/dtcum
          END DO
       END DO
       DO i = 1, klon
          pyv1(i) = pyv1(i)/dtcum
          pyu1(i) = pyu1(i)/dtcum
       END DO
       DO k = 1, nbsrf
          DO i = 1, klon
             pftsol(i, k) = pftsol(i, k)/dtcum
             pftsol1(i) = pftsol(i, 1)
             pftsol2(i) = pftsol(i, 2)
             pftsol3(i) = pftsol(i, 3)
             pftsol4(i) = pftsol(i, 4)

             ppsrf(i, k) = ppsrf(i, k)/dtcum
             ppsrf1(i) = ppsrf(i, 1)
             ppsrf2(i) = ppsrf(i, 2)
             ppsrf3(i) = ppsrf(i, 3)
             ppsrf4(i) = ppsrf(i, 4)
          END DO
       END DO

       ! ecriture des champs 

       irec = irec + 1

       CALL gr_fi_ecrit(klev, klon, iim, jjm+1, t, zx_tmp_3d)
       CALL histwrite(physid, 't', itap, zx_tmp_3d)

       CALL gr_fi_ecrit(klev, klon, iim, jjm+1, mfu, zx_tmp_3d)
       CALL histwrite(physid, 'mfu', itap, zx_tmp_3d)
       CALL gr_fi_ecrit(klev, klon, iim, jjm+1, mfd, zx_tmp_3d)
       CALL histwrite(physid, 'mfd', itap, zx_tmp_3d)
       CALL gr_fi_ecrit(klev, klon, iim, jjm+1, en_u, zx_tmp_3d)
       CALL histwrite(physid, 'en_u', itap, zx_tmp_3d)
       CALL gr_fi_ecrit(klev, klon, iim, jjm+1, de_u, zx_tmp_3d)
       CALL histwrite(physid, 'de_u', itap, zx_tmp_3d)
       CALL gr_fi_ecrit(klev, klon, iim, jjm+1, en_d, zx_tmp_3d)
       CALL histwrite(physid, 'en_d', itap, zx_tmp_3d)
       CALL gr_fi_ecrit(klev, klon, iim, jjm+1, de_d, zx_tmp_3d)
       CALL histwrite(physid, 'de_d', itap, zx_tmp_3d)
       CALL gr_fi_ecrit(klev, klon, iim, jjm+1, coefh, zx_tmp_3d)
       CALL histwrite(physid, 'coefh', itap, zx_tmp_3d)

       DO k = 1, klev
          DO i = 1, klon
             fm_therm1(i, k) = fm_therm(i, k)
          END DO
       END DO

       CALL gr_fi_ecrit(klev, klon, iim, jjm+1, fm_therm1, zx_tmp_3d)
       CALL histwrite(physid, 'fm_th', itap, zx_tmp_3d)

       CALL gr_fi_ecrit(klev, klon, iim, jjm+1, entr_therm, zx_tmp_3d)
       CALL histwrite(physid, 'en_th', itap, zx_tmp_3d)
       !ccc 
       CALL gr_fi_ecrit(klev, klon, iim, jjm+1, frac_impa, zx_tmp_3d)
       CALL histwrite(physid, 'frac_impa', itap, zx_tmp_3d)

       CALL gr_fi_ecrit(klev, klon, iim, jjm+1, frac_nucl, zx_tmp_3d)
       CALL histwrite(physid, 'frac_nucl', itap, zx_tmp_3d)

       CALL gr_fi_ecrit(1, klon, iim, jjm+1, pyu1, zx_tmp_2d)
       CALL histwrite(physid, 'pyu1', itap, zx_tmp_2d)

       CALL gr_fi_ecrit(1, klon, iim, jjm+1, pyv1, zx_tmp_2d)
       CALL histwrite(physid, 'pyv1', itap, zx_tmp_2d)

       CALL gr_fi_ecrit(1, klon, iim, jjm+1, pftsol1, zx_tmp_2d)
       CALL histwrite(physid, 'ftsol1', itap, zx_tmp_2d)
       CALL gr_fi_ecrit(1, klon, iim, jjm+1, pftsol2, zx_tmp_2d)
       CALL histwrite(physid, 'ftsol2', itap, zx_tmp_2d)
       CALL gr_fi_ecrit(1, klon, iim, jjm+1, pftsol3, zx_tmp_2d)
       CALL histwrite(physid, 'ftsol3', itap, zx_tmp_2d)
       CALL gr_fi_ecrit(1, klon, iim, jjm+1, pftsol4, zx_tmp_2d)
       CALL histwrite(physid, 'ftsol4', itap, zx_tmp_2d)

       CALL gr_fi_ecrit(1, klon, iim, jjm+1, ppsrf1, zx_tmp_2d)
       CALL histwrite(physid, 'psrf1', itap, zx_tmp_2d)
       CALL gr_fi_ecrit(1, klon, iim, jjm+1, ppsrf2, zx_tmp_2d)
       CALL histwrite(physid, 'psrf2', itap, zx_tmp_2d)
       CALL gr_fi_ecrit(1, klon, iim, jjm+1, ppsrf3, zx_tmp_2d)
       CALL histwrite(physid, 'psrf3', itap, zx_tmp_2d)
       CALL gr_fi_ecrit(1, klon, iim, jjm+1, ppsrf4, zx_tmp_2d)
       CALL histwrite(physid, 'psrf4', itap, zx_tmp_2d)

       IF (ok_sync) CALL histsync(physid)

       ! Test sur la valeur des coefficients de lessivage 

       zmin = 1E33
       zmax = -1E33
       DO k = 1, klev
          DO i = 1, klon
             zmax = max(zmax, frac_nucl(i, k))
             zmin = min(zmin, frac_nucl(i, k))
          END DO
       END DO
       PRINT *, 'coefs de lessivage (min et max)'
       PRINT *, 'facteur de nucleation ', zmin, zmax
       zmin = 1E33
       zmax = -1E33
       DO k = 1, klev
          DO i = 1, klon
             zmax = max(zmax, frac_impa(i, k))
             zmin = min(zmin, frac_impa(i, k))
          END DO
       END DO
       PRINT *, 'facteur d impaction ', zmin, zmax
    END IF

  END SUBROUTINE phystokenc

end module phystokenc_m
