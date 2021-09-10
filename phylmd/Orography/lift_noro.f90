module lift_noro_m

  IMPLICIT NONE

contains

  SUBROUTINE lift_noro(paprs, pplay, pmea, pstd, ppic, t, u, v, pulow, pvlow, &
       pustr, pvstr, d_t, d_u, d_v, ktest)

    ! Author: F.Lott (LMD/CNRS) date: 1995/02/01
    ! Objet: Frottement de la montagne, interface

    use conf_gcm_m, only: dtphys
    USE dimphy, only: klon, klev
    use orolift_m, only: orolift
    use phyetat0_m, only: rlat
    USE suphec_m, only: rd, rg
    
    REAL, INTENT (IN) :: paprs(klon, klev + 1)
    ! paprs---input-R-pression pour chaque inter-couche (en Pa)
    REAL, INTENT (IN) :: pplay(klon, klev)
    ! pplay---input-R-pression pour le mileu de chaque couche (en Pa)
    REAL pmea(klon)
    REAL, INTENT (IN):: pstd(klon)
    REAL ppic(klon)
    REAL, INTENT (IN):: t(klon, klev)
    ! t-------input-R-temperature (K)
    real, INTENT (IN):: u(klon, klev), v(klon, klev)
    ! u-------input-R-vitesse horizontale (m / s)
    ! v-------input-R-vitesse horizontale (m / s)
    REAL pulow(klon), pvlow(klon), pustr(klon), pvstr(klon)
    REAL d_t(klon, klev), d_u(klon, klev), d_v(klon, klev)
    ! d_t-----output-R-increment de la temperature
    ! d_u-----output-R-increment de la vitesse u
    ! d_v-----output-R-increment de la vitesse v

    integer, INTENT(IN):: ktest(klon)

    ! Local:
    INTEGER i, k
    REAL zgeom(klon, klev)
    REAL pdtdt(klon, klev), pdudt(klon, klev), pdvdt(klon, klev)
    REAL pt(klon, klev), pu(klon, klev), pv(klon, klev)
    REAL papmf(klon, klev), papmh(klon, klev + 1)

    !----------------------------------------------------------------------

    ! initialiser les variables de sortie (pour securite)

    DO i = 1, klon
       pulow(i) = 0.0
       pvlow(i) = 0.0
       pustr(i) = 0.0
       pvstr(i) = 0.0
    END DO
    DO k = 1, klev
       DO i = 1, klon
          d_t(i, k) = 0.0
          d_u(i, k) = 0.0
          d_v(i, k) = 0.0
          pdudt(i, k) = 0.0
          pdvdt(i, k) = 0.0
          pdtdt(i, k) = 0.0
       END DO
    END DO

    ! preparer les variables d'entree (attention: l'ordre des niveaux
    ! verticaux augmente du haut vers le bas)

    DO k = 1, klev
       DO i = 1, klon
          pt(i, k) = t(i, klev-k + 1)
          pu(i, k) = u(i, klev-k + 1)
          pv(i, k) = v(i, klev-k + 1)
          papmf(i, k) = pplay(i, klev-k + 1)
       END DO
    END DO
    DO k = 1, klev + 1
       DO i = 1, klon
          papmh(i, k) = paprs(i, klev-k + 2)
       END DO
    END DO
    DO i = 1, klon
       zgeom(i, klev) = rd * pt(i, klev) &
            * log(papmh(i, klev + 1) / papmf(i, klev))
    END DO
    DO k = klev - 1, 1, -1
       DO i = 1, klon
          zgeom(i, k) = zgeom(i, k + 1) + rd * (pt(i, k) + pt(i, k + 1)) &
               / 2. * log(papmf(i, k + 1) / papmf(i, k))
       END DO
    END DO

    ! appeler la routine principale

    CALL orolift(klon, klev, ktest, dtphys, papmh, zgeom, pt, pu, pv, rlat, &
         pmea, pstd, ppic, pulow, pvlow, pdudt, pdvdt, pdtdt)

    DO k = 1, klev
       DO i = 1, klon
          d_u(i, klev + 1-k) = dtphys * pdudt(i, k)
          d_v(i, klev + 1-k) = dtphys * pdvdt(i, k)
          d_t(i, klev + 1-k) = dtphys * pdtdt(i, k)
          pustr(i) = pustr(i) &
               + pdudt(i, k) * (papmh(i, k + 1)-papmh(i, k)) / rg
          pvstr(i) = pvstr(i) &
               + pdvdt(i, k) * (papmh(i, k + 1)-papmh(i, k)) / rg
       END DO
    END DO

  END SUBROUTINE lift_noro

end module lift_noro_m
