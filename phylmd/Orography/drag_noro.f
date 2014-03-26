module drag_noro_m

  IMPLICIT NONE

contains

  SUBROUTINE drag_noro(nlon, nlev, dtime, paprs, pplay, pmea, pstd, psig, &
       pgam, pthe, ppic, pval, kgwd, kdx, ktest, t, u, v, pulow, pvlow, &
       pustr, pvstr, d_t, d_u, d_v)

    ! From LMDZ4/libf/phylmd/orografi.F, version 1.4 2005/12/01 11:27:29

    USE dimphy, ONLY : klev, klon
    USE suphec_m, ONLY : rd, rg

    ! Auteur(s): F.Lott (LMD/CNRS) date: 19950201
    ! Objet: Frottement de la montagne Interface
    !======================================================================
    ! Arguments:
    ! dtime---input-R- pas d'integration (s)
    ! paprs---input-R-pression pour chaque inter-couche (en Pa)
    ! pplay---input-R-pression pour le mileu de chaque couche (en Pa)
    ! t-------input-R-temperature (K)
    ! u-------input-R-vitesse horizontale (m/s)
    ! v-------input-R-vitesse horizontale (m/s)

    ! d_t-----output-R-increment de la temperature
    ! d_u-----output-R-increment de la vitesse u
    ! d_v-----output-R-increment de la vitesse v
    !======================================================================

    ! ARGUMENTS

    INTEGER nlon, nlev
    REAL, INTENT (IN) :: dtime
    REAL, INTENT (IN) :: paprs(klon, klev+1)
    REAL, INTENT (IN) :: pplay(klon, klev)
    REAL pmea(nlon)
    REAL, INTENT (IN):: pstd(nlon), psig(nlon)
    REAL pgam(nlon), pthe(nlon)
    REAL ppic(nlon), pval(nlon)
    REAL pulow(nlon), pvlow(nlon), pustr(nlon), pvstr(nlon)
    REAL, INTENT (IN):: t(nlon, nlev)
    real, INTENT (IN):: u(nlon, nlev), v(nlon, nlev)
    REAL d_t(nlon, nlev), d_u(nlon, nlev), d_v(nlon, nlev)

    INTEGER i, k, kgwd, kdx(nlon), ktest(nlon)

    ! Variables locales:

    REAL zgeom(klon, klev)
    REAL pdtdt(klon, klev), pdudt(klon, klev), pdvdt(klon, klev)
    REAL pt(klon, klev), pu(klon, klev), pv(klon, klev)
    REAL papmf(klon, klev), papmh(klon, klev+1)

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
          pt(i, k) = t(i, klev-k+1)
          pu(i, k) = u(i, klev-k+1)
          pv(i, k) = v(i, klev-k+1)
          papmf(i, k) = pplay(i, klev-k+1)
       END DO
    END DO
    DO k = 1, klev + 1
       DO i = 1, klon
          papmh(i, k) = paprs(i, klev-k+2)
       END DO
    END DO
    DO i = 1, klon
       zgeom(i, klev) = rd*pt(i, klev)*log(papmh(i, klev+1)/papmf(i, klev))
    END DO
    DO k = klev - 1, 1, -1
       DO i = 1, klon
          zgeom(i, k) = zgeom(i, k + 1) + rd * (pt(i, k) + pt(i, k + 1)) / 2. &
               * log(papmf(i, k + 1) / papmf(i, k))
       END DO
    END DO

    ! appeler la routine principale

    CALL orodrag(klon, klev, kgwd, kdx, ktest, dtime, papmh, papmf, zgeom, &
         pt, pu, pv, pmea, pstd, psig, pgam, pthe, ppic, pval, pulow, pvlow, &
         pdudt, pdvdt, pdtdt)

    DO k = 1, klev
       DO i = 1, klon
          d_u(i, klev+1-k) = dtime*pdudt(i, k)
          d_v(i, klev+1-k) = dtime*pdvdt(i, k)
          d_t(i, klev+1-k) = dtime*pdtdt(i, k)
          pustr(i) = pustr(i) & 
               + pdudt(i, k)*(papmh(i, k+1)-papmh(i, k))/rg
          pvstr(i) = pvstr(i) & 
               + pdvdt(i, k)*(papmh(i, k+1)-papmh(i, k))/rg
       END DO
    END DO

  END SUBROUTINE drag_noro

end module drag_noro_m
