module lift_noro_m

  IMPLICIT NONE

contains

  SUBROUTINE lift_noro(paprs, play, zmea, zstd, zpic, t_seri, u_seri, v_seri, &
       ustrli, vstrli, d_t_lif, d_u_lif, d_v_lif)

    ! Author: F.Lott (LMD/CNRS) date: 1995/02/01
    ! Objet: Frottement de la montagne, interface

    use conf_gcm_m, only: dtphys
    USE dimphy, only: klon, klev
    use orolift_m, only: orolift
    use phyetat0_m, only: rlat
    USE suphec_m, only: rd, rg

    REAL, INTENT (IN) :: paprs(klon, klev + 1)
    ! pression pour chaque inter-couche (en Pa)

    REAL, INTENT (IN) :: play(klon, klev)
    ! pression pour le mileu de chaque couche (en Pa)

    REAL, INTENT (IN):: zmea(klon)
    REAL, INTENT (IN):: zstd(klon)
    REAL, INTENT (IN):: zpic(klon)

    REAL, INTENT (IN):: t_seri(klon, klev) ! temperature (K)

    real, INTENT (IN):: u_seri(klon, klev), v_seri(klon, klev)
    ! vitesse horizontale (m / s)

    REAL, intent(out):: ustrli(klon), vstrli(klon)
    REAL, intent(out):: d_t_lif(klon, klev) ! incr\'ement de la temperature

    REAL, intent(out):: d_u_lif(klon, klev), d_v_lif(klon, klev)
    ! incr\'ement de la vitesse u_seri, v_seri

    ! Local:
    INTEGER i, k
    REAL zgeom(klon, klev)
    REAL pdtdt(klon, klev), pdudt(klon, klev), pdvdt(klon, klev)
    REAL pt(klon, klev), pu(klon, klev), pv(klon, klev)
    REAL papmf(klon, klev), papmh(klon, klev + 1)

    !----------------------------------------------------------------------

    ! initialiser les variables de sortie (pour securite)
    pdudt = 0.
    pdvdt = 0.
    pdtdt = 0.

    ! preparer les variables d'entree (attention: l'ordre des niveaux
    ! verticaux augmente du haut vers le bas)

    DO k = 1, klev
       DO i = 1, klon
          pt(i, k) = t_seri(i, klev-k + 1)
          pu(i, k) = u_seri(i, klev-k + 1)
          pv(i, k) = v_seri(i, klev-k + 1)
          papmf(i, k) = play(i, klev-k + 1)
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
    CALL orolift(klon, klev, dtphys, papmh, zgeom, pt, pu, pv, rlat, zmea, &
         zstd, zpic, pdudt, pdvdt, pdtdt)

    ustrli = 0.
    vstrli = 0.

    DO k = 1, klev
       DO i = 1, klon
          d_u_lif(i, klev + 1-k) = dtphys * pdudt(i, k)
          d_v_lif(i, klev + 1-k) = dtphys * pdvdt(i, k)
          d_t_lif(i, klev + 1-k) = dtphys * pdtdt(i, k)
          ustrli(i) = ustrli(i) &
               + pdudt(i, k) * (papmh(i, k + 1)-papmh(i, k)) / rg
          vstrli(i) = vstrli(i) &
               + pdvdt(i, k) * (papmh(i, k + 1)-papmh(i, k)) / rg
       END DO
    END DO

  END SUBROUTINE lift_noro

end module lift_noro_m
