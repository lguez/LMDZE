module drag_noro_m

  IMPLICIT NONE

contains

  SUBROUTINE drag_noro(paprs, play, zmea, zstd, zsig, zgam, zthe, zpic, zval, &
       t_seri, u_seri, v_seri, ustrdr, vstrdr, d_t_oro, d_u_oro, d_v_oro)

    ! From LMDZ4/libf/phylmd/orografi.F, version 1.4 2005/12/01 11:27:29

    ! Author: F. Lott (LMD/CNRS). Date: 1995/02/01.
    ! Objet : frottement de la montagne, interface.

    use conf_gcm_m, only: dtphys
    USE dimphy, ONLY: klev, klon
    use orodrag_m, only: orodrag
    USE suphec_m, ONLY: rd, rg

    REAL, INTENT(IN):: paprs(klon, klev + 1)
    ! pression pour chaque inter-couche (en Pa)

    REAL, INTENT(IN):: play(klon, klev)
    ! pression pour le mileu de chaque couche (en Pa)

    REAL, INTENT(IN):: zmea(klon)
    REAL, INTENT(IN):: zstd(klon), zsig(klon)
    REAL, INTENT(INout):: zgam(klon)
    real, INTENT(IN):: zthe(klon), zpic(klon), zval(klon)

    REAL, INTENT(IN):: t_seri(klon, klev) ! temperature (K)

    real, INTENT(IN):: u_seri(klon, klev), v_seri(klon, klev)
    ! vitesse horizontale (m / s)

    REAL, intent(out):: ustrdr(klon), vstrdr(klon)
    REAL, intent(out):: d_t_oro(klon, klev) ! increment de la temperature

    REAL, intent(out):: d_u_oro(klon, klev), d_v_oro(klon, klev)
    ! incr\'ement de la vitesse

    ! Local:
    INTEGER i, k
    REAL zgeom(klon, klev) ! variation of geopotential
    REAL pdtdt(klon, klev), pdudt(klon, klev), pdvdt(klon, klev)
    REAL tm1(klon, klev), um1(klon, klev), vm1(klon, klev)
    REAL papm1(klon, klev), paphm1(klon, klev + 1)

    !--------------------------------------------------------------------

    ! Initialiser les variables de sortie (pour securite)
    pdudt = 0.
    pdvdt = 0.
    pdtdt = 0.

    ! Preparer les variables d'entree (attention: l'indice des niveaux
    ! verticaux augmente du haut vers le bas)

    DO k = 1, klev
       DO i = 1, klon
          tm1(i, k) = t_seri(i, klev - k + 1)
          um1(i, k) = u_seri(i, klev - k + 1)
          vm1(i, k) = v_seri(i, klev - k + 1)
          papm1(i, k) = play(i, klev - k + 1)
       END DO
    END DO

    DO k = 1, klev + 1
       DO i = 1, klon
          paphm1(i, k) = paprs(i, klev - k + 2)
       END DO
    END DO

    DO i = 1, klon
       zgeom(i, klev) = rd * tm1(i, klev) &
            * log(paphm1(i, klev + 1) / papm1(i, klev))
    END DO

    DO k = klev - 1, 1, - 1
       DO i = 1, klon
          zgeom(i, k) = zgeom(i, k + 1) + rd * (tm1(i, k) + tm1(i, k + 1)) &
               / 2. * log(papm1(i, k + 1) / papm1(i, k))
       END DO
    END DO

    ! Appeler la routine principale
    CALL orodrag(paphm1, papm1, zgeom, tm1, um1, vm1, zmea, zstd, zsig, zgam, &
         zthe, zpic, zval, pdudt, pdvdt, pdtdt)

    ustrdr = 0.
    vstrdr = 0.

    DO k = 1, klev
       DO i = 1, klon
          d_u_oro(i, klev + 1 - k) = dtphys * pdudt(i, k)
          d_v_oro(i, klev + 1 - k) = dtphys * pdvdt(i, k)
          d_t_oro(i, klev + 1 - k) = dtphys * pdtdt(i, k)
          ustrdr(i) = ustrdr(i) &
               + pdudt(i, k) * (paphm1(i, k + 1) - paphm1(i, k)) / rg
          vstrdr(i) = vstrdr(i) &
               + pdvdt(i, k) * (paphm1(i, k + 1) - paphm1(i, k)) / rg
       END DO
    END DO

  END SUBROUTINE drag_noro

end module drag_noro_m
