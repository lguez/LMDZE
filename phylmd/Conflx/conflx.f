module conflx_m

  IMPLICIT none

contains

  SUBROUTINE conflx (dtime, pres_h, pres_f, t, q, con_t, con_q, qhfl, w, &
       d_t, d_q, rain, snow, mfu, mfd, pen_u, pde_u, pen_d, pde_d, kcbot, &
       kctop, kdtop, pmflxr, pmflxs)

    ! From LMDZ4/libf/phylmd/conflx.F, version 1.1.1.1 2004/05/19 12:53:08

    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1994/10/14

    ! Objet: schéma en flux de masse pour la convection (schéma de
    ! Tiedtke avec quelques modifications mineures)

    ! Décembre 1997 : prise en compte des modifications introduites
    ! par Olivier Boucher et Alexandre Armengaud pour le mélange et le
    ! lessivage des traceurs passifs.

    use flxmain_m, only: flxmain
    USE dimphy, ONLY: klev, klon
    USE suphec_m, ONLY: rd, retv, rtt
    USE yoethf_m, ONLY: r2es
    USE fcttre, ONLY: foeew

    REAL, intent(in):: dtime ! pas d'integration (s)
    REAL, intent(in):: pres_h(:, :) ! (klon, klev + 1) pression half-level (Pa)
    REAL, intent(in):: pres_f(:, :) ! (klon, klev) pression full-level (Pa)
    REAL, intent(in):: t(:, :) ! (klon, klev) temperature (K)
    REAL, intent(in):: q(:, :) ! (klon, klev) humidité spécifique (no dimension)

    REAL, intent(in):: con_t(:, :)
    ! (klon, klev) convergence de temperature (K/s)

    REAL, intent(in):: con_q(:, :) 
    ! (klon, klev) convergence de l'eau vapeur (g/g/s)

    REAL, intent(in):: qhfl(:) ! (klon) evaporation (negative vers haut) mm/s
    REAL, intent(in):: w(:, :) ! (klon, klev) vitesse verticale (Pa/s)

    REAL, intent(out):: d_t(:, :) ! (klon, klev) incrementation de temperature
    REAL, intent(out):: d_q(:, :) ! (klon, klev) incrementation d'humidite
    REAL, intent(out):: rain(:) ! (klon) pluie (mm/s)
    REAL, intent(out):: snow(:) ! (klon) neige (mm/s)

    REAL, intent(out):: mfu(:, :) ! (klon, klev)
    ! flux de masse (kg/m2/s) panache ascendant

    REAL, intent(out):: mfd(:, :) ! (klon, klev)
    ! flux de masse (kg/m2/s) panache descendant

    REAL, intent(out):: pen_u(:, :) ! (klon, klev)
    REAL, intent(out):: pde_u(:, :) ! (klon, klev)
    REAL, intent(out):: pen_d(:, :) ! (klon, klev)
    REAL, intent(out):: pde_d(:, :) ! (klon, klev)
    INTEGER, intent(out):: kcbot(:) ! (klon) niveau du bas de la convection
    INTEGER, intent(out):: kctop(:) ! (klon) niveau du haut de la convection
    INTEGER, intent(out):: kdtop(:) ! (klon) niveau du haut des downdrafts
    REAL, intent(out):: pmflxr(:, :) ! (klon, klev + 1)
    REAL, intent(out):: pmflxs(:, :) ! (klon, klev + 1)

    ! Local:

    REAL qsen(klon, klev)
    REAL pvervel(klon, klev)
    LOGICAL land(klon)

    REAL d_t_bis(klon, klev)
    REAL d_q_bis(klon, klev)
    REAL paprs(klon, klev + 1)
    REAL paprsf(klon, klev)
    REAL zgeom(klon, klev)
    REAL zcvgq(klon, klev)
    REAL zcvgt(klon, klev)

    REAL zen_u(klon, klev)
    REAL zen_d(klon, klev)
    REAL zde_u(klon, klev)
    REAL zde_d(klon, klev)
    REAL zmflxr(klon, klev + 1)
    REAL zmflxs(klon, klev + 1)

    INTEGER i, k
    REAL zqsat

    !--------------------------------------------------------------------

    ! Initialiser les variables de sortie (pour securité):
    rain = 0.
    snow = 0.
    kcbot = 0
    kctop = 0
    kdtop = 0
    d_t = 0.
    d_q = 0.

    zen_u = 0.
    zde_u = 0.
    zen_d = 0.
    zde_d = 0.
    zmflxr = 0.
    zmflxs = 0.

    ! Calculer la nature du sol (pour l'instant, océan partout):
    land = .FALSE.

    ! Préparer les variables d'entrée (attention: l'indice des niveaux
    ! verticaux augmente du haut vers le bas) :
    DO k = 1, klev
       DO i = 1, klon
          paprsf(i, k) = pres_f(i, klev-k + 1)
          paprs(i, k) = pres_h(i, klev + 1-k + 1)
          pvervel(i, k) = w(i, klev + 1-k)
          zcvgt(i, k) = con_t(i, klev-k + 1)
          zcvgq(i, k) = con_q(i, klev-k + 1)

          zqsat = MIN(0.5, R2ES * FOEEW(t(i, k), &
               merge(0., 1., rtt < t(i, k))) / paprsf(i, k))
          qsen(i, k) = zqsat / (1. - RETV * zqsat)
       ENDDO
    ENDDO
    DO i = 1, klon
       paprs(i, klev + 1) = pres_h(i, 1)
       zgeom(i, klev) = RD * t(i, klev) &
            / (0.5*(paprs(i, klev + 1) + paprsf(i, klev))) &
            * (paprs(i, klev + 1)-paprsf(i, klev))
    ENDDO
    DO k = klev-1, 1, -1
       DO i = 1, klon
          zgeom(i, k) = zgeom(i, k + 1) &
               + RD * 0.5*(t(i, k + 1) + t(i, k)) / paprs(i, k + 1) &
               * (paprsf(i, k + 1)-paprsf(i, k))
       ENDDO
    ENDDO

    ! Appeler la routine principale :
    CALL flxmain(dtime, t, q, qsen, qhfl, paprsf, paprs, zgeom, land, &
         zcvgt, zcvgq, pvervel, rain, snow, kcbot, kctop, kdtop, mfu, mfd, &
         zen_u, zde_u, zen_d, zde_d, d_t_bis, d_q_bis, zmflxr, zmflxs)

    ! De la même façon que l'on effectue le réindiçage pour la
    ! température t et le champ q, on réindice les flux nécessaires à
    ! la convection des traceurs.
    DO k = 1, klev
       DO i = 1, klon
          d_q(i, klev + 1-k) = dtime*d_q_bis(i, k)
          d_t(i, klev + 1-k) = dtime*d_t_bis(i, k)
       ENDDO
    ENDDO

    mfu = eoshift(mfu, shift=1, dim=2)
    mfd = eoshift(mfd, shift=1, dim=2)
    pen_d(:, 1)= 0.
    pde_d(:, 1)= 0.

    DO k = 1, klev
       DO i = 1, klon
          pen_u(i, klev + 1-k)= zen_u(i, k)
          pde_u(i, klev + 1-k)= zde_u(i, k)
       ENDDO
    ENDDO

    DO k = 1, klev-1
       DO i = 1, klon
          pen_d(i, klev + 1-k)= -zen_d(i, k + 1)
          pde_d(i, klev + 1-k)= -zde_d(i, k + 1)
       ENDDO
    ENDDO

    DO k = 1, klev + 1
       DO i = 1, klon
          pmflxr(i, klev + 2-k)= zmflxr(i, k)
          pmflxs(i, klev + 2-k)= zmflxs(i, k)
       ENDDO
    ENDDO

  END SUBROUTINE conflx

end module conflx_m
