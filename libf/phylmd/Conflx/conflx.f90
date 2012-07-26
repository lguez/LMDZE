module conflx_m

  IMPLICIT none

contains

  SUBROUTINE conflx (dtime, pres_h, pres_f, t, q, con_t, con_q, pqhfl, w, &
       d_t, d_q, rain, snow, pmfu, pmfd, pen_u, pde_u, pen_d, pde_d, kcbot, &
       kctop, kdtop, pmflxr, pmflxs)

    ! From LMDZ4/libf/phylmd/conflx.F, version 1.1.1.1 2004/05/19 12:53:08

    ! Author: Z. X. Li (LMD/CNRS)
    ! date: 1994/10/14

    ! Objet: schéma flux de masse pour la convection (schéma de
    ! Tiedtke avec quelques modifications mineures)

    ! Décembre 1997 : prise en compte des modifications introduites
    ! par Olivier Boucher et Alexandre Armengaud pour mélange et
    ! lessivage des traceurs passifs.

    use flxmain_m, only: flxmain
    USE dimphy, ONLY: klev, klon
    USE suphec_m, ONLY: rd, retv, rtt
    USE yoethf_m, ONLY: r2es
    USE fcttre, ONLY: foeew

    ! Entree:
    REAL, intent(in):: dtime            ! pas d'integration (s)
    REAL, intent(in):: pres_h(klon, klev+1) ! pression half-level (Pa)
    REAL, intent(in):: pres_f(klon, klev)! pression full-level (Pa)
    REAL, intent(in):: t(klon, klev)     ! temperature (K)
    REAL q(klon, klev)     ! humidite specifique (g/g)
    REAL w(klon, klev)     ! vitesse verticale (Pa/s)
    REAL con_t(klon, klev) ! convergence de temperature (K/s)
    REAL con_q(klon, klev) ! convergence de l'eau vapeur (g/g/s)
    REAL pqhfl(klon)      ! evaporation (negative vers haut) mm/s

    ! Sortie:
    REAL d_t(klon, klev)   ! incrementation de temperature
    REAL d_q(klon, klev)   ! incrementation d'humidite

    REAL, intent(out):: pmfu(:, :) ! (klon, klev)
    ! flux masse (kg/m2/s) panache ascendant
    
    REAL, intent(out):: pmfd(:, :) ! (klon, klev)
    ! flux masse (kg/m2/s) panache descendant

    REAL pen_u(klon, klev)
    REAL pen_d(klon, klev)
    REAL pde_u(klon, klev)
    REAL pde_d(klon, klev)
    REAL rain(klon)       ! pluie (mm/s)
    REAL snow(klon)       ! neige (mm/s)
    REAL pmflxr(klon, klev+1)
    REAL pmflxs(klon, klev+1)
    INTEGER kcbot(klon)  ! niveau du bas de la convection
    INTEGER kctop(klon)  ! niveau du haut de la convection
    INTEGER kdtop(klon)  ! niveau du haut des downdrafts

    ! Local:

    REAL pt(klon, klev)
    REAL pq(klon, klev)
    REAL pqs(klon, klev)
    REAL pvervel(klon, klev)
    LOGICAL land(klon)

    REAL d_t_bis(klon, klev)
    REAL d_q_bis(klon, klev)
    REAL paprs(klon, klev+1)
    REAL paprsf(klon, klev)
    REAL zgeom(klon, klev)
    REAL zcvgq(klon, klev)
    REAL zcvgt(klon, klev)

    REAL zmfu(klon, klev)
    REAL zmfd(klon, klev)
    REAL zen_u(klon, klev)
    REAL zen_d(klon, klev)
    REAL zde_u(klon, klev)
    REAL zde_d(klon, klev)
    REAL zmflxr(klon, klev+1)
    REAL zmflxs(klon, klev+1)

    INTEGER i, k
    REAL zdelta, zqsat

    !--------------------------------------------------------------------

    ! initialiser les variables de sortie (pour securite)
    DO i = 1, klon
       rain(i) = 0.0
       snow(i) = 0.0
       kcbot(i) = 0
       kctop(i) = 0
       kdtop(i) = 0
    ENDDO
    DO k = 1, klev
       DO i = 1, klon
          d_t(i, k) = 0.0
          d_q(i, k) = 0.0
          pmfu(i, k) = 0.0
          pmfd(i, k) = 0.0
          pen_u(i, k) = 0.0
          pde_u(i, k) = 0.0
          pen_d(i, k) = 0.0
          pde_d(i, k) = 0.0
          zmfu(i, k) = 0.0
          zmfd(i, k) = 0.0
          zen_u(i, k) = 0.0
          zde_u(i, k) = 0.0
          zen_d(i, k) = 0.0
          zde_d(i, k) = 0.0
       ENDDO
    ENDDO
    DO k = 1, klev+1
       DO i = 1, klon
          zmflxr(i, k) = 0.0
          zmflxs(i, k) = 0.0
       ENDDO
    ENDDO

    ! calculer la nature du sol (pour l'instant, ocean partout)
    DO i = 1, klon
       land(i) = .FALSE.
    ENDDO

    ! preparer les variables d'entree (attention: l'ordre des niveaux
    ! verticaux augmente du haut vers le bas)
    DO k = 1, klev
       DO i = 1, klon
          pt(i, k) = t(i, klev-k+1)
          pq(i, k) = q(i, klev-k+1)
          paprsf(i, k) = pres_f(i, klev-k+1)
          paprs(i, k) = pres_h(i, klev+1-k+1)
          pvervel(i, k) = w(i, klev+1-k)
          zcvgt(i, k) = con_t(i, klev-k+1)
          zcvgq(i, k) = con_q(i, klev-k+1)

          zdelta=MAX(0., SIGN(1., RTT-pt(i, k)))
          zqsat=R2ES*FOEEW ( pt(i, k), zdelta ) / paprsf(i, k)
          zqsat=MIN(0.5, zqsat)
          zqsat=zqsat/(1.-RETV  *zqsat)
          pqs(i, k) = zqsat
       ENDDO
    ENDDO
    DO i = 1, klon
       paprs(i, klev+1) = pres_h(i, 1)
       zgeom(i, klev) = RD * pt(i, klev) &
            / (0.5*(paprs(i, klev+1)+paprsf(i, klev))) &
            * (paprs(i, klev+1)-paprsf(i, klev))
    ENDDO
    DO k = klev-1, 1, -1
       DO i = 1, klon
          zgeom(i, k) = zgeom(i, k+1) &
               + RD * 0.5*(pt(i, k+1)+pt(i, k)) / paprs(i, k+1) &
               * (paprsf(i, k+1)-paprsf(i, k))
       ENDDO
    ENDDO

    ! appeler la routine principale

    CALL flxmain(dtime, pt, pq, pqs, pqhfl, paprsf, paprs, zgeom, land, &
         zcvgt, zcvgq, pvervel, rain, snow, kcbot, kctop, kdtop, zmfu, zmfd, &
         zen_u, zde_u, zen_d, zde_d, d_t_bis, d_q_bis, zmflxr, zmflxs)

    ! De la même façon que l'on effectue le réindiçage pour la
    ! température t et le champ q, on réindice les flux nécessaires à
    ! la convection des traceurs.
    DO k = 1, klev
       DO i = 1, klon
          d_q(i, klev+1-k) = dtime*d_q_bis(i, k)
          d_t(i, klev+1-k) = dtime*d_t_bis(i, k)
       ENDDO
    ENDDO

    DO i = 1, klon
       pmfu(i, 1)= 0.
       pmfd(i, 1)= 0.
       pen_d(i, 1)= 0.
       pde_d(i, 1)= 0.
    ENDDO

    DO k = 2, klev
       DO i = 1, klon
          pmfu(i, klev+2-k)= zmfu(i, k)
          pmfd(i, klev+2-k)= zmfd(i, k)
       ENDDO
    ENDDO

    DO k = 1, klev
       DO i = 1, klon
          pen_u(i, klev+1-k)=  zen_u(i, k)
          pde_u(i, klev+1-k)=  zde_u(i, k)
       ENDDO
    ENDDO

    DO k = 1, klev-1
       DO i = 1, klon
          pen_d(i, klev+1-k)= -zen_d(i, k+1)
          pde_d(i, klev+1-k)= -zde_d(i, k+1)
       ENDDO
    ENDDO

    DO k = 1, klev+1
       DO i = 1, klon
          pmflxr(i, klev+2-k)= zmflxr(i, k)
          pmflxs(i, klev+2-k)= zmflxs(i, k)
       ENDDO
    ENDDO

  END SUBROUTINE conflx

end module conflx_m
