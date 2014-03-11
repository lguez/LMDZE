module pvtheta_m

  IMPLICIT NONE

contains

  SUBROUTINE pvtheta(ilon, ilev, pucov, pvcov, pteta, ztfi, zplay, zplev, &
       nbteta, theta, pvteta)

    USE comconst, ONLY: cpp, g, r
    USE comgeom, ONLY: aire_2d, alpha1_2d, alpha2_2d, alpha3_2d, alpha4_2d, &
         apoln, apols
    USE dimens_m, ONLY: iim, jjm, llm
    USE disvert_m, ONLY: preff
    USE paramet_m, ONLY: iip1, ip1jm, ip1jmp1, jjp1
    USE tourabs_m, ONLY: tourabs

    ! Author: I. Musat

    ! Calcul de la vorticite potentielle PVteta sur des iso-theta selon
    ! la methodologie du NCEP/NCAR :

    ! 1) on calcule n2

    ! 2) on interpole les vents, la temperature et le N2 sur des
    ! isentropes (en fait sur des iso-theta) lineairement en
    ! log(theta), d'où ucovteta, vcovteta, N2teta

    ! 3) on calcule vorateta
    ! 4) on calcule rhoteta
    ! 5) on calcule PVteta

    INTEGER, INTENT(IN):: ilon
    INTEGER, INTENT(IN):: ilev

    ! pucov, pvcov, pteta, ztfi, zplay, zplev sur la grille dynamique
    REAL, INTENT(IN):: pucov(iip1, jjp1, ilev)
    REAL, INTENT(IN):: pvcov(iip1, jjm, ilev)
    REAL, INTENT(IN):: pteta(iip1, jjp1, ilev)
    REAL, INTENT(IN):: ztfi(ilon, ilev)
    REAL, INTENT(IN):: zplay(ilon, ilev), zplev(ilon, ilev + 1)

    INTEGER, INTENT(IN):: nbteta
    REAL, INTENT(IN):: theta(nbteta) ! sur la grille dynamique

    REAL, INTENT(out):: pvteta(ilon, nbteta) 
    ! vorticité potentielle sur des iso-theta, sur la grille physique,
    ! en Pa-1 s-1
    ! PVteta = vorateta * N2 / (g**2 * rhoteta)

    ! Local:

    INTEGER i, j, l, ig0
    REAL ssum
    REAL teta(ilon, ilev)
    REAL ptetau(ip1jmp1, ilev), ptetav(ip1jm, ilev)
    REAL ucovteta(ip1jmp1, ilev), vcovteta(ip1jm, ilev)

    REAL n2(ilon, ilev - 1) ! stabilité statique sur les niveaux du modèle
    ! N**2 = g / T * (dT / dz + g / cp)

    real n2teta(ilon, nbteta) ! N**2 sur une iso-theta
    REAL ztfiteta(ilon, nbteta)

    REAL rhoteta(ilon, nbteta) ! densite sur des iso-theta
    ! rhoteta = (T / theta)**(cp / R) * p0 / (R * T)

    REAL vorateta(iip1, jjm, nbteta) ! vorticite absolue sur des iso-theta
    REAL voratetafi(ilon, nbteta), vorpol(iim)

    !-------------------------------------------------------------------------

    ! projection teta sur la grille physique
    DO l = 1, llm
       teta(1, l) = pteta(1, 1, l)
       ig0 = 2
       DO j = 2, jjm
          DO i = 1, iim
             teta(ig0, l) = pteta(i, j, l)
             ig0 = ig0 + 1
          END DO
       END DO
       teta(ig0, l) = pteta(1, jjp1, l)
    END DO

    ! calcul pteta sur les grilles U et V
    DO l = 1, llm
       DO j = 1, jjp1
          DO i = 1, iip1
             ig0 = i + (j - 1) * iip1
             ptetau(ig0, l) = pteta(i, j, l)
          END DO
       END DO
       DO j = 1, jjm
          DO i = 1, iip1
             ig0 = i + (j - 1) * iip1
             ptetav(ig0, l) = 0.5 * (pteta(i, j, l) + pteta(i, j + 1, l))
          END DO
       END DO
    END DO

    ! projection pucov, pvcov sur une surface de theta constante
    DO l = 1, nbteta
       CALL tetalevel(ip1jmp1, llm, .TRUE., ptetau, theta(l), pucov, &
            ucovteta(:, l))
       CALL tetalevel(ip1jm, llm, .TRUE., ptetav, theta(l), pvcov, &
            vcovteta(:, l))
    END DO

    CALL tourabs(nbteta, vcovteta, ucovteta, vorateta)

    ! projection vorateta sur la grille physique => voratetafi

    DO l = 1, nbteta
       DO j = 2, jjm
          ig0 = 1 + (j - 2) * iim
          DO i = 1, iim
             voratetafi(ig0 + i + 1, l) = vorateta(i, j - 1, l) &
                  * alpha4_2d(i + 1, j) + vorateta(i + 1, j - 1, l) &
                  * alpha1_2d(i + 1, j) + vorateta(i, j, l) &
                  * alpha3_2d(i + 1, j) + vorateta(i + 1, j, l) &
                  * alpha2_2d(i + 1, j)
          END DO
          voratetafi(ig0 + 1, l) = voratetafi(ig0 + 1 + iim, l)
       END DO
    END DO

    DO l = 1, nbteta
       DO i = 1, iim
          vorpol(i) = vorateta(i, 1, l) * aire_2d(i, 1)
       END DO
       voratetafi(1, l) = ssum(iim, vorpol, 1) / apoln
    END DO

    DO l = 1, nbteta
       DO i = 1, iim
          vorpol(i) = vorateta(i, jjm, l) * aire_2d(i, jjm + 1)
       END DO
       voratetafi(ilon, l) = ssum(iim, vorpol, 1) / apols
    END DO

    DO l = 1, llm - 1
       DO i = 1, ilon
          n2(i, l) = (g**2 * zplay(i, l) * (ztfi(i, l + 1) - ztfi(i, l))) &
               / (r * ztfi(i, l) * ztfi(i, l) &
               * (zplev(i, l) - zplev(i, l + 1))) + (g**2) / (ztfi(i, l) * cpp)
       END DO
    END DO

    ! calcul N2teta
    DO l = 1, nbteta
       CALL tetalevel(ilon, llm - 1, .TRUE., teta, theta(l), n2, n2teta(:, l))
       CALL tetalevel(ilon, llm, .TRUE., teta, theta(l), ztfi, ztfiteta(:, l))
    END DO

    DO l = 1, nbteta
       DO i = 1, ilon
          rhoteta(i, l) = (ztfiteta(i, l) / theta(l))**(cpp / r) &
               * (preff / (r * ztfiteta(i, l)))
          pvteta(i, l) = (voratetafi(i, l) * n2teta(i, l)) &
               / (g**2 * rhoteta(i, l))
       END DO
    END DO

  END SUBROUTINE pvtheta

end module pvtheta_m
