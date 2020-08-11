module cltracrn_m

  IMPLICIT none

contains 

  SUBROUTINE cltracrn(itr, dtime, u1lay, v1lay, coef, cdragh, t, ftsol, &
       pctsrf, tr, trs, paprs, pplay, delp, masktr, fshtr, hsoltr, tautr, &
       vdeptr, lat, d_tr, d_trs)

    ! From phylmd/cltracrn.F, version 1.2 2005/05/25 13:10:09

    ! Author: Alexandre ARMENGAUD
    ! Date: February 1999
    ! Inspiré de clqh et clvent

    ! Objet: diffusion verticale de traceurs avec quantité de traceur
    ! dans le sol (réservoir de sol de radon)

    ! Note : pour l'instant le traceur dans le sol et le flux sont
    ! calculés mais ils ne servent que de diagnostics. Seule la
    ! tendance sur le traceur est sortie (d_tr).

    use indicesol, only: nbsrf
    use dimphy, only: klon, klev
    use SUPHEC_M, only: RD, rg

    INTEGER itr
    ! itr--- -input-R- le type de traceur 1- Rn 2 - Pb

    REAL, intent(in):: dtime
    ! dtime----input-R- intervalle de temps (en second)
    REAL, intent(in):: u1lay(klon), v1lay(klon) ! vent de la premiere
                                                ! couche (m/s)

    REAL, intent(in):: coef(:, 2:) ! (klon, 2:klev)
    ! coefficient d'echange (m^2 / s)

    real, intent(in):: cdragh(:) ! klon
    REAL, intent(in):: t(klon, klev) ! temperature (K)
    real, intent(in):: ftsol(klon, nbsrf), pctsrf(klon, nbsrf) 
    ! ftsol----input-R- temperature du sol (en Kelvin)
    REAL, intent(in):: tr(klon, klev) ! traceur
    REAL, intent(in):: trs(:) ! (klon) traceur dans le sol
    REAL, intent(in):: paprs(klon, klev+1)
    ! paprs----input-R- pression a l'inter-couche (Pa)
    real, intent(in):: pplay(klon, klev)
    ! pplay----input-R- pression au milieu de couche (Pa)
    real delp(klon, klev)
    ! delp-----input-R- epaisseur de couche (Pa)
    REAL masktr(klon) 
    ! masktr---input-R- Masque reservoir de sol traceur (1 = reservoir)
    REAL fshtr(klon) 
    ! fshtr----input-R- Flux surfacique de production dans le sol
    REAL hsoltr
    ! hsoltr---input-R- Epaisseur equivalente du reservoir de sol
    REAL tautr
    ! tautr----input-R- Constante de decroissance du traceur

    REAL, intent(in):: vdeptr
    ! vitesse de d\'ep\^ot sec dans la couche brownienne

    REAL, intent(in):: lat(klon) 
    ! lat-----input-R- latitude en degree
    REAL d_tr(klon, klev)
    ! d_tr-----output-R- le changement de "tr"

    REAL, intent(out):: d_trs(:) ! (klon) (diagnostic) changement de "trs"

    ! Local:
    INTEGER i, k, n, l
    REAL rotrhi(klon)
    REAL zx_coef(klon, klev)
    REAL zx_buf(klon)
    REAL zx_ctr(klon, klev)
    REAL zx_dtr(klon, klev)
    REAL zx_a, zx_b

    REAL local_tr(klon, klev)
    REAL local_trs(klon)
    REAL zts(klon)
    REAL zx_alpha1(klon), zx_alpha2(klon)

    !------------------------------------------------------------------

    ! Pour l'instant les quatre types de surface ne sont pas pris en
    ! compte. On fabrique avec zts un champ de température de sol que
    ! l'on pondère par la fraction de sol.

    DO i = 1, klon
       zts(i) = 0. 
    ENDDO

    DO n=1, nbsrf
       DO i = 1, klon
          zts(i) = zts(i) + ftsol(i, n)*pctsrf(i, n)
       ENDDO
    ENDDO

    DO i = 1, klon
       rotrhi(i) = RD * zts(i) / hsoltr 
    END DO

    DO k = 1, klev
       DO i = 1, klon
          local_tr(i, k) = tr(i, k)
       ENDDO
    ENDDO

    local_trs = trs

    ! Attention si dans pbl_surface zx_alf1(i) = 1.
    ! Il doit y avoir coherence (donc la meme chose ici)

    DO i = 1, klon
       zx_alpha1(i) = 1.0
       zx_alpha2(i) = 1.0 - zx_alpha1(i)
    ENDDO

    DO i = 1, klon
       zx_coef(i, 1) = cdragh(i) &
            * (1.0+SQRT(u1lay(i)**2+v1lay(i)**2)) &
            * pplay(i, 1)/(RD*t(i, 1))
       zx_coef(i, 1) = zx_coef(i, 1) * dtime*RG
    ENDDO

    DO k = 2, klev
       DO i = 1, klon
          zx_coef(i, k) = coef(i, k)*RG/(pplay(i, k-1)-pplay(i, k)) &
               *(paprs(i, k)*2/(t(i, k)+t(i, k-1))/RD)**2
          zx_coef(i, k) = zx_coef(i, k) * dtime*RG
       ENDDO
    ENDDO

    DO i = 1, klon
       zx_buf(i) = delp(i, klev) + zx_coef(i, klev)
       zx_ctr(i, klev) = local_tr(i, klev)*delp(i, klev)/zx_buf(i)
       zx_dtr(i, klev) = zx_coef(i, klev) / zx_buf(i)
    ENDDO

    DO l = klev-1, 2 , -1
       DO i = 1, klon
          zx_buf(i) = delp(i, l)+zx_coef(i, l) &
               +zx_coef(i, l+1)*(1.-zx_dtr(i, l+1))
          zx_ctr(i, l) = (local_tr(i, l)*delp(i, l) &
               + zx_coef(i, l+1)*zx_ctr(i, l+1))/zx_buf(i)
          zx_dtr(i, l) = zx_coef(i, l) / zx_buf(i)
       ENDDO
    ENDDO

    DO i = 1, klon
       zx_buf(i) = delp(i, 1) + zx_coef(i, 2)*(1.-zx_dtr(i, 2)) &
            + masktr(i) * zx_coef(i, 1) &
            *(zx_alpha1(i)+zx_alpha2(i)*zx_dtr(i, 2))
       zx_ctr(i, 1) = (local_tr(i, 1)*delp(i, 1) &
            + zx_ctr(i, 2) &
            *(zx_coef(i, 2) &
            - masktr(i) * zx_coef(i, 1) &
            *zx_alpha2(i))) / zx_buf(i)
       zx_dtr(i, 1) = masktr(i) * zx_coef(i, 1) / zx_buf(i)
    ENDDO

    ! Calculer d'abord local_trs nouvelle quantite dans le reservoir
    ! de sol

    ! Au dessus des continents
    ! Le pb peut se deposer partout : vdeptr = 10-3 m/s
    ! Le Rn est traité commme une couche Brownienne puisque vdeptr = 0.

    DO i = 1, klon
       IF (NINT(masktr(i)) .EQ. 1) THEN
          zx_a = local_trs(i) &
               +fshtr(i)*dtime*rotrhi(i) &
               +rotrhi(i)*masktr(i)*zx_coef(i, 1)/RG &
               *(zx_ctr(i, 1)*(zx_alpha1(i)+zx_alpha2(i)*zx_dtr(i, 2)) &
               +zx_alpha2(i)*zx_ctr(i, 2))
          ! Pour l'instant, pour aller vite, le d\'ep\^ot sec est trait\'e
          ! comme une d\'ecroissance :
          zx_b = 1. + rotrhi(i)*masktr(i)*zx_coef(i, 1)/RG &
               * (1.-zx_dtr(i, 1) &
               *(zx_alpha1(i)+zx_alpha2(i)*zx_dtr(i, 2))) &
               + dtime / tautr &
               + dtime * vdeptr / hsoltr
          local_trs(i) = zx_a / zx_b
       ENDIF

       ! Si on est entre 60N et 70N on divise par 2 l'emanation

       IF ((itr.eq.1.AND.NINT(masktr(i)).EQ.1.AND.lat(i).GE.60. &
            .AND.lat(i).LE.70.) .OR. &
            (itr.eq.2.AND.NINT(masktr(i)).EQ.1.AND.lat(i).GE.60. &
            .AND.lat(i).LE.70.)) THEN
          zx_a = local_trs(i) &
               +(fshtr(i)/2.)*dtime*rotrhi(i) &
               +rotrhi(i)*masktr(i)*zx_coef(i, 1)/RG &
               *(zx_ctr(i, 1)*(zx_alpha1(i)+zx_alpha2(i)*zx_dtr(i, 2)) &
               +zx_alpha2(i)*zx_ctr(i, 2))
          zx_b = 1. + rotrhi(i)*masktr(i)*zx_coef(i, 1)/RG &
               * (1.-zx_dtr(i, 1) &
               *(zx_alpha1(i)+zx_alpha2(i)*zx_dtr(i, 2))) &
               + dtime / tautr &
               + dtime * vdeptr / hsoltr
          local_trs(i) = zx_a / zx_b
       ENDIF

       ! Au dessus des oceans et aux hautes latitudes

       ! au dessous de -60S pas d'emission de radon au dessus 
       ! des oceans et des continents

       IF ((itr.EQ.1.AND.NINT(masktr(i)).EQ.0) .OR. &
            (itr.EQ.1.AND.NINT(masktr(i)).EQ.1.AND.lat(i).LT.-60.)) THEN
          local_trs(i) = 0.
       END IF

       ! au dessus de 70 N pas d'emission de radon au dessus 
       ! des oceans et des continents

       IF ((itr.EQ.1.AND.NINT(masktr(i)).EQ.0) .OR. &
            (itr.EQ.1.AND.NINT(masktr(i)).EQ.1.AND.lat(i).GT.70.)) THEN
          local_trs(i) = 0.
       END IF

       ! Au dessus des oceans la source est nulle

       IF (itr.eq.1.AND.NINT(masktr(i)).EQ.0) THEN
          local_trs(i) = 0.
       END IF
    ENDDO ! sur le i=1, klon

    ! une fois qu'on a local_trs, on peut faire l'iteration

    DO i = 1, klon
       local_tr(i, 1) = zx_ctr(i, 1)+zx_dtr(i, 1)*local_trs(i)
    ENDDO
    DO l = 2, klev
       DO i = 1, klon
          local_tr(i, l) &
               = zx_ctr(i, l) + zx_dtr(i, l)*local_tr(i, l-1)
       ENDDO
    ENDDO

    ! Calcul des tendances du traceur dans le sol et dans l'atmosphere

    DO l = 1, klev
       DO i = 1, klon
          d_tr(i, l) = local_tr(i, l) - tr(i, l)
       ENDDO
    ENDDO
    d_trs = local_trs - trs

  END SUBROUTINE cltracrn

end module cltracrn_m
