module clvent_m

  IMPLICIT none

contains

  SUBROUTINE clvent(knon, dtime, u1lay, v1lay, coef, t, ven, paprs, pplay, &
       delp, d_ven, flux_v)

    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1993/08/18
    ! Objet : diffusion verticale de la vitesse

    USE dimphy, ONLY: klev, klon
    USE suphec_m, ONLY: rd, rg

    INTEGER knon
    REAL, intent(in):: dtime
    ! dtime----input-R- intervalle du temps (en second)

    REAL u1lay(klon), v1lay(klon)
    ! u1lay----input-R- vent u de la premiere couche (m/s)
    ! v1lay----input-R- vent v de la premiere couche (m/s)

    REAL, intent(in):: coef(:, :) ! (knon, klev)
    ! Coefficient d'echange (m**2/s) multiplié par le cisaillement du
    ! vent (dV/dz). La première valeur indique la valeur de Cdrag (sans
    ! unité).

    REAL t(klon, klev), ven(klon, klev)
    ! t--------input-R- temperature (K)
    ! ven------input-R- vitesse horizontale (m/s)
    REAL paprs(klon, klev+1), pplay(klon, klev), delp(klon, klev)
    ! paprs----input-R- pression a inter-couche (Pa)
    ! pplay----input-R- pression au milieu de couche (Pa)
    ! delp-----input-R- epaisseur de couche (Pa)
    REAL d_ven(klon, klev)
    ! d_ven----output-R- le changement de "ven"
    REAL flux_v(klon, klev)
    ! flux_v---output-R- (diagnostic) flux du vent: (kg m/s)/(m**2 s)

    ! Local:
    INTEGER i, k
    REAL zx_cv(klon, 2:klev)
    REAL zx_dv(klon, 2:klev)
    REAL zx_buf(klon)
    REAL zx_coef(klon, klev)
    REAL local_ven(klon, klev)
    REAL zx_alf1(klon), zx_alf2(klon)

    !------------------------------------------------------------------

    DO k = 1, klev
       DO i = 1, knon
          local_ven(i, k) = ven(i, k)
       ENDDO
    ENDDO

    DO i = 1, knon
       zx_alf1(i) = 1.0
       zx_alf2(i) = 1.0 - zx_alf1(i)
       zx_coef(i, 1) = coef(i, 1) * (1. + SQRT(u1lay(i)**2 + v1lay(i)**2)) &
            * pplay(i, 1) / (RD * t(i, 1))
       zx_coef(i, 1) = zx_coef(i, 1) * dtime * RG
    ENDDO

    DO k = 2, klev
       DO i = 1, knon
          zx_coef(i, k) = coef(i, k) * RG / (pplay(i, k-1) - pplay(i, k)) &
               * (paprs(i, k) * 2 / (t(i, k) + t(i, k - 1)) / RD)**2
          zx_coef(i, k) = zx_coef(i, k) * dtime * RG
       ENDDO
    ENDDO

    DO i = 1, knon
       zx_buf(i) = delp(i, 1) + zx_coef(i, 1)*zx_alf1(i)+zx_coef(i, 2)
       zx_cv(i, 2) = local_ven(i, 1)*delp(i, 1) / zx_buf(i)
       zx_dv(i, 2) = (zx_coef(i, 2)-zx_alf2(i)*zx_coef(i, 1)) &
            /zx_buf(i)
    ENDDO
    DO k = 3, klev
       DO i = 1, knon
          zx_buf(i) = delp(i, k-1) + zx_coef(i, k) &
               + zx_coef(i, k-1)*(1.-zx_dv(i, k-1))
          zx_cv(i, k) = (local_ven(i, k-1)*delp(i, k-1) &
               +zx_coef(i, k-1)*zx_cv(i, k-1) )/zx_buf(i)
          zx_dv(i, k) = zx_coef(i, k)/zx_buf(i)
       ENDDO
    ENDDO
    DO i = 1, knon
       local_ven(i, klev) = ( local_ven(i, klev)*delp(i, klev) &
            +zx_coef(i, klev)*zx_cv(i, klev) ) &
            / ( delp(i, klev) + zx_coef(i, klev) &
            -zx_coef(i, klev)*zx_dv(i, klev) )
    ENDDO
    DO k = klev-1, 1, -1
       DO i = 1, knon
          local_ven(i, k) = zx_cv(i, k+1) + zx_dv(i, k+1)*local_ven(i, k+1)
       ENDDO
    ENDDO

    ! flux_v est le flux de moment angulaire (positif vers bas) dont
    ! l'unite est: (kg m/s)/(m**2 s)
    DO i = 1, knon
       flux_v(i, 1) = zx_coef(i, 1)/(RG*dtime) &
            *(local_ven(i, 1)*zx_alf1(i) &
            +local_ven(i, 2)*zx_alf2(i))
    ENDDO
    DO k = 2, klev
       DO i = 1, knon
          flux_v(i, k) = zx_coef(i, k)/(RG*dtime) &
               * (local_ven(i, k)-local_ven(i, k-1))
       ENDDO
    ENDDO

    DO k = 1, klev
       DO i = 1, knon
          d_ven(i, k) = local_ven(i, k) - ven(i, k)
       ENDDO
    ENDDO

  END SUBROUTINE clvent

end module clvent_m
