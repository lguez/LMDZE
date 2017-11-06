module clvent_m

  IMPLICIT none

contains

  SUBROUTINE clvent(dtime, u1lay, v1lay, coef, cdrag, t, ven, paprs, pplay, &
       delp, d_ven, flux_v)

    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1993/08/18
    ! Objet : diffusion verticale de la vitesse

    USE dimphy, ONLY: klev
    use nr_util, only: assert_eq
    USE suphec_m, ONLY: rd, rg

    REAL, intent(in):: dtime ! intervalle de temps (en s)

    REAL, intent(in):: u1lay(:), v1lay(:) ! (knon)
    ! vent de la premiere couche (m/s)

    REAL, intent(in):: coef(:, 2:) ! (knon, 2:klev)
    ! Coefficient d'echange (m**2/s) multiplié par le cisaillement du
    ! vent (dV/dz)

    REAL, intent(in):: cdrag(:) ! (knon) sans unité

    REAL, intent(in):: t(:, :) ! (knon, klev) ! temperature (K)
    REAL, intent(in):: ven(:, :) ! (knon, klev) vitesse horizontale (m/s)
    REAL, intent(in):: paprs(:, :) ! (knon, klev+1) pression a inter-couche (Pa)
    real, intent(in):: pplay(:, :) ! (knon, klev) pression au milieu
                                   ! de couche (Pa)
    real, intent(in):: delp(:, :) ! (knon, klev) epaisseur de couche (Pa)
    REAL, intent(out):: d_ven(:, :) ! (knon, klev) ! le changement de "ven"

    REAL, intent(out):: flux_v(:) ! (knon)
    ! (diagnostic) flux du vent à la surface, en (kg m/s)/(m**2 s)
    ! flux_v est le flux de moment angulaire (positif vers bas)

    ! Local:
    INTEGER knon, i, k
    REAL zx_cv(size(u1lay), 2:klev) ! (knon, 2:klev)
    REAL zx_dv(size(u1lay), 2:klev) ! (knon, 2:klev)
    REAL zx_buf(size(u1lay)) ! (knon)
    REAL zx_coef(size(u1lay), klev) ! (knon, klev)
    REAL local_ven(size(u1lay), klev) ! (knon, klev)

    !------------------------------------------------------------------

    knon = assert_eq([size(u1lay), size(v1lay), size(coef, 1), size(t, 1), &
         size(ven, 1), size(paprs, 1), size(pplay, 1), size(delp, 1), &
         size(d_ven, 1), size(flux_v)], "clvent knon")
    local_ven = ven

    DO i = 1, knon
       zx_coef(i, 1) = cdrag(i) * (1. + SQRT(u1lay(i)**2 + v1lay(i)**2)) &
            * pplay(i, 1) / (RD * t(i, 1)) * dtime * RG
    ENDDO

    DO k = 2, klev
       DO i = 1, knon
          zx_coef(i, k) = coef(i, k) * RG / (pplay(i, k-1) - pplay(i, k)) &
               * (paprs(i, k) * 2 / (t(i, k) + t(i, k - 1)) / RD)**2
          zx_coef(i, k) = zx_coef(i, k) * dtime * RG
       ENDDO
    ENDDO

    DO i = 1, knon
       zx_buf(i) = delp(i, 1) + zx_coef(i, 1)+zx_coef(i, 2)
       zx_cv(i, 2) = local_ven(i, 1)*delp(i, 1) / zx_buf(i)
       zx_dv(i, 2) = zx_coef(i, 2) &
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

    DO i = 1, knon
       flux_v(i) = zx_coef(i, 1)/(RG*dtime) &
            *local_ven(i, 1)
    ENDDO

    DO k = 1, klev
       DO i = 1, knon
          d_ven(i, k) = local_ven(i, k) - ven(i, k)
       ENDDO
    ENDDO

  END SUBROUTINE clvent

end module clvent_m
