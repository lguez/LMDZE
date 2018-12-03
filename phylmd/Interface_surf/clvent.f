module clvent_m

  IMPLICIT none

contains

  SUBROUTINE clvent(u1lay, v1lay, coef, cdrag, t, ven, paprs, pplay, delp, &
       d_ven, flux_v)

    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1993/08/18
    ! Objet : diffusion verticale de la vitesse

    ! Library:
    use nr_util, only: assert

    use comconst, only: dtphys
    USE dimphy, ONLY: klev
    USE suphec_m, ONLY: rd, rg

    REAL, intent(in):: u1lay(:), v1lay(:) ! (knon)
    ! vent de la premiere couche (m / s)

    REAL, intent(in):: coef(:, 2:) ! (knon, 2:klev)
    ! Coefficient d'echange (m**2 / s) multiplié par le cisaillement du
    ! vent (dV / dz)

    REAL, intent(in):: cdrag(:) ! (knon) sans dimension
    REAL, intent(in):: t(:, :) ! (knon, klev) ! temperature (K)
    REAL, intent(in):: ven(:, :) ! (knon, klev) vitesse horizontale (m / s)
    REAL, intent(in):: paprs(:, :) ! (knon, klev + 1) pression a
    ! inter-couche (Pa)
    real, intent(in):: pplay(:, :) ! (knon, klev) pression au milieu
    ! de couche (Pa)
    real, intent(in):: delp(:, :) ! (knon, klev) epaisseur de couche (Pa)

    REAL, intent(out):: d_ven(:, :) ! (knon, klev)
    ! le changement de "ven", en m s-1

    REAL, intent(out):: flux_v(:) ! (knon)
    ! (diagnostic) stress du vent à la surface, en Pa
    ! flux_v est le flux de quantité de mouvement (positif vers bas)

    ! Local:
    INTEGER k
    REAL zx_cv(size(u1lay), 2:klev) ! (knon, 2:klev)
    REAL zx_dv(size(u1lay), 2:klev) ! (knon, 2:klev)
    REAL zx_buf(size(u1lay)) ! (knon)
    REAL zx_coef(size(u1lay), klev) ! (knon, klev) in Pa
    REAL local_ven(size(u1lay), klev) ! (knon, klev) in m s-1

    !------------------------------------------------------------------

    call assert(size(u1lay) == [size(v1lay), size(coef, 1), size(t, 1), &
         size(ven, 1), size(paprs, 1), size(pplay, 1), size(delp, 1), &
         size(d_ven, 1), size(flux_v)], "clvent knon")

    zx_coef(:, 1) = cdrag * (1. + SQRT(u1lay**2 + v1lay**2)) * pplay(:, 1) &
         / (RD * t(:, 1)) * dtphys * RG

    forall (k = 2:klev) zx_coef(:, k) = coef(:, k) * RG / (pplay(:, k - 1) &
         - pplay(:, k)) * (paprs(:, k) * 2 / (t(:, k) + t(:, k - 1)) / RD)**2 &
         * dtphys * RG

    zx_buf = delp(:, 1) + zx_coef(:, 1) + zx_coef(:, 2)
    zx_cv(:, 2) = ven(:, 1) * delp(:, 1) / zx_buf
    zx_dv(:, 2) = zx_coef(:, 2) / zx_buf

    DO k = 3, klev
       zx_buf = delp(:, k - 1) + zx_coef(:, k) &
            + zx_coef(:, k - 1) * (1. - zx_dv(:, k - 1))
       zx_cv(:, k) = (ven(:, k - 1) * delp(:, k - 1) &
            + zx_coef(:, k - 1) * zx_cv(:, k - 1)) / zx_buf
       zx_dv(:, k) = zx_coef(:, k) / zx_buf
    ENDDO

    local_ven(:, klev) = (ven(:, klev) * delp(:, klev) &
         + zx_coef(:, klev) * zx_cv(:, klev)) &
         / (delp(:, klev) + zx_coef(:, klev) &
         - zx_coef(:, klev) * zx_dv(:, klev))

    DO k = klev - 1, 1, - 1
       local_ven(:, k) = zx_cv(:, k + 1) + zx_dv(:, k + 1) * local_ven(:, k + 1)
    ENDDO

    flux_v = zx_coef(:, 1) / (RG * dtphys) * local_ven(:, 1)
    d_ven = local_ven - ven

  END SUBROUTINE clvent

end module clvent_m
