module climb_hq_down_m

  implicit none

contains

  subroutine climb_hq_down(pkf, cq, dq, ch, dh, paprs, pplay, t, coefh, delp, q)

    use calc_coef_m, only: calc_coef
    use comconst, only: dtphys
    USE conf_phys_m, ONLY: iflag_pbl
    USE dimphy, ONLY: klev
    USE suphec_m, ONLY: rcpd, rd, rg, rkappa

    REAL, intent(in):: pkf(:, :) ! (knon, klev)
    REAL, intent(out), dimension(:, :):: cq, dq, ch, dh ! (knon, klev)

    REAL, intent(in):: paprs(:, :) ! (knon, klev + 1)
    ! pression \`a inter-couche (Pa)

    REAL, intent(in):: pplay(:, :) ! (knon, klev)
    ! pression au milieu de couche (Pa)

    REAL, intent(in):: t(:, :) ! (knon, klev) temperature (K)

    REAL, intent(in):: coefh(:, 2:) ! (knon, 2:klev)
    ! diffusion coefficient at layer interface, for heat and humidity, in m2 s-1

    REAL, intent(in):: delp(:, :) ! (knon, klev)
    ! epaisseur de couche en pression (Pa)

    REAL, intent(in):: q(:, :) ! (knon, klev) humidite specifique (kg / kg)

    ! Local:
    
    INTEGER k
    REAL h(size(paprs, 1), klev) ! (knon, klev) enthalpie potentielle
    REAL zx_coef(size(paprs, 1), 2:klev) ! (knon, 2:klev)

    REAL gamt(size(paprs, 1), 2:klev) ! (knon, 2:klev)
    ! contre-gradient pour la chaleur sensible, en J kg-1 m-1

    REAL gamma(size(paprs, 1), 2:klev) ! (knon, 2:klev)

    real rho(size(paprs, 1), 2:klev) ! (knon, 2:klev)
    ! mass density of air at layer interface

    !----------------------------------------------------------------

    h = RCPD * t * pkf

    forall (k = 2:klev) &
         rho(:, k) = paprs(:, k) * 2 / (t(:, k) + t(:, k - 1)) / RD

    ! Convertir les coefficients en variables convenables au calcul :
    forall (k = 2:klev) zx_coef(:, k) = coefh(:, k) &
         / (pplay(:, k - 1) - pplay(:, k)) * rho(:, k)**2 * dtphys * RG**2

    ! Pr\'eparer les flux li\'es aux contre-gradients :

    if (iflag_pbl == 1) then
       gamt(:, 2) = 2.5e-3 * RCPD * (paprs(:, 1) / paprs(:, 2))**RKAPPA
       forall (k = 3:klev) gamt(:, k)= 1e-3 * RCPD * (paprs(:, 1) &
            / paprs(:, k))**RKAPPA
       forall (k = 2:klev) gamma(:, k) = gamt(:, k) / rho(:, k) / RG &
            * (pplay(:, k - 1) - pplay(:, k))
    else
       gamma = 0.
    endif

    call calc_coef(ch, dh, h, gamma, delp, zx_coef)
    if (iflag_pbl == 1) gamma = 0.
    call calc_coef(cq, dq, q, gamma, delp, zx_coef)
    
  end subroutine climb_hq_down
  
end module climb_hq_down_m
