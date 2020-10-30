module climb_hq_down_m

  implicit none

contains

  subroutine climb_hq_down(pkf, cq, dq, ch, dh, aq, bq, ah, bh, paprs, pplay, &
       t, coefh, delp, q)

    use calc_coef_m, only: calc_coef
    use comconst, only: dtphys
    USE conf_phys_m, ONLY: iflag_pbl
    USE dimphy, ONLY: klev
    USE suphec_m, ONLY: rcpd, rd, rg, rkappa

    REAL, intent(in):: pkf(:, :) ! (knon, klev)
    REAL, intent(out), dimension(:, 2:):: cq, dq, ch, dh ! (knon, 2:klev)
    REAL, intent(out):: aq(:), bq(:), ah(:), bh(:) ! (knon)

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
    REAL zx_coef(size(paprs, 1), 2:klev) ! (knon, 2:klev)
    REAL gamma(size(paprs, 1), 2:klev) ! (knon, 2:klev)

    real rho(size(paprs, 1), 2:klev) ! (knon, 2:klev)
    ! mass density of air at layer interface

    !----------------------------------------------------------------

    forall (k = 2:klev) &
         rho(:, k) = paprs(:, k) * 2 / (t(:, k) + t(:, k - 1)) / RD

    ! Convertir les coefficients en variables convenables au calcul :
    forall (k = 2:klev) zx_coef(:, k) = coefh(:, k) &
         / (pplay(:, k - 1) - pplay(:, k)) * rho(:, k)**2 * dtphys * RG**2

    ! Counter-gradient for potential enthalpy:
    if (iflag_pbl == 1) then
       gamma(:, 2) = 2.5e-3
       gamma(:, 3:klev)= 1e-3
       forall (k = 2:klev) gamma(:, k) = gamma(:, k) * RCPD * (paprs(:, 1) &
            / paprs(:, k))**RKAPPA / rho(:, k) / RG * (pplay(:, k - 1) &
            - pplay(:, k))
    else
       gamma = 0.
    endif

    call calc_coef(ch, dh, ah, bh, RCPD * t * pkf, gamma, delp, zx_coef)

    ! Counter-gradient for humidity:
    if (iflag_pbl == 1) gamma = 0.

    call calc_coef(cq, dq, aq, bq, q, gamma, delp, zx_coef)
    
  end subroutine climb_hq_down
  
end module climb_hq_down_m
