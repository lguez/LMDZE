module climb_hq_up_m

  IMPLICIT none

contains

  subroutine climb_hq_up(d_t, d_q, cq, dq, ch, dh, ah, aq, bh, bq, flux_t, &
       flux_q, pkf, t, q)

    use conf_gcm_m, only: dtphys
    USE dimphy, ONLY: klev
    USE suphec_m, ONLY: rcpd

    REAL, intent(out):: d_t(:, :) ! (knon, klev) variation of air temperature t
    REAL, intent(out):: d_q(:, :) ! (knon, klev) incrementation de "q"
    REAL, intent(in), dimension(:, 2:):: cq, dq, ch, dh ! (knon, 2:klev)
    REAL, intent(in):: ah(:), aq(:), bh(:), bq(:) ! (knon)

    REAL, intent(in):: flux_t(:) ! (knon)
    ! (diagnostic) flux de chaleur sensible (Cp T) à la surface,
    ! positif vers le bas, W / m2

    REAL, intent(in):: flux_q(:) ! (knon)
    ! flux de la vapeur d'eau à la surface, en kg / (m**2 s)

    REAL, intent(in):: pkf(:, :) ! (knon, klev)
    REAL, intent(in):: t(:, :) ! (knon, klev) temperature (K)
    REAL, intent(in):: q(:, :) ! (knon, klev) humidite specifique (kg / kg)

    ! Local:
    REAL h(size(flux_t), klev) ! (knon, klev) enthalpie potentielle
    INTEGER k
    REAL local_q(size(flux_t), klev) ! (knon, klev)

    !----------------------------------------------------------------------

    h(:, 1) = ah + bh * flux_t * dtphys
    local_q(:, 1) = aq + bq * flux_q * dtphys

    DO k = 2, klev
       h(:, k) = ch(:, k) + dh(:, k) * h(:, k - 1)
       local_q(:, k) = cq(:, k) + dq(:, k) * local_q(:, k - 1)
    ENDDO

    d_t = h / pkf / RCPD - t
    d_q = local_q - q

  end subroutine climb_hq_up
  
end module climb_hq_up_m
