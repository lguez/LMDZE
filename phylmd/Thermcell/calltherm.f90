module calltherm_m

  implicit none

contains

  subroutine calltherm(pplay, paprs, pphi, u_seri, v_seri, t_seri, q_seri, &
       fm_therm, entr_therm)

    ! From LMDZ4/libf/phylmd/calltherm.F, version 1.2 2004/12/10 11:27:46

    ! Thermiques.

    use conf_gcm_m, only: dtphys
    USE dimphy, ONLY: klev, klon
    use thermcell_m, only: thermcell

    REAL, intent(in):: pplay(klon, klev)
    REAL, intent(in):: paprs(klon, klev + 1)
    REAL, intent(in):: pphi(klon, klev)
    REAL, intent(inout):: u_seri(klon, klev), v_seri(klon, klev)
    REAL, intent(inout):: t_seri(klon, klev)
    real, intent(inout):: q_seri(klon, klev)
    real, intent(out):: fm_therm(klon, klev + 1), entr_therm(klon, klev)

    ! Local:
    REAL d_t_the(klon, klev), d_q_the(klon, klev)
    REAL d_u_the(klon, klev), d_v_the(klon, klev)
    integer i, k

    !----------------------------------------------------------------

    ! tests sur les valeurs negatives de l'eau
    do k = 1, klev
       do i = 1, klon
          if (q_seri(i, k) < 0.) then
             print *, 'Warning: eau < 0 avant thermcell, i =', i, ', k =', k, &
                  ', q =', q_seri(i, k)
             q_seri(i, k) = 1e-15
          endif
       enddo
    enddo

    CALL thermcell(klev, dtphys, pplay, paprs, pphi, u_seri, v_seri, t_seri, &
         q_seri, d_u_the, d_v_the, d_t_the, d_q_the, fm_therm, entr_therm)
    fm_therm(:, klev + 1) = 0.

    ! incrementation des variables meteo
    t_seri = t_seri + d_t_the * dtphys
    u_seri = u_seri + d_u_the * dtphys
    v_seri = v_seri + d_v_the * dtphys
    q_seri = q_seri + d_q_the * dtphys

    ! tests sur les valeurs negatives de l'eau
    DO k = 1, klev
       DO i = 1, klon
          if (q_seri(i, k) < 0.) then
             print *, 'Warning: eau < 0 apr`es thermcell, i =', i, ' k =', &
                  k, ' dq, q', d_q_the(i, k), q_seri(i, k)
             q_seri(i, k) = 1e-15
          endif
       ENDDO
    ENDDO

  end subroutine calltherm

end module calltherm_m
