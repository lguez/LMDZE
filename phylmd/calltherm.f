module calltherm_m

  implicit none

contains

  subroutine calltherm(pplay, paprs, pphi, u_seri, v_seri, t_seri, q_seri, &
       d_u_ajs, d_v_ajs, d_t_ajs, d_q_ajs, fm_therm, entr_therm)

    ! From LMDZ4/libf/phylmd/calltherm.F, version 1.2 2004/12/10 11:27:46

    ! Thermiques.

    use comconst, only: dtphys
    USE ctherm, ONLY: l_mix_thermals, nsplit_thermals, r_aspect_thermals, &
         tho_thermals, w2di_thermals
    USE dimphy, ONLY: klev, klon
    use thermcell_m, only: thermcell

    REAL, intent(in):: pplay(klon, klev)
    REAL, intent(in):: paprs(klon, klev+1)
    REAL, intent(in):: pphi(klon, klev)
    REAL, intent(inout):: u_seri(klon, klev), v_seri(klon, klev)
    REAL, intent(inout):: t_seri(klon, klev)
    real, intent(inout):: q_seri(klon, klev)

    ! Update thermiques
    REAL d_u_ajs(klon, klev), d_v_ajs(klon, klev)
    REAL d_t_ajs(klon, klev), d_q_ajs(klon, klev)
    real fm_therm(klon, klev+1), entr_therm(klon, klev)

    ! Local:

    REAL d_t_the(klon, klev), d_q_the(klon, klev)
    REAL d_u_the(klon, klev), d_v_the(klon, klev)

    real, save:: zfm_therm(klon, klev+1), zentr_therm(klon, klev)
    real zdt

    integer i, k, isplit

    !----------------------------------------------------------------

    ! Modele du thermique
    print*, 'avant isplit ', nsplit_thermals

    fm_therm=0.
    entr_therm=0.

    ! tests sur les valeurs negatives de l'eau
    do k=1, klev
       do i=1, klon
          if (.not.q_seri(i, k) >= 0.) then
             print*, 'WARN eau<0 avant therm i=', i, ' k=', k, ' dq, q', &
                  d_q_the(i, k), q_seri(i, k)
             q_seri(i, k)=1.e-15
          endif
       enddo
    enddo

    zdt=dtphys/real(nsplit_thermals)
    do isplit = 1, nsplit_thermals
       CALL thermcell(klon, klev, zdt, pplay, paprs, pphi, u_seri, v_seri, &
            t_seri, q_seri, d_u_the, d_v_the, d_t_the, d_q_the, zfm_therm, &
            zentr_therm, r_aspect_thermals, l_mix_thermals, w2di_thermals, &
            tho_thermals)

       ! transformation de la derivee en tendance
       d_t_the=d_t_the*dtphys/real(nsplit_thermals)
       d_u_the=d_u_the*dtphys/real(nsplit_thermals)
       d_v_the=d_v_the*dtphys/real(nsplit_thermals)
       d_q_the=d_q_the*dtphys/real(nsplit_thermals)
       fm_therm=fm_therm +zfm_therm/real(nsplit_thermals)
       entr_therm=entr_therm +zentr_therm/real(nsplit_thermals)
       fm_therm(:, klev+1)=0.

       ! accumulation de la tendance
       d_t_ajs=d_t_ajs+d_t_the
       d_u_ajs=d_u_ajs+d_u_the
       d_v_ajs=d_v_ajs+d_v_the
       d_q_ajs=d_q_ajs+d_q_the

       ! incrementation des variables meteo
       t_seri = t_seri + d_t_the
       u_seri = u_seri + d_u_the
       v_seri = v_seri + d_v_the
       q_seri = q_seri + d_q_the

       ! tests sur les valeurs negatives de l'eau
       DO k = 1, klev
          DO i = 1, klon
             if (.not.q_seri(i, k) >= 0.) then
                print*, 'WARN eau<0 apres therm i=', i, ' k=', k, ' dq, q', &
                     d_q_the(i, k), q_seri(i, k)
                q_seri(i, k)=1.e-15
             endif
          ENDDO
       ENDDO
    enddo

  end subroutine calltherm

end module calltherm_m
