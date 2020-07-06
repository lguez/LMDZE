module ajsec_m

  IMPLICIT none

contains

  SUBROUTINE ajsec(paprs, pplay, t, q, d_t, d_q)

    ! From LMDZ4/libf/phylmd/ajsec.F, version 1.1.1.1 2004/05/19 12:53:08

    ! Author: Z. X. Li (LMD/CNRS) date: 1993/08/18
    ! Objet : ajustement sec

    USE dimphy, ONLY : klev, klon
    USE suphec_m, ONLY : rcpd, rkappa

    REAL, intent(in):: paprs(:, :) ! (klon, klev+1)
    real, intent(in):: pplay(:, :) ! (klon, klev)
    REAL, intent(in):: t(:, :) ! (klon, klev) temperature
    real, intent(in):: q(:, :) ! (klon, klev)

    REAL, intent(out):: d_t(:, :) ! (klon, klev)
    ! incrémentation de la température

    REAL, intent(out):: d_q(:, :) ! (klon, klev)

    ! Local:
    INTEGER, PARAMETER:: limbas=1 ! les couches à ajuster
    REAL, dimension(klon, limbas: klev):: zh, zq, zpk, zpkdp
    REAL hm, sm, qm
    LOGICAL down
    INTEGER i, k, k1, k2

    !--------------------------------------------------------------------

    zpk = pplay(:, limbas:)**RKAPPA
    zh = RCPD * t(:, limbas:) / zpk
    zq = q(:, limbas:)
    forall (k = limbas:klev) &
         zpkdp(:, k) = zpk(:, k) * (paprs(:, k) - paprs(:, k+1))

    ! Correction des profils instables :
    DO i = 1, klon
       IF (any((/(zh(i, k) < zh(i, k - 1), k = limbas + 1, klev)/))) THEN
          ! Profil instable, à modifier
          k2 = limbas
          do while (k2 <= klev - 1)
             k2 = k2 + 1
             IF (zh(i, k2) < zh(i, k2-1)) THEN
                k1 = k2 - 1
                k = k1
                sm = zpkdp(i, k2)
                hm = zh(i, k2)
                qm = zq(i, k2)
                do
                   sm = sm + zpkdp(i, k)
                   hm = hm + zpkdp(i, k) * (zh(i, k)-hm) / sm
                   qm = qm + zpkdp(i, k) * (zq(i, k)-qm) / sm
                   down = .FALSE.
                   IF (k1 /= limbas) THEN
                      IF (hm < zh(i, k1-1)) down = .TRUE.
                   ENDIF
                   IF (down) THEN
                      k1 = k1 - 1
                      k = k1
                   ELSE
                      IF (k2 == klev) exit
                      IF (zh(i, k2 + 1) >= hm) exit
                      k2 = k2 + 1
                      k = k2
                   ENDIF
                end do

                ! nouveau profil : constant (valeur moyenne)
                DO k = k1, k2
                   zh(i, k) = hm
                   zq(i, k) = qm
                ENDDO
                k2 = k2 + 1
             ENDIF
          end do
       ENDIF
    end DO

    d_t(:, : limbas - 1) = 0.
    d_t(:, limbas: klev) = zh * zpk / RCPD - t(:, limbas: klev)
    d_t(:, klev + 1:) = 0.
    d_q = 0.

  END SUBROUTINE ajsec

end module ajsec_m
