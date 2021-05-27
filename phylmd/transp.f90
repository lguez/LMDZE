module transp_m

  IMPLICIT NONE

contains

  SUBROUTINE transp(paprs, t_seri, q_seri, u_seri, v_seri, zphi, ve, vq, ue, uq)

    ! From LMDZ4/libf/phylmd/transp.F, version 1.1.1.1 2004/05/19 12:53:09

    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: April, 25th 1994

    ! Objet: calculer le transport total de l'\'energie et de la
    ! vapeur d'eau (diagnostique)

    USE dimphy, only: klon, klev
    USE suphec_m, only: rcpd, rg, rlvtt
    USE histwrite_phy_m, ONLY: histwrite_phy

    REAL, INTENT(IN):: paprs(:, :) ! (klon, klev + 1)
    ! pression pour chaque inter-couche, en Pa

    REAL, INTENT(IN):: t_seri(:, :) ! (klon, klev)
    REAL, INTENT(IN):: q_seri(:, :), u_seri(:, :), v_seri(:, :) ! (klon, klev)
    REAL, INTENT(IN):: zphi(:, :) ! (klon, klev)
    REAL, INTENT(out):: ve(:), vq(:), ue(:), uq(:) ! (klon)

    ! Local:
    INTEGER i, l
    real e

    !------------------------------------------------------------------

    ue = 0.
    uq = 0.
    ve = 0.
    vq = 0.

    DO l = 1, klev
       DO i = 1, klon
          e = rcpd * t_seri(i, l) + rlvtt * q_seri(i, l) + zphi(i, l)
          ue(i) = ue(i) + u_seri(i, l) * e * (paprs(i, l) - paprs(i, l + 1)) &
               / rg
          uq(i) = uq(i) + u_seri(i, l) * q_seri(i, l) * (paprs(i, l) &
               - paprs(i, l + 1)) / rg
          ve(i) = ve(i) + v_seri(i, l) * e * (paprs(i, l) - paprs(i, l + 1)) &
               / rg
          vq(i) = vq(i) + v_seri(i, l) * q_seri(i, l) * (paprs(i, l) &
               - paprs(i, l + 1)) / rg
       END DO
    END DO

    CALL histwrite_phy("ue", ue)
    CALL histwrite_phy("ve", ve)
    CALL histwrite_phy("uq", uq)
    CALL histwrite_phy("vq", vq)
    
  END SUBROUTINE transp

end module transp_m
