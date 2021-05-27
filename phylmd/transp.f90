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

    INTEGER l

    real esh(klon, klev)
    ! moist static energy per unit surface in a 3D grid cell, in J m-2

    real sigma(klon, klev) ! mass per unit surface in a 3D grid cell, in kg m-2

    real sigma_w(klon, klev)
    ! mass of water vapor per unit surface in a 3D grid cell, in kg m-2

    !------------------------------------------------------------------

    forall (l = 1:klev) sigma(:, l) = (paprs(:, l) - paprs(:, l + 1)) / rg
    esh = (rcpd * t_seri + rlvtt * q_seri + zphi) * sigma
    sigma_w = q_seri * sigma
    ue = sum(u_seri * esh, dim = 2)
    uq = sum(u_seri * sigma_w, dim = 2)
    ve = sum(v_seri * esh, dim = 2)
    vq = sum(v_seri * sigma_w, dim = 2)

    CALL histwrite_phy("ue", ue)
    CALL histwrite_phy("ve", ve)
    CALL histwrite_phy("uq", uq)
    CALL histwrite_phy("vq", vq)
    
  END SUBROUTINE transp

end module transp_m
