module transp_lay_m

  IMPLICIT NONE

contains

  SUBROUTINE transp_lay(paprs, t_seri, q_seri, u_seri, v_seri, zphi, ve_lay, &
       vq_lay, ue_lay, uq_lay)

    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: April, 25th 1994
    ! Objet : calculer le transport de l'\'energie et de la vapeur d'eau

    USE dimphy, only: klon, klev
    USE suphec_m, only: rcpd, rg, rlvtt
    USE histwrite_phy_m, ONLY: histwrite_phy

    REAL, INTENT(IN):: paprs(:, :) ! (klon, klev + 1)
    ! pression pour chaque inter-couche, en Pa

    REAL, INTENT(IN):: t_seri(:, :) ! (klon, klev)
    REAL, INTENT(IN):: q_seri(:, :), u_seri(:, :), v_seri(:, :) ! (klon, klev)
    REAL, INTENT(IN):: zphi(:, :)
    REAL, INTENT(out):: ve_lay(:, :), vq_lay(:, :) ! (klon, klev)
    REAL, INTENT(out):: ue_lay(:, :), uq_lay(:, :) ! (klon, klev)

    ! Local:
    INTEGER l
    real esh(klon, klev) ! moist static energy, in J kg-1
    real sigma(klon, klev) ! mass per unit surface in a 3D grid cell, in kg m-2

    !------------------------------------------------------------------

    forall (l = 1:klev) sigma(:, l) = (paprs(:, l) - paprs(:, l + 1)) / rg
    esh = rcpd * t_seri + rlvtt * q_seri + zphi
    ue_lay = u_seri * esh * sigma
    uq_lay = u_seri * q_seri * sigma
    ve_lay = v_seri * esh * sigma
    vq_lay = v_seri * q_seri * sigma
    CALL histwrite_phy("ue_lay", ue_lay)
    CALL histwrite_phy("ve_lay", ve_lay)
    CALL histwrite_phy("uq_lay", uq_lay)
    CALL histwrite_phy("vq_lay", vq_lay)

  END SUBROUTINE transp_lay

end module transp_lay_m
