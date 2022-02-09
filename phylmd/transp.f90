module transp_m

  IMPLICIT NONE

contains

  SUBROUTINE transp(zmasse, t_seri, q_seri, u_seri, v_seri, zphi)

    ! From LMDZ4/libf/phylmd/transp.F, version 1.1.1.1 2004/05/19 12:53:09

    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: April, 25th 1994

    ! Objet: calculer le transport total de l'\'energie et de la
    ! vapeur d'eau (diagnostique)

    USE dimphy, only: klon, klev
    USE suphec_m, only: rcpd
    USE histwrite_phy_m, ONLY: histwrite_phy

    real, INTENT(IN):: zmasse(:, :) ! (klon, klev)
    ! column-density of mass of air in a cell, in kg m-2

    REAL, INTENT(IN):: t_seri(:, :) ! (klon, klev)
    REAL, INTENT(IN):: q_seri(:, :), u_seri(:, :), v_seri(:, :) ! (klon, klev)

    REAL, INTENT(IN):: zphi(:, :) ! (klon, klev)
    ! geopotential at mid-layer, in m2 s-2

    ! Local:

    real dse(klon, klev)
    ! dry static energy per unit surface in a 3D grid cell, in J m-2

    real sigma_w(klon, klev)
    ! mass of water vapor per unit surface in a 3D grid cell, in kg m-2

    !------------------------------------------------------------------

    dse = (rcpd * t_seri + zphi) * zmasse
    sigma_w = q_seri * zmasse

    ! Intégrales verticales:
    CALL histwrite_phy("ue", sum(u_seri * dse, dim = 2))
    CALL histwrite_phy("ve", sum(v_seri * dse, dim = 2))
    CALL histwrite_phy("uq", sum(u_seri * sigma_w, dim = 2))
    CALL histwrite_phy("vq", sum(v_seri * sigma_w, dim = 2))
    
  END SUBROUTINE transp

end module transp_m
