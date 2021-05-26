module transp_lay_m

  IMPLICIT NONE

contains

  SUBROUTINE transp_lay(paprs, t_seri, q_seri, u_seri, v_seri, zphi, ve_lay, vq_lay, ue_lay, &
       uq_lay)

    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: April, 25th 1994
    ! Objet : calculer le transport de l'\'energie et de la vapeur d'eau

    USE dimphy, only: klon, klev
    USE suphec_m, only: rcpd, rg, rlvtt
    USE histwrite_phy_m, ONLY: histwrite_phy

    REAL, INTENT(IN):: paprs(klon, klev+1)
    REAL, INTENT(IN):: t_seri(klon, klev)
    REAL, INTENT(IN):: q_seri(klon, klev), u_seri(klon, klev), v_seri(klon, klev)
    REAL, INTENT(IN):: zphi(klon, klev)
    REAL, INTENT(out):: ve_lay(klon, klev), vq_lay(klon, klev)
    REAL, INTENT(out):: ue_lay(klon, klev), uq_lay(klon, klev)

    ! Local:
    INTEGER i, l
    real esh
    
    !------------------------------------------------------------------
    
    DO l = 1, klev
       DO i = 1, klon
          ue_lay(i, l) = 0.
          uq_lay(i, l) = 0.
          ve_lay(i, l) = 0.
          vq_lay(i, l) = 0.
       END DO
    END DO

    DO l = 1, klev
       DO i = 1, klon
          esh = rcpd * t_seri(i, l) + rlvtt * q_seri(i, l) + zphi(i, l)
          ue_lay(i, l) = ue_lay(i, l) + u_seri(i, l) * esh &
               * (paprs(i, l) - paprs(i, l+1)) / rg
          uq_lay(i, l) = uq_lay(i, l) + u_seri(i, l) * q_seri(i, l) &
               * (paprs(i, l) - paprs(i, l+1)) / rg
          ve_lay(i, l) = ve_lay(i, l) + v_seri(i, l) * esh &
               * (paprs(i, l) - paprs(i, l+1)) / rg
          vq_lay(i, l) = vq_lay(i, l) + v_seri(i, l) * q_seri(i, l) &
               * (paprs(i, l) - paprs(i, l+1)) / rg
       END DO
    END DO

    CALL histwrite_phy("ue_lay", ue_lay)
    CALL histwrite_phy("ve_lay", ve_lay)
    CALL histwrite_phy("uq_lay", uq_lay)
    CALL histwrite_phy("vq_lay", vq_lay)
    
  END SUBROUTINE transp_lay

end module transp_lay_m
