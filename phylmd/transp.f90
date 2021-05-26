module transp_m

  IMPLICIT NONE

contains

  SUBROUTINE transp(paprs, t_seri, q_seri, u_seri, v_seri, zphi, ve, vq, ue, uq)

    ! From LMDZ4/libf/phylmd/transp.F,v 1.1.1.1 2004/05/19 12:53:09

    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 25 avril 1994
    ! Objet: calculer le transport total de l'\'energie et de la vapeur d'eau

    USE dimensions
    USE dimphy
    USE suphec_m
    USE histwrite_phy_m, ONLY: histwrite_phy

    REAL, INTENT (IN) :: paprs(klon, klev+1)
    REAL, INTENT (IN) :: t_seri(klon, klev)
    REAL, INTENT (IN) :: q_seri(klon, klev), u_seri(klon, klev), v_seri(klon, klev)
    REAL ue(klon), uq(klon), ve(klon), vq(klon)

    INTEGER i, l
    ! ------------------------------------------------------------------
    REAL zphi(klon, klev), e
    ! ------------------------------------------------------------------
    DO i = 1, klon
       ue(i) = 0.0
       uq(i) = 0.0
       ve(i) = 0.0
       vq(i) = 0.0
    END DO

    DO l = 1, klev
       DO i = 1, klon
          e = rcpd*t_seri(i, l) + rlvtt*q_seri(i, l) + zphi(i, l)
          ue(i) = ue(i) + u_seri(i, l)*e*(paprs(i,l)-paprs(i,l+1))/rg
          uq(i) = uq(i) + u_seri(i, l)*q_seri(i, l)*(paprs(i,l)-paprs(i,l+1))/rg
          ve(i) = ve(i) + v_seri(i, l)*e*(paprs(i,l)-paprs(i,l+1))/rg
          vq(i) = vq(i) + v_seri(i, l)*q_seri(i, l)*(paprs(i,l)-paprs(i,l+1))/rg
       END DO
    END DO

    CALL histwrite_phy("ue", ue)
    CALL histwrite_phy("ve", ve)
    CALL histwrite_phy("uq", uq)
    CALL histwrite_phy("vq", vq)
    
  END SUBROUTINE transp

end module transp_m
