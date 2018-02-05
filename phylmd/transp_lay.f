module transp_lay_m

  IMPLICIT NONE

contains

  SUBROUTINE transp_lay(paprs, t, q, u, v, geom, vtran_e, vtran_q, utran_e, &
       utran_q)

    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: April, 25th 1994
    ! Objet : calculer le transport de l'\'energie et de la vapeur d'eau

    USE dimphy, only: klon, klev
    USE suphec_m, only: rcpd, rg, rlvtt

    REAL, INTENT(IN):: paprs(klon, klev+1)
    REAL, INTENT(IN):: t(klon, klev)
    REAL, INTENT(IN):: q(klon, klev), u(klon, klev), v(klon, klev)
    REAL, INTENT(IN):: geom(klon, klev)
    REAL, INTENT(out):: vtran_e(klon, klev), vtran_q(klon, klev)
    REAL, INTENT(out):: utran_e(klon, klev), utran_q(klon, klev)

    ! Local:
    INTEGER i, l
    real esh
    
    !------------------------------------------------------------------
    
    DO l = 1, klev
       DO i = 1, klon
          utran_e(i, l) = 0.
          utran_q(i, l) = 0.
          vtran_e(i, l) = 0.
          vtran_q(i, l) = 0.
       END DO
    END DO

    DO l = 1, klev
       DO i = 1, klon
          esh = rcpd * t(i, l) + rlvtt * q(i, l) + geom(i, l)
          utran_e(i, l) = utran_e(i, l) + u(i, l) * esh &
               * (paprs(i, l) - paprs(i, l+1)) / rg
          utran_q(i, l) = utran_q(i, l) + u(i, l) * q(i, l) &
               * (paprs(i, l) - paprs(i, l+1)) / rg
          vtran_e(i, l) = vtran_e(i, l) + v(i, l) * esh &
               * (paprs(i, l) - paprs(i, l+1)) / rg
          vtran_q(i, l) = vtran_q(i, l) + v(i, l) * q(i, l) &
               * (paprs(i, l) - paprs(i, l+1)) / rg
       END DO
    END DO

  END SUBROUTINE transp_lay

end module transp_lay_m
