module transp_lay_m

  IMPLICIT NONE

contains

  SUBROUTINE transp_lay(paprs, t, q, u, v, geom, vtran_e, vtran_q, &
       utran_e, utran_q)

    USE dimens_m
    USE dimphy
    USE suphec_m
    ! ======================================================================
    ! Auteur(s): Z.X.Li (LMD/CNRS)
    ! Date: le 25 avril 1994
    ! Objet: Calculer le transport de l'energie et de la vapeur d'eau
    ! ======================================================================


    REAL, INTENT (IN):: paprs(klon, klev+1)
    REAL, INTENT (IN):: t(klon, klev)
    REAL, INTENT (IN):: q(klon, klev), u(klon, klev), v(klon, klev)
    REAL utran_e(klon, klev), utran_q(klon, klev)
    REAL vtran_e(klon, klev), vtran_q(klon, klev)

    INTEGER i, l
    ! ------------------------------------------------------------------
    REAL geom(klon, klev), esh
    ! ------------------------------------------------------------------
    DO l = 1, klev
       DO i = 1, klon
          utran_e(i, l) = 0.0
          utran_q(i, l) = 0.0
          vtran_e(i, l) = 0.0
          vtran_q(i, l) = 0.0
       END DO
    END DO

    DO l = 1, klev
       DO i = 1, klon
          esh = rcpd*t(i, l) + rlvtt*q(i, l) + geom(i, l)
          utran_e(i, l) = utran_e(i, l) + u(i, l)*esh*(paprs(i,l)-paprs(i,l+1))/ &
               rg
          utran_q(i, l) = utran_q(i, l) + u(i, l)*q(i, l)*(paprs(i,l)-paprs(i,l+1 &
               ))/rg
          vtran_e(i, l) = vtran_e(i, l) + v(i, l)*esh*(paprs(i,l)-paprs(i,l+1))/ &
               rg
          vtran_q(i, l) = vtran_q(i, l) + v(i, l)*q(i, l)*(paprs(i,l)-paprs(i,l+1 &
               ))/rg
       END DO
    END DO

  END SUBROUTINE transp_lay

end module transp_lay_m
