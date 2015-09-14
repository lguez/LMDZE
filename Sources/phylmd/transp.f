module transp_m

  IMPLICIT NONE

contains

  SUBROUTINE transp(paprs, tsol, t, q, u, v, geom, vtran_e, vtran_q, utran_e, &
       utran_q)

    ! From LMDZ4/libf/phylmd/transp.F,v 1.1.1.1 2004/05/19 12:53:09

    USE dimens_m
    USE dimphy
    USE suphec_m
    ! ======================================================================
    ! Auteur(s): Z.X.Li (LMD/CNRS)
    ! Date: le 25 avril 1994
    ! Objet: Calculer le transport total de l'energie et de la vapeur d'eau
    ! ======================================================================


    REAL, INTENT (IN) :: paprs(klon, klev+1)
    REAL tsol(klon)
    REAL, INTENT (IN) :: t(klon, klev)
    REAL, INTENT (IN) :: q(klon, klev), u(klon, klev), v(klon, klev)
    REAL utran_e(klon), utran_q(klon), vtran_e(klon), vtran_q(klon)

    INTEGER i, l
    ! ------------------------------------------------------------------
    REAL geom(klon, klev), e
    ! ------------------------------------------------------------------
    DO i = 1, klon
       utran_e(i) = 0.0
       utran_q(i) = 0.0
       vtran_e(i) = 0.0
       vtran_q(i) = 0.0
    END DO

    DO l = 1, klev
       DO i = 1, klon
          e = rcpd*t(i, l) + rlvtt*q(i, l) + geom(i, l)
          utran_e(i) = utran_e(i) + u(i, l)*e*(paprs(i,l)-paprs(i,l+1))/rg
          utran_q(i) = utran_q(i) + u(i, l)*q(i, l)*(paprs(i,l)-paprs(i,l+1))/rg
          vtran_e(i) = vtran_e(i) + v(i, l)*e*(paprs(i,l)-paprs(i,l+1))/rg
          vtran_q(i) = vtran_q(i) + v(i, l)*q(i, l)*(paprs(i,l)-paprs(i,l+1))/rg
       END DO
    END DO

  END SUBROUTINE transp

end module transp_m
