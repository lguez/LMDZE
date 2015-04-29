
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/interpost.F,v 1.1.1.1 2004/05/19
! 12:53:06 lmdzadmin Exp $

SUBROUTINE interpost(q, qppm)

  USE dimens_m
  USE paramet_m
  USE comconst
  USE disvert_m
  USE comgeom
  IMPLICIT NONE



  ! Arguments
  REAL q(iip1, jjp1, llm)
  REAL qppm(iim, jjp1, llm)
  ! Local
  INTEGER l, i, j

  ! RE-INVERSION DES NIVEAUX
  ! le programme ppm3d travaille avec une 3ème coordonnée inversée par
  ! rapport
  ! de celle du LMDZ: z=1<=>niveau max, z=llm+1<=>surface
  ! On passe donc des niveaux de Lin à ceux du LMDZ

  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iim
        q(i, j, l) = qppm(i, j, llm-l+1)
      END DO
    END DO
  END DO

  ! BOUCLAGE EN LONGITUDE PAS EFFECTUE DANS PPM3D

  DO l = 1, llm
    DO j = 1, jjp1
      q(iip1, j, l) = q(1, j, l)
    END DO
  END DO


  RETURN

END SUBROUTINE interpost
