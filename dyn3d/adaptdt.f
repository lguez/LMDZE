
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/adaptdt.F,v 1.1.1.1 2004/05/19
! 12:53:05 lmdzadmin Exp $

SUBROUTINE adaptdt(dtbon, n, pbaru, masse)

  USE dimens_m
  USE paramet_m
  USE comconst
  USE disvert_m
  USE conf_gcm_m
  USE conf_gcm_m
  USE comgeom
  USE temps
  IMPLICIT NONE


  ! ----------------------------------------------------------
  ! Arguments
  ! ----------------------------------------------------------
  INTEGER n
  REAL dtbon
  REAL, INTENT (IN) :: pbaru(iip1, jjp1, llm)
  REAL masse(iip1, jjp1, llm)
  ! ----------------------------------------------------------
  ! Local
  ! ----------------------------------------------------------
  INTEGER i, j, l
  REAL cflmax, aaa, bbb

  cflmax = 0.
  DO l = 1, llm
    DO j = 2, jjm
      DO i = 1, iim
        aaa = pbaru(i, j, l)*dtvr/masse(i, j, l)
        cflmax = max(cflmax, aaa)
        bbb = -pbaru(i, j, l)*dtvr/masse(i+1, j, l)
        cflmax = max(cflmax, bbb)
      END DO
    END DO
  END DO
  n = int(cflmax) + 1
  dtbon = dtvr/n

  RETURN
END SUBROUTINE adaptdt







