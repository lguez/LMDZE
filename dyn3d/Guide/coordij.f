
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/coordij.F,v 1.1.1.1 2004/05/19
! 12:53:05 lmdzadmin Exp $

SUBROUTINE coordij(lon, lat, ilon, jlat)

  ! =======================================================================

  ! calcul des coordonnees i et j de la maille scalaire dans
  ! laquelle se trouve le point (lon,lat) en radian

  ! =======================================================================

  USE dimens_m
  USE paramet_m
  USE comconst
  USE comgeom
  USE serre
  USE nr_util, ONLY: pi
  IMPLICIT NONE
  REAL, intent(in):: lon, lat
  INTEGER ilon, jlat
  INTEGER i, j


  REAL zlon, zlat

  zlon = lon*pi/180.
  zlat = lat*pi/180.

  DO i = 1, iim + 1
    IF (rlonu(i)>zlon) THEN
      ilon = i
      GO TO 10
    END IF
  END DO
10 CONTINUE

  j = 0
  DO j = 1, jjm
    IF (rlatv(j)<zlat) THEN
      jlat = j
      GO TO 20
    END IF
  END DO
20 CONTINUE
  IF (j==0) j = jjm + 1

  RETURN
END SUBROUTINE coordij
