module coordij_m

  IMPLICIT NONE

contains

  SUBROUTINE coordij(lon, lat, ilon, jlat)

    ! From LMDZ4/libf/dyn3d/coordij.F, version 1.1.1.1 2004/05/19 12:53:05

    ! calcul des coordonnees i et j de la maille scalaire dans
    ! laquelle se trouve le point (lon, lat) en radian

    USE dimens_m, only: iim, jjm
    USE comgeom, only: rlonu, rlatv
    USE nr_util, ONLY: pi

    REAL, intent(in):: lon, lat
    INTEGER ilon, jlat
    INTEGER i, j

    DO i = 1, iim + 1
       IF (rlonu(i)>lon) THEN
          ilon = i
          exit
       END IF
    END DO

    j = 0
    DO j = 1, jjm
       IF (rlatv(j)<lat) THEN
          jlat = j
          exit
       END IF
    END DO

    IF (j==0) j = jjm + 1

  END SUBROUTINE coordij

end module coordij_m
