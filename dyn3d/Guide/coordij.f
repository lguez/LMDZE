module coordij_m

  IMPLICIT NONE

contains

  SUBROUTINE coordij(lon, lat, ilon, jlat)

    ! From LMDZ4/libf/dyn3d/coordij.F, version 1.1.1.1 2004/05/19 12:53:05

    ! Calcul des coordonnées ilon et jlat de la maille scalaire dans
    ! laquelle se trouve le point (lon, lat).

    USE dimens_m, only: iim, jjm
    USE dynetat0_m, only: rlonu, rlatv

    REAL, intent(in):: lon, lat ! in rad
    INTEGER, intent(out):: ilon, jlat

    !----------------------------------------------------------

    ilon = 1
    do while (ilon <= iim .and. rlonu(ilon) <= lon)
       ilon = ilon + 1
    end do

    jlat = 1
    do while (jlat <= jjm - 1 .and. rlatv(jlat) >= lat)
       jlat = jlat + 1
    end do

  END SUBROUTINE coordij

end module coordij_m
