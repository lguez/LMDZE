module dist_sphe_m

  ! From grid_atob.F, v 1.1.1.1 2004/05/19 12:53:05

  IMPLICIT none

contains

  SUBROUTINE dist_sphe(rf_lon, rf_lat, rlon, rlat, im, jm, distance)

    ! Auteur: Laurent Li (le 30 decembre 1996)

    ! Ce programme calcule la distance minimale (selon le grand cercle)
    ! entre deux points sur la terre

    use nr_util, only: pi

    INTEGER, intent(in):: im, jm ! dimensions
    REAL, intent(in):: rf_lon ! longitude du point de reference (degres)
    REAL, intent(in):: rf_lat ! latitude du point de reference (degres)
    REAL, intent(in):: rlon(im), rlat(jm) ! longitude et latitude des points

    REAL, intent(out):: distance(im, jm) ! distances en metre

    REAL rlon1, rlat1
    REAL rlon2, rlat2
    REAL dist
    REAL pa, pb, p

    REAL radius
    PARAMETER (radius=6371229.)
    integer i, j

    !---------------------------------------------------------------------

    DO j = 1, jm
       DO i = 1, im
          rlon1=rf_lon
          rlat1=rf_lat
          rlon2=rlon(i)
          rlat2=rlat(j)
          pa = pi/2.0 - rlat1*pi/180.0 ! dist. entre pole n et point a
          pb = pi/2.0 - rlat2*pi/180.0 ! dist. entre pole n et point b
          p = (rlon1-rlon2)*pi/180.0 ! angle entre a et b (leurs meridiens)

          dist = ACOS(COS(pa)*COS(pb) + SIN(pa)*SIN(pb)*COS(p))
          dist = radius * dist
          distance(i, j) = dist
       end DO
    end DO

  END SUBROUTINE dist_sphe

end module dist_sphe_m
