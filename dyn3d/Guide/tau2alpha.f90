module tau2alpha_m

  IMPLICIT NONE

contains

  SUBROUTINE tau2alpha(lat_min_guide, lat_max_guide, factt, dxdy, rlat, &
       taumin, taumax, alpha)

    USE jumble, ONLY: assert_eq

    use init_tau2alpha_m, only: dxdy_min, dxdy_max, gamma

    ! Dans le cas où on n'a les analyses que sur une bande de latitudes :
    REAL, intent(in):: lat_min_guide ! minimum latitude for nudging, in rad
    real, intent(in):: lat_max_guide ! maximum latitude for nudging, in rad

    REAL, intent(in):: factt
    ! pas de temps entre deux appels au guidage, en jours
    
    REAL, intent(in):: dxdy(:, :) ! (n_lon, n_lat)
    REAL, intent(in):: rlat(:) ! (n_lat)
    REAL, intent(in):: taumin, taumax
    real, intent(out):: alpha(:, :) ! (n_lon, n_lat)

    ! Local:
    REAL a_min, a_max, xi
    INTEGER i, j, n_lon, n_lat

    !------------------------------------------------------------

    PRINT *, 'Call sequence information: tau2alpha'

    n_lon = assert_eq(size(alpha, 1), size(dxdy, 1), "tau2alpha n_lon")
    n_lat = assert_eq(size(alpha, 2), size(dxdy, 2), size(rlat), &
         "tau2alpha n_lat")

    a_min = factt / taumax
    a_max = factt / taumin

    DO j = 1, n_lat
       IF (lat_min_guide <= rlat(j) .AND. rlat(j) <= lat_max_guide) THEN
          DO i = 1, n_lon
             xi = min(((dxdy_max - dxdy(i, j)) &
                  / (dxdy_max - dxdy_min))**gamma, 1.)
             alpha(i, j) = 1. - exp(- xi * a_min - (1. - xi) * a_max)
          END DO
       ELSE
          alpha(:, j) = 0.
       END IF
    END DO

  END SUBROUTINE tau2alpha

end module tau2alpha_m
