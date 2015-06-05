module tau2alpha_m

  IMPLICIT NONE

contains

  SUBROUTINE tau2alpha(dxdy, rlat, taumin, taumax, alpha)

    use conf_guide_m, only: lat_min_guide, lat_max_guide, factt
    use init_tau2alpha_m, only: dxdy_min, dxdy_max, gamma
    USE nr_util, ONLY: pi, assert_eq

    REAL, intent(in):: dxdy(:, :) ! (n_lon, n_lat)
    REAL, intent(in):: rlat(:) ! (n_lat)
    REAL, intent(in):: taumin, taumax
    real, intent(out):: alpha(:, :) ! (n_lon, n_lat)

    ! Local:
    REAL alphamin, alphamax, xi
    INTEGER i, j, n_lon, n_lat

    !------------------------------------------------------------

    PRINT *, 'Call sequence information: tau2alpha'

    n_lon = assert_eq(size(alpha, 1), size(dxdy, 1), "tau2alpha n_lon")
    n_lat = assert_eq(size(alpha, 2), size(dxdy, 2), size(rlat), &
         "tau2alpha n_lat")

    alphamin = factt / taumax
    alphamax = factt / taumin

    DO j = 1, n_lat
       IF (lat_min_guide <= rlat(j) .AND. rlat(j) <= lat_max_guide) THEN
          DO i = 1, n_lon
             xi = min(((dxdy_max - dxdy(i, j)) &
                  / (dxdy_max - dxdy_min))**gamma, 1.)
             alpha(i, j) = xi * alphamin + (1. - xi) * alphamax
          END DO
       ELSE
          alpha(:, j) = 0.
       END IF
    END DO

  END SUBROUTINE tau2alpha

end module tau2alpha_m
