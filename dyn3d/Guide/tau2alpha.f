module tau2alpha_m

  IMPLICIT NONE

contains

  SUBROUTINE tau2alpha(type, factt, taumin, taumax, alpha)

    USE comgeom, ONLY: cu_2d, cv_2d, rlatu, rlatv
    use conf_guide_m, only: lat_min_guide, lat_max_guide
    USE dimens_m, ONLY: iim, jjm
    USE nr_util, ONLY: pi
    USE paramet_m, ONLY: iip1, jjp1
    USE serre, ONLY: clat, clon, grossismx, grossismy
    use writefield_m, only: writefield

    INTEGER, intent(in):: type
    REAL, intent(in):: factt, taumin, taumax
    real, intent(out):: alpha(:, :)

    ! Local:
    REAL dxdy
    REAL, save:: dxdy_min, dxdy_max
    REAL alphamin, alphamax, xi
    REAL, SAVE:: gamma
    INTEGER i, j, ilon, ilat
    LOGICAL:: first = .TRUE.
    REAL dx(iip1, jjp1), dy(iip1, jjp1)
    REAL zlat
    REAL, save:: dxdys(iip1, jjp1), dxdyu(iip1, jjp1), dxdyv(iip1, jjm)

    !------------------------------------------------------------

    IF (first) THEN
       DO j = 2, jjm
          DO i = 2, iip1
             dx(i, j) = 0.5 * (cu_2d(i - 1, j) + cu_2d(i, j)) / cos(rlatu(j))
          END DO
          dx(1, j) = dx(iip1, j)
       END DO
       DO j = 2, jjm
          DO i = 1, iip1
             dy(i, j) = 0.5 * (cv_2d(i, j - 1) + cv_2d(i, j))
          END DO
       END DO
       DO i = 1, iip1
          dx(i, 1) = dx(i, 2)
          dx(i, jjp1) = dx(i, jjm)
          dy(i, 1) = dy(i, 2)
          dy(i, jjp1) = dy(i, jjm)
       END DO

       DO j = 1, jjp1
          DO i = 1, iip1
             dxdys(i, j) = sqrt(dx(i, j)**2 + dy(i, j)**2)
          END DO
       END DO
       CALL writefield("dxdys", dxdys)

       DO j = 1, jjp1
          DO i = 1, iim
             dxdyu(i, j) = 0.5 * (dxdys(i, j) + dxdys(i + 1, j))
          END DO
          dxdyu(iip1, j) = dxdyu(1, j)
       END DO

       DO j = 1, jjm
          DO i = 1, iip1
             dxdyv(i, j) = 0.5 * (dxdys(i, j) + dxdys(i, j + 1))
          END DO
       END DO

       ! coordonnees du centre du zoom
       CALL coordij(clon, clat, ilon, ilat)
       ! aire de la maille au centre du zoom
       dxdy_min = dxdys(ilon, ilat)

       ! dxdy maximal de la maille :
       dxdy_max = 0.
       DO j = 1, jjp1
          DO i = 1, iip1
             dxdy_max = max(dxdy_max, dxdys(i, j))
          END DO
       END DO

       IF (abs(grossismx - 1.) < 0.1 .OR. abs(grossismy - 1.) < 0.1) THEN
          PRINT *, 'Attention : modèle peu zoomé.'
          PRINT *, 'On prend une constante de guidage constante.'
       ELSE
          gamma = (dxdy_max - 2. * dxdy_min) / (dxdy_max - dxdy_min)
          IF (gamma < 1E-5) THEN
             PRINT *, '(dxdy_max - 2. * dxdy_min) / (dxdy_max - dxdy_min) ' &
                  // '< 1e-5'
             STOP 1
          END IF
          gamma = log(0.5) / log(gamma)
          PRINT *, 'gamma=', gamma
       END IF
       first = .false.
    END IF

    alphamin = factt / taumax
    alphamax = factt / taumin

    DO j = 1, size(alpha, 2)
       DO i = 1, size(alpha, 1)
          IF (type==1) THEN
             dxdy = dxdys(i, j)
             zlat = rlatu(j) * 180. / pi
          ELSE IF (type==2) THEN
             dxdy = dxdyu(i, j)
             zlat = rlatu(j) * 180. / pi
          ELSE IF (type==3) THEN
             dxdy = dxdyv(i, j)
             zlat = rlatv(j) * 180. / pi
          END IF
          IF (abs(grossismx - 1.) < 0.1 .OR. abs(grossismy - 1.) < 0.1) THEN
             ! grille regulière
             alpha(i, j) = alphamin
          ELSE
             xi = ((dxdy_max - dxdy) / (dxdy_max - dxdy_min))**gamma
             xi = min(xi, 1.)
             IF (lat_min_guide <= zlat .AND. zlat <= lat_max_guide) THEN
                alpha(i, j) = xi * alphamin + (1. - xi) * alphamax
             ELSE
                alpha(i, j) = 0.
             END IF
          END IF
       END DO
    END DO

  END SUBROUTINE tau2alpha

end module tau2alpha_m
