module init_tau2alpha_m

  IMPLICIT NONE

  REAL dxdy_min, dxdy_max, gamma

contains

  SUBROUTINE init_tau2alpha(dxdys, dxdyu, dxdyv)

    USE comgeom, ONLY: cu_2d, cv_2d
    use conf_guide_m, only: guide_u, guide_v
    use coordij_m, only: coordij
    USE dimens_m, ONLY: iim, jjm
    USE dynetat0_m, ONLY: clat, clon, grossismx, grossismy, rlatu
    USE paramet_m, ONLY: iip1, jjp1
    use writefield_m, only: writefield

    REAL, intent(out):: dxdys(iip1, jjp1), dxdyu(iip1, jjp1), dxdyv(iip1, jjm)

    ! Local:
    INTEGER i, j, ilon, ilat
    REAL dx(iip1, jjp1), dy(iip1, jjp1)

    !------------------------------------------------------------

    PRINT *, 'Call sequence information: init_tau2alpha'

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

    if (guide_u) then
       DO j = 1, jjp1
          DO i = 1, iim
             dxdyu(i, j) = 0.5 * (dxdys(i, j) + dxdys(i + 1, j))
          END DO
          dxdyu(iip1, j) = dxdyu(1, j)
       END DO
    end if

    if (guide_v) then
       DO j = 1, jjm
          DO i = 1, iip1
             dxdyv(i, j) = 0.5 * (dxdys(i, j) + dxdys(i, j + 1))
          END DO
       END DO
    end if

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

  END SUBROUTINE init_tau2alpha

end module init_tau2alpha_m
