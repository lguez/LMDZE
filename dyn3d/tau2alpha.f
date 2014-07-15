module tau2alpha_m

    USE paramet_m, ONLY : iip1, jjp1
    USE dimens_m, ONLY : jjm

    IMPLICIT NONE

    private iip1, jjp1, jjm

    REAL dxdys(iip1, jjp1), dxdyu(iip1, jjp1), dxdyv(iip1, jjm)

contains

  SUBROUTINE tau2alpha(type, pim, pjm, factt, taumin, taumax, alpha)

    USE comgeom, ONLY : cu_2d, cv_2d, rlatu, rlatv
    use conf_guide_m, only: lat_min_guide, lat_max_guide
    use dump2d_m, only: dump2d
    USE dimens_m, ONLY : iim
    USE nr_util, ONLY : pi
    USE serre, ONLY : clat, clon, grossismx, grossismy

    !   arguments :
    INTEGER type
    INTEGER pim, pjm
    REAL, intent(in):: factt, taumin, taumax
    REAL dxdy_, alpha(pim, pjm)
    REAL dxdy_min, dxdy_max

    !  local :
    REAL alphamin, alphamax, gamma, xi
    SAVE gamma
    INTEGER i, j, ilon, ilat

    LOGICAL first
    SAVE first
    DATA first/ .TRUE./

    REAL zdx(iip1, jjp1), zdy(iip1, jjp1)
    REAL zlat

    !------------------------------------------------------------

    IF (first) THEN
       DO j = 2, jjm
          DO i = 2, iip1
             zdx(i, j) = 0.5*(cu_2d(i-1, j)+cu_2d(i, j))/cos(rlatu(j))
          END DO
          zdx(1, j) = zdx(iip1, j)
       END DO
       DO j = 2, jjm
          DO i = 1, iip1
             zdy(i, j) = 0.5*(cv_2d(i, j-1)+cv_2d(i, j))
          END DO
       END DO
       DO i = 1, iip1
          zdx(i, 1) = zdx(i, 2)
          zdx(i, jjp1) = zdx(i, jjm)
          zdy(i, 1) = zdy(i, 2)
          zdy(i, jjp1) = zdy(i, jjm)
       END DO
       DO j = 1, jjp1
          DO i = 1, iip1
             dxdys(i, j) = sqrt(zdx(i, j)*zdx(i, j)+zdy(i, j)*zdy(i, j))
          END DO
       END DO
       DO j = 1, jjp1
          DO i = 1, iim
             dxdyu(i, j) = 0.5*(dxdys(i, j)+dxdys(i+1, j))
          END DO
          dxdyu(iip1, j) = dxdyu(1, j)
       END DO
       DO j = 1, jjm
          DO i = 1, iip1
             dxdyv(i, j) = 0.5*(dxdys(i, j)+dxdys(i, j + 1))
          END DO
       END DO

       CALL dump2d(iip1, jjp1, dxdys, 'DX2DY2 SCAL  ')
       CALL dump2d(iip1, jjp1, dxdyu, 'DX2DY2 U     ')
       CALL dump2d(iip1, jjp1, dxdyv, 'DX2DY2 v     ')

       !   coordonnees du centre du zoom
       CALL coordij(clon, clat, ilon, ilat)
       !   aire de la maille au centre du zoom
       dxdy_min = dxdys(ilon, ilat)
       !   dxdy maximale de la maille
       dxdy_max = 0.
       DO j = 1, jjp1
          DO i = 1, iip1
             dxdy_max = max(dxdy_max, dxdys(i, j))
          END DO
       END DO

       IF (abs(grossismx-1.)<0.1 .OR. abs(grossismy-1.)<0.1) THEN
          PRINT *, 'ATTENTION modele peu zoome'
          PRINT *, 'ATTENTION on prend une constante de guidage cste'
          gamma = 0.
       ELSE
          gamma = (dxdy_max-2.*dxdy_min)/(dxdy_max-dxdy_min)
          PRINT *, 'gamma=', gamma
          IF (gamma<1.E-5) THEN
             PRINT *, 'gamma =', gamma, '<1e-5'
             STOP
          END IF
          PRINT *, 'gamma=', gamma
          gamma = log(0.5)/log(gamma)
       END IF
    END IF

    alphamin = factt/taumax
    alphamax = factt/taumin

    DO j = 1, pjm
       DO i = 1, pim
          IF (type==1) THEN
             dxdy_ = dxdys(i, j)
             zlat = rlatu(j)*180./pi
          ELSE IF (type==2) THEN
             dxdy_ = dxdyu(i, j)
             zlat = rlatu(j)*180./pi
          ELSE IF (type==3) THEN
             dxdy_ = dxdyv(i, j)
             zlat = rlatv(j)*180./pi
          END IF
          IF (abs(grossismx-1.)<0.1 .OR. abs(grossismy-1.)<0.1) THEN
             !  pour une grille reguliere, xi=xxx**0=1 -> alpha=alphamin
             alpha(i, j) = alphamin
          ELSE
             xi = ((dxdy_max-dxdy_)/(dxdy_max-dxdy_min))**gamma
             xi = min(xi, 1.)
             IF (lat_min_guide<=zlat .AND. zlat<=lat_max_guide) THEN
                alpha(i, j) = xi*alphamin + (1.-xi)*alphamax
             ELSE
                alpha(i, j) = 0.
             END IF
          END IF
       END DO
    END DO

  END SUBROUTINE tau2alpha

end module tau2alpha_m
