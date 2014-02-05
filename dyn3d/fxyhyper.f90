module fxyhyper_m

  IMPLICIT NONE

contains

  SUBROUTINE fxyhyper(yzoom, grossy, dzoomy, tauy, xzoom, grossx, dzoomx, &
       taux, rlatu, yprimu, rlatv, yprimv, rlatu1, yprimu1, rlatu2, yprimu2, &
       rlonu, xprimu, rlonv, xprimv, rlonm025, xprimm025, rlonp025, xprimp025)

    ! From dyn3d/fxyhyper.F, version 1.1.1.1 2004/05/19 12:53:06

    USE dimens_m, ONLY: jjm
    use fxhyp_m, only: fxhyp
    USE paramet_m, ONLY: iip1, jjp1

    ! Auteur : P. Le Van d'après formulations de R. Sadourny

    ! Cette procédure calcule les latitudes (routine fyhyp) et
    ! longitudes (fxhyp) par des fonctions à tangente hyperbolique.

    ! Il y a trois paramètres, en plus des coordonnées du centre du
    ! zoom (xzoom et yzoom) :

    ! a) le grossissement du zoom : grossy (en y) et grossx (en x)
    ! b) l' extension du zoom : dzoomy (en y) et dzoomx (en x)
    ! c) la raideur de la transition du zoom : taux et tauy 

    ! Nota bene : il vaut mieux avoir : grossx * dzoomx < pi (radians)
    ! et grossy * dzoomy < pi/2 (radians)

    REAL yzoom, grossy, dzoomy, tauy, xzoom, grossx, dzoomx, taux
    REAL rlatu(jjp1), yprimu(jjp1), rlatv(jjm), yprimv(jjm)
    real rlatu1(jjm), yprimu1(jjm), rlatu2(jjm), yprimu2(jjm)
    REAL rlonu(iip1), xprimu(iip1), rlonv(iip1), xprimv(iip1)
    REAL rlonm025(iip1), xprimm025(iip1), rlonp025(iip1), xprimp025(iip1)

    ! Variables locales :

    double precision dxmin, dxmax, dymin, dymax
    INTEGER i, j

    !----------------------------------------------------------

    CALL fyhyp(yzoom, grossy, dzoomy, tauy, rlatu, yprimu, rlatv, yprimv, &
         rlatu2, yprimu2, rlatu1, yprimu1, dymin, dymax)
    CALL fxhyp(xzoom, grossx, dzoomx, taux, rlonm025, xprimm025, rlonv, &
         xprimv, rlonu, xprimu, rlonp025, xprimp025, dxmin, dxmax)

    DO i = 1, iip1
       IF(rlonp025(i).LT.rlonv(i)) THEN
          print *, ' Attention ! rlonp025 < rlonv', i
          STOP 1
       ENDIF

       IF(rlonv(i).LT.rlonm025(i)) THEN 
          print *, ' Attention ! rlonm025 > rlonv', i
          STOP 1
       ENDIF

       IF(rlonp025(i).GT.rlonu(i)) THEN
          print *, ' Attention ! rlonp025 > rlonu', i
          STOP 1
       ENDIF
    ENDDO

    print *, 'Test de coherence ok pour fx'

    DO j = 1, jjm
       IF(rlatu1(j).LE.rlatu2(j)) THEN
          print *, 'Attention ! rlatu1 < rlatu2 ', rlatu1(j), rlatu2(j), j
          STOP 13
       ENDIF

       IF(rlatu2(j).LE.rlatu(j+1)) THEN
          print *, 'Attention ! rlatu2 < rlatup1 ', rlatu2(j), rlatu(j+1), j
          STOP 14
       ENDIF

       IF(rlatu(j).LE.rlatu1(j)) THEN
          print *, ' Attention ! rlatu < rlatu1 ', rlatu(j), rlatu1(j), j
          STOP 15
       ENDIF

       IF(rlatv(j).LE.rlatu2(j)) THEN
          print *, ' Attention ! rlatv < rlatu2 ', rlatv(j), rlatu2(j), j
          STOP 16
       ENDIF

       IF(rlatv(j).ge.rlatu1(j)) THEN
          print *, ' Attention ! rlatv > rlatu1 ', rlatv(j), rlatu1(j), j
          STOP 17
       ENDIF

       IF(rlatv(j).ge.rlatu(j)) THEN
          print *, ' Attention ! rlatv > rlatu ', rlatv(j), rlatu(j), j
          STOP 18
       ENDIF
    ENDDO

    print *, 'Test de coherence ok pour fy'

    print *, 'Latitudes'
    print 3, dymin, dymax
    print *, 'Si cette derniere est trop lache, modifiez les parametres'
    print *, 'grossism, tau, dzoom pour Y et repasser ! '

    print *, ' Longitudes '
    print 3, dxmin, dxmax
    print *, 'Si cette derniere est trop lache, modifiez les parametres'
    print *, 'grossism, tau, dzoom pour Y et repasser ! '

3   Format(1x, ' Au centre du zoom, la longueur de la maille est', &
         ' d environ ', f0.2, ' degres ', /, &
         ' alors que la maille en dehors de la zone du zoom est ', &
         "d'environ", f0.2, ' degres ')

  END SUBROUTINE fxyhyper

end module fxyhyper_m
