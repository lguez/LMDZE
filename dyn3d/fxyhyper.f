module fxyhyper_m

  IMPLICIT NONE

contains

  SUBROUTINE fxyhyper(rlatu, yprimu, rlatv, yprimv, rlatu1, yprimu1, rlatu2, &
       yprimu2, rlonu, xprimu, rlonv, xprimv, rlonm025, xprimm025, rlonp025, &
       xprimp025)

    ! From dyn3d/fxyhyper.F, version 1.1.1.1, 2004/05/19 12:53:06

    USE dimens_m, ONLY: jjm
    use fxhyp_m, only: fxhyp
    use fyhyp_m, only: fyhyp
    USE paramet_m, ONLY: iip1
    use serre, only: clat, grossismy, dzoomy, tauy, clon, grossismx, dzoomx, &
         taux

    ! Auteur : P. Le Van d'après les formulations de R. Sadourny

    ! f(x, y) à dérivée tangente hyperbolique

    ! Cette procédure calcule les latitudes (routine fyhyp) et
    ! longitudes (fxhyp) par des fonctions tangente hyperbolique.

    ! Il y a trois paramètres, en plus des coordonnées du centre du
    ! zoom (clon et clat) :

    ! a) le grossissement du zoom : grossismy (en y) et grossismx (en x)
    ! b) l'extension du zoom : dzoomy (en y) et dzoomx (en x)
    ! c) la raideur de la transition du zoom : taux et tauy 

    ! Nota bene : il vaut mieux avoir : grossismx * dzoomx < pi (radians)
    ! et grossismy * dzoomy < pi/2 (radians)

    REAL rlatu(:), yprimu(:) ! (jjm + 1)
    real rlatv(:), yprimv(:) ! (jjm)
    real rlatu1(:), yprimu1(:), rlatu2(:), yprimu2(:) ! (jjm)
    REAL rlonu(:), xprimu(:), rlonv(:), xprimv(:) ! (iim + 1)
    REAL rlonm025(:), xprimm025(:), rlonp025(:), xprimp025(:) ! (iim + 1)

    ! Local:

    double precision dxmin, dxmax, dymin, dymax
    INTEGER i, j

    !----------------------------------------------------------

    CALL fyhyp(clat, grossismy, dzoomy, tauy, rlatu, yprimu, rlatv, yprimv, &
         rlatu2, yprimu2, rlatu1, yprimu1, dymin, dymax)
    CALL fxhyp(clon, grossismx, dzoomx, taux, rlonm025, xprimm025, rlonv, &
         xprimv, rlonu, xprimu, rlonp025, xprimp025, dxmin, dxmax)

    DO i = 1, iip1
       IF (rlonp025(i).LT.rlonv(i)) THEN
          print *, ' Attention ! rlonp025 < rlonv', i
          STOP 1
       ENDIF

       IF (rlonv(i).LT.rlonm025(i)) THEN 
          print *, ' Attention ! rlonm025 > rlonv', i
          STOP 1
       ENDIF

       IF (rlonp025(i).GT.rlonu(i)) THEN
          print *, ' Attention ! rlonp025 > rlonu', i
          STOP 1
       ENDIF
    ENDDO

    print *, 'Test de coherence ok pour fx'

    DO j = 1, jjm
       IF (rlatu1(j).LE.rlatu2(j)) THEN
          print *, 'Attention ! rlatu1 < rlatu2 ', rlatu1(j), rlatu2(j), j
          STOP 13
       ENDIF

       IF (rlatu2(j).LE.rlatu(j+1)) THEN
          print *, 'Attention ! rlatu2 < rlatup1 ', rlatu2(j), rlatu(j+1), j
          STOP 14
       ENDIF

       IF (rlatu(j).LE.rlatu1(j)) THEN
          print *, ' Attention ! rlatu < rlatu1 ', rlatu(j), rlatu1(j), j
          STOP 15
       ENDIF

       IF (rlatv(j).LE.rlatu2(j)) THEN
          print *, ' Attention ! rlatv < rlatu2 ', rlatv(j), rlatu2(j), j
          STOP 16
       ENDIF

       IF (rlatv(j).ge.rlatu1(j)) THEN
          print *, ' Attention ! rlatv > rlatu1 ', rlatv(j), rlatu1(j), j
          STOP 17
       ENDIF

       IF (rlatv(j).ge.rlatu(j)) THEN
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
