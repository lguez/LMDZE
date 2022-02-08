module inifilr_m

  IMPLICIT NONE

  INTEGER jfiltnu, jfiltnv
  ! index of the last line filtered in northern hemisphere at rlat[uv]
  ! latitudes

  integer jfiltsu, jfiltsv
  ! index of the first line filtered in southern hemisphere at
  ! rlat[uv] latitudes

  ! Filtre pour les champs situes sur la grille scalaire (longitudes
  ! rlonv, latitudes rlatu) :
  real, pointer:: matriceun(:, :, :) ! (iim, iim, jfiltnu - 1)
  real, pointer:: matriceus(:, :, :) ! (iim, iim, jjm - jfiltsu + 1)

  ! Filtre pour les champs situes sur la grille scalaire (longitudes
  ! rlonv, latitudes rlatu), pour le filtre inverse :
  real, pointer:: matrinvn(:, :, :) ! (iim, iim, jfiltnu - 1)
  real, pointer:: matrinvs(:, :, :) ! (iim, iim, jjm - jfiltsu + 1)

  ! Filtre pour les champs situes sur la grille de la vorticit\'e
  ! (longitudes rlonu, latitudes rlatv)
  real, pointer:: matricevn(:, :, :) ! (iim, iim, jfiltnv) matrice
  real, pointer:: matricevs(:, :, :) ! (iim, iim, jjm - jfiltsv + 1)

contains

  SUBROUTINE inifilr

    ! From filtrez/inifilr.F, version 1.1.1.1, 2004/05/19 12:53:09
    ! H. Upadhyaya, O. Sharma

    ! This procedure computes the filtering coefficients for scalar
    ! lines and meridional wind v lines. The modes are filtered from
    ! modfrst to iim. We filter all those latitude lines where coefil
    ! < 1. No filtering at poles. colat0 is to be used when alpha
    ! (stretching coefficient) is set equal to zero for the regular
    ! grid case.

    ! Libraries:
    use jumble, only: new_unit, pi, ifirstloc, assert

    USE dimensions, ONLY: iim, jjm
    USE dynetat0_m, ONLY: rlatu, rlatv, xprimu
    USE dynetat0_chosen_m, ONLY: grossismx
    use inifgn_m, only: inifgn
    use inifilr_hemisph_m, only: inifilr_hemisph

    ! Local:

    REAL dlatu(jjm)
    REAL rlamda(2:iim) ! > 0, in descending order
    real eignvl(iim) ! eigenvalues (<= 0) sorted in descending order
    INTEGER j, unit
    REAL colat0 ! > 0
    integer j1 ! index of negative latitude closest to the equator

    real eignfnu(iim, iim), eignfnv(iim, iim)
    ! eigenvectors of the discrete second derivative with respect to
    ! longitude, at rlon[uv] longitudes

    !-----------------------------------------------------------

    print *, "Call sequence information: inifilr"

    CALL inifgn(eignvl, eignfnu, eignfnv)

    ! Calcul de colat0
    forall (j = 1:jjm) dlatu(j) = rlatu(j) - rlatu(j + 1)
    colat0 = min(0.5, minval(dlatu) / minval(xprimu(:iim)))
    PRINT *, 'colat0 = ', colat0

    rlamda = iim / pi / colat0 * grossismx / sqrt(- eignvl(2: iim))
    print *, "1 / rlamda(iim) = ", 1. / rlamda(iim)
    ! This is demonstrated in the notes but just to be sure:
    call assert(rlamda(iim) * colat0 >= 1. - 2. * epsilon(0.), &
         "inifilr rlamda(iim) * colat0")

    call new_unit(unit)
    open(unit, file = "modfrst.csv", status = "replace", action = "write") 
    write(unit, fmt = *) '"rlat (degrees)" modfrst' ! title line

    j1 = ifirstloc(rlatu <= 0.)

    call inifilr_hemisph(rlatu(j1 - 1:2:- 1), rlamda, unit, eignfnv, jfiltnu, &
         matriceun, matrinvn)
    jfiltnu = j1 - jfiltnu
    matriceun = matriceun(:, :, jfiltnu - 1:1:- 1)
    matrinvn = matrinvn(:, :, jfiltnu - 1:1:- 1)

    call inifilr_hemisph(- rlatu(j1:jjm), rlamda, unit, eignfnv, jfiltsu, &
         matriceus, matrinvs)
    jfiltsu = j1 - 1 + jfiltsu

    j1 = ifirstloc(rlatv <= 0.)

    call inifilr_hemisph(rlatv(j1 - 1:1:- 1), rlamda, unit, eignfnu, jfiltnv, &
         matricevn)
    jfiltnv = j1 - jfiltnv
    matricevn = matricevn(:, :, jfiltnv:1:- 1)

    call inifilr_hemisph(- rlatv(j1:jjm), rlamda, unit, eignfnu, jfiltsv, &
         matricevs)
    jfiltsv = j1 - 1 + jfiltsv

    close(unit)
    PRINT *, 'jfiltnu =', jfiltnu
    PRINT *, 'jfiltsu =', jfiltsu
    PRINT *, 'jfiltnv =', jfiltnv
    PRINT *, 'jfiltsv =', jfiltsv

  END SUBROUTINE inifilr

end module inifilr_m
