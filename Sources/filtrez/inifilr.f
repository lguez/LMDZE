module inifilr_m

  IMPLICIT NONE

  ! North:

  INTEGER jfiltnu, jfiltnv
  ! index of the last scalar line filtered in northern hemisphere

  real, pointer:: matriceun(:, :, :) ! (iim, iim, jfiltnu - 1)
  ! matrice filtre pour les champs situes sur la grille scalaire

  real, pointer:: matrinvn(:, :, :) ! (iim, iim, jfiltnu - 1)
  ! matrice filtre pour les champs situes sur la grille scalaire, pour
  ! le filtre inverse

  real, pointer:: matricevn(:, :, :) ! (iim, iim, jfiltnv)
  ! matrice filtre pour les champs situes sur la grille de V ou de Z

  ! South:

  integer jfiltsu, jfiltsv
  ! index of the first line filtered in southern hemisphere

  real, pointer:: matriceus(:, :, :) ! (iim, iim, jjm - jfiltsu + 1)
  ! matrice filtre pour les champs situes sur la grille scalaire

  real, pointer:: matrinvs(:, :, :) ! (iim, iim, jjm - jfiltsu + 1)
  ! matrice filtre pour les champs situes sur la grille scalaire, pour
  ! le filtre inverse

  real, pointer:: matricevs(:, :, :) ! (iim, iim, jjm - jfiltsv + 1)
  ! matrice filtre pour les champs situes sur la grille de V ou de Z

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

    USE dimens_m, ONLY : iim, jjm
    USE dynetat0_m, ONLY : rlatu, rlatv, xprimu, grossismx
    use inifgn_m, only: inifgn
    use inifilr_hemisph_m, only: inifilr_hemisph
    use jumble, only: new_unit
    use nr_util, only: pi, ifirstloc

    ! Local:

    REAL dlatu(jjm)
    REAL rlamda(2:iim) ! > 0, in descending order
    real eignvl(iim) ! eigenvalues sorted in descending order (<= 0)
    INTEGER j, unit
    REAL colat0 ! > 0
    integer j1 ! index of smallest positive latitude

    real eignfnu(iim, iim), eignfnv(iim, iim)
    ! eigenvectors of the discrete second derivative with respect to longitude

    !-----------------------------------------------------------

    print *, "Call sequence information: inifilr"

    CALL inifgn(eignvl, eignfnu, eignfnv)

    ! Calcul de colat0
    forall (j = 1:jjm) dlatu(j) = rlatu(j) - rlatu(j + 1)
    colat0 = min(0.5, minval(dlatu) / minval(xprimu(:iim)))
    PRINT *, 'colat0 = ', colat0

    rlamda = iim / (pi * colat0 / grossismx) / sqrt(- eignvl(2: iim))
    call new_unit(unit)
    open(unit, file = "modfrst.csv", status = "replace", action = "write") 
    write(unit, fmt = *) '"rlat (degrees)" modfrst' ! title line

   ! D\'etermination de jfilt[ns][uv] :

    j1 = ifirstloc(rlatu <= 0.)

    call inifilr_hemisph(rlatu(j1 - 1:2:- 1), colat0, rlamda, unit, eignfnv, &
         jfiltnu, matriceun, matrinvn)
    jfiltnu = j1 - jfiltnu
    matriceun = matriceun(:, :, jfiltnu - 1:1:- 1)
    matrinvn = matrinvn(:, :, jfiltnu - 1:1:- 1)

    call inifilr_hemisph(- rlatu(j1:jjm), colat0, rlamda, unit, eignfnv, &
         jfiltsu, matriceus, matrinvs)
    jfiltsu = j1 - 1 + jfiltsu

    j1 = ifirstloc(rlatv <= 0.)

    call inifilr_hemisph(rlatv(j1 - 1:1:- 1), colat0, rlamda, unit, eignfnu, &
         jfiltnv, matricevn)
    jfiltnv = j1 - jfiltnv
    matricevn = matricevn(:, :, jfiltnv:1:- 1)

    call inifilr_hemisph(- rlatv(j1:jjm), colat0, rlamda, unit, eignfnu, &
         jfiltsv, matricevs)
    jfiltsv = j1 - 1 + jfiltsv

    close(unit)
    PRINT *, 'jfiltnu =', jfiltnu
    PRINT *, 'jfiltsu =', jfiltsu
    PRINT *, 'jfiltnv =', jfiltnv
    PRINT *, 'jfiltsv =', jfiltsv

  END SUBROUTINE inifilr

end module inifilr_m
