module inifilr_m

  IMPLICIT NONE

  INTEGER jfiltnu, jfiltsu, jfiltnv, jfiltsv
  ! jfiltn index of the last scalar line filtered in NH
  ! jfilts index of the first line filtered in SH

  ! North:
  real, allocatable:: matriceun(:, :, :), matrinvn(:, :, :) 
  ! (iim, iim, 2:jfiltnu)

  real, allocatable:: matricevn(:, :, :) ! (iim, iim, jfiltnv)

  ! South:
  real, allocatable:: matriceus(:, :, :), matrinvs(:, :, :)
  ! (iim, iim, jfiltsu:jjm)

  real, allocatable:: matricevs(:, :, :) ! (iim, iim, jfiltsv:jjm)

contains

  SUBROUTINE inifilr

    ! From filtrez/inifilr.F, version 1.1.1.1 2004/05/19 12:53:09
    ! H. Upadhyaya, O. Sharma

    ! This routine computes the eigenfunctions of the laplacian on the
    ! stretched grid, and the filtering coefficients. The modes are
    ! filtered from modfrst to iim.

    USE dimens_m, ONLY : iim, jjm
    USE dynetat0_m, ONLY : rlatu, rlatv, xprimu, grossismx
    use inifgn_m, only: inifgn, eignfnu, eignfnv
    use jumble, only: new_unit
    use nr_util, only: pi

    ! Local:
    REAL dlatu(jjm)
    REAL rlamda(2: iim)
    real eignvl(iim) ! eigenvalues
    REAL cof
    INTEGER i, j, k, kf, unit
    REAL colat0
    REAL eignft(iim, iim), coff

    ! Filtering coefficients (lamda_max * cos(rlat) / lamda):
    real coefilu(iim, jjm), coefilv(iim, jjm)
    real coefilu2(iim, jjm), coefilv2(iim, jjm)

    integer modfrstu(jjm), modfrstv(jjm)
    ! index of the mode from where modes are filtered

    !-----------------------------------------------------------

    print *, "Call sequence information: inifilr"

    CALL inifgn(eignvl)

    call new_unit(unit)
    open(unit, file = "eignvl.txt", status = "replace", action = "write") 
    write(unit, fmt = *) EIGNVL
    close(unit)

    ! compute eigenvalues and eigenfunctions
    ! compute the filtering coefficients for scalar lines and
    ! meridional wind v-lines
    ! we filter all those latitude lines where coefil < 1
    ! NO FILTERING AT POLES
    ! colat0 is to be used when alpha (stretching coefficient)
    ! is set equal to zero for the regular grid case

    ! Calcul de colat0
    forall (j = 1:jjm) dlatu(j) = rlatu(j) - rlatu(j + 1)
    colat0 = min(0.5, minval(dlatu) / minval(xprimu(:iim)))
    PRINT *, 'colat0 = ', colat0

    rlamda = iim / (pi * colat0 / grossismx) / sqrt(abs(eignvl(2: iim)))
    coefilu = 0.
    coefilv = 0.
    coefilu2 = 0.
    coefilv2 = 0.

    ! Determination de jfiltnu, jfiltnv, jfiltsu, jfiltsv

    DO j = 2, jjm / 2 + 1
       IF (cos(rlatu(j)) / colat0 < 1. &
            .and. rlamda(iim) * cos(rlatu(j)) < 1.) jfiltnu = j

       IF (cos(rlatu(jjm - j + 2)) / colat0 < 1. &
            .and. rlamda(iim) * cos(rlatu(jjm - j + 2)) < 1.) &
            jfiltsu = jjm - j + 2
    END DO

    DO j = 1, jjm / 2
       IF (cos(rlatv(j)) / colat0 < 1. .and. rlamda(iim) * cos(rlatv(j)) < 1.) &
            jfiltnv = j

       IF (cos(rlatv(jjm - j + 1)) / colat0 < 1. .and. rlamda(iim) &
            * cos(rlatv(jjm - j + 1)) < 1.) jfiltsv = jjm - j + 1
    END DO

    IF (jfiltnu <= 0) jfiltnu = 1
    IF (jfiltnu > jjm / 2 + 1) THEN
       PRINT *, 'jfiltnu en dehors des valeurs acceptables ', jfiltnu
       STOP 1
    END IF

    IF (jfiltsu <= 0) jfiltsu = 1
    IF (jfiltsu > jjm + 1) THEN
       PRINT *, 'jfiltsu en dehors des valeurs acceptables ', jfiltsu
       STOP 1
    END IF

    IF (jfiltnv <= 0) jfiltnv = 1
    IF (jfiltnv > jjm / 2) THEN
       PRINT *, 'jfiltnv en dehors des valeurs acceptables ', jfiltnv
       STOP 1
    END IF

    IF (jfiltsv <= 0) jfiltsv = 1
    IF (jfiltsv > jjm) THEN
       PRINT *, 'jfiltsv en dehors des valeurs acceptables ', jfiltsv
       STOP 1
    END IF

    PRINT *, 'jfiltnv jfiltsv jfiltnu jfiltsu ', jfiltnv, jfiltsv, jfiltnu, &
         jfiltsu

    ! Determination de coefilu, coefilv, n=modfrstu, modfrstv

    DO j = 1, jjm
       modfrstu(j) = iim
       modfrstv(j) = iim
    END DO

    DO j = 2, jfiltnu
       DO k = 2, iim
          IF (rlamda(k) * cos(rlatu(j)) < 1.) exit
       end DO
       if (k == iim + 1) cycle
       modfrstu(j) = k

       kf = modfrstu(j)
       DO k = kf, iim
          cof = rlamda(k) * cos(rlatu(j))
          coefilu(k, j) = cof - 1.
          coefilu2(k, j) = cof**2 - 1.
       end DO
    END DO

    DO j = 1, jfiltnv
       DO k = 2, iim
          IF (rlamda(k) * cos(rlatv(j)) < 1.) exit
       end DO
       if (k == iim + 1) cycle
       modfrstv(j) = k

       kf = modfrstv(j)
       DO k = kf, iim
          cof = rlamda(k) * cos(rlatv(j))
          coefilv(k, j) = cof - 1.
          coefilv2(k, j) = cof**2 - 1.
       end DO
    end DO

    DO j = jfiltsu, jjm
       DO k = 2, iim
          IF (rlamda(k) * cos(rlatu(j)) < 1.) exit
       end DO
       if (k == iim + 1) cycle
       modfrstu(j) = k

       kf = modfrstu(j)
       DO k = kf, iim
          cof = rlamda(k) * cos(rlatu(j))
          coefilu(k, j) = cof - 1.
          coefilu2(k, j) = cof**2 - 1.
       end DO
    end DO

    DO j = jfiltsv, jjm
       DO k = 2, iim
          IF (rlamda(k) * cos(rlatv(j)) < 1.) exit
       end DO
       if (k == iim + 1) cycle
       modfrstv(j) = k

       kf = modfrstv(j)
       DO k = kf, iim
          cof = rlamda(k) * cos(rlatv(j))
          coefilv(k, j) = cof - 1.
          coefilv2(k, j) = cof**2 - 1.
       end DO
    END DO

    IF (jfiltnv>=jjm / 2 .OR. jfiltnu>=jjm / 2) THEN
       IF (jfiltnv == jfiltsv) jfiltsv = 1 + jfiltnv
       IF (jfiltnu == jfiltsu) jfiltsu = 1 + jfiltnu

       PRINT *, 'jfiltnv jfiltsv jfiltnu jfiltsu', jfiltnv, jfiltsv, jfiltnu, &
            jfiltsu
    END IF

    PRINT *, 'Modes premiers v '
    PRINT 334, modfrstv
    PRINT *, 'Modes premiers u '
    PRINT 334, modfrstu

    allocate(matriceun(iim, iim, 2:jfiltnu), matrinvn(iim, iim, 2:jfiltnu))
    allocate(matricevn(iim, iim, jfiltnv))
    allocate(matricevs(iim, iim, jfiltsv:jjm))
    allocate(matriceus(iim, iim, jfiltsu:jjm), matrinvs(iim, iim, jfiltsu:jjm))

    ! Calcul de la matrice filtre 'matriceu' pour les champs situes
    ! sur la grille scalaire

    DO j = 2, jfiltnu
       DO i = 1, iim
          IF (i < modfrstu(j)) then
             coff = 0.
          else
             coff = coefilu(i, j)
          end IF
          eignft(i, :) = eignfnv(:, i) * coff
       END DO
       matriceun(:, :, j) = matmul(eignfnv, eignft)
    END DO

    DO j = jfiltsu, jjm
       DO i = 1, iim
          IF (i < modfrstu(j)) then
             coff = 0.
          else
             coff = coefilu(i, j)
          end IF
          eignft(i, :) = eignfnv(:, i) * coff
       END DO
       matriceus(:, :, j) = matmul(eignfnv, eignft)
    END DO

    ! Calcul de la matrice filtre 'matricev' pour les champs situes
    ! sur la grille de V ou de Z

    DO j = 1, jfiltnv
       DO i = 1, iim
          IF (i < modfrstv(j)) then
             coff = 0.
          else
             coff = coefilv(i, j)
          end IF
          eignft(i, :) = eignfnu(:, i) * coff
       END DO
       matricevn(:, :, j) = matmul(eignfnu, eignft)
    END DO

    DO j = jfiltsv, jjm
       DO i = 1, iim
          IF (i < modfrstv(j)) then
             coff = 0.
          else
             coff = coefilv(i, j)
          end IF
          eignft(i, :) = eignfnu(:, i) * coff
       END DO
       matricevs(:, :, j) = matmul(eignfnu, eignft)
    END DO

    ! Calcul de la matrice filtre 'matrinv' pour les champs situes
    ! sur la grille scalaire , pour le filtre inverse

    DO j = 2, jfiltnu
       DO i = 1, iim
          IF (i < modfrstu(j)) then
             coff = 0.
          else
             coff = coefilu(i, j) / (1. + coefilu(i, j))
          end IF
          eignft(i, :) = eignfnv(:, i) * coff
       END DO
       matrinvn(:, :, j) = matmul(eignfnv, eignft)
    END DO

    DO j = jfiltsu, jjm
       DO i = 1, iim
          IF (i < modfrstu(j)) then
             coff = 0.
          else
             coff = coefilu(i, j) / (1. + coefilu(i, j))
          end IF
          eignft(i, :) = eignfnv(:, i) * coff
       END DO
       matrinvs(:, :, j) = matmul(eignfnv, eignft)
    END DO

334 FORMAT (1X, 24I3)

  END SUBROUTINE inifilr

end module inifilr_m
