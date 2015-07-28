module inifilr_m

  IMPLICIT NONE

  ! North:

  INTEGER jfiltnu, jfiltnv
  ! index of the last scalar line filtered in northern hemisphere

  real, allocatable:: matriceun(:, :, :), matrinvn(:, :, :) 
  ! (iim, iim, 2:jfiltnu)

  real, allocatable:: matricevn(:, :, :) ! (iim, iim, jfiltnv)

  ! South:

  integer jfiltsu, jfiltsv
  ! index of the first line filtered in southern hemisphere

  real, allocatable:: matriceus(:, :, :), matrinvs(:, :, :)
  ! (iim, iim, jfiltsu:jjm)

  real, allocatable:: matricevs(:, :, :) ! (iim, iim, jfiltsv:jjm)

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
    use jumble, only: new_unit
    use nr_util, only: pi

    ! Local:

    REAL dlatu(jjm)
    REAL rlamda(2: iim)
    real eignvl(iim) ! eigenvalues sorted in descending order (<= 0)
    INTEGER i, j, unit
    REAL colat0 ! > 0
    REAL eignft(iim, iim)

    real eignfnu(iim, iim), eignfnv(iim, iim)
    ! eigenvectors of the discrete second derivative with respect to longitude

    ! Filtering coefficients (lamda_max * cos(rlat) / lamda):
    real, allocatable:: coefilnu(:, :) ! (iim, 2:jfiltnu)
    real, allocatable:: coefilsu(:, :) ! (iim, jfiltsu:jjm)
    real, allocatable:: coefilnv(:, :) ! (iim, jfiltnv)
    real, allocatable:: coefilsv(:, :) ! (iim, jfiltsv:jjm)

    ! Index of the mode from where modes are filtered:
    integer, allocatable:: modfrstnu(:) ! (2:jfiltnu)
    integer, allocatable:: modfrstsu(:) ! (jfiltsu:jjm)
    integer, allocatable:: modfrstnv(:) ! (jfiltnv)
    integer, allocatable:: modfrstsv(:) ! (jfiltsv:jjm)

    !-----------------------------------------------------------

    print *, "Call sequence information: inifilr"

    CALL inifgn(eignvl, eignfnu, eignfnv)

    ! Calcul de colat0
    forall (j = 1:jjm) dlatu(j) = rlatu(j) - rlatu(j + 1)
    colat0 = min(0.5, minval(dlatu) / minval(xprimu(:iim)))
    PRINT *, 'colat0 = ', colat0

    rlamda = iim / (pi * colat0 / grossismx) / sqrt(- eignvl(2: iim))

    ! Determination de jfiltnu, jfiltsu, jfiltnv, jfiltsv

    jfiltnu = (jjm + 1) / 2
    do while (cos(rlatu(jfiltnu)) >= colat0 &
         .or. rlamda(iim) * cos(rlatu(jfiltnu)) >= 1.)
       jfiltnu = jfiltnu - 1
    end do

    jfiltsu = jjm / 2 + 2
    do while (cos(rlatu(jfiltsu)) >= colat0 &
         .or. rlamda(iim) * cos(rlatu(jfiltsu)) >= 1.)
       jfiltsu = jfiltsu + 1
    end do

    jfiltnv = jjm / 2
    do while ((cos(rlatv(jfiltnv)) >= colat0 &
         .or. rlamda(iim) * cos(rlatv(jfiltnv)) >= 1.) .and. jfiltnv >= 2)
       jfiltnv = jfiltnv - 1
    end do

    if (cos(rlatv(jfiltnv)) >= colat0 &
         .or. rlamda(iim) * cos(rlatv(jfiltnv)) >= 1.) then
       ! {jfiltnv == 1}
       PRINT *, 'Could not find jfiltnv.'
       STOP 1
    END IF

    jfiltsv = (jjm + 1)/ 2 + 1
    do while ((cos(rlatv(jfiltsv)) >= colat0 &
         .or. rlamda(iim) * cos(rlatv(jfiltsv)) >= 1.) .and. jfiltsv <= jjm - 1)
       jfiltsv = jfiltsv + 1
    end do

    IF (cos(rlatv(jfiltsv)) >= colat0 &
         .or. rlamda(iim) * cos(rlatv(jfiltsv)) >= 1.) THEN
       ! {jfiltsv == jjm}
       PRINT *, 'Could not find jfiltsv.'
       STOP 1
    END IF

    PRINT *, 'jfiltnu =', jfiltnu
    PRINT *, 'jfiltsu =', jfiltsu
    PRINT *, 'jfiltnv =', jfiltnv
    PRINT *, 'jfiltsv =', jfiltsv

    ! D\'etermination de coefil[ns][uv], modfrst[ns][uv]:

    allocate(modfrstnu(2:jfiltnu), modfrstsu(jfiltsu:jjm))
    allocate(modfrstnv(jfiltnv), modfrstsv(jfiltsv:jjm))
    allocate(coefilnu(iim, 2:jfiltnu), coefilsu(iim, jfiltsu:jjm))
    allocate(coefilnv(iim, jfiltnv), coefilsv(iim, jfiltsv:jjm))

    coefilnu = 0.
    coefilnv = 0.
    coefilsu = 0.
    coefilsv = 0.

    DO j = 2, jfiltnu
       modfrstnu(j) = 2
       do while (rlamda(modfrstnu(j)) * cos(rlatu(j)) >= 1. &
            .and. modfrstnu(j) <= iim - 1)
          modfrstnu(j) = modfrstnu(j) + 1
       end do

       if (rlamda(modfrstnu(j)) * cos(rlatu(j)) < 1.) then
          DO i = modfrstnu(j), iim
             coefilnu(i, j) = rlamda(i) * cos(rlatu(j)) - 1.
          end DO
       end if
    END DO

    DO j = 1, jfiltnv
       modfrstnv(j) = 2
       do while (rlamda(modfrstnv(j)) * cos(rlatv(j)) >= 1. &
            .and. modfrstnv(j) <= iim - 1)
          modfrstnv(j) = modfrstnv(j) + 1
       end do

       if (rlamda(modfrstnv(j)) * cos(rlatv(j)) < 1.) then
          DO i = modfrstnv(j), iim
             coefilnv(i, j) = rlamda(i) * cos(rlatv(j)) - 1.
          end DO
       end if
    end DO

    DO j = jfiltsu, jjm
       modfrstsu(j) = 2
       do while (rlamda(modfrstsu(j)) * cos(rlatu(j)) >= 1. &
            .and. modfrstsu(j) <= iim - 1)
          modfrstsu(j) = modfrstsu(j) + 1
       end do

       if (rlamda(modfrstsu(j)) * cos(rlatu(j)) < 1.) then
          DO i = modfrstsu(j), iim
             coefilsu(i, j) = rlamda(i) * cos(rlatu(j)) - 1.
          end DO
       end if
    end DO

    DO j = jfiltsv, jjm
       modfrstsv(j) = 2
       do while (rlamda(modfrstsv(j)) * cos(rlatv(j)) >= 1. &
            .and. modfrstsv(j) <= iim - 1)
          modfrstsv(j) = modfrstsv(j) + 1
       end do

       if (rlamda(modfrstsv(j)) * cos(rlatv(j)) < 1.) then
          DO i = modfrstsv(j), iim
             coefilsv(i, j) = rlamda(i) * cos(rlatv(j)) - 1.
          end DO
       end if
    END DO

    call new_unit(unit)
    open(unit, file = "inifilr_out.txt", status = "replace", action = "write") 
    write(unit, fmt = *) '"EIGNVL"', eignvl
    write(unit, fmt = *) '"modfrstnu"', modfrstnu
    write(unit, fmt = *) '"modfrstsu"', modfrstsu
    write(unit, fmt = *) '"modfrstnv"', modfrstnv
    write(unit, fmt = *) '"modfrstsv"', modfrstsv
    close(unit)

    allocate(matriceun(iim, iim, 2:jfiltnu), matrinvn(iim, iim, 2:jfiltnu))
    allocate(matricevn(iim, iim, jfiltnv))
    allocate(matricevs(iim, iim, jfiltsv:jjm))
    allocate(matriceus(iim, iim, jfiltsu:jjm), matrinvs(iim, iim, jfiltsu:jjm))

    ! Calcul de la matrice filtre 'matriceu' pour les champs situes
    ! sur la grille scalaire

    DO j = 2, jfiltnu
       eignft(:modfrstnu(j) - 1, :) = 0.
       forall (i = modfrstnu(j):iim) eignft(i, :) = eignfnv(:, i) &
            * coefilnu(i, j)
       matriceun(:, :, j) = matmul(eignfnv, eignft)
    END DO

    DO j = jfiltsu, jjm
       eignft(:modfrstsu(j) - 1, :) = 0.
       forall (i = modfrstsu(j):iim) eignft(i, :) = eignfnv(:, i) &
            * coefilsu(i, j)
       matriceus(:, :, j) = matmul(eignfnv, eignft)
    END DO

    ! Calcul de la matrice filtre 'matricev' pour les champs situes
    ! sur la grille de V ou de Z

    DO j = 1, jfiltnv
       eignft(:modfrstnv(j) - 1, :) = 0.
       forall (i = modfrstnv(j): iim) eignft(i, :) = eignfnu(:, i) &
            * coefilnv(i, j)
       matricevn(:, :, j) = matmul(eignfnu, eignft)
    END DO

    DO j = jfiltsv, jjm
       eignft(:modfrstsv(j) - 1, :) = 0.
       forall (i = modfrstsv(j):iim) eignft(i, :) = eignfnu(:, i) &
            * coefilsv(i, j)
       matricevs(:, :, j) = matmul(eignfnu, eignft)
    END DO

    ! Calcul de la matrice filtre 'matrinv' pour les champs situes
    ! sur la grille scalaire , pour le filtre inverse

    DO j = 2, jfiltnu
       eignft(:modfrstnu(j) - 1, :) = 0.
       forall (i = modfrstnu(j):iim) eignft(i, :) = eignfnv(:, i) &
            * coefilnu(i, j) / (1. + coefilnu(i, j))
       matrinvn(:, :, j) = matmul(eignfnv, eignft)
    END DO

    DO j = jfiltsu, jjm
       eignft(:modfrstsu(j) - 1, :) = 0.
       forall (i = modfrstsu(j):iim) eignft(i, :) = eignfnv(:, i) &
            * coefilsu(i, j) / (1. + coefilsu(i, j))
       matrinvs(:, :, j) = matmul(eignfnv, eignft)
    END DO

  END SUBROUTINE inifilr

end module inifilr_m
