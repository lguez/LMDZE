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
    use inifgn_m, only: inifgn
    use jumble, only: new_unit
    use nr_util, only: pi

    ! Local:
    REAL dlatu(jjm)
    REAL rlamda(2: iim)
    real eignvl(iim) ! eigenvalues sorted in descending order
    REAL cof
    INTEGER i, j, k, unit
    REAL colat0 ! > 0
    REAL eignft(iim, iim), coff

    real eignfnu(iim, iim), eignfnv(iim, iim)
    ! eigenfunctions of the discrete laplacian

    ! Filtering coefficients (lamda_max * cos(rlat) / lamda):
    real coefilu(iim, jjm), coefilv(iim, jjm)
    real coefilu2(iim, jjm), coefilv2(iim, jjm)

    ! Index of the mode from where modes are filtered:
    integer, allocatable:: modfrstnu(:), modfrstsu(:)
    integer, allocatable:: modfrstnv(:), modfrstsv(:)

    !-----------------------------------------------------------

    print *, "Call sequence information: inifilr"

    CALL inifgn(eignvl, eignfnu, eignfnv)

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

    ! Determination de coefilu, coefilv, modfrst[ns][uv]:

    allocate(modfrstnu(2:jfiltnu), modfrstsu(jfiltsu:jjm))
    allocate(modfrstnv(jfiltnv), modfrstsv(jfiltsv:jjm))
    coefilu = 0.
    coefilv = 0.
    coefilu2 = 0.
    coefilv2 = 0.

    DO j = 2, jfiltnu
       modfrstnu(j) = 2
       do while (rlamda(modfrstnu(j)) * cos(rlatu(j)) >= 1. &
            .and. modfrstnu(j) <= iim - 1)
          modfrstnu(j) = modfrstnu(j) + 1
       end do

       if (rlamda(modfrstnu(j)) * cos(rlatu(j)) < 1.) then
          DO k = modfrstnu(j), iim
             cof = rlamda(k) * cos(rlatu(j))
             coefilu(k, j) = cof - 1.
             coefilu2(k, j) = cof**2 - 1.
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
          DO k = modfrstnv(j), iim
             cof = rlamda(k) * cos(rlatv(j))
             coefilv(k, j) = cof - 1.
             coefilv2(k, j) = cof**2 - 1.
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
          DO k = modfrstsu(j), iim
             cof = rlamda(k) * cos(rlatu(j))
             coefilu(k, j) = cof - 1.
             coefilu2(k, j) = cof**2 - 1.
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
          DO k = modfrstsv(j), iim
             cof = rlamda(k) * cos(rlatv(j))
             coefilv(k, j) = cof - 1.
             coefilv2(k, j) = cof**2 - 1.
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
       DO i = 1, iim
          IF (i < modfrstnu(j)) then
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
          IF (i < modfrstsu(j)) then
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
          IF (i < modfrstnv(j)) then
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
          IF (i < modfrstsv(j)) then
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
          IF (i < modfrstnu(j)) then
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
          IF (i < modfrstsu(j)) then
             coff = 0.
          else
             coff = coefilu(i, j) / (1. + coefilu(i, j))
          end IF
          eignft(i, :) = eignfnv(:, i) * coff
       END DO
       matrinvs(:, :, j) = matmul(eignfnv, eignft)
    END DO

  END SUBROUTINE inifilr

end module inifilr_m
