module inifilr_hemisph_m

  implicit none

contains

  subroutine inifilr_hemisph(rlat, colat0, rlamda, unit, eignfn, jfilt, &
       matrice, matrinv)

    ! See notes.

    USE dimens_m, ONLY : iim
    use nr_util, only: pi, ifirstloc

    real, intent(in):: rlat(:) ! (n)
    ! latitudes, in rad, in [0, pi / 2[, in strictly ascending order

    REAL, intent(in):: colat0 ! > 0
    REAL, intent(in):: rlamda(2:) ! (2:iim) > 0, in descending order
    integer, intent(in):: unit

    real, intent(in):: eignfn(:, :) ! (iim, iim)
    ! eigenvectors of the discrete second derivative with respect to longitude

    integer, intent(out):: jfilt
    real, pointer:: matrice(:, :, :) ! (iim, iim, n - jfilt + 1) matrice filtre

    real, pointer, optional:: matrinv(:, :, :) ! (iim, iim, n - jfilt + 1) 
    ! matrice pour le filtre inverse

    ! Local:

    integer n, i, j
    REAL eignft(iim, iim)

    ! Index of the mode from where modes are filtered:
    integer, allocatable:: modfrst(:) ! (jfilt:n) in {2, ..., iim}

    ! Filtering coefficients (lamda_max * cos(rlat) / lamda):
    real coefil(iim)

    !-----------------------------------------------------------

    n = size(rlat)
    jfilt = ifirstloc(cos(rlat) < min(colat0, 1. / rlamda(iim)))
    allocate(modfrst(jfilt:n))

    DO j = jfilt, n
       modfrst(j) = ifirstloc(rlamda < 1. / cos(rlat(j)), my_lbound = 2)
       write(unit, fmt = *) rlat(j) / pi * 180., modfrst(j)
    END DO

    allocate(matrice(iim, iim, n - jfilt + 1))
    if (present(matrinv)) allocate(matrinv(iim, iim, n - jfilt + 1))

    DO j = jfilt, n
       DO i = modfrst(j), iim
          coefil(i) = rlamda(i) * cos(rlat(j)) - 1.
       end DO

       eignft(:modfrst(j) - 1, :) = 0.

       forall (i = modfrst(j):iim) eignft(i, :) = eignfn(:, i) * coefil(i)
       matrice(:, :, j - jfilt + 1) = matmul(eignfn, eignft)

       if (present(matrinv)) then
          forall (i = modfrst(j):iim) eignft(i, :) = eignfn(:, i) * coefil(i) &
               / (1. + coefil(i))
          matrinv(:, :, j - jfilt + 1) = matmul(eignfn, eignft)
       end if
    END DO

  end subroutine inifilr_hemisph

end module inifilr_hemisph_m
