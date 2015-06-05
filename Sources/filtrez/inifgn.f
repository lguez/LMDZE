module inifgn_m

  use dimens_m, only: iim

  IMPLICIT NONE

  private iim

  real sddu(iim), sddv(iim) ! SQRT(dx)
  real unsddu(iim), unsddv(iim)

  real eignfnu(iim, iim), eignfnv(iim, iim)
  ! eignfn eigenfunctions of the discrete laplacian

contains

  SUBROUTINE inifgn(dv)

    ! From LMDZ4/libf/filtrez/inifgn.F, v 1.1.1.1 2004/05/19 12:53:09

    ! H.Upadyaya, O.Sharma 

    USE dimens_m, ONLY: iim
    USE dynetat0_m, ONLY: xprimu, xprimv
    use nr_util, only: pi
    use numer_rec_95, only: jacobi

    real, intent(out):: dv(iim)

    ! Local:
    REAL vec(iim, iim), vec1(iim, iim)
    REAL du(iim)
    real d(iim)
    INTEGER i, j, k, nrot

    EXTERNAL acc

    !----------------------------------------------------------------

    sddv = sqrt(xprimv(:iim))
    sddu = sqrt(xprimu(:iim))
    unsddu = 1. / sddu
    unsddv = 1. / sddv

    DO j = 1, iim
       DO i = 1, iim
          vec(i, j) = 0.
          vec1(i, j) = 0.
          eignfnv(i, j) = 0.
          eignfnu(i, j) = 0.
       END DO
    END DO

    eignfnv(1, 1) = - 1.
    eignfnv(iim, 1) = 1.
    DO i = 1, iim - 1
       eignfnv(i+1, i+1) = - 1.
       eignfnv(i, i+1) = 1.
    END DO

    DO j = 1, iim
       DO i = 1, iim
          eignfnv(i, j) = eignfnv(i, j) / (sddu(i) * sddv(j))
       END DO
    END DO

    DO j = 1, iim
       DO i = 1, iim
          eignfnu(i, j) = - eignfnv(j, i)
       END DO
    END DO

    DO j = 1, iim
       DO i = 1, iim
          vec(i, j) = 0.0
          vec1(i, j) = 0.0
          DO k = 1, iim
             vec(i, j) = vec(i, j) + eignfnu(i, k) * eignfnv(k, j)
             vec1(i, j) = vec1(i, j) + eignfnv(i, k) * eignfnu(k, j)
          END DO
       END DO
    END DO

    CALL jacobi(vec, dv, eignfnv, nrot)
    CALL acc(eignfnv, d, iim)
    CALL eigen_sort(dv, eignfnv, iim, iim)

    CALL jacobi(vec1, du, eignfnu, nrot)
    CALL acc(eignfnu, d, iim)
    CALL eigen_sort(du, eignfnu, iim, iim)

  END SUBROUTINE inifgn

end module inifgn_m
