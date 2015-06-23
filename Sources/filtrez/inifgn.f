module inifgn_m

  use dimens_m, only: iim

  IMPLICIT NONE

  private iim

  real sddu(iim), sddv(iim) ! SQRT(dx / di)
  real unsddu(iim), unsddv(iim)

  real eignfnu(iim, iim), eignfnv(iim, iim)
  ! eigenfunctions of the discrete laplacian

contains

  SUBROUTINE inifgn(dv)

    ! From LMDZ4/libf/filtrez/inifgn.F, v 1.1.1.1 2004/05/19 12:53:09

    ! H. Upadyaya, O. Sharma 

    use acc_m, only: acc 
    USE dimens_m, ONLY: iim
    USE dynetat0_m, ONLY: xprimu, xprimv
    use numer_rec_95, only: jacobi, eigsrt

    real, intent(out):: dv(:) ! (iim) eigenvalues sorted in descending order

    ! Local:
    REAL, dimension(iim, iim):: a, b, c
    REAL du(iim)
    INTEGER i

    !----------------------------------------------------------------

    print *, "Call sequence information: inifgn"

    sddv = sqrt(xprimv(:iim))
    sddu = sqrt(xprimu(:iim))
    unsddu = 1. / sddu
    unsddv = 1. / sddv

    b = 0.
    b(iim, 1) = unsddu(iim) * unsddv(1)
    forall (i = 1:iim) b(i, i) = - unsddu(i) * unsddv(i)
    forall (i = 1:iim - 1) b(i, i + 1) = unsddu(i) * unsddv(i + 1)

    c = - transpose(b)

    a = matmul(c, b)
    CALL jacobi(a, dv, eignfnv)
    CALL acc(eignfnv)
    CALL eigsrt(dv, eignfnv)

    a = matmul(b, c)
    CALL jacobi(a, du, eignfnu)
    CALL acc(eignfnu)
    CALL eigsrt(du, eignfnu)

  END SUBROUTINE inifgn

end module inifgn_m
