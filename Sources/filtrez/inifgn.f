module inifgn_m

  use dimens_m, only: iim

  IMPLICIT NONE

  private iim

  real sddu(iim), sddv(iim) 
  ! sdd[uv] = sqrt(2 pi / iim * (derivative of the longitudinal zoom
  ! function)(rlon[uv]))

  real unsddu(iim), unsddv(iim)

contains

  SUBROUTINE inifgn(eignval_v, eignfnu, eignfnv)

    ! From LMDZ4/libf/filtrez/inifgn.F, v 1.1.1.1 2004/05/19 12:53:09

    ! Authors: H. Upadyaya, O. Sharma

    ! Computes the eigenvalues and eigenvectors of the discrete analog
    ! of the second derivative with respect to longitude.

    use acc_m, only: acc 
    USE dimens_m, ONLY: iim
    USE dynetat0_m, ONLY: xprimu, xprimv
    use numer_rec_95, only: jacobi, eigsrt

    real, intent(out):: eignval_v(:) ! (iim)
    ! eigenvalues sorted in descending order

    real, intent(out):: eignfnu(:, :), eignfnv(:, :) ! (iim, iim) eigenvectors

    ! Local:

    REAL a(iim, iim) ! second derivative, symmetric, elements are angle^{-2}

    REAL deriv_u(iim, iim), deriv_v(iim, iim) 
    ! first derivative at u and v longitudes, elements are angle^{-1}

    REAL eignval_u(iim)
    INTEGER i

    !----------------------------------------------------------------

    print *, "Call sequence information: inifgn"

    sddv = sqrt(xprimv(:iim))
    sddu = sqrt(xprimu(:iim))
    unsddu = 1. / sddu
    unsddv = 1. / sddv

    deriv_u = 0.
    deriv_u(iim, 1) = unsddu(iim) * unsddv(1)
    forall (i = 1:iim) deriv_u(i, i) = - unsddu(i) * unsddv(i)
    forall (i = 1:iim - 1) deriv_u(i, i + 1) = unsddu(i) * unsddv(i + 1)

    deriv_v = - transpose(deriv_u)

    a = matmul(deriv_v, deriv_u) ! second derivative at v longitudes
    CALL jacobi(a, eignval_v, eignfnv)
    CALL acc(eignfnv)
    CALL eigsrt(eignval_v, eignfnv)

    a = matmul(deriv_u, deriv_v) ! second derivative at u longitudes
    CALL jacobi(a, eignval_u, eignfnu)
    CALL acc(eignfnu)
    CALL eigsrt(eignval_u, eignfnu)

  END SUBROUTINE inifgn

end module inifgn_m
