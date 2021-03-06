module inifgn_m

  IMPLICIT NONE

  real, allocatable, protected:: sddu(:), sddv(:) ! (iim) 
  ! sdd[uv] = sqrt(2 pi / iim * (derivative of the longitudinal zoom
  ! function)(rlon[uv]))

  real, allocatable, protected:: unsddu(:), unsddv(:) ! (iim)

contains

  SUBROUTINE inifgn(eignval_v, eignfnu, eignfnv)

    ! From LMDZ4/libf/filtrez/inifgn.F, v 1.1.1.1 2004/05/19 12:53:09

    ! Authors: H. Upadyaya, O. Sharma

    ! Computes the eigenvalues and eigenvectors of the discrete analog
    ! of the second derivative with respect to longitude.

    USE dimensions, ONLY: iim
    USE dynetat0_m, ONLY: xprimu, xprimv
    use jumble, only: new_unit
    use numer_rec_95, only: jacobi, eigsrt

    real, intent(out):: eignval_v(:) ! (iim)
    ! eigenvalues sorted in descending order

    real, intent(out):: eignfnu(:, :), eignfnv(:, :) ! (iim, iim) eigenvectors

    ! Local:

    REAL delta(iim, iim) ! second derivative, symmetric, elements are angle^{-2}

    REAL deriv_u(iim, iim), deriv_v(iim, iim) 
    ! first derivative at u and v longitudes, elements are angle^{-1}

    REAL eignval_u(iim)
    INTEGER i, unit

    !----------------------------------------------------------------

    print *, "Call sequence information: inifgn"

    allocate(sddu(iim), sddv(iim))
    allocate(unsddu(iim), unsddv(iim))

    sddv = sqrt(xprimv(:iim))
    sddu = sqrt(xprimu(:iim))
    unsddu = 1. / sddu
    unsddv = 1. / sddv

    deriv_u = 0.
    deriv_u(iim, 1) = unsddu(iim) * unsddv(1)
    forall (i = 1:iim) deriv_u(i, i) = - unsddu(i) * unsddv(i)
    forall (i = 1:iim - 1) deriv_u(i, i + 1) = unsddu(i) * unsddv(i + 1)

    deriv_v = - transpose(deriv_u)

    delta = matmul(deriv_v, deriv_u) ! second derivative at v longitudes
    CALL jacobi(delta, eignval_v, eignfnv)
    CALL eigsrt(eignval_v, eignfnv)

    delta = matmul(deriv_u, deriv_v) ! second derivative at u longitudes
    CALL jacobi(delta, eignval_u, eignfnu)
    CALL eigsrt(eignval_u, eignfnu)

    call new_unit(unit)
    open(unit, file = "inifgn_out.txt", status = "replace", action = "write") 
    write(unit, fmt = *) '"eignval_v"', eignval_v
    close(unit)

  END SUBROUTINE inifgn

end module inifgn_m
