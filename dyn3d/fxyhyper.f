module fxyhyper_m

  IMPLICIT NONE

contains

  SUBROUTINE fxyhyper(rlatu, yprimu, rlatv, rlatu1, yprimu1, rlatu2, yprimu2, &
       rlonu, xprimu, rlonv, xprimv, xprimm025, xprimp025)

    ! From dyn3d/fxyhyper.F, version 1.1.1.1, 2004/05/19 12:53:06

    use fxhyp_m, only: fxhyp
    use fyhyp_m, only: fyhyp

    REAL, intent(out):: rlatu(:), yprimu(:) ! (jjm + 1)
    real, intent(out):: rlatv(:) ! (jjm)
    real, intent(out):: rlatu1(:), yprimu1(:), rlatu2(:), yprimu2(:) ! (jjm)
    REAL, intent(out):: rlonu(:), xprimu(:), rlonv(:), xprimv(:) ! (iim + 1)
    REAL, intent(out):: xprimm025(:) ! (iim + 1)
    REAL, intent(out):: xprimp025(:) ! (iim + 1)

    !----------------------------------------------------------

    CALL fyhyp(rlatu, yprimu, rlatv, rlatu2, yprimu2, rlatu1, yprimu1)
    CALL fxhyp(xprimm025, rlonv, xprimv, rlonu, xprimu, xprimp025)

  END SUBROUTINE fxyhyper

end module fxyhyper_m
