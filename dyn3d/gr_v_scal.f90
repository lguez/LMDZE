
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/gr_v_scal.F,v 1.1.1.1 2004/05/19
! 12:53:06 lmdzadmin Exp $

SUBROUTINE gr_v_scal(nx, x_v, x_scal)
  ! %W%    %G%
  ! =======================================================================

  ! Author:    Frederic Hourdin      original: 11/11/92
  ! -------

  ! Subject:
  ! ------

  ! Method:
  ! --------

  ! Interface:
  ! ----------

  ! Input:
  ! ------

  ! Output:
  ! -------

  ! =======================================================================
  USE dimens_m
  USE paramet_m
  USE comgeom
  IMPLICIT NONE
  ! -----------------------------------------------------------------------
  ! Declararations:
  ! ---------------


  ! Arguments:
  ! ----------

  INTEGER nx
  REAL x_v(ip1jm, nx), x_scal(ip1jmp1, nx)

  ! Local:
  ! ------

  INTEGER l, ij

  ! -----------------------------------------------------------------------

  DO l = 1, nx
    DO ij = iip2, ip1jm
      x_scal(ij, l) = (airev(ij-iip1)*x_v(ij-iip1,l)+airev(ij)*x_v(ij,l))/ &
        (airev(ij-iip1)+airev(ij))
    END DO
    DO ij = 1, iip1
      x_scal(ij, l) = 0.
    END DO
    DO ij = ip1jm + 1, ip1jmp1
      x_scal(ij, l) = 0.
    END DO
  END DO

  RETURN
END SUBROUTINE gr_v_scal
