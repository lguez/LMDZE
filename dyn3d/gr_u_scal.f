
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/gr_u_scal.F,v 1.1.1.1 2004/05/19
! 12:53:06 lmdzadmin Exp $

SUBROUTINE gr_u_scal(nx, x_u, x_scal)
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
  USE dimensions
  USE paramet_m
  USE comgeom
  IMPLICIT NONE
  ! -----------------------------------------------------------------------
  ! Declararations:
  ! ---------------


  ! Arguments:
  ! ----------

  INTEGER nx
  REAL x_u(ip1jmp1, nx), x_scal(ip1jmp1, nx)

  ! Local:
  ! ------

  INTEGER l, ij

  ! -----------------------------------------------------------------------

  DO l = 1, nx
    DO ij = ip1jmp1, 2, -1
      x_scal(ij, l) = (aireu(ij)*x_u(ij,l)+aireu(ij-1)*x_u(ij-1,l))/ &
        (aireu(ij)+aireu(ij-1))
    END DO
  END DO

  CALL scopy(nx*jjp1, x_scal(iip1,1), iip1, x_scal(1,1), iip1)

  RETURN
END SUBROUTINE gr_u_scal
