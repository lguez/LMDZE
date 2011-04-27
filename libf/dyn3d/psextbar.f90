SUBROUTINE psextbar(ps, psexbarxy)

  ! From LMDZ4/libf/dyn3d/psextbar.F, version 1.1.1.1 2004/05/19 12:53:06
  ! Author: P. Le Van

  ! Objet : calcul des moyennes en x et en y de (pression au sol * aire
  ! variable)
  ! Cf. "psextbar.txt".

  use dimens_m
  use paramet_m
  use comgeom

  IMPLICIT NONE

  REAL, intent(in):: ps(ip1jmp1)
  real, intent(out):: psexbarxy(ip1jm)

  ! Local variables:
  real pext(ip1jmp1)
  INTEGER l, ij

  !--------------------------------------------------------

  DO ij = 1, ip1jmp1
     pext(ij) = ps(ij) * aire(ij)
  ENDDO

  DO ij = 1, ip1jm - 1
     psexbarxy(ij) = pext(ij) * alpha2(ij) + pext(ij+1) * alpha3(ij+1) &
          + pext(ij+iip1) * alpha1(ij+iip1) + pext(ij+iip2) * alpha4(ij+iip2)
  end DO

  ! Correction pour psexbarxy(iip1,j) :
  DO ij = iip1, ip1jm, iip1
     psexbarxy(ij) = psexbarxy(ij - iim)
  end DO

END SUBROUTINE psextbar
