!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/gr_u_scal.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE gr_u_scal(nx,x_u,x_scal)
c%W%    %G%
c=======================================================================
c
c   Author:    Frederic Hourdin      original: 11/11/92
c   -------
c
c   Subject:
c   ------
c
c   Method:
c   --------
c
c   Interface:
c   ----------
c
c      Input:
c      ------
c
c      Output:
c      -------
c
c=======================================================================
      use dimens_m
      use paramet_m
      use comgeom
      IMPLICIT NONE
c-----------------------------------------------------------------------
c   Declararations:
c   ---------------


c   Arguments:
c   ----------

      INTEGER nx
      REAL x_u(ip1jmp1,nx),x_scal(ip1jmp1,nx)

c   Local:
c   ------

      INTEGER l,ij

c-----------------------------------------------------------------------

      DO l=1,nx
         DO ij=ip1jmp1,2,-1
            x_scal(ij,l)=
     s      (aireu(ij)*x_u(ij,l)+aireu(ij-1)*x_u(ij-1,l))
     s      /(aireu(ij)+aireu(ij-1))
         ENDDO
      ENDDO

      CALL SCOPY(nx*jjp1,x_scal(iip1,1),iip1,x_scal(1,1),iip1)

      RETURN
      END
