!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/gr_v_scal.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE gr_v_scal(nx,x_v,x_scal)
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
      REAL x_v(ip1jm,nx),x_scal(ip1jmp1,nx)

c   Local:
c   ------

      INTEGER l,ij

c-----------------------------------------------------------------------

      DO l=1,nx
         DO ij=iip2,ip1jm
            x_scal(ij,l)=
     s      (airev(ij-iip1)*x_v(ij-iip1,l)+airev(ij)*x_v(ij,l))
     s      /(airev(ij-iip1)+airev(ij))
         ENDDO
         DO ij=1,iip1
            x_scal(ij,l)=0.
         ENDDO
         DO ij=ip1jm+1,ip1jmp1
            x_scal(ij,l)=0.
         ENDDO
      ENDDO

      RETURN
      END
