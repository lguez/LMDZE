!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/initial0.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
      SUBROUTINE initial0(n,x)
      IMPLICIT NONE
      INTEGER n,i
      REAL x(n)
      DO 10 i=1,n
         x(i)=0.
10    CONTINUE
      RETURN
      END
