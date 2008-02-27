!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/ismin.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
      FUNCTION ismin(n,sx,incx)
c
      IMPLICIT NONE
c
      integer n,i,incx,ismin,ix
      real sx((n-1)*incx+1),sxmin
c
      ix=1
      ismin=1
      sxmin=sx(1)
      DO i=1,n-1
         ix=ix+incx
         if(sx(ix).lt.sxmin) then
             sxmin=sx(ix)
             ismin=i+1
         endif
      ENDDO
c
      return
      end
C
