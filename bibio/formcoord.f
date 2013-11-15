!
! $Header: /home/cvsroot/LMDZ4/libf/bibio/formcoord.F,v 1.1.1.1 2004/05/19 12:53:05 lmdzadmin Exp $
!
      subroutine formcoord(unit,n,x,a,rev,text)
      implicit none
      integer n,unit,ndec
      logical rev
      real x(n),a
      character*4 text

      integer i,id,i1,i2,in
      real dx,dxmin

      if(rev) then
         id=-1
         i1=n
         i2=n-1
         in=1
         write(unit,3000) text(1:1)
      else
         id=1
         i1=1
         i2=2
         in=n
      endif

      if (n.lt.2) then
         ndec=1
         write(unit,1000) text,n,x(1)*a
      else
         dxmin=abs(x(2)-x(1))
         do i=2,n-1
            dx=abs(x(i+1)-x(i))
            if (dx.lt.dxmin) dxmin=dx
         enddo

         ndec=-log10(dxmin)+2
         if(mod(n,6).eq.1) then
            write(unit,1000) text,n,x(i1)*a
            write(unit,2000) (x(i)*a,i=i2,in,id)
         else
            write(unit,1000) text,n
            write(unit,2000) (x(i)*a,i=i1,in,id)
         endif
      endif

1000  format(a4,2x,i4,' LEVELS',43x,f12.2)
2000  format(6f12.2)
c1000  format(a4,2x,i4,' LEVELS',43x,f12.<ndec>)
c2000  format(6f12.<ndec>)
3000  format('FORMAT ',a1,'REV')
      return

      end
