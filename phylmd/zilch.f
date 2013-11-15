!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/zilch.F,v 1.1.1.1 2004/05/19 12:53:09 lmdzadmin Exp $
!
      subroutine zilch(x,m)
c
c Zero the real array x dimensioned m.
c
      implicit none
c
      integer, intent(in):: m
      integer i
      real x(m)
      do 1 i=1,m
      x(i)= 0.0  
    1 continue
      return
      end
