!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/heavyside.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
c
c
       FUNCTION heavyside(a)

c      ...   P. Le Van  ....
c
       IMPLICIT NONE

       REAL*8 heavyside , a

       IF ( a.LE.0. )  THEN
         heavyside = 0.
       ELSE
         heavyside = 1.
       ENDIF

       RETURN
       END


