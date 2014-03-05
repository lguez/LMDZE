
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/heavyside.F,v 1.1.1.1 2004/05/19
! 12:53:06 lmdzadmin Exp $



FUNCTION heavyside(a)

  ! ...   P. Le Van  ....

  IMPLICIT NONE

  DOUBLE PRECISION heavyside, a

  IF (a<=0.) THEN
    heavyside = 0.
  ELSE
    heavyside = 1.
  END IF

  RETURN
END FUNCTION heavyside


