!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/coordij.F,v 1.1.1.1 2004/05/19 12:53:05 lmdzadmin Exp $
!
      SUBROUTINE coordij(lon,lat,ilon,jlat)

c=======================================================================
c
c   calcul des coordonnees i et j de la maille scalaire dans
c   laquelle se trouve le point (lon,lat) en radian
c
c=======================================================================

      use dimens_m
      use paramet_m
      use comconst
      use comgeom
      use serre
      IMPLICIT NONE
      REAL lon,lat
      INTEGER ilon,jlat
      INTEGER i,j


      real zlon,zlat

      zlon=lon*pi/180.
      zlat=lat*pi/180.

      DO i=1,iim+1
         IF (rlonu(i).GT.zlon) THEN
            ilon=i
            GOTO 10
         ENDIF
      ENDDO
10    CONTINUE

      j=0
      DO j=1,jjm
         IF(rlatv(j).LT.zlat) THEN
            jlat=j
            GOTO 20
         ENDIF
      ENDDO
20    CONTINUE
      IF(j.EQ.0) j=jjm+1

      RETURN
      END
