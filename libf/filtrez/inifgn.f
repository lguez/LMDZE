!
! $Header: /home/cvsroot/LMDZ4/libf/filtrez/inifgn.F,v 1.1.1.1 2004/05/19 12:53:09 lmdzadmin Exp $
!
      SUBROUTINE inifgn(dv)
c  
c    ...  H.Upadyaya , O.Sharma  ... 
c
      use dimens_m
      use paramet_m
      use comgeom
      use serre
      IMPLICIT NONE
c

c
      REAL vec(iim,iim),vec1(iim,iim)
      REAL dlonu(iim),dlonv(iim)
      REAL du(iim),dv(iim),d(iim)
      REAL pi
      INTEGER i,j,k,imm1,nrot
C
      include "coefils.h"
c
      EXTERNAL SSUM, acc, jacobi
CC      EXTERNAL eigen
      REAL SSUM
c

      imm1  = iim -1
      pi = 2.* ASIN(1.)
C
      DO 5 i=1,iim
       dlonu(i)=  xprimu( i )
       dlonv(i)=  xprimv( i )
   5  CONTINUE

      DO 12 i=1,iim
      sddv(i)   = SQRT(dlonv(i))
      sddu(i)   = SQRT(dlonu(i))
      unsddu(i) = 1./sddu(i)
      unsddv(i) = 1./sddv(i)
  12  CONTINUE
C
      DO 17 j=1,iim
      DO 17 i=1,iim
      vec(i,j)     = 0.
      vec1(i,j)    = 0.
      eignfnv(i,j) = 0.
      eignfnu(i,j) = 0.
  17  CONTINUE
c
c
      eignfnv(1,1)    = -1.
      eignfnv(iim,1)  =  1.
      DO 20 i=1,imm1
      eignfnv(i+1,i+1)= -1.
      eignfnv(i,i+1)  =  1.
  20  CONTINUE
      DO 25 j=1,iim
      DO 25 i=1,iim
      eignfnv(i,j) = eignfnv(i,j)/(sddu(i)*sddv(j))
  25  CONTINUE
      DO 30 j=1,iim
      DO 30 i=1,iim
      eignfnu(i,j) = -eignfnv(j,i)
  30  CONTINUE
c
      DO j = 1, iim
      DO i = 1, iim
        vec (i,j) = 0.0
        vec1(i,j) = 0.0
       DO k = 1, iim
        vec (i,j) = vec(i,j)  + eignfnu(i,k) * eignfnv(k,j)
        vec1(i,j) = vec1(i,j) + eignfnv(i,k) * eignfnu(k,j)
       ENDDO
      ENDDO
      ENDDO

c
      CALL jacobi(vec,iim,iim,dv,eignfnv,nrot)
      CALL acc(eignfnv,d,iim)
      CALL eigen_sort(dv,eignfnv,iim,iim)
c
      CALL jacobi(vec1,iim,iim,du,eignfnu,nrot)
      CALL acc(eignfnu,d,iim)
      CALL eigen_sort(du,eignfnu,iim,iim)

cc   ancienne version avec appels IMSL
c
c     CALL MXM(eignfnu,iim,eignfnv,iim,vec,iim)
c     CALL MXM(eignfnv,iim,eignfnu,iim,vec1,iim)
c     CALL EVCSF(iim,vec,iim,dv,eignfnv,iim)
c     CALL acc(eignfnv,d,iim)
c     CALL eigen(eignfnv,dv)
c
c     CALL EVCSF(iim,vec1,iim,du,eignfnu,iim)
c     CALL acc(eignfnu,d,iim)
c     CALL eigen(eignfnu,du)

      RETURN
      END

