!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/dump2d.F,v 1.1.1.1 2004/05/19 12:53:05 lmdzadmin Exp $
!
      SUBROUTINE dump2d(im,jm,z,nom_z)
      IMPLICIT NONE
      INTEGER im,jm
      REAL z(im,jm)
      CHARACTER*80 nom_z

      INTEGER i,j,imin,illm,jmin,jllm
      REAL zmin,zllm

      PRINT*,nom_z

      zmin=z(1,1)
      zllm=z(1,1)
      imin=1
      illm=1
      jmin=1
      jllm=1

      DO j=1,jm
         DO i=1,im
            IF(z(i,j).GT.zllm) THEN
               illm=i
               jllm=j
               zllm=z(i,j)
            ENDIF
            IF(z(i,j).LT.zmin) THEN
               imin=i
               jmin=j
               zmin=z(i,j)
            ENDIF
         ENDDO
      ENDDO

      PRINT*,'MIN: ',zmin
      PRINT*,'MAX: ',zllm

      IF(zllm.GT.zmin) THEN
      DO j=1,jm
      WRITE(*,'(72i1)') (NINT(10.*(z(i,j)-zmin)/(zllm-zmin)),i=1,im)
      ENDDO
      ENDIF
      RETURN
      END
